#include <enoki/stl.h>
#include <mitsuba/core/plugin.h>
#include <mitsuba/core/properties.h>
#include <mitsuba/render/bsdf.h>
#include <mitsuba/render/fermatNEE.h>
#include <mitsuba/render/integrator.h>
#include <mitsuba/render/kdtree.h>
#include <mitsuba/render/medium.h>
#include <mitsuba/render/scene.h>
#include <tbb/parallel_reduce.h>

#if defined(MTS_ENABLE_EMBREE)
#include "scene_embree.inl"
#else
#include "scene_native.inl"
#endif

#if defined(MTS_ENABLE_OPTIX)
#include "scene_optix.inl"
#endif

NAMESPACE_BEGIN(mitsuba)

MTS_VARIANT Scene<Float, Spectrum>::Scene(const Properties &props) {
    for (auto &kv : props.objects()) {
        m_children.push_back(kv.second.get());

        Shape *shape           = dynamic_cast<Shape *>(kv.second.get());
        Emitter *emitter       = dynamic_cast<Emitter *>(kv.second.get());
        Sensor *sensor         = dynamic_cast<Sensor *>(kv.second.get());
        Integrator *integrator = dynamic_cast<Integrator *>(kv.second.get());

        if (shape) {
            if (shape->is_emitter())
                m_emitters.push_back(shape->emitter());
            if (shape->is_sensor())
                m_sensors.push_back(shape->sensor());

            m_bbox.expand(shape->bbox());
            m_shapes.push_back(shape);
        } else if (emitter) {
            // Surface emitters will be added to the list when attached to a
            // shape
            if (!has_flag(emitter->flags(), EmitterFlags::Surface))
                m_emitters.push_back(emitter);

            if (emitter->is_environment()) {
                if (m_environment)
                    Throw("Only one environment emitter can be specified per "
                          "scene.");
                m_environment = emitter;
            }
        } else if (sensor) {
            m_sensors.push_back(sensor);
        } else if (integrator) {
            if (m_integrator)
                Throw("Only one integrator can be specified per scene.");
            m_integrator = integrator;
        }
    }

    if (m_sensors.empty()) {
        Log(Warn, "No sensors found! Instantiating a perspective camera..");
        Properties sensor_props("perspective");
        sensor_props.set_float("fov", 45.0f);

        /* Create a perspective camera with a 45 deg. field of view
           and positioned so that it can see the entire scene */
        if (m_bbox.valid()) {
            ScalarPoint3f center   = m_bbox.center();
            ScalarVector3f extents = m_bbox.extents();

            ScalarFloat distance = hmax(extents) / (2.f * std::tan(45.f * .5f * math::Pi<ScalarFloat> / 180.f));

            sensor_props.set_float("far_clip", hmax(extents) * 5 + distance);
            sensor_props.set_float("near_clip", distance / 100);

            sensor_props.set_float("focus_distance", distance + extents.z() / 2);
            sensor_props.set_transform("to_world",
                                       ScalarTransform4f::translate(ScalarVector3f(center.x(), center.y(), m_bbox.min.z() - distance)));
        }

        m_sensors.push_back(PluginManager::instance()->create_object<Sensor>(sensor_props));
    }

    if (!m_integrator) {
        Log(Warn, "No integrator found! Instantiating a path tracer..");
        m_integrator = PluginManager::instance()->create_object<Integrator>(Properties("path"));
    }

    if constexpr (is_cuda_array_v<Float>)
        accel_init_gpu(props);
    else
        accel_init_cpu(props);

    // Create emitters' shapes (environment luminaires)
    for (Emitter *emitter : m_emitters)
        emitter->set_scene(this);

    for (Shape *shape : m_shapes) {
        if (shape->is_caustic_caster_single_scatter())
            m_caustic_casters_single.push_back(shape);
        if (shape->is_caustic_caster_multi_scatter())
            m_caustic_casters_multi.push_back(shape);
        if (shape->is_caustic_caster_double_refraction())
            m_caustic_casters_double_refraction.push_back(shape);
    }

    for (Emitter *emitter : m_emitters) {
        if (emitter->is_caustic_emitter_single_scatter())
            m_caustic_emitters_single.push_back(emitter);
        if (emitter->is_caustic_emitter_multi_scatter())
            m_caustic_emitters_multi.push_back(emitter);
    }
}

MTS_VARIANT Scene<Float, Spectrum>::~Scene() {
    if constexpr (is_cuda_array_v<Float>)
        accel_release_gpu();
    else
        accel_release_cpu();
}

MTS_VARIANT typename Scene<Float, Spectrum>::SurfaceInteraction3f Scene<Float, Spectrum>::ray_intersect(const Ray3f &ray,
                                                                                                        Mask active) const {
    MTS_MASKED_FUNCTION(ProfilerPhase::RayIntersect, active);

    if constexpr (is_cuda_array_v<Float>)
        return ray_intersect_gpu(ray, active);
    else
        return ray_intersect_cpu(ray, active);
}

MTS_VARIANT typename Scene<Float, Spectrum>::SurfaceInteraction3f Scene<Float, Spectrum>::ray_intersect_naive(const Ray3f &ray,
                                                                                                              Mask active) const {
    MTS_MASKED_FUNCTION(ProfilerPhase::RayIntersect, active);

#if !defined(MTS_ENABLE_EMBREE)
    if constexpr (!is_cuda_array_v<Float>)
        return ray_intersect_naive_cpu(ray, active);
#endif
    ENOKI_MARK_USED(ray);
    ENOKI_MARK_USED(active);
    NotImplementedError("ray_intersect_naive");
}

MTS_VARIANT typename Scene<Float, Spectrum>::Mask Scene<Float, Spectrum>::ray_test(const Ray3f &ray, Mask active) const {
    MTS_MASKED_FUNCTION(ProfilerPhase::RayTest, active);

    if constexpr (is_cuda_array_v<Float>)
        return ray_test_gpu(ray, active);
    else
        return ray_test_cpu(ray, active);
}

MTS_VARIANT
std::pair<typename Scene<Float, Spectrum>::DirectionSample3f, Spectrum>
Scene<Float, Spectrum>::sample_emitter_direction(const Interaction3f &ref, const Point2f &sample_, bool test_visibility,
                                                 Mask active) const {
    MTS_MASKED_FUNCTION(ProfilerPhase::SampleEmitterDirection, active);

    using EmitterPtr = replace_scalar_t<Float, Emitter *>;

    Point2f sample(sample_);
    DirectionSample3f ds;
    Spectrum spec;

    if (likely(!m_emitters.empty())) {
        if (m_emitters.size() == 1) {
            // Fast path if there is only one emitter
            std::tie(ds, spec) = m_emitters[0]->sample_direction(ref, sample, active);
        } else {
            ScalarFloat emitter_pdf = 1.f / m_emitters.size();

            // Randomly pick an emitter
            UInt32 index = min(UInt32(sample.x() * (ScalarFloat) m_emitters.size()), (uint32_t) m_emitters.size() - 1);

            // Rescale sample.x() to lie in [0,1) again
            sample.x() = (sample.x() - index * emitter_pdf) * m_emitters.size();

            EmitterPtr emitter = gather<EmitterPtr>(m_emitters.data(), index, active);

            // Sample a direction towards the emitter
            std::tie(ds, spec) = emitter->sample_direction(ref, sample, active);

            // Account for the discrete probability of sampling this emitter
            ds.pdf *= emitter_pdf;
            spec *= rcp(emitter_pdf);
        }

        active &= neq(ds.pdf, 0.f);

        // Perform a visibility test if requested
        if (test_visibility && any_or<true>(active)) {
            Ray3f ray(ref.p, ds.d, math::RayEpsilon<Float> * (1.f + hmax(abs(ref.p))), ds.dist * (1.f - math::ShadowEpsilon<Float>),
                      ref.time, ref.wavelengths);
            spec[ray_test(ray, active)] = 0.f;
        }
    } else {
        ds   = zero<DirectionSample3f>();
        spec = 0.f;
    }

    return { ds, spec };
}

// TODO direction sampling
// TODO move sampling in fermaNEE.cpp
MTS_VARIANT
std::pair<typename Scene<Float, Spectrum>::DirectionSample3f, Spectrum> Scene<Float, Spectrum>::sample_emitter_direction_double_refraction(
    FermatNEE<Float, Spectrum> *NEE_sampler, const BSDFContext &ctx, const SurfaceInteraction3f &ref, const Point2f &sample_,
    const Point2f &sample2_, const Float &sample3_, bool useSMS, bool *foundDoubleRefraction, ShapePtr *heightfield1,
    ShapePtr *heightfield2, Mask active) const {
    MTS_MASKED_FUNCTION(ProfilerPhase::SampleEmitterDirection, active);

    using EmitterPtr = replace_scalar_t<Float, Emitter *>;

    Point2f sample(sample_);
    Spectrum spec(1.0f);

    if (likely(!m_emitters.empty())) {

        // emitter sampling:
        // Uniformly sample an emitter, same as in Scene::sample_emitter_direction
        Float emitter_sample     = sample_.x();
        Float emitter_pdf        = 1.f / m_emitters.size();
        UInt32 index             = min(UInt32(emitter_sample * (ScalarFloat) m_emitters.size()), (uint32_t) m_emitters.size() - 1);
        const EmitterPtr emitter = gather<EmitterPtr>(m_emitters.data(), index);

        Spectrum ei_weight(0.f);
        DirectionSample3f ps;
        // surface
        if (emitter->shape()) {
            const ShapePtr shape = emitter->shape();
            ps                   = shape->sample_position(ref.time, sample_);

            SurfaceInteraction3f si_emitter;
            if (ps.pdf > 0) {
                si_emitter.p           = ps.p;
                si_emitter.wi          = Vector3f(0.f, 0.f, 1.f);
                si_emitter.wavelengths = ref.wavelengths;
                si_emitter.time        = ref.time;
                si_emitter.n           = ps.n;
                ei_weight              = rcp(ps.pdf) * rcp(emitter_pdf) * emitter->eval(si_emitter);
                ps.d = Vector3f(0); // point sample
            }

        }
        // environnement map
        else {
            std::tie(ps, spec) = emitter->sample_direction(ref, sample, active);
            if (ps.pdf > 0) {
                ei_weight = rcp(ps.pdf) * rcp(emitter_pdf) * spec;
            }
        }

        *foundDoubleRefraction = false;
        if (any_or<true>(active)) {
            // guided_path_AG test first
            Vector3f wo;
            if(ps.d == Vector3f(0)){
                wo = normalize(ps.p - ref.p);
            }else{
                wo = ps.d;
            }
            Ray3f ray(ref.p, wo, math::RayEpsilon<Float> * (1.f + hmax(abs(ref.p))),
                      norm(ps.p - ref.p) * (1.f - math::ShadowEpsilon<Float>), ref.time, ref.wavelengths);

            // check if 2 refractive heighfields are present:
            SurfaceInteraction3f si_heightfieldTest = ray_intersect(ray, active);
            // check for first heighfield
            if (si_heightfieldTest.shape) {
                if (si_heightfieldTest.shape->is_caustic_caster_double_refraction()) {
                    // check for 2nd heighfield:
                    auto [bs, bsdf_val] = si_heightfieldTest.bsdf()->sample(ctx, si_heightfieldTest, sample3_, sample2_, active);
                    Ray3f secondary_ray = si_heightfieldTest.spawn_ray(si_heightfieldTest.to_world(bs.wo));
                    SurfaceInteraction3f si_heightfieldTest_2 = ray_intersect(secondary_ray, active);

                    if (si_heightfieldTest_2.shape) {
                        if (si_heightfieldTest_2.shape->is_caustic_caster_double_refraction()) {
                            // both refractive interface found:
                            *foundDoubleRefraction = true;
                            if (useSMS) {
                                return { ps, 0 };
                            } else {
                                *heightfield1 = si_heightfieldTest.shape;
                                *heightfield2 = si_heightfieldTest_2.shape;
                                spec          = ei_weight;
                                return { ps, ei_weight };
                            }
                        }
                    }
                }
                if (*foundDoubleRefraction == false) {
                    return { ps, Spectrum(0.0f) };
                }
            }
        } else {
            return { ps, Spectrum(0.0f) };
        }

        return { ps, ei_weight };
    }
}

MTS_VARIANT
bool Scene<Float, Spectrum>::test_double_refraction_NEE(const BSDFContext &ctx, const SurfaceInteraction3f &ref, Vector3f wo, Point3f S, ShapePtr *H1, ShapePtr *H2, Mask active, const Point2f &sample2_, const Float &sample3_ ) const
{

  MTS_MASKED_FUNCTION(ProfilerPhase::SampleEmitterDirection, active);
  if (any_or<true>(active)) {
    // guided_path_AG test first
    Vector3f wo;
      wo = normalize(S - ref.p);
    Ray3f ray(ref.p, wo, math::RayEpsilon<Float> * (1.f + hmax(abs(ref.p))),
	norm(S - ref.p) * (1.f - math::ShadowEpsilon<Float>), ref.time, ref.wavelengths);

    // check if 2 refractive heighfields are present:
    SurfaceInteraction3f si_heightfieldTest = ray_intersect(ray, active);
    // check for first heighfield
    if (si_heightfieldTest.shape) {
      if (si_heightfieldTest.shape->is_caustic_caster_double_refraction()) {
	// check for 2nd heighfield:
	auto [bs, bsdf_val] = si_heightfieldTest.bsdf()->sample(ctx, si_heightfieldTest, sample3_, sample2_, active);
	Ray3f secondary_ray = si_heightfieldTest.spawn_ray(si_heightfieldTest.to_world(bs.wo));
	SurfaceInteraction3f si_heightfieldTest_2 = ray_intersect(secondary_ray, active);

	if (si_heightfieldTest_2.shape) {
	  if (si_heightfieldTest_2.shape->is_caustic_caster_double_refraction()) {
	    // both refractive interface found:
	      *H1 = si_heightfieldTest.shape;
	      *H2 = si_heightfieldTest_2.shape;
	      return true;
	    }
	  }
	}
      }
    }

  return false;
}



MTS_VARIANT Float Scene<Float, Spectrum>::pdf_emitter_direction(const Interaction3f &ref, const DirectionSample3f &ds, Mask active) const {
    MTS_MASK_ARGUMENT(active);
    using EmitterPtr = replace_scalar_t<Float, const Emitter *>;

    if (m_emitters.size() == 1) {
        // Fast path if there is only one emitter
        return m_emitters[0]->pdf_direction(ref, ds, active);
    } else {
        return reinterpret_array<EmitterPtr>(ds.object)->pdf_direction(ref, ds, active) * (1.f / m_emitters.size());
    }
}

MTS_VARIANT void Scene<Float, Spectrum>::traverse(TraversalCallback *callback) {
    for (auto &child : m_children) {
        std::string id = child->id();
        if (id.empty() || string::starts_with(id, "_unnamed_"))
            id = child->class_()->name();
        callback->put_object(id, child.get());
    }
}

MTS_VARIANT void Scene<Float, Spectrum>::parameters_changed(const std::vector<std::string> & /*keys*/) {
    if (m_environment)
        m_environment->set_scene(this); // TODO use
                                        // parameters_changed({"scene"})
}

MTS_VARIANT std::string Scene<Float, Spectrum>::to_string() const {
    std::ostringstream oss;
    oss << "Scene[" << std::endl << "  children = [" << std::endl;
    for (size_t i = 0; i < m_children.size(); ++i) {
        oss << "    " << string::indent(m_children[i], 4);
        if (i + 1 < m_children.size())
            oss << ",";
        oss << std::endl;
    }
    oss << "  ]" << std::endl << "]";
    return oss.str();
}

void librender_nop() {}

MTS_IMPLEMENT_CLASS_VARIANT(Scene, Object, "scene")
MTS_INSTANTIATE_CLASS(Scene)
NAMESPACE_END(mitsuba)
