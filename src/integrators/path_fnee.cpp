#include <enoki/stl.h>
#include <mitsuba/core/properties.h>
#include <mitsuba/core/ray.h>
#include <mitsuba/render/bsdf.h>
#include <mitsuba/render/emitter.h>
#include <mitsuba/render/integrator.h>
#include <mitsuba/render/records.h>
#include <random>

#include <mitsuba/render/OpenGL_viewer_client.h>
#include <mitsuba/render/fermatNEE.h>
#include <mitsuba/render/manifold_ms_AG.h>
#include <mitsuba/render/manifold_ss_AG.h>

NAMESPACE_BEGIN(mitsuba)

/**!

.. _integrator-path:

Path tracer (:monosp:`path`)
-------------------------------------------

.. pluginparameters::

 * - max_depth
   - |int|
   - Specifies the longest path depth in the generated output image (where -1
corresponds to :math:`\infty`). A value of 1 will only render directly visible
light sources. 2 will lead to single-bounce (direct-only) illumination, and so
on. (Default: -1)
 * - rr_depth
   - |int|
   - Specifies the minimum path depth, after which the implementation will start
to use the *russian roulette* path termination criterion. (Default: 5)
 * - hide_emitters
   - |bool|
   - Hide directly visible emitters. (Default: no, i.e. |false|)

This integrator implements a basic path tracer and is a **good default choice**
when there is no strong reason to prefer another method.

To use the path tracer appropriately, it is instructive to know roughly how
it works: its main operation is to trace many light paths using *random walks*
starting from the sensor. A single random walk is shown below, which entails
casting a ray associated with a pixel in the output image and searching for
the first visible intersection. A new direction is then chosen at the
intersection, and the ray-casting step repeats over and over again (until one of
several stopping criteria applies).

.. image:: ../images/integrator_path_figure.png
    :width: 95%
    :align: center

At every intersection, the path tracer tries to create a connection to
the light source in an attempt to find a *complete* path along which
light can flow from the emitter to the sensor. This of course only works
when there is no occluding object between the intersection and the emitter.

This directly translates into a category of scenes where
a path tracer can be expected to produce reasonable results: this is the case
when the emitters are easily "accessible" by the contents of the scene. For
instance, an interior scene that is lit by an area light will be considerably
harder to render when this area light is inside a glass enclosure (which
effectively counts as an occluder).

Like the :ref:`direct <integrator-direct>` plugin, the path tracer internally
relies on multiple importance sampling to combine BSDF and emitter samples. The
main difference in comparison to the former plugin is that it considers light
paths of arbitrary length to compute both direct and indirect illumination.

.. _sec-path-strictnormals:

.. Commented out for now
.. Strict normals
   --------------

.. Triangle meshes often rely on interpolated shading normals
   to suppress the inherently faceted appearance of the underlying geometry.
These "fake" normals are not without problems, however. They can lead to
paradoxical situations where a light ray impinges on an object from a direction
that is classified as "outside" according to the shading normal, and "inside"
according to the true geometric normal.

.. The :paramtype:`strict_normals` parameter specifies the intended behavior
when such cases arise. The default (|false|, i.e. "carry on") gives precedence
to information given by the shading normal and considers such light paths to be
valid. This can theoretically cause light "leaks" through boundaries, but it is
not much of a problem in practice.

.. When set to |true|, the path tracer detects inconsistencies and ignores these
paths. When objects are poorly tesselated, this latter option may cause them to
lose a significant amount of the incident radiation (or, in other words, they
will look dark).

.. note:: This integrator does not handle participating media

 */

template <typename Float, typename Spectrum>
class path_fnee : public MonteCarloIntegrator<Float, Spectrum> {
public:
    MTS_IMPORT_BASE(MonteCarloIntegrator, m_max_depth, m_rr_depth)
    MTS_IMPORT_TYPES(Scene, Sampler, Medium, Emitter, EmitterPtr, BSDF, BSDFPtr,
                     ShapePtr)
    using SpecularManifold = SpecularManifold<Float, Spectrum>;
    using SpecularManifoldMultiScatter_AG =
        SpecularManifoldMultiScatter_AG<Float, Spectrum>;
    using SpecularManifoldSingleScatter_AG =
        SpecularManifoldSingleScatter_AG<Float, Spectrum>;
    using FermatNEE = FermatNEE<Float, Spectrum>;

    static inline ThreadLocal<SpecularManifoldMultiScatter_AG> tl_manifold{};
    static inline ThreadLocal<FermatNEE> tl_fermaNEE{};

    path_fnee(const Properties &props) : Base(props) {
        m_sms_config = SMSConfig();
        m_sms_config.biased = props.bool_("biased", false);
        m_sms_config.twostage = props.bool_("twostage", false);
        m_sms_config.halfvector_constraints =
            props.bool_("halfvector_constraints", false);
        m_sms_config.mnee_init = props.bool_("mnee_init", false);
        m_sms_config.step_scale = props.float_("step_scale", 1.f);
        m_sms_config.max_iterations = props.int_("max_iterations", 20);
        m_sms_config.solver_threshold = props.float_("solver_threshold", 1e-5f);
        m_sms_config.uniqueness_threshold =
            props.float_("uniqueness_threshold", 1e-4f);
        m_sms_config.max_trials = props.int_("max_trials", -1);

        m_sms_config.bounces = props.int_("bounces", 2);

        alpha = props.float_("alpha", 0.5f);
        beta = props.float_("beta", 0.8f);

        crop_caustic = props.bool_("crop_caustic", false);
        use_SMS = props.bool_("use_SMS", false);
    }

    void Specular_nee(const Scene *scene, Sampler *sampler, BSDFContext &ctx,
                      SurfaceInteraction3f &si, Spectrum &result,
                      FermatNEE &NEE_sampler, Spectrum &throughput,
                      Mask active) const {

        ShapePtr H1;
        ShapePtr H2;
        Spectrum NEE_specularOutput;

        typename SpecularManifold::EmitterInteraction vy =
            SpecularManifold::sample_emitter_interaction(
                si, scene->caustic_emitters_multi_scatter(), sampler);

        /* +--------------------+
        // |// double refraction|
        // +--------------------+ */
        bool sample_heighfield_success = scene->test_double_refraction_NEE(
            ctx, si, Vector3f(0), vy.p, &H1, &H2, active, sampler->next_2d(),
            sampler->next_1d());

        if (sample_heighfield_success) {
            Point3f origin_H1 =
                H1->get_world_transform().transform_affine(Point3f(0.f));
            Point3f origin_H2 =
                H2->get_world_transform().transform_affine(Point3f(0.f));
            Point3f origin_H2_local_H1 =
                H1->get_local_transform().transform_affine(origin_H2);
            Point3f O_local_H1 =
                H1->get_local_transform().transform_affine(si.p);
            Float e = origin_H2_local_H1.x();

            if (squared_norm(O_local_H1) >
                squared_norm(O_local_H1 - origin_H2_local_H1)) {
                std::swap(H1, H2);
                e = -e;
            }
            Float NEE_invpdf = 1.0f;
            Point3f S = H1->get_local_transform().transform_affine(vy.p);
            Point3f O = H1->get_local_transform().transform_affine(si.p);
            // std::cout<< Vector3f(1,1,1) << "|" <<
            // H1->get_local_transform()*(Vector3f(1,1,1))
            // <<std::endl;
            Vector3f wo_r = normalize(vy.p - si.p);
            RayDifferential3f r(
                si.p, wo_r, math::RayEpsilon<Float> * (1.f + hmax(abs(si.p))),
                norm(vy.p - si.p) * (1.f - math::ShadowEpsilon<Float>), si.time,
                si.wavelengths);

            // find ray proposal
            Float partition_size = H1->eval_attribute_1("partition_size", si);
            Float inv_partition_size =
                H1->eval_attribute_1("inv_partition_size", si);
            Float inv_partition_area =
                H1->eval_attribute_1("inv_partition_area", si);
            typename FermatNEE::L_data_float fermat_connection_data = {
                H1,
                H2,
                H1->eval_attribute_1("H", si),
                H2->eval_attribute_1("H", si),
                H1->eval_attribute_1("L", si),
                H2->eval_attribute_1("L", si),
                H1->eval_attribute_1("W", si),
                H2->eval_attribute_1("W", si),
                O,
                S,
                1.00f,
                1.54f,
                e,
                partition_size,
                inv_partition_size,
                inv_partition_area
            }; // needs local coordinates

            auto [FNEE_success, fermat_connection_direction] =
                NEE_sampler.specular_connection(si, vy, &fermat_connection_data);
            r.d = fermat_connection_direction;
            r.update();

            if (FNEE_success) {
                Vector3f wo_bsdf = (si.to_local(r.d));

                Spectrum bsdf_shading_point = si.bsdf()->eval(ctx, si, wo_bsdf);
                bsdf_shading_point =
                    si.to_world_mueller(bsdf_shading_point, -wo_bsdf, si.wi);

                NEE_specularOutput = NEE_sampler.compute_ray_contribution(
                    si, ctx, r, vy, 1.0f, 1.54f, H1, H2, active);

                if (NEE_specularOutput != Spectrum(0)) {
                    if (!use_SMS) {
                        NEE_invpdf = NEE_sampler.eval_invPDF(
                            &fermat_connection_data,
                            fermat_connection_direction);

                    } else {
                        NEE_invpdf = NEE_sampler.SMS_eval_invPDF(
                            si, vy, &fermat_connection_data,
                            fermat_connection_direction);
                    }

                    result[active] += throughput * bsdf_shading_point *
                                      vy.weight * NEE_specularOutput *
                                      NEE_invpdf;
                }
            }
        }

        // /* +-------------+
        // // |// reflection|
        // // +-------------+ */
        // H1                        = nullptr;
        // H2                        = nullptr;
        // Float invPdf_H_sampling = 0;
        // sample_heighfield_success =
        // NEE_sampler.sample_heightfield(&H1, &invPdf_H_sampling,
        // active); if (sample_heighfield_success) {
        //     // TODO gerer direction sample
        //     typename SpecularManifold::EmitterInteraction vy =
        //     SpecularManifold::sample_emitter_interaction(si,
        //     scene->caustic_emitters_multi_scatter(), sampler);
        //     Float NEE_invpdf                      = 1.0f;
        //     Point3f S                             =
        //     H1->get_local_transform().transform_affine(vy.p);
        //     Point3f O                             =
        //     H1->get_local_transform().transform_affine(si.p);

        //     Vector3f wo_r = normalize(vy.p - si.p);
        //     // Ray3f r(si.p, wo_r, math::RayEpsilon<Float> * (1.f
        //     + hmax(abs(si.p))),
        //     //                     math::Infinity<Float>,
        //     si.time, si.wavelengths); Ray3f r; r.o    = si.p;
        //     r.mint = math::RayEpsilon<Float> * (1.f +
        //     hmax(abs(si.p))); r.time = si.time;

        //     // find ray proposal
        //     typename FermatNEE::L_data_float
        //     fermat_connection_data = { H1,
        //                                                                 nullptr,
        //                                                                 H1->eval_attribute_1("H",
        //                                                                 si),
        //                                                                 0.f,
        //                                                                 H1->eval_attribute_1("L",
        //                                                                 si),
        //                                                                 0.f,
        //                                                                 H1->eval_attribute_1("W",
        //                                                                 si),
        //                                                                 0.f,
        //                                                                 O,
        //                                                                 S,
        //                                                                 1.00f,
        //                                                                 0.0f,
        //                                                                 0.0f
        //                                                                 };

        //     Vector3f fermat_connection_direction;
        //     bool FNEE_success   = false;
        //     bool use_two_stages = false;

        //     FNEE_success =
        //     NEE_sampler.fermat_connection_reflection(&fermat_connection_data,
        //     &fermat_connection_direction); if (FNEE_success) {

        //         r.d =
        //         H1->get_world_transform().transform_affine(fermat_connection_direction);
        //         r.d = normalize(r.d);
        //         r.update();

        //         Vector3f wo_bsdf            = (si.to_local(r.d));
        //         Spectrum bsdf_shading_point =
        //         si.bsdf()->eval(ctx, si, wo_bsdf);
        //         bsdf_shading_point          =
        //         si.to_world_mueller(bsdf_shading_point, -wo_bsdf,
        //         si.wi);

        //         NEE_specularOutput =
        //             NEE_sampler.compute_ray_contribution_reflection(si,
        //             ctx, r, vy, 1.0f, 1.54f, H1, nullptr,
        //             active);

        //         if (NEE_specularOutput != Spectrum(0)) {

        //             NEE_invpdf     =
        //             NEE_sampler.eval_invPDF_reflection(&fermat_connection_data,
        //             fermat_connection_direction); Vector3f wo =
        //             (si.to_local(r.d)); Float bsdf_pdf =
        //             bsdf->pdf(ctx, si, wo, active); Float mis =
        //             mis_weight(vy.pdf, bsdf_pdf);

        //             result[active] += throughput *
        //             bsdf_shading_point * vy.weight *
        //             NEE_specularOutput * NEE_invpdf *
        //             invPdf_H_sampling; // *NEE_invpdf
        //         }
        //     }
        // }
    }

    std::pair<Spectrum, Mask> sample(const Scene *scene, Sampler *sampler,
                                     const RayDifferential3f &ray_,
                                     const Medium * /* medium */,
                                     Float * /* aovs */,
                                     Mask active) const override {
        MTS_MASKED_FUNCTION(ProfilerPhase::SamplingIntegratorSample, active);
        // auto &mf_ms       = (SpecularManifoldMultiScatter_AG &) tl_manifold;
        FermatNEE &NEE_sampler = (FermatNEE &) tl_fermaNEE;
        // mf_ms.init(scene, m_sms_config);
        NEE_sampler.init(scene, sampler, m_sms_config.uniqueness_threshold,
                         m_sms_config.max_trials, alpha, beta, crop_caustic,
                         m_sms_config, use_SMS);

        RayDifferential3f ray = ray_;
        // Tracks radiance scaling due to index of refraction changes
        Float eta(1.f);

        // MIS weight for intersected emitters (set by prev. iteration)
        Float emission_weight(1.f);

        Spectrum throughput(1.f), result(0.f);
        bool specular_camera_path =
            true; // To capture emitters visible direcly through purely specular
                  // reflection/refractions

        // ---------------------- First intersection ----------------------

        SurfaceInteraction3f si = scene->ray_intersect(ray);
        Mask valid_ray = si.is_valid();
        EmitterPtr emitter = si.emitter(scene);

        if (emitter) {
            result += emitter->eval(si);
        }

        for (int depth = 1;; ++depth) {

            // ------------------ Possibly terminate path -----------------

            if (!si.is_valid())
                break;
            si.compute_partials(ray);

            if (depth > m_rr_depth) {
                Float q = min(hmax(depolarize(throughput)) * sqr(eta), .95f);
                if (sampler->next_1d() > q)
                    break;
                throughput *= rcp(q);
            }

            if (uint32_t(depth) >= uint32_t(m_max_depth))
                break;

            BSDFContext ctx;
            ctx.sampler = sampler;
            BSDFPtr bsdf = si.bsdf(ray);

            // TODO have two integrator: <launch solver for shadow only> /
            // <global>
            bool doubleRefraction = false;
            bool on_caustic_receiver = si.shape->is_caustic_receiver();
            bool on_caustic_caster =
                si.shape->is_caustic_caster_single_scatter() ||
                si.shape->is_caustic_caster_multi_scatter() ||
                si.shape->is_caustic_bouncer();

            if (on_caustic_receiver && !on_caustic_caster) {
                Specular_nee(scene, sampler, ctx, si, result, NEE_sampler,
                             throughput, active);
            }

            // --------------------- Emitter sampling ---------------------

            Mask active_e =
                active && has_flag(bsdf->flags(), BSDFFlags::Smooth);

            /* As usual, emitter sampling only makes sense on Smooth BSDFs
               that can be evaluated.
               Additionally, filter out:
                - paths that we could previously sample with SMS
                - paths that are even harder to sample, e.g. paths bouncing
                  off several caustic casters before hitting the light.
               As a result, we only do emitter sampling on non-caustic
               casters---with the exception of the first bounce where we might
               see a direct (glossy) reflection of a light source this way.

               Note: of course, SMS might not always be the optimal sampling
               strategy. For example, when rough surfaces are involved it
               would be still better to do emitter sampling.
               A way of incoorporating MIS with all of this would be super
               useful. */
            if (has_flag(bsdf->flags(), BSDFFlags::Smooth) &&
                !on_caustic_caster) {
                // if (NEE_specularOutput == Spectrum(0.f)) {
                auto [ds, emitter_val] = scene->sample_emitter_direction(
                    si, sampler->next_2d(active_e), true, active_e);

                active_e &= neq(ds.pdf, 0.f);

                // Query the BSDF for that emitter-sampled direction
                Vector3f wo = si.to_local(ds.d);
                Spectrum bsdf_val = bsdf->eval(ctx, si, wo, active_e);
                bsdf_val = si.to_world_mueller(bsdf_val, -wo, si.wi);

                // Determine density of sampling that same direction using BSDF
                // sampling
                Float bsdf_pdf = bsdf->pdf(ctx, si, wo, active_e);
                Float mis = select(ds.delta, 1.f, mis_weight(ds.pdf, bsdf_pdf));
                result[active_e] += mis * throughput * bsdf_val * emitter_val;
            }

            // ----------------------- BSDF sampling ----------------------

            // Sample BSDF * cos(theta)
            auto [bs, bsdf_weight] =
                bsdf->sample(ctx, si, sampler->next_1d(), sampler->next_2d());
            bsdf_weight = si.to_world_mueller(bsdf_weight, -bs.wo, si.wi);

            throughput = throughput * bsdf_weight;
            eta *= bs.eta;
            if (!has_flag(bs.sampled_type, BSDFFlags::Delta)) {
                specular_camera_path = false;
            }

            if (all(eq(throughput, 0.f)))
                break;

            // Intersect the BSDF ray against the scene geometry
            ray = si.spawn_ray(si.to_world(bs.wo));
            SurfaceInteraction3f si_bsdf = scene->ray_intersect(ray);
            emitter = si_bsdf.emitter(scene);

            // Hit emitter after BSDF sampling
            if (emitter) {
                /* With the same reasoning as in the emitter sampling case,
                   filter out some of the light paths here.
                   Again, this is unfortunately not robust in all cases,
                   for large light sources, BSDF sampling would be more
                   appropriate than relying purely on SMS. */
                if (!on_caustic_caster || specular_camera_path) {
                    /* Only do BSDF sampling in usual way if we don't interact
                       with a caustic caster now. */

                    // Evaluate the emitter for that direction
                    Spectrum emitter_val = emitter->eval(si_bsdf);

                    /* Determine probability of having sampled that same
                       direction using emitter sampling. */
                    DirectionSample3f ds(si_bsdf, si);
                    ds.object = emitter;
                    Float emitter_pdf =
                        select(!has_flag(bs.sampled_type, BSDFFlags::Delta),
                               scene->pdf_emitter_direction(si, ds), 0.f);
                    Float mis = mis_weight(bs.pdf, emitter_pdf);
                    result += mis * throughput * emitter_val;
                }
            }

            si = std::move(si_bsdf);
        }

        return { result, valid_ray };
    }

    //! @}
    // =============================================================

    std::string to_string() const override {
        return tfm::format("path_fnee[\n"
                           "  max_depth = %i,\n"
                           "  rr_depth = %i\n"
                           "]",
                           m_max_depth, m_rr_depth);
    }

    Float mis_weight(Float pdf_a, Float pdf_b) const {
        pdf_a *= pdf_a;
        pdf_b *= pdf_b;
        return select(pdf_a > 0.f, pdf_a / (pdf_a + pdf_b), 0.f);
    }
    MTS_DECLARE_CLASS()
protected:
    SMSConfig m_sms_config;
    float alpha, beta;
    bool init = false;
    bool crop_caustic;
    bool use_SMS;
};

MTS_IMPLEMENT_CLASS_VARIANT(path_fnee, MonteCarloIntegrator)
MTS_EXPORT_PLUGIN(path_fnee, "Path Tracer integrator");
NAMESPACE_END(mitsuba)
