#include <mitsuba/render/bsdf.h>
#include <mitsuba/core/properties.h>
#include <mitsuba/render/shape.h>

NAMESPACE_BEGIN(mitsuba)

MTS_VARIANT BSDF<Float, Spectrum>::BSDF(const Properties &props)
    : m_flags(+BSDFFlags::None), m_id(props.id()) { }

MTS_VARIANT BSDF<Float, Spectrum>::~BSDF() { }

MTS_VARIANT Spectrum BSDF<Float, Spectrum>::eval_null_transmission(
    const SurfaceInteraction3f & /* si */, Mask /* active */) const {
    return 0.f;
}

MTS_VARIANT typename BSDF<Float, Spectrum>::Frame3f
BSDF<Float, Spectrum>::frame(const SurfaceInteraction3f &si,
                             Float /* smoothing */,
                             Mask /* active */) const {
    return compute_shading_frame(si.sh_frame.n, si.dp_du);
}

MTS_VARIANT std::pair<typename BSDF<Float, Spectrum>::Frame3f,
                      typename BSDF<Float, Spectrum>::Frame3f>
BSDF<Float, Spectrum>::frame_derivative(const SurfaceInteraction3f &si,
                                        Float /* smoothing */,
                                        Mask /* active */) const {
    auto [dn_du, dn_dv] = (si.instance ? si.instance : si.shape)->normal_derivative(si);
        // std::cout << "frame_derivative" << dn_du << dn_dv << std::endl;


    return compute_shading_frame_derivative(si.sh_frame.n, si.dp_du, dn_du, dn_dv);
}

MTS_VARIANT std::pair<typename BSDF<Float, Spectrum>::Point2f,
                      typename BSDF<Float, Spectrum>::Matrix2f>
BSDF<Float, Spectrum>::lean(const SurfaceInteraction3f & /* si */,
                            Mask /* active */) const {
    return { 0.f, 0.f };
}

MTS_VARIANT typename BSDF<Float, Spectrum>::Point2f
BSDF<Float, Spectrum>::slope(const Point2f & /* uv */,
                             Mask /* active */) const {
    return Point2f(0.f);
}

MTS_VARIANT std::pair<typename BSDF<Float, Spectrum>::Vector2f,
                      typename BSDF<Float, Spectrum>::Vector2f>
BSDF<Float, Spectrum>::slope_derivative(const Point2f & /* uv */,
                                        Mask /* active */) const {
    return { Vector2f(0.f), Vector2f(0.f) };
}

MTS_VARIANT Float BSDF<Float, Spectrum>::roughness() const {
    return math::Infinity<Float>;
}

MTS_VARIANT Complex<Spectrum> BSDF<Float, Spectrum>::ior(const SurfaceInteraction3f & /* si */,
                                                         Mask /* active */) const {
    return {0.f, 1.f};
}

MTS_VARIANT Float BSDF<Float, Spectrum>::glint_component_weight(const SurfaceInteraction3f & /* si */,
                                                                Mask /* active */) const {
    return 1.f;
}

MTS_VARIANT std::tuple<const Float *, size_t, Float>
BSDF<Float, Spectrum>::lean_pointer() const {
    Throw("Vectorized glints is trying to access LEAN pointer for non-normalmap BSDF!");
    return std::make_tuple(nullptr, 0, 0.f);
}

MTS_VARIANT std::string BSDF<Float, Spectrum>::id() const { return m_id; }

template <typename Index>
std::string type_mask_to_string(Index type_mask) {
    std::ostringstream oss;
    oss << "{ ";

#define is_set(flag) has_flag(type_mask, flag)
    if (is_set(BSDFFlags::All)) {
        oss << "all ";
        type_mask = type_mask & ~BSDFFlags::All;
    }
    if (is_set(BSDFFlags::Reflection)) {
        oss << "reflection ";
        type_mask = type_mask & ~BSDFFlags::Reflection;
    }
    if (is_set(BSDFFlags::Transmission)) {
        oss << "transmission ";
        type_mask = type_mask & ~BSDFFlags::Transmission;
    }
    if (is_set(BSDFFlags::Smooth)) {
        oss << "smooth ";
        type_mask = type_mask & ~BSDFFlags::Smooth;
    }
    if (is_set(BSDFFlags::Diffuse)) {
        oss << "diffuse ";
        type_mask = type_mask & ~BSDFFlags::Diffuse;
    }
    if (is_set(BSDFFlags::Glossy)) {
        oss << "glossy ";
        type_mask = type_mask & ~BSDFFlags::Glossy;
    }
    if (is_set(BSDFFlags::Delta)) {
        oss << "delta";
        type_mask = type_mask & ~BSDFFlags::Delta;
    }
    if (is_set(BSDFFlags::Delta1D)) {
        oss << "delta_1d ";
        type_mask = type_mask & ~BSDFFlags::Delta1D;
    }
    if (is_set(BSDFFlags::DiffuseReflection)) {
        oss << "diffuse_reflection ";
        type_mask = type_mask & ~BSDFFlags::DiffuseReflection;
    }
    if (is_set(BSDFFlags::DiffuseTransmission)) {
        oss << "diffuse_transmission ";
        type_mask = type_mask & ~BSDFFlags::DiffuseTransmission;
    }
    if (is_set(BSDFFlags::GlossyReflection)) {
        oss << "glossy_reflection ";
        type_mask = type_mask & ~BSDFFlags::GlossyReflection;
    }
    if (is_set(BSDFFlags::GlossyTransmission)) {
        oss << "glossy_transmission ";
        type_mask = type_mask & ~BSDFFlags::GlossyTransmission;
    }
    if (is_set(BSDFFlags::DeltaReflection)) {
        oss << "delta_reflection ";
        type_mask = type_mask & ~BSDFFlags::DeltaReflection;
    }
    if (is_set(BSDFFlags::DeltaTransmission)) {
        oss << "delta_transmission ";
        type_mask = type_mask & ~BSDFFlags::DeltaTransmission;
    }
    if (is_set(BSDFFlags::Delta1DReflection)) {
        oss << "delta_1d_reflection ";
        type_mask = type_mask & ~BSDFFlags::Delta1DReflection;
    }
    if (is_set(BSDFFlags::Delta1DTransmission)) {
        oss << "delta_1d_transmission ";
        type_mask = type_mask & ~BSDFFlags::Delta1DTransmission;
    }
    if (is_set(BSDFFlags::Null)) {
        oss << "null ";
        type_mask = type_mask & ~BSDFFlags::Null;
    }
    if (is_set(BSDFFlags::Anisotropic)) {
        oss << "anisotropic ";
        type_mask = type_mask & ~BSDFFlags::Anisotropic;
    }
    if (is_set(BSDFFlags::FrontSide)) {
        oss << "front_side ";
        type_mask = type_mask & ~BSDFFlags::FrontSide;
    }
    if (is_set(BSDFFlags::BackSide)) {
        oss << "back_side ";
        type_mask = type_mask & ~BSDFFlags::BackSide;
    }
    if (is_set(BSDFFlags::SpatiallyVarying)) {
        oss << "spatially_varying ";
        type_mask = type_mask & ~BSDFFlags::SpatiallyVarying;
    }
    if (is_set(BSDFFlags::NonSymmetric)) {
        oss << "non_symmetric ";
        type_mask = type_mask & ~BSDFFlags::NonSymmetric;
    }
#undef is_set

    Assert(type_mask == 0);
    oss << "}";
    return oss.str();
}

std::ostream &operator<<(std::ostream &os, const BSDFContext& ctx) {
    os << "BSDFContext[" << std::endl
        << "  mode = " << ctx.mode << "," << std::endl
        << "  type_mask = " << type_mask_to_string(ctx.type_mask) << "," << std::endl
        << "  component = ";
    if (ctx.component == (uint32_t) -1)
        os << "all";
    else
        os << ctx.component;
    os  << std::endl << "]";
    return os;
}

std::ostream &operator<<(std::ostream &os, const TransportMode &mode) {
    switch (mode) {
        case TransportMode::Radiance:   os << "radiance"; break;
        case TransportMode::Importance: os << "importance"; break;
        default:                        os << "invalid"; break;
    }
    return os;
}

MTS_IMPLEMENT_CLASS_VARIANT(BSDF, Object, "bsdf")
MTS_INSTANTIATE_CLASS(BSDF)
NAMESPACE_END(mitsuba)
