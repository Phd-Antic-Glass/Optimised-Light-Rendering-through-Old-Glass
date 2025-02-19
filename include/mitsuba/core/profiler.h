/* Builtin sampling profiler based on that of PBRTv3 */

#pragma once

#include <mitsuba/core/object.h>

#if !defined(MTS_PROFILE_HASH_SIZE)
#define MTS_PROFILE_HASH_SIZE 256
#endif

NAMESPACE_BEGIN(mitsuba)

/**
 * List of 'phases' that are handled by the profiler. Note that a partial order
 * is assumed -- if a method "B" can occur in a call graph of another method
 * "A", then "B" must occur after "A" in the list below.
 */
enum class ProfilerPhase : int {
    InitScene = 0,            /* Scene initialization */
    LoadGeometry,             /* Geometry loading */
    LoadTexture,              /* Texture loading */
    InitKDTree,               /* kd-tree construction */
    Render,                   /* Integrator::render() */
    SamplingIntegratorSample, /* SamplingIntegrator::sample() */
    SMSCaustics,
    SMSCausticsBernoulliTrials,
    SMSGlints,
    SMSGlintsBernoulliTrials,
    SampleEmitterRay,         /* Scene::sample_emitter_ray() */
    SampleEmitterDirection,   /* Scene::sample_emitter_direction() */
    RayTest,                  /* Scene::ray_test() */
    RayIntersect,             /* Scene::ray_intersect() */
    CreateSurfaceInteraction, /* KDTree::create_surface_interaction() */
    ImageBlockPut,            /* ImageBlock::put() */
    BSDFEvaluate,             /* BSDF::eval() and BSDF::pdf() */
    BSDFSample,               /* BSDF::sample() */
    PhaseFunctionEvaluate,    /* PhaseFunction::eval() and PhaseFunction::pdf() */
    PhaseFunctionSample,      /* PhaseFunction::sample() */
    MediumEvaluate,           /* Medium::eval() and Medium::pdf() */
    MediumSample,             /* Medium::sample() */
    EndpointEvaluate,         /* Endpoint::eval() and Endpoint::pdf() */
    EndpointSampleRay,        /* Endpoint::sample_ray() */
    EndpointSampleDirection,  /* Endpoint::sample_direction() */
    TextureSample,            /* Texture::sample() */
    TextureEvaluate,          /* Texture::eval() and Texture::pdf() */
    SMS_NewtonSolver,
    SMS_computeDirection,
    FermaNEE_principal,
    FermaNEE_NewtonSolver,
    FermaNEE_computeDirection,
    Eval_heighfield,

    ProfilerPhaseCount
};

constexpr const char *profiler_phase_id[int(ProfilerPhase::ProfilerPhaseCount)] = { "Scene initialization",
                                                                                    "Geometry loading",
                                                                                    "Texture loading",
                                                                                    "kd-tree construction",
                                                                                    "Integrator::render()",
                                                                                    "SamplingIntegrator::sample()",
                                                                                    "SMS::Caustics",
                                                                                    "SMS::CausticsBernoulliTrials",
                                                                                    "SMS::Glints",
                                                                                    "SMS::GlintsBernoulliTrials",
                                                                                    "Scene::sample_emitter_ray()",
                                                                                    "Scene::sample_emitter_direction()",
                                                                                    "Scene::ray_test()",
                                                                                    "Scene::ray_intersect()",
                                                                                    "KDTree::create_surface_interaction()",
                                                                                    "ImageBlock::put()",
                                                                                    "BSDF::eval(), pdf()",
                                                                                    "BSDF::sample()",
                                                                                    "PhaseFunction::eval(), pdf()",
                                                                                    "PhaseFunction::sample()",
                                                                                    "Medium::eval(), pdf()",
                                                                                    "Medium::sample()",
                                                                                    "Endpoint::eval(), pdf()",
                                                                                    "Endpoint::sample_ray()",
                                                                                    "Endpoint::sample_direction()",
                                                                                    "Texture::sample()",
                                                                                    "Texture::eval()",
                                                                                    "SMS_NewtonSolver",
                                                                                    "SMS_computeDirection",
                                                                                    "FermaNEE_principal",
                                                                                    "FermaNEE_NewtonSolver",
                                                                                    "FermaNEE_fonction",
                                                                                    "Eval_heighfield"
                                                                                     };

static_assert(int(ProfilerPhase::ProfilerPhaseCount) <= 64, "List of profiler phases is limited to 64 entries");

static_assert(std::extent_v<decltype(profiler_phase_id)> == int(ProfilerPhase::ProfilerPhaseCount),
              "Profiler phases and descriptions don't have matching length!");

#if defined(MTS_ENABLE_PROFILER)
/* Inlining the access to a thread_local variable produces *awful* machine code
   with Clang on OSX. The combination of weak and noinline is needed to prevent
   the compiler from inlining it (just noinline does not seem to be enough). It
   is marked as 'const' because separate function calls always produce the same
   pointer. */
extern MTS_EXPORT_CORE uint64_t *profiler_flags() __attribute__((noinline, weak, const));

struct ScopedPhase {
    ScopedPhase(ProfilerPhase phase) : m_target(profiler_flags()), m_flag(1ull << int(phase)) {
        if ((*m_target & m_flag) == 0)
            *m_target |= m_flag;
        else
            m_flag = 0;
    }

    ~ScopedPhase() { *m_target &= ~m_flag; }

    ScopedPhase(const ScopedPhase &) = delete;
    ScopedPhase &operator=(const ScopedPhase &) = delete;

private:
    uint64_t *m_target;
    uint64_t m_flag;
};

class MTS_EXPORT_CORE Profiler : public Object {
public:
    static void static_initialization();
    static void static_shutdown();
    static void print_report();
    MTS_DECLARE_CLASS()
private:
    Profiler() = delete;
};

#else

/* Profiler not supported on this platform */
struct ScopedPhase {
    ScopedPhase(ProfilerPhase) {}
};
class Profiler {
public:
    static void static_initialization() {}
    static void static_shutdown() {}
    static void print_report() {}
};

#endif

NAMESPACE_END(mitsuba)
