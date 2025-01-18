#pragma once

#include "fwd.h"
#include <mitsuba/core/fwd.h>
#include <mitsuba/core/warp.h>
#include <mitsuba/render/OpenGL_viewer_client.h>
#include <mitsuba/render/bsdf.h>
#include <mitsuba/render/emitter.h>
#include <mitsuba/render/fwd.h>
#include <mitsuba/render/manifold.h>
#include <mitsuba/render/shape.h>

NAMESPACE_BEGIN(mitsuba)

template <typename Float_, typename Spectrum_>
class MTS_EXPORT_RENDER FermatNEE {
public:
    using Float = Float_;
    using Spectrum = Spectrum_;
    MTS_IMPORT_TYPES(Scene, Sampler, Emitter) // Sampler, Scene, Emitter, Shape,
                                              // Heightfield, MeshAttribute
    using EmitterPtr = typename RenderAliases::EmitterPtr;
    using ShapePtr = typename RenderAliases::ShapePtr;
    using ManifoldVertex = ManifoldVertex<Float, Spectrum>;
    using EmitterInteraction = EmitterInteraction<Float, Spectrum>;
    using SpecularManifold = SpecularManifold<Float, Spectrum>;

    FermatNEE() {}

    void init(const Scene *scene, Sampler *sampler,
              Float solutionIdentical_threshold, int maxBernouilliTrial,
              double alpha1, double beta, bool crop_caustic,
              const SMSConfig &config, bool use_SMS) {
        m_scene = scene;
        m_sampler = sampler;

        this->solutionIdentical_threshold = solutionIdentical_threshold;
        this->maxBernouilliTrial = maxBernouilliTrial;
        this->alpha1 = alpha1;
        this->beta = beta;

        this->crop_caustic = crop_caustic;
        this->m_config = config;
        m_use_SMS = use_SMS;
    }

    struct L_data_float {
        ShapePtr pH1, pH2;
        Float H1, H2, L1, L2, W1, W2;
        Point3f O, S;
        Float n1, n2;
        Float e;
        Float partition_size; // length of partition side
        Float inv_partition_size;
        Float inv_partition_area;
        Vector2f partition; // partition index
    };

    /************
     *  Solver  *
     ************/

    bool newton_solver_double(Array<Float, 4> *X0_, L_data_float *objfn_data,
                              Float *f_out) const;

    bool newton_solver_reflection(Vector2f *X0_, L_data_float *objfn_data,
                                  Float *f_out, Float *out = nullptr,
                                  Matrix<double, 2> *H_out = nullptr) const;

    Array<Float, 4> inline solve_4x4_enoki(const Matrix4f &m,
                                           const Array<Float, 4> b) const;

    Array<double, 4> inline solve_4x4_enoki_double(
        const Matrix<double, 4> &m, const Array<double, 4> b) const;

    bool newton_solver_SMS(const SurfaceInteraction3f &si,
                           const EmitterInteraction &ei, L_data_float *data);

    Mask invert_tridiagonal_step(std::vector<ManifoldVertex> &v);

    ManifoldVertex manifoldVertex_from_uv(Point2f uv, ShapePtr heightfield,
                                          SurfaceInteraction3f &si_);

    /* +===============================================================================================+
    // |                              OBJECTIVE FUNCTION: OPTICAL LENGTH |
    //
    +===============================================================================================+
  */

    inline Float static L_flat(Array<double, 4> x_in,
                               Array<double, 4> *grad_out,
                               Matrix<double, 4> *hess_out,
                               L_data_float *objfn_data, Float *f_out);
    inline double static L_double(Array<double, 4> x_in,
                                  Array<double, 4> *grad_out,
                                  Matrix<double, 4> *hess_out,
                                  L_data_float *objfn_data, Float *f_out);
    inline double static L_reflection(Vector2f x_in, Vector2f *grad_out,
                                      Matrix<double, 2> *hess_out,
                                      L_data_float *objfn_data, Float *f_out);

    // copied from manifold_ms to have everything selfcontained
    Mask compute_step_halfvector(const Point3f &x0,
                                 const EmitterInteraction &ei) {
        ScopedPhase scope_phase(ProfilerPhase::SMS_computeDirection);
        std::vector<ManifoldVertex> &v = m_current_path;

        size_t k = v.size();
        for (size_t i = 0; i < k; ++i) {
            v[i].C = Vector2f(0.f);
            v[i].dC_dx_prev = Matrix2f(0.f);
            v[i].dC_dx_cur = Matrix2f(0.f);
            v[i].dC_dx_next = Matrix2f(0.f);

            Point3f x_prev = (i == 0) ? x0 : v[i - 1].p;
            Point3f x_next = (i == k - 1) ? ei.p : v[i + 1].p;
            Point3f x_cur = v[i].p;

            bool at_endpoint_with_fixed_direction =
                (i == (k - 1) && ei.is_directional());

            // Setup wi / wo
            Vector3f wo;
            if (at_endpoint_with_fixed_direction) {
                // Case of fixed 'wo' direction
                wo = ei.d;
            } else {
                // Standard case for fixed emitter position
                wo = x_next - x_cur;
            }
            Float ilo = norm(wo);
            if (ilo < 1e-3f) {
                return false;
            }
            ilo = rcp(ilo);
            wo *= ilo;

            Vector3f wi = x_prev - x_cur;
            Float ili = norm(wi);
            if (ili < 1e-3f) {
                return false;
            }
            ili = rcp(ili);
            wi *= ili;

            // Setup generalized half-vector
            Float eta = v[i].eta;
            if (dot(wi, v[i].gn) < 0.f) {
                eta = rcp(eta);
            }
            Vector3f h = wi + eta * wo;
            if (eta != 1.f)
                h *= -1.f;
            Float ilh = rcp(norm(h));
            h *= ilh;

            ilo *= eta * ilh;
            ili *= ilh;

            Vector3f dh_du, dh_dv;

            // Derivative of specular constraint w.r.t. x_{i-1}
            if (i > 0) {
                dh_du = ili * (v[i - 1].dp_du - wi * dot(wi, v[i - 1].dp_du));
                dh_dv = ili * (v[i - 1].dp_dv - wi * dot(wi, v[i - 1].dp_dv));

                dh_du -= h * dot(dh_du, h);
                dh_dv -= h * dot(dh_dv, h);
                if (eta != 1.f) {
                    dh_du *= -1.f;
                    dh_dv *= -1.f;
                }

                v[i].dC_dx_prev =
                    Matrix2f(dot(v[i].s, dh_du), dot(v[i].s, dh_dv),
                             dot(v[i].t, dh_du), dot(v[i].t, dh_dv));
            }

            // Derivative of specular constraint w.r.t. x_{i}
            if (at_endpoint_with_fixed_direction) {
                // When the 'wo' direction is fixed, the derivative here
                // simplifies.
                dh_du = ili * (-v[i].dp_du + wi * dot(wi, v[i].dp_du));
                dh_dv = ili * (-v[i].dp_dv + wi * dot(wi, v[i].dp_dv));
            } else {
                // Standard case for fixed emitter position
                dh_du = -v[i].dp_du * (ili + ilo) +
                        wi * (dot(wi, v[i].dp_du) * ili) +
                        wo * (dot(wo, v[i].dp_du) * ilo);
                dh_dv = -v[i].dp_dv * (ili + ilo) +
                        wi * (dot(wi, v[i].dp_dv) * ili) +
                        wo * (dot(wo, v[i].dp_dv) * ilo);
            }
            dh_du -= h * dot(dh_du, h);
            dh_dv -= h * dot(dh_dv, h);
            if (eta != 1.f) {
                dh_du *= -1.f;
                dh_dv *= -1.f;
            }

            v[i].dC_dx_cur = Matrix2f(dot(v[i].ds_du, h) + dot(v[i].s, dh_du),
                                      dot(v[i].ds_dv, h) + dot(v[i].s, dh_dv),
                                      dot(v[i].dt_du, h) + dot(v[i].t, dh_du),
                                      dot(v[i].dt_dv, h) + dot(v[i].t, dh_dv));

            // Derivative of specular constraint w.r.t. x_{i+1}
            if (i < k - 1) {
                dh_du = ilo * (v[i + 1].dp_du - wo * dot(wo, v[i + 1].dp_du));
                dh_dv = ilo * (v[i + 1].dp_dv - wo * dot(wo, v[i + 1].dp_dv));

                dh_du -= h * dot(dh_du, h);
                dh_dv -= h * dot(dh_dv, h);
                if (eta != 1.f) {
                    dh_du *= -1.f;
                    dh_dv *= -1.f;
                }

                v[i].dC_dx_next =
                    Matrix2f(dot(v[i].s, dh_du), dot(v[i].s, dh_dv),
                             dot(v[i].t, dh_du), dot(v[i].t, dh_dv));
            }

            // Evaluate specular constraint
            Vector2f H(dot(v[i].s, h), dot(v[i].t, h));
            Vector3f n_offset = m_offset_normals[i];

            Vector2f N(n_offset[0], n_offset[1]);
            v[i].C = H - N;
        }

        if (!invert_tridiagonal_step(v)) {
            return false;
        }

        return true;
    }

    // copied from manifold_ms to have everything selfcontained
    Mask compute_step_anglediff(const Point3f &x0,
                                const EmitterInteraction &ei) {
        std::vector<ManifoldVertex> &v = m_current_path;
        bool success = true;
        size_t k = v.size();
        for (size_t i = 0; i < k; ++i) {
            v[i].C = Vector2f(0.f);
            v[i].dC_dx_prev = Matrix2f(0.f);
            v[i].dC_dx_cur = Matrix2f(0.f);
            v[i].dC_dx_next = Matrix2f(0.f);

            Point3f x_prev = (i == 0) ? x0 : v[i - 1].p;
            Point3f x_next = (i == k - 1) ? ei.p : v[i + 1].p;
            Point3f x_cur = v[i].p;

            bool at_endpoint_with_fixed_direction =
                (i == (k - 1) && ei.is_directional());

            // Setup wi / wo
            Vector3f wo;
            if (at_endpoint_with_fixed_direction) {
                // Case of fixed 'wo' direction
                wo = ei.d;
            } else {
                // Standard case for fixed emitter position
                wo = x_next - x_cur;
            }
            Float ilo = norm(wo);
            if (ilo < 1e-3f) {
                return false;
            }
            ilo = rcp(ilo);
            wo *= ilo;

            Vector3f dwo_du_cur, dwo_dv_cur;
            if (at_endpoint_with_fixed_direction) {
                // Fixed 'wo' direction means its derivative must be zero
                dwo_du_cur = Vector3f(0.f);
                dwo_dv_cur = Vector3f(0.f);
            } else {
                // Standard case for fixed emitter position
                dwo_du_cur = -ilo * (v[i].dp_du - wo * dot(wo, v[i].dp_du));
                dwo_dv_cur = -ilo * (v[i].dp_dv - wo * dot(wo, v[i].dp_dv));
            }

            Vector3f wi = x_prev - x_cur;
            Float ili = norm(wi);
            if (ili < 1e-3f) {
                return false;
            }
            ili = rcp(ili);
            wi *= ili;

            Vector3f dwi_du_cur =
                         -ili * (v[i].dp_du - wi * dot(wi, v[i].dp_du)),
                     dwi_dv_cur =
                         -ili * (v[i].dp_dv - wi * dot(wi, v[i].dp_dv));

            // Set up constraint function and its derivatives
            bool success_i = false;

            auto transform = [&](const Vector3f &w, const Vector3f &n,
                                 Float eta) {
                if (eta == 1.f) {
                    return SpecularManifold::reflect(w, n);
                } else {
                    return SpecularManifold::refract(w, n, eta);
                }
            };
            auto d_transform = [&](const Vector3f &w, const Vector3f &dw_du,
                                   const Vector3f &dw_dv, const Vector3f &n,
                                   const Vector3f &dn_du, const Vector3f &dn_dv,
                                   Float eta) {
                if (eta == 1.f) {
                    return SpecularManifold::d_reflect(w, dw_du, dw_dv, n,
                                                       dn_du, dn_dv);
                } else {
                    return SpecularManifold::d_refract(w, dw_du, dw_dv, n,
                                                       dn_du, dn_dv, eta);
                }
            };

            // Handle offset normal. These are no-ops in case n_offset=[0,0,1]
            Vector3f n_offset = m_offset_normals[i];
            Normal3f n = v[i].s * n_offset[0] + v[i].t * n_offset[1] +
                         v[i].n * n_offset[2];
            Vector3f dn_du = v[i].ds_du * n_offset[0] +
                             v[i].dt_du * n_offset[1] +
                             v[i].dn_du * n_offset[2];
            Vector3f dn_dv = v[i].ds_dv * n_offset[0] +
                             v[i].dt_dv * n_offset[1] +
                             v[i].dn_dv * n_offset[2];

            auto [valid_i_refr_i, wio] = transform(wi, n, v[i].eta);
            if (valid_i_refr_i) {
                auto [to, po] = SpecularManifold::sphcoords(wo);
                auto [tio, pio] = SpecularManifold::sphcoords(wio);

                Float dt = to - tio, dp = po - pio;
                if (dp < -math::Pi<Float>) {
                    dp += 2.f * math::Pi<Float>;
                } else if (dp > math::Pi<Float>) {
                    dp -= 2.f * math::Pi<Float>;
                }
                v[i].C = Vector2f(dt, dp);

                Float dto_du, dpo_du, dto_dv, dpo_dv;
                Float dtio_du, dpio_du, dtio_dv, dpio_dv;

                // Derivative of specular constraint w.r.t. x_{i-1}
                if (i > 0) {
                    Vector3f dwi_du_prev = ili * (v[i - 1].dp_du -
                                                  wi * dot(wi, v[i - 1].dp_du)),
                             dwi_dv_prev = ili * (v[i - 1].dp_dv -
                                                  wi * dot(wi, v[i - 1].dp_dv));
                    // Vector3f dwo_du_prev = ilo * (v[i-1].dp_du - wo*dot(wo,
                    // v[i-1].dp_du)),  // = 0
                    //          dwo_dv_prev = ilo * (v[i-1].dp_dv - wo*dot(wo,
                    //          v[i-1].dp_dv));  // = 0
                    auto [dwio_du_prev, dwio_dv_prev] = d_transform(
                        wi, dwi_du_prev, dwi_dv_prev, n, Vector3f(0.f),
                        Vector3f(0.f),
                        v[i].eta); // Possible optimization: specific
                                   // implementation here that already knows
                                   // some of these are 0.

                    // std::tie(dto_du, dpo_du, dto_dv, dpo_dv)     =
                    // SpecularManifold::d_sphcoords(wo, dwo_du_prev,
                    // dwo_dv_prev);  // = 0
                    std::tie(dtio_du, dpio_du, dtio_dv, dpio_dv) =
                        SpecularManifold::d_sphcoords(wio, dwio_du_prev,
                                                      dwio_dv_prev);

                    v[i].dC_dx_prev(0, 0) = -dtio_du;
                    v[i].dC_dx_prev(1, 0) = -dpio_du;
                    v[i].dC_dx_prev(0, 1) = -dtio_dv;
                    v[i].dC_dx_prev(1, 1) = -dpio_dv;
                }

                // Derivative of specular constraint w.r.t. x_{i}
                auto [dwio_du_cur, dwio_dv_cur] = d_transform(
                    wi, dwi_du_cur, dwi_dv_cur, n, dn_du, dn_dv, v[i].eta);

                std::tie(dto_du, dpo_du, dto_dv, dpo_dv) =
                    SpecularManifold::d_sphcoords(wo, dwo_du_cur, dwo_dv_cur);
                std::tie(dtio_du, dpio_du, dtio_dv, dpio_dv) =
                    SpecularManifold::d_sphcoords(wio, dwio_du_cur,
                                                  dwio_dv_cur);

                v[i].dC_dx_cur(0, 0) = dto_du - dtio_du;
                v[i].dC_dx_cur(1, 0) = dpo_du - dpio_du;
                v[i].dC_dx_cur(0, 1) = dto_dv - dtio_dv;
                v[i].dC_dx_cur(1, 1) = dpo_dv - dpio_dv;

                // Derivative of specular constraint w.r.t. x_{i+1}
                if (i < k - 1) {
                    // Vector3f dwi_du_next = ili * (v[i+1].dp_du - wi*dot(wi,
                    // v[i+1].dp_du)),  // = 0
                    //          dwi_dv_next = ili * (v[i+1].dp_dv - wi*dot(wi,
                    //          v[i+1].dp_dv));  // = 0
                    Vector3f dwo_du_next = ilo * (v[i + 1].dp_du -
                                                  wo * dot(wo, v[i + 1].dp_du)),
                             dwo_dv_next = ilo * (v[i + 1].dp_dv -
                                                  wo * dot(wo, v[i + 1].dp_dv));
                    // auto [dwio_du_next, dwio_dv_next] = d_transform(wi,
                    // dwi_du_next, dwi_dv_next, n, Vector3f(0.f),
                    // Vector3f(0.f), v[i].eta); // = 0

                    std::tie(dto_du, dpo_du, dto_dv, dpo_dv) =
                        SpecularManifold::d_sphcoords(wo, dwo_du_next,
                                                      dwo_dv_next);
                    // std::tie(dtio_du, dpio_du, dtio_dv, dpio_dv) =
                    // SpecularManifold::d_sphcoords(wio, dwio_du_next,
                    // dwio_dv_next);   // = 0

                    v[i].dC_dx_next(0, 0) = dto_du;
                    v[i].dC_dx_next(1, 0) = dpo_du;
                    v[i].dC_dx_next(0, 1) = dto_dv;
                    v[i].dC_dx_next(1, 1) = dpo_dv;
                }

                success_i = true;
            }

            auto [valid_o_refr_o, woi] = transform(wo, n, v[i].eta);
            if (valid_o_refr_o && !success_i) {
                auto [ti, pi] = SpecularManifold::sphcoords(wi);
                auto [toi, poi] = SpecularManifold::sphcoords(woi);

                Float dt = ti - toi, dp = pi - poi;
                if (dp < -math::Pi<Float>) {
                    dp += 2.f * math::Pi<Float>;
                } else if (dp > math::Pi<Float>) {
                    dp -= 2.f * math::Pi<Float>;
                }
                v[i].C = Vector2f(dt, dp);

                Float dti_du, dpi_du, dti_dv, dpi_dv;
                Float dtoi_du, dpoi_du, dtoi_dv, dpoi_dv;

                // Derivative of specular constraint w.r.t. x_{i-1}
                if (i > 0) {
                    Vector3f dwi_du_prev = ili * (v[i - 1].dp_du -
                                                  wi * dot(wi, v[i - 1].dp_du)),
                             dwi_dv_prev = ili * (v[i - 1].dp_dv -
                                                  wi * dot(wi, v[i - 1].dp_dv));
                    // Vector3f dwo_du_prev = ilo * (v[i-1].dp_du - wo*dot(wo,
                    // v[i-1].dp_du)),  // = 0 dwo_dv_prev = ilo * (v[i-1].dp_dv
                    // - wo*dot(wo, v[i-1].dp_dv));  // = 0 auto [dwoi_du_prev,
                    // dwoi_dv_prev] = d_transform(wo, dwo_du_prev, dwo_dv_prev,
                    // n, Vector3f(0.f), Vector3f(0.f), v[i].eta);   // = 0

                    std::tie(dti_du, dpi_du, dti_dv, dpi_dv) =
                        SpecularManifold::d_sphcoords(wi, dwi_du_prev,
                                                      dwi_dv_prev);
                    // std::tie(dtoi_du, dpoi_du, dtoi_dv, dpoi_dv) =
                    // SpecularManifold::d_sphcoords(woi, dwoi_du_prev,
                    // dwoi_dv_prev);   // = 0

                    v[i].dC_dx_prev(0, 0) = dti_du;
                    v[i].dC_dx_prev(1, 0) = dpi_du;
                    v[i].dC_dx_prev(0, 1) = dti_dv;
                    v[i].dC_dx_prev(1, 1) = dpi_dv;
                }

                // Derivative of specular constraint w.r.t. x_{i}
                auto [dwoi_du_cur, dwoi_dv_cur] = d_transform(
                    wo, dwo_du_cur, dwo_dv_cur, n, dn_du, dn_dv, v[i].eta);

                std::tie(dti_du, dpi_du, dti_dv, dpi_dv) =
                    SpecularManifold::d_sphcoords(wi, dwi_du_cur, dwi_dv_cur);
                std::tie(dtoi_du, dpoi_du, dtoi_dv, dpoi_dv) =
                    SpecularManifold::d_sphcoords(woi, dwoi_du_cur,
                                                  dwoi_dv_cur);

                v[i].dC_dx_cur(0, 0) = dti_du - dtoi_du;
                v[i].dC_dx_cur(1, 0) = dpi_du - dpoi_du;
                v[i].dC_dx_cur(0, 1) = dti_dv - dtoi_dv;
                v[i].dC_dx_cur(1, 1) = dpi_dv - dpoi_dv;

                // Derivative of specular constraint w.r.t. x_{i+1}
                if (i < k - 1) {
                    // Vector3f dwi_du_next = ili * (v[i+1].dp_du - wi*dot(wi,
                    // v[i+1].dp_du)),  // = 0 dwi_dv_next = ili * (v[i+1].dp_dv
                    // - wi*dot(wi, v[i+1].dp_dv));  // = 0
                    Vector3f dwo_du_next = ilo * (v[i + 1].dp_du -
                                                  wo * dot(wo, v[i + 1].dp_du)),
                             dwo_dv_next = ilo * (v[i + 1].dp_dv -
                                                  wo * dot(wo, v[i + 1].dp_dv));
                    auto [dwoi_du_next, dwoi_dv_next] = d_transform(
                        wo, dwo_du_next, dwo_dv_next, n, Vector3f(0.f),
                        Vector3f(0.f),
                        v[i].eta); // Possible optimization: specific
                                   // implementation here that already knows
                                   // some of these are 0.

                    // std::tie(dti_du, dpi_du, dti_dv, dpi_dv)  =
                    // SpecularManifold::d_sphcoords(wi, dwi_du_next,
                    // dwi_dv_next);  // = 0
                    std::tie(dtoi_du, dpoi_du, dtoi_dv, dpoi_dv) =
                        SpecularManifold::d_sphcoords(woi, dwoi_du_next,
                                                      dwoi_dv_next);

                    v[i].dC_dx_next(0, 0) = -dtoi_du;
                    v[i].dC_dx_next(1, 0) = -dpoi_du;
                    v[i].dC_dx_next(0, 1) = -dtoi_dv;
                    v[i].dC_dx_next(1, 1) = -dpoi_dv;
                }

                success_i = true;
            }

            success &= success_i;
        }

        if (!success || !invert_tridiagonal_step(v)) {
            return false;
        }
        return true;
    }

    // next event estimation through double refraction event
    std::tuple<bool, Vector3f> specular_connection(SurfaceInteraction3f &si,
                                                   const EmitterInteraction &ei,
                                                   L_data_float *data);

    std::tuple<bool, Vector3f> fermat_connection(SurfaceInteraction3f &si,
                                                 const EmitterInteraction &ei,
                                                 L_data_float *data) const;

    Float eval_invPDF(L_data_float *data, Vector3f &proposal) const;

    std::tuple<bool, Vector3f> SMS_connection(SurfaceInteraction3f &si,
                                              const EmitterInteraction &ei,
                                              L_data_float *data);

    Float SMS_eval_invPDF(SurfaceInteraction3f &si,
                          const EmitterInteraction &ei, L_data_float *data,
                          Vector3f &proposal);

    // python debug
    std::pair<bool, Vector3f> fermat_connection_(SurfaceInteraction3f &si,
                                                 const EmitterInteraction &vy,
                                                 Point2f init = Point2f(-1, -1),
                                                 ShapePtr shape = nullptr,
                                                 Float H_out = 0) const;

    // naive connection with a straight line
    Ray3f straighline_approx(ShapePtr H1, ShapePtr H2, Point3f O, Point3f S,
                             Float n1, Float n2, const Float *invpdf,
                             bool bComputePDF) const;

    // next event estimation through a single reflection event
    bool fermat_connection_reflection(L_data_float *data, Vector3f *result,
                                      Float *out = nullptr,
                                      Point2f init = Point2f(-1, -1),
                                      Matrix<double, 2> *H_out = nullptr) const;
    Float eval_invPDF_reflection(L_data_float *data, Vector3f &proposal) const;

    Float invert_tridiagonal_geo(std::vector<ManifoldVertex> &v) const;

    // adapted from manifold_ms
    Float generalizedGeometryFactor(const SurfaceInteraction3f &si0,
                                    const SurfaceInteraction3f &si1,
                                    const SurfaceInteraction3f &si2,
                                    const SurfaceInteraction3f &si3,
                                    const bool reflection = false) const;

    // adapted from manifold_ms
    Spectrum compute_ray_contribution(SurfaceInteraction3f si_O,
                                      const BSDFContext &ctx, Ray3f ray,
                                      EmitterInteraction vy, Float n1, Float n2,
                                      ShapePtr H1, ShapePtr H2,
                                      Mask active) const;

    // adapted from manifold_ss
    Spectrum compute_ray_contribution_reflection(
        SurfaceInteraction3f si_O, const BSDFContext &ctx,
        RayDifferential3f ray, EmitterInteraction vy, Float n1, Float n2,
        ShapePtr H1, ShapePtr H2, Mask active) const;

    // +==============================================================================+
    // |                                   UTILITY |
    // +==============================================================================+

    static MTS_INLINE std::pair<Mask, Vector3f>
    refract(const Vector3f &w, const Normal3f &n_, Float eta_) {
        Normal3f n = n_;
        Float eta = rcp(eta_);
        if (dot(w, n) < 0) {
            // Coming from the "inside"
            eta = rcp(eta);
            n *= -1.f;
        }
        Float dot_w_n = dot(w, n);
        Float root_term = 1.f - eta * eta * (1.f - dot_w_n * dot_w_n);
        if (root_term < 0.f) {
            return std::make_pair(false, Vector3f(0.f));
        }
        Vector3f wt = -eta * (w - dot_w_n * n) - n * sqrt(root_term);
        return std::make_pair(true, wt);
    }

    Mask reproject_raytrace(const SurfaceInteraction3f &si_);
    Mask reproject(const SurfaceInteraction3f &si_, L_data_float *data);

    // +==============================================================================+
    // |                                   Sampling |
    // +==============================================================================+
    EmitterInteraction sample_emitter(SurfaceInteraction3f &si,
                                      Mask active) const;

    bool sample_heightfield_pair(ShapePtr *H1, ShapePtr *H2, Float *invPdf,
                                 Mask active = true) const;
    bool sample_heightfield(ShapePtr *H1, Float *invPdf,
                            Mask active = true) const;

    // +==============================================================================+
    // |                                   MEMBERS |
    // +==============================================================================+

    Float solutionIdentical_threshold;
    int maxBernouilliTrial;
    double alpha1;
    double beta;
    bool crop_caustic;
    bool m_use_SMS;

protected:
    const Scene *m_scene = nullptr;
    Sampler *m_sampler = nullptr;
    SMSConfig m_config;

    std::vector<ManifoldVertex> m_seed_path, m_current_path, m_proposed_path;
    std::vector<Point3f> m_proposed_positions;
    std::vector<Vector3f> m_offset_normals;
};
MTS_EXTERN_CLASS_RENDER(FermatNEE)
NAMESPACE_END(mitsuba)
