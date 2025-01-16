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

template <typename Float_, typename Spectrum_> class MTS_EXPORT_RENDER fermatNEE {
public:
    using Float    = Float_;
    using Spectrum = Spectrum_;
    MTS_IMPORT_TYPES(Scene, Sampler, Emitter) // Sampler, Scene, Emitter, Shape, Heightfield, MeshAttribute
    using EmitterPtr         = typename RenderAliases::EmitterPtr;
    using ShapePtr           = typename RenderAliases::ShapePtr;
    using ManifoldVertex     = ManifoldVertex<Float, Spectrum>;
    using EmitterInteraction = EmitterInteraction<Float, Spectrum>;
    using SpecularManifold   = SpecularManifold<Float, Spectrum>;

    fermatNEE() {}

    void init(const Scene *scene, Sampler *sampler, Float solutionIdentical_threshold, int maxBernouilliTrial, double alpha1, double beta,
              bool crop_caustic, const SMSConfig &config) {
        m_scene   = scene;
        m_sampler = sampler;

        this->solutionIdentical_threshold = solutionIdentical_threshold;
        this->maxBernouilliTrial          = maxBernouilliTrial;
        this->alpha1                      = alpha1;
        this->beta                        = beta;

        this->crop_caustic = crop_caustic;
        this->m_config         = config;
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

    // struct ManifoldVertex {
    //     using Float    = Float_;
    //     using Spectrum = Spectrum_;
    //     MTS_IMPORT_RENDER_BASIC_TYPES()
    //     using ShapePtr             = typename RenderAliases::ShapePtr;
    //     using SurfaceInteraction3f = typename RenderAliases::SurfaceInteraction3f;

    //     // Position and partials
    //     Point3f p;
    //     Vector3f dp_du, dp_dv;

    //     // Normal and partials
    //     Normal3f n, gn;
    //     Vector3f dn_du, dn_dv;

    //     // Tangents and partials
    //     Vector3f s, t;
    //     Vector3f ds_du, ds_dv;
    //     Vector3f dt_du, dt_dv;

    //     // Further information
    //     Float eta;
    //     Vector2f uv;
    //     ShapePtr shape;
    //     Mask fixed_direction;

    //     // Used in multi-bounce version
    //     Vector2f C;
    //     Matrix2f dC_dx_prev, dC_dx_cur, dC_dx_next;
    //     Matrix2f tmp, inv_lambda;
    //     Vector2f dx;

    //     ManifoldVertex(const Point3f &p = Point3f(0.f))
    //         : p(p), dp_du(0.f), dp_dv(0.f), n(0.f), gn(0.f), dn_du(0.f), dn_dv(0.f), s(0.f), t(0.f), ds_du(0.f), ds_dv(0.f), dt_du(0.f),
    //           dt_dv(0.f), eta(1.f), uv(0.f), shape(nullptr), fixed_direction(false) {}

    //     ManifoldVertex(const SurfaceInteraction3f &si, Float smoothing = 0.f)
    //         : p(si.p), dp_du(si.dp_du), dp_dv(si.dp_dv), gn(si.n), uv(si.uv), shape(si.shape), fixed_direction(false) {

    //         // Encode conductors with eta=1.0, and dielectrics with their relative IOR
    //         Complex<Spectrum> ior = si.bsdf()->ior(si);
    //         eta                   = select(all(eq(0.f, imag(ior))), hmean(real(ior)),
    //                                        1.f); // Assumption here is that real (dielectric) IOR is not spectrally varying.

    //         // Compute frame and its derivative
    //         Frame3f frame = si.bsdf()->frame(si, smoothing);
    //         n             = frame.n;
    //         s             = frame.s;
    //         t             = frame.t;

    //         auto [dframe_du, dframe_dv] = si.bsdf()->frame_derivative(si, smoothing);
    //         dn_du                       = dframe_du.n;
    //         dn_dv                       = dframe_dv.n;
    //         ds_du                       = dframe_du.s;
    //         ds_dv                       = dframe_dv.s;
    //         dt_du                       = dframe_du.t;
    //         dt_dv                       = dframe_dv.t;

    //         // In rare cases, e.g. 'twosided' materials, the geometric normal needs to be flipped
    //         masked(gn, dot(n, gn) < 0.f) *= -1.f;
    //     }

    //     void make_orthonormal() {
    //         // Turn into orthonormal parameterization at 'p'
    //         Float inv_norm = rcp(norm(dp_du));
    //         dp_du *= inv_norm;
    //         dn_du *= inv_norm;
    //         Float dp           = dot(dp_du, dp_dv);
    //         Vector3f dp_dv_tmp = dp_dv - dp * dp_du;
    //         Vector3f dn_dv_tmp = dn_dv - dp * dn_du;
    //         inv_norm           = rcp(norm(dp_dv_tmp));
    //         dp_dv              = dp_dv_tmp * inv_norm;
    //         dn_dv              = dn_dv_tmp * inv_norm;
    //     }

    //     std::string to_string() const {
    //         std::ostringstream oss;
    //         oss << "ManifoldVertex[" << std::endl
    //             << "  p = " << p << "," << std::endl
    //             << "  n = " << n << "," << std::endl
    //             << "  gn = " << gn << "," << std::endl
    //             << "  dp_du = " << dp_du << "," << std::endl
    //             << "  dp_dv = " << dp_dv << "," << std::endl
    //             << "  dn_du = " << dn_du << "," << std::endl
    //             << "  dn_dv = " << dn_dv << "," << std::endl
    //             << "  eta = " << eta << "," << std::endl
    //             << "  uv = " << uv << std::endl
    //             << "]";
    //         return oss.str();
    //     }
    // };

    /* +===============================================================================================+
    // |                                            SOLVER                                             |
    // +===============================================================================================+ */

    bool newton_solver_double(Array<Float, 4> *X0_, L_data_float *objfn_data, Float *f_out) const;
    bool newton_solver_double_(Array<Float, 4> *X0_, L_data_float *objfn_data, Float *f_out, Matrix<double, 4> *H_out=nullptr) const;
    bool newton_solver_reflection(Vector2f *X0_, L_data_float *objfn_data, Float *f_out, Float *out = nullptr, Matrix<double,2> *H_out=nullptr) const;
    Array<Float, 4> inline solve_4x4_enoki(const Matrix4f &m, const Array<Float, 4> b) const;
    Array<double, 4> inline solve_4x4_enoki_double(const Matrix<double, 4> &m, const Array<double, 4> b) const;

    bool newton_solver_SMS(const SurfaceInteraction3f &si, const EmitterInteraction &ei, L_data_float *data);
    bool newton_solver_SMS_test(const SurfaceInteraction3f &si, const EmitterInteraction &ei, L_data_float *data);

    Mask invert_tridiagonal_step(std::vector<ManifoldVertex> &v) {
        // Solve block tri-diagonal linear system with full RHS vector

        // From "The Natural-Constraint Representation of the Path Space for Efficient Light Transport Simulation"
        // by Kaplanyan et al. 2014 Supplemental material, Figure 2.

        auto invert = [](const Matrix2f &A, Matrix2f &Ainv) {
            Float determinant = det(A);
            if (abs(determinant) == 0) {
                return false;
            }
            Ainv = inverse(A);
            return true;
        };

        int n = int(v.size());
        if (n == 0)
            return true;

        v[0].tmp   = v[0].dC_dx_prev;
        Matrix2f m = v[0].dC_dx_cur;
        if (!invert(m, v[0].inv_lambda))
            return false;

        for (int i = 1; i < n; ++i) {
            v[i].tmp   = v[i].dC_dx_prev * v[i - 1].inv_lambda;
            Matrix2f m = v[i].dC_dx_cur - v[i].tmp * v[i - 1].dC_dx_next;
            if (!invert(m, v[i].inv_lambda))
                return false;
        }

        v[0].dx = v[0].C;
        for (int i = 1; i < n; ++i) {
            v[i].dx = v[i].C - v[i].tmp * v[i - 1].dx;
        }

        v[n - 1].dx = v[n - 1].inv_lambda * v[n - 1].dx;
        for (int i = n - 2; i >= 0; --i) {
            v[i].dx = v[i].inv_lambda * (v[i].dx - v[i].dC_dx_next * v[i + 1].dx);
        }

        return true;
    }

    ManifoldVertex manifoldVertex_from_uv(Point2f uv, ShapePtr heightfield, SurfaceInteraction3f &si_) ;

    /* +===============================================================================================+
    // |                              OBJECTIVE FUNCTION: OPTICAL LENGTH                               |
    // +===============================================================================================+ */

    inline Float static L_flat(Array<double, 4> x_in, Array<double, 4> *grad_out, Matrix<double, 4> *hess_out, L_data_float *objfn_data,
                               Float *f_out);
    inline double static L_double(Array<double, 4> x_in, Array<double, 4> *grad_out, Matrix<double, 4> *hess_out, L_data_float *objfn_data,
                                  Float *f_out);
    inline double static L_double_mesh(Array<double, 4> x_in, Array<double, 4> *grad_out, Matrix<double, 4> *hess_out, L_data_float *objfn_data,
                                  Float *f_out);
    inline double static L_reflection(Vector2f x_in, Vector2f *grad_out, Matrix<double, 2> *hess_out, L_data_float *objfn_data,
                                      Float *f_out);

    Mask compute_step_halfvector(const Point3f &x0, const EmitterInteraction &ei) {
        ScopedPhase scope_phase(ProfilerPhase::SMS_computeDirection);
        std::vector<ManifoldVertex> &v = m_current_path;

        size_t k = v.size();
        for (size_t i = 0; i < k; ++i) {
            v[i].C          = Vector2f(0.f);
            v[i].dC_dx_prev = Matrix2f(0.f);
            v[i].dC_dx_cur  = Matrix2f(0.f);
            v[i].dC_dx_next = Matrix2f(0.f);

            Point3f x_prev = (i == 0) ? x0 : v[i - 1].p;
            Point3f x_next = (i == k - 1) ? ei.p : v[i + 1].p;
            Point3f x_cur  = v[i].p;

            bool at_endpoint_with_fixed_direction = (i == (k - 1) && ei.is_directional());

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
            Float ili   = norm(wi);
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

                v[i].dC_dx_prev = Matrix2f(dot(v[i].s, dh_du), dot(v[i].s, dh_dv), dot(v[i].t, dh_du), dot(v[i].t, dh_dv));
            }

            // Derivative of specular constraint w.r.t. x_{i}
            if (at_endpoint_with_fixed_direction) {
                // When the 'wo' direction is fixed, the derivative here simplifies.
                dh_du = ili * (-v[i].dp_du + wi * dot(wi, v[i].dp_du));
                dh_dv = ili * (-v[i].dp_dv + wi * dot(wi, v[i].dp_dv));
            } else {
                // Standard case for fixed emitter position
                dh_du = -v[i].dp_du * (ili + ilo) + wi * (dot(wi, v[i].dp_du) * ili) + wo * (dot(wo, v[i].dp_du) * ilo);
                dh_dv = -v[i].dp_dv * (ili + ilo) + wi * (dot(wi, v[i].dp_dv) * ili) + wo * (dot(wo, v[i].dp_dv) * ilo);
            }
            dh_du -= h * dot(dh_du, h);
            dh_dv -= h * dot(dh_dv, h);
            if (eta != 1.f) {
                dh_du *= -1.f;
                dh_dv *= -1.f;
            }

            v[i].dC_dx_cur = Matrix2f(dot(v[i].ds_du, h) + dot(v[i].s, dh_du), dot(v[i].ds_dv, h) + dot(v[i].s, dh_dv),
                                      dot(v[i].dt_du, h) + dot(v[i].t, dh_du), dot(v[i].dt_dv, h) + dot(v[i].t, dh_dv));

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

                v[i].dC_dx_next = Matrix2f(dot(v[i].s, dh_du), dot(v[i].s, dh_dv), dot(v[i].t, dh_du), dot(v[i].t, dh_dv));
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

    Mask compute_step_anglediff(const Point3f &x0, const EmitterInteraction &ei) {
        std::vector<ManifoldVertex> &v = m_current_path;
        bool success = true;
        size_t k = v.size();
        for (size_t i = 0; i < k; ++i) {
            v[i].C = Vector2f(0.f);
            v[i].dC_dx_prev = Matrix2f(0.f);
            v[i].dC_dx_cur  = Matrix2f(0.f);
            v[i].dC_dx_next = Matrix2f(0.f);

            Point3f x_prev = (i == 0)   ? x0   : v[i-1].p;
            Point3f x_next = (i == k-1) ? ei.p : v[i+1].p;
            Point3f x_cur  = v[i].p;

            bool at_endpoint_with_fixed_direction = (i == (k-1) && ei.is_directional());

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
                dwo_du_cur = -ilo * (v[i].dp_du - wo*dot(wo, v[i].dp_du));
                dwo_dv_cur = -ilo * (v[i].dp_dv - wo*dot(wo, v[i].dp_dv));
            }

            Vector3f wi = x_prev - x_cur;
            Float ili = norm(wi);
            if (ili < 1e-3f) {
                return false;
            }
            ili = rcp(ili);
            wi *= ili;

            Vector3f dwi_du_cur = -ili * (v[i].dp_du - wi*dot(wi, v[i].dp_du)),
                     dwi_dv_cur = -ili * (v[i].dp_dv - wi*dot(wi, v[i].dp_dv));

            // Set up constraint function and its derivatives
            bool success_i = false;

            auto transform = [&](const Vector3f &w, const Vector3f &n, Float eta) {
                if (eta == 1.f) {
                    return SpecularManifold::reflect(w, n);
                } else {
                    return SpecularManifold::refract(w, n, eta);
                }
            };
            auto d_transform = [&](const Vector3f &w, const Vector3f &dw_du, const Vector3f &dw_dv,
                    const Vector3f &n, const Vector3f &dn_du, const Vector3f &dn_dv,
                    Float eta) {
                if (eta == 1.f) {
                    return SpecularManifold::d_reflect(w, dw_du, dw_dv, n, dn_du, dn_dv);
                } else {
                    return SpecularManifold::d_refract(w, dw_du, dw_dv, n, dn_du, dn_dv, eta);
                }
            };

            // Handle offset normal. These are no-ops in case n_offset=[0,0,1]
            Vector3f n_offset = m_offset_normals[i];
            Normal3f n = v[i].s * n_offset[0] +
                v[i].t * n_offset[1] +
                v[i].n * n_offset[2];
            Vector3f dn_du = v[i].ds_du * n_offset[0] +
                v[i].dt_du * n_offset[1] +
                v[i].dn_du * n_offset[2];
            Vector3f dn_dv = v[i].ds_dv * n_offset[0] +
                v[i].dt_dv * n_offset[1] +
                v[i].dn_dv * n_offset[2];

            auto [valid_i_refr_i, wio] = transform(wi, n, v[i].eta);
            if (valid_i_refr_i) {
                auto [to, po]   = SpecularManifold::sphcoords(wo);
                auto [tio, pio] = SpecularManifold::sphcoords(wio);

                Float dt = to - tio,
                      dp = po - pio;
                if (dp < -math::Pi<Float>) {
                    dp += 2.f*math::Pi<Float>;
                } else if (dp > math::Pi<Float>) {
                    dp -= 2.f*math::Pi<Float>;
                }
                v[i].C = Vector2f(dt, dp);

                Float dto_du, dpo_du, dto_dv, dpo_dv;
                Float dtio_du, dpio_du, dtio_dv, dpio_dv;

                // Derivative of specular constraint w.r.t. x_{i-1}
                if (i > 0) {
                    Vector3f dwi_du_prev = ili * (v[i-1].dp_du - wi*dot(wi, v[i-1].dp_du)),
                             dwi_dv_prev = ili * (v[i-1].dp_dv - wi*dot(wi, v[i-1].dp_dv));
                    // Vector3f dwo_du_prev = ilo * (v[i-1].dp_du - wo*dot(wo, v[i-1].dp_du)),  // = 0
                    //          dwo_dv_prev = ilo * (v[i-1].dp_dv - wo*dot(wo, v[i-1].dp_dv));  // = 0
                    auto [dwio_du_prev, dwio_dv_prev] = d_transform(wi, dwi_du_prev, dwi_dv_prev, n, Vector3f(0.f), Vector3f(0.f), v[i].eta);   // Possible optimization: specific implementation here that already knows some of these are 0.

                    // std::tie(dto_du, dpo_du, dto_dv, dpo_dv)     = SpecularManifold::d_sphcoords(wo, dwo_du_prev, dwo_dv_prev);  // = 0
                    std::tie(dtio_du, dpio_du, dtio_dv, dpio_dv) = SpecularManifold::d_sphcoords(wio, dwio_du_prev, dwio_dv_prev);

                    v[i].dC_dx_prev(0,0) = -dtio_du;
                    v[i].dC_dx_prev(1,0) = -dpio_du;
                    v[i].dC_dx_prev(0,1) = -dtio_dv;
                    v[i].dC_dx_prev(1,1) = -dpio_dv;
                }

                // Derivative of specular constraint w.r.t. x_{i}
                auto [dwio_du_cur, dwio_dv_cur] = d_transform(wi, dwi_du_cur, dwi_dv_cur, n, dn_du, dn_dv, v[i].eta);

                std::tie(dto_du, dpo_du, dto_dv, dpo_dv)     = SpecularManifold::d_sphcoords(wo, dwo_du_cur, dwo_dv_cur);
                std::tie(dtio_du, dpio_du, dtio_dv, dpio_dv) = SpecularManifold::d_sphcoords(wio, dwio_du_cur, dwio_dv_cur);

                v[i].dC_dx_cur(0,0) = dto_du - dtio_du;
                v[i].dC_dx_cur(1,0) = dpo_du - dpio_du;
                v[i].dC_dx_cur(0,1) = dto_dv - dtio_dv;
                v[i].dC_dx_cur(1,1) = dpo_dv - dpio_dv;

                // Derivative of specular constraint w.r.t. x_{i+1}
                if (i < k-1) {
                    // Vector3f dwi_du_next = ili * (v[i+1].dp_du - wi*dot(wi, v[i+1].dp_du)),  // = 0
                    //          dwi_dv_next = ili * (v[i+1].dp_dv - wi*dot(wi, v[i+1].dp_dv));  // = 0
                    Vector3f dwo_du_next = ilo * (v[i+1].dp_du - wo*dot(wo, v[i+1].dp_du)),
                             dwo_dv_next = ilo * (v[i+1].dp_dv - wo*dot(wo, v[i+1].dp_dv));
                    // auto [dwio_du_next, dwio_dv_next] = d_transform(wi, dwi_du_next, dwi_dv_next, n, Vector3f(0.f), Vector3f(0.f), v[i].eta); // = 0

                    std::tie(dto_du, dpo_du, dto_dv, dpo_dv) = SpecularManifold::d_sphcoords(wo, dwo_du_next, dwo_dv_next);
                    // std::tie(dtio_du, dpio_du, dtio_dv, dpio_dv) = SpecularManifold::d_sphcoords(wio, dwio_du_next, dwio_dv_next);   // = 0

                    v[i].dC_dx_next(0,0) = dto_du;
                    v[i].dC_dx_next(1,0) = dpo_du;
                    v[i].dC_dx_next(0,1) = dto_dv;
                    v[i].dC_dx_next(1,1) = dpo_dv;
                }

                success_i = true;
            }

            auto [valid_o_refr_o, woi] = transform(wo, n, v[i].eta);
            if (valid_o_refr_o && !success_i) {
                auto [ti, pi]   = SpecularManifold::sphcoords(wi);
                auto [toi, poi] = SpecularManifold::sphcoords(woi);

                Float dt = ti - toi,
                      dp = pi - poi;
                if (dp < -math::Pi<Float>) {
                    dp += 2.f*math::Pi<Float>;
                } else if (dp > math::Pi<Float>) {
                    dp -= 2.f*math::Pi<Float>;
                }
                v[i].C = Vector2f(dt, dp);

                Float dti_du, dpi_du, dti_dv, dpi_dv;
                Float dtoi_du, dpoi_du, dtoi_dv, dpoi_dv;

                // Derivative of specular constraint w.r.t. x_{i-1}
                if (i > 0) {
                    Vector3f dwi_du_prev = ili * (v[i-1].dp_du - wi*dot(wi, v[i-1].dp_du)),
                             dwi_dv_prev = ili * (v[i-1].dp_dv - wi*dot(wi, v[i-1].dp_dv));
                    // Vector3f dwo_du_prev = ilo * (v[i-1].dp_du - wo*dot(wo, v[i-1].dp_du)),  // = 0
                    // dwo_dv_prev = ilo * (v[i-1].dp_dv - wo*dot(wo, v[i-1].dp_dv));  // = 0
                    // auto [dwoi_du_prev, dwoi_dv_prev] = d_transform(wo, dwo_du_prev, dwo_dv_prev, n, Vector3f(0.f), Vector3f(0.f), v[i].eta);   // = 0

                    std::tie(dti_du, dpi_du, dti_dv, dpi_dv) = SpecularManifold::d_sphcoords(wi, dwi_du_prev, dwi_dv_prev);
                    // std::tie(dtoi_du, dpoi_du, dtoi_dv, dpoi_dv) = SpecularManifold::d_sphcoords(woi, dwoi_du_prev, dwoi_dv_prev);   // = 0

                    v[i].dC_dx_prev(0,0) = dti_du;
                    v[i].dC_dx_prev(1,0) = dpi_du;
                    v[i].dC_dx_prev(0,1) = dti_dv;
                    v[i].dC_dx_prev(1,1) = dpi_dv;
                }


                // Derivative of specular constraint w.r.t. x_{i}
                auto [dwoi_du_cur, dwoi_dv_cur] = d_transform(wo, dwo_du_cur, dwo_dv_cur, n, dn_du, dn_dv, v[i].eta);

                std::tie(dti_du, dpi_du, dti_dv, dpi_dv)     = SpecularManifold::d_sphcoords(wi, dwi_du_cur, dwi_dv_cur);
                std::tie(dtoi_du, dpoi_du, dtoi_dv, dpoi_dv) = SpecularManifold::d_sphcoords(woi, dwoi_du_cur, dwoi_dv_cur);

                v[i].dC_dx_cur(0,0) = dti_du - dtoi_du;
                v[i].dC_dx_cur(1,0) = dpi_du - dpoi_du;
                v[i].dC_dx_cur(0,1) = dti_dv - dtoi_dv;
                v[i].dC_dx_cur(1,1) = dpi_dv - dpoi_dv;

                // Derivative of specular constraint w.r.t. x_{i+1}
                if (i < k-1) {
                    // Vector3f dwi_du_next = ili * (v[i+1].dp_du - wi*dot(wi, v[i+1].dp_du)),  // = 0
                    // dwi_dv_next = ili * (v[i+1].dp_dv - wi*dot(wi, v[i+1].dp_dv));  // = 0
                    Vector3f dwo_du_next = ilo * (v[i+1].dp_du - wo*dot(wo, v[i+1].dp_du)),
                             dwo_dv_next = ilo * (v[i+1].dp_dv - wo*dot(wo, v[i+1].dp_dv));
                    auto [dwoi_du_next, dwoi_dv_next] = d_transform(wo, dwo_du_next, dwo_dv_next, n, Vector3f(0.f), Vector3f(0.f), v[i].eta);   // Possible optimization: specific implementation here that already knows some of these are 0.

                    // std::tie(dti_du, dpi_du, dti_dv, dpi_dv)  = SpecularManifold::d_sphcoords(wi, dwi_du_next, dwi_dv_next);  // = 0
                    std::tie(dtoi_du, dpoi_du, dtoi_dv, dpoi_dv) = SpecularManifold::d_sphcoords(woi, dwoi_du_next, dwoi_dv_next);

                    v[i].dC_dx_next(0,0) = -dtoi_du;
                    v[i].dC_dx_next(1,0) = -dpoi_du;
                    v[i].dC_dx_next(0,1) = -dtoi_dv;
                    v[i].dC_dx_next(1,1) = -dpoi_dv;
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


    /* +===============================================================================================+
    // |                                           GUIDE RAY                                           |
    // +===============================================================================================+ */

    // next event estimation through double refraction event
    bool fermat_connection(L_data_float *data, Vector3f *result) const;
    Float eval_invPDF(L_data_float *data, Vector3f &proposal) const;

    bool SMS_connection(SurfaceInteraction3f &si, const EmitterInteraction &ei, L_data_float *data, Vector3f *result,
                        ref<Sampler> sampler);
    Float SMS_eval_invPDF(SurfaceInteraction3f &si, const EmitterInteraction &ei, L_data_float *data, Vector3f &proposal,
                          ref<Sampler> sampler);

    // python debug
    std::pair<bool, Vector3f> fermat_connection_(SurfaceInteraction3f &si, const EmitterInteraction &vy, Point2f init = Point2f(-1, -1),
                                                 ShapePtr shape = nullptr, Float H_out = 0) const;


    // naive connection with a straight line
    Ray3f straighline_approx(ShapePtr H1, ShapePtr H2, Point3f O, Point3f S, Float n1, Float n2, const Float *invpdf,
                             bool bComputePDF) const;

    // next event estimation through a single reflection event
    bool fermat_connection_reflection(L_data_float *data, Vector3f *result, Float *out = nullptr, Point2f init = Point2f(-1, -1), Matrix<double, 2> *H_out = nullptr) const;
    Float eval_invPDF_reflection(L_data_float *data, Vector3f &proposal) const;

    /* +===============================================================================================+
    // |                                        SPECULAR OUTPUT                                        |
    // +===============================================================================================+ */
    Float invert_tridiagonal_geo(std::vector<ManifoldVertex> &v) const;

    Float generalizedGeometryFactor(const SurfaceInteraction3f &si0, const SurfaceInteraction3f &si1, const SurfaceInteraction3f &si2,
                                    const SurfaceInteraction3f &si3, const bool reflection = false) const;

    Spectrum compute_ray_contribution(SurfaceInteraction3f si_O, const BSDFContext &ctx, Ray3f ray, EmitterInteraction vy, Float n1,
                                      Float n2, ShapePtr H1, ShapePtr H2, Mask active) const;

    Spectrum compute_ray_contribution_reflection(SurfaceInteraction3f si_O, const BSDFContext &ctx, RayDifferential3f ray,
                                                 EmitterInteraction vy, Float n1, Float n2, ShapePtr H1, ShapePtr H2, Mask active) const;

    // +==============================================================================+
    // |                                   UTILITY                                    |
    // +==============================================================================+

    void OpenGL_draw_ray_diff(RayDifferential3f r, Vector3f color, SurfaceInteraction3f si) const;

    void OpenGL_draw_frame(Vector3f n, Vector3f t, Vector3f s, Vector3f color, SurfaceInteraction3f &si) const;

    static MTS_INLINE std::pair<Mask, Vector3f> refract(const Vector3f &w, const Normal3f &n_, Float eta_) {
        Normal3f n = n_;
        Float eta  = rcp(eta_);
        if (dot(w, n) < 0) {
            // Coming from the "inside"
            eta = rcp(eta);
            n *= -1.f;
        }
        Float dot_w_n   = dot(w, n);
        Float root_term = 1.f - eta * eta * (1.f - dot_w_n * dot_w_n);
        if (root_term < 0.f) {
            return std::make_pair(false, Vector3f(0.f));
        }
        Vector3f wt = -eta * (w - dot_w_n * n) - n * sqrt(root_term);
        return std::make_pair(true, wt);
    }

    bool get_partition(L_data_float *data) const;

    MTS_INLINE
    Array<Float, 4> sample_partition(Vector2f partition, L_data_float *data) const {
        return Array<Float, 4>(partition.x() * data->partition_size + m_sampler->next_1d() * data->partition_size,
                               partition.y() * data->partition_size + m_sampler->next_1d() * data->partition_size,
                               partition.x() * data->partition_size + m_sampler->next_1d() * data->partition_size,
                               partition.y() * data->partition_size + m_sampler->next_1d() * data->partition_size);
    }

    MTS_INLINE
    bool is_inside_partition(Array<Float, 4> X, Vector2f partition, L_data_float *data) const {
        // std::cout<<partition<<std::endl;
        if (X.x() > 1.0f || X.x() < 0.0f||
            X.y() > 1.0f || X.y() < 0.0f||
            X.z() > 1.0f || X.z() < 0.0f||
            X.w() > 1.0f || X.w() < 0.0f) {
            return false;
        } else {
            return true;
        }
    }

    Mask reproject_raytrace(const SurfaceInteraction3f &si_);
    Mask reproject(const SurfaceInteraction3f &si_, L_data_float *data);

    // +==============================================================================+
    // |                                   Sampling                                   |
    // +==============================================================================+
    EmitterInteraction sample_emitter(SurfaceInteraction3f &si, Mask active) const;

    bool sample_heightfield_pair(ShapePtr *H1, ShapePtr *H2, Float * invPdf , Mask active = true) const;
    bool sample_heightfield(ShapePtr *H1, Float *invPdf, Mask active = true) const;
    // +==============================================================================+
    // |                                   MEMBERS                                    |
    // +==============================================================================+

    Float solutionIdentical_threshold;
    int maxBernouilliTrial;
    double alpha1;
    double beta;
    bool crop_caustic;

protected:
    const Scene *m_scene = nullptr;
    Sampler *m_sampler   = nullptr;
    SMSConfig m_config;

    std::vector<ManifoldVertex> m_seed_path, m_current_path, m_proposed_path;
    std::vector<Point3f> m_proposed_positions;
    std::vector<Vector3f> m_offset_normals;
};
MTS_EXTERN_CLASS_RENDER(fermatNEE)
NAMESPACE_END(mitsuba)
