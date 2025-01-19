#include <mitsuba/render/fermatNEE.h>

#include "stdlib.h"
#include <cstddef>
#include <enoki/stl.h>
#include <fstream>
#include <iterator>
#include <mitsuba/core/plugin.h>
#include <mitsuba/core/properties.h>
#include <mitsuba/render/bsdf.h>
#include <mitsuba/render/fermatNEE.h>
#include <mitsuba/render/integrator.h>
#include <mitsuba/render/kdtree.h>
#include <mitsuba/render/medium.h>
#include <mitsuba/render/scene.h>
#include <random>
#include <type_traits>

#if defined(MTS_ENABLE_EMBREE)
#include "scene_embree.inl"
#else
#include "scene_native.inl"
#endif

NAMESPACE_BEGIN(mitsuba)

template <typename Float, typename Spectrum>
bool FermatNEE<Float, Spectrum>::newton_solver_double(Array<Float, 4> *X0_,
                                                      L_data_float *objfn_data,
                                                      Float *f_out) const {
    ScopedPhase scope_phase(ProfilerPhase::FermaNEE_NewtonSolver);

    double t = 1.0;

    Array<double, 4> G(0);
    Matrix<double, 4> H(0);
    Array<double, 4> X1, X0 = Array<double, 4>(*(X0_));

    for (int i = 0; i < 20; i++) {

        double f = L_double(X0, &G, &H, objfn_data, f_out);

        Array<double, 4> v = -solve_4x4_enoki_double(H, G);

        // line search Armijo
        X1 = X0 + t * v;

        // backtracking Armijo strategy
        t = 1.0;
        while (L_double(X1, nullptr, nullptr, objfn_data, nullptr) >
               f + alpha1 * t * dot(G, v)) {
            t = beta * t;
            X1 = X0 + t * v;

            // TODO param
            if (t < 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5)
                break;
        }
        X0 = X1;

        // gradient norm small enough, return
        if (norm(G) < 0.00001) {
            *X0_ = X0;
            L_double(X0, nullptr, nullptr, objfn_data, f_out);
            return true;
        }
    }
    *X0_ = X0;
    return false;
}

template <typename Float, typename Spectrum>
bool FermatNEE<Float, Spectrum>::newton_solver_reflection(
    Vector2f *X0_, L_data_float *objfn_data, Float *f_out, Float *out,
    Matrix<double, 2> *H_out) const {
    ScopedPhase scope_phase(ProfilerPhase::FermaNEE_NewtonSolver);

    double t = 1.0;

    Vector2f G(0);
    Matrix<double, 2> H(0);
    Vector2f X1, X0 = Vector2f(*(X0_));

    for (int i = 0; i < 200; i++) {

        double f = L_reflection(X0, &G, &H, objfn_data, f_out);

        // Array<double, 4> v = -solve_4x4_enoki_double(H, G);
        Vector2f v = -inverse(H) * G;
        // std::cout << v << std::endl;

        // line search Armijo
        t = 1.0;
        X1 = X0 + t * v;
        while (L_reflection(X1, nullptr, nullptr, objfn_data, nullptr) >
               f + alpha1 * t * dot(G, v)) {
            t = beta * t;
            X1 = X0 + t * v;
            if (t < 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5)
                break;
        }

        X0 = X1;

        // gradient norm small enough, return
        if (norm(G) < 0.00001) {
            *X0_ = X0;
            if (H_out != nullptr)
                *H_out = H;
            L_reflection(X0, nullptr, nullptr, objfn_data, f_out);
            return true;
        }
    }
    *X0_ = X0;
    return false;
}
template <typename Float, typename Spectrum>
Array<Float, 4>
FermatNEE<Float, Spectrum>::solve_4x4_enoki(const Matrix4f &m,
                                            const Array<Float, 4> b) const {
    using Vector = Array<Float, 4>;

    Vector col0 = m.coeff(0), col1 = m.coeff(1), col2 = m.coeff(2),
           col3 = m.coeff(3);

    col1 = shuffle<2, 3, 0, 1>(col1);
    col3 = shuffle<2, 3, 0, 1>(col3);

    Vector tmp, row0, row1, row2, row3;

    tmp = shuffle<1, 0, 3, 2>(col2 * col3);
    row0 = col1 * tmp;
    row1 = col0 * tmp;
    tmp = shuffle<2, 3, 0, 1>(tmp);
    row0 = fmsub(col1, tmp, row0);
    row1 = shuffle<2, 3, 0, 1>(fmsub(col0, tmp, row1));

    tmp = shuffle<1, 0, 3, 2>(col1 * col2);
    row0 = fmadd(col3, tmp, row0);
    row3 = col0 * tmp;
    tmp = shuffle<2, 3, 0, 1>(tmp);
    row0 = fnmadd(col3, tmp, row0);
    row3 = shuffle<2, 3, 0, 1>(fmsub(col0, tmp, row3));

    tmp = shuffle<1, 0, 3, 2>(shuffle<2, 3, 0, 1>(col1) * col3);
    col2 = shuffle<2, 3, 0, 1>(col2);
    row0 = fmadd(col2, tmp, row0);
    row2 = col0 * tmp;
    tmp = shuffle<2, 3, 0, 1>(tmp);
    row0 = fnmadd(col2, tmp, row0);
    row2 = shuffle<2, 3, 0, 1>(fmsub(col0, tmp, row2));

    tmp = shuffle<1, 0, 3, 2>(col0 * col1);
    row2 = fmadd(col3, tmp, row2);
    row3 = fmsub(col2, tmp, row3);
    tmp = shuffle<2, 3, 0, 1>(tmp);
    row2 = fmsub(col3, tmp, row2);
    row3 = fnmadd(col2, tmp, row3);

    tmp = shuffle<1, 0, 3, 2>(col0 * col3);
    row1 = fnmadd(col2, tmp, row1);
    row2 = fmadd(col1, tmp, row2);
    tmp = shuffle<2, 3, 0, 1>(tmp);
    row1 = fmadd(col2, tmp, row1);
    row2 = fnmadd(col1, tmp, row2);

    tmp = shuffle<1, 0, 3, 2>(col0 * col2);
    row1 = fmadd(col3, tmp, row1);
    row3 = fnmadd(col1, tmp, row3);
    tmp = shuffle<2, 3, 0, 1>(tmp);
    row1 = fnmadd(col3, tmp, row1);
    row3 = fmadd(col1, tmp, row3);

    Matrix inv_times_det = Matrix4f(row0, row1, row2, row3);

    Float inv_det = (rcp(dot(col0, row0)));

    return inv_det * (inv_times_det * b);
}

template <typename Float, typename Spectrum>
Array<double, 4> FermatNEE<Float, Spectrum>::solve_4x4_enoki_double(
    const Matrix<double, 4> &m, const Array<double, 4> b) const {
    using Vector = Array<double, 4>;

    Vector col0 = m.coeff(0), col1 = m.coeff(1), col2 = m.coeff(2),
           col3 = m.coeff(3);

    col1 = shuffle<2, 3, 0, 1>(col1);
    col3 = shuffle<2, 3, 0, 1>(col3);

    Vector tmp, row0, row1, row2, row3;

    tmp = shuffle<1, 0, 3, 2>(col2 * col3);
    row0 = col1 * tmp;
    row1 = col0 * tmp;
    tmp = shuffle<2, 3, 0, 1>(tmp);
    row0 = fmsub(col1, tmp, row0);
    row1 = shuffle<2, 3, 0, 1>(fmsub(col0, tmp, row1));

    tmp = shuffle<1, 0, 3, 2>(col1 * col2);
    row0 = fmadd(col3, tmp, row0);
    row3 = col0 * tmp;
    tmp = shuffle<2, 3, 0, 1>(tmp);
    row0 = fnmadd(col3, tmp, row0);
    row3 = shuffle<2, 3, 0, 1>(fmsub(col0, tmp, row3));

    tmp = shuffle<1, 0, 3, 2>(shuffle<2, 3, 0, 1>(col1) * col3);
    col2 = shuffle<2, 3, 0, 1>(col2);
    row0 = fmadd(col2, tmp, row0);
    row2 = col0 * tmp;
    tmp = shuffle<2, 3, 0, 1>(tmp);
    row0 = fnmadd(col2, tmp, row0);
    row2 = shuffle<2, 3, 0, 1>(fmsub(col0, tmp, row2));

    tmp = shuffle<1, 0, 3, 2>(col0 * col1);
    row2 = fmadd(col3, tmp, row2);
    row3 = fmsub(col2, tmp, row3);
    tmp = shuffle<2, 3, 0, 1>(tmp);
    row2 = fmsub(col3, tmp, row2);
    row3 = fnmadd(col2, tmp, row3);

    tmp = shuffle<1, 0, 3, 2>(col0 * col3);
    row1 = fnmadd(col2, tmp, row1);
    row2 = fmadd(col1, tmp, row2);
    tmp = shuffle<2, 3, 0, 1>(tmp);
    row1 = fmadd(col2, tmp, row1);
    row2 = fnmadd(col1, tmp, row2);

    tmp = shuffle<1, 0, 3, 2>(col0 * col2);
    row1 = fmadd(col3, tmp, row1);
    row3 = fnmadd(col1, tmp, row3);
    tmp = shuffle<2, 3, 0, 1>(tmp);
    row1 = fnmadd(col3, tmp, row1);
    row3 = fmadd(col1, tmp, row3);

    Matrix<double, 4> inv_times_det = Matrix<double, 4>(row0, row1, row2, row3);

    double inv_det = (rcp(dot(col0, row0)));

    return inv_det * (inv_times_det * b);
}

template <typename Float, typename Spectrum>
bool FermatNEE<Float, Spectrum>::newton_solver_SMS(
    const SurfaceInteraction3f &si, const EmitterInteraction &ei,
    L_data_float *data) {
    ScopedPhase scope_phase(ProfilerPhase::SMS_NewtonSolver);

    // Newton iterations..
    bool success = false;
    size_t iterations = 0;
    Float beta = 1.0f;

    bool needs_step_update = true;

    // init
    m_proposed_positions.clear();
    for (int i = 0; i < m_seed_path.size(); i++) {
        m_proposed_positions.push_back(m_seed_path[i].p);
    }

    while (iterations < m_config.max_iterations) {
        bool step_success = true;
        if (needs_step_update) {
            if (m_config.halfvector_constraints) {
                // Use standard manifold formulation using half-vector
                // constraints
                step_success = compute_step_halfvector(si.p, ei);
            } else {
                // Use angle-difference constraint formulation
                step_success = compute_step_anglediff(si.p, ei);
            }
        }
        if (!step_success) {
            break;
        }

        // Check for success
        bool converged = true;
        for (size_t i = 0; i < m_current_path.size(); ++i) {
            const ManifoldVertex &v = m_current_path[i];
            if (norm(v.C) > m_config.solver_threshold) {
                converged = false;
                break;
            }
        }
        if (converged) {
            success = true;
            break;
        }

        // Make a proposal
        m_proposed_positions.clear();
        for (size_t i = 0; i < m_current_path.size(); ++i) {
            const ManifoldVertex &v = m_current_path[i];
            Point3f p_prop = v.p - m_config.step_scale * beta *
                                       (v.dp_du * v.dx[0] + v.dp_dv * v.dx[1]);
            m_proposed_positions.push_back(p_prop);
        }

        // Project back to surfaces
        bool project_success = reproject_raytrace(si);

        // orthogonal reprojection, just for the sake of showing that it doesn't
        // work bool project_success = reproject(si, data);

        if (!project_success) {
            beta = 0.5f * beta;
            needs_step_update = false;
        } else {
            beta = std::min(Float(1), Float(2) * beta);
            m_current_path = m_proposed_path;
            needs_step_update = true;
        }

        iterations++;
    }

    if (!success) {
        return false;
    }

    /// from manifold_ms
    /* "In the refraction case, the half-vector formulation of Manifold
       walks will often converge to invalid solutions that are actually
       reflections. Here we need to reject those". */
    size_t n = m_current_path.size();
    for (size_t i = 0; i < n; ++i) {
        Point3f x_prev = (i == 0) ? si.p : m_current_path[i - 1].p;
        Point3f x_next = (i == n - 1) ? ei.p : m_current_path[i + 1].p;
        Point3f x_cur = m_current_path[i].p;

        bool at_endpoint_with_fixed_direction =
            (i == (n - 1) && ei.is_directional());
        Vector3f wi = normalize(x_prev - x_cur);
        Vector3f wo =
            at_endpoint_with_fixed_direction ? ei.d : normalize(x_next - x_cur);

        Float cos_theta_i = dot(m_current_path[i].gn, wi),
              cos_theta_o = dot(m_current_path[i].gn, wo);
        bool refraction = cos_theta_i * cos_theta_o < 0.f,
             reflection = !refraction;
        if ((m_current_path[i].eta == 1.f && !reflection) ||
            (m_current_path[i].eta != 1.f && !refraction)) {
            return false;
        }
    }

    return true;
}

template <typename Float, typename Spectrum>
typename FermatNEE<Float, Spectrum>::Mask
FermatNEE<Float, Spectrum>::invert_tridiagonal_step(
    std::vector<ManifoldVertex> &v) {
    // Solve block tri-diagonal linear system with full RHS vector

    // From "The Natural-Constraint Representation of the Path Space for
    // Efficient Light Transport Simulation" by Kaplanyan et al. 2014
    // Supplemental material, Figure 2.

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

    v[0].tmp = v[0].dC_dx_prev;
    Matrix2f m = v[0].dC_dx_cur;
    if (!invert(m, v[0].inv_lambda))
        return false;

    for (int i = 1; i < n; ++i) {
        v[i].tmp = v[i].dC_dx_prev * v[i - 1].inv_lambda;
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

template <typename Float, typename Spectrum>
inline Float FermatNEE<Float, Spectrum>::L_flat(Array<double, 4> x_in,
                                                Array<double, 4> *grad_out,
                                                Matrix<double, 4> *hess_out,
                                                L_data_float *objfn_data,
                                                Float *f_out) {
    ScopedPhase scope_phase(ProfilerPhase::FermaNEE_computeDirection);

    /* +---------------------+
    // |extract geometry data|
    // +---------------------+ */
    double e = objfn_data->e;
    double n1 = objfn_data->n1;
    double n2 = objfn_data->n2;
    Point3f O = objfn_data->O;
    Point3f S = objfn_data->S;

    // local parametrisation [0,1]x[0,1]
    Point2f x1(x_in.x(), x_in.y());
    Point2f x2(x_in.z(), x_in.w());

    // local parametrisation [0,L]x[0,H]
    double u1 = x_in.x() * objfn_data->L1;
    double v1 = x_in.y() * objfn_data->W1;
    double u2 = x_in.z() * objfn_data->L2;
    double v2 = x_in.w() * objfn_data->W2;

    double f = 0.f;
    double g = 0.f;
    if (f_out)
        *f_out = f;

    // vertex positions
    Point3f P1(f, u1, v1);
    Point3f P2(g + e, u2, v2);

    // norms
    double normP1P2 = norm(P1 - P2);
    double normP2S = norm(S - P2);
    double normOP1 = norm(P1 - O);

    /* +--------------------+
    // |Gradient and Hessian|
    // +--------------------+ */
    if (grad_out || hess_out) {
        // first order curvature

        // precompute some terms
        double D = rcp(normP2S);
        double E = rcp(normP1P2);
        double F = rcp(normOP1);

        double A = E * E * E;
        double B = F * F * F;
        double C = D * D * D;

        double G = e;
        double H = S.x() - e;
        double I = O.x();
        double J = v1 - v2;
        double K = u1 - u2;

        /* +--------------------+
        // |Gradient computation|
        // +--------------------+ */
        if (grad_out) {

            (*grad_out)[0] = E * n2 * (K) + F * n1 * (u1 - O.y());
            (*grad_out)[1] = E * n2 * (J) + F * n1 * (v1 - O.z());
            (*grad_out)[2] = D * n1 * (u2 - S.y()) + E * n2 * (-K);
            (*grad_out)[3] = D * n1 * (v2 - S.z()) + E * n2 * (-J);
        }
        /* +-------------------+
        // |Hessian computation|
        // +-------------------+ */
        if (hess_out) {
            // Hessian coeficients
            (*hess_out)[0][0] = A * n2 * (K) * (-K) +
                                B * n1 * (u1 - O.y()) * (O.y() - u1) + E * n2 +
                                F * n1;
            (*hess_out)[0][1] =
                A * n2 * (K) * (-J) + B * n1 * (-O.y() + u1) * (O.z() - v1);
            (*hess_out)[0][2] = A * n2 * (K) * (K) -E * n2;
            (*hess_out)[0][3] = A * n2 * (K) * (J);

            (*hess_out)[1][1] = A * n2 * (J) * (-J) +
                                B * n1 * (-O.z() + v1) * (O.z() - v1) + E * n2 +
                                F * n1;
            (*hess_out)[1][2] = A * n2 * (K) * (J);
            (*hess_out)[1][3] = A * n2 * (J) * (J) -E * n2;

            (*hess_out)[2][2] = A * n2 * (K) * (-K) +
                                C * n1 * (u2 - S.y()) * (S.y() - u2) + D * n1 +
                                E * n2;
            (*hess_out)[2][3] =
                A * n2 * (-K) * (J) + C * n1 * (u2 - S.y()) * (S.z() - v2);
            (*hess_out)[3][0] = A * n2 * (-K) * (-J);

            (*hess_out)[3][3] = A * n2 * (J) * (-J) +
                                C * n1 * (v2 - S.z()) * (S.z() - v2) + D * n1 +
                                E * n2;

            // symmetric
            (*hess_out)[1][0] = (*hess_out)[0][1];
            (*hess_out)[2][0] = (*hess_out)[0][2];
            (*hess_out)[2][1] = (*hess_out)[1][2];
            (*hess_out)[3][1] = (*hess_out)[1][3];
            (*hess_out)[3][2] = (*hess_out)[2][3];
        }
    }
    return n1 * (normOP1 + normP2S) + n2 * normP1P2;
}

template <typename Float, typename Spectrum>
inline double FermatNEE<Float, Spectrum>::L_double(Array<double, 4> x_in,
                                                   Array<double, 4> *grad_out,
                                                   Matrix<double, 4> *hess_out,
                                                   L_data_float *objfn_data,
                                                   Float *f_out) {
    ScopedPhase scope_phase(ProfilerPhase::FermaNEE_computeDirection);

    // if(objfn_data->H1 == 0 && objfn_data->H2 ==0)
    //     return L_flat(x_in, grad_out, hess_out, objfn_data, f_out);
    double h = 1.0;

    /* +---------------------+
    // |extract geometry data|
    // +---------------------+ */
    double e = objfn_data->e;
    double n1 = objfn_data->n1;
    double n2 = objfn_data->n2;
    Array<double, 3> O = Array<double, 3>(objfn_data->O);
    Array<double, 3> S = Array<double, 3>(objfn_data->S);

    // local parametrisation [0,1]x[0,1]
    Point2f x1(x_in.x(), x_in.y());
    Point2f x2(x_in.z(), x_in.w());

    // local parametrisation [0,L]x[0,H]
    double u1 = x_in.x() * objfn_data->L1;
    double v1 = x_in.y() * objfn_data->W1;
    double u2 = x_in.z() * objfn_data->L2;
    double v2 = x_in.w() * objfn_data->W2;

    double f = h * objfn_data->pH1->eval_texture(0, x1);
    double g = h * objfn_data->pH2->eval_texture(0, x2);
    if (f_out)
        *f_out = f;

    // vertex positions
    Array<double, 3> P1(f, u1, v1);
    Array<double, 3> P2(g + e, u2, v2);

    // norms
    double normP1P2 = norm(P1 - P2);
    double normP2S = norm(S - P2);
    double normOP1 = norm(P1 - O);

    // +--------------------+
    // |Gradient and Hessian|
    // +--------------------+
    if (grad_out || hess_out) {
        // first order curvature
        double duf = h * ((objfn_data->pH1->eval_texture(1, x1)));
        double dvf = h * ((objfn_data->pH1->eval_texture(2, x1)));
        double dug = h * ((objfn_data->pH2->eval_texture(1, x2)));
        double dvg = h * ((objfn_data->pH2->eval_texture(2, x2)));
        // std::cout << x1 << "|" << x2 <<"  ||"<< duf <<"|"<< dvf <<"|"<< dug
        // <<"|"<< dvg <<"|" << std::endl; precompute some terms
        double D = rsqrt(squared_norm(S - P2));
        double E = rsqrt(squared_norm(P1 - P2));
        double F = rsqrt(squared_norm(P1 - O));

        double A = sqr(E) * E;
        double B = sqr(F) * F;
        double C = sqr(D) * D;

        double G = (e - f + g);
        double H = (S.x() - e - g);
        double I = O.x() - f;
        double J = v1 - v2;
        double K = u1 - u2;

        double duf2 = sqr(duf);
        double dug2 = sqr(dug);
        double dvf2 = sqr(dvf);
        double dvg2 = sqr(dvg);

        /* +--------------------+
        // |Gradient computation|
        // +--------------------+ */
        if (grad_out) {

            (*grad_out)[0] = E * n2 * fnmadd(G, duf, K) +
                             F * n1 * fnmadd(I, duf, u1 - O.y());
            (*grad_out)[1] = E * n2 * fnmadd(G, dvf, J) +
                             F * n1 * fnmadd(I, dvf, v1 - O.z());
            (*grad_out)[2] =
                D * n1 * fnmadd(H, dug, u2 - S.y()) + E * n2 * fmsub(G, dug, K);
            (*grad_out)[3] =
                D * n1 * fnmadd(H, dvg, v2 - S.z()) + E * n2 * fmsub(G, dvg, J);
        }
        /* +-------------------+
        // |Hessian computation|
        // +-------------------+ */
        if (hess_out) {

            // second order curvature information
            double duuf = h * ((objfn_data->pH1->eval_texture(3, x1)));
            double dvvf = h * ((objfn_data->pH1->eval_texture(4, x1)));
            double duvf = h * ((objfn_data->pH1->eval_texture(5, x1)));
            double duug = h * ((objfn_data->pH2->eval_texture(3, x2)));
            double dvvg = h * ((objfn_data->pH2->eval_texture(4, x2)));
            double duvg = h * ((objfn_data->pH2->eval_texture(5, x2)));
            // Hessian coeficients

            double a = fnmadd(G, duf, K);
            double b = fnmadd(G, dvg, J);
            double c = fnmadd(G, dug, K);
            double d = fnmadd(G, dvf, J);
            double aa = fmsub(G, dug, K);
            double bb = fnmadd(H, dug, u2 - S.y());
            double cc = fmadd(H, dvg, S.z() - v2);

            (*hess_out)[0][0] = A * n2 * (-sqr(a)) +
                                B * n1 * fnmadd(I, duf, u1 - O.y()) *
                                    fnmadd(I, duf, O.y() - u1) +
                                E * n2 * fnmadd(G, duuf, duf2 + 1.0) +
                                F * n1 * fnmadd(I, duuf, duf2 + 1.0);

            (*hess_out)[0][1] = A * n2 * a * fmsub(G, dvf, J) +
                                B * n1 * fnmadd(I, duf, u1 - O.y()) *
                                    fnmadd(I, dvf, O.z() - v1) +
                                E * n2 * fnmadd(G, duvf, duf * dvf) +
                                F * n1 * fnmadd(I, duvf, duf * dvf);

            (*hess_out)[0][3] = A * n2 * a * b - E * duf * dvg * n2;

            (*hess_out)[0][2] = A * n2 * a * c + E * n2 * fnmsub(duf, dug, 1.0);

            (*hess_out)[1][1] = A * n2 * d * fmsub(G, dvf, J) +
                                B * n1 * fnmadd(I, dvf, v1 - O.z()) *
                                    fnmadd(I, dvf, O.z() - v1) +
                                E * n2 * fnmadd(G, dvvf, dvf2 + 1.0) +
                                F * n1 * fnmadd(I, dvvf, dvf2 + 1.0);

            (*hess_out)[1][2] = A * n2 * c * d - E * dug * dvf * n2;

            (*hess_out)[1][3] = A * n2 * d * b + E * n2 * fnmsub(dvf, dvg, 1.0);

            (*hess_out)[2][2] = A * n2 * c * aa +
                                C * n1 * bb * fmadd(H, dug, S.y() - u2) +
                                D * n1 * fnmadd(H, duug, dug2 + 1.0) +
                                E * n2 * fmadd(G, duug, dug2 + 1.0);

            (*hess_out)[2][3] = A * n2 * aa * b + C * n1 * bb * cc +
                                D * n1 * fnmadd(H, duvg, dug * dvg) +
                                E * n2 * (G * duvg + dug * dvg);

            (*hess_out)[3][0] = A * n2 * (a) * (b) -E * duf * dvg * n2;

            (*hess_out)[3][3] = A * n2 * (-sqr(b)) + C * n1 * (-sqr(cc)) +
                                D * n1 * fnmadd(H, dvvg, dvg2 + 1.0) +
                                E * n2 * fmadd(G, dvvg, dvg2 + 1.0);

            // symmetric
            (*hess_out)[1][0] = (*hess_out)[0][1];
            (*hess_out)[2][0] = (*hess_out)[0][2];
            (*hess_out)[2][1] = (*hess_out)[1][2];
            (*hess_out)[3][1] = (*hess_out)[1][3];
            (*hess_out)[3][2] = (*hess_out)[2][3];
        }

        // std::cout << "L=" << n1 * (normOP1 + normP2S) + n2 * normP1P2
        // <<std::endl; std::cout << "L_grad=" << *grad_out <<std::endl;
        // std::cout << "L_HESS=" << *hess_out <<std::endl;
    }
    return n1 * (normOP1 + normP2S) + n2 * normP1P2;
}

template <typename Float, typename Spectrum>
inline double static FermatNEE<Float, Spectrum>::L_reflection(
    Vector2f x_in, Vector2f *grad_out, Matrix<double, 2> *hess_out,
    L_data_float *objfn_data, Float *f_out) {
    ScopedPhase scope_phase(ProfilerPhase::FermaNEE_computeDirection);

    double h = 1.0;

    /* +---------------------+
    // |extract geometry data|
    // +---------------------+ */
    double n1 = objfn_data->n1;
    Array<double, 3> O = Array<double, 3>(objfn_data->O);
    Array<double, 3> S = Array<double, 3>(objfn_data->S);

    // local parametrisation [0,1]x[0,1]
    Point2f x1(x_in.x(), x_in.y());

    // local parametrisation [0,L]x[0,H]
    double u1 = x_in.x() * objfn_data->L1;
    double v1 = x_in.y() * objfn_data->W1;

    double f = h * objfn_data->pH1->eval_texture(0, x1);
    if (f_out)
        *f_out = f;

    // vertex positions
    Array<double, 3> P1(f, u1, v1);

    // norms
    double normP1S = norm(S - P1);
    double normOP1 = norm(P1 - O);

    // +--------------------+
    // |Gradient and Hessian|
    // +--------------------+
    if (grad_out || hess_out) {
        // first order curvature
        double duf = h * ((objfn_data->pH1->eval_texture(1, x1)));
        double dvf = h * ((objfn_data->pH1->eval_texture(2, x1)));

        // std::cout << x1 << "|" << x2 <<"  ||"<< duf <<"|"<< dvf <<"|"<< dug
        // <<"|"<< dvg <<"|" << std::endl; precompute some terms
        double D = rsqrt(squared_norm(S - P1));
        double F = rsqrt(squared_norm(P1 - O));

        double B = sqr(F) * F;
        double C = sqr(D) * D;

        double I = O.x() - f;

        double duf2 = sqr(duf);
        double dvf2 = sqr(dvf);

        /* +--------------------+
        // |Gradient computation|
        // +--------------------+ */
        if (grad_out) {
            (*grad_out)[0] = n1 * ((-O.y() - duf * I + u1) * F +
                                   (-S.y() - duf * (S.x() - f) + u1) * D);
            (*grad_out)[1] = n1 * ((-O.z() - dvf * I + v1) * F +
                                   (-S.z() - dvf * (S.x() - f) + v1) * D);
        }
        /* +-------------------+
        // |Hessian computation|
        // +-------------------+ */
        if (hess_out) {

            // second order curvature information
            double duuf = h * ((objfn_data->pH1->eval_texture(3, x1)));
            double dvvf = h * ((objfn_data->pH1->eval_texture(4, x1)));
            double duvf = h * ((objfn_data->pH1->eval_texture(5, x1)));
            // Hessian coeficients

            (*hess_out)[0][0] =
                n1 * (B * (-I * duf - O.y() + u1) * (I * duf + O.y() - u1) +
                      F * (-I * duuf + duf2 + 1) +
                      (-S.y() - duf * (S.x() - f) + u1) *
                          (S.y() + duf * (S.x() - f) - u1) /
                          pow(pow(S.x() - f, 2) + pow(S.y() - u1, 2) +
                                  pow(S.z() - v1, 2),
                              3.0 / 2.0) +
                      (duf2 + duuf * (-S.x() + f) + 1) /
                          sqrt(pow(S.x() - f, 2) + pow(S.y() - u1, 2) +
                               pow(S.z() - v1, 2)));
            (*hess_out)[0][1] =
                n1 * (B * (I * duf + O.y() - u1) * (-I * dvf - O.z() + v1) +
                      F * (-I * duvf + duf * dvf) +
                      (duf * dvf + duvf * (-S.x() + f)) /
                          sqrt(pow(S.x() - f, 2) + pow(S.y() - u1, 2) +
                               pow(S.z() - v1, 2)) +
                      (S.y() + duf * (S.x() - f) - u1) *
                          (-S.z() - dvf * (S.x() - f) + v1) /
                          pow(pow(S.x() - f, 2) + pow(S.y() - u1, 2) +
                                  pow(S.z() - v1, 2),
                              3.0 / 2.0));
            (*hess_out)[1][0] =
                n1 * (B * (-I * duf - O.y() + u1) * (I * dvf + O.z() - v1) +
                      F * (-I * duvf + duf * dvf) +
                      (duf * dvf + duvf * (-S.x() + f)) /
                          sqrt(pow(S.x() - f, 2) + pow(S.y() - u1, 2) +
                               pow(S.z() - v1, 2)) +
                      (-S.y() - duf * (S.x() - f) + u1) *
                          (S.z() + dvf * (S.x() - f) - v1) /
                          pow(pow(S.x() - f, 2) + pow(S.y() - u1, 2) +
                                  pow(S.z() - v1, 2),
                              3.0 / 2.0));
            (*hess_out)[1][1] =
                n1 * (B * (-I * dvf - O.z() + v1) * (I * dvf + O.z() - v1) +
                      F * (-I * dvvf + dvf2 + 1) +
                      (-S.z() - dvf * (S.x() - f) + v1) *
                          (S.z() + dvf * (S.x() - f) - v1) /
                          pow(pow(S.x() - f, 2) + pow(S.y() - u1, 2) +
                                  pow(S.z() - v1, 2),
                              3.0 / 2.0) +
                      (dvf2 + dvvf * (-S.x() + f) + 1) /
                          sqrt(pow(S.x() - f, 2) + pow(S.y() - u1, 2) +
                               pow(S.z() - v1, 2)));
        }
    }
    return n1 * (normOP1 + normP1S);
}

/* +===============================================================================================+
// |                                        SPECULAR OUTPUT |
//
+===============================================================================================+
*/

template <typename Float, typename Spectrum>
Float FermatNEE<Float, Spectrum>::invert_tridiagonal_geo(
    std::vector<ManifoldVertex> &v) const {
    // Solve block tri-diagonal linear system with RHS vector where only last
    // element in non-zero

    // Procedure as outlined in original "Manifold Exploration" by Jakob and
    // Marschner 2012. Based on the implementation in Mitsuba 0.6: manifold.cpp
    // line 382

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
        return 0.f;

    Matrix2f Li;
    if (!invert(v[0].dC_dx_cur, Li))
        return 0.f;

    for (int i = 0; i < n - 2; ++i) {
        v[i].tmp = Li * v[i].dC_dx_next;
        Matrix2f m = v[i + 1].dC_dx_cur - v[i + 1].dC_dx_prev * v[i].tmp;
        if (!invert(m, Li))
            return 0.f;
    }

    v[n - 2].inv_lambda = -Li * v[n - 2].dC_dx_next;
    for (int i = n - 3; i >= 0; --i) {
        v[i].inv_lambda = -v[i].tmp * v[i + 1].inv_lambda;
    }

    return abs(det(-v[0].inv_lambda));
}

template <typename Float, typename Spectrum>
Float FermatNEE<Float, Spectrum>::generalizedGeometryFactor(
    const SurfaceInteraction3f &si0, const SurfaceInteraction3f &si1,
    const SurfaceInteraction3f &si2, const SurfaceInteraction3f &si3,
    const bool reflection) const {
    Float n1 = 1.0f, n2 = 1.54f;
    bool smoothing = false;
    std::vector<ManifoldVertex> v;
    ManifoldVertex vx(si0);
    vx.make_orthonormal();
    ManifoldVertex vy(si3);
    vy.make_orthonormal();

    if (!reflection) {

        // X1
        ManifoldVertex x1(si1);
        x1.make_orthonormal();

        // X2
        ManifoldVertex x2(si2);
        x2.make_orthonormal();

        // v.push_back(vx);
        v.push_back(x1);
        v.push_back(x2);
        v.push_back(vy);
    } else {
        // vertex:

        // X1
        ManifoldVertex x1(si1);
        x1.make_orthonormal();

        v.push_back(x1);
        v.push_back(vy);
    }

    size_t k = v.size();

    for (size_t i = 0; i < k - 1; ++i) {
        v[i].dC_dx_prev = Matrix2f(0.f);
        v[i].dC_dx_cur = Matrix2f(0.f);
        v[i].dC_dx_next = Matrix2f(0.f);

        Point3f x_cur = v[i].p;
        Point3f x_next = v[i + 1].p;

        Vector3f wo = x_next - x_cur;
        Float ilo = norm(wo);
        if (ilo < 1e-3f) {
            return false;
        }
        ilo = rcp(ilo);
        wo *= ilo;

        if (v[i].fixed_direction) {
            // Derivative of directional constraint w.r.t. x_{i}
            Vector3f dc_du_cur = -ilo * (v[i].dp_du - wo * dot(wo, v[i].dp_du)),
                     dc_dv_cur = -ilo * (v[i].dp_dv - wo * dot(wo, v[i].dp_dv));
            v[i].dC_dx_cur = Matrix2f(
                dot(dc_du_cur, v[i].dp_du), dot(dc_dv_cur, v[i].dp_du),
                dot(dc_du_cur, v[i].dp_dv), dot(dc_dv_cur, v[i].dp_dv));

            // Derivative of directional constraint w.r.t. x_{i+1}
            Vector3f dc_du_next =
                         ilo * (v[i + 1].dp_du - wo * dot(wo, v[i + 1].dp_du)),
                     dc_dv_next =
                         ilo * (v[i + 1].dp_dv - wo * dot(wo, v[i + 1].dp_dv));
            v[i].dC_dx_next = Matrix2f(
                dot(dc_du_next, v[i].dp_du), dot(dc_dv_next, v[i].dp_du),
                dot(dc_du_next, v[i].dp_dv), dot(dc_dv_next, v[i].dp_dv));
            continue;
        }

        Point3f x_prev =
            (i == 0) ? vx.p
                     : v[i - 1].p; // Note that we only end up here for
                                   // positionally fixed endpoints, thus x is
                                   // not part of the path array directly.

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

        // Local shading tangent frame
        Float dot_dpdu_n = dot(v[i].dp_du, v[i].n),
              dot_dpdv_n = dot(v[i].dp_dv, v[i].n);
        Vector3f s = v[i].dp_du - dot_dpdu_n * v[i].n,
                 t = v[i].dp_dv - dot_dpdv_n * v[i].n;

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

            v[i].dC_dx_prev = Matrix2f(dot(s, dh_du), dot(s, dh_dv),
                                       dot(t, dh_du), dot(t, dh_dv));
        }

        // Derivative of specular constraint w.r.t. x_{i}
        dh_du = -v[i].dp_du * (ili + ilo) + wi * (dot(wi, v[i].dp_du) * ili) +
                wo * (dot(wo, v[i].dp_du) * ilo);
        dh_dv = -v[i].dp_dv * (ili + ilo) + wi * (dot(wi, v[i].dp_dv) * ili) +
                wo * (dot(wo, v[i].dp_dv) * ilo);
        dh_du -= h * dot(dh_du, h);
        dh_dv -= h * dot(dh_dv, h);
        if (eta != 1.f) {
            dh_du *= -1.f;
            dh_dv *= -1.f;
        }

        Float dot_h_n = dot(h, v[i].n), dot_h_dndu = dot(h, v[i].dn_du),
              dot_h_dndv = dot(h, v[i].dn_dv);

        v[i].dC_dx_cur =
            Matrix2f(dot(dh_du, s) - dot(v[i].dp_du, v[i].dn_du) * dot_h_n -
                         dot_dpdu_n * dot_h_dndu,
                     dot(dh_dv, s) - dot(v[i].dp_du, v[i].dn_dv) * dot_h_n -
                         dot_dpdu_n * dot_h_dndv,
                     dot(dh_du, t) - dot(v[i].dp_dv, v[i].dn_du) * dot_h_n -
                         dot_dpdv_n * dot_h_dndu,
                     dot(dh_dv, t) - dot(v[i].dp_dv, v[i].dn_dv) * dot_h_n -
                         dot_dpdv_n * dot_h_dndv);

        // Derivative of specular constraint w.r.t. x_{i+1}
        dh_du = ilo * (v[i + 1].dp_du - wo * dot(wo, v[i + 1].dp_du));
        dh_dv = ilo * (v[i + 1].dp_dv - wo * dot(wo, v[i + 1].dp_dv));
        dh_du -= h * dot(dh_du, h);
        dh_dv -= h * dot(dh_dv, h);
        if (eta != 1.f) {
            dh_du *= -1.f;
            dh_dv *= -1.f;
        }

        v[i].dC_dx_next = Matrix2f(dot(s, dh_du), dot(s, dh_dv), dot(t, dh_du),
                                   dot(t, dh_dv));
    }

    if (vy.fixed_direction) {
        Float G = invert_tridiagonal_geo(v);
        // Cancel out cosine term that will be added during lighting integral in
        // integrator
        Vector3f d = normalize(v[k - 1].p - v[k - 2].p);
        G /= abs_dot(d, v[k - 1].n);
        return G;
    } else {
        Float dx1_dxend = invert_tridiagonal_geo(v);
        // std::cout<<"det = "<<dx1_dxend<<std::endl;

        /* Unfortunately, these geometric terms can be unstable, so to avoid
           severe variance we can clamp here. */
        dx1_dxend = min(dx1_dxend, Float(2.f));

        Vector3f d = vx.p - v[0].p;
        Float inv_r2 = rcp(squared_norm(d));
        d *= sqrt(inv_r2);
        Float dw0_dx1 = abs_dot(d, v[0].gn) * inv_r2;
        Float G = dw0_dx1 * dx1_dxend;
        // std::cout<<"geometry factor = "<<G<<std::endl;
        return G;
    }
}

template <typename Float, typename Spectrum>
Spectrum FermatNEE<Float, Spectrum>::compute_ray_contribution(
    SurfaceInteraction3f si_O, const BSDFContext &ctx, Ray3f ray,
    EmitterInteraction vy, Float n1, Float n2, ShapePtr H1, ShapePtr H2,
    Mask active) const {
    MTS_MASK_ARGUMENT(active);

    Point3f S = vy.p;
    Spectrum specular_out(0.0f);
    SurfaceInteraction3f si_heightfield_1 = m_scene->ray_intersect(ray, active);

    // +------------------------------+
    // |first Heightfield intersection:|
    // +------------------------------+
    if (si_heightfield_1.shape == H1 &&
        si_heightfield_1.shape->is_caustic_caster_double_refraction()) {

        //  std::tie(scatter_success, wo) = SpecularManifold::refract(wi, m,
        //  vertex.eta);
        auto [scatter_success_1, wo_2] =
            refract(-ray.d, si_heightfield_1.n, n2);
        if (!scatter_success_1)
            return specular_out;
        // Ray3f ray_2 =
        //     Ray3f(si_heightfield_1.p, wo_2, math::RayEpsilon<Float> * (1.f +
        //     hmax(abs(si_heightfield_1.p))),
        //           norm(vy.p - si_heightfield_1.p) * (1.f -
        //           math::RayEpsilon<Float>), si_heightfield_1.time,
        //           si_heightfield_1.wavelengths);

        auto [bs, bsdf_val] = si_heightfield_1.bsdf()->sample(
            ctx, si_heightfield_1, m_sampler->next_1d(), m_sampler->next_2d(),
            active);
        RayDifferential3f ray_2 =
            si_heightfield_1.spawn_ray(si_heightfield_1.to_world(bs.wo));
        SurfaceInteraction3f si_heightfield_2 =
            m_scene->ray_intersect(ray_2, active);

        // +-------------------------------+
        // |second heightfield intersection:|
        // +-------------------------------+
        if (si_heightfield_2.shape == H2 &&
            si_heightfield_2.shape->is_caustic_caster_double_refraction()) {

            auto [scatter_success_2, wo_final] =
                refract(-ray_2.d, si_heightfield_2.n, n2);
            if (!scatter_success_2)
                return Spectrum(0.0f);

            Vector3f d_tmp = normalize(vy.p - si_heightfield_2.p);
            // Ray3f ray_final(si_heightfield_2.p, d_tmp, ray.time,
            // ray.wavelengths);

            Ray3f ray_final = Ray3f(
                si_heightfield_2.p, normalize(vy.p - si_heightfield_2.p),
                math::RayEpsilon<Float> * (1.f + hmax(abs(si_heightfield_2.p))),
                norm(vy.p - si_heightfield_2.p) *
                    (1.f + math::RayEpsilon<Float>),
                ray.time, ray.wavelengths);
            /* +-------------------+
            // |emitter interaction|
            // +-------------------+ */
            // check if we reached the sampled point (e.g: is the sampled point
            // on the backface of the emitter?)
            Ray3f ray_final_test = Ray3f(
                si_heightfield_2.p, normalize(vy.p - si_heightfield_2.p),
                math::RayEpsilon<Float> * (1.f + hmax(abs(si_heightfield_2.p))),
                norm(vy.p - si_heightfield_2.p) *
                    (1.f - math::RayEpsilon<Float>),
                ray.time, ray.wavelengths);
            SurfaceInteraction3f si_test_final =
                m_scene->ray_intersect(ray_final_test, active);
            if (si_test_final.shape)
                return Spectrum(0.0f);

            // compute emitter interaction
            SurfaceInteraction3f si_testFinal =
                m_scene->ray_intersect(ray_final, active);
            if (!si_testFinal.shape)
                return Spectrum(0.0f);

            // +-------------------+
            // |Specular throughput|
            // +-------------------+
            // fresnel P1
            Spectrum f1(0.0f);
            Float eta, cos_theta;
            Complex<Spectrum> ior =
                si_heightfield_1.bsdf()->ior(si_heightfield_1);
            Mask reflection = any(neq(0.f, imag(ior)));
            Frame3f frame =
                si_heightfield_1.bsdf()->frame(si_heightfield_1, 0.f);
            cos_theta = dot(frame.n, ray_2.d);
            eta = hmean(real(ior));

            eta = hmean(real(ior));
            if (cos_theta < 0.f) {
                eta = rcp(eta);
            }
            auto [F, unused_0, unused_1, unused_2] = fresnel(cos_theta, eta);
            f1 = 1.f - F;
            f1 *= sqr(eta);

            // fresnel P2
            Spectrum f2(0.0f);
            ior = si_heightfield_2.bsdf()->ior(si_heightfield_2);
            reflection = any(neq(0.f, imag(ior)));
            frame = si_heightfield_2.bsdf()->frame(si_heightfield_2, 0.f);
            cos_theta = dot(frame.n, ray_final.d);
            eta = hmean(real(ior));

            eta = hmean(real(ior));
            if (cos_theta < 0.f) {
                eta = rcp(eta);
            }
            auto [F2, unused_4, unused_5, unused_6] = fresnel(cos_theta, eta);
            f2 = 1.f - F2;
            f2 *= sqr(eta);

            // geometry factor from contraint derivative (see W. Jakob, manifold
            // exploration, 2012)
            Float geometry_factor = generalizedGeometryFactor(
                si_O, si_heightfield_1, si_heightfield_2, si_testFinal);

            specular_out = f1 * f2 * geometry_factor;
        }
    }
    return specular_out;
}

template <typename Float, typename Spectrum>
Spectrum FermatNEE<Float, Spectrum>::compute_ray_contribution_reflection(
    SurfaceInteraction3f si_O, const BSDFContext &ctx, RayDifferential3f ray,
    EmitterInteraction vy, Float n1, Float n2, ShapePtr H1, ShapePtr H2,
    Mask active) const {
    MTS_MASK_ARGUMENT(active);

    Point3f S = vy.p;
    Spectrum specular_out(0.0f);
    SurfaceInteraction3f si_heightfield_1 = m_scene->ray_intersect(ray, active);

    // +------------------------------+
    // |first heightfield intersection:|
    // +------------------------------+
    // OpenGL_viewer_client<Float, Spectrum>::OpenGL_command_drawline(ray.o,
    // ray.o + ray.d*0.1 , Vector3f(1,1,1));

    if (si_heightfield_1.shape &&
        si_heightfield_1.shape->is_caustic_caster_double_refraction() &&
        si_heightfield_1.shape == H1) {

        // std::cout << "success:" << std::endl;
        // std::cout << ray << std::endl;

        // OpenGL_viewer_client<Float, Spectrum>::OpenGL_command_drawline(ray.o,
        // ray.o + ray.d * si_heightfield_1.t, Vector3f(0,0,1)); specular_out =
        // Spectrum(1);

        auto [bs, bsdf_val] = si_heightfield_1.bsdf()->sample(
            ctx, si_heightfield_1, m_sampler->next_1d(), m_sampler->next_2d(),
            active);
        RayDifferential3f ray_2 =
            si_heightfield_1.spawn_ray(si_heightfield_1.to_world(bs.wo));

        Vector3f d_tmp = normalize(vy.p - si_heightfield_1.p);
        RayDifferential3f ray_final(si_heightfield_1.p +
                                        math::ShadowEpsilon<Float> * d_tmp,
                                    d_tmp, ray.time, ray.wavelengths);
        // +----------------------+
        // |check light visibility|
        // +----------------------+

        SurfaceInteraction3f si_testFinal =
            m_scene->ray_intersect(ray_final, active);
        if (vy.is_directional()) {
            if (si_testFinal.shape) // something in between envmap and X2
                return specular_out * 0;
        } else {
            if (!si_testFinal.shape) // no emitter is visible
                return specular_out * 0;
        }

        // check if we've hit the point we intended to reach (e.g: is the
        // sampled point on the backface of the emitter?)
        ray_final = Ray3f(
            si_heightfield_1.p, normalize(vy.p - si_heightfield_1.p),
            math::RayEpsilon<Float> * (1.f + hmax(abs(si_heightfield_1.p))),
            norm(vy.p - si_heightfield_1.p) * (1.f - math::RayEpsilon<Float>),
            si_heightfield_1.time, si_heightfield_1.wavelengths);

        if (m_scene->ray_test(ray_final))
            return 0.f;

        // OpenGL_viewer_client<Float,
        // Spectrum>::OpenGL_command_drawline(ray_final.o, ray_final.o +
        // ray_final.d * si_testFinal.t, Vector3f(1,0,0));

        // // +-------------------+
        // // |Specular throughput|
        // // +-------------------+
        // fresnel P1
        Spectrum f1(0.0f);
        Float eta, cos_theta;
        Complex<Spectrum> ior = si_heightfield_1.bsdf()->ior(si_heightfield_1);
        Mask reflection = any(neq(0.f, imag(ior)));
        Frame3f frame = si_heightfield_1.bsdf()->frame(si_heightfield_1, 0.f);
        cos_theta = dot(frame.n, ray_2.d);
        eta = hmean(real(ior));

        auto [F_, cos_theta_t, eta_it, eta_ti] =
            fresnel(Spectrum(abs(cos_theta)), real(ior));
        f1 = F_;

        // // geometry factor from contraint derivative (see W. Jakob, manifold
        // exploration, 2012)
        Float geometry_factor = generalizedGeometryFactor(
            si_O, si_heightfield_1, si_heightfield_1, si_testFinal, true);

        specular_out = f1 * geometry_factor;
    }
    // if (si_heightfield_1.shape){
    //     OpenGL_viewer_client<Float, Spectrum>::OpenGL_command_drawline(ray.o,
    //     ray.o + ray.d * si_heightfield_1.t, Vector3f(0,1,1));

    // }
    // else{
    //     std::cout << "failure:" << std::endl;
    //     OpenGL_viewer_client<Float, Spectrum>::OpenGL_command_drawline(ray.o,
    //     ray.o + ray.d * 0.25, Vector3f(1,0,0));
    //     // std::cout << si_heightfield_1 << std::endl;
    //     // std::cout << H1 << std::endl;
    //     std::cout << ray << std::endl;
    // }
    // else{
    //     std::cout << ray << std::endl;
    //     std::cout << si_heightfield_1 << std::endl;

    // }
    return specular_out;
}

// +===============================================================================================+
// |                                       FIND RAY PROPOSAL |
// +===============================================================================================+

template <typename Float, typename Spectrum>
std::tuple<bool, typename FermatNEE<Float, Spectrum>::Vector3f>
FermatNEE<Float, Spectrum>::specular_connection(SurfaceInteraction3f &si,
                                                const EmitterInteraction &ei,
                                                L_data_float *data,
                                                Vector4f init) {
    if (m_use_SMS) {
        return SMS_connection(si, ei, data, init);
    } else {
        return fermat_connection(si, ei, data, init);
    }
}

template <typename Float, typename Spectrum>
std::tuple<bool, typename FermatNEE<Float, Spectrum>::Vector3f>
FermatNEE<Float, Spectrum>::fermat_connection(SurfaceInteraction3f &si,
                                              const EmitterInteraction &ei,
                                              L_data_float *data,
                                              Vector4f init) const {
    ScopedPhase scope_phase(ProfilerPhase::FermaNEE_principal);

    // make data
    // initial value
    Array<Float, 4> X(0);
    if (init == Vector4f(-1)) {
        X = Array<Float, 4>(m_sampler->next_1d(), m_sampler->next_1d(),
                            m_sampler->next_1d(), m_sampler->next_1d());
    } else {
        X = init;
    }
    Float f_out = 0.0f;

    bool success = newton_solver_double(&X, data, &f_out);

    if (!success)
        return { false, Vector3f(0) };

    Point3f P1(f_out, X.x() * data->L1, X.y() * data->W1);

    Vector3f result = data->pH1->get_world_transform().transform_affine(
        normalize(P1 - data->O));

    return { true, result };
}

template <typename Float, typename Spectrum>
Float FermatNEE<Float, Spectrum>::eval_invPDF(L_data_float *data,
                                              Vector3f &proposal) const {
    ScopedPhase scope_phase(ProfilerPhase::FermaNEE_principal);

    Float inv_pdf = 0.f;
    while (1) {
        inv_pdf += 1.f;

        Array<Float, 4> X;
        X = Array<Float, 4>(m_sampler->next_1d(), m_sampler->next_1d(),
                            m_sampler->next_1d(), m_sampler->next_1d());

        Float f_out = 0.0f;
        bool success = newton_solver_double(&X, data, &f_out);
        if (success) {
            Point3f P1(f_out, X.x() * data->L1, X.y() * data->W1);
            Vector3f bernoulli_trial =
                data->pH1->get_world_transform().transform_affine(
                    normalize(P1 - data->O));
            Float diff = abs(dot(proposal, bernoulli_trial) - 1.f);

            if (diff < m_config.uniqueness_threshold)
                return inv_pdf;
        }
        if (inv_pdf > m_config.max_trials)
            return 0.f;
    }
    return inv_pdf;
}

template <typename Float, typename Spectrum>
std::tuple<bool, typename FermatNEE<Float, Spectrum>::Vector3f>
FermatNEE<Float, Spectrum>::SMS_connection(SurfaceInteraction3f &si,
                                           const EmitterInteraction &ei,
                                           L_data_float *data, Vector4f init) {

    ScopedPhase scope_phase(ProfilerPhase::FermaNEE_principal);

    /// specular path initialisation
    m_current_path.clear();
    m_seed_path.clear();
    m_offset_normals.clear();

    Point2f uv_1(0);
    Point2f uv_2(0);
    // random init
    if(init == Vector4f(-1)){
        uv_1= Point2f(m_sampler->next_1d(), m_sampler->next_1d());
        uv_2= Point2f(m_sampler->next_1d(), m_sampler->next_1d());
    }else{
        uv_1=Point2f(init.x(), init.y());
        uv_2=Point2f(init.z(), init.w());

    }
    ManifoldVertex x1 = manifoldVertex_from_uv(uv_1, data->pH1, si);
    ManifoldVertex x2 = manifoldVertex_from_uv(uv_2, data->pH2, si);

    // point X1
    m_current_path.push_back(x1);
    m_seed_path.push_back(x1);

    // point X2
    m_current_path.push_back(x2);
    m_seed_path.push_back(x2);

    // null normal offset
    Vector3f n_offset = Vector3f(0.f, 0.f, 1.f);
    m_offset_normals.push_back(n_offset);
    m_offset_normals.push_back(n_offset);

    /// newton_solver
    bool success = newton_solver_SMS(si, ei, data);

    if (!success)
        return { false, Vector3f(0) };

    /// output stage
    Point3f P1 = m_current_path[0].p;
    Vector3f result = normalize(P1 - si.p); // already in world coord
    return { true, result };
}

template <typename Float, typename Spectrum>
Float FermatNEE<Float, Spectrum>::SMS_eval_invPDF(SurfaceInteraction3f &si,
                                                  const EmitterInteraction &ei,
                                                  L_data_float *data,
                                                  Vector3f &proposal) {
    ScopedPhase scope_phase(ProfilerPhase::FermaNEE_principal);

    Float inv_pdf = 0.f;
    while (1) {
        inv_pdf += 1.f;
        auto [success, result] = SMS_connection(si, ei, data);
        if (success) {
            float diff = abs(dot(proposal, result) - 1.f);
            if (diff < m_config.uniqueness_threshold)
                return inv_pdf;
        }
        if (inv_pdf > m_config.max_trials)
            return 0.f;
    }
    return inv_pdf;
}

template <typename Float, typename Spectrum>
std::pair<bool, typename FermatNEE<Float, Spectrum>::Vector3f>
FermatNEE<Float, Spectrum>::specular_connection_debug(
    SurfaceInteraction3f &si, const EmitterInteraction &vy, Vector4f init) {
        
    ShapePtr H1;
    ShapePtr H2;
    Float e = 0;
    Float sample_heightfield_invpdf = 0.0f;
    bool sample_heighfield_success =
        sample_heightfield_pair(&H1, &H2, &sample_heightfield_invpdf, true);

    Point3f origin_H1 =
        H1->get_world_transform().transform_affine(Point3f(0.f));
    Point3f origin_H2 =
        H2->get_world_transform().transform_affine(Point3f(0.f));
    Point3f origin_H2_local_H1 =
        H1->get_local_transform().transform_affine(origin_H2);
    Point3f O_local_H1 = H1->get_local_transform().transform_affine(si.p);
    e = origin_H2_local_H1.x();

    if (squared_norm(O_local_H1) >
        squared_norm(O_local_H1 - origin_H2_local_H1)) {
        std::swap(H1, H2);
        e = -e;
    }
    if (!sample_heighfield_success)
        return { false, Vector3f(0) };

    Float NEE_invpdf = 1.0f;
    Point3f S = H1->get_local_transform().transform_affine(vy.p);
    Point3f O = H1->get_local_transform().transform_affine(si.p);

    Vector3f wo_r = normalize(vy.p - si.p);
    RayDifferential3f r(si.p, wo_r,
                        math::RayEpsilon<Float> * (1.f + hmax(abs(si.p))),
                        norm(vy.p - si.p) * (1.f - math::ShadowEpsilon<Float>),
                        si.time, si.wavelengths);

    // find ray proposal
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
        e
    };

    return specular_connection(si, vy, &fermat_connection_data, init);
}

template <typename Float, typename Spectrum>
bool FermatNEE<Float, Spectrum>::fermat_connection_reflection(
    L_data_float *data, Vector3f *result, Float *out, Point2f init,
    Matrix<double, 2> *H_out) const {
    // +-----------------------+
    // |Setup solver input data|
    // +-----------------------+
    ScopedPhase scope_phase(ProfilerPhase::FermaNEE_principal);
    SurfaceInteraction3f si;

    // make data
    // initial value
    // Float x1 = m_sampler->next_1d();
    // Float x2 = m_sampler->next_1d();
    // Array<Float, 4> X(x1,x2,x1,x2);
    Vector2f X(m_sampler->next_1d(), m_sampler->next_1d());
    if (init != Point2f(-1, -1)) {
        X = Array<Float, 2>(init.x(), init.y());
        // std::cout << "init" << X << std::endl;
    }
    // +--------------------+
    // |proposal computation|
    // +--------------------+
    Float f_out = 0.0f;

    bool success = newton_solver_reflection(&X, data, &f_out, nullptr, H_out);
    if (!success) {
        // std::cout<<"FAILLURE: ";
        // std::cout<<X<<std::endl;
        return false;
    }

    // std::cout<<"SUCCESS: ";
    // std::cout<<X<<std::endl;
    // +------+
    // |Output|
    // +------+
    if (X.x() < 1 && X.x() > 0 && X.y() > 0 && X.y() < 1) {
        Point3f P1(f_out, X.x() * data->L1, X.y() * data->W1);
        *result = normalize(P1 - data->O);
        return true;
    }
    return false;
}

template <typename Float, typename Spectrum>
Float FermatNEE<Float, Spectrum>::eval_invPDF_reflection(
    L_data_float *data, Vector3f &proposal) const {
    // +-----------------------+
    // |Setup solver input data|
    // +-----------------------+
    ScopedPhase scope_phase(ProfilerPhase::FermaNEE_principal);
    SurfaceInteraction3f si;

    // FIXME parameters
    Float solutionIdentical_threshold = 0.000001f;
    int maxBernouilliTrial = 1000;
    /* +----------------------+
    // |inverse pdf evaluation|
    // +----------------------+ */

    Float inv_pdf = 0.f;
    while (1) {
        inv_pdf += 1.f;
        // Float x1 = m_sampler->next_1d();
        // Float x2 = m_sampler->next_1d();
        // Array<Float, 4> X(x1,x2,x1,x2);
        Vector2f X(m_sampler->next_1d(), m_sampler->next_1d());
        // Array<Float, 4> X(m_sampler->next_1d(), m_sampler->next_1d(),
        // m_sampler->next_1d(), m_sampler->next_1d());

        Float f_out = 0.0f;
        Float g_out = 0.0f;
        bool success = newton_solver_reflection(&X, data, &f_out);
        if (success) {
            Point3f P1(f_out, X.x() * data->L1, X.y() * data->W1);
            Vector3f bernoulli_trial = normalize(P1 - data->O);
            Float diff = abs(dot(proposal, bernoulli_trial) - 1.f);

            if (inv_pdf > maxBernouilliTrial) {
                inv_pdf = 0.0f;
                break;
            }
            if (diff < solutionIdentical_threshold)
                return inv_pdf;
        }
    }

    return inv_pdf;
}

template <typename Float, typename Spectrum>
typename FermatNEE<Float, Spectrum>::Ray3f
FermatNEE<Float, Spectrum>::straighline_approx(ShapePtr H1, ShapePtr H2,
                                               Point3f O, Point3f S, Float n1,
                                               Float n2, const Float *invpdf,
                                               bool bComputePDF) const {
    Ray3f r = Ray3f();
    r.o = O;
    r.d = normalize(S - O);

    return r;
}

template <typename Float, typename Spectrum>
typename FermatNEE<Float, Spectrum>::Mask
FermatNEE<Float, Spectrum>::reproject_raytrace(
    const SurfaceInteraction3f &si_) {
    m_proposed_path.clear();
    SurfaceInteraction3f si(si_);

    // Start ray-tracing towards the first specular vertex along the chain
    Point3f x0 = si.p, x1 = m_proposed_positions[0];
    Vector3f wo = normalize(x1 - x0);
    Ray3f ray(x0, wo, si.time, si.wavelengths);

    // std::cout << "Reproject raytrace :" << std::endl;

    while (true) {
        int bounce = m_proposed_path.size();

        if (bounce >= m_config.bounces) {
            /* We reached the number of specular bounces that was requested.
               (Implicitly) connect to the light source now by terminating. */
            break;
        }

        si = m_scene->ray_intersect(ray);
        if (!si.is_valid()) {
            return false;
        }
        const ShapePtr shape = si.shape;
        if (shape != m_seed_path[bounce].shape) {
            // We hit some other shape than previously
            return false;
        }

        // Create the path vertex
        ManifoldVertex vertex = ManifoldVertex(si, 0.f);
        m_proposed_path.push_back(vertex);

        // std::cout << vertex.shape << std::endl;

        // Get current (potentially offset) normal in world space
        Vector3f n_offset = m_offset_normals[bounce];
        // Vector3f n_offset(0);
        Vector3f m = vertex.s * n_offset[0] + vertex.t * n_offset[1] +
                     vertex.n * n_offset[2];

        // Perform scattering at vertex
        Vector3f wi = -wo;
        bool scatter_success = false;
        if (vertex.eta == 1.f) {
            std::tie(scatter_success, wo) = SpecularManifold::reflect(wi, m);
        } else {
            std::tie(scatter_success, wo) =
                SpecularManifold::refract(wi, m, vertex.eta);
        }

        if (!scatter_success) {
            // We must have hit total internal reflection. Abort.
            return false;
        }

        ray = si.spawn_ray(wo);
    }

    return m_proposed_path.size() == m_seed_path.size();
}

// FIXME si_ useless ?
template <typename Float, typename Spectrum>
typename FermatNEE<Float, Spectrum>::Mask FermatNEE<Float, Spectrum>::reproject(
    const SurfaceInteraction3f &si_,
    L_data_float *data) { // m_proposed_path.clear();
    SurfaceInteraction3f si(si_);
    // ManifoldVertex vy = m_proposed_path[2];

    // point X1
    Point3f X1_local = data->pH1->get_local_transform().transform_affine(
        m_proposed_positions[0]);
    Point2f uv_1(X1_local.y() * rcp(data->L1), X1_local.z() * rcp(data->W1));
    ManifoldVertex x1 = manifoldVertex_from_uv(uv_1, data->pH1, si);

    // point X2
    Point3f X2_local = data->pH2->get_local_transform().transform_affine(
        m_proposed_positions[1]);
    Point2f uv_2(X2_local.y() * rcp(data->L2), X2_local.z() * rcp(data->W2));
    ManifoldVertex x2 = manifoldVertex_from_uv(uv_2, data->pH2, si);

    if (uv_1.x() > data->L1 || uv_2.x() > data->L1 || uv_1.x() < 0 ||
        uv_2.x() < 0 || uv_1.y() > data->L1 || uv_2.y() > data->L1 ||
        uv_1.y() < 0 || uv_2.y() < 0) {
        return false;
    }

    m_proposed_path.clear();
    m_proposed_path.push_back(x1);
    m_proposed_path.push_back(x2);

    return true;
}

template <typename Float, typename Spectrum>
bool FermatNEE<Float, Spectrum>::sample_heightfield_pair(ShapePtr *H1,
                                                         ShapePtr *H2,
                                                         Float *invPdf,
                                                         Mask active) const {
    MTS_MASK_ARGUMENT(active);

    auto shapes = m_scene->caustic_casters_double_refraction();
    uint32_t size = shapes.size();
    *invPdf = ((ScalarFloat) (size)) * 0.5f; // (* 0.5) as we are sampling pairs

    if (unlikely(shapes.size() == 0)) {
        return 0.f;
    }
    // select a random pair:
    UInt32 index =
        min(UInt32(m_sampler->next_1d() * (ScalarFloat) shapes.size()),
            (uint32_t) shapes.size() - 1);
    *H1 = shapes[index];

    // TODO direct reference
    for (size_t i = 0; i < shapes.size(); i++) {
        if (shapes[i]->get_self_id(active) == (*H1)->get_pair_id(active)) {
            *H2 = shapes[i];
            *invPdf = shapes.size() * 0.5f;
            return true;
        }
    }
    return false;
}

template <typename Float, typename Spectrum>
bool FermatNEE<Float, Spectrum>::sample_heightfield(ShapePtr *H1, Float *invPdf,
                                                    Mask active) const {
    MTS_MASK_ARGUMENT(active);

    auto shapes = m_scene->caustic_casters_double_refraction();
    uint32_t size = shapes.size();
    if (unlikely(size == 0)) {
        return false;
    }
    // select a random pair:
    UInt32 index = min(UInt32(m_sampler->next_1d() * (ScalarFloat) size),
                       (uint32_t) size - 1);
    *H1 = shapes[index];
    *invPdf = ((ScalarFloat) (size)) * 0.5f; // (* 0.5) as we are sampling pairs

    return true;
}

template <typename Float, typename Spectrum>
typename FermatNEE<Float, Spectrum>::ManifoldVertex
FermatNEE<Float, Spectrum>::manifoldVertex_from_uv(Point2f uv,
                                                   ShapePtr heightfield,
                                                   SurfaceInteraction3f &si_) {

    SurfaceInteraction3f si =
        heightfield->eval_surfaceInteraction_from_uv(uv, true, si_);
    si.shape = heightfield;
    ManifoldVertex x0(si, 0.0f);
    x0.make_orthonormal();
    // std::cout<<"si"<<std::endl;
    return x0;
}

MTS_INSTANTIATE_CLASS(FermatNEE)
NAMESPACE_END(mitsuba)
