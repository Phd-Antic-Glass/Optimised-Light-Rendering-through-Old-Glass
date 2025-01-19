#include <iomanip>
#include <mitsuba/core/fwd.h>
#include <mitsuba/core/math.h>
#include <mitsuba/core/properties.h>
#include <mitsuba/core/spectrum.h>
#include <mitsuba/core/transform.h>
#include <mitsuba/core/util.h>
#include <mitsuba/core/warp.h>
#include <mitsuba/render/bsdf.h>
#include <mitsuba/render/fwd.h>
#include <mitsuba/render/interaction.h>
#include <mitsuba/render/shape.h>
#include <mitsuba/render/texture.h>
#include <sys/select.h>

// #include <optim.hpp>

#if defined(MTS_ENABLE_OPTIX)
#include "optix/heightfield.cuh"
#include <mitsuba/render/optix_api.h>
#endif

NAMESPACE_BEGIN(mitsuba)

/**!

.. _shape-heightfield:

Heightfield (:monosp:`heightfield`)
-------------------------------------------------

.. pluginparameters::

 * - origin
   - |point|
   - Origin of the heightfield (Default: (0, 0, 0))
* - L
   - |float|
   - length (|) of the heightfield (Default: 0)
* - W
   - |float|
   - width (---) of the heightfield (Default: 0)
 * - H
   - |float|
   - maximum elevation of the heightfield (Default: 0). Extends in normals
direction
 * - flip_normals
   - |bool|
   - Is the heightfield inverted, i.e. should the normal vectors be flipped?
(Default:|false|, i.e. the normals point outside)
 * - to_world
   - |transform|
   -  Specifies an optional linear object-to-world transformation.
      Note that non-uniform scales and shears are not permitted!
      (Default: none, i.e. object space = world space)

This shape plugin describes a heightfield intersection primitive: given an .exr
elevation profile, a surface is reconstructed using bicubic interpolation

A heightfield can either be configured using a linear :monosp:`to_world`
transformation or the :monosp:`center` and :monosp:`radius` parameters (or
both). The two declarations below are equivalent.


.. warning:: This plugin is currently CPU only.

 */

template <typename Float, typename Spectrum>
class MTS_EXPORT_RENDER Heightfield final : public Shape<Float, Spectrum> {
public:
    MTS_IMPORT_BASE(Shape, m_to_world, m_to_object, set_children,
                    get_children_string)
    MTS_IMPORT_TYPES(Texture)

    using typename Base::ScalarSize;

    void generate_bicubic_coefs(ref<Texture> heightmap,
                                Array<Float, 16> **coefs) {

        Float scale = rcp((Float) (N));

        // precompute bicubic interpolator coefficients:
        // interpolator coefficients
        *coefs = new Array<Float, 16>[N * N];
        // for each pixels:
        for (size_t i = 0; i < N; i++) {
            for (size_t j = 0; j < N; j++) {

                // account for border cells
                size_t X = i, Y = j;

                X = min(N - 3, max(1, X));
                Y = min(N - 3, max(1, Y));

                // Fill FXY: get the 16 reference texture points
                SurfaceInteraction3f si;
                si.uv = Point2f(X - 1, Y - 1) * scale;
                Float F00 = heightmap->eval_1(si, "bilinear");
                si.uv = Point2f(X - 1, Y + 0) * scale;
                Float F01 = heightmap->eval_1(si, "bilinear");
                si.uv = Point2f(X - 1, Y + 1) * scale;
                Float F02 = heightmap->eval_1(si, "bilinear");
                si.uv = Point2f(X - 1, Y + 2) * scale;
                Float F03 = heightmap->eval_1(si, "bilinear");

                si.uv = Point2f(X + 0, Y - 1) * scale;
                Float F10 = heightmap->eval_1(si, "bilinear");
                si.uv = Point2f(X + 0, Y + 0) * scale;
                Float F11 = heightmap->eval_1(si, "bilinear");
                si.uv = Point2f(X + 0, Y + 1) * scale;
                Float F12 = heightmap->eval_1(si, "bilinear");
                si.uv = Point2f(X + 0, Y + 2) * scale;
                Float F13 = heightmap->eval_1(si, "bilinear");
                si.uv = Point2f(X + 1, Y - 1) * scale;

                Float F20 = heightmap->eval_1(si, "bilinear");
                si.uv = Point2f(X + 1, Y + 0) * scale;
                Float F21 = heightmap->eval_1(si, "bilinear");
                si.uv = Point2f(X + 1, Y + 1) * scale;
                Float F22 = heightmap->eval_1(si, "bilinear");
                si.uv = Point2f(X + 1, Y + 2) * scale;
                Float F23 = heightmap->eval_1(si, "bilinear");

                si.uv = Point2f(X + 2, Y - 1) * scale;
                Float F30 = heightmap->eval_1(si, "bilinear");
                si.uv = Point2f(X + 2, Y + 0) * scale;
                Float F31 = heightmap->eval_1(si, "bilinear");
                si.uv = Point2f(X + 2, Y + 1) * scale;
                Float F32 = heightmap->eval_1(si, "bilinear");
                si.uv = Point2f(X + 2, Y + 2) * scale;
                Float F33 = heightmap->eval_1(si, "bilinear");

                // C1
                Array<Float, 16> coef(
                    1.0F * F00,
                    -1.83333333F * F00 + 3.0F * F01 - 1.5F * F02 +
                        0.333333333F * F03,
                    1.0F * F00 - 2.5F * F01 + 2.0F * F02 - 0.5F * F03,
                    -0.166666667F * F00 + 0.5F * F01 - 0.5F * F02 +
                        0.166666667F * F03,
                    -1.83333333F * F00 + 3.0F * F10 - 1.5F * F20 +
                        0.333333333F * F30,
                    3.36111111F * F00 - 5.5F * F01 + 2.75F * F02 -
                        0.611111111F * F03 - 5.5F * F10 + 9.0F * F11 -
                        4.5F * F12 + 1.0F * F13 + 2.75F * F20 - 4.5F * F21 +
                        2.25F * F22 - 0.5F * F23 - 0.611111111F * F30 +
                        1.0F * F31 - 0.5F * F32 + 0.111111111F * F33,
                    -1.83333333F * F00 + 4.58333333F * F01 - 3.66666667F * F02 +
                        0.916666667F * F03 + 3.0F * F10 - 7.5F * F11 +
                        6.0F * F12 - 1.5F * F13 - 1.5F * F20 + 3.75F * F21 -
                        3.0F * F22 + 0.75F * F23 + 0.333333333F * F30 -
                        0.833333333F * F31 + 0.666666667F * F32 -
                        0.166666667F * F33,
                    0.305555556F * F00 - 0.916666667F * F01 +
                        0.916666667F * F02 - 0.305555556F * F03 - 0.5F * F10 +
                        1.5F * F11 - 1.5F * F12 + 0.5F * F13 + 0.25F * F20 -
                        0.75F * F21 + 0.75F * F22 - 0.25F * F23 -
                        0.0555555556F * F30 + 0.166666667F * F31 -
                        0.166666667F * F32 + 0.0555555556F * F33,
                    1.0F * F00 - 2.5F * F10 + 2.0F * F20 - 0.5F * F30,
                    -1.83333333F * F00 + 3.0F * F01 - 1.5F * F02 +
                        0.333333333F * F03 + 4.58333333F * F10 - 7.5F * F11 +
                        3.75F * F12 - 0.833333333F * F13 - 3.66666667F * F20 +
                        6.0F * F21 - 3.0F * F22 + 0.666666667F * F23 +
                        0.916666667F * F30 - 1.5F * F31 + 0.75F * F32 -
                        0.166666667F * F33,
                    1.0F * F00 - 2.5F * F01 + 2.0F * F02 - 0.5F * F03 -
                        2.5F * F10 + 6.25F * F11 - 5.0F * F12 + 1.25F * F13 +
                        2.0F * F20 - 5.0F * F21 + 4.0F * F22 - 1.0F * F23 -
                        0.5F * F30 + 1.25F * F31 - 1.0F * F32 + 0.25F * F33,
                    -0.166666667F * F00 + 0.5F * F01 - 0.5F * F02 +
                        0.166666667F * F03 + 0.416666667F * F10 - 1.25F * F11 +
                        1.25F * F12 - 0.416666667F * F13 - 0.333333333F * F20 +
                        1.0F * F21 - 1.0F * F22 + 0.333333333F * F23 +
                        0.0833333333F * F30 - 0.25F * F31 + 0.25F * F32 -
                        0.0833333333F * F33,
                    -0.166666667F * F00 + 0.5F * F10 - 0.5F * F20 +
                        0.166666667F * F30,
                    0.305555556F * F00 - 0.5F * F01 + 0.25F * F02 -
                        0.0555555556F * F03 - 0.916666667F * F10 + 1.5F * F11 -
                        0.75F * F12 + 0.166666667F * F13 + 0.916666667F * F20 -
                        1.5F * F21 + 0.75F * F22 - 0.166666667F * F23 -
                        0.305555556F * F30 + 0.5F * F31 - 0.25F * F32 +
                        0.0555555556F * F33,
                    -0.166666667F * F00 + 0.416666667F * F01 -
                        0.333333333F * F02 + 0.0833333333F * F03 + 0.5F * F10 -
                        1.25F * F11 + 1.0F * F12 - 0.25F * F13 - 0.5F * F20 +
                        1.25F * F21 - 1.0F * F22 + 0.25F * F23 +
                        0.166666667F * F30 - 0.416666667F * F31 +
                        0.333333333F * F32 - 0.0833333333F * F33,
                    0.0277777778F * F00 - 0.0833333333F * F01 +
                        0.0833333333F * F02 - 0.0277777778F * F03 -
                        0.0833333333F * F10 + 0.25F * F11 - 0.25F * F12 +
                        0.0833333333F * F13 + 0.0833333333F * F20 -
                        0.25F * F21 + 0.25F * F22 - 0.0833333333F * F23 -
                        0.0277777778F * F30 + 0.0833333333F * F31 -
                        0.0833333333F * F32 + 0.0277777778F * F33);
                // store everything:
                for (size_t k = 0; k < 16; k++) {
                    (*(*coefs + i + N * j))[k] = coef[k];
                }
            }
        }
        return;
    }

    Heightfield(const Properties &props) : Base(props) {
        /// Are the heightfield normals pointing inwards? default: no
        m_flip_normals = props.bool_("flip_normals", false);
        flip = select(!m_flip_normals, 1.f, -1.f);
        // dimensions
        L = props.float_("L", 1);
        W = props.float_("W", 1);
        H = props.float_("H", 1);

        // load texture
        heightmap = props.texture<Texture>("heightmap", .5f);
        N = heightmap->resolution().x();
        // precompute scaling factors
        invW = rcp(W);
        invL = rcp(L);
        dx_scale = invL * (Float) (N);
        dy_scale = invW * (Float) (N);
        dx_scale_sqr = dx_scale * dx_scale;
        dy_scale_sqr = dy_scale * dy_scale;
        du_scale = L * rcp(Float(N));
        dv_scale = W * rcp(Float(N));
        du_scale_sqr = du_scale * du_scale;
        // select surface reconstruction mode (bilinear or bicubic)
        reconstructionMode = props.string("reconstructionMode", "bilinear");

        Float scale = rcp((Float) (N));
        scale_w = scale; // same since heightmap is a N*N texture
        scale_l = scale; // same since heightmap is a N*N texture

        generate_bicubic_coefs(heightmap, &cubic_coefs);

        // define partition size (for multi-piece glass panels)
        partition_size = rcp(props.float_("partition_number", 1.0));
        inv_partition_size = rcp(partition_size);
        inv_partition_area = sqr(sqr(inv_partition_size));

        // reference to heightfield pair
        m_self_id = props.int_("self_id", 0);
        m_pair_id = props.int_("pair_id", 0);

        update();
        set_children();
    }

    void update() {
        m_to_object = m_to_world.inverse();

        auto [S, Q, T] = transform_decompose(m_to_world.matrix);

        if (abs(S[0][1]) > 1e-6f || abs(S[0][2]) > 1e-6f ||
            abs(S[1][0]) > 1e-6f || abs(S[1][2]) > 1e-6f ||
            abs(S[2][0]) > 1e-6f || abs(S[2][1]) > 1e-6f)
            Log(Warn, "'to_world' transform shouldn't contain any shearing!");

        m_center = T;

        ScalarVector3f dp_du = m_to_world * ScalarVector3f(0.f, 2.f, 0.f);
        ScalarVector3f dp_dv = m_to_world * ScalarVector3f(0.f, 0.f, 2.f);
        ScalarNormal3f normal =
            flip * normalize(m_to_world * ScalarNormal3f(-1.f, 0.f, 0.f));

        m_frame = ScalarFrame3f(dp_du, dp_dv, normal);

        m_inv_surface_area = rcp(surface_area());
    }

    ScalarBoundingBox3f bbox() const override {
        ScalarBoundingBox3f bbox;
        bbox.expand(m_to_world.transform_affine(ScalarPoint3f(0, L, 0)));
        bbox.expand(m_to_world.transform_affine(ScalarPoint3f(0, L, W)));
        bbox.expand(m_to_world.transform_affine(ScalarPoint3f(0, 0, W)));
        bbox.expand(m_to_world.transform_affine(ScalarPoint3f(0, 0, 0)));
        bbox.expand(m_to_world.transform_affine(ScalarPoint3f(-H, L, 0)));
        bbox.expand(m_to_world.transform_affine(ScalarPoint3f(-H, L, W)));
        bbox.expand(m_to_world.transform_affine(ScalarPoint3f(-H, 0, W)));
        bbox.expand(m_to_world.transform_affine(ScalarPoint3f(-H, 0, 0)));

        return bbox;
    }

    // since H is generally very small, surface_area() could be aproximated as
    // L*W
    ScalarFloat surface_area() const override {
        Float S = L * W;
        return S;
    }

    // =============================================================
    //! @{ \name Sampling routines
    // =============================================================

    PositionSample3f sample_position(Float time, const Point2f &sample,
                                     Mask active) const override {
        MTS_MASK_ARGUMENT(active);

        PositionSample3f ps;
        Float S_x = h_f(sample, active);

        ps.p = m_to_world.transform_affine(
            Point3f(S_x, sample.x() * L,
                    sample.y() * W)); // H being quite small, this is considered
                                      // uniform sampling
        ps.n = (m_to_world * normalize(flip * grad_h(sample, active)));
        ps.pdf = m_inv_surface_area;
        ps.uv = sample;
        ps.time = time;
        ps.delta = false;

        return ps;
    }

    Float pdf_position(const PositionSample3f & /*ps*/,
                       Mask active) const override {
        MTS_MASK_ARGUMENT(active);
        return m_inv_surface_area;
    }

    // =============================================================
    //! @{ \name attribute get (messy...)
    // =============================================================
    Float eval_attribute_1(const std::string &name,
                           const SurfaceInteraction3f &si,
                           Mask active = true) const override {
        MTS_MASK_ARGUMENT(active);
        if (name == "H")
            return H;
        else if (name == "L")
            return L;
        else if (name == "W")
            return W;
        else if (name == "m_flip_normals")
            return m_flip_normals ? -1.0f : 1.0f;
        else if (name == "partition_size")
            return partition_size;
        else if (name == "inv_partition_size")
            return inv_partition_size;
        else if (name == "inv_partition_area")
            return inv_partition_area;
        else if (name == "N")
            return Float(N);
        else if (name == "nbEval")
            return Float(nbEval);
        else
            return 0.0f;
    }

    Array<Float, 6> eval_curvature(Point2f p, Mask active, bool full) const {

        // eval point is outside of heightfield definition domain
        if (unlikely(p.x() > 1 || p.y() > 1 || p.x() < 0 || p.y() < 0)) {
            return 0.0f;
        } else {

            Float x = p.x();
            Float y = p.y();
            Float x_ = x * Float(N);
            Float y_ = y * Float(N);
            UInt32 X = (floor(x_));
            UInt32 Y = (floor(y_));
            Float dx = (x_ - X);
            Float dy = (y_ - Y);

            Float d2y = dy * dy, d3y = d2y * dy, d2x = dx * dx, d3x = d2x * dx;

            X = min(max(X, 1), N - 3);
            Y = min(max(Y, 1), N - 3);

            Array<Float, 4> vY(1.f, dy, d2y, d3y);
            Array<Float, 4> vX(1.f, dx, d2x, d3x);

            Array<Float, 16> coefs(0);

            coefs = load<Array<Float, 16>>(cubic_coefs + X + N * Y);
            Matrix<Float, 4> A(coefs[0], coefs[1], coefs[2], coefs[3], coefs[4],
                               coefs[5], coefs[6], coefs[7], coefs[8], coefs[9],
                               coefs[10], coefs[11], coefs[12], coefs[13],
                               coefs[14], coefs[15]);

            Array<Float, 4> A_t_vY = A * vY;
            Float f = -H * dot(vX, A_t_vY);
            if (full) {
                Float du =
                    -H * dx_scale *
                    dot(Array<Float, 4>(0.f, 1.f, 2.f * dx, 3.f * d2x), A_t_vY);
                Float dv =
                    -H * dx_scale *
                    dot(vX, A * Array<Float, 4>(0.f, 1.f, 2.f * dy, 3.f * d2y));
                Float duu =
                    -H * dx_scale_sqr *
                    dot(Array<Float, 4>(0.f, 0.f, 2.f, 6.f * dx), A_t_vY);
                Float dvv =
                    -H * dx_scale_sqr *
                    dot(vX, A * Array<Float, 4>(0.f, 0.f, 2.f, 6.f * dy));
                Float duv =
                    -H * dx_scale_sqr *
                    dot(Array<Float, 4>(0.f, 1.f, 2.f * dx, 3.f * d2x),
                        A * Array<Float, 4>(0.f, 1.f, 2.f * dy, 3.f * d2y));

                return Array<Float, 6>(f, du, dv, duu, dvv, duv);
            } else {
                return Array<Float, 6>(f, 0, 0, 0, 0, 0);
            }
        }
    }

    // =============================================================
    //! @{ \name Window heightfield pair reference
    // =============================================================
    int get_pair_id(Mask /*active*/) const override { return m_pair_id; }

    int get_self_id(Mask /*active*/) const override { return m_self_id; }

    Float eval_texture(int index, Point2f p,
                       Mask active = true) const override {
        MTS_MASK_ARGUMENT(active);


        // nbEval++;
        if (index == 0) {
            // return h_f(p, active);
            return -H * bicubic(p, active, index);

        } else if (index == 1) {
            // return duf(p, active);
            return -H * dx_scale * bicubic(p, active, index);

        } else if (index == 2) {
            // return dvf(p, active);1
            return -H * dy_scale * bicubic(p, active, index);

        } else if (index == 3) {
            // return duuf(p, active);
            return -H * dx_scale * dx_scale * bicubic(p, active, index);

        } else if (index == 4) {
            // return dvvf(p, active);
            return -H * dy_scale * dy_scale * bicubic(p, active, index);

        } else if (index == 5) {
            // return dvuf(p, active);
            return -H * dx_scale * dy_scale * bicubic(p, active, index);
        } else {
            std::cout << "error, Heighfield eval_texture() <name> " << index
                      << "is not a valid parameter" << std::endl;
        }
        return 0.0f;
    }

    // =============================================================
    //! @{ \name Ray tracing routines
    // =============================================================

    // Heighfield function description
    Float h_f(Point2f p, Mask active) const {
        MTS_MASK_ARGUMENT(active);

        return select(H != 0, -H * bicubic(p, active), 0.f);
    }

    // bicubic interpolator
    Float bicubic(Point2f p, Mask active, int partial = 0) const {
        MTS_MASK_ARGUMENT(active);
        ScopedPhase scope_phase(ProfilerPhase::Eval_heighfield);

        // eval point is outside of heightfield definition domain
        if (p.x() > 1 || p.y() > 1 || p.x() < 0 || p.y() < 0) {
            return 0.0f;
        }

        Float x = p.x();
        Float y = p.y();
        Float x_ = x * Float(N);
        Float y_ = y * Float(N);
        UInt32 X = (floor(x_));
        UInt32 Y = (floor(y_));
        Float dx = (x_ - X);
        Float dy = (y_ - Y);

        Float d2y = dy * dy, d3y = d2y * dy, d2x = dx * dx, d3x = d2x * dx;

        X = min(max(X, 1), N - 3);
        Y = min(max(Y, 1), N - 3);

        Array<Float, 4> vY(1.f, dy, d2y, d3y);
        Array<Float, 4> vX(1.f, dx, d2x, d3x);

        Array<Float, 16> coefs(0);

        coefs = load<Array<Float, 16>>(cubic_coefs + X + N * Y);

        if (partial == 0) {
            // nothing to do
        } else if (partial == 1) {
            vX = Array<Float, 4>(0.f, 1.f, 2.f * dx, 3.f * d2x);
        } else if (partial == 2) {
            vY = Array<Float, 4>(0.f, 1.f, 2.f * dy, 3.f * d2y);
        } else if (partial == 3) {
            vX = Array<Float, 4>(0.f, 0.f, 2.f, 6.f * dx);
        } else if (partial == 4) {
            vY = Array<Float, 4>(0.f, 0.f, 2.f, 6.f * dy);
        } else if (partial == 5) {
            vX = Array<Float, 4>(0.f, 1.f, 2.f * dx, 3.f * d2x);
            vY = Array<Float, 4>(0.f, 1.f, 2.f * dy, 3.f * d2y);
        } else {
            std::cout << "error, Heighfield bicubic() <partial> " << partial
                      << "is not a valid parameter" << std::endl;
            return 0;
        }

        Matrix<Float, 4> A(coefs[0], coefs[1], coefs[2], coefs[3], coefs[4],
                           coefs[5], coefs[6], coefs[7], coefs[8], coefs[9],
                           coefs[10], coefs[11], coefs[12], coefs[13],
                           coefs[14], coefs[15]);
        return dot(vX, A * vY);
    }

    // Heighfield distance field description
    Float h(Point3f P, Mask active) const {
        MTS_MASK_ARGUMENT(active);
        Float x = P.x();
        Point2f p(P.y() * invL, P.z() * invW);
        return h_f(p, active) - x;
    }

    // gradient computation
    Normal3f grad_h(Point2f P, Mask active) const {
        MTS_MASK_ARGUMENT(active);
        if (H != 0) {
            Float df_du = -H * bicubic(P, active, 1);
            Float df_dv = -H * bicubic(P, active, 2);
            Normal3f dp_du = Normal3f(df_du, du_scale, 0);
            Normal3f dp_dv = Normal3f(df_dv, 0, du_scale);

            // return normalize(Normal3f(du_scale_sqr, df_du*du_scale,
            // -df_dv*du_scale));
            return normalize(cross(dp_du, dp_dv));
        } else {
            return Normal3f(1.0f, 0, 0);
        }
    }

    // Intersect ray_ (local space) and Heightfield surface by sphere tracing
    Float sphereTrace(const Ray3f &ray_, Mask active) const {
        MTS_MASK_ARGUMENT(active);
        Float epsilon = 1e-6f;

        Float t = 0;
        Float d = h(ray_(t), active);

        // prevent unwanted backface culling
        Float minus = select(d >= 0, 1.0f, -1.0f);
        d = minus * d;

        // main loop
        int c = 0;
        while (d > epsilon) {
            t = t + d;
            d = minus * h(ray_(t), active);
            c++;

            if (t > ray_.maxt || t < 0 || c > 1000)
                return std::numeric_limits<ScalarFloat>::infinity();
        }
        return t;
    }
    /** \brief Ray-triangle intersection test
     *
     * Uses the algorithm by Moeller and Trumbore discussed at
     * <tt>http://www.acm.org/jgt/papers/MollerTrumbore97/code.html</tt>.
     *
     * \param index
     *    Index of the triangle to be intersected.
     * \param ray
     *    The ray segment to be used for the intersection query.
     * \return
     *    Returns an ordered tuple <tt>(mask, u, v, t)</tt>, where \c mask
     *    indicates whether an intersection was found, \c t contains the
     *    distance from the ray origin to the intersection point, and \c u and
     *    \c v contains the first two components of the intersection in
     *    barycentric coordinates
     */
    MTS_INLINE std::tuple<Mask, Float, Float, Float>
    ray_intersect_triangle(const Ray3f &ray, Point3f p0, Point3f p1, Point3f p2,
                           identity_t<Mask> active = true) const {

        Vector3f e1 = p1 - p0, e2 = p2 - p0;

        Vector3f pvec = cross(ray.d, e2);
        Float inv_det = rcp(dot(e1, pvec));

        Vector3f tvec = ray.o - p0;
        Float u = dot(tvec, pvec) * inv_det;
        active &= u >= 0.f && u <= 1.f;

        Vector3f qvec = cross(tvec, e1);
        Float v = dot(ray.d, qvec) * inv_det;
        active &= v >= 0.f && u + v <= 1.f;

        Float t = dot(e2, qvec) * inv_det;
        active &= t >= ray.mint && t <= ray.maxt;

        return { active, u, v, t };
    }
    std::pair<Mask, Float> ray_intersect(const Ray3f &ray_, Float *cache,
                                         Mask active) const override {
        MTS_MASK_ARGUMENT(active);
        Ray3f ray = m_to_object.transform_affine(ray_);
        if (H != 0) {
            Float t = sphereTrace(ray, active);
            Point3f local = ray(t);

            // Is intersection within ray segment and rectangle?
            active = active && t >= ray.mint && t <= ray.maxt &&
                     local.z() <= W && local.y() <= L && local.y() >= 0.f &&
                     local.z() >= 0.f;

            t = select(active, t, Float(math::Infinity<Float>));

            if (cache) {
                masked(cache[0], active) = t;
            }

            return { active, t };

        } else {
            Point3f p0 = Point3f(0, 0, 0), p1 = Point3f(0, L, 0),
                    p2 = Point3f(0, 0, W), p3 = Point3f(0, L, W);

            Point2f uv0 = Point2f(0, 0), uv1 = Point2f(1, 0),
                    uv2 = Point2f(0, 1), uv3 = Point2f(1, 1);

            auto [sucess_tri_1, tri1_u, tri1_v, tri1_t] =
                ray_intersect_triangle(ray, p0, p2, p3);
            if (sucess_tri_1) {
                cache[0] = tri1_t;
                return { true, tri1_t };
            }
            auto [sucess_tri_2, tri2_u, tri2_v, tri2_t] =
                ray_intersect_triangle(ray, p0, p3, p1);
            if (sucess_tri_2) {
                cache[0] = tri2_t;
                return { true, tri2_t };
            }
            return { false, 0 };
        }
    }

    Mask ray_test(const Ray3f &ray_, Mask active) const override {
        MTS_MASK_ARGUMENT(active);

        Ray3f ray = m_to_object.transform_affine(ray_);

        Point3f p0 = Point3f(0, 0, 0), p1 = Point3f(0, L, 0),
                p2 = Point3f(0, 0, W), p3 = Point3f(0, L, W);

        Point2f uv0 = Point2f(0, 0), uv1 = Point2f(1, 0), uv2 = Point2f(0, 1),
                uv3 = Point2f(1, 1);

        auto [sucess_tri_1, tri1_u, tri1_v, tri1_t] =
            ray_intersect_triangle(ray, p0, p2, p3);
        if (sucess_tri_1) {
            return { true };
        }
        auto [sucess_tri_2, tri2_u, tri2_v, tri2_t] =
            ray_intersect_triangle(ray, p0, p3, p1);
        if (sucess_tri_2) {
            return { true };
        }
        return { false };
    }

    void fill_surface_interaction(const Ray3f &ray_, const Float *cache,
                                  SurfaceInteraction3f &si_out,
                                  Mask active) const override {
        MTS_MASK_ARGUMENT(active);

        // get intersection info
        Float t = cache[0];
        SurfaceInteraction3f si(si_out);
        si.p = ray_(t);
        si.t = t;
        si.time = ray_.time;
        Point3f local = m_to_object * si.p;
        local =
            Point3f(-H * bicubic(Point2f(local.y() * invL, local.z() * invW),
                                 active, 0),
                    local.y(), local.z());

        // surfaceInteraction fill
        si.uv = Point2f(local.y() * invL, local.z() * invW);
        Vector3f grad_local = grad_h(si.uv, active);
        si.n = m_to_world * normalize(flip * grad_local);

        si.sh_frame.n = si.n;

        Normal3f dp_du =
            m_to_world *
            Normal3f(-H * N * bicubic(si.uv, active, 1), N * du_scale, 0);
        Normal3f dp_dv =
            m_to_world *
            Normal3f(-H * N * bicubic(si.uv, active, 2), 0, N * du_scale);
        si.dp_du = dp_du;
        si.dp_dv = dp_dv;

        si_out[active] = si;
    }

    SurfaceInteraction3f
    eval_surfaceInteraction_from_uv(Point2f p, Mask active,
                                    SurfaceInteraction3f &si_) const override {
        MTS_MASK_ARGUMENT(active);
        SurfaceInteraction3f si;
        Float f = eval_texture(0, p);
        si.p = m_to_world.transform_affine(Point3f(f, p.x() * L, p.y() * W));

        Vector3f grad_local = grad_h(p, active);
        si.n = m_to_world * normalize(flip * grad_local);
        si.dp_du = m_to_world *
                   Normal3f(-H * N * bicubic(p, active, 1), N * du_scale, 0);
        si.dp_dv = m_to_world *
                   Normal3f(-H * N * bicubic(p, active, 2), 0, N * du_scale);
        si.uv = p;

        return si;
    }

    Array<Float, 6> eval_curvature_du(Point2f p, Mask active, bool full) const {
        // eval point is outside of heightfield definition domain
        if (p.x() > 1 || p.y() > 1 || p.x() < 0 || p.y() < 0) {
            return 0.0f;
        }

        Float x = p.x();
        Float y = p.y();
        Float x_ = x * Float(N);
        Float y_ = y * Float(N);
        UInt32 X = (floor(x_));
        UInt32 Y = (floor(y_));
        Float dx = (x_ - X);
        Float dy = (y_ - Y);

        Float d2y = dy * dy, d3y = d2y * dy, d2x = dx * dx, d3x = d2x * dx;

        X = min(max(X, 1), N - 3);
        Y = min(max(Y, 1), N - 3);

        Array<Float, 4> vY(1.f, dy, d2y, d3y);
        Array<Float, 4> vX(1.f, dx, d2x, d3x);

        Array<Float, 16> coefs(0);

        coefs = load<Array<Float, 16>>(cubic_coefs + X + N * Y);
        Matrix<Float, 4> A(coefs[0], coefs[1], coefs[2], coefs[3], coefs[4],
                           coefs[5], coefs[6], coefs[7], coefs[8], coefs[9],
                           coefs[10], coefs[11], coefs[12], coefs[13],
                           coefs[14], coefs[15]);

        Float f = -H * dot(vX, A * vY);

        Float du =
            -H * dot(Array<Float, 4>(0.f, 1.f, 2 * dx, 3.f * d2x), A * vY);
        Float dv =
            -H * dot(vX, A * Array<Float, 4>(0.f, 1.f, 2 * dy, 3.f * d2y));

        Float duu = -H * dot(Array<Float, 4>(0.f, 0.f, 2.f, 6.f * dx), A * vY);
        Float dvv = -H * dot(vX, A * Array<Float, 4>(0.f, 0.f, 2.f, 6.f * dy));
        Float duv =
            -H * dot(Array<Float, 4>(0.f, 1.f, 2.f * dx, 3.f * d2x),
                     A * Array<Float, 4>(0.f, 1.f, 2.f * dy, 3.f * d2y));

        return Array<Float, 6>(0, du, dv, duu, dvv, duv);
    }

    std::pair<Vector3f, Vector3f>
    normal_derivative(const SurfaceInteraction3f &si, bool /*shading_frame*/,
                      Mask active) const override {
        MTS_MASK_ARGUMENT(active);
        if (H != 0) {
            Normal3f dn_du, dn_dv;

            Point3f p = m_to_object.transform_affine(si.p);

            // Vector2f uv(p.y() * invL, p.z() * invW);

            Array<Float, 6> curvature = eval_curvature_du(si.uv, active, true);

            Float du = curvature[1];
            Float dv = curvature[2];
            Float dudu = curvature[3];
            Float dvdv = curvature[4];
            Float dvdu = curvature[5];

            // gradient
            Normal3f V_du = Normal3f(N * du, N * du_scale, 0.f);
            Normal3f V_dv = Normal3f(N * dv, 0.f, N * du_scale);

            Vector3f n = flip * normalize(cross(V_du, V_dv));

            // gradient 2nd order
            Vector3f V_dudu = N * N * Vector3f(dudu, 0.f, 0.f);
            Vector3f V_dvdu = N * N * Vector3f(dvdu, 0.f, 0.f);
            Vector3f V_dvdv = N * N * Vector3f(dvdv, 0.f, 0.f);
            // n = (cross(V_du, V_dv));

            // first fondamental form
            Float E = dot((V_du), V_du);
            Float F = dot((V_du), V_dv);
            Float G = dot((V_dv), V_dv);

            // second fondamental form
            Float e = dot(n, V_dudu);
            Float f = dot(n, V_dvdu);
            Float g = dot(n, V_dvdv);

            // Weingarten eq:
            Float inv_EGtimesFF = rcp(E * G - F * F);
            dn_du = (((f * F - e * G) * inv_EGtimesFF * (V_du)) +
                     ((e * F - f * E) * inv_EGtimesFF * (V_dv)));
            dn_dv = (((g * F - f * G) * inv_EGtimesFF * (V_du)) +
                     ((f * F - g * E) * inv_EGtimesFF * (V_dv)));

            return { m_to_world * dn_du, m_to_world * dn_dv };
        } else {
            return { Normal3f(0), Normal3f(0) };
        }
    }

    //! @}
    // =============================================================

    ScalarSize primitive_count() const override { return 1; }

    ScalarSize effective_primitive_count() const override { return 1; }

    void traverse(TraversalCallback *callback) override {
        Base::traverse(callback);
    }

    void
    parameters_changed(const std::vector<std::string> & /*keys*/) override {
        update();
        Base::parameters_changed();
#if defined(MTS_ENABLE_OPTIX)
        optix_prepare_geometry();
#endif
    }

#if defined(MTS_ENABLE_OPTIX)
    using Base::m_optix_data_ptr;

    void optix_prepare_geometry() override {
        if constexpr (is_cuda_array_v<Float>) {
            if (!m_optix_data_ptr)
                m_optix_data_ptr = cuda_malloc(sizeof(OptixHeightfieldData));

            OptixHeightfieldData data = {
                bbox(),   m_to_world, m_to_object,
                m_center, m_radius,   m_flip_normals
            };

            cuda_memcpy_to_device(m_optix_data_ptr, &data,
                                  sizeof(OptixHeightfieldData));
        }
    }
#endif

    std::string to_string() const override {
        std::ostringstream oss;
        oss << "Heightfield[" << std::endl
            << "  to_world = " << string::indent(m_to_world) << "," << std::endl
            << "  center = " << m_center << "," << std::endl
            << "  H = " << H << "," << std::endl
            << "  W = " << W << "," << std::endl
            << "  L = " << L << "," << std::endl
            << "  flip = " << flip << "," << std::endl
            << "  surface_area = " << surface_area() << "," << std::endl
            << "  " << string::indent(get_children_string()) << std::endl
            << "]";
        return oss.str();
    }

    MTS_DECLARE_CLASS()

private:
    /// Center in world-space
    ScalarPoint3f m_center;
    // dimensions
    ScalarFloat L, W, H;
    // frame
    ScalarFrame3f m_frame;
    ScalarFloat m_inv_surface_area;

    // pair reference
    int m_self_id; // name of the heightfield (ex: window front)
    int m_pair_id; // id of the associated heightfield (ex: window back)

    // flip normal
    bool m_flip_normals;
    ScalarFloat flip;
    size_t N;
    ScalarFloat scale_l, scale_w; // scale factor to normalize local coordinates
                                  // ([0,L] --> [0,1])
    ScalarFloat invL, invW;
    Float dx_scale, dy_scale, dx_scale_sqr, dy_scale_sqr, du_scale, dv_scale,
        du_scale_sqr;

    std::string reconstructionMode;
    // heightmap
    ref<Texture> heightmap;
    Array<Float, 16> *cubic_coefs;
    Float partition_size, inv_partition_size, inv_partition_area;

    // number of heighfield evaluation
    mutable long nbEval;
};

MTS_IMPLEMENT_CLASS_VARIANT(Heightfield, Shape)
// MTS_EXTERN_CLASS_RENDER(Heightfield)
MTS_EXPORT_PLUGIN(Heightfield, "Heightfield intersection primitive");
// MTS_INSTANTIATE_CLASS(Heightfield)
NAMESPACE_END(mitsuba)
