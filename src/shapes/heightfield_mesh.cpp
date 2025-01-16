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
   - |Float|
   - length (|) of the heightfield (Default: 0)
* - W
   - |Float|
   - width (---) of the heightfield (Default: 0)
 * - H
   - |Float|
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

template <typename Float, typename Spectrum> class MTS_EXPORT_RENDER Heightfield final : public Shape<Float, Spectrum> {
public:
    MTS_IMPORT_BASE(Shape, m_to_world, m_to_object, set_children, get_children_string)
    MTS_IMPORT_TYPES(Texture)

    using typename Base::ScalarSize;

    void generate_coef(ref<Texture> heightmap, Array<Float, 1> **height_data, Array<Normal3f, 1> **normal_data) {

        Float scale = rcp((Float) (N));

        // precompute bicubic interpolator coefficients:
        // interpolator coefficients
        *height_data = new Array<Float, 1>[N * N];
        *normal_data = new Array<Normal3f, 1>[N * N];
        // for each pixels:
        for (size_t i = 0; i < N; i++) {
            for (size_t j = 0; j < N; j++) {

                // account for border cells
                int X = i, Y = j;

                Y = max(min(Y, N - 2), 1);
                X = max(min(X, N - 2), 1);
                // Fill FXY: get the 16 reference texture points
                SurfaceInteraction3f si;
                si.uv     = Point2f(X + 0, Y + 0) * scale;
                Float F00 = (-H * 0.99f) * heightmap->eval_1(si, "bilinear");
                // store everything:
                (*(*height_data + i + N * j))[0] = F00;
            }
        }

        for (size_t i = 0; i < N; i++) {
            for (size_t j = 0; j < N; j++) {

                // account for border cells
                int X = i, Y = j;

                X = max(min(X, N - 2), 1);
                Y = max(min(Y, N - 2), 1);
                // Fill FXY: get the 16 reference texture points
                Normal3f vn = vertexNormal(X, Y);

                // store everything:

                (*(*normal_data + i + N * j))[0] = Normal3f(vn.x(), vn.y(), vn.z());
            }
        }

        return;
    }

    Heightfield(const Properties &props) : Base(props) {
        /// Are the heightfield normals pointing inwards? default: no
        m_flip_normals = props.bool_("flip_normals", false);
        flip           = select(!m_flip_normals, 1.f, -1.f);
        // dimensions
        L = props.float_("L", 1);
        W = props.float_("W", 1);
        H = props.float_("H", 1);

        // load texture
        heightmap = props.texture<Texture>("heightmap", .5f);
        N         = heightmap->resolution().x();
        // precompute scaling factors
        invW     = rcp(W);
        invL     = rcp(L);
        dv_scale = invW * (Float) (N);
        du_scale = invL * (Float) (N);
        dh_scale = rcp(H);
        // select surface reconstruction mode (bilinear or bicubic)
        reconstructionMode = props.string("reconstructionMode", "bilinear");

        Float scale = rcp((Float) (N));
        scale_w     = scale; // same since heightmap is a N*N texture
        scale_l     = scale; // same since heightmap is a N*N texture

        generate_coef(heightmap, &height_data, &normal_data);

        // define partition size (for multi-piece glass panels)
        partition_size     = rcp(props.float_("partition_number", 1.0));
        inv_partition_size = rcp(partition_size);
        inv_partition_area = sqr(sqr(inv_partition_size));

        // reference to heightfield pair
        m_self_id = props.int_("self_id", 0);
        m_pair_id = props.int_("pair_id", 0);

        update();
        set_children();
        // std::cout << this << std::endl;
    }

    // ~Heightfield() { std::cout << "destroy heigfield: eval number = " << nbEvall << std::endl; }

    void update() {
        m_to_object = m_to_world.inverse();

        auto [S, Q, T] = transform_decompose(m_to_world.matrix);

        if (abs(S[0][1]) > 1e-6f || abs(S[0][2]) > 1e-6f || abs(S[1][0]) > 1e-6f || abs(S[1][2]) > 1e-6f || abs(S[2][0]) > 1e-6f ||
            abs(S[2][1]) > 1e-6f)
            Log(Warn, "'to_world' transform shouldn't contain any shearing!");

        m_center = T;

        ScalarVector3f dp_du  = m_to_world * ScalarVector3f(0.f, 2.f, 0.f);
        ScalarVector3f dp_dv  = m_to_world * ScalarVector3f(0.f, 0.f, 2.f);
        ScalarNormal3f normal = -flip * normalize(m_to_world * ScalarNormal3f(-1.f, 0.f, 0.f));

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

    // since H is generally very small, surface_area() could be aproximated as L*W
    ScalarFloat surface_area() const override {

        Float S = L * W;

        // // more accurate surface estimation
        // int Nr   = N;
        // Float dy = rcp(Float(Nr));
        // Float dx = rcp(Float(Nr));
        // Float S  = 0;
        // for (int i = 0; i < Nr; i++) {
        //     for (int j = 0; j < Nr; j++) {
        //         Float du = H * bicubic(Point2f(i * dx, j * dy), true, 1);
        //         Float dv = H * bicubic(Point2f(i * dx, j * dy), true, 2);
        //         S        = S + L * dx * W * dy * sqrt(du * du + dv * dv + 1);
        //     }
        // }
        return S;
    }

    // =============================================================
    //! @{ \name Sampling routines
    // =============================================================

    PositionSample3f sample_position(Float time, const Point2f &sample, Mask active) const override {
        MTS_MASK_ARGUMENT(active);

        PositionSample3f ps;

        Float scale_x = L * scale_l;
        Float scale_y = W * scale_w;

        Float x_ = sample.x() * L * rcp(scale_x);
        Float y_ = sample.y() * W * rcp(scale_y);
        UInt32 Y = (floor(x_));
        UInt32 Z = (floor(y_));

        Float u_proj = x_ - Y; // projected yz coordinate in object space, inside triangle
        Float v_proj = y_ - Z;

        // triangle contruction
        Y = max(min(Y, N - 1), 1);
        Z = max(min(Z, N - 1), 1);
        Float H00 = load<Array<Float, 1>>(height_data + Y + N * Z)[0];
        Float H10 = load<Array<Float, 1>>(height_data + Y + 1 + N * Z)[0];
        Float H01 = load<Array<Float, 1>>(height_data + Y + N * (Z + 1))[0];
        Float H11 = load<Array<Float, 1>>(height_data + (Y + 1) + N * (Z + 1))[0];

        Float v_ = Y / du_scale, v1 = (Y + 1) / du_scale, w_ = Z / dv_scale,
              w1   = (Z + 1) / dv_scale; // voxel coordinates into object local coordinates
        Point3f P0 = m_to_world.transform_affine(Point3f((H00), v_, w_)), P1 = m_to_world.transform_affine(Point3f((H10), v1, w_)),
                P2 = m_to_world.transform_affine(Point3f((H01), v_, w1)), P3 = m_to_world.transform_affine(Point3f((H11), v1, w1));

        Normal3f N0 = vertexNormal_2(Y, Z), N1 = vertexNormal_2(Y + 1, Z), N2 = vertexNormal_2(Y, Z + 1), N3 = vertexNormal_2(Y + 1, Z + 1);

        Point3f p0, p1, p2;
        Normal3f n0, n1, n2;
        if (v_proj < 1 - u_proj) { // upper triangle
            p0 = P0, p1 = P2, p2 = P3;
            n0 = N0, n1 = N2, n2 = N3;
        } else { // lower triangle
            p0 = P0, p1 = P3, p2 = P1;
            n0 = N0, n1 = N3, n2 = N1;
        }

        // postion sample contruction
        Float b1 = u_proj, b2 = v_proj;

        Float b0 = 1.f - b1 - b2;

        Vector3f dp0 = p1 - p0, dp1 = p2 - p0;

        // Re-interpolate intersection using barycentric coordinates
        ps.p = (p0 * b0 + p1 * b1 + p2 * b2); // H being quite small, this is considered
                                              // uniform sampling

        // ps.p =  m_to_world.transform_affine(Point3f(0, 0.5, 0.5));
        // std::cout << ps.p  << "H = " << H << std::endl;
        ps.n = m_to_world.transform_affine(-flip * normalize(n0 * b1 + n1 * b2 + n2 * b0));

        // ps.p     = m_to_world.transform_affine(Vector3f(0, sample.x() * L, sample.y() * W));
        ps.n      = m_to_world.transform_affine(Vector3f(-1, 0, 0));
        ps.pdf    = m_inv_surface_area;
        ps.uv     = Vector2f(0, 0);
        ps.time   = time;
        ps.delta  = false;
        ps.object = this;

        // std::cout << H << std::endl;

        return ps;

        // //--------------------------------------------------------
        // PositionSample3f ps;
        // Float S_x = 0;
        // ps.p      = m_to_world.transform_affine(Point3f(S_x, sample.x() * L, sample.y() * W)); // H being quite small, this is considered
        //                                                                                        // uniform sampling
        // ps.n     = (m_to_world * normalize(-flip * grad_h(Point3f(S_x, sample.x() * L, sample.y() * W), active)));
        // ps.pdf   = m_inv_surface_area;
        // ps.uv    = sample;
        // ps.time  = time;
        // ps.delta = false;
        // return ps;
    }

    Float pdf_position(const PositionSample3f & /*ps*/, Mask active) const override {
        MTS_MASK_ARGUMENT(active);
        return m_inv_surface_area;
    }

    // =============================================================
    //! @{ \name attribute get (messy...)
    // =============================================================
    Float eval_attribute_1(const std::string &name, const SurfaceInteraction3f &si, Mask active = true) const override {
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
        else
            return 0.0f;
    }

    // =============================================================
    //! @{ \name Window heightfield pair reference
    // =============================================================
    int get_pair_id(Mask /*active*/) const override { return m_pair_id; }

    int get_self_id(Mask /*active*/) const override { return m_self_id; }

    // TODO refactor
    Float eval_texture(int index, Point2f p, Mask active = true) const override {
        MTS_MASK_ARGUMENT(active);
        // if (index == 0) {
        //     // return h_f(p, active);
        //     return -H * bicubic(p, active, index);

        // } else if (index == 1) {
        //     // return duf(p, active);
        //     return -H * du_scale * bicubic(p, active, index);

        // } else if (index == 2) {
        //     // return dvf(p, active);1
        //     return -H * dv_scale * bicubic(p, active, index);

        // } else if (index == 3) {
        //     // return duuf(p, active);
        //     return -H * du_scale * du_scale * bicubic(p, active, index);

        // } else if (index == 4) {
        //     // return dvvf(p, active);
        //     return -H * dv_scale * dv_scale * bicubic(p, active, index);

        // } else if (index == 5) {
        //     // return dvuf(p, active);
        //     return -H * du_scale * dv_scale * bicubic(p, active, index);
        // } else {
        //     std::cout << "error, Heighfield eval_texture() <name> " << index << "is not a valid parameter" << std::endl;
        // }
        return 0.0f;
    }

    // =============================================================
    //! @{ \name Ray tracing routines
    // =============================================================

#define SIGN(x) (x > 0 ? 1 : (x < 0 ? -1 : 0))
#define FRAC0(x) (x - std::floorf(x))
#define FRAC1(x) (1 - x + std::floorf(x))
    std::pair<Mask, Float> ray_intersect(const Ray3f &ray_, Float *cache, Mask active) const override {
        MTS_MASK_ARGUMENT(active);
        // Ray3f ray     = m_to_object.transform_affine(ray_);
        // Float t       = sphereTrace(ray, active);
        // Point3f local = ray(t);

        // // Is intersection within ray segment and rectangle?
        // active = active && t >= ray.mint && t <= ray.maxt && local.z() <= W && local.y() <= L && local.y() >= 0.f && local.z() >= 0.f;

        // t = select(active, t, Float(math::Infinity<Float>));

        // if (cache) {
        //     masked(cache[0], active) = t;
        // }

        // return { active, t };

        // ----------------------------------------------------------------
        // Ray in local space
        Ray rayLocal = m_to_object.transform_affine(ray_);

        // triangle intersection structure
        SurfaceInteraction3f triIsect1_HF, triIsect2_HF;
        SurfaceInteraction3f intr;

        // test against bounding box: voxel traversal starting point
        Point3f lb = Point3f(0.1, 0, 0), ub = Point3f(-H, L, W);
        Float t_In = 0, t_Out = 0;
        ScalarBoundingBox3f bounds;
        bounds.expand(ub);
        bounds.expand(lb);
        auto [b_hit, mint, maxt] = bounds.ray_intersect(rayLocal);
        t_In                     = mint;
        t_Out                    = maxt;
        // std::cout << "valid box ?" << bounds<< std::endl;
        // std::cout << mint << " " << maxt << std::endl;
        if (!b_hit)
            return {};

        if (H != 0.0f) {
            //----------------------------------------------
            // voxel traversal
            //----------------------------------------------

            Float tMaxX = 0, tMaxY = 0, tMaxZ = 0, tDeltaX = 0, tDeltaY = 0, tDeltaZ = 0;
            int X = 0, Y = 0, Z = 0;
            Point3f In  = rayLocal(t_In);
            Point3f Out = rayLocal(t_Out);
            // account for voxel length != 1
            In  = Point3f(In.x() * dh_scale, In.y() * du_scale, In.z() * dv_scale);
            Out = Point3f(Out.x() * dh_scale, Out.y() * du_scale, Out.z() * dv_scale);

            // X init
            int dx = SIGN(Out.x() - In.x());
            if (dx != 0)
                tDeltaX = std::fmin(dx / (Out.x() - In.x()), 10000000.0f);
            else
                tDeltaX = 10000000.0f;
            if (dx > 0)
                tMaxX = tDeltaX * FRAC1(In.x());
            else
                tMaxX = tDeltaX * FRAC0(In.x());
            X = (int) std::floorf(In.x());
            // Y init
            int dy = SIGN(Out.y() - In.y());
            if (dy != 0)
                tDeltaY = std::fmin(dy / (Out.y() - In.y()), 10000000.0f);
            else
                tDeltaY = 10000000.0f;
            if (dy > 0)
                tMaxY = tDeltaY * FRAC1(In.y());
            else
                tMaxY = tDeltaY * FRAC0(In.y());
            Y = (int) (std::floorf(In.y()));
            // Z init
            int dz = SIGN(Out.z() - In.z());
            if (dz != 0)
                tDeltaZ = std::fmin(dz / (Out.z() - In.z()), 1000000000.0f);
            else
                tDeltaZ = 1000000000.0f;
            if (dz > 0)
                tMaxZ = tDeltaZ * FRAC1(In.z());
            else
                tMaxZ = tDeltaZ * FRAC0(In.z());
            Z = (int) (std::floorf(In.z()));

            // std::cout << "tMaxX = " << tMaxX << std::endl;
            // std::cout << "tMaxY = " << tMaxY << std::endl;
            // std::cout << "tMaxZ = " << tMaxZ << std::endl;

            // traversal
            while (true) { // while
                if (tMaxX < tMaxY) {
                    if (tMaxX < tMaxZ) {
                        X += dx;
                        tMaxX += tDeltaX;
                    } else {
                        Z += dz;
                        tMaxZ += tDeltaZ;
                    }
                } else {
                    if (tMaxY < tMaxZ) {
                        Y += dy;
                        tMaxY += tDeltaY;
                    } else {
                        Z += dz;
                        tMaxZ += tDeltaZ;
                    }
                }
                if ((tMaxX) > 1 && (tMaxY) > 1 && (tMaxZ) > 1)
                    break; // outside of voxel grid

                // for (int Y = 0; Y < N; Y++) {
                //     for (int Z = 0; Z < N; Z++) {

                SurfaceInteraction3f si;

                Z         = max(min(Z, N - 1), 1);
                Y         = max(min(Y, N - 1), 1);
                Float H00 = load<Array<Float, 1>>(height_data + Y + N * Z)[0];
                Float H10 = load<Array<Float, 1>>(height_data + Y + 1 + N * Z)[0];
                Float H01 = load<Array<Float, 1>>(height_data + Y + N * (Z + 1))[0];
                Float H11 = load<Array<Float, 1>>(height_data + (Y + 1) + N * (Z + 1))[0];

                Float v = Y / du_scale, v1 = (Y + 1) / du_scale, w = Z / dv_scale,
                      w1 = (Z + 1) / dv_scale; // voxel coordinates into object local coordinates

                Float voxelHeight    = std::fmax(std::fmax(std::fmax(H00, H10), H01), H11); // maximum height in cell
                Float minvoxelHeight = std::fmin(std::fmin(std::fmin(H00, H10), H01), H11); // minimum height in cell

                lb = Point3f(minvoxelHeight, (Y + 1) / du_scale, (Z + 1) / dv_scale);
                ub = Point3f(voxelHeight, (Y + 2) / du_scale, (Z + 2) / dv_scale);
                Float dummy;
                ScalarBoundingBox3f bounds_2;
                bounds_2.expand(lb);
                bounds_2.expand(ub);
                // auto [b_hit_2, mint_2, maxt_2] = bounds_2.ray_intersect(rayLocal);
                // if (!b_hit_2) // test against voxel bounding box
                // {
                Point3f p0 = m_to_world.transform_affine(Point3f((H00), v, w)), p1 = m_to_world.transform_affine(Point3f((H10), v1, w)),
                        p2 = m_to_world.transform_affine(Point3f((H01), v, w1)), p3 = m_to_world.transform_affine(Point3f((H11), v1, w1));

                // Point3f p0 = m_to_world.transform_affine(Point3f(0, 0, 0)), p1 = m_to_world.transform_affine(Point3f(0, L, 0)),
                // p2 = m_to_world.transform_affine(Point3f(0, 0, W)), p3 = m_to_world.transform_affine(Point3f(0, L, W));

                Point2f uv0 = Point2f(0, 0), uv1 = Point2f(1, 0), uv2 = Point2f(0, 1), uv3 = Point2f(1, 1);

                // triangle patch intersection
                auto [sucess_tri_1, tri1_u, tri1_v, tri1_t] = ray_intersect_triangle(ray_, p0, p2, p3);
                if (sucess_tri_1) {
                    // std::cout << "test " << ray_(tri1_t) << std::endl;
                    cache[0] = tri1_t;
                    cache[1] = Y;
                    cache[2] = Z;
                    cache[3] = 0.0f; // low triangle
                    cache[4] = tri1_u;
                    cache[5] = tri1_v;
                    return { true, tri1_t };
                }
                auto [sucess_tri_2, tri2_u, tri2_v, tri2_t] = ray_intersect_triangle(ray_, p0, p3, p1);
                if (sucess_tri_2) {
                    cache[0] = tri2_t;
                    cache[1] = Y;
                    cache[2] = Z;
                    cache[3] = 1.0f; // high triangle
                    cache[4] = tri2_u;
                    cache[5] = tri2_v;
                    return { true, tri2_t };
                }
                // }
                // }
            } // if in voxel
            // } // while
        } else {
            Point3f p0 = m_to_world.transform_affine(Point3f(0, 0, 0)), p1 = m_to_world.transform_affine(Point3f(0, L, 0)),
                    p2 = m_to_world.transform_affine(Point3f(0, 0, W)), p3 = m_to_world.transform_affine(Point3f(0, L, W));

            Point2f uv0 = Point2f(0, 0), uv1 = Point2f(1, 0), uv2 = Point2f(0, 1), uv3 = Point2f(1, 1);

            auto [sucess_tri_1, tri1_u, tri1_v, tri1_t] = ray_intersect_triangle(ray_, p0, p2, p3);
            if (sucess_tri_1) {
                // std::cout << "test " << tri1_t << std::endl;
                cache[0] = tri1_t;
                cache[1] = 0;
                cache[2] = 0;
                cache[3] = 0.0f; // low triangle
                cache[4] = tri1_u;
                cache[5] = tri1_v;

                return { true, tri1_t };
            }
            auto [sucess_tri_2, tri2_u, tri2_v, tri2_t] = ray_intersect_triangle(ray_, p0, p3, p1);
            if (sucess_tri_2) {
                cache[0] = tri2_t;
                cache[1] = 0;
                cache[2] = 0;
                cache[3] = 1.0f; // high triangle
                cache[4] = tri2_u;
                cache[5] = tri2_v;
                return { true, tri2_t };
            }
            return { false, 0 };
        }

        return { false, 0 };
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
    MTS_INLINE std::tuple<Mask, Float, Float, Float> ray_intersect_triangle(const Ray3f &ray, Point3f p0, Point3f p1, Point3f p2,
                                                                            identity_t<Mask> active = true) const {

        Vector3f e1 = p1 - p0, e2 = p2 - p0;

        Vector3f pvec = cross(ray.d, e2);
        Float inv_det = rcp(dot(e1, pvec));

        Vector3f tvec = ray.o - p0;
        Float u       = dot(tvec, pvec) * inv_det;
        active &= u >= 0.f && u <= 1.f;

        Vector3f qvec = cross(tvec, e1);
        Float v       = dot(ray.d, qvec) * inv_det;
        active &= v >= 0.f && u + v <= 1.f;

        Float t = dot(e2, qvec) * inv_det;
        active &= t >= ray.mint && t <= ray.maxt;

        return { active, u, v, t };
    }

    Normal3f vertexNormal(int Y, int Z) const {
        SurfaceInteraction3f si;
        Y = max(min(Y, N - 2), 1);
        Z = max(min(Z, N - 2), 1);
        si.uv       = Point2f(Y - 1, Z) * scale_w;
        Float Zleft = load<Array<Float, 1>>(height_data + Y - 1 + N * Z)[0];
        ;

        si.uv        = Point2f(Y + 1, Z) * scale_w;
        Float Zright = load<Array<Float, 1>>(height_data + Y + 1 + N * Z)[0];
        ;

        si.uv          = Point2f(Y - 1, Z + 1) * scale_w; //
        Float Zupright = load<Array<Float, 1>>(height_data + Y - 1 + N * (Z + 1))[0];
        ;

        si.uv           = Point2f(Y + 1, Z - 1) * scale_w; //
        Float Zdownleft = load<Array<Float, 1>>(height_data + Y + 1 + N * (Z - 1))[0];
        ;

        si.uv       = Point2f(Y, Z - 1) * scale_w;
        Float Zdown = load<Array<Float, 1>>(height_data + Y + N * (Z - 1))[0];
        ;

        si.uv     = Point2f(Y, Z + 1) * scale_w;
        Float Zup = load<Array<Float, 1>>(height_data + Y + N * (Z + 1))[0];
        ;

        si.uv    = Point2f(Y, Z) * scale_w;
        Float Z0 = load<Array<Float, 1>>(height_data + Y + N * Z)[0];
        ;

        Point3f P_upright  = Point3f(Zupright, (Y - 1) / du_scale, (Z + 1) / dv_scale);
        Point3f P_up       = Point3f(Zup, (Y + 0) / du_scale, (Z + 1) / dv_scale);
        Point3f P_right    = Point3f(Zright, (Y + 1) / du_scale, (Z + 0) / dv_scale);
        Point3f P_downleft = Point3f(Zdownleft, (Y + 1) / du_scale, (Z - 1) / dv_scale);
        Point3f P_down     = Point3f(Zdown, (Y + 0) / du_scale, (Z - 1) / dv_scale);
        Point3f P_left     = Point3f(Zleft, (Y - 1) / du_scale, (Z + 0) / dv_scale);
        Point3f P0         = Point3f(Z0, (Y + 0) / du_scale, (Z + 0) / dv_scale);

        Normal3f N1 = (cross(P_upright - P0, P_up - P0));
        Normal3f N2 = (cross(P_up - P0, P_right - P0));
        Normal3f N3 = (cross(P_right - P0, P_downleft - P0));
        Normal3f N4 = (cross(P_downleft - P0, P_down - P0));
        Normal3f N5 = (cross(P_down - P0, P_left - P0));
        Normal3f N6 = (cross(P_left - P0, P_upright - P0));

        Normal3f N = (N1 + N2 + N3 + N4 + N5 + N6);

        if (H == 0.0f) {
            N = Normal3f(-1, 0, 0);
        }

        return N;
    }

    Normal3f vertexNormal_2(int Y, int Z) const {

        Y = max(min(Y, N - 1), 1);
        Z = max(min(Z, N - 1), 1);

        Array<Normal3f, 1> array = load<Array<Normal3f, 1>>(normal_data + Y + N * Z);
        Normal3f N0              = array[0];
        return N0;
    }

    Array<Float, 6> eval_curvature(Point2f p, Mask active, bool full) const {

        Float scale_x = L * scale_l;
        Float scale_y = W * scale_w;

        Float x_ = p.x() * rcp(scale_x);
        Float y_ = p.y() * rcp(scale_y);
        UInt32 Y = (floor(x_));
        UInt32 Z = (floor(y_));

        Float u_proj = x_ - Y; // projected yz coordinate in object space, inside triangle
        Float v_proj = y_ - Z;

        // triangle contruction
        Y = max(min(Y, N - 1), 1);
        Z = max(min(Z, N - 1), 1);

        Float H00 = load<Array<Float, 1>>(height_data + Y + N * Z)[0];
        Float H10 = load<Array<Float, 1>>(height_data + Y + 1 + N * Z)[0];
        Float H01 = load<Array<Float, 1>>(height_data + Y + N * (Z + 1))[0];
        Float H11 = load<Array<Float, 1>>(height_data + (Y + 1) + N * (Z + 1))[0];

        Float v0 = Y / du_scale, v1 = (Y + 1) / du_scale, w0 = Z / dv_scale,
              w1   = (Z + 1) / dv_scale; // voxel coordinates into object local coordinates
        Point3f P0 = m_to_world.transform_affine(Point3f((H00), v0, w0)), P1 = m_to_world.transform_affine(Point3f((H10), v1, w0)),
                P2 = m_to_world.transform_affine(Point3f((H01), v0, w1)), P3 = m_to_world.transform_affine(Point3f((H11), v1, w1));
        Normal3f N0 = vertexNormal_2(Y, Z), N1 = vertexNormal_2(Y + 1, Z), N2 = vertexNormal_2(Y, Z + 1), N3 = vertexNormal_2(Y + 1, Z + 1);

        Point3f p0, p1, p2;
        Normal3f n0, n1, n2;
        if (v_proj > u_proj) { // upper triangle
            p0 = P0, p1 = P2, p2 = P3;
            n0 = N0, n1 = N2, n2 = N3;
        } else { // lower triangle
            p0 = P0, p1 = P3, p2 = P1;
            n0 = N0, n1 = N3, n2 = N1;
        }

        Float u_ = u_proj, v_ = v_proj;
        Float w_ = 1.f - u_ - v_;
        Point3f localP = u_ * p1 + v_ * p2 + w_ * p0;

        Vector3f rel = localP - p0, du = p1 - p0, dv = p2 - p0;

        /* Solve a least squares problem to determine
           the UV coordinates within the current triangle */
        Float b1 = dot(du, rel), b2 = dot(dv, rel), a11 = dot(du, du), a12 = dot(du, dv), a22 = dot(dv, dv),
              inv_det = rcp(a11 * a22 - a12 * a12);

        Float u = fmsub(a22, b1, a12 * b2) * inv_det, v = fnmadd(a12, b1, a11 * b2) * inv_det, w = 1.f - u - v;

        /* Now compute the derivative of "normalize(u*n1 + v*n2 + (1-u-v)*n0)"
           with respect to [u, v] in the local triangle parameterization.

           Since d/du [f(u)/|f(u)|] = [d/du f(u)]/|f(u)|
             - f(u)/|f(u)|^3 <f(u), d/du f(u)>, this results in
        */
        if (full) {
            Normal3f N(u * n1 + v * n2 + w * n0);
            Float il = rsqrt(squared_norm(N));
            N *= il;

            Vector3f dndu = (n1 - n0) * il;
            Vector3f dndv = (n2 - n0) * il;

            dndu = fnmadd(N, dot(N, dndu), dndu);
            dndv = fnmadd(N, dot(N, dndv), dndv);

            Point2f uv0, uv1, uv2;
            uv0 = Point2f(0.f, 0.f);
            uv1 = Point2f(1.f, 0.f);
            uv2 = Point2f(0.f, 1.f);

            Vector3f dp0 = p1 - p0, dp1 = p2 - p0;
            Vector2f duv1 = uv1 - uv0, duv2 = uv2 - uv0;
            inv_det = rcp(duv1.x() * duv2.y() - duv1.y() * duv2.x());

            Vector3f dndu_ = (duv2.y() * dndu - duv1.y() * dndv) * inv_det, dndv_ = (-duv2.x() * dndu + duv1.x() * dndv) * inv_det;
            dndu = dndu_;
            dndv = dndv_;

            /*
            Compute dp_du and dp_dv
            */
            auto [dp_du, dp_dv] = coordinate_system(N);


            dp_du               = fmsub(duv2.y(), dp0, duv1.y() * dp1) * inv_det;
            dp_dv               = fnmadd(duv2.x(), dp0, duv1.x() * dp1) * inv_det;


            Float dp_duu = -dot(dndu, dp_du) * rcp(N.x());
            Float dp_dvv = -dot(dndv, dp_dv) * rcp(N.x());
            Float dp_dvu = -dot(dndv, dp_du) * rcp(N.x());

            return Array<Float, 6>(m_to_object.transform_affine(localP).x(), dp_du.x(), dp_dv.x(), dp_duu, dp_dvv, dp_dvu);
        } else {
            return Array<Float, 6>(m_to_object.transform_affine(localP).x(), 0, 0, 0, 0, 0);
        }
    }

    void fill_surface_interaction_triangle(Float u, Float v, Point3f p0, Point3f p1, Point3f p2, Normal3f n0, Normal3f n1, Normal3f n2,
                                           SurfaceInteraction3f &si, Mask active, Ray3f r, Vector2f uv0, Vector2f uv1 ,Vector2f uv2) const {
        if (H != 0) {
            // Barycentric coordinates within triangle
            Float b1 = u, b2 = v;

            Float b0 = 1.f - b1 - b2;

            Vector3f dp0 = p1 - p0, dp1 = p2 - p0;

            // Re-interpolate intersection using barycentric coordinates
            si.p = (p0 * b0 + p1 * b1 + p2 * b2);

            // Face normal
            Normal3f n = normalize(-flip * cross(dp0, dp1));
            si.n       = n;

            // Shading normal (if available)

            n    = m_to_world * (-flip * normalize(n0 * b0 + n1 * b1 + n2 * b2));
            si.n = n;
            // si.n       = m_to_world.transform_affine(-flip * Vector3f(-1, 0, 0));
            si.sh_frame.n = si.n;

            // Texture coordinates (if available)
            auto [dp_du, dp_dv] = coordinate_system(si.n);
            Point2f uv(0, 0);


        //    uv0 = Point2f(0.f, 0.f);
        //    uv1 = Point2f(1.f, 0.f);
        //    uv2 = Point2f(0.f, 1.f);

            uv = uv0 * b0 + uv1 * b1 + uv2 * b2;

            Vector2f duv0 = uv1 - uv0, duv1 = uv2 - uv0;

            Float det = fmsub(duv0.x(), duv1.y(), duv0.y() * duv1.x()), inv_det = rcp(det);

            Mask valid = neq(det, 0.f);

            dp_du[valid] = fmsub(duv1.y(), dp0, duv0.y() * dp1) * inv_det;
            dp_dv[valid] = fnmadd(duv1.x(), dp0, duv0.x() * dp1) * inv_det;

            si.uv = uv;

            // Tangents
            si.dp_du = dp_du;
            si.dp_dv = dp_dv;
        

	    //si = eval_surfaceInteraction_from_uv(si.uv, active, si);
	    //std::cout << si << std::endl;
	    // std::cout << si.dp_du << " || " << si.dp_dv << std::endl;

            //             if(this->is_caustic_caster_multi_scatter()){
            //     std::cout << "heightfield :" << si << std::endl;
            //                 std::cout << r << std::endl;
            // }
            // si.dp_du      = m_to_world * Vector3f(0, 1, 0);
            // si.dp_dv      = m_to_world * Vector3f(0, 0, 1);
        } else {
            // Barycentric coordinates within triangle
            Float b1 = u, b2 = v;

            Float b0 = 1.f - b1 - b2;

            Vector3f dp0 = p1 - p0, dp1 = p2 - p0;

            // Re-interpolate intersection using barycentric coordinates
            Point3f p_local = p0 * b0 + p1 * b1 + p2 * b2;

            si.p = p_local;
            // si.p            = m_to_world.transform_affine(p_local);

            Float scale_x = L * scale_l;
            Float scale_y = W * scale_w;

            Float x_ = p_local.y() * rcp(scale_x);
            Float y_ = p_local.z() * rcp(scale_y);
            UInt32 Y = (floor(x_));
            UInt32 Z = (floor(y_));

            // Face normal
            // // flat shading :
            // Normal3f n   = m_to_world.transform_affine(flip * normalize(cross(dp0, dp1)));
            // si.n[active] = -n;

            // smooth shading :
            Normal3f n = normalize(n0 * b0 + n1 * b1 + n2 * b2);
            si.n       = m_to_world.transform_affine(-flip * n);
            si.n       = m_to_world.transform_affine(-flip * Vector3f(-1, 0, 0));

            // Texture coordinates (if available)
            auto [dp_du, dp_dv] = coordinate_system(si.n);
            // si.sh_frame.s = dp_du;
            // si.sh_frame.t = dp_dv;
            Point2f uv(x_, y_);

            Point2f uv0, uv1, uv2;
            // if (has_vertex_texcoords()) {
            //     uv0 = vertex_texcoord(fi[0], active);
            //     uv1 = vertex_texcoord(fi[1], active);
            //     uv2 = vertex_texcoord(fi[2], active);
            // } else {
            uv0 = Point2f(Y * scale_x, Z * scale_y);
            uv1 = Point2f((Y + 1) * scale_x, Z * scale_y);
            uv2 = Point2f(Y * scale_x, (Z + 1) * scale_y);
            // }

            // // if (has_vertex_texcoords() || m_use_default_uv_parameterization) {
            // uv = uv0 * b0 + uv1 * b1 + uv2 * b2;

            // Vector2f duv0 = uv1 - uv0, duv1 = uv2 - uv0;

            // Float det = fmsub(duv0.x(), duv1.y(), duv0.y() * duv1.x()), inv_det = rcp(det);

            // Mask valid = neq(det, 0.f);

            // dp_du = fmsub(duv1.y(), dp0, duv0.y() * dp1) * inv_det;
            // dp_dv = fnmadd(duv1.x(), dp0, duv0.x() * dp1) * inv_det;
            // // }
            // si.uv = Point2f(0,0);

            // Tangents
            si.sh_frame.n = si.n;
            si.dp_du      = m_to_world * dp_du;
            si.dp_dv      = m_to_world * dp_dv;
            si.dp_du      = m_to_world * Vector3f(0, 1, 0);
            si.dp_dv      = m_to_world * Vector3f(0, 0, 1);
            // std::cout << si << std::endl;
        }
    }

    // NOT IMPLEMENTED
    Mask ray_test(const Ray3f &ray_, Mask active) const override {
        MTS_MASK_ARGUMENT(active);
        
//	if(H!=0.0f){
//	  
//	 auto [hit, t] = ray_intersect(ray_);
//	 return hit;
//
//	}
//	else{
	Ray rayLocal = m_to_object.transform_affine(ray_);

        Point3f p0 = m_to_world.transform_affine(Point3f(0, 0, 0)), p1 = m_to_world.transform_affine(Point3f(0, L, 0)),
                p2 = m_to_world.transform_affine(Point3f(0, 0, W)), p3 = m_to_world.transform_affine(Point3f(0, L, W));
        Point2f uv0 = Point2f(0, 0), uv1 = Point2f(1, 0), uv2 = Point2f(0, 1), uv3 = Point2f(1, 1);

        auto [sucess_tri_1, tri1_u, tri1_v, tri1_t] = ray_intersect_triangle(ray_, p0, p2, p3);
        if (sucess_tri_1) {
            return true;
        }

        auto [sucess_tri_2, tri2_u, tri2_v, tri2_t] = ray_intersect_triangle(ray_, p0, p3, p1);
        if (sucess_tri_2) {
            return true;
        }
        return false;

        // Ray3f ray = m_to_object.transform_affine(ray_);

        // Float t       = sphereTrace(ray, active);
        // Point3f local = ray(t);

        // // Is intersection within ray segment and rectangle?
        // return active && t >= ray.mint && t <= ray.maxt && local.z() <= W && local.y() <= L && local.y() >= 0.f && local.z() >= 0.f;
//	}
    }

    void fill_surface_interaction(const Ray3f &ray_, const Float *cache, SurfaceInteraction3f &si_out, Mask active) const override {
        MTS_MASK_ARGUMENT(active);

        // get intersection info
        Float t = cache[0];
        Float Y = cache[1];
        Float Z = cache[2];

        Float cache_u = cache[4];
        Float cache_v = cache[5];
        si_out.t      = t;
        si_out.time   = ray_.time;
        if (H != 0.0f) {
            Point3f local_p = m_to_object.transform_affine(ray_(t));
            Float scale_x   = L * scale_l;
            Float scale_y   = W * scale_w;

            Float x_ = local_p.y() * rcp(scale_x);
            Float y_ = local_p.z() * rcp(scale_y);
            UInt32 Y = (floor(x_));
            UInt32 Z = (floor(y_));

            Float u_proj = x_ - Y; // projected yz coordinate in object space, inside triangle
            Float v_proj = y_ - Z;
        Y = max(min(Y, N - 1), 1);
        Z = max(min(Z, N - 1), 1);
            Float H00 = load<Array<Float, 1>>(height_data + Y + N * Z)[0];
            Float H10 = load<Array<Float, 1>>(height_data + Y + 1 + N * Z)[0];
            Float H01 = load<Array<Float, 1>>(height_data + Y + N * (Z + 1))[0];
            Float H11 = load<Array<Float, 1>>(height_data + (Y + 1) + N * (Z + 1))[0];
            // std::cout << si.uv << std::endl;

            Float v = Y / du_scale, v1 = (Y + 1) / du_scale, w = Z / dv_scale,
                  w1 = (Z + 1) / dv_scale; // voxel coordinates into object local coordinates
            // Point3f p0 = (Point3f((H00), v, w)), p1 = (Point3f((H10), v1, w)), p2 = (Point3f((H01), v, w1)), p3 = (Point3f((H11), v1,
            // w1));
            Point3f p0 = m_to_world.transform_affine(Point3f((H00), v, w)), p1 = m_to_world.transform_affine(Point3f((H10), v1, w)),
                    p2 = m_to_world.transform_affine(Point3f((H01), v, w1)), p3 = m_to_world.transform_affine(Point3f((H11), v1, w1));
            // Point3f p0 = m_to_world.transform_affine(Point3f(0, 0, 0)), p1 = m_to_world.transform_affine(Point3f(0, L, 0)),
            //         p2 = m_to_world.transform_affine(Point3f(0, 0, W)), p3 = m_to_world.transform_affine(Point3f(0, L, W));
            Normal3f n0 = (vertexNormal_2(Y, Z)), n1 = (vertexNormal_2(Y + 1, Z)), n2 = (vertexNormal_2(Y, Z + 1)),
                     n3 = (vertexNormal_2(Y + 1, Z + 1));
            // Normal3f n0 = m_to_world * Normal3f(-1, 0, 0), n1 = m_to_world * Normal3f(-1, 0, 0), n2 = m_to_world * Normal3f(-1, 0, 0),
            //          n3 = m_to_world * Normal3f(-1, 0, 0);

	    Float v0_ = v * rcp(L);
	    Float w0_ = w * rcp(W);
	    Float v1_ = v1 * rcp(L);
	    Float w1_ = w1 * rcp(W);
            Vector2f uv0 = Vector2f(v0_, w0_), uv1 = Vector2f( v1_, w0_),
                    uv2 = Vector2f( v0_, w1_), uv3 = Vector2f( v1_, w1_);

            if (v_proj > u_proj) {
                // if (cache[3]==0.0f) {
                fill_surface_interaction_triangle(cache_u, cache_v, p0, p2, p3, n0, n2, n3, si_out, active, ray_, uv0, uv2, uv3);
            } else {
                fill_surface_interaction_triangle(cache_u, cache_v, p0, p3, p1, n0, n3, n1, si_out, active, ray_,uv0,uv3,uv1);
            }

            // std::cout << si_out << std::endl;
            return;
        } else {

            Point3f p0 = m_to_world.transform_affine(Point3f(0, 0, 0)), p1 = m_to_world.transform_affine(Point3f(0, L, 0)),
                    p2 = m_to_world.transform_affine(Point3f(0, 0, W)), p3 = m_to_world.transform_affine(Point3f(0, L, W));

            Point2f uv0 = Point2f(0, 0), uv1 = Point2f(1, 0), uv2 = Point2f(0, 1), uv3 = Point2f(1, 1);
            Normal3f n0 = Normal3f(-1, 0, 0), n1 = Normal3f(-1, 0, 0), n2 = Normal3f(-1, 0, 0), n3 = Normal3f(-1, 0, 0);

            if (cache[3] == 0.0f) {
                fill_surface_interaction_triangle(cache_u, cache_v, p0, p2, p3, n0, n2, n3, si_out, active, ray_,Vector2f(0),Vector2f(0),Vector2f(0));
            } else if (cache[3] == 1.0f) {
                fill_surface_interaction_triangle(cache_u, cache_v, p0, p3, p1, n0, n3, n1, si_out, active, ray_,Vector2f(0),Vector2f(0),Vector2f(0));
            }

            // std::cout << si_out << std::endl;
            return;
        }
    }

    SurfaceInteraction3f eval_surfaceInteraction_from_uv(Point2f p, Mask active, SurfaceInteraction3f &si_) const {
        MTS_MASK_ARGUMENT(active);
        SurfaceInteraction3f si = si_;
//////////////////////////////////////////////////////////////////////////////////
        Float scale_x = L * scale_l;
        Float scale_y = W * scale_w;

        Float x_ = p.x() * L * rcp(scale_x);
        Float y_ = p.y() * W * rcp(scale_y);
        UInt32 Y = (floor(x_));
        UInt32 Z = (floor(y_));

        Float u_proj = x_ - Y; // projected yz coordinate in object space, inside triangle
        Float v_proj = y_ - Z;

        // triangle contruction
        Y = max(min(Y, N - 1), 1);
        Z = max(min(Z, N - 1), 1);
        Float H00 = load<Array<Float, 1>>(height_data + Y + N * Z)[0];
        Float H10 = load<Array<Float, 1>>(height_data + Y + 1 + N * Z)[0];
        Float H01 = load<Array<Float, 1>>(height_data + Y + N * (Z + 1))[0];
        Float H11 = load<Array<Float, 1>>(height_data + (Y + 1) + N * (Z + 1))[0];

        Float v_ = Y / du_scale, v1 = (Y + 1) / du_scale, w_ = Z / dv_scale,
              w1   = (Z + 1) / dv_scale; // voxel coordinates into object local coordinates
        Point3f P0 = m_to_world.transform_affine(Point3f((H00), v_, w_)), P1 = m_to_world.transform_affine(Point3f((H10), v1, w_)),
                P2 = m_to_world.transform_affine(Point3f((H01), v_, w1)), P3 = m_to_world.transform_affine(Point3f((H11), v1, w1));

        Normal3f N0 = vertexNormal_2(Y, Z), N1 = vertexNormal_2(Y + 1, Z), N2 = vertexNormal_2(Y, Z + 1), N3 = vertexNormal_2(Y + 1, Z + 1);

	    Float v0_ = v_ * rcp(L);
	    Float w0_ = w_ * rcp(W);
	    Float v1_ = v1 * rcp(L);
	    Float w1_ = w1 * rcp(W);
            Vector2f UV0 = Vector2f(v0_, w0_), UV1 = Vector2f( v1_, w0_),
                    UV2 = Vector2f( v0_, w1_), UV3 = Vector2f( v1_, w1_);

        Point3f p0, p1, p2;
        Normal3f n0, n1, n2;
	Point2f uv0, uv1, uv2;
        if (v_proj > u_proj) { // upper triangle
            p0 = P0, p1 = P2, p2 = P3;
            n0 = N0, n1 = N2, n2 = N3;
	    uv0 = UV0,uv1=UV2,uv2=UV3;
        } else { // lower triangle
            p0 = P0, p1 = P3, p2 = P1;
            n0 = N0, n1 = N3, n2 = N1;
	    uv0 = UV0,uv1=UV3,uv2=UV1;
        }

        // postion sample contruction
        Float b1 = u_proj, b2 = v_proj;

        Float b0 = 1.f - b1 - b2;

        Vector3f dp0 = p1 - p0, dp1 = p2 - p0;

        // Re-interpolate intersection using barycentric coordinates
        si.p = (p0 * b0 + p1 * b1 + p2 * b2); // H being quite small, this is considered
                                              // uniform sampling

        // ps.p =  m_to_world.transform_affine(Point3f(0, 0.5, 0.5));
        // std::cout << ps.p  << "H = " << H << std::endl;
        si.n = m_to_world * (-flip * normalize(n0 * b0 + n1 * b1 + n2 * b2));
	si.sh_frame.n = si.n;
        // std::cout << H << std::endl;

            auto [dp_du, dp_dv] = coordinate_system(si.n);
            Point2f uv(0, 0);

            uv = uv0 * b0 + uv1 * b1 + uv2 * b2;

            Vector2f duv0 = uv1 - uv0, duv1 = uv2 - uv0;

            Float det = fmsub(duv0.x(), duv1.y(), duv0.y() * duv1.x()), inv_det = rcp(det);

            Mask valid = neq(det, 0.f);

            dp_du = fmsub(duv1.y(), dp0, duv0.y() * dp1) * inv_det;
            dp_dv = fnmadd(duv1.x(), dp0, duv0.x() * dp1) * inv_det;

            si.uv = uv;

            // Tangents
            si.dp_du = dp_du;
            si.dp_dv = dp_dv;
	    //std::cout << si  << std::endl;
	    return si;
    }

    std::pair<Vector3f, Vector3f> normal_derivative(const SurfaceInteraction3f &si, bool shading_frame, Mask active) const override {
        MTS_MASK_ARGUMENT(active);

        // if (!shading_frame)
        //     return { zero<Vector3f>(), zero<Vector3f>() };
        if (H != 0) {
            Point3f p_local = m_to_object.transform_affine(si.p);

            Float scale_x = L * scale_l;
            Float scale_y = W * scale_w;

            Float x_ = p_local.y() * rcp(scale_x);
            Float y_ = p_local.z() * rcp(scale_y);
            UInt32 Y = (floor(x_));
            UInt32 Z = (floor(y_));

            Float u_proj = x_ - Y; // projected yz coordinate in object space, inside triangle
            Float v_proj = y_ - Z;

            // triangle contruction
        Y = max(min(Y, N - 1), 1);
        Z = max(min(Z, N - 1), 1);
            Float H00 = load<Array<Float, 1>>(height_data + Y + N * Z)[0];
            Float H10 = load<Array<Float, 1>>(height_data + Y + 1 + N * Z)[0];
            Float H01 = load<Array<Float, 1>>(height_data + Y + N * (Z + 1))[0];
            Float H11 = load<Array<Float, 1>>(height_data + (Y + 1) + N * (Z + 1))[0];

            Float v0 = Y / du_scale, v1 = (Y + 1) / du_scale, w0 = Z / dv_scale,
                  w1   = (Z + 1) / dv_scale; // voxel coordinates into object local coordinates
            Point3f P0 = m_to_world.transform_affine(Point3f((H00), v0, w0)), P1 = m_to_world.transform_affine(Point3f((H10), v1, w0)),
                    P2 = m_to_world.transform_affine(Point3f((H01), v0, w1)), P3 = m_to_world.transform_affine(Point3f((H11), v1, w1));
            Normal3f N0 = vertexNormal_2(Y, Z), N1 = vertexNormal_2(Y + 1, Z), N2 = vertexNormal_2(Y, Z + 1),
                     N3 = vertexNormal_2(Y + 1, Z + 1);

	    Float v0_ = v0 * rcp(L);
	    Float w0_ = w0 * rcp(W);
	    Float v1_ = v1 * rcp(L);
	    Float w1_ = w1 * rcp(W);
            Vector2f UV0 = Vector2f(v0_, w0_), UV1 = Vector2f( v1_, w0_),
                    UV2 = Vector2f( v0_, w1_), UV3 = Vector2f( v1_, w1_);
            Point3f p0, p1, p2;
            Normal3f n0, n1, n2;
            Point2f uv0, uv1, uv2;
            if (v_proj > u_proj) { // upper triangle
                p0 = P0, p1 = P2, p2 = P3;
                n0 = N0, n1 = N2, n2 = N3;
		uv0= UV0,uv1=UV2,uv2=UV3;
            } else { // lower triangle
                p0 = P0, p1 = P3, p2 = P1;
                n0 = N0, n1 = N3, n2 = N1;
		uv0= UV0,uv1=UV3,uv2=UV1;
            }

            Vector3f rel = si.p - p0, du = p1 - p0, dv = p2 - p0;

            /* Solve a least squares problem to determine
               the UV coordinates within the current triangle */
            Float b1 = dot(du, rel), b2 = dot(dv, rel), a11 = dot(du, du), a12 = dot(du, dv), a22 = dot(dv, dv),
                  inv_det = rcp(a11 * a22 - a12 * a12);

            Float u = fmsub(a22, b1, a12 * b2) * inv_det, v = fnmadd(a12, b1, a11 * b2) * inv_det, w = 1.f - u - v;

            /* Now compute the derivative of "normalize(u*n1 + v*n2 + (1-u-v)*n0)"
               with respect to [u, v] in the local triangle parameterization.

               Since d/du [f(u)/|f(u)|] = [d/du f(u)]/|f(u)|
                 - f(u)/|f(u)|^3 <f(u), d/du f(u)>, this results in
            */
            // Normal3f N = m_to_world*(-flip*normalize(u * n1 + v * n2 + w * n0));
            // Normal3f N(u * n0 + v * n1 + w * n2);
            Normal3f N(u * n1 + v * n2 + w * n0);
            Float il = rsqrt(squared_norm(N));
            N *= il;

            Vector3f dndu = (n1 - n0) * il;
            Vector3f dndv = (n2 - n0) * il;

            dndu = fnmadd(N, dot(N, dndu), dndu);
            dndv = fnmadd(N, dot(N, dndv), dndv);

            // if (has_vertex_texcoords()) {
            //     uv0 = vertex_texcoord(fi[0], active);
            //     uv1 = vertex_texcoord(fi[1], active);
            //     uv2 = vertex_texcoord(fi[2], active);
            // } else {
            // }

            // if (has_vertex_texcoords() || m_use_default_uv_parameterization) {
            Vector2f duv1 = uv1 - uv0, duv2 = uv2 - uv0;
            inv_det = rcp(duv1.x() * duv2.y() - duv1.y() * duv2.x());

            Vector3f dndu_ = (duv2.y() * dndu - duv1.y() * dndv) * inv_det, dndv_ = (-duv2.x() * dndu + duv1.x() * dndv) * inv_det;
            dndu =dndu_;
            dndv =dndv_;
            // }
            //

            // return { Vector3f(0), Vector3f(0) };
            // std::cout << si << std::endl;
            // std::cout << dndu << dndv << std::endl;
	    if(!this->is_caustic_caster_double_refraction()){
	      std::cout << "bug" << std::endl;
	    }
            return { m_to_world * dndu, m_to_world * dndv };
        } else {
            return { Vector3f(0), Vector3f(0) };
        }
    }

    //! @}
    // =============================================================

    ScalarSize primitive_count() const override { return 1; }

    ScalarSize effective_primitive_count() const override { return 1; }

    void traverse(TraversalCallback *callback) override { Base::traverse(callback); }

    void parameters_changed(const std::vector<std::string> & /*keys*/) override {
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

            OptixHeightfieldData data = { bbox(), m_to_world, m_to_object, m_center, m_radius, m_flip_normals };

            cuda_memcpy_to_device(m_optix_data_ptr, &data, sizeof(OptixHeightfieldData));
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
    ScalarFloat scale_l, scale_w; // scale factor to normalize local coordinates ([0,L] --> [0,1])
    ScalarFloat invL, invW;
    Float du_scale, dv_scale, dh_scale;

    std::string reconstructionMode;
    // heightmap
    ref<Texture> heightmap;
    Array<Float, 1> *height_data;    // height samples
    Array<Normal3f, 1> *normal_data; // vertex normals
    Float partition_size, inv_partition_size, inv_partition_area;
    mutable long nbEvall;
};

MTS_IMPLEMENT_CLASS_VARIANT(Heightfield, Shape)
// MTS_EXTERN_CLASS_RENDER(Heightfield)
MTS_EXPORT_PLUGIN(Heightfield, "Heightfield intersection primitive");
// MTS_INSTANTIATE_CLASS(Heightfield)
NAMESPACE_END(mitsuba)
