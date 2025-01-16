#include <mitsuba/render/OpenGL_viewer_client.h>

#include <cstddef>
#include <enoki/stl.h>
#include <mitsuba/core/plugin.h>
#include <mitsuba/core/properties.h>
#include <mitsuba/render/bsdf.h>
#include <mitsuba/render/fermatNEE.h>
#include <mitsuba/render/integrator.h>
#include <mitsuba/render/kdtree.h>
#include <mitsuba/render/medium.h>
#include <mitsuba/render/scene.h>

#if defined(MTS_ENABLE_EMBREE)
#include "scene_embree.inl"
#else
#include "scene_native.inl"
#endif

NAMESPACE_BEGIN(mitsuba)

template <typename Float, typename Spectrum>
static int OpenGL_viewer_client<Float, Spectrum>::OpenGL_command_drawline(Point3f start, Point3f end, Point3f color) {
    try {
        Connection conn("127.0.0.1", 8080);

        std::string msg = "drawline:";
        msg             = msg + std::to_string(start.x()) + ",";
        msg             = msg + std::to_string(start.y()) + ",";
        msg             = msg + std::to_string(start.z()) + ";";

        msg = msg + std::to_string(end.x()) + ",";
        msg = msg + std::to_string(end.y()) + ",";
        msg = msg + std::to_string(end.z()) + ";";

        msg = msg + std::to_string(color.x()) + ",";
        msg = msg + std::to_string(color.y()) + ",";
        msg = msg + std::to_string(color.z()) + ";";

        conn.tx(msg);

        // cout << msg << endl;
        // cout << i << endl;
        close(conn.getSocket());

    } catch (std::exception &e) {
        std::cerr << "OpenGL_viewer_client: " << e.what() << std::endl;
        return EXIT_FAILURE;
    }
    return 0;
}

MTS_INSTANTIATE_CLASS(OpenGL_viewer_client)
NAMESPACE_END(mitsuba)