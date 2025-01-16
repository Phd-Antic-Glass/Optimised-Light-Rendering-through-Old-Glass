#include <mitsuba/render/fermatNEE.h>
#include <mitsuba/render/scene.h>
#include <mitsuba/python/python.h>

MTS_PY_EXPORT(fermatNEE) {
    MTS_PY_IMPORT_TYPES()
    using fermatNEE = fermatNEE<Float, Spectrum>;
    using Matrix2d = Matrix<double, 2>;

    MTS_PY_STRUCT(fermatNEE)
        .def(py::init<>())
        .def("init", &fermatNEE::init,
            "scene"_a, "sampler"_a,  "solutionIdentical_threshold"_a, "maxBernouilliTrial"_a, "alpha1"_a, "beta"_a, "multi_piece_mode"_a, "m_sms_config"_a)
        .def("fermat_connection_", &fermatNEE::fermat_connection_,
             "si"_a, "vy"_a, "init"_a=Point2f(-1,-1),"shape"_a=nullptr, "H_out"_a=0);
        // .def("sample_emitter", &fermatNEE::sample_emitter, "si"_a, "active"_a=true);
	//
}
