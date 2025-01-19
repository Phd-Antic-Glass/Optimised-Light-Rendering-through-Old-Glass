#include <mitsuba/render/fermatNEE.h>
#include <mitsuba/render/scene.h>
#include <mitsuba/python/python.h>

MTS_PY_EXPORT(FermatNEE) {
    MTS_PY_IMPORT_TYPES()
    using FermatNEE = FermatNEE<Float, Spectrum>;
    using Matrix2d = Matrix<double, 2>;

    MTS_PY_STRUCT(FermatNEE)
        .def(py::init<>())
        .def("init", &FermatNEE::init,
            "scene"_a, "sampler"_a,  "solutionIdentical_threshold"_a, "maxBernouilliTrial"_a, "alpha1"_a, "beta"_a, "multi_piece_mode"_a, "m_sms_config"_a, "use_SMS"_a)
        .def("specular_connection_debug", &FermatNEE::specular_connection_debug,
             "si"_a, "vy"_a, "init"_a=Vector3f(0));
        // .def("sample_emitter", &FermatNEE::sample_emitter, "si"_a, "active"_a=true);
	//
}
