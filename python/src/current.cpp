// SPDX-License-Identifier: LGPL-3.0-or-other
// Copyright (C) 2021 Contributors to the SLS Detector Package
#include <pybind11/chrono.h>
#include <pybind11/numpy.h>
#include <pybind11/operators.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

// #include "sls/Pattern.h"
#include "sls/ToString.h"
#include "sls/sls_detector_defs.h"

namespace py = pybind11;
void init_source(py::module &m) {

    using src = slsDetectorDefs::currentSrcParameters;
    py::class_<src> currentSrcParameters(m, "currentSrcParameters");

    currentSrcParameters.def(py::init());
    currentSrcParameters.def_readwrite("enable", &src::enable);
    currentSrcParameters.def_readwrite("fix", &src::fix);
    currentSrcParameters.def_readwrite("normal", &src::normal);
    currentSrcParameters.def_readwrite("select", &src::select);
    currentSrcParameters.def(pybind11::self == pybind11::self);

    currentSrcParameters.def("__repr__",
                             [](const src &a) { return sls::ToString(a); });
}
