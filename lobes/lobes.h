#include <Eigen/Core>
#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>

namespace py = pybind11;


class BeamModel {

public:

    BeamModel() {}

    virtual py::array_t<std::complex<double>> eval(py::EigenDRef<const Eigen::ArrayXd> theta, py::EigenDRef<const Eigen::ArrayXd> phi)=0;
};


class LobesBeamModel : public BeamModel {

public:

    LobesBeamModel(const std::string &data_file_name);

    py::array_t<std::complex<double>> eval(py::EigenDRef<const Eigen::ArrayXd> theta, py::EigenDRef<const Eigen::ArrayXd> phi) override;

private:
    Eigen::ArrayXXcd m_coefficients;
    std::vector<size_t> m_coefficients_shape;

    std::vector<double> m_frequencies;
    py::array_t<int> m_nms;
};

class HamakerBeamModel : public BeamModel {

public:

    HamakerBeamModel() {}

    py::array_t<std::complex<double>> eval(py::EigenDRef<const Eigen::ArrayXd> theta, py::EigenDRef<const Eigen::ArrayXd> phi) override {}
};
