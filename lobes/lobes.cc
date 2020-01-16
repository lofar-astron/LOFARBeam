#define __STDCPP_WANT_MATH_SPEC_FUNCS__ 1
#include <cmath>

#include <iostream>
#include <complex>

#include "ElementResponse.h"

#include <H5Cpp.h>

#include "lobes.h"

using namespace std::complex_literals;


#include <tuple>

// python style enumerate function
// To make it possible to write:
//   for (auto [i, v] : enumerate(iterable)) {...}
template <typename T,
          typename TIter = decltype(std::begin(std::declval<T>())),
          typename = decltype(std::end(std::declval<T>()))>
constexpr auto enumerate(T && iterable)
{
    struct iterator
    {
        size_t i;
        TIter iter;
        bool operator != (const iterator & other) const { return iter != other.iter; }
        void operator ++ () { ++i; ++iter; }
        auto operator * () const { return std::tie(i, *iter); }
    };
    struct iterable_wrapper
    {
        T iterable;
        auto begin() { return iterator{ 0, std::begin(iterable) }; }
        auto end() { return iterator{ 0, std::end(iterable) }; }
    };
    return iterable_wrapper{ std::forward<T>(iterable) };
}


double P20(double x) { return 0.5*(3*x*x-1); }
double P21(double x) { return 3.0*x*std::sqrt(1-x*x); }
double P22(double x) { return 3*(1-x*x); }

int plustwo(int a)
{
    std::cout << std::assoc_legendre(2, 0, 0.5) << '=' << P20(0.5) << '\n'
              << std::assoc_legendre(2, 1, 0.5) << '=' << P21(0.5) << '\n'
              << std::assoc_legendre(2, 2, 0.5) << '=' << P22(0.5) << '\n';

    return a+2;
}


Eigen::MatrixXd big_mat()
{
    return Eigen::MatrixXd::Zero(10000, 10000);
}


Eigen::ArrayXd P(int m, int n, py::EigenDRef<Eigen::ArrayXd> x)
{
    auto N = x.rows();

    Eigen::ArrayXd result(N); // = Eigen::VectorXd::Zero(10, 1);
    for(int i = 0; i<N; i++)
    {
        result[i] = std::assoc_legendre(n, std::abs(m), x[i]);
    }

    if (m<0)
    {
        result *= std::pow(-1, -m) * std::tgamma(n+m+1) / std::tgamma(n-m+1);
    }

    return std::move(result);
}

Eigen::ArrayXd Pacc(int m, int n, py::EigenDRef<Eigen::ArrayXd> x)
{
    auto N = x.rows();

    return (-(n+m) * (n-m+1.0) * Eigen::sqrt(1.0 - x*x) * P(m-1, n, x) - m * x * P(m, n, x))/(x*x-1.0);
}


/** Compute Spherical Wave Function
  *
  * returns std::pair(Etheta, Ephi)
  */
std::pair<Eigen::VectorXcd,Eigen::VectorXcd> F4far_new(int s, int m, int n,  py::EigenDRef<const Eigen::ArrayXd> theta, py::EigenDRef<const Eigen::ArrayXd> phi, double beta)
{
    int N = theta.rows();
    Eigen::VectorXd result(N);

    double C;
    if (m) {
        C = beta * std::sqrt(60.0) * 1.0/std::sqrt(n*(n+1.0)) * std::pow(-m/std::abs(m), m);
    }
    else {
        C = beta * std::sqrt(60.0) * 1.0/std::sqrt(n * (n+1.0));
    }

    Eigen::ArrayXcd q2;
    Eigen::ArrayXcd q3;


    Eigen::ArrayXd cos_theta = Eigen::cos(theta);
    Eigen::ArrayXd sin_theta = Eigen::sin(theta);
    Eigen::ArrayXd P_cos_theta = P(std::abs(m),n,cos_theta);
    Eigen::ArrayXd Pacc_cos_theta = Pacc(std::abs(m),n,cos_theta);
    Eigen::ArrayXcd exp_i_m_phi = Eigen::exp(1.0i*double(m)*phi);

    if (s == 1) {

        q2 = C * std::pow(-1.0i, -n-1)/beta * 1.0i * double(m) / (sin_theta) *
             std::sqrt((2.*n+1)/2.0 * std::tgamma(n-std::abs(m)+1) / std::tgamma(n+std::abs(m)+1)) *
             P_cos_theta * exp_i_m_phi;

        q3 = C * std::pow(-1.0i, -n-1)/beta * std::sqrt((2.*n+1) / 2.0 * std::tgamma(n-abs(m)+1) / std::tgamma(n+abs(m)+1)) *
             Pacc_cos_theta * sin_theta * exp_i_m_phi;
    }
    else if (s==2) {
        q2 = -C * std::pow(-1.0i, -n) / beta * std::sqrt((2.*n+1)/2.0 * std::tgamma(n-abs(m)+1) / std::tgamma(n+abs(m)+1)) *
             Pacc_cos_theta * sin_theta * exp_i_m_phi;

        q3 = C * std::pow(-1.0i, -n) / beta * 1.0i * double(m)/sin_theta * std::sqrt((2.*n+1) / 2.0 * std::tgamma(n-abs(m)+1) / std::tgamma(n+abs(m)+1)) *
             P_cos_theta * exp_i_m_phi;
    }

    return std::make_pair(q2,q3);
}

Eigen::Array<std::complex<double>, 2, 2> element_response_lba(double freq, double theta, double phi)
{
    Eigen::Array<std::complex<double>, 2, 2> response;

    typedef std::complex<double> Response[2][2];
    Response &response_ = *((Response *) response.data());

    LOFAR::element_response_lba(freq, theta, phi, response_);
    return std::move(response);
}


LobesBeamModel::LobesBeamModel(const std::string &data_file_name)
{
    H5::H5File h5file( data_file_name.c_str(), H5F_ACC_RDONLY );

    H5::DataSet dataset = h5file.openDataSet( "coefficients" );
    H5::DataSpace dataspace = dataset.getSpace();
    int nr_elements = dataspace.getSimpleExtentNpoints();

    const std::string REAL( "r" );
    const std::string IMAG( "i" );

    H5::CompType h5_dcomplex( sizeof(std::complex<double>) );
    h5_dcomplex.insertMember( REAL, 0, H5::PredType::NATIVE_DOUBLE);
    h5_dcomplex.insertMember( IMAG, sizeof(double), H5::PredType::NATIVE_DOUBLE);

    /*
        * Get the number of dimensions in the dataspace.
        */
    int ndims_coefficients = dataspace.getSimpleExtentNdims();
    /*
        * Get the dimension size of each dimension in the dataspace and
        * display them.
        */
    hsize_t dims_coefficients[ndims_coefficients];
    dataspace.getSimpleExtentDims( dims_coefficients, NULL);

    std::cout << "Dimensions of coefficients:" << std::endl << "    ";
    for(int i = 0; i<ndims_coefficients; i++)
    {
        std::cout << (size_t)(dims_coefficients[i]);
        if (i < ndims_coefficients-1)
            std::cout  << " x ";
        else
            std::cout << std::endl;
    }

    // Reshape coefficients
    size_t nr_elements_first_dims = 1;
    for(int i = 0; i<ndims_coefficients-1; i++) nr_elements_first_dims *= dims_coefficients[i];
    std::cout << nr_elements_first_dims << " x " << (size_t)(dims_coefficients[ndims_coefficients-1]) << std::endl;

    m_coefficients = Eigen::ArrayXXcd(nr_elements_first_dims, (size_t)dims_coefficients[ndims_coefficients-1]);
    m_coefficients_shape = std::vector<size_t>(dims_coefficients, dims_coefficients + ndims_coefficients);
    dataset.read(m_coefficients.data(), h5_dcomplex);

//=====================================
    dataset = h5file.openDataSet( "frequencies" );
    dataspace = dataset.getSpace();
    nr_elements = dataspace.getSimpleExtentNpoints();

    m_frequencies = std::vector<double> (nr_elements);

    dataset.read(m_frequencies.data(), H5::PredType::NATIVE_DOUBLE);
    std::cout << "Frequencies:" << std::endl;
    for(auto freq : m_frequencies) {
        std::cout << freq << " ";
    }
    std::cout << std::endl << std::endl;


    dataset = h5file.openDataSet( "nms" );
    dataspace = dataset.getSpace();
    nr_elements = dataspace.getSimpleExtentNpoints();


    /*
        * Get the number of dimensions in the dataspace.
        */
    int ndims_nms = dataspace.getSimpleExtentNdims();
    /*
        * Get the dimension size of each dimension in the dataspace and
        * display them.
        */
    hsize_t dims_nms[ndims_nms];
    dataspace.getSimpleExtentDims( dims_nms, NULL);

    std::cout << "Dimensions of nms:" << ndims_nms << std::endl << "    ";

    for(int i = 0; i<ndims_nms; i++)
    {
        std::cout << (size_t)(dims_nms[i]);
        if (i < ndims_nms-1)
                std::cout  << " x ";
    }
    std::cout << std::endl;

    m_nms = py::array_t<int>({dims_nms[0], dims_nms[1]});
    dataset.read(m_nms.mutable_data(), H5::PredType::NATIVE_INT);
}

py::array_t<std::complex<double>> LobesBeamModel::eval(py::EigenDRef<const Eigen::ArrayXd> theta, py::EigenDRef<const Eigen::ArrayXd> phi)
{

    double beta;
    for (auto [i, freq] : enumerate(m_frequencies)) {
        beta = 2.0 * M_PI * freq / 2.99792458e8;
    }


    beta = 2.0 * M_PI * m_frequencies[0] / 2.99792458e8;

    int s,m,n;

    auto nms = m_nms.unchecked<2>();

//     auto result = py::array_t<std::complex<double>>({2*theta.size(), m_nms.shape(0)});
//     Eigen::Map<Eigen::MatrixXcd> Base(result.mutable_data(), 2*theta.size(), m_nms.shape(0));

    Eigen::MatrixXcd Base(2*theta.size(), m_nms.shape(0));

    for(int i = 0; i < m_nms.shape(0); i++) {
        int n = nms(i,0);
        int m = nms(i,1);
        int s = nms(i,2);
//         std::cout << i << ": (" << n << ", " << m << ", " << s << ")" << std::endl;
        Eigen::VectorXcd e_phi, e_theta;
        // Compute SWF
        std::tie(e_theta, e_phi) = F4far_new(s, m, n,  theta, phi, beta);
        // Interleave
        Base.col(i)(Eigen::seq(0,Eigen::last,2)) = e_theta;
        Base.col(i)(Eigen::seq(1,Eigen::last,2)) = e_phi;
    }

    std::cout << Base.rows() << ", " << m_coefficients.rows() << std::endl;
    std::cout << Base.cols() << ", " << m_coefficients.cols() << std::endl;

    // Create array to return
    auto result = py::array_t<std::complex<double>>({Base.rows(), m_coefficients.rows()});

    // Map the result array to an Eigen array, to be able to assign it the result of an Eigen matrix multiplication
    Eigen::Map<Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>>
        result_matrix(result.mutable_data(), Base.rows(), m_coefficients.rows());

    // Multiply base functions by coefficients
    result_matrix = Base *  m_coefficients.matrix().transpose();

    return result;

}


// py::array_t<std::complex<double>> LobesBeamModel::eval(py::array_t<double> theta, py::array_t<double> phi)
// {
//     py::buffer_info theta_buffer = theta.request();
//     py::buffer_info phi_buffer = phi.request();
//
//     if (theta_buffer.ndim != 1 || phi_buffer.ndim != 1)
//         throw std::runtime_error("Number of dimensions must be one");
//
//     if (theta_buffer.size != phi_buffer.size)
//         throw std::runtime_error("Input shapes must match");
//
//     /* No pointer is passed, so NumPy will allocate the buffer */
//     int size = theta_buffer.size;
//     auto result = py::array_t<std::complex<double>>({size, 2, 2});
//
//     double beta = 1.0;
//
//     Eigen::Map<Eigen::ArrayXd> theta_((double *) theta_buffer.ptr, size);
//     Eigen::Map<Eigen::ArrayXd> phi_(size, (double *) phi_buffer.ptr);
//
// //     std::pair<Eigen::VectorXcd,Eigen::VectorXcd> F4far_new(s, m, n,  theta, phi, beta)
//
//     int s,m,n;
//
// //     F4far_new(s, m, n,  theta_, phi_, beta);
//
// //     py::buffer_info buf3 = result.request();
// //
// //     double *ptr1 = (double *) buf1.ptr,
// //            *ptr2 = (double *) buf2.ptr,
// //            *ptr3 = (double *) buf3.ptr;
// //
// //     for (size_t idx = 0; idx < buf1.shape[0]; idx++)
// //         ptr3[idx] = ptr1[idx] + ptr2[idx];
//
//     return result;
//
// }


PYBIND11_MODULE(lobes, m) {

    m.doc() = "LOBES module"; // optional module docstring

    m.def("plustwo", &plustwo, "A function which adds two to a number");

    m.def("assoc_legendre", &std::assoc_legendre<double>, "Associated legendre function");

    m.def("scale", [](py::EigenDRef<Eigen::MatrixXd> m, double c) { m *= c; });

    m.def("big_mat", &big_mat, "Big matrix");

    m.def("P", &P, "Legendre");
    m.def("Pacc", &Pacc, "Legendre");
    m.def("F4far_new", &F4far_new, "F4far_new");
    m.def("element_response_lba", &element_response_lba, "element_response_lba");

    py::class_<LobesBeamModel>(m, "LobesBeamModel")
        .def(py::init<const std::string &>())
        .def("eval", &LobesBeamModel::eval);

}
