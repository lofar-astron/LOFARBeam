#include "OskarSphericalWaveCoeff.h"

H5::CompType get_complex_double_type()
{
    H5::CompType complex_type(sizeof(std::complex<double>));
    complex_type.insertMember("real", 0, H5::PredType::NATIVE_DOUBLE);
    complex_type.insertMember("imag", sizeof(double), H5::PredType::NATIVE_DOUBLE);
    return complex_type;
}

// Constructor for reading coeff from file
OskarSphericalWaveCoefficients::OskarSphericalWaveCoefficients(
    std::string& filename)
{
    read_coeffs(filename);
};

size_t OskarSphericalWaveCoefficients::get_l_max() const
{
    return 0;
}

void OskarSphericalWaveCoefficients::read_coeffs(
    std::string& filename)
{
    // Open file
    m_h5_file.reset(new H5::H5File(filename, H5F_ACC_RDONLY));

    // Disable HDF5 error prints
    H5::Exception::dontPrint();
}

bool OskarSphericalWaveCoefficients::read_frequency(
    const unsigned int freq)
{
    if (m_freq != freq) {
        std::string dataset_name = std::to_string((int) (freq / 1e6));
        std::cout << "read oskar dataset: " << dataset_name << std::endl;
        try {
            H5::DataSet dataset = m_h5_file->openDataSet(dataset_name);
            m_dataset_available = true;
        } catch (H5::FileIException& e) {
            std::cerr << "Could not load dataset for frequency " << dataset_name << " Mhz" << std::endl;
            m_dataset_available = false;
        }
    }

    m_freq = freq;
    return m_dataset_available;
}

#if 0
    // Read dataset
    H5::DataSet dataset = file.openDataSet(m_dataset_name);

    // Open attribute and read its contents
    H5::Attribute freq_center_attr = dataset.openAttribute("freq_center");
    H5::Attribute freq_range_attr = dataset.openAttribute("freq_range");
    freq_center_attr.read(H5::PredType::NATIVE_DOUBLE, &m_freq_center);
    freq_range_attr.read(H5::PredType::NATIVE_DOUBLE, &m_freq_range);

    // Read dataset dimensions
    H5::DataSpace dataspace = dataset.getSpace();
    unsigned int rank = dataspace.getSimpleExtentNdims();
    assert(rank == m_dataset_rank);
    hsize_t dims[rank];
    dataspace.getSimpleExtentDims(dims, NULL);
    m_nHarmonics  = dims[0];
    m_nPowerTheta = dims[1];
    m_nPowerFreq  = dims[2];

    // Read coeff
    std::complex<double> coeff[m_nHarmonics][m_nPowerTheta][m_nPowerFreq][m_nInner];
    H5::CompType complex_type = get_complex_double_type();
    dataset.read(coeff, complex_type, dataspace);

    // Set m_coeff
    m_coeff.resize(get_nr_coeffs());
    set_coeffs((std::complex<double> *) coeff);
#endif


void OskarSphericalWaveCoefficients::print_coeffs()
{
    const int l_max = m_l_max;

    for (int l = 1; l <= l_max; ++l) {
        for (int abs_m = l; abs_m >=0; --abs_m) {
            // TODO
        }
    }
}