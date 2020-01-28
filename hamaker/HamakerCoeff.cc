#include "HamakerCoeff.h"

H5::CompType get_complex_double_type()
{
    H5::CompType complex_type(sizeof(std::complex<double>));
    complex_type.insertMember("r", 0, H5::PredType::NATIVE_DOUBLE);
    complex_type.insertMember("i", sizeof(double), H5::PredType::NATIVE_DOUBLE);
    return complex_type;
}

size_t HamakerCoefficients::get_index(
    const unsigned int n,
    const unsigned int t,
    const unsigned int f)
{
    return n * m_nPowerTheta * m_nPowerFreq * m_nInner +
                           t * m_nPowerFreq * m_nInner +
                                          f * m_nInner;
}

// Constructor for reading coeff from file
HamakerCoefficients::HamakerCoefficients(
    std::string& filename)
{
    read_coeffs(filename);
};

// Constructor for writing coeff to file
HamakerCoefficients::HamakerCoefficients(
    const double freq_center,
    const double freq_range,
    const unsigned int nHarmonics,
    const unsigned int nPowerTheta,
    const unsigned int nPowerFreq) :
    m_freq_center(freq_center),
    m_freq_range(freq_range),
    m_nHarmonics(nHarmonics),
    m_nPowerTheta(nPowerTheta),
    m_nPowerFreq(nPowerFreq),
    m_coeff(get_nr_coeffs())
{}

size_t HamakerCoefficients::get_nr_coeffs() const
{
    return m_nHarmonics * m_nPowerTheta * m_nPowerFreq * m_nInner;
}

void HamakerCoefficients::set_coeff(
    const unsigned int n,
    const unsigned int t,
    const unsigned int f,
    std::pair<std::complex<double>, std::complex<double>> value)
{
    size_t index = get_index(n, t, f);
    m_coeff[index] = value.first;
    m_coeff[index + 1] = value.second;
}

void HamakerCoefficients::set_coeffs(
    const std::complex<double>* coeff) {
    memcpy(m_coeff.data(), coeff, m_coeff.size() * sizeof(std::complex<double>));
}

void HamakerCoefficients::set_coeffs(
    const std::vector<std::complex<double>> coeff)
{
    assert(coeff.size() == m_coeff.size());
    std::copy(coeff.begin(), coeff.end(), m_coeff.begin());
}

std::pair<std::complex<double>, std::complex<double>> HamakerCoefficients::get_coeff(
    const unsigned int n,
    const unsigned int t,
    const unsigned int f)
{
    size_t index = get_index(n, t, f);
    std::pair<std::complex<double>, std::complex<double>> value;
    value.first = m_coeff[index];
    value.second = m_coeff[index + 1];
    return value;
}

void HamakerCoefficients::read_coeffs(
    std::string& filename)
{
    // Open file
    H5::H5File file(filename, H5F_ACC_RDONLY);

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

    // Read coeffs
    m_coeff.resize(get_nr_coeffs());
    H5::DataType data_type = dataset.getDataType();
    assert(data_type.getSize() == sizeof(std::complex<double>));
    dataset.read(m_coeff.data(), data_type, dataspace);
}

void HamakerCoefficients::write_coeffs(
    std::string& filename)
{
    // Open file
    H5::H5File file(filename, H5F_ACC_TRUNC);

    // Create dataspace
    const unsigned int rank = 4;
    hsize_t dims[rank] = { m_nHarmonics, m_nPowerTheta, m_nPowerFreq, m_nInner };
    H5::DataSpace coeff_dataspace(rank, dims);

    // Create pointer to coeff
    typedef std::complex<double> coeff_type[m_nHarmonics][m_nPowerTheta][m_nPowerFreq][m_nInner];
    coeff_type *coeff_ptr = (coeff_type *) m_coeff.data();

    // Create complex type
    H5::CompType complex_type = get_complex_double_type();

    // Write dataset
    H5::DataSet dataset = file.createDataSet("coeff", complex_type, coeff_dataspace);
    dataset.write(coeff_ptr, complex_type);

    // Create new dataspace for attribute
    H5::DataSpace attr_dataspace(H5S_SCALAR);

    // Write frequency center attribute
    H5::Attribute freq_center_attr = dataset.createAttribute("freq_center", H5::PredType::NATIVE_DOUBLE, attr_dataspace);
    freq_center_attr.write(H5::PredType::NATIVE_DOUBLE, &m_freq_center);

    // Write frequency range attribute
    H5::Attribute freq_range_attr = dataset.createAttribute("freq_range", H5::PredType::NATIVE_DOUBLE, attr_dataspace);
    freq_range_attr.write(H5::PredType::NATIVE_DOUBLE, &m_freq_range);

    file.flush(H5F_SCOPE_LOCAL);
}

void HamakerCoefficients::print_coeffs()
{
    for (unsigned int h = 0; h < m_nHarmonics; h++) {
        for (unsigned int t = 0; t < m_nPowerTheta; t++) {
            for (unsigned int f = 0; f < m_nPowerFreq; f++) {
                auto coeff = get_coeff(h, t, f);
                std::cout << coeff.first << ", "
                          << coeff.second << std::endl;
            }
        }
    }
    std::cout << std::endl;
}
