#include "OskarSphericalWaveCoeff.h"

size_t OskarSphericalWaveCoefficients::get_index(
    const unsigned int element,
    const unsigned int l,
    const unsigned int m) const
{
    return element * m_nr_pols * m_nr_tetm * m_l_max * m_l_max +
                                                   l * m_l_max +
                                                       m;
}

// Constructor for reading coeff from file
OskarSphericalWaveCoefficients::OskarSphericalWaveCoefficients(
    std::string& filename)
{
    read_coeffs(filename);
};

size_t OskarSphericalWaveCoefficients::get_nr_elements() const
{
    return m_nr_elements;
}

size_t OskarSphericalWaveCoefficients::get_l_max() const
{
    return m_l_max;
}

size_t OskarSphericalWaveCoefficients::get_nr_coeffs() const
{
    return m_nr_pols * m_nr_tetm * m_l_max * m_l_max;
}

std::complex<double>* OskarSphericalWaveCoefficients::get_alpha_ptr(
    const unsigned int element,
    const unsigned int l,
    const unsigned int m)
{
    assert(element < get_nr_elements());
    auto l_max = get_l_max();
    //assert(m <= l);
    assert(l <= l_max);
    size_t index = get_index(element, l, m);
    return m_coeff.data() + index;
}

void OskarSphericalWaveCoefficients::read_coeffs(
    std::string& filename)
{
    // Open file
    std::cout << "read oskar datafile: " << filename << std::endl;
    m_h5_file.reset(new H5::H5File(filename, H5F_ACC_RDONLY));

    // Disable HDF5 error prints
    H5::Exception::dontPrint();
}

bool OskarSphericalWaveCoefficients::read_frequency(
    const unsigned int freq)
{
    if (m_freq != freq) {
        std::string dataset_name = std::to_string((int) (freq / 1e6));

        try {
            // Open dataset
            H5::DataSet dataset = m_h5_file->openDataSet(dataset_name);
            m_dataset_available = true;

            // Read dataset dimensions
            H5::DataSpace dataspace = dataset.getSpace();
            unsigned int rank = dataspace.getSimpleExtentNdims();
            assert(rank == m_dataset_rank);

            // Get dimensions
            hsize_t dims[rank];
            dataspace.getSimpleExtentDims(dims, NULL);
            m_nr_elements = dims[0];
            assert(dims[1] == m_nr_pols); // pola, polb
            assert(dims[2] == m_nr_tetm); // te, tm
            m_l_max = dims[3]; // l_max
            assert(dims[4] == m_l_max); // m_abs

            // Read coeff
            m_coeff.resize(get_nr_elements() * get_nr_coeffs());
            assert(dims[0]*dims[1]*dims[2]*dims[3]*dims[4]==m_coeff.size());
            H5::DataType data_type = dataset.getDataType();
            assert(data_type.getSize() == sizeof(std::complex<double>));
            dataset.read(m_coeff.data(), data_type, dataspace);
        } catch (H5::FileIException& e) {
            std::cerr << "Could not load dataset for frequency " << dataset_name << " Mhz" << std::endl;
            m_dataset_available = false;
        }
    }

    m_freq = freq;
    return m_dataset_available;
}

void OskarSphericalWaveCoefficients::print_alpha(
    const unsigned int element)
{
    if (!m_dataset_available) {
        return;
    }

    const int l_max = get_l_max();

    for (int l = 1; l <= l_max; ++l) {
        for (int abs_m = l; abs_m >=0; --abs_m) {
            auto alpha_ptr = get_alpha_ptr(element, l, abs_m);
            std::cout << *alpha_ptr;
            if (abs_m > 0) {
                std::cout << ", ";
            }
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
}