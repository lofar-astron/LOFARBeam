#ifndef OSKAR_SPHERICAL_WAVE_COEFF_H
#define OSKAR_SPHERICAL_WAVE_COEFF_H

#include <iostream>
#include <string>
#include <complex>
#include <cassert>
#include <vector>
#include <cstring>
#include <memory>
#include <mutex>

#include <H5Cpp.h>

class OskarSphericalWaveCoefficients {
    public:
        // Constructor for reading coeff from file
        OskarSphericalWaveCoefficients(
            std::string& filename);

        // Get
        size_t get_nr_elements() const;
        size_t get_l_max() const;
        size_t get_nr_coeffs() const;

        std::complex<double>* get_alpha_ptr(
            const unsigned int element,
            const unsigned int l = 0,
            const unsigned int m = 0);

        // HDF5 I/O
        void read_coeffs(
            std::string& filename);

        bool read_frequency(
            const unsigned int freq);

        // Debugging
        void print_alpha(unsigned int element = 0);

    private:
        // Methods
        size_t get_index(
            const unsigned int element,
            const unsigned int l,
            const unsigned int m) const;

        // Parameters
        unsigned int m_nr_elements;
        const unsigned int m_nr_pols = 2;
        const unsigned int m_nr_tetm = 2;
        unsigned int m_l_max;
        unsigned int m_freq = 0;

        // Data
        std::vector<std::complex<double>> m_coeff;

        // HDF5
        std::string m_filename;
        const unsigned int m_dataset_rank = 5;
        std::unique_ptr<H5::H5File> m_h5_file;
        bool m_dataset_available = false;
};

#endif
