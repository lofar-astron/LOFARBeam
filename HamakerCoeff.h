#include <iostream>
#include <string>
#include <complex>
#include <cassert>
#include <vector>
#include <cstring>

#include <H5Cpp.h>

class HamakerCoefficients {
    public:
        // Constructor for reading coeff from file
        HamakerCoefficients(
            std::string& filename);

        // Constructor for writing coeff to file
        HamakerCoefficients(
            const double freq_center,
            const double freq_range,
            const unsigned int nHarmonics,
            const unsigned int nPowerTheta,
            const unsigned int nPowerFreq);

        size_t get_nr_coeffs() const;

        // Set
        void set_coeff(
            const unsigned int n,
            const unsigned int t,
            const unsigned int f,
            std::pair<std::complex<double>, std::complex<double>> value);

        void set_coeffs(
            const std::complex<double>* coeff);

        void set_coeffs(
            const std::vector<std::complex<double>> coeff);

        // Get
        std::pair<std::complex<double>, std::complex<double>> get_coeff(
            const unsigned int n,
            const unsigned int t,
            const unsigned int f);

        // HDF5 I/O
        void read_coeffs(
            std::string& filename);

        void write_coeffs(
            std::string& filename);

        // Debugging
        void print_coeffs();

    private:
        // Methods
        size_t get_index(
            const unsigned int n,
            const unsigned int t,
            const unsigned int f);

        // Parameters
        double m_freq_center;
        double m_freq_range;
        unsigned int m_nHarmonics;
        unsigned int m_nPowerTheta;
        unsigned int m_nPowerFreq;
        const unsigned int m_nInner = 2;

        // Data
        std::vector<std::complex<double>> m_coeff;

        // HDF5
        std::string m_dataset_name = "coeff";
        const unsigned int m_dataset_rank = 4;
};
