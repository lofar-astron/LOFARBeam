#include <string>

#include "DefaultCoeffHBA.cc"
#include "DefaultCoeffLBA.cc"

#include "HamakerCoeff.h"

void write_hba_coeff(
    std::string filename)
{
    HamakerCoefficients hba_coeff(
        default_hba_freq_center,
        default_hba_freq_range,
        default_hba_coeff_shape[0],
        default_hba_coeff_shape[1],
        default_hba_coeff_shape[2]);

    hba_coeff.set_coeffs(default_hba_coeff);

    hba_coeff.write_coeffs(filename);
}

void write_lba_coeff(
    std::string filename)
{
    HamakerCoefficients lba_coeff(
        default_lba_freq_center,
        default_lba_freq_range,
        default_lba_coeff_shape[0],
        default_lba_coeff_shape[1],
        default_lba_coeff_shape[2]);

    lba_coeff.set_coeffs(default_lba_coeff);

    lba_coeff.write_coeffs(filename);
}

int main (void)
{
    // Write hba coeff from DefaultCoeffHBA.cc
    std::string hba_coeff_filename("HamakerHBACoeff.h5");
    write_hba_coeff(hba_coeff_filename);

    // Write lba coeff from DefaultCoeffLBA.cc
    std::string lba_coeff_filename("HamakerLBACoeff.h5");
    write_lba_coeff(lba_coeff_filename);
}
