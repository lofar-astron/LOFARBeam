#include "ElementResponse.h"

namespace LOFAR {
namespace StationResponse {

std::ostream& operator<<(std::ostream& os, ElementResponseModel model )
{
    switch (model)
    {
        case Unknown: os << "Unknown"; break;
        case Hamaker: os << "Hamaker"; break;
        case LOBES  : os << "LOBES"; break;
        case OSKARDipole : os << "OSKARDipole"; break;
        case OSKARSphericalWave : os << "OSKARSphericalWave"; break;
        default     : os.setstate(std::ios_base::failbit);
    }
    return os;
}

} // namespace StationResponse
} // namespace LOFAR
