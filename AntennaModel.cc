#include "AntennaModel.h"

std::ostream& operator<<(std::ostream& os, ElementResponseModel model)
{
    switch (model)
    {
        case Unknown: os << "Unknown"; break;
        case Hamaker: os << "Hamaker"; break;
        case LOBES  : os << "LOBES"; break;
        case OSKAR  : os << "OSKAR"; break;
        default     : os.setstate(std::ios_base::failbit);
    }
    return os;
}
