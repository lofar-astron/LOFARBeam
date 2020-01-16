#ifndef ELEMENT_RESPONSE_MODEL_H
#define ELEMENT_RESPONSE_MODEL_H

#include <ostream>

enum ElementResponseModel {
    Unknown,
    Hamaker,
    LOBES,
    OSKAR
};

std::ostream& operator<<(std::ostream& os, ElementResponseModel model);

#endif
