#ifndef ELEMENT_H
#define ELEMENT_H

#include <complex>
#include <memory>

#include "Antenna.h"
#include "ElementResponse.h"
#include "Types.h"

namespace LOFAR {
namespace StationResponse {

class Element : public Antenna
{
public:

    typedef std::shared_ptr<Element> Ptr;

    Element(int id, ElementResponse::Ptr element_response) :
        m_id(id),
        m_element_response(element_response)
    {}

    virtual matrix22c_t response(
        real_t time,
        real_t freq,
        const vector3r_t &direction,
        const Options &options) const final override;

private:
    int m_id;
    double m_orientation;
    ElementResponse::Ptr m_element_response;
};

} // namespace StationResponse
} // namespace LOFAR

#endif
