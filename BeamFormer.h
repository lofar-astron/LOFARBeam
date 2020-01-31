#ifndef BEAMFORMER_H
#define BEAMFORMER_H

#include <complex>

class BeamFormer : public Element
{
public:

    virtual void response(
        double time,
        double freq,
        double theta,
        double phi,
        std::complex<double> (&result)[2][2]) override;

    add_element(Element::Ptr element) {m_elements->insert(element);}

private:
    std::vector<Element::Ptr> m_elements;

};

#endif
