#include <iostream>

#include "../Station.h"

#include "../BeamFormer.h"
#include "../LofarMetaDataUtil.h"

using namespace LOFAR::StationResponse;

int main()
{
    std::string name = "CS001LBA";
    vector3r_t position;
    ElementResponseModel model = ElementResponseModel::Hamaker;

    Station station(name, position, model);

    double time;
    double freq;
    vector3r_t direction = {0.0, 0.0, 1.0};
    matrix22c_t response;

    auto antenna0 = Element::Ptr(new Element(0, station.get_element_response()));
    auto antenna1 = Element::Ptr(new Element(1, station.get_element_response()));

    station.set_antenna(antenna0);

    std::cout << response[0][0] << std::endl;

    response = station.response(time, freq, direction);

    std::cout << response[0][0] << std::endl;

    auto beam_former = BeamFormer::Ptr(new BeamFormer());
    beam_former->add_antenna(antenna0);
    beam_former->add_antenna(antenna1);
    station.set_antenna(beam_former);
    response = station.response(time, freq, direction);
    std::cout << response[0][0] << std::endl;

    casacore::MeasurementSet ms("/home/vdtol/data/imagtest96.MS");

    Station::Ptr station0 = readStation(ms, 0, model);
    response = station0->response(time, freq, direction);
    response = station0->response(time, freq, direction);
    std::cout << response[0][0] << std::endl;

    return 0;
}

