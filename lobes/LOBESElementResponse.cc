#include "LOBESElementResponse.h"

#include <map>

namespace LOFAR {
namespace StationResponse {

LOBESElementResponse::LOBESElementResponse(std::string name)
{
}

void LOBESElementResponse::response(
        double freq,
        double theta,
        double phi,
        std::complex<double> (&response)[2][2]) const
{
}

std::shared_ptr<LOBESElementResponse> LOBESElementResponse::getInstance(std::string name)
{
    static std::map<std::string, std::shared_ptr<LOBESElementResponse>> name_response_map;

    auto entry = name_response_map.find(name);
    if (entry == name_response_map.end()) {
        entry = name_response_map.insert(entry, {name, std::make_shared<LOBESElementResponse>(name)});
    }
    return entry->second;

};

} // namespace StationResponse
} // namespace LOFAR

