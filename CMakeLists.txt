# $Id$

lofar_package(StationResponse 0.1 DEPENDS Common ElementResponse)

include(LofarFindPackage)
lofar_find_package(Boost REQUIRED)
lofar_find_package(Casacore REQUIRED COMPONENTS casa measures ms tables images)

# Uncomment to check for unsafe conversions (gcc), for example conversion of
# size_t to unsigned int (truncation).
#add_definitions(-Wconversion)

add_subdirectory(include/StationResponse)
add_subdirectory(src)
