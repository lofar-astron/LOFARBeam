# $Id$

lofar_package(ElementResponse 0.1)

# Uncomment to check for unsafe conversions (gcc), for example conversion of
# size_t to unsigned int (truncation).
#add_definitions(-Wconversion)

add_subdirectory(include/ElementResponse)
add_subdirectory(src)
