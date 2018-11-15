# $Id: CMakeLists.txt 18775 2011-09-06 13:36:45Z zwieten $

include(LofarPackageVersion)

lofar_add_library(elementresponse
  Package__Version.cc
  ElementResponse.cc)
