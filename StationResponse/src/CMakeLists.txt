# $Id$

include(LofarPackageVersion)

lofar_add_library(stationresponse
  Package__Version.cc
  AntennaField.cc
  AntennaFieldHBA.cc
  AntennaFieldLBA.cc
  AntennaModelHBA.cc
  AntennaModelLBA.cc
  DualDipoleAntenna.cc
  ITRFDirection.cc
  LofarMetaDataUtil.cc
  MathUtil.cc
  Station.cc
  TileAntenna.cc
  Types.cc)

lofar_add_bin_program(makeresponseimage makeresponseimage.cc)
