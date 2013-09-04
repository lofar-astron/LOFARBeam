# $Id$

# Create symbolic link to include directory.
execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_BINARY_DIR}/include/${PACKAGE_NAME})

# Install header files.
install(FILES
  AntennaField.h
  AntennaFieldHBA.h
  AntennaFieldLBA.h
  AntennaModelHBA.h
  AntennaModelLBA.h
  Constants.h
  DualDipoleAntenna.h
  ITRFDirection.h
  LofarMetaDataUtil.h
  MathUtil.h
  Station.h
  TileAntenna.h
  Types.h
  DESTINATION include/${PACKAGE_NAME})
