
add_executable(tStation tStation.cc)

target_link_libraries(tStation stationresponse)


add_test(station tStation)
