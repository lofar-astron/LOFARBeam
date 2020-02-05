find_package(pybind11 REQUIRED)
find_package (Eigen3 REQUIRED NO_MODULE)
find_package(HDF5 COMPONENTS C CXX REQUIRED)

# pybind11_add_module(pylobes SHARED lobes.cc ElementResponse.cc LOBESElementResponse.cc)
# target_link_libraries (pylobes PUBLIC Eigen3::Eigen ${HDF5_LIBRARIES})

add_library(lobes SHARED LOBESElementResponse.cc)
install(TARGETS lobes)

configure_file(test.py test.py COPYONLY)

configure_file(F4far_new.py F4far_new.py COPYONLY)
configure_file(calc_modes_1deg_noback.py calc_modes_1deg_noback.py COPYONLY)
configure_file(hamaker_vs_lobes.py hamaker_vs_lobes.py COPYONLY)
configure_file(test_beam_model.py test_beam_model.py COPYONLY)

configure_file(CS302_coords.mat CS302_coords.mat COPYONLY)

# TODO Too large for git repo, get file from somewhere else
# configure_file(LBA_CS302_fine.mat LBA_CS302_fine.mat COPYONLY)



