External dependencies included in source tree

Subdirectories were created as described below.


Install pybind11 into LOFARBeam source tree


    git clone https://github.com/pybind/pybind11.git
    cd pybind11
    git checkout v2.4.3
    mkdir build
    cd build
    cmake .. -DCMAKE_INSTALL_PREFIX=<srcdir>/LOFARBeam/external/pybind11
    make install


Install Eigen3 into LOFARBeam source tree

    git clone https://gitlab.com/libeigen/eigen.git
    cd eigen
    git checkout 'master@{2020-01-01 00:00:00}'
    mkdir build
    cd build
    cmake .. -DCMAKE_INSTALL_PREFIX=<srcdir>/LOFARBeam/external/eigen -DEIGEN_BUILD_PKGCONFIG=Off
    make install
