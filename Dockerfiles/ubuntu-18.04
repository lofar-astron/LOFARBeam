FROM ubuntu:18.04

RUN apt-get update
RUN apt-get install -y libboost-all-dev \
		               build-essential \
                       libpython3-dev \
                       python3-dev \
                       python3-numpy \
		               gfortran \
		               wget \
                       libsofa1-dev \
                       flex \
                       bison \
                       libbison-dev \
                       libatlas-base-dev \
                       liblapack-dev \
                       libcfitsio-dev \
                       wcslib-dev \
                       cmake 

#####################################################################
## BUILD CASACORE FROM SOURCE
#####################################################################
RUN mkdir /src
WORKDIR /src
RUN wget https://github.com/casacore/casacore/archive/v3.1.1.tar.gz
RUN tar xvf v3.1.1.tar.gz
RUN mkdir casacore-3.1.1/build
WORKDIR /src/casacore-3.1.1/build
RUN cmake -DCMAKE_INSTALL_PREFIX=/usr -DCMAKE_BUILD_TYPE=Release -DBUILD_DEPRECATED=ON -DBUILD_PYTHON=OFF -DBUILD_PYTHON3=ON ../
RUN make -j 4
RUN make install
RUN ldconfig

RUN apt-get install -y python3-pip
WORKDIR /src
RUN rm v3.1.1.tar.gz
RUN wget https://github.com/casacore/python-casacore/archive/v3.1.1.tar.gz
RUN tar xvf v3.1.1.tar.gz
WORKDIR /src/python-casacore-3.1.1
RUN pip3 install .
WORKDIR /
RUN python3 -c "from pyrap.tables import table as tbl"

#####################################################################
## BUILD STATION RESPONSE FROM SOURCE
#####################################################################
ADD . /src
WORKDIR /src
RUN mkdir build
RUN cd build && cmake -DCMAKE_INSTALL_PREFIX=/usr -DPYTHON_EXECUTABLE=$(which python3) ../ && make -j16 && make install
RUN touch "/usr/lib/python3.6/site-packages/lofar/__init__.py"
# installs in the site-packages for some reason, need to add this to the search path
ENV PYTHONPATH "/usr/lib/python3.6/site-packages:$PYTHONPATH"
RUN python3 -c "import lofar.stationresponse as lsr"

