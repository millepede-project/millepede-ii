FROM almalinux:9.6
RUN dnf --enablerepo=crb install -y cmake gcc-c++ gcc lapack  openblas zlib-devel git openblas-openmp gfortran blas-devel openblas-devel lapack-devel doxygen texlive-latex texlive-dvips ghostscript 
WORKDIR /pede
# RUN git clone --recursive https://gitlab.desy.de/millepede/millepede-ii.git src 
# RUN pwd 
# RUN ls 
# RUN mkdir build install
# COPY CMakeLists.txt src/CMakeLists.txt
# WORKDIR /pede/build
# RUN cmake -DCMAKE_INSTALL_PREFIX=/pede/install ../src 
# RUN cmake --build . 
# RUN cmake --install .
# WORKDIR /pede
# RUN source install/mp2setup.sh
