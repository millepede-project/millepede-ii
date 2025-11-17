FROM almalinux:9.6
RUN dnf --enablerepo=crb install -y cmake gcc-c++ gcc git lapack  openblas zlib-devel git openblas-openmp gfortran blas-devel openblas-devel lapack-devel texlive-latex texlive-dvips ghostscript diffutils make  flex bison 
RUN git clone -b Release_1_9_8 https://github.com/doxygen/doxygen.git
WORKDIR doxygen
WORKDIR build
RUN pwd 
RUN cmake ../ -G "Unix Makefiles" 
RUN make 
RUN make install 
