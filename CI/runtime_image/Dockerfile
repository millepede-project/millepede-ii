FROM gitlab.desy.de:5555/millepede/millepede-ii/mp2dev_base
WORKDIR /pede
RUN git clone --recursive https://gitlab.desy.de/millepede/millepede-ii.git src 
RUN pwd 
RUN ls 
RUN mkdir build install
COPY CMakeLists.txt src/CMakeLists.txt
WORKDIR /pede/build
RUN cmake -DCMAKE_INSTALL_PREFIX=/pede/install ../src 
RUN cmake --build . 
RUN cmake --install .
WORKDIR /pede
RUN cp install/mp2setup.sh /pede/mp2setup.sh
