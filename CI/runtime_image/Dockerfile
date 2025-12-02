FROM gitlab.desy.de:5555/millepede/millepede-ii/mp2dev_base
WORKDIR /pede
RUN git clone --recursive https://gitlab.desy.de/millepede/millepede-ii.git src 
RUN pwd 
RUN ls 
RUN mkdir build 
WORKDIR /pede/build
RUN cmake ../src 
RUN cmake --build . 
RUN cmake --install .
ENV LD_LIBRARY_PATH=/usr/local/lib
