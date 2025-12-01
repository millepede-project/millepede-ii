FROM gitlab-registry.cern.ch/linuxsupport/alma9-base

ARG ROOT_TARBALL=root_v6.36.06.Linux-almalinux9.6-x86_64-gcc11.5.tar.gz

WORKDIR /opt

COPY dependencies dependencies

RUN dnf update -q -y \
 && dnf install -y epel-release \
 && dnf --enablerepo=crb  install -y  $(cat dependencies)\
 && rm -f dependencies \
 && curl -O https://root.cern/download/${ROOT_TARBALL} \
 && tar -xzvf ${ROOT_TARBALL} \
 && rm -f ${ROOT_TARBALL} \
 && echo /opt/root/lib >> /etc/ld.so.conf \
 && ldconfig

ENV ROOTSYS /opt/root
ENV PATH $ROOTSYS/bin:$PATH
ENV PYTHONPATH $ROOTSYS/lib:$PYTHONPATH
ENV CLING_STANDARD_PCH none

RUN git clone -b Release_1_9_8 https://github.com/doxygen/doxygen.git
WORKDIR doxygen
WORKDIR build
RUN pwd 
RUN cmake ../ -G "Unix Makefiles" 
RUN make 
RUN make install 
