FROM quay.io/fenicsproject/stable
MAINTAINER M. Parno

USER root

RUN apt-get update && \
    apt-get install -yy pwgen npm nodejs-legacy python3-pip && \
    npm install npm@latest -g && \
    npm install -g configurable-http-proxy && \
    pip3 install --upgrade pip && \
    pip3 install --upgrade numpy && \
    pip3 install jupyter pandas jupyterlab && \
    apt-get clean && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

USER fenics

# Install MUQ
RUN cd /home/fenics && \
    mkdir Installations; mkdir Installations/MUQ_INSTALL && \
    git clone --depth 1 https://bitbucket.org/mituq/muq2.git && \
    cd muq2/; mkdir build; cd build;  \
    cmake -DCMAKE_INSTALL_PREFIX=/home/fenics/Installations/MUQ_INSTALL -DMUQ_USE_PYTHON=ON -DMUQ_HDF5_DIR=/usr/local/lib/python3.5/dist-packages/h5py/.libs/ ../ && \
    make -j2 install && \
    cd /home/fenics && \
    rm -r muq2

# Install hIPPYlib
RUN cd /home/fenics/Installations && \
    git clone https://github.com/hippylib/hippylib.git && \
    chmod -R o+rx hippylib


# Configure jupyter notebooks
COPY python3_config.json /usr/local/share/jupyter/kernels/python3/kernel.json
ENV LD_LIBRARY_PATH /home/fenics/Installations/MUQ_INSTALL/lib:/home/fenics/Installations/MUQ_INSTALL/muq_external/lib
ENV PYTHONPATH /home/fenics/Installations/MUQ_INSTALL/lib:/home/fenics/Installations/hippylib

USER root

# Install h5py with the MUQ-installed HDF5 library
RUN cd /home/fenics/Installations && \
    git clone --depth=1 https://github.com/h5py/h5py.git && \
    cd h5py && \
    python3 setup.py configure --hdf5=/home/fenics/Installations/MUQ_INSTALL/muq_external && \
    python3 setup.py install

RUN mkdir -p /home/fenics/.jupyter
COPY jupyter_notebook_config.py /home/fenics/.jupyter/jupyter_notebook_config.py


WORKDIR /home/fenics/
CMD ["bash"]
