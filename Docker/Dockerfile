#################################################################################################################
# Modified from the Edward Docker image: https://github.com/blei-lab/edward/tree/master/docker
#################################################################################################################

#################################################################################################################
# Base Images
#################################################################################################################
FROM ubuntu:14.04

#################################################################################################################
#           ENV Setting
#################################################################################################################
ENV CONDA_DIR /opt/conda
ENV PATH $CONDA_DIR/bin:$PATH

#################################################################################################################
#           Initial Setting
#################################################################################################################
RUN mkdir -p $CONDA_DIR && \
    echo export PATH=$CONDA_DIR/bin:'$PATH' > /etc/profile.d/conda.sh && \
    apt-get update && \
    apt-get install -y wget git libhdf5-dev g++ graphviz unzip && \
    wget --quiet https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh && \
    /bin/bash /Miniconda3-latest-Linux-x86_64.sh -f -b -p $CONDA_DIR && \
    rm Miniconda3-latest-Linux-x86_64.sh

#################################################################################################################
#           User Setting
#################################################################################################################
ENV NB_USER edward
ENV NB_UID 1000

RUN useradd -m -s /bin/bash -N -u $NB_UID $NB_USER && \
    mkdir -p $CONDA_DIR && \
    chown edward $CONDA_DIR -R && \
    mkdir -p /src && \
    chown edward /src

USER edward

#################################################################################################################
#           Python Setting
#################################################################################################################
# Python
ARG python_version=3.5.3-0
ARG python_qt_version=4
RUN mkdir /tmp/taigapy && \
    cd /tmp/taigapy && \
    wget 'https://stash.broadinstitute.org/plugins/servlet/archive/projects/CPDS/repos/taigapy?at=refs%2Fheads%2Fmaster' -O taigapy-master.zip && \
    unzip taigapy-master.zip && \
    conda install -y python=${python_version} && \
    pip install numpy six tensorflow keras prettytensor && \
    pip install ipdb pytest pytest-cov python-coveralls coverage==3.7.1 pytest-xdist pep8 pytest-pep8 pydot_ng && \
    pip install google-cloud attrs && \
    conda install Pillow scikit-learn notebook pandas matplotlib seaborn pyyaml h5py && \
    conda install -y pyqt=${python_qt_version} && \
    pip install edward && \
    conda clean -yt

ENV PYTHONPATH='/src/:$PYTHONPATH'

# install taigapy
RUN cd /tmp/taigapy && \
    python setup.py install


# install kubeque
RUN echo v3 && mkdir /tmp/kubeque && \
    cd /tmp/kubeque && \
    wget 'https://github.com/broadinstitute/kubeque/archive/master.zip' -O master.zip && \
    unzip master.zip && \
    cd kubeque-master && \
    python setup.py install 

#RUN pip install simplejson
RUN pip install json-tricks

# install dev version of kubeque
#USER root
#COPY kubeque-0.2.tar.gz /tmp
#USER edward
#RUN mkdir /tmp/kubeque && \
#    cd /tmp/kubeque && \
#    tar xzf /tmp/kubeque-0.2.tar.gz && \
#    cd kubeque-0.2 && \
#    python setup.py install
    
#################################################################################################################
#           WORK Jupyter
#################################################################################################################
WORKDIR /src

EXPOSE 8888

CMD jupyter notebook --port=8888 --ip=0.0.0.0

USER root
# give edward account sudo access.  Needed for accessing /host-var which kubeque uses for caching
RUN echo 'edward ALL = NOPASSWD : ALL' >> /etc/sudoers

# USER edward


