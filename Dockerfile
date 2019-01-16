## INSTALL R 3.5
FROM r-base


## INSTALL SEURAT ##
ADD install.R /software/scripts/install.R
ADD cluster.R /software/scripts/cluster.R
ADD cluster.py /software/scripts/cluster.py

ENV PATH "$PATH:/software/scripts"

# Install ssl, curl, and ssh
RUN apt-get update -qq
RUN apt-get install -t unstable -y libssl-dev libcurl4-openssl-dev libssh2-1-dev libhdf5-dev
RUN apt-get clean

RUN apt-get install -t unstable -y gfortran
RUN Rscript /software/scripts/install.R

RUN chmod -R a+x /software/scripts

## INSTALL PYTHON

RUN wget https://repo.continuum.io/miniconda/Miniconda3-4.5.4-Linux-x86_64.sh -O ~/miniconda.sh

RUN bash ~/miniconda.sh -b
ENV PATH /root/miniconda3/bin:$PATH
RUN apt-get update && apt-get install --no-install-recommends -y \
		build-essential \
		automake \
		zlib1g-dev \
		libxml2-dev \
		cmake \
		gnupg \
		lsb-release \
		libfftw3-dev \
		libboost-iostreams-dev \
		default-jdk\ 
		curl && \
	export CLOUD_SDK_REPO="cloud-sdk-$(lsb_release -c -s)" && \
	echo "deb http://packages.cloud.google.com/apt $CLOUD_SDK_REPO main" | tee -a /etc/apt/sources.list.d/google-cloud-sdk.list && \
	curl https://packages.cloud.google.com/apt/doc/apt-key.gpg | apt-key add - && \
	apt-get update && apt-get install -y google-cloud-sdk

RUN pip install python-igraph
RUN pip install louvain

ADD 1M_neurons_filtered_gene_bc_matrices_h5.h5 /home/1M_neurons_filtered_gene_bc_matrices_h5.h5

RUN apt-get -y install libhdf5-dev
RUN apt-get -y install git-all

RUN pip install --upgrade pip && \
	pip install numpy && \
	pip install Cython && \
	pip install pybind11 && \
	git clone https://github.com/nmslib/hnsw.git /software/hnswlib && \
	cd /software/hnswlib/python_bindings && python setup.py install && cd ../../..

COPY scCloud-0.8.0.tar.gz /software
COPY extract_ADT_UMIs /software/extract_ADT_UMIs

ADD https://raw.githubusercontent.com/klarman-cell-observatory/KCO/master/docker/monitor_script.sh /software

RUN tar -xzf /software/scCloud-0.8.0.tar.gz -C /software && \
	cd /software/scCloud-0.8.0 && pip install . && cd ../.. && \
	cd /software/extract_ADT_UMIs && make all && cd ../.. && \
	chmod a+rx /software/monitor_script.sh

RUN pip install scanpy
ADD 1M_neurons_matrix_subsampled_1M.h5 /home/1M_neurons_matrix_subsampled_1M.h5
ADD 1M_neurons_matrix_subsampled_500K.h5 /home/1M_neurons_matrix_subsampled_500K.h5
ADD 1M_neurons_matrix_subsampled_250K.h5 /home/1M_neurons_matrix_subsampled_250K.h5
ADD 1M_neurons_matrix_subsampled_100K.h5 /home/1M_neurons_matrix_subsampled_100K.h5
ADD 1M_neurons_matrix_subsampled_25K.h5 /home/1M_neurons_matrix_subsampled_25K.h5
ADD 1M_neurons_filtered_gene_bc_matrices_h5.h5 /home/1M_neurons_filtered_gene_bc_matrices_h5.h5

CMD ["bash"]
