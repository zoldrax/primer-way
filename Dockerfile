FROM python:2-slim
RUN apt-get update && apt-get install -y build-essential libgomp1 wget git && \
        pip install --no-cache-dir pysam pyvcf primer3-py && \
        cd /tmp && wget http://public.lanl.gov/jgans/tntblast/tntblast-2.04.tar.gz && \
        tar zfx tntblast*.tar.gz && rm tntblast*.tar.gz* && cd /tmp/tntblast* && \
        ./configure --enable-openmp && make CXXFLAGS='-std=c++03 -O1' install && rm -R /tmp/tntblast* && \
        apt-get purge -y build-essential && apt-get autoremove -y && apt-get clean && rm -rf /var/lib/apt/lists/*
RUN cd /opt && git clone https://github.com/zoldrax/primer-way.git
ENV PATH="/opt/primer-way:${PATH}"
WORKDIR /data
CMD ["start.sh"]
