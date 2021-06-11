FROM python:3.8-slim
RUN apt-get update && apt-get install -y build-essential libgomp1 wget git && \
        pip install --no-cache-dir pysam pyvcf primer3-py && \
        cd /opt && git clone https://github.com/zoldrax/primer-way.git && \
        cd /tmp && tar zfx /opt/primer-way/tntblast*.tar.gz && rm /opt/primer-way/tntblast*.tar.gz* && cd /tmp/tntblast* && \
        ./configure --enable-openmp && make CXXFLAGS='-std=c++03 -O1' install && rm -R /tmp/tntblast* && \
        apt-get purge -y build-essential && apt-get autoremove -y && apt-get clean && rm -rf /var/lib/apt/lists/*
ENV PATH="/opt/primer-way:${PATH}"
WORKDIR /data
CMD ["start.sh"]
