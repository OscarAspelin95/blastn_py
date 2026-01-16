FROM python:3.14.2-bookworm

ENV BLAST_VERSION="2.17.0"

RUN mkdir -p /usr/src/deps \
	&& wget https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/${BLAST_VERSION}/ncbi-blast-${BLAST_VERSION}+-x64-linux.tar.gz \
	&& tar -xf ncbi-blast-${BLAST_VERSION}+-x64-linux.tar.gz \
	&& chmod +x ncbi-blast-${BLAST_VERSION}+/bin/blastn \
	&& mv ncbi-blast-${BLAST_VERSION}+/bin/blastn /usr/local/bin/ \
	&& chmod +x ncbi-blast-${BLAST_VERSION}+/bin/makeblastdb \
	&& mv ncbi-blast-${BLAST_VERSION}+/bin/makeblastdb /usr/local/bin/ \
	#
	&& rm -rf /usr/src/deps



COPY requirements.txt .

RUN apt-get update \
	&& pip install -r requirements.txt

WORKDIR /app
ENTRYPOINT ["/bin/bash"]
