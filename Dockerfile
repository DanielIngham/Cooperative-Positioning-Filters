FROM ubuntu:24.04

RUN apt-get update && \
	apt-get install -y --no-install-recommends \
		gcc \
		g++ \
		gnuplot \
		make && \
	apt-get clean && \ 
	rm -rf /var/lib/apt/lists/*

RUN mkdir -p ./app/output

# Copy files maintaining the structure
COPY . ./app
WORKDIR ./app

CMD ["make", "run"]
