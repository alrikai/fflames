FROM ubuntu:latest

LABEL maintainer="Alrik Firl afirlortwo@gmail.com" \
      version="1.0" \
      description="FractalFlame Dockerfile"

ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update --fix-missing && apt-get --fix-missing -y install \
    curl \
    gcc \ 
	g++ \ 
    build-essential \
    cmake \
    libopencv-dev \
    ffmpeg \
    libboost-all-dev \
    git

COPY ./fflames /fflames
COPY ./tests /tests
COPY main.cpp /.
COPY ./CMakeLists.txt* /

RUN mkdir -p /build
WORKDIR /build
RUN cmake ../. && make -j4
WORKDIR /

CMD ["bash", "-c", "./build/fflame_gen", "--help"]
