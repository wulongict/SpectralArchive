# syntax=docker/dockerfile:1
FROM ubuntu:18.04

# the default gcc version is 7.5.0, lower than the one required by spectroscape binary package.
# So let's install from source code

# First install the dependencies
RUN apt update
RUN apt install wget git cmake build-essential liblapack-dev  spawn-fcgi nginx -y
RUN wget https://github.com/Kitware/CMake/releases/download/v3.26.4/cmake-3.26.4-linux-x86_64.sh
RUN chmod +x cmake-3.26.4-linux-x86_64.sh
RUN ./cmake-3.26.4-linux-x86_64.sh --skip-license --prefix=/usr/local

RUN  git clone --recurse-submodules  https://github.com/wulongict/SpectralArchive.git --branch release
RUN cd SpectralArchive && ./cleanMake.bash
RUN cd SpectralArchive && ./compile.bash

