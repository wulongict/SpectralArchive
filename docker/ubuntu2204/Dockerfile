# syntax=docker/dockerfile:1
FROM ubuntu:22.04

# install app dependencies
RUN apt update && apt install -y wget 
RUN wget https://github.com/wulongict/SpectralArchive/releases/download/v1.2.0/Spectroscape_CPU-1.2.0.deb
RUN apt install -y ./Spectroscape_CPU-1.2.0.deb

