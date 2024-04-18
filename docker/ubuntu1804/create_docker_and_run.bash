#!/bin/bash

# create the docker
docker build -t spectroscape-ubuntu1804 .

# run spectroscape
docker run -it spectroscape-ubuntu1804:latest SpectralArchive/build/bin/spectroscape 

# run the docker with bash
docker run -it spectroscape-ubuntu1804:latest bash

