#!/bin/bash

# create the docker
docker build -t spectroscape-v1.1.3-ubuntu1804 .

# run spectroscape
docker run -it spectroscape-v1.1.3-ubuntu1804:latest SpectralArchive/build/bin/spectroscape 

# run the docker with bash
docker run -it spectroscape-v1.1.3-ubuntu1804:latest bash

