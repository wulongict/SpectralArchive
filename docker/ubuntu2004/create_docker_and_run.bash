#!/bin/bash

# create the docker
docker build -t spectroscape-v1.1.3 .

# run spectroscape
docker run -it spectroscape-v1.1.3:latest spectroscape 

# run the docker with bash
docker run -it spectroscape-v1.1.3:latest bash

