#!/bin/bash

# create the docker
docker build -t spectroscape-v1.1.2 .

# run spectroscape
docker run -it spectroscape-v1.1.2:latest spectroscape 

# run the docker with bash
docker run -it spectroscape-v1.1.2:latest bash

