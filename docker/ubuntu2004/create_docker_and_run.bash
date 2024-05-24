#!/bin/bash

# create the docker
docker build --no-cache -t spectroscape-ubuntu2004 .

# docker run -it spectroscape-ubuntu2004:latest  gcc --version 
# run spectroscape
docker run -it   spectroscape-ubuntu2004:latest spectroscape 
output=$?
status="compiled successfully;"
if [ $output -eq 0 ]; then
    echo "spectroscape run successfully"
    # run the docker with bash
    echo docker run -it --rm spectroscape-ubuntu2004:latest bash
else
    status="compiled failed;"
    echo "spectroscape run failed"
    echo `date` `pwd` $status tested by Long > ../log.txt
fi




