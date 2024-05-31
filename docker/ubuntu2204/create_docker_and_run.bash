#!/bin/bash

# create the docker
docker build -t spectroscape-v1.2.0 .

# run spectroscape
docker run -it --rm  spectroscape-v1.2.0:latest spectroscape 
output=$?
status="compiled successfully;"
if [ $output -eq 0 ]; then
    echo "spectroscape run successfully"
    # run the docker with bash
    echo docker run -it --rm spectroscape-v1.2.0:latest bash
else
    status="compiled failed;"
    echo "spectroscape run failed"
fi

echo `date` `pwd` $status tested by Long > ../log.txt
