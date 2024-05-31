#!/bin/bash

# create the docker
docker build --no-cache -t spectroscape-ubuntu1804 .

# run spectroscape
docker run -it  spectroscape-ubuntu1804:latest SpectralArchive/build/bin/spectroscape 
container_name=$(docker ps -a --format "{{.Names}}" | head -n 1)
echo $container_name
docker cp $container_name:/SpectralArchive/package ./
echo ".deb and .sh package copied to ./package"


output=$?
status="compiled successfully;"
if [ $output -eq 0 ]; then
    echo "spectroscape run successfully"

    # run the docker with bash
    echo docker run -it --rm spectroscape-ubuntu1804:latest bash
    
else
    status="compiled failed;"
    echo "spectroscape run failed"
fi

echo `date` `pwd` $status tested by Long > ../log.txt

# clear the cache, release the storage of docker build 
# to check the space used by docker, run 
# docker system df 
# reference: https://depot.dev/blog/docker-clear-cache 
docker container prune -f
