# check disk usage of docker. 
docker system df

# if it is not for CI, remove unused images and layers er safe. 
# remove all stoped containners.
docker container prun -f 

# -f options means no prompt before removing a artifacts. 

# list all th unused containsers.
docker ps --filter status=exited --filter status=dead -q 

# remove all the containers.
docker stop $(docker ps -q )
docker container prune

# another method to remove all the containers. 
docker rum $(docker ps -a -q)

# remove dangling docker images
docker image prun -f 

# remove all the iamges. 
docker image prun -a -f

# remove anonymous volumes
docker volume prune -f 

# remove all the volumes
docker volume rm -a -f 

# finally, the most important step. remove all the build cache.
docker buildx prune -f 

# remove everything
docker system prune -f 

# remove volumes as well
docker system prune --volumes -af




