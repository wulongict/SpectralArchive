for i in `find -type d | grep ubuntu`; 
do 
cd $i && ./create_docker_and_run.bash; 
 cd ..; 
done