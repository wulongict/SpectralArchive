# part one, install GSL 1.16
# sudo apt-get install libgsl-dev
# the command above install a different version.
# use command below:
GSL_VERSION=$(gsl-config --version)
if [ "$GSL_VERSION" == "1.16" ]; then
  echo "GSL version 1.16 already installed!"
else
  echo GSL version is $GSL_VERSION
  wget -nc https://mirror-hk.koddos.net/gnu/gsl/gsl-1.16.tar.gz
  tar xvf gsl-1.16.tar.gz
  # shellcheck disable=SC2164
  (
    cd gsl-1.16 || exit
    ./configure
    make all -j 10
    sudo make install
  )
fi

# part two install other packages.
sudo apt-get install sqlite3 libsqlite3-dev libcgicc-dev liblog4cpp5-dev libboost-all-dev libfcgi-dev spawn-fcgi nginx curl  gnuplot gnuplot-qt  

# part three: install faiss
sudo apt-get install python3.8-dev swig python-numpy python-all-dev python3-all-dev liblapack-dev

# install faiss from soruce code
# download source code, unzip and run the following folder in top level folder.
# cmake -B build  .
# make -C build 
