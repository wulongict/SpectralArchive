# cd ~/code/faiss-1.7.1/ 
# Options to control GPU fetures.
# eneral options:
# -DFAISS_ENABLE_GPU=OFF in order to disable building GPU indices (possible values are ON and OFF),
# -DFAISS_ENABLE_PYTHON=OFF in order to disable building python bindings (possible values are ON and OFF),
# -DBUILD_TESTING=OFF in order to disable building C++ tests,
# -DBUILD_SHARED_LIBS=ON in order to build a shared library (possible values are ON and OFF),
# optimization-related options:
# -DCMAKE_BUILD_TYPE=Release in order to enable generic compiler optimization options (enables -O3 on gcc for instance),
# -DFAISS_OPT_LEVEL=avx2 in order to enable the required compiler flags to generate code using optimized SIMD instructions (possible values are generic, sse4, and avx2, by increasing order of optimization),



cmake -DBUILD_TESTING=ON -DFAISS_ENABLE_PYTHON=OFF -DBUILD_SHARED_LIBS=ON  -DFAISS_ENABLE_GPU=ON -DCMAKE_BUILD_TYPE=Release -B build .
# cmake -DBUILD_TESTING=ON -DFAISS_ENABLE_PYTHON=OFF -DBUILD_SHARED_LIBS=ON  -DFAISS_ENABLE_GPU=ON -DCMAKE_BUILD_TYPE=Release -B build . 
make -C build clean
make  --directory=build -j faiss 
# remove first
 sudo rm -rf /usr/local/include/faiss/
  sudo rm /usr/local/lib/libfaiss.a /usr/local/lib/libfaiss.so

sudo make -C build install

build/tests/faiss_test
