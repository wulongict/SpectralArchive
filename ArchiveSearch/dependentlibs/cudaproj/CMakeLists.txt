cmake_minimum_required(VERSION 3.8)
project(cuda_dotproduct LANGUAGES NONE)

#get_property(languages GLOBAL PROPERTY ENABLED_LANGUAGES)
#MESSAGE("LANGUAGES ENABLES " ${languages})

# warning: every time you change the cmakelist file, you should use
# tool-> cmake-> reset cmake cache and reload project
# according to this page: https://alessandrogar.io/post/osquery-development-with-qtcreator-and-clion/

#unset(CMAKE_CUDA_COMPILER)

# Note: our server is using the Pascal architecture.
#
# Fri Aug  6 18:59:39 2021
# +-----------------------------------------------------------------------------+
# | NVIDIA-SMI 384.130                Driver Version: 384.130                   |
# |-------------------------------+----------------------+----------------------+
# | GPU  Name        Persistence-M| Bus-Id        Disp.A | Volatile Uncorr. ECC |
# | Fan  Temp  Perf  Pwr:Usage/Cap|         Memory-Usage | GPU-Util  Compute M. |
# |===============================+======================+======================|
# |   0  TITAN X (Pascal)    Off  | 00000000:02:00.0 Off |                  N/A |
# | 23%   41C    P8    17W / 250W |   4119MiB / 12189MiB |      0%      Default |
# +-------------------------------+----------------------+----------------------+
# the gencode can be found in the following link:
# https://arnon.dk/matching-sm-architectures-arch-and-gencode-for-various-nvidia-cards/
#

# Fermi†	Kepler†	Maxwell‡	Pascal	Volta	Turing	Ampere	Lovelace*	Hopper**
# sm_20	    sm_30	sm_50	    sm_60	sm_70	sm_75	sm_80	sm_90?	    sm_100c?
# ----      sm_35   sm_52	    sm_61	sm_72	-----   sm_86
# ----      sm_37 	sm_53	    sm_62
# choose 60 61 62 is good for our GPU
# set(CMAKE_CUDA_ARCHITECTURES 60 61 62 70 72 80 86 75 87)



add_library(dpgpu SHARED dpcuda.cu)
set_property(TARGET dpgpu PROPERTY CUDA_STANDARD 11)
# add -fPIC to shared library
set_target_properties(dpgpu PROPERTIES POSITION_INDEPENDENT_CODE ON)
get_target_property(LIBB_INCLUDES dpgpu INCLUDE_DIRECTORIES)
#message("target libdpgpu.so has following include directory: " ${LIBB_INCLUDES})
target_link_libraries(dpgpu msms mstoolkit)


add_library(dpgpu_static STATIC dpcuda.cu)
set_property(TARGET dpgpu_static PROPERTY CUDA_STANDARD 11)
target_link_libraries(dpgpu_static msms_static mstoolkit)


add_executable(dpcudatest dpcuda_main.cu)
set_property(TARGET dpcudatest PROPERTY CUDA_STANDARD 11)
target_link_libraries(dpcudatest dpgpu msms_static mstoolkit)





