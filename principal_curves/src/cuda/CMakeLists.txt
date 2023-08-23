
FIND_PACKAGE(CUDA REQUIRED)

#set (CUDA_NVCC_FLAGS "-lineinfo" CACHE STRING "nvcc flags" FORCE)
set (CUDA_NVCC_FLAGS "-g -G -lineinfo" CACHE STRING "nvcc flags" FORCE)

SET (CUDA_VERBOSE_BUILD ON CACHE BOOL "nvcc verbose" FORCE) 
#HW capability 1.3

SET(LIB_TYPE STATIC)    #set the lib type
CUDA_ADD_LIBRARY(skelGPU ${LIB_TYPE} skelft.h skelft.cu )

