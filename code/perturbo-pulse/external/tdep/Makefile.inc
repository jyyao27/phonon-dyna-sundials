
# Specify fortran compiler
#FC = mpifort -fpp  #-DGPx11 -Dclusterprogressbar
#F77C = mpifort

# And the flag that tells the compiler where to put modules
#MODULE_FLAG = -module

# Specify optimization
#OPT = -O0
#OPTLOW = 

# BLAS/LAPACK paths
#blaslapackLPATH = -L/opt/intel/oneapi/mkl/2022.2.1/lib/intel64
#blaslapackIPATH = -I/opt/intel/oneapi/mkl/2022.2.1/include
#blaslapackLIBS = -qmkl=cluster
# MPI
#incLPATHmpi = -L/opt/apps/mpi/mvapich2-2.3.7-1_intel-2021.7.1/lib
#incIPATHmpi = -I/opt/apps/mpi/mvapich2-2.3.7-1_intel-2021.7.1/include
#incLIBSmpi = -lmpi -lmpifort
# fftw
#incLPATHfft = 
#incIPATHfft = -I/opt/intel/oneapi/mkl/2022.2.1/include/fftw
#incLIBSfft = 
# HDF5
#incLPATHhdf = #-L/opt/apps/hdf5-1.12.2_intel-2021.7.1/lib
#incIPATHhdf = #-I/opt/apps/hdf5-1.12.2_intel-2021.7.1/include
#incLIBShdf = #-lhdf5 -lhdf5_fortran
# c-compiler
#CC = icc 

# paths to things.
incbasedir = ./
bindir = ./
libdir = ./
