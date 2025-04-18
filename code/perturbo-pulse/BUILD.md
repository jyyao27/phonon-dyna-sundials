NOTE:
-----

All build operations assume that the user's environment has been initialized
in the way specified by the NVIDIA HPC 22.5 documentation.  Here is a link to
the relevant details:

https://docs.nvidia.com/hpc-sdk/hpc-sdk-install-guide/index.html#install-linux-end-usr-env-settings

HDF5 Build:
-----------

With MPI:

    CC=mpicc FC=mpifort CXX=mpic++ FCFLAGS=-fPIC ./configure --prefix=/home/donnie/hdf5 --enable-fortran --enable-parallel --enable-symbols --enable-trace --enable-codestack

Without MPI:

    CC=nvc FC=nvfortran CXX=nvc++ CFLAGS="-mcmodel=medium" FCFLAGS="-fPIC -Mlarge_arrays -mcmodel=medium" ./configure --prefix=/home/donnie/hdf5 --enable-fortran --enable-parallel=no --enable-symbols --enable-trace --enable-codestack

Using GCC toolchain:

    ./configure --prefix=/home/donnie/hdf5 --enable-fortran --enable-parallel=no --enable-symbols --enable-trace --enable-codestack

Debug settings:

    CC=nvc FC=nvfortran CXX=nvc++ CFLAGS="-tp=px -fPIC -mcmodel=medium" FCFLAGS="-tp=px -fPIC -Mlarge_arrays -mcmodel=medium" ./configure --prefix=/home/donnie/hdf5 --enable-fortran --enable-parallel=no --enable-symbols --enable-trace --enable-codestack --enable-asserts --enable-internal-debug=all --enable-diags=yes --enable-memory-alloc-sanity-check --enable-build-mode=debug

Quantum Espresso Build:
-----------------------

No MPI, no OpenMP:

    FC=nvfortran CC=nvc FCFLAGS="-Mlarge_arrays -mcmodel=medium" CFLAGS="-mcmodel=medium" ./configure --enable-parallel=no --enable-openmp=no

No MPI, yes OpenMP:

    FC=nvfortran CC=nvc FCFLAGS="-Mlarge_arrays -mcmodel=medium" CFLAGS="-mcmodel=medium" ./configure --enable-parallel=no --enable-openmp=yes

Maybe this?
-----------

-acc -ta=tesla -Minfo=accel -fast -Mfree -mp -Mlarge_arrays

*   -acc - Enable OpenACC directives
*   -mp - Enable OpenMP directives

*   -Minfo=accel - Emit information about accelerator region targeting.

*   -Mfree - Process Fortran source using free form specifications.



