## Perturbo compilation

In this folder you can find the files that are needed to compile PERTURBO for different platforms and compilation scripts in particular:

1. GCC
2. Intel compilers classic with(out) MPI
3. Intel compilers with(out) MPI

You can take one of the provided file as a basis, and modify it to suit your case by adding, for example, OpenMP compilation options or the path to the HDF5 library (it can be omitted if it is specified in _make.inc_ of QE).

Copy your chosen sample to the root folder _"perturbo"_:

```bash
$ cp config/make_gcc_serial.sys make.sys
```

and make the necessary changes to it. You will then be ready to compile PERTURBO:

```bash
$ make
```

## Configuration of Quantum Espresso

Even if you already have Quantum Espresso compiled, you need to configure the _make.inc_ file for it, since this file is used to compile PERTURBO. 

Below are Quantum Espresso compilation options for several cases:

1. GCC: 
```bash
$ ./configure --with-hdf5="/path/to/hdf5" CFLAGS=-O3 FFLAGS=-O3
```
2. Intel Compilers Classic:
```bash
$  ./configure CFLAGS=-O3 FFLAGS=-O3 --disable-parallel --with-hdf5="/path/to/hdf5"
```
3. Intel Compilers Classic with MPI:
```bash
$ ./configure CFLAGS=-O3 FFLAGS=-O3 F90=ifort MPIF90=mpiifort --with-hdf5="/path/to/hdf5"
```
4. Intel Compilers:
```bash
$ ./configure F90=ifx CC=icx CFLAGS=-O3 FFLAGS=-O3 --disable-parallel --with-hdf5="/path/to/hdf5"
```
5. Intel Compilers with MPI:
```bash
$ ./configure F90=ifx CC=mpiicx CFLAGS=-O3 FFLAGS=-O3 MPIF90=mpiifx --with-hdf5="/path/to/hdf5"
```
