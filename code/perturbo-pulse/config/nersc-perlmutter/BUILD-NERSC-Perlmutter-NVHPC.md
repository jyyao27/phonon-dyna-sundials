# Build and Run Perturbo on NERSC Perlmutter with MPI, OpenMP and OpenACC

These directions show how to set up the Perturbo program to run on the
NERSC Perlmutter compute cluster, using multiple compute nodes with MPI,
multiple CPUs with OpenMP, and GPU acceleration with OpenACC.

# TL;DR Summary

If you want to just give it a go, here are instructions to follow on
NERSC Perlmutter for a GPU build:

```sh
# Set up the environment for NVHPC Toolkit version 23.1 on NERSC Perlmutter
module unload cudatoolkit/12.2
module load nvhpc/23.1
module load cray-libsci/23.02.1.1
module load cray-hdf5-parallel
module load python/3.11

# Get the source code for Quantum ESPRESSO version 7.2 and Perturbo.  Then,
# replace the QE7.2 configure script with the version included in Perturbo.
git clone https://gitlab.com/QEF/q-e.git
cd q-e
git checkout qe-7.2
git clone git@github.com:perturbo-code/perturbo-dev.git
cp perturbo-dev/config/nersc-perlmutter/qe72-configure-fixed install/configure

# Configure Quantum ESPRESSO to use the NVHPC compiler + OpenMP + OpenACC.
# Note this is one long command that spans multiple lines.
FC=ftn F90=ftn MPIF90=ftn CC=cc \
  FFLAGS="-fast -Mlarge_arrays -mcmodel=medium" CFLAGS="-fast -mcmodel=medium" \
  LDFLAGS="-L/opt/nvidia/hpc_sdk/Linux_x86_64/23.1/math_libs/lib64" \
  BLAS_LIBS="-lblas" LAPACK_LIBS="-llapack" \
  ./configure --enable-parallel=yes --enable-openmp=yes \
  --with-cuda=/opt/nvidia/hpc_sdk/Linux_x86_64/23.1/cuda \
  --with-cuda-cc=80 --with-cuda-runtime=12.0

# Build the necessary parts of Quantum ESPRESSO 7.2.
make clean
make pw ph pp w90

# Build Perturbo using the provided make.sys file for NERSC Perlmutter.
cd perturbo-dev
cp config/make_nersc_nvhpc.sys ./make.sys
make clean
make
```

# Detailed Build Notes

Building Perturbo on NERSC Perlmutter is somewhat fussy because the NVIDIA
compiler toolchain is itself a bit fussy.  As of this writing (May 2024),
the `nvfortran` compiler is the only compiler we are aware of that can use
_both OpenMP and OpenACC at the same time_, with OpenMP generating parallel
CPU code, and OpenACC generating NVIDIA GPU code.  At the same time, Perturbo
invariably seems to expose `nvfortran` compiler bugs and other dependency
issues.  Thus, moving to a new version of NVHPC is usually a challenge.

**These instructions have been verified to work with NVHPC Toolkit 23.1.**
Specifically, we have identified [an `nvfortran` compiler bug](https://forums.developer.nvidia.com/t/nvhpc-23-9-vectorization-issue/292768)
in version 23.9 that prevents us from using that version of the Toolkit. This bug still exists in the latest NVHPC 24.5.

>   **NOTE:**  Currently, NERSC doesn't keep up with the absolute latest
>   versions of the NVIDIA HPC Toolkit (NVHPC) on Perlmutter.  The reason
>   is that NERSC also tries to provide a CUDA-aware MPI implementation
>   usable across multiple compilers; upgrading the NVHPC Toolkit requires
>   the rebuilding of many other dependencies on the cluster.

>   **NOTE:**  The NVIDIA compilers are fussy in general, but the fussiness
>   is exacerbated by using OpenBLAS.  The NVIDIA HPC Toolkit includes
>   GPU-aware libraries like cublas, and using the NVIDIA-provided libraries
>   seems to produce a much more stable program.

## Step 1:  Set up Perlmutter environment.

Make sure to switch to the NVIDIA HPC Toolkit version 23.1.

```sh
module load nvhpc/23.1
module load cray-hdf5-parallel
```

The details of what modules are loaded may vary depending on what you are
doing, but this is what the output of `module list` should _roughly_ look
like after everything is configured properly:

```
Currently Loaded Modules:
  1) craype-x86-milan
  2) libfabric/1.15.2.0
  3) craype-network-ofi
  4) xpmem/2.6.2-2.5_2.38__gd067c3f.shasta
  5) perftools-base/23.12.0
  6) cpe/23.12
  7) craype-accel-nvidia80
  8) gpu/1.0
  9) craype/2.7.30                         (c)
 10) cray-dsmml/0.2.2
 11) cray-mpich/8.1.28                     (mpi)
 12) PrgEnv-nvhpc/8.5.0
 13) nvhpc/23.1                            (c)
 14) cray-libsci/23.02.1.1                 (math)
 15) cray-hdf5-parallel/1.12.2.3           (io)
 16) evp-patch
 17) python/3.11                           (dev)
 18) conda/Miniconda3-py311_23.11.0-2

  Where:
   mpi:   MPI Providers
   math:  Mathematical libraries
   io:    Input/output software
   c:     Compiler
   dev:   Development Tools and Programming Languages
```

## Step 2:  Download QuantumESPRESSO 7.2 and Perturbo.

Quantum ESPRESSO (QE) has a generally straightforward build process, but
there are some wrinkles to be aware of.  Specifically there is a bug in
the QE7.2 `configure` script that must be patched before building it on
Perlmutter.  This is solved by copying in a patched version of the script
from Perturbo's codebase.

First, clone the QE repository and switch to version 7.2.

```sh
# Go to the home directory to download Quantum ESPRESSO
cd ~
git clone https://gitlab.com/QEF/q-e.git

# Switch to the version of Quantum ESPRESSO tagged as 7.2
cd q-e
git checkout qe-7.2
```

Next, clone the Perturbo repository _as a subdirectory of the Quantum
ESPRESSO codebase_.  We do this step now so we may copy in the fixed
`configure` script.

>   **NOTE:**  The QuantumEspresso 7.2 release has
>   [a bug](https://gitlab.com/QEF/q-e/-/issues/632) in its configuration
>   script that prevents it from recognizing the NVIDIA `nvfortran` compiler
>   through the NERSC `ftn` wrapper.  This bug can be fixed by overwriting
>   the QE7.2 `install/configure` script with [this file](./qe72-configure-fixed).
>   (Note there is also a `configure` script in the top-level QE directory;
>   _you must overwrite the one in the `install` subdirectory_.)
>
>   This bug has been fixed in QE7.3, so the extra step of copying in this
>   file can be removed when Perturbo migrates to QE7.3.

```sh
# Clone Perturbo's codebase inside Quantum ESPRESSO.
# NOTE:  Still in the q-e directory!!
git clone git@github.com:perturbo-code/perturbo-dev.git

# Patch the QE7.2 install/configure script
cp perturbo-dev/config/nersc-perlmutter/qe72-configure-fixed install/configure
```

## Step 3:  Configure and build QuantumESPRESSO 7.2.

Next, configure QE to build with the NVIDIA HPC compilers, using MPI,
OpenMP and OpenACC.  (Examples are given later of how to disable various
technologies.)  NERSC strongly recommends the use of their compiler-wrapper
programs, to ensure that all essential parameters are passed to the
compilers.

>   **NOTE:**  Specifying optimization levels above `-O2` and `-fast` cause
>   nondeterminism bugs in Perturbo!  Thus, we strongly recommend sticking
>   with `-fast` for optimized builds using the NVIDIA compilers.

```sh
# Configure Quantum ESPRESSO with compiler switches required by Perturbo.
# Note that this is one single command wrapped across multiple lines.
FC=ftn F90=ftn MPIF90=ftn CC=cc \
  FFLAGS="-fast -Mlarge_arrays -mcmodel=medium" CFLAGS="-fast -mcmodel=medium" \
  LDFLAGS="-L/opt/nvidia/hpc_sdk/Linux_x86_64/23.1/math_libs/lib64" \
  BLAS_LIBS="-lblas" LAPACK_LIBS="-llapack" \
  ./configure --enable-parallel=yes --enable-openmp=yes \
  --with-cuda=/opt/nvidia/hpc_sdk/Linux_x86_64/23.1/cuda \
  --with-cuda-cc=80 --with-cuda-runtime=12.0
```

Finally, build the necessary parts of Quantum ESPRESSO 7.2.

>   **NOTE:**  Using `-j` to enable parallel compilation (e.g.
>   `make pw ph pp w90 -j`) is broken on QE 7.2, so this step takes some
>   time.  Go get a cup of coffee and read a nice paper about condensed
>   matter physics while you wait.

```sh
# Build the necessary parts of QE 7.2.
make clean
make pw ph pp w90
```

>   **NOTE:**  If you see the C preprocessor `cpp` explicitly invoked during
>   the Quantum ESPRESSO build process then something has gone wrong with
>   the QE `configure` step.  The `nvfortran` compiler has a built-in
>   preprocessor, and QE knows to use this if it detects the NVHPC compiler.

### OPTIONAL:  Run QE 7.2 Tests

Once the Quantum ESPRESSO build process finishes, you may want to run the
QE 7.2 test suite.  However, you should note that it's common for at least
a few tests to fail, so this may be of questionable value.  It can be useful
to see if the QE 7.2 build fails catastrophically; if _many_ tests fail then
the QE7.2 build is probably not reliable.

```sh
cd test-suite
make run-tests
# [lots of test output!]
cd ..
```

## Step 4:  Configure and build Perturbo.

Recall that we already cloned Perturbo in Step 2 above.  Thus, there should
already be a `perturbo-dev` subdirectory within the `q-e` directory.

The Perturbo build process must be configured by putting a `make.sys` file
in the `perturbo-dev` directory; this file includes Perturbo-specific build
settings.  Multiple example files are included in the `perturbo-dev/config`
directory.  The `perturbo-dev/config/make_nersc_nvhpc.sys` file is provided
specifically for NERSC Perlmutter, and can simply be copied to the path
`perturbo-dev/make.sys`.

```sh
# NOTE:  Start in the ~/q-e directory!!
cd perturbo-dev
cp config/make_nersc_nvhpc.sys ./make.sys
```

Next it should be straightforward to build Perturbo.

```sh
# The "make clean" is just in case any build artifacts are hanging around.
make clean
make
```
>   **NOTE:**  Certain flags must be passed to the NVIDIA compilers to ensure
>   that Perturbo will run correctly.  If you see bad behavior from Perturbo,
>   a good first step is to rebuild it from scratch (`make clean && make`),
>   and ensure that the necessary flags are present.  These are as follows:
>
>   *   `-mpreprocess` to enable the C preprocessor in `nvfortran`, which is
>       off by default
>
>   *   `-Mlarge_arrays -mcmodel=medium` to allow for arrays and allocations
>       larger than 2GB in size
>       ([see here](https://forums.developer.nvidia.com/t/segmentation-fault-with-openmp/213527/2) for more details)
>
>   *   `-fast` (or `-O2`) should be the **maximum** optimization level;
>       `-O3` or `-O4` will generate spurious results
>
>   If the above instructions are carefully followed, these flags will be
>   picked up from the Quantum ESPRESSO configuration step.

## Step 5:  Testing Perturbo with `perturbopy` test suite.

Perturbo includes some regression tests for verifying that the program
operates correctly; it's important to perform this verification after
building it, especially due to the fussiness of the NVIDIA compilers.

### Install `perturbopy` Package

Here are some specific instructions for setting up `perturobpy` on NERSC
Perlmutter.  (For reference, here are some general instructions for
[installing the `perturbopy` package](https://perturbopy.readthedocs.io/en/latest/index.html).)

>   **NOTE:**  I don't know if this is the best approach to use, but it
>   certainly works.  Please let me know if you have a better approach.

I set up the `perturbopy` project as a separate codebase outside of the
Quantum ESPRESSO and Perturbo codebases, since it's not required for
`perturbopy` to live within QE, and it just keeps things cleaner and
more well-separated.

From my home directory:

```sh
# Go back to the user's top-level directory.
cd ~

# Set up Python and install virtualenv project.
module load python/3.11

# Check out perturbopy sources, and set up a virtual environment
# inside the project.
git clone https://github.com/perturbo-code/perturbopy.git
cd perturbopy

# Run the "venv" module from the Python standard library
# ("python -m venv ..."), to generate a new virtual environment
# in the "venv" directory ("... venv").  Then activate that
# environment ("source venv/bin/activate").
python -m venv venv
source venv/bin/activate

# Install perturbopy and all of its dependencies.
pip install .
```

>   **NOTE:**  If you make any edits to the `perturbopy` sources after
>   completing these instructions, or if you `git pull` on the `perturbopy`
>   repository, you will need to run `pip install -U .` again in the main
>   `perturbopy` directory to propagate these changes into the Python
>   virtual environment.

### Configure the Perturbo Test Suite

Once `perturbopy` is installed, the Perturbo test suite must be configured
so it knows how to run Perturbo.  This topic is very complex, so we will
start with a very simple configuration that just runs Perturbo on a single
machine without MPI or SLURM.

Back to the `perturbo-dev` directory:

```sh
# Make a config_machine.xml file from config_machine_serial.yml,
# which has only basic invocations for various Perturbo programs.
cd ~/q-e/perturbo-dev/tests/config_machine
cp cp config_machine_serial.yml config_machine.yml
```

To make sure that the test suite can find the `perturbo.x` and `qe2pert.x`
binaries, update the `config_machine.yml` file with the full paths to your
`perturbo-dev/bin/perturbo.x` and `perturbo-dev/bin/qe2pert.x` files.
**Failure to do this will result in all tests failing.**

```sh
# Go back to the perturbo-dev/tests directory
cd ..

# Run the test suite
python run_tests.py
```

### Running Perturbo Tests with MPI

TODO

# Disabling Various Technologies

This section describes how to disable the use of MPI, OpenMP or OpenACC.
Keep in mind that there are two distinct build steps - building QE and
building Perturbo - and some technologies are disabled only in one build
stage, whereas others can be disabled in one or both stages.

**NOTE:  EVERYTHING BEYOND THIS POINT HAS NOT BEEN TESTED FOR
NVHPC 23.1 COMPILER.**

## Disabling MPI

Disabling MPI support must be done during the QuantumESPRESSO build step.
Instead of specifying `--enable-parallel=yes`, specify `--enable-parallel=no`
or `--disable-parallel` to disable MPI.  Also, `MPIF90` no longer needs to
be specified.  Like this:

```sh
# Configure Quantum ESPRESSO with compiler switches required by Perturbo.
# Note that this is one single command wrapped across multiple lines.
FC=ftn F90=ftn CC=cc \
    FFLAGS="-fast -Mlarge_arrays -mcmodel=medium" CFLAGS="-fast -mcmodel=medium" \
    LDFLAGS="-L/opt/nvidia/hpc_sdk/Linux_x86_64/23.1/math_libs/lib64" \
    BLAS_LIBS="-lblas" LAPACK_LIBS="-llapack" \
    ./configure --enable-parallel=no --enable-openmp=yes \
    --with-cuda=/opt/nvidia/hpc_sdk/Linux_x86_64/23.1/cuda \
    --with-cuda-cc=80 --with-cuda-runtime=12.0
```

## Disabling OpenMP

OpenMP can be disabled in the QuantumESPRESSO build step, or in the Perturbo
build step, or in both steps.

To disable OpenMP in the QE build step, specify `--enable-openmp=no` or
`--disable-openmp` to the QE `configure` script.

```sh
# Configure Quantum ESPRESSO with compiler switches required by Perturbo.
# Note that this is one single command wrapped across multiple lines.
FC=ftn F90=ftn MPIF90=ftn CC=cc \
  FFLAGS="-fast -Mlarge_arrays -mcmodel=medium" CFLAGS="-fast -mcmodel=medium" \
  LDFLAGS="-L/opt/nvidia/hpc_sdk/Linux_x86_64/23.1/math_libs/lib64" \
  BLAS_LIBS="-lblas" LAPACK_LIBS="-llapack" \
  ./configure --enable-parallel=yes --enable-openmp=no \
  --with-cuda=/opt/nvidia/hpc_sdk/Linux_x86_64/23.1/cuda \
  --with-cuda-cc=80 --with-cuda-runtime=12.0
```

To disable OpenMP in the Perturbo build step, edit the `make.sys` file:
in the section marked "OpenMP" make sure to comment out all lines other
than the ones that say `-nomp`.

## Disabling OpenACC

OpenACC can be disabled in the QuantumESPRESSO build step, or in the Perturbo
build step, or in both steps.

Use of OpenACC is _enabled_ by default in QuantumESPRESSO.  To disable
OpenACC, make sure to specify `--enable-openacc=no` or `--disable-openacc`
to the QE `configure` script.  Also, the CUDA-related options should also
be removed.

```sh
# Configure Quantum ESPRESSO with compiler switches required by Perturbo.
# Note that this is one single command wrapped across multiple lines.
FC=ftn F90=ftn MPIF90=ftn CC=cc \
    FFLAGS="-fast -Mlarge_arrays -mcmodel=medium" CFLAGS="-fast -mcmodel=medium" \
    LDFLAGS="-L/opt/nvidia/hpc_sdk/Linux_x86_64/22.7/math_libs/lib64" \
    BLAS_LIBS="-lblas" LAPACK_LIBS="-llapack" \
    ./configure --enable-parallel=yes --enable-openmp=yes --enable-openacc=no

FC=ftn F90=ftn MPIF90=ftn CC=cc \
  FFLAGS="-fast -Mlarge_arrays -mcmodel=medium" CFLAGS="-fast -mcmodel=medium" \
  LDFLAGS="-L/opt/nvidia/hpc_sdk/Linux_x86_64/23.1/math_libs/lib64" \
  BLAS_LIBS="-lblas" LAPACK_LIBS="-llapack" \
  ./configure --enable-parallel=yes --enable-openmp=yes --enable-openacc=no
```

# Error in Message Passing (mp) module 9100

This error can occur when running Perturbo on NERSC Perlmutter CPU nodes.
It is caused by building Perturbo on a GPU environment and then trying to
run it in a CPU environment.  **Note this seems to have nothing to do with
whether Perturbo is built with GPU acceleration;** rather, it has to do with
whether Pertubo is built in the CPU or GPU environment on Perlmutter.  This
can happen easily, since the Perlmutter login nodes have a GPU.

If you want to run Perturbo on CPU nodes and you encounter this error 9100,
you need to rebuild the program in a CPU-only environment so it will run
correctly.  These modifications must be made to the above build steps:

*   _After_ completing all other `module unload` / `module load` steps, but
    _before_ beginning to build QE7.2, run this:  `module swap gpu cpu`

*   When configuring QE7.2, you should exclude all `"--with-cuda..."`
    arguments.  Your configure command will look like this:

    ```sh
    FC=ftn F90=ftn MPIF90=ftn CC=cc \
      FFLAGS="-fast -Mlarge_arrays -mcmodel=medium" CFLAGS="-fast -mcmodel=medium" \
      LDFLAGS="-L/opt/nvidia/hpc_sdk/Linux_x86_64/23.1/math_libs/lib64" \
      BLAS_LIBS="-lblas" LAPACK_LIBS="-llapack" \
      ./configure --enable-parallel=yes --enable-openmp=yes
    ```

*   _After_ you have copied the `perturbo-dev/config/make_nersc_nvhpc.sys`
    file to `perturbo-dev/make.sys`, and _before_ you compile Perturbo,
    you must edit the `make.sys` file:  Comment out the two lines that
    enable OpenACC, and uncomment the two lines that disable OpenACC (i.e.
    the lines that add `“--noacc”` to the arguments.)

This should allow you to build a CPU-only version of Perturbo that will
run correctly on the Perlmutter CPU nodes.
