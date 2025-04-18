# Parallel Quicksort

The files `parallel_qsort.c` and `parallel_qsort.h` contain a quicksort
implementation parallelized with OpenMP tasks.

This implementation follows the Tony Hoare quicksort algorithm since it
generates fewer swaps than the more widely known Lamuto quicksort.  Wikipedia
has an excellent article on quicksort.

To understand the basics of how OpenMP is used to parallelize the quicksort,
see [this PDF](../../docs/sc16-openmp-booth-tasking-ruud.pdf) in the Perturbo
`docs` directory.  As suggested by the slide-deck, once the subproblem gets
below a certain size, the implementation just uses the C Standard Library
`qsort()` to complete the sort, since the overhead of threading and task
management exceeds the cost of just doing the rest of the task sequentially.

Compiling and running with GCC:

    gcc -O2 -fopenmp -Wall -Werror parallel_qsort.c test_parallel_qsort.c -o test_parallel_qsort

Compiling and running with NVIDIA C compiler:

    nvc -fast -mcmodel=medium -Mlarge_arrays -Minfo=all -mp=multicore -Wall -Werror parallel_qsort.c test_parallel_qsort.c -o test_parallel_qsort

Running:

    export OMP_NUM_THREADS=32
    ./test_parallel_qsort
