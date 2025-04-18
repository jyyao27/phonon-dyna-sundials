# Contributing to Perturbo: code and documentation guidelines

## Code guidelines

Perturbo is written in modern Fortran (f90). To ensure code readability and minimize errors, we seek to practice a high standard of coding. 

Here is a list of the Perturbo coding guidelines:

### Indentation

* Use 3-space indentation throughout the code.

* Indent all subroutines, `if` blocks, loops, etc.

### Style 

* Put white spaces between operations:
```fortran
! Bad
a=b+c

! Good
a = b + c
```

* Separate the subroutine (function) arguments and array indices with spaces after comas, do not put comas between a paranthesis and a variable. Example:
```fortran
! Bad
call my_function(var1,var2,var3 )

! Good
call my_function(var1, var2, var3)
```

* Put blank lines between blocks of subroutines (functions), especially, before and after the `if` blocks, `do` loops, etc.

* Separate the local variables of a subroutine (function) from the input parameters. Example:
```fortran
subroutine init_dist_gaussian(kgrid, dist, e0)
   implicit none
   type(grid), intent(in) :: kgrid                     !< k-point grid
   real(dp), intent(in) :: e0                          !< Energy center of the Gaussian
   real(dp), intent(out) :: dist(kgrid%numb, kgrid%nk) !< Occupation function
   ! local variables
   integer :: ik, ib
   real(dp) :: enk
```

### Naming

* Avoid short one- or two- letter variable names.

* Variable names should not be very long, but descriptive. Example: `kgrid`.

* We recommend start the loop variable names with "i" and number of iterations with "n". Example: 
```fortran
do istep = 1, nstep
   ...
end do
```

### Profiling

* Set up the start and stop clocks around your most time consuming blocks of code. The timing will be automatically output by Perturbo. Example:
```fortran
call start_clock('dynamics_tot')
do istep = 1, nstep
...
end do
call stop_clock('dynamics_tot')
```

* Memory profiling will be added in future.

### Performance

* Fortran is column-major, which means that the first index of array should be the fastest. Example:
```fortran
do i = 1, 1000
   do j = 1, 1000
      array(j, i)
   end do
end do
```

* The largest loop must be the inner one:
```fortran
do i = 1, 5
   do j = 1, 1000
      array(j, i)
   end do
end do
```

* Perturbo supports the hybrid OpenMP+MPI parallelization. Please, follow the examples from the source code to implement parallelization of new functions/subroutines.

### Output formats

* Relatively small amounts of data (up to 1 MB) should be output as YAML files using `YAML_utils.f90` or to ASCII text files (optional).

* Large arrays should be output to HDF5 files using `HDF5_utils.f90`.

Examples of all the three ways of output (YAML, ASCII, and HDF5) can be found in the code.

### Miscellaneous

* Use the double precision for Fortran reals. Example: `real(dp) :: x`

* Specify the double precision for reals. Example: 
```fortran
! Bad
x = 0.0

! Good
x = 0.0_dp
```

* Always use `implicit none`

* Use `intent(in)`, `intent(out)`, or `intent(inout)` for **EVERY** input parameter for **EVERY** function or subroutine you implement.

* Avoid using the global variables. Global variables within a module are allowed if they are appropriate for a given module.

* Avoid putting multiple Fortran instructions in one lines using `;`. Example:
```fortran
var1 = 1.0_dp
var2 = 2.0_dp
```
is better than:
```fortran
var1 = 1.0_dp; var2 = 2.0_dp;
```

* When using modules, specify the names of the imported subroutines, functions, or variables. Example: 
```fortran
! Bad
use boltz_grid

! Good
use boltz_grid, only: grid, init_boltz_grid, init_boltz_qgrid
```

* Never assume eager or short-circuit evaluation of logical statements in Fortran. For example, treat `.and.` syntax to be **NOT** short-circuit,
```fortran
! Bad, should be avoided!
if(present(step_max) .and. step_max > 0)

! Good
if(present(step_max)) then
    if(step_max > 0) ...
endif
```
If the variable `step_max` is not defined, the first condition in the `.and.` logical statements is false and the compilers (such as `ifort` and `nvfortran`) will still continue to evaluate the second condition because the short-circuit is not guaranteed by the Fortran standard. Consequently, the compilers (such as `ifort` and `nvfortran`) will probably raise `segment fault` error even if the failure of the first condition. 

* Group your functions and subroutines that are dedicated to a certain computation into modules. Example: `module selfenergy`.

* As a general rule, new code should follow the style of the existing core Perturbo source code.

## Testing Perturbo

When doing new implementations, please ensure that they do not affect the results of the existing calculation modes. To do so, run the [Perturbo testsuite](https://perturbopy.readthedocs.io/en/latest/testsuite.html). Create new tests for the new features you implemented.

## Adding a new input parameter

* There is only one place where you should introduce new variables: `./docs/input_parameters_perturbo.yml` for `perturbo.x` and `./docs/input_parameters_qe2pert.yml` for `qe2pert.x`. Please follow the format for other variable and choose a 'variable family' that fits. For example, if you would like to introduce a new string variable that affects the behavior of the ultrafast dynamics, it might look like so:
```
new_variable:
   family: dynamics
   type: string
   default: "'text'"
   options:
      "'text'": Description 1.
      "'text2'": Description 2.
   len: 80
   attributes: save
   description: Description of the variable.
```
Please follow examples for other parameter types.

* Run the `./utils/input_parameters.py` script to generate the Fortran files.

* Add a new Fortran variable in the `pert_param.f90` file, follow the procedure for existing variables in the file (including the declaration, initialization, reading from the input file, and MPI broadcasting).

The dedicated page on the [Perturbo website](https://perturbo-code.github.io/mydoc_param_perturbo.html) is automatically generated based on these `YAML` files.

## Documenting the source code

Source code documentation is indispensible for the collaboration and sustainability of the project. We use [Doxygen](https://www.doxygen.nl) for the generation of the documentation and follow the Doxygen Fortran comments format.

### Reading the existing Perturbo source code documentation

Open the `index.html` file located in the [path-to-perturbo]/docs/html with a browser.

### Documentation guidelines

* Before every subroutine or function, put a comment section consisting of a few lines. It should not be very long, but explaining the purpose of a given function or subroutine and containing any additional information the coder finds appropriate. Usually, 2-3 sentences is enough. But for the key subroutines or functions it can reach up to 20-30 lines of text with (sub-)sections, latex equations, etc.

Start the comment block before a subroutine or a function with `!>` and continue with `!!`. Example for a subroutine (similar for function):
```fortran
!> Initialize the electron occupation function with
!! the Guassian distribution around a given energy.
subroutine init_dist_gaussian(kgrid, dist, e0, smear, lhole)
```

* Similarly, put a header before every module, like so:
```fortran
!> Description of what this module does. The common purpose
!! of the subroutines and functions in this module.
module test
...
end module test
```
Module headers might be larger than the subroutine ones.

* Comment on a variable after the `!<` symbol and after the variable declaration:
```fortran
integer :: nstep !< number of steps of the ultrafast simulation
```

Explain **all** the input and output variables of a function or subroutine. Comments on the local variables are *optional* but appreciated. If multiple variables are declared on the same line and you would like to comment on them, write the declaration as multiple lines. Example:

```fortran
real(dp), intent(in) :: e0, smear
```
split onto two lines and comment:

```fortran
real(dp), intent(in) :: e0 !< Energy center of the Gaussian
real(dp), intent(in) :: smear !< Guassian smearing
```

**Align** the comment blocks for variables. E.g., replace this:
```fortran
subroutine init_dist_gaussian(kgrid, dist, e0, smear, lhole)
   implicit none
   type(grid), intent(in) :: kgrid  !< k-point grid
   real(dp), intent(in) :: e0 !< Energy center of the Gaussian
   real(dp), intent(in) :: smear !< Guassian smearing
   real(dp), intent(out) :: dist(kgrid%numb, kgrid%nk) !< Occupation function
```
with this:
```fortran
subroutine init_dist_gaussian(kgrid, dist, e0, smear, lhole)
   implicit none
   type(grid), intent(in) :: kgrid                     !< k-point grid
   real(dp), intent(in) :: e0                          !< Energy center of the Gaussian
   real(dp), intent(in) :: smear                       !< Guassian smearing
   real(dp), intent(out) :: dist(kgrid%numb, kgrid%nk) !< Occupation function
```

* Thoughout the code, put the necessary amount of comments to explain all important code blocks and non-trivial decisions (at least a few words).

### Doxygen syntax

Enhance your comments with stylized text and cross-references to other parts of the code. Doxygen uses standard markdown syntax and provides additional tools for referencing code.

### Doxygen uses basic markdown syntax
As in markdown, use `#` for sections, `##` for subsections, etc. Use the (sub-)sections **only** if you write a large comment block (more than five lines).
In the module, function, or subroutine headers, the (sub-)sections are not required if they are short.

Use `>` for quote blocks.
> quote

Lists can be done with `*`. Here is an example of a nested list:

```fortran
!! * one
!!   + nested one
!!   + nested two
!! * two
```

*Italicize* text with one `*` around the text and **bold** with two `*`:
```fortran
!! *italic* **bold**
```

If you would like to ~~strike~~ text, you can do so with:
```fortran
!! ~~text~~
```

Code `spans` should be surrounded by the backquotes: "`". Use those, when you refer to Fortran variables in the code.
```fortran
!! `var`
```

Latex in-line equations inside the comment blocks should be surounded by `\f$`:  
```fortran
!! \f$ a = b \f$
```

Separate line Latex equations are designated by `\f[` and `\f]`:
```fortran
!! \f[
!! f(x) = \sum_i g_i(x)
!! \f]
```

Links to websites:
```fortran
!! [perturbo tutorial](https://perturbo-code.github.io/mydoc_running_perturbo.html) 
```

Notes and todo sections can be done with the `\note` and `\todo` commands, which will be very nicely rendered by Doxygen. Example:
```fortran
!! \todo
!!   - reading and writing ( real, strings )
!!   - hdf_get_*
!!     - hdf_get_obj_name  (h5iget_name_f)
!!     - hdf_get_obj_id (not needed)
!!   - error checking,
!!     - check dims when reading
!!     - check group name when reading/writing
!! \note In this subroutine, we assume ...
```

### Links to other subroutines, functions, modules, etc.

If you reference a file of the source code, e.g.:
```fortran
!! see boltz_grid.f90 for more details
```
a link to the file will be automatically generated.

To link a subroutine within the same module or global, type `::` in front of the subroutine name without spaces:
```fortran
!! see ::init_dist_lorentz
```

To link a subroutine (function) from a different module, type `<module name>::` in front of the subroutine (function) name without spaces:
```fortran
!! see boltz_grid::init_boltz_grid
```

The same works for references to types from different modules:
```fortran
!! see boltz_grid::grid
```

## Generating documentation
If you have modified the Perturbo source code and created the Doxygen comment blocks, you can run `doxygen` to re-generate the documentation.

### Doxygen installation:

On Linux OS, it can be installed as follows:
```bash
sudo apt-get install doxygen
```

On mac, Doxygen can be installed with brew:
```bash
brew install doxygen
```

More information [here](https://www.doxygen.nl/manual/install.html).

### Generating documentation

To generate documentation, run `doxygen` from your terminal in the root folder of Perturbo source code.

The documentation in the HTML format will be generated in the folder: [path-to-perturbo]/docs/html .
Other output formats are available with Doxygen, which can be configured in the `Doxyfile`. This file can be found in the root folder of the source code.

## Technical notes

The following parameters were changed from their default values in the Doxygen configure file:

```bash
EXTRACT_ALL            = YES
OPTIMIZE_FOR_FORTRAN   = YES
PROJECT_NAME           = "PERTURBO"
OUTPUT_DIRECTORY       = docs
PROJECT_LOGO           = docs/perturbo_logo.pdf
EXTENSION_MAPPING      = f90
```
