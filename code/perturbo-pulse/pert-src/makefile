
include ../../make.inc
include ../make.sys


PERTMOD = \
HDF5_utils.f90 \
wigner_seitz_cell.f90 \
triplet_cell.f90 \
polar_correction.f90 \
electron_wannier.f90 \
force_constant.f90 \
force_constant_thirdorder.f90 \
elph_matrix_wannier.f90 \
pert_const.f90 \
qe_mpi_mod.f90 \
epr_hdf5_io.f90 \
pert_utils.f90 \
autogen_param_perturbo.f90 \
pert_param.f90  \
sys_utils.f90 \
YAML_utils.f90 \
param_to_yaml_perturbo.f90 \
pert_data.f90   \
YAML_utils_perturbo.f90 \
pert_output.f90 \
band_structure.f90 \
phonon_dispersion.f90 \
elphon_coupling_matrix.f90 \
boltz_utils.f90 \
vector_list.f90 \
boltz_grid.f90  \
bfield_utils.f90 \
imaginary_selfenergy.f90 \
transition_dipole.f90 \
phonon_anharm_prop.f90 \
phonon_anharm2.f90 \
boltz_trans_output.f90 \
boltz_trans_mod.f90 \
boltz_scatter_sizest.f90 \
ph-ph-svd_mod.f90 \
boltz_scatter.f90 \
boltz_scatter_integral.f90 \
sorting.f90 \
boltz_scatter_integral_tgt.f90 \
boltz_grid_neighbors.f90 \
boltz_external_field.f90 \
boltz_dynamics_solver.f90 \
boltz_dynamics_mod.f90 \
selfenergy.f90 \
cumulant_expansion.f90 \
cumulant_utility.f90 \
sundials.f90
# TODO:  Prepend nvgpu.f90 to above list for NVML output
# PERTMOD := nvgpu.f90 $(PERTMOD)


PERTSRC = \
environment_setup.f90 \
weight_dos.f90 \
load_data_pert.f90 \
calc_bands.f90 \
calc_ephmat.f90 \
calc_electron_mfp.f90 \
real_selfenergy.f90 \
electron_imsigma.f90 \
electron_dipole.f90 \
phonon_imsigma.f90 \
electron_resigma.f90 \
ph-ph-svd.f90 \
boltz_setup.f90 \
transport.f90 \
trans_postproc.f90 \
carrier_dynamics_run.f90 \
carrier_dynamics_postproc.f90 \
calc_selfenergy.f90 \
process_spectral.f90 \
calc_spectral_cumulant.f90 \
calc_trans_cumulant.f90 \
calc_trans_cumulant_opcond.f90 \
perturbo.f90


OBJ = $(PERTMOD:.f90=.o) $(PERTSRC:.f90=.o)

# TODO:  Append nvgpu_info.o to OBJ for NVML output
# OBJ := $(OBJ) nvgpu_info.o


EXE = perturbo.x

default: $(EXE)
 
$(EXE): $(OBJ)
	$(LD) $(LDFLAGS) -o $(EXE) $(OBJ) $(PHMODS) $(LRMODS) $(PWOBJS) $(QEMODS) $(LIBOBJS) $(QELIBS) $(LINKFLAGS)

%.o: %.f90
	$(LD) $(F90FLAGS) $(MODFLAG) -c $<

%.o: %.c
	$(CC) $(CFLAGS) -c $<

boltz_scatter.f90: cdyn_config.h

boltz_scatter_integral.f90: cdyn_config.h

nvgpu.f90: nvgpu_info.h nvgpu_info.c

nvgpu_info.c: nvgpu_info.h

clean:
	$(RM) *.o *~ *.d *.mod *.i *.optrpt

