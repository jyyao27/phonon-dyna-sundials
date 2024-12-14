
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
phonon_anharm_prop.f90 \
phonon_anharm2.f90 \
boltz_trans_output.f90 \
boltz_trans_mod.f90 \
boltz_scatter.f90 \
boltz_scatter_integral.f90 \
boltz_grid_neighbors.f90 \
boltz_external_field.f90 \
boltz_dynamics_solver.f90 \
boltz_dynamics_mod.f90 \
selfenergy.f90 \
cumulant_expansion.f90 \
cumulant_utility.f90 \
sundials.f90


PERTSRC = \
environment_setup.f90 \
weight_dos.f90 \
load_data_pert.f90 \
calc_bands.f90 \
calc_ephmat.f90 \
calc_electron_mfp.f90 \
real_selfenergy.f90 \
electron_imsigma.f90 \
phonon_imsigma.f90 \
electron_resigma.f90 \
phonon_gruneisen.f90 \
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

EXE = perturbo.x

default: $(OBJ)
	$(LD) $(LDFLAGS) -o $(EXE) $(OBJ) $(PHMODS) $(LRMODS) $(PWOBJS) $(QEMODS) $(LIBOBJS) $(QELIBS) $(LINKFLAGS)


%.o: %.f90 
	$(LD) $(F90FLAGS) $(MODFLAG) -c $<

boltz_grid.o: boltz_grid.f90 boltz_grid_generate.f90 boltz_grid_load.f90
	$(LD) $(F90FLAGS) $(MODFLAG) -c $<

boltz_scatter.o: boltz_scatter.f90 setup_kq_pair.f90 compute_scatter_eph_g2.f90 load_scatter_eph_g2.f90 setup_qq_pair.f90 compute_scatter_ph3_psi2.f90 
	$(LD) $(F90FLAGS) $(MODFLAG) -c $<

clean:
	$(RM) *.o *~ *.d *.mod *.i
