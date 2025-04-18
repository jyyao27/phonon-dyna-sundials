
include ../../make.inc
include ../make.sys


PERTDIR=../pert-src

PERTOBJ = \
HDF5_utils.f90 \
qe_mpi_mod.f90 \
pert_const.f90 \
pert_utils.f90 \
sys_utils.f90 \
YAML_utils.f90 \
environment_setup.f90 \
wigner_seitz_cell.f90 \
triplet_cell.f90 \
polar_correction.f90 \
electron_wannier.f90 \
force_constant.f90 \
force_constant_thirdorder.f90 \
elph_matrix_wannier.f90 \
epr_hdf5_io.f90 \

# TODO:  Prepend nvgpu.f90 to above list for NVML output
# PERTMOD := nvgpu.f90 $(PERTOBJ)


Q2PTMOD = \
autogen_param_qe2pert.f90 \
input_param.f90 \
param_to_yaml_qe2pert.f90 \
YAML_utils_qe2pert.f90 \
electronic_data.f90 \
lattice_data.f90 \
lattice_tdep.f90 \
dvscf_mod.f90 \
dvhub_mod.f90 \
elph_matrix.f90 \
perturbo_data.f90



Q2PTSRC = \
get_igk_kqg.f90 \
pw_read_file.f90 \
pw_read_collected_to_evc.f90 \
pw_init_us_2.f90 \
pw_atomic_wfc.f90 \
pw_s_psi.f90 \
ph_transform_int.f90 \
ph_transform_int3.f90 \
ph_dvqpsi_us_local.f90 \
ph_dvqpsi_us_noloc.f90 \
ph_dvanqq.f90 \
ph_adddvscf.f90 \
ph_newdq.f90 \
ph_adddvhubscf.f90 \
ph_dwfc.f90 \
ph_swfc.f90 \
ph_delta_sphi.f90 \
ph_lr_orthoUwfc.f90 \
ph_dvqhub_barepsi_us.f90 \
calc_electron_wann.f90 \
calc_lattice_ifc.f90 \
tdep_lattice_ifc.f90 \
save_elph_mat_wann.f90 \
collect_elph_matrix.f90 \
qe2pert.f90


OBJ = $(PERTOBJ:.f90=.o) $(Q2PTMOD:.f90=.o) $(Q2PTSRC:.f90=.o)

# TODO:  Append nvgpu_info.o to OBJ for NVML output
# OBJ := $(OBJ) nvgpu_info.o


EXE = qe2pert.x

default: $(EXE)

$(EXE): $(OBJ)
	$(LD) $(LDFLAGS) -o $(EXE) $(OBJ) $(PHMODS) $(LRMODS) $(PWOBJS) $(QEMODS) $(LIBOBJS) $(QELIBS) $(LIBRARY_DIRS) $(LINKFLAGS)

%.o: ${PERTDIR}/%.f90
	$(LD) $(F90FLAGS) $(MODFLAG) -c $<

%.o: %.f90 
	$(LD) $(F90FLAGS) $(MODFLAG) -c $<

%.o: ${PERTDIR}/%.c
	$(CC) $(CFLAGS) -c $<

${PERTDIR}/nvgpu.f90: ${PERTDIR}/nvgpu_info.h ${PERTDIR}/nvgpu_info.c

${PERTDIR}/nvgpu_info.c: ${PERTDIR}/nvgpu_info.h

clean:
	$(RM) *.o *~ *.d *.mod *.i

