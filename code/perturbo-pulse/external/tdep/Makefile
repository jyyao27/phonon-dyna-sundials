#
# Makefile for Ollelib
#
# Computer dependent stuff
include ../../../make.inc
include ../../make.sys
include Makefile.inc


# Where to put things
OBJECT_PATH=./
MODULE_PATH=./
DEST_PATH=./

#LPATH = $(blaslapackLPATH) $(incLPATHhdf) $(incLPATHfft) $(incLPATHmpi)
#IPATH = $(blaslapackIPATH) $(incIPATHhdf) $(incIPATHfft) $(incIPATHmpi)
#LIBS = $(blaslapackLIBS) $(incLIBShdf) $(incLIBSfft) $(incLIBSmpi)

#OPT = -O0 -fbacktrace -fcheck=all -finit-real=nan -finit-derived -fmax-errors=10 --warn-all # -ffpe-trap=zero,overflow,underflow
#OPT = -O0 -fbacktrace -fcheck=all -finit-real=nan -finit-derived -fmax-errors=10 --pedantic -Wall -Wextra -pedantic -Wcast-align -Wdisabled-optimization -Wmissing-include-dirs -Wshadow -Wunused -fdiagnostics-show-option -Wstrict-overflow=0 -Wno-maybe-uninitialized -Wno-compare-reals -Wno-character-truncation -Werror #-Wrealloc-lhs
#OPT = -O0 -check bounds -check uninit -check pointers -traceback -g -fpe0
#OPT = -traceback -fpe0 #-check bounds -check uninit -check pointers -traceback -g -fpe0


# For debugging/throwing warnings and all that.
#FLG0 = -O0 -fbacktrace -fcheck=all -finit-real=nan -finit-derived -fmax-errors=10 --pedantic -Wall -Wextra -pedantic -Wcast-align -Wdisabled-optimization -Wmissing-include-dirs -Wshadow -Wunused -fdiagnostics-show-option -Wstrict-overflow=0 -Wno-maybe-uninitialized -Werror #-Wrealloc-lhs
#FLG1 = -O0 -fbacktrace -fcheck=all -finit-real=nan -finit-derived -fmax-errors=10 --pedantic -Wall -Wextra -pedantic -Wcast-align -Wdisabled-optimization -Wmissing-include-dirs -Wshadow -Wunused -fdiagnostics-show-option -Wstrict-overflow=0 -Wno-maybe-uninitialized -Wno-compare-reals -Werror #-Wrealloc-lhs

#FFLAGS0 =       $(OPT) $(LPATH) $(LIBS) $(IPATH) $(MODULE_FLAG) $(MODULE_PATH) $(FLG0)
#FFLAGS1 =       $(OPT) $(LPATH) $(LIBS) $(IPATH) $(MODULE_FLAG) $(MODULE_PATH) $(FLG1)
#FFLAGS_LOW = $(OPTLOW) $(LPATH) $(LIBS) $(IPATH) $(MODULE_FLAG) $(MODULE_PATH) $(FLG0)

#LPATH = $(HDF5_LIB) $(QELIBS) 

LPATH = $(HDF5_LIB) $(SCALAPACK_LIBS) $(LAPACK_LIBS) $(FFT_LIBS) $(BLAS_LIBS) $(MPI_LIBS)
IPATH = $(IFLAGS_TDEP)  

FC = $(LD) $(F90FLAGS_TDEP) $(FDFLAGS)   
#$(FDFLAGS) $(MODFLAGS) $(IFLAGS)
FFLAGS0 = $(OPT) $(LPATH) $(IPATH) $(MODULE_FLAG) $(MODULE_PATH) 
FFLAGS1 = $(OPT) $(LPATH) $(IPATH) $(MODULE_FLAG) $(MODULE_PATH)
FFLAGS_LOW = $(OPTLOW) $(LPATH) $(IPATH) $(MODULE_FLAG) $(MODULE_PATH)


OBJSf = \
$(OBJECT_PATH)konstanter.o \
$(OBJECT_PATH)gottochblandat.o \
$(OBJECT_PATH)gottochblandat_tensors.o \
$(OBJECT_PATH)gottochblandat_sorting.o \
$(OBJECT_PATH)gottochblandat_linalg.o \
$(OBJECT_PATH)gottochblandat_calculus.o \
$(OBJECT_PATH)gottochblandat_physics.o \
$(OBJECT_PATH)gottochblandat_boxes.o \
$(OBJECT_PATH)lo_memtracker.o \
$(OBJECT_PATH)lo_sorting.o \
$(OBJECT_PATH)lo_verletboxes.o \
$(OBJECT_PATH)lo_voronoi.o \
$(OBJECT_PATH)lo_randomnumbers.o \
$(OBJECT_PATH)lo_brents_method.o \
$(OBJECT_PATH)lo_longrange_electrostatics.o \
$(OBJECT_PATH)lo_longrange_electrostatics_optz.o \
$(OBJECT_PATH)lo_longrange_electrostatics_dynmat.o \
$(OBJECT_PATH)quadratures_stencils.o \
$(OBJECT_PATH)geometryfunctions.o \
$(OBJECT_PATH)geometryfunctions_slicingroutines.o \
$(OBJECT_PATH)geometryfunctions_objectprocedures.o \
$(OBJECT_PATH)geometryfunctions_generatingfunctions.o \
$(OBJECT_PATH)cgal_wrappers.o \
$(OBJECT_PATH)hdf5_wrappers.o \
$(OBJECT_PATH)mpi_wrappers.o \
$(OBJECT_PATH)type_distancetable.o \
$(OBJECT_PATH)type_graph_distancetable.o \
$(OBJECT_PATH)type_voronoi.o \
$(OBJECT_PATH)type_voronoi_distancetable.o \
$(OBJECT_PATH)type_blas_lapack_wrappers.o \
$(OBJECT_PATH)type_crystalstructure.o \
$(OBJECT_PATH)type_qpointmesh.o \
$(OBJECT_PATH)type_qpointmesh_wedgegeneration.o \
$(OBJECT_PATH)type_qpointmesh_gridgeneration.o \
$(OBJECT_PATH)type_qpointmesh_integrationweights.o \
$(OBJECT_PATH)type_qpointmesh_bandstructure.o \
$(OBJECT_PATH)type_qpointmesh_commensurate.o \
$(OBJECT_PATH)type_qpointmesh_io.o \
$(OBJECT_PATH)type_symmetryoperation.o \
$(OBJECT_PATH)type_forceconstant_firstorder.o \
$(OBJECT_PATH)type_linkedlist.o \
$(OBJECT_PATH)type_jij_secondorder.o \
$(OBJECT_PATH)type_forceconstant_thirdorder.o \
$(OBJECT_PATH)type_forceconstant_fourthorder.o \
$(OBJECT_PATH)type_forceconstant_secondorder.o \


#$(OBJECT_PATH)type_forceconstant_secondorder_io.o \
#$(OBJECT_PATH)type_forceconstant_secondorder_dynamicalmatrix.o \
#$(OBJECT_PATH)type_forceconstant_secondorder_loto.o \
#$(OBJECT_PATH)type_forceconstant_secondorder_aux.o \
#$(OBJECT_PATH)type_forceconstant_secondorder.o \
#$(OBJECT_PATH)type_forceconstant_firstorder_io.o \
#$(OBJECT_PATH)type_forceconstant_secondorder_aux.o \
#$(OBJECT_PATH)type_forceconstant_secondorder_io.o \
#$(OBJECT_PATH)type_forceconstant_secondorder_loto.o \
#$(OBJECT_PATH)type_forceconstant_secondorder_dynamicalmatrix.o \

OBJS = $(OBJSf) $(OBJSa) $(OBJScg)

# rules to make the actual library
libolle.a: $(OBJS)
	-rm $(DEST_PATH)libolle.a
	ar rc $(DEST_PATH)libolle.a $(OBJS)

# The general stuff, does not depend on too complicated things
$(OBJECT_PATH)konstanter.o:
	$(FC) $(FFLAGS0) -c konstanter.f90 -o $@

$(OBJECT_PATH)lo_memtracker.o:\
    $(OBJECT_PATH)konstanter.o
	$(FC) $(FFLAGS0) -c lo_memtracker.f90 -o $@

$(OBJECT_PATH)lo_sorting.o:\
    $(OBJECT_PATH)konstanter.o\
    $(OBJECT_PATH)lo_memtracker.o
	$(FC) $(FFLAGS0) -c lo_sorting.f90 -o $@

$(OBJECT_PATH)lo_voronoi.o:\
    $(OBJECT_PATH)konstanter.o\
    $(OBJECT_PATH)lo_memtracker.o\
	$(OBJECT_PATH)gottochblandat.o\
	$(OBJECT_PATH)geometryfunctions.o\
	$(OBJECT_PATH)type_distancetable.o\
	$(OBJECT_PATH)mpi_wrappers.o
	$(FC) $(FFLAGS0) -c lo_voronoi.f90 -o $@

$(OBJECT_PATH)gottochblandat.o:\
	gottochblandat.f90 \
    $(OBJECT_PATH)konstanter.o
	$(FC) $(FFLAGS0) -c gottochblandat.f90 -o $@
$(OBJECT_PATH)gottochblandat_tensors.o:\
    $(OBJECT_PATH)gottochblandat.o\
    $(OBJECT_PATH)konstanter.o
	$(FC) $(FFLAGS0) -c gottochblandat_tensors.f90 -o $@
$(OBJECT_PATH)gottochblandat_sorting.o:\
    $(OBJECT_PATH)gottochblandat.o\
    $(OBJECT_PATH)konstanter.o
	$(FC) $(FFLAGS0) -c gottochblandat_sorting.f90 -o $@
$(OBJECT_PATH)gottochblandat_linalg.o:\
	gottochblandat_linalg.f90 \
    $(OBJECT_PATH)gottochblandat.o\
    $(OBJECT_PATH)konstanter.o
	$(FC) $(FFLAGS1) -c gottochblandat_linalg.f90 -o $@
$(OBJECT_PATH)gottochblandat_calculus.o:\
    $(OBJECT_PATH)gottochblandat.o\
    $(OBJECT_PATH)konstanter.o
	$(FC) $(FFLAGS0) -c gottochblandat_calculus.f90 -o $@
$(OBJECT_PATH)gottochblandat_physics.o:\
    $(OBJECT_PATH)gottochblandat.o\
    $(OBJECT_PATH)konstanter.o
	$(FC) $(FFLAGS0) -c gottochblandat_physics.f90 -o $@
$(OBJECT_PATH)gottochblandat_boxes.o:\
    $(OBJECT_PATH)gottochblandat.o\
    $(OBJECT_PATH)konstanter.o
	$(FC) $(FFLAGS0) -c gottochblandat_boxes.f90 -o $@

$(OBJECT_PATH)mpi_wrappers.o:\
    $(OBJECT_PATH)konstanter.o\
    $(OBJECT_PATH)gottochblandat.o
	$(FC) $(FFLAGS0) -c mpi_wrappers.f90 -o $@

$(OBJECT_PATH)lo_randomnumbers.o:\
    $(OBJECT_PATH)konstanter.o
	$(FC) $(FFLAGS0) -c lo_randomnumbers.f90 -o $@

$(OBJECT_PATH)lo_brents_method.o:\
    $(OBJECT_PATH)konstanter.o
	$(FC) $(FFLAGS0) -c lo_brents_method.f90 -o $@

$(OBJECT_PATH)lo_longrange_electrostatics.o:\
	lo_longrange_electrostatics.f90 \
    $(OBJECT_PATH)konstanter.o \
    $(OBJECT_PATH)mpi_wrappers.o \
    $(OBJECT_PATH)lo_memtracker.o \
	$(OBJECT_PATH)lo_brents_method.o \
	$(OBJECT_PATH)type_crystalstructure.o \
    $(OBJECT_PATH)gottochblandat.o
	$(FC) $(FFLAGS0) -c lo_longrange_electrostatics.f90 -o $@
$(OBJECT_PATH)lo_longrange_electrostatics_optz.o:\
	lo_longrange_electrostatics_optz.f90 \
    $(OBJECT_PATH)lo_longrange_electrostatics.o \
	$(OBJECT_PATH)quadratures_stencils.o
	$(FC) $(FFLAGS0) -c lo_longrange_electrostatics_optz.f90 -o $@
$(OBJECT_PATH)lo_longrange_electrostatics_dynmat.o:\
	lo_longrange_electrostatics_dynmat.f90 \
    $(OBJECT_PATH)lo_longrange_electrostatics.o
	$(FC) $(FFLAGS0) -c lo_longrange_electrostatics_dynmat.f90 -o $@

$(OBJECT_PATH)quadratures_stencils.o:\
    $(OBJECT_PATH)konstanter.o\
    $(OBJECT_PATH)gottochblandat.o
	$(FC) $(FFLAGS0) -c quadratures_stencils.f90 -o $@

$(OBJECT_PATH)type_linkedlist.o:\
    $(OBJECT_PATH)konstanter.o
	$(FC) $(FFLAGS_LOW) -c type_linkedlist.f90 -o $@

$(OBJECT_PATH)geometryfunctions.o: \
    $(OBJECT_PATH)konstanter.o\
    $(OBJECT_PATH)gottochblandat.o
	$(FC) $(FFLAGS0) -c geometryfunctions.f90 -o $@
$(OBJECT_PATH)geometryfunctions_slicingroutines.o: \
    $(OBJECT_PATH)geometryfunctions.o\
    $(OBJECT_PATH)konstanter.o\
    $(OBJECT_PATH)gottochblandat.o
	$(FC) $(FFLAGS_LOW) -c geometryfunctions_slicingroutines.f90 -o $@
$(OBJECT_PATH)geometryfunctions_objectprocedures.o: \
    $(OBJECT_PATH)geometryfunctions.o\
    $(OBJECT_PATH)konstanter.o\
    $(OBJECT_PATH)gottochblandat.o
	$(FC) $(FFLAGS0) -c geometryfunctions_objectprocedures.f90 -o $@
$(OBJECT_PATH)geometryfunctions_generatingfunctions.o: \
    $(OBJECT_PATH)geometryfunctions.o\
    $(OBJECT_PATH)konstanter.o\
    $(OBJECT_PATH)gottochblandat.o
	$(FC) $(FFLAGS0) -c geometryfunctions_generatingfunctions.f90 -o $@

$(OBJECT_PATH)hdf5_wrappers.o:\
    $(OBJECT_PATH)gottochblandat.o\
    $(OBJECT_PATH)konstanter.o
	$(FC) $(FFLAGS0) -O0 -c hdf5_wrappers.f90 -o $@

$(OBJECT_PATH)lo_verletboxes.o:\
    $(OBJECT_PATH)konstanter.o \
    $(OBJECT_PATH)lo_memtracker.o \
    $(OBJECT_PATH)mpi_wrappers.o
	$(FC) $(FFLAGS0) -c lo_verletboxes.f90 -o $@

$(OBJECT_PATH)type_blas_lapack_wrappers.o:\
    $(OBJECT_PATH)konstanter.o
	$(FC) $(FFLAGS0) -c type_blas_lapack_wrappers.f90 -o $@

$(OBJECT_PATH)type_distancetable.o:\
    $(OBJECT_PATH)konstanter.o\
    $(OBJECT_PATH)gottochblandat.o\
    $(OBJECT_PATH)type_blas_lapack_wrappers.o\
    $(OBJECT_PATH)geometryfunctions.o\
    $(OBJECT_PATH)mpi_wrappers.o
	$(FC) $(FFLAGS0) -c type_distancetable.f90 -o $@

$(OBJECT_PATH)type_voronoi.o:\
    $(OBJECT_PATH)konstanter.o\
    $(OBJECT_PATH)gottochblandat.o\
    $(OBJECT_PATH)geometryfunctions.o\
    $(OBJECT_PATH)type_distancetable.o\
    $(OBJECT_PATH)type_linkedlist.o
	$(FC) $(FFLAGS_LOW) -c type_voronoi.f90 -o $@
$(OBJECT_PATH)type_graph_distancetable.o:\
    $(OBJECT_PATH)konstanter.o\
    $(OBJECT_PATH)gottochblandat.o\
    $(OBJECT_PATH)geometryfunctions.o\
    $(OBJECT_PATH)type_voronoi.o\
    $(OBJECT_PATH)type_distancetable.o
	$(FC) $(FFLAGS0) -c type_graph_distancetable.f90 -o $@
$(OBJECT_PATH)type_voronoi_distancetable.o:\
    $(OBJECT_PATH)konstanter.o\
    $(OBJECT_PATH)gottochblandat.o\
    $(OBJECT_PATH)geometryfunctions.o\
    $(OBJECT_PATH)type_distancetable.o\
    $(OBJECT_PATH)type_voronoi.o
	$(FC) $(FFLAGS0) -c type_voronoi_distancetable.f90 -o $@
$(OBJECT_PATH)type_symmetryoperation.o:\
    $(OBJECT_PATH)konstanter.o\
    $(OBJECT_PATH)gottochblandat.o\
    $(OBJECT_PATH)geometryfunctions.o
	$(FC) $(FFLAGS0) -c type_symmetryoperation.f90 -o $@
#
# If I have CGAL installed
#
ifeq ($(USECGAL),yes)
$(OBJECT_PATH)cgal_wrappers.o:\
    $(OBJECT_PATH)konstanter.o\
    $(OBJECT_PATH)gottochblandat.o\
    $(OBJECT_PATH)cgal_cleanup.o\
    $(OBJECT_PATH)cgal_chull2.o\
    $(OBJECT_PATH)cgal_chull3.o\
    $(OBJECT_PATH)cgal_deltri2.o\
    $(OBJECT_PATH)cgal_deltri3.o\
    $(OBJECT_PATH)cgal_tesselate_polyhedron.o
	$(FC) $(FFLAGS0) -c cgal_wrappers.f90 -o $@
else
$(OBJECT_PATH)cgal_wrappers.o:\
    $(OBJECT_PATH)konstanter.o\
    $(OBJECT_PATH)gottochblandat.o
	$(FC) $(FFLAGS0) -c cgal_wrappers.f90 -o $@
endif
#
# Sockets for I-PI
#
#$(OBJECT_PATH)type_ipi.o:\
#    $(OBJECT_PATH)type_ipi_sockets.o
#	$(FC) $(FFLAGS) -c type_ipi.f90 -o $@
#
# A bit more advanced
#
$(OBJECT_PATH)type_crystalstructure.o:\
	type_crystalstructure.f90\
    $(OBJECT_PATH)konstanter.o\
    $(OBJECT_PATH)hdf5_wrappers.o\
    $(OBJECT_PATH)gottochblandat.o\
    $(OBJECT_PATH)geometryfunctions.o\
    $(OBJECT_PATH)type_voronoi.o\
    $(OBJECT_PATH)type_distancetable.o\
    $(OBJECT_PATH)type_symmetryoperation.o
	$(FC) $(FFLAGS0) -c type_crystalstructure.f90 -o $@

$(OBJECT_PATH)type_qpointmesh.o:\
    $(OBJECT_PATH)konstanter.o\
    $(OBJECT_PATH)mpi_wrappers.o\
    $(OBJECT_PATH)lo_memtracker.o\
    $(OBJECT_PATH)gottochblandat.o\
	$(OBJECT_PATH)geometryfunctions.o\
    $(OBJECT_PATH)type_crystalstructure.o\
    $(OBJECT_PATH)type_symmetryoperation.o
	$(FC) $(FFLAGS0) -c type_qpointmesh.f90 -o $@
$(OBJECT_PATH)type_qpointmesh_wedgegeneration.o:\
    $(OBJECT_PATH)cgal_wrappers.o\
    $(OBJECT_PATH)lo_sorting.o\
    $(OBJECT_PATH)lo_verletboxes.o\
    $(OBJECT_PATH)type_qpointmesh.o
	$(FC) $(FFLAGS0) -c type_qpointmesh_wedgegeneration.f90 -o $@
$(OBJECT_PATH)type_qpointmesh_gridgeneration.o:\
    type_qpointmesh_gridgeneration.f90\
    $(OBJECT_PATH)type_qpointmesh.o
	$(FC) $(FFLAGS0) -c type_qpointmesh_gridgeneration.f90 -o $@
$(OBJECT_PATH)type_qpointmesh_integrationweights.o:\
    $(OBJECT_PATH)lo_voronoi.o\
	$(OBJECT_PATH)lo_sorting.o\
	$(OBJECT_PATH)lo_verletboxes.o\
	$(OBJECT_PATH)type_distancetable.o\
    $(OBJECT_PATH)type_qpointmesh.o
	$(FC) $(FFLAGS0) -c type_qpointmesh_integrationweights.f90 -o $@
$(OBJECT_PATH)type_qpointmesh_bandstructure.o:\
    type_qpointmesh_bandstructure.f90\
	$(OBJECT_PATH)lo_sorting.o\
    $(OBJECT_PATH)type_qpointmesh.o
	$(FC) $(FFLAGS0) -c type_qpointmesh_bandstructure.f90 -o $@
$(OBJECT_PATH)type_qpointmesh_commensurate.o:\
    $(OBJECT_PATH)type_qpointmesh.o
	$(FC) $(FFLAGS0) -c type_qpointmesh_commensurate.f90 -o $@
$(OBJECT_PATH)type_qpointmesh_io.o:\
    $(OBJECT_PATH)hdf5_wrappers.o\
    $(OBJECT_PATH)type_qpointmesh.o
	$(FC) $(FFLAGS_LOW) -c type_qpointmesh_io.f90 -o $@

$(OBJECT_PATH)type_forceconstant_firstorder.o:\
	type_forceconstant_firstorder.f90\
    $(OBJECT_PATH)konstanter.o\
    $(OBJECT_PATH)type_crystalstructure.o
	$(FC) $(FFLAGS0) -c type_forceconstant_firstorder.f90 -o $@


$(OBJECT_PATH)type_forceconstant_secondorder.o:\
	type_forceconstant_secondorder.f90\
    $(OBJECT_PATH)konstanter.o\
    $(OBJECT_PATH)gottochblandat.o\
    $(OBJECT_PATH)mpi_wrappers.o\
	$(OBJECT_PATH)lo_memtracker.o\
	$(OBJECT_PATH)lo_longrange_electrostatics.o\
    $(OBJECT_PATH)type_crystalstructure.o\
    $(OBJECT_PATH)type_qpointmesh.o
	$(FC) $(FFLAGS0) -c type_forceconstant_secondorder.f90 -o $@
$(OBJECT_PATH)type_forceconstant_secondorder_dynamicalmatrix.o:\
	type_forceconstant_secondorder_dynamicalmatrix.f90\
    $(OBJECT_PATH)type_forceconstant_secondorder.o\
    $(OBJECT_PATH)konstanter.o\
    $(OBJECT_PATH)gottochblandat.o\
	$(OBJECT_PATH)lo_sorting.o\
    $(OBJECT_PATH)type_blas_lapack_wrappers.o
	$(FC) $(FFLAGS0) -c type_forceconstant_secondorder_dynamicalmatrix.f90 -o $@
$(OBJECT_PATH)type_forceconstant_secondorder_io.o:\
	type_forceconstant_secondorder_io.f90\
    $(OBJECT_PATH)type_forceconstant_secondorder.o
	$(FC) $(FFLAGS0) -c type_forceconstant_secondorder_io.f90 -o $@
$(OBJECT_PATH)type_forceconstant_secondorder_loto.o:\
	type_forceconstant_secondorder_loto.f90\
    $(OBJECT_PATH)type_forceconstant_secondorder.o
	$(FC) $(FFLAGS0) -c type_forceconstant_secondorder_loto.f90 -o $@


#$(OBJECT_PATH)type_forceconstant_secondorder_aux.o:\
#	type_forceconstant_secondorder_aux.f90\
#    $(OBJECT_PATH)type_forceconstant_secondorder.o\
#    $(OBJECT_PATH)type_distancetable.o\
#    $(OBJECT_PATH)lo_randomnumbers.o\
#	$(OBJECT_PATH)lo_sorting.o\
#    $(OBJECT_PATH)type_blas_lapack_wrappers.o
#	$(FC) $(FFLAGS0) -c type_forceconstant_secondorder_aux.f90 -o $@
#$(OBJECT_PATH)type_forceconstant_secondorder.o:\
#	type_forceconstant_secondorder.f90\
#    $(OBJECT_PATH)konstanter.o\
#    $(OBJECT_PATH)gottochblandat.o\
#    $(OBJECT_PATH)mpi_wrappers.o\
#	$(OBJECT_PATH)lo_memtracker.o\
#	$(OBJECT_PATH)lo_longrange_electrostatics.o\
#    $(OBJECT_PATH)type_crystalstructure.o\
#    $(OBJECT_PATH)type_qpointmesh.o\
#    $(OBJECT_PATH)type_distancetable.o\
#    $(OBJECT_PATH)lo_randomnumbers.o\
#	$(OBJECT_PATH)lo_sorting.o\
#    $(OBJECT_PATH)type_blas_lapack_wrappers.o
#	$(FC) $(FFLAGS0) -c type_forceconstant_secondorder.f90 -o $@

$(OBJECT_PATH)type_forceconstant_thirdorder.o:\
    type_forceconstant_thirdorder.f90\
    $(OBJECT_PATH)konstanter.o\
    $(OBJECT_PATH)gottochblandat.o\
    $(OBJECT_PATH)type_crystalstructure.o\
    $(OBJECT_PATH)type_qpointmesh.o\
    $(OBJECT_PATH)type_distancetable.o
	$(FC) $(FFLAGS0) -c type_forceconstant_thirdorder.f90 -o $@
$(OBJECT_PATH)type_forceconstant_fourthorder.o:\
    $(OBJECT_PATH)konstanter.o\
    $(OBJECT_PATH)gottochblandat.o\
    $(OBJECT_PATH)type_crystalstructure.o\
    $(OBJECT_PATH)type_distancetable.o
	$(FC) $(FFLAGS0) -c type_forceconstant_fourthorder.f90 -o $@
$(OBJECT_PATH)type_jij_secondorder.o:\
    $(OBJECT_PATH)konstanter.o\
    $(OBJECT_PATH)gottochblandat.o\
    $(OBJECT_PATH)type_distancetable.o\
    $(OBJECT_PATH)type_crystalstructure.o
	$(FC) $(FFLAGS0) -c type_jij_secondorder.f90 -o $@



clean:
	rm -rf $(OBJECT_PATH)*.o $(MODULE_PATH)*mod 
