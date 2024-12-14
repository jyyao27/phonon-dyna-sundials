subroutine phonon_gruneisen()
   use pert_const, only: dp
   use qe_mpi_mod, only: ionode, stdout
   use pert_output, only: output_gruneisen
   use pert_data,  only:  epr_fid, kc_dim, qc_dim
   use vector_list, only: vlist, load_vector_list, rand_vector_list
   use pert_param, only: fklist, fqlist, ntemper, temper, efermi, &
      cauchy_scale, sampling, nsamples
   use phonon_anharm_prop, only: calc_gruneisen
   !
   use phonon_dispersion, only: lattice_ifc, lattice_ifct, &
                 init_lattice_ifc, init_lattice_ifct
   implicit none
   !integer :: ntpr !ntpr defines phonon distribution
   integer :: i
   !load k- and q- list
   type(vlist) :: ql
   real(dp), allocatable:: wq(:,:), grun(:,:)
   !
   type(lattice_ifc)   :: phon
   type(lattice_ifct)  :: phont
   !check if temperatures have been readin 
   !if(ntemper .eq. 0) call errore('phonon_gruneisen',&
   !'ftemper is not specified, missing temperatures and chemical potential',1)
   !ntpr = size(temper)

   !load k- and q-grid
   call load_vector_list(fqlist, ql)
   if(ionode) then
      write(stdout,'(5x,a,i11)') 'total number of k-points:', ql%nvec
   endif
   !init
   call init_lattice_ifc(epr_fid, qc_dim, phon)
   call init_lattice_ifct(epr_fid, qc_dim, phont)

   allocate(grun(phon%nm, ql%nvec), wq(phon%nm, ql%nvec))
   grun = 0.0_dp; wq = 0.0_dp   

   call calc_gruneisen(phon, phont, ql%vec, wq, grun)

   !output results
   if(ionode) then
      call output_gruneisen(wq, grun)
   endif
   !release memory
   deallocate(wq, grun)
   
   if(ionode) write(stdout,'(5x, a)') '>Gruneisen done.'
end subroutine phonon_gruneisen
