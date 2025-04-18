
module third_order_fc  
   use hdf5_utils
   use pert_const, only: dp, ci, twopi 
   implicit none 
   type, public :: fct_matrix 
      !
      integer, public :: nmod, nr 
      !
      integer, public :: nsize 
      !
      integer, public, allocatable :: im_index(:)
      !
      real(dp), public, allocatable :: rset(:,:)
      !rset(:,ir(i)) = the real space vector of the i-th entry 
      integer, public, allocatable :: ir(:)
      !
      real(dp), public, allocatable :: phi3_matrix(:,:,:)
      !
      real(dp), public, allocatable :: phi3_matrix_dou(:,:,:), phi3_matrix_tri(:,:,:)
   end type 

   contains 

   subroutine allocate_fct(fct)
      type(fct_matrix) :: fct 
      allocate(fct%rset(3, fct%nr), fct%ir(fct%nsize), fct%im_index(fct%nsize), fct%phi3_matrix(fct%nsize, fct%nsize, fct%nmod))
      allocate(fct%phi3_matrix_tri(fct%nsize, fct%nsize, fct%nmod))
      allocate(fct%phi3_matrix_dou(fct%nsize, fct%nsize, fct%nmod))
   end subroutine 

   subroutine deallocate_fct(fct)
      type(fct_matrix) :: fct 
      deallocate(fct%rset, fct%ir, fct%im_index, fct%phi3_matrix, fct%phi3_matrix_dou, fct%phi3_matrix_tri)
   end subroutine 

   subroutine load_fct(prefix, fct) 
      character(len=80) :: prefix
      type(fct_matrix) :: fct 
      !local 
      integer(HID_T) :: group_id, fct_id 
      character(len=120) :: dset_name,fname
      
      !read to hdf5: prefix_fct.h5  
      fname = trim(prefix)//"_fct.h5"
      call hdf_open_file(fct_id, trim(fname), status='OLD', action='READ')
      call hdf_open_group(fct_id, "fct", group_id)
      dset_name = "nsize"
      call hdf_read_dataset(group_id, trim(dset_name), fct%nsize)
      dset_name = "nr"
      call hdf_read_dataset(group_id, trim(dset_name), fct%nr)
      dset_name = "nmod"
      call hdf_read_dataset(group_id, trim(dset_name), fct%nmod)
      !
      call allocate_fct(fct)
      dset_name = "mode_index"
      call hdf_read_dataset(group_id, trim(dset_name), fct%im_index)
      dset_name = "ir"
      call hdf_read_dataset(group_id, trim(dset_name), fct%ir)   
      dset_name = "rset"
      call hdf_read_dataset(group_id, trim(dset_name), fct%rset)   
      dset_name = "phi3_matrix"
      call hdf_read_dataset(group_id, trim(dset_name), fct%phi3_matrix)
      dset_name = "phi3_matrix_2"
      call hdf_read_dataset(group_id, trim(dset_name), fct%phi3_matrix_dou)
      dset_name = "phi3_matrix_3"
      call hdf_read_dataset(group_id, trim(dset_name), fct%phi3_matrix_tri)
      call hdf_close_group(group_id)
      !
      call hdf_close_file(fct_id)

      !write(*,'(A20,3i10)') 'nsize, nr, nm = ',fct%nsize, fct%nr, fct%nmod 
      !write(*,'(i20)') fct%ir
      !write(*,'(E20.10)') maxval(abs(fct%phi3_matrix))
   end subroutine

   subroutine write_fct(prefix, fct)
      character(len=80) :: prefix
      type(fct_matrix) :: fct 
      !local 
      integer(HID_T) :: group_id, fct_id 
      character(len=120) :: dset_name,fname

      !write to hdf5: prefix_fct.h5  
      fname = trim(prefix)//"_fct.h5"
      call hdf_open_file(fct_id, trim(fname), status='NEW', action='WRITE')
      call hdf_create_group(fct_id, 'fct')
      call hdf_open_group(fct_id, "fct", group_id)
      dset_name = "nsize"
      call hdf_write_dataset(group_id, trim(dset_name), fct%nsize)
      dset_name = "nr"
      call hdf_write_dataset(group_id, trim(dset_name), fct%nr)
      dset_name = "nmod"
      call hdf_write_dataset(group_id, trim(dset_name), fct%nmod)
      dset_name = "mode_index"
      call hdf_write_dataset(group_id, trim(dset_name), fct%im_index)
      dset_name = "ir"
      call hdf_write_dataset(group_id, trim(dset_name), fct%ir) 
      dset_name = "rset"
      call hdf_write_dataset(group_id, trim(dset_name), fct%rset)   
      dset_name = "phi3_matrix"
      call hdf_write_dataset(group_id, trim(dset_name), fct%phi3_matrix)
      dset_name = "phi3_matrix_3"
      call hdf_write_dataset(group_id, trim(dset_name), fct%phi3_matrix_tri)
      dset_name = "phi3_matrix_2"
      call hdf_write_dataset(group_id, trim(dset_name), fct%phi3_matrix_dou)
      call hdf_close_group(group_id)
      !
      call hdf_close_file(fct_id)
   end subroutine 

   subroutine body_order_decomp(fct)
      implicit none 
      type(fct_matrix) :: fct 
      integer :: ia, ia1, ia2, ir, i,j, im, im1, im2  
      integer :: nsize, nmod 
      real(dp) :: r1(3), r2(3) 
      !local 
      integer(HID_T) :: group_id, fct_id 
      character(len=120) :: dset_name,fname
      
      nsize = fct%nsize 
      nmod = fct%nmod 

      fct%phi3_matrix_tri = fct%phi3_matrix
      do im = 1,nmod
         ia = (im-1) / 3  
         do i = 1, nsize 
            im1 = fct%im_index(i)
            ia1 = (im1-1)/3
            r1 = fct%rset(:,fct%ir(i))
            do j = 1,nsize
               im2 = fct%im_index(j)
               ia2 = (im2-1)/3
               r2 = fct%rset(:, fct%ir(j))
               if(ia.eq.ia1 .and. maxval(abs(r1))<=1e-5 ) fct%phi3_matrix_tri(i,j,im) = 0.d0 
               if(ia.eq.ia2 .and. maxval(abs(r2))<=1e-5 ) fct%phi3_matrix_tri(i,j,im) = 0.d0 
               if(ia1.eq.ia2 .and. maxval(abs(r1-r2))<=1e-5 ) fct%phi3_matrix_tri(i,j,im) = 0.d0 
            enddo
         enddo  
      enddo 
      !one and two body term 
      fct%phi3_matrix_dou = fct%phi3_matrix - fct%phi3_matrix_tri

   end subroutine 

   !phi(r2=0 im2; r1 im1; q im) = \sum_{r} phirr(r1,r) e^{iqr}
   subroutine phirr_to_phiqr(fct, qpt, phiqr)
      implicit none 
      type(fct_matrix) :: fct 
      real(dp) :: qpt(3)
      complex(dp) :: phiqr(fct%nr, fct%nmod, fct%nmod, fct%nmod)
      !local 
      integer :: il, jl, im2, im1, im, ir2, ir1, ir  
      real(dp) :: r1(3)
      complex(dp) :: exp_qr(fct%nsize)
      do il = 1,fct%nsize  
         exp_qr(il) = exp(ci*sum(qpt*fct%rset(:,fct%ir(il)))*twopi)
      enddo 

      phiqr = 0.d0 
      do im2 = 1,fct%nmod 
         do il = 1,fct%nsize
            im1 = fct%im_index(il) 
            ir1 = fct%ir(il)
            do jl = 1,fct%nsize 
               im = fct%im_index(jl) 
               phiqr(ir1, im2, im1, im) = phiqr(ir1, im2, im1, im) + exp_qr(jl) * fct%phi3_matrix(il, jl, im2)
            enddo 
         enddo 
      enddo 
   end subroutine

   !phi(q2=(-q1-q2) im2; q1 im1; q im) = \sum_{r1} phirr(r1,q) e^{iq_1r_1}
   subroutine phiqr_to_phiqq(fct, qpt, phiqr, phiqq)
      implicit none 
      type(fct_matrix) :: fct 
      complex(dp) :: phiqr(fct%nr, fct%nmod, fct%nmod, fct%nmod)
      complex(dp) :: phiqq(fct%nmod, fct%nmod, fct%nmod)
      
      real(dp) :: qpt(3)

      !local 
      integer ::  im2, im1, im, ir2, ir1, ir  
      real(dp) :: r1(3)
      complex(dp) :: exp_qr(fct%nr)
      do ir1 = 1,fct%nr  
         exp_qr(ir1) = exp( ci*sum(qpt * fct%rset(:,ir1))*twopi )
      enddo 

      do im2 = 1,fct%nmod 
         do im1 = 1,fct%nmod
            do im = 1,fct%nmod 
               phiqq(im2,im1,im) = sum(exp_qr*phiqr(:,im2,im1,im) )
            enddo 
         enddo 
      enddo 

   end subroutine

   !transform ph-ph matrix from cartesian to eigenmode
   ! \sum_{la}u_1(la, \nu_1)\sum_{l1b}u_2(l1b, \nu_2)\sum_{l2c} u_2(l2c, \nu_2) phi(la,l1b,l2c)
   subroutine phi3_transform(nmod, u1, u2, u3, phi3_qq)
      implicit none 
      integer,intent(in) :: nmod 
      complex(dp) :: phi3_qq(nmod, nmod, nmod)
      complex(dp) :: u1(nmod,nmod), u2(nmod,nmod), u3(nmod,nmod),xt 
      !
      integer :: im1, im2, im3, im1p, im2p, im3p 
      complex(dp) :: phi3_tmp(nmod, nmod, nmod)
      
      phi3_tmp = 0.d0 
      do im1=1,nmod; do im2=1,nmod; do im3=1,nmod
         xt = 0.d0 
         do im1p=1,nmod; do im2p=1,nmod; do im3p=1,nmod
            xt = xt + phi3_qq(im1p, im2p, im3p) * u1(im1p, im1) * u2(im2p, im2) * u3(im3p, im3)
         enddo;enddo;enddo 
         phi3_tmp(im1,im2,im3) = xt 
      enddo;enddo;enddo 
      phi3_qq = phi3_tmp
      return 
   end subroutine

end module 
