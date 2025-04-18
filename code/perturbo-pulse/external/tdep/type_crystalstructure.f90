#include "precompilerdefinitions"
module type_crystalstructure
!! Some text about the crystalstructure module
use konstanter, only: flyt,i8,lo_huge,lo_tiny,lo_pi,lo_hugeint,lo_tol,lo_sqtol,lo_radiantol,lo_status,lo_bohr_to_A,&
                      lo_A_to_bohr,lo_velocity_Afs_to_au,lo_amu_to_emu,lo_velocity_au_to_Afs,lo_exitcode_param,&
                      lo_exitcode_symmetry,lo_exitcode_io
use gottochblandat, only: open_file,tochar,qsort,walltime,lo_clean_fractional_coordinates,lo_chop,lo_determ,lo_choplarge,&
                   lo_frobnorm,lo_sqnorm,lo_reciprocal_basis,lo_mean,lo_unsigned_tetrahedron_volume,&
                   lo_index_in_periodic_array,lo_cross,lo_get_axis_angles,lo_invert3x3matrix,lo_stop_gracefully,&
                   lo_permutations,lo_return_unique
use geometryfunctions, only: lo_plane,lo_inscribed_sphere_in_box,lo_bounding_sphere_of_box
use type_distancetable, only: lo_distancetable
use type_voronoi, only: lo_voronoi_diagram,lo_voronoi_diagram_cell
use type_symmetryoperation, only: lo_symset, lo_operate_on_vector
use hdf5

implicit none
private
public :: lo_crystalstructure

!> define an atom with partial occupancy of different species, as in a random alloy
type lo_alloyatom
    !> how many components
    integer :: n=-lo_hugeint
    !> what are the components
    character(len=40), dimension(:), allocatable :: atomic_symbol
    !> concentration
    real(flyt), dimension(:), allocatable :: concentration
    !> atomic number of component
    integer, dimension(:), allocatable :: atomic_number
    !> mass of the components
    real(flyt), dimension(:), allocatable :: mass
end type

!> some shorthand to keep track of how disordered magnetic states are to be created
type lo_maginfo
    !> keep track of which atoms have magnetic moments
    logical, dimension(:), allocatable :: atom_has_moment
    !> collinear moment on each atom
    real(flyt), dimension(:), allocatable :: collinear_moment
    !> noncollinear moment on each atom
    real(flyt), dimension(:,:), allocatable :: noncollinear_moment
    !> species list is a little bit different in case of magnetism
    integer, dimension(:), allocatable :: magspecies
    !> number of (magnetic) species
    integer :: nmagspecies=-lo_hugeint
    !> counter for the magnetic species
    integer, dimension(:), allocatable :: magspeciescounter
end type

!> Distribution of isotopes for an atom.
type lo_isotope_distribution
    !> how many isotopes
    integer :: n=-lo_hugeint
    !> concentration of isotopes
    real(flyt), dimension(:), allocatable :: conc
    !> mass of isotope
    real(flyt), dimension(:), allocatable :: mass
    !> average mass
    real(flyt) :: mean_mass=lo_huge
    !> disorderparameter
    real(flyt) :: disorderparameter=lo_huge
    contains
        !> fetch the natural distribution
        procedure :: naturaldistribution
        !> calculate the mass disorder parameter
        procedure :: mass_disorder_parameter
end type

!> Brilluoin zone, stored as a polyhedron.
type, extends(lo_voronoi_diagram_cell) :: lo_brillouin_zone
    !> how many high symmetry points
    integer :: nhighsymmetrypoints=-lo_hugeint
    !> high-symmetry points, all of them
    real(flyt), dimension(:,:), allocatable :: highsymmetrypoints
    !> labels for the high-symmetry points
    character(len=10), dimension(:), allocatable :: label
    contains
        !> reciprocal lattice vector that shifts a point to the first bz
        procedure :: gshift
        !> calculates the distance to the zone edge.
        procedure :: distance_to_zone_edge
end type

!> A face in the wedge in the Brilluoin zone.
type lo_brillouin_zone_wedge_face
    !> how many points on this face?
    integer :: n=-lo_hugeint
    !> which nodes make up this face
    integer, dimension(:), allocatable :: ind
    !> what plane defines this face
    type(lo_plane) :: plane
end type

!> Irreducible wedge in the Brilluoin zone, polyhedron in reciprocal space.
type :: lo_brillouin_zone_irreducible_wedge
    !> how many irreducible points
    integer :: nnodes=-lo_hugeint
    !> how many faces in the polyhedron
    integer :: nfaces=-lo_hugeint
    !> faces
    type(lo_brillouin_zone_wedge_face), dimension(:), allocatable :: face
    !> coordinates to the high symmetry points
    real(flyt), dimension(:,:), allocatable :: r
    !> the labels of these points
    character(len=10), dimension(:), allocatable :: label
end type

!> Subcontainer with classification, none of which is particularly important.
type lo_crystalstructure_classification
    !> Just a title for this crystalstructure
    character(len=1000) :: title='empty header'
    ! Heuristics

    !> Is this an alloy?
    logical :: alloy=.false.
    !> have I specified collinear magnetic moments?
    logical :: collmag=.false.
    !> have I specified noncollinear magnetic moments?
    logical :: noncollmag=.false.
    !> have I calculated the symmetry operations?
    logical :: havespacegroup=.false.
    !> have I calculated the Brillouin zone
    logical :: havebz=.false.
    !> have I calculated the irreducible wedge
    logical :: havewedge=.false.
    !> have I figured out which Bravais lattice it is?
    logical :: havebravais=.false.
    !> Are all high symmetry points labelled ok?
    logical :: pointslabelled=.false.
    !> Have I deced on time-reversal symmetry?
    logical :: decidedtimereversal=.false.
    !> Have I calculated the character table?
    logical :: havecharactertable=.false.

    ! Classification things

    !> unique representation of the primitive lattice
    real(flyt), dimension(3,3) :: unique_primitive_basis=lo_huge
    !> unique representation of the conventional lattice
    real(flyt), dimension(3,3) :: unique_conventional_basis=lo_huge
    !> Permutation matrix from the basis to the unique primitive basis
    real(flyt), dimension(3,3) :: permutation_to_unique=lo_huge
    !> Transformation from the basis to the unique primitive basis
    real(flyt), dimension(3,3) :: transformation_to_unique=lo_huge
    !> the "a" lattice parameter, for pretty output.
    real(flyt) :: unitcell_lattice_parameter=lo_huge
    !> What Bravais lattice is it?
    character(len=10) :: bravaislattice='nothing'
    character(len=50) :: bravaislatticelongname='nothing'

    ! Supercell stuff

    !> Is this a supercell?
    logical :: supercell=.false.
    !> What are the dimensions?
    integer, dimension(3,3) :: supercellmatrix=-lo_hugeint
    !> For Fourier transforms and things like that it's good to know the index to the unitcell in the supercell
    integer, dimension(:,:), allocatable :: cellindex
    !> Might also be nice to have an index in the unit cell that it can be related to
    integer, dimension(:), allocatable :: index_in_unitcell
    !> How much should it talk?
    integer :: verbosity=-lo_hugeint
end type

!> Crystal structure class. Contains where the atoms are, how many and so on, basically everything about the crystal structure. When initialized by reading a structure from file, a whole bunch of stuff is calculated and gathered from tables, such as atomic masses, atomic numbers, space group, Brillouin zone and so on.
type lo_crystalstructure
    !> Isotope distributions
    type(lo_isotope_distribution), dimension(:), allocatable :: isotope
    !> Info about a species defined as a random alloy
    type(lo_alloyatom), dimension(:), allocatable :: alloyspecies
    !> Brillouin zone
    type(lo_brillouin_zone) :: bz
    !> irreducible wedge of the Brillouin zone
    type(lo_brillouin_zone_irreducible_wedge) :: irrw
    !> The space group for this structure
    type(lo_symset) :: sym
    !> Some extra classification information, if that is useful
    type(lo_crystalstructure_classification) :: info
    !> Information about the magnetic state
    type(lo_maginfo) :: mag
    !> number of atoms in the cell
    integer :: na=-lo_hugeint
    !> basis vectors
    real(flyt), dimension(3,3) :: latticevectors=lo_huge
    !> inverse of basis vectors
    real(flyt), dimension(3,3) :: inv_latticevectors=lo_huge
    !> reciprocal lattice vectors
    real(flyt), dimension(3,3) :: reciprocal_latticevectors=lo_huge
    !> inverse reciprocal lattice vectors
    real(flyt), dimension(3,3) :: inv_reciprocal_latticevectors=lo_huge
    !> Volume of the cell
    real(flyt) :: volume=lo_huge
    !> How many different elements
    integer :: nelements=-lo_hugeint
    !> Counter for each element type
    integer, dimension(:), allocatable :: element_counter
    !> What are the atomic symbols of these elements, e.g. Li,Fe,O
    character(len=8), dimension(:), allocatable :: atomic_symbol
    !> Atomic number for each atom
    integer, dimension(:), allocatable :: atomic_number
    !> What species is atom i?
    integer, allocatable, dimension(:) :: species
    !> Positions, fractional coordinates
    real(flyt), allocatable, dimension(:,:) :: r
    !> Positions, cartesian coordinates
    real(flyt), allocatable, dimension(:,:) :: rcart
    !> Velocities
    real(flyt), allocatable, dimension(:,:) :: v
    !> Forces
    real(flyt), allocatable, dimension(:,:) :: f
    !> Displacements
    real(flyt), allocatable, dimension(:,:) :: u
    !> Mass of each atom
    real(flyt), allocatable, dimension(:) :: mass
    !> Inverse square root of mass of each atom
    real(flyt), allocatable, dimension(:) :: invsqrtmass
    !> The inelastic neutron cross-section
    real(flyt), allocatable, dimension(:) :: inelastic_neutron_cross_section
    !> Flavor on each atom, that is to be preserved
    integer, allocatable, dimension(:) :: flavor
    contains
        !> create the structure
        procedure :: generate
        !> classify it in different ways
        procedure :: classify
        !> coordinate conversion
        procedure :: fractional_to_cartesian
        procedure :: cartesian_to_fractional
        procedure :: displacement_fractional_to_cartesian
        !> get the kinetic energy from velocities
        procedure :: kinetic_energy
        !> build a supercell
        procedure :: build_supercell
        !> largest pair cutoff in cell
        procedure :: maxcutoff
        !> smallest pair cutoff in cell
        procedure :: mincutoff
        !> nearest neighbour distance
        procedure :: nearest_neighbour_distance
        !> convert a coordinate to a high symmetry point
        procedure :: coordinate_from_high_symmetry_point_label
        !> return the permutation of an SQS
        procedure :: alloy_site_permutation
        !> permute the positions in a structure
        procedure :: permute_positions
        !> return the names of the unique atoms
        procedure :: unique_atom_label

        !> Initialize from file
        procedure :: readfromfile
        !> write to file, various formats
        procedure :: writetofile
        !> read isotope distribution from file
        procedure :: readisotopefromfile
        !> write the Brillouin zone to hdf5
        procedure :: write_bz_to_hdf5
        !> write structure information to hdf5
        procedure :: write_structure_to_hdf5

        !> measure size in memory
        procedure :: size_in_mem=>structure_size_in_mem
end type

! Global max for the max number of components in an alloy
integer, parameter :: max_n_components=6

! Interfaces to things in type_crystalstructure_atomdata
interface
    module function z_to_symbol(z_nucleus) result(symbol)
        integer, intent(in) :: z_nucleus
        character(len=2) :: symbol
    end function
    module function symbol_to_z(symbol) result(z_nucleus)
        character(len=*), intent(in) :: symbol
        integer :: z_nucleus
    end function
    module function neutron_cross_section(atomic_number) result(xs)
        integer, intent(in) :: atomic_number
        real(flyt) :: xs
    end function
end interface
! Interfaces to things in type_crystalstructure_io
interface
    module subroutine readfromfile(p,filename,verbosity)
        class(lo_crystalstructure), intent(out) :: p
        character(len=*), intent(in) :: filename
        integer, intent(in), optional :: verbosity
    end subroutine
    module subroutine writetofile(p,filename,output_format,write_velocities,transformationmatrix)
        class(lo_crystalstructure), intent(in) :: p
        character(len=*), intent(in) :: filename
        integer, intent(in) :: output_format
        logical, intent(in), optional :: write_velocities
        real(flyt), dimension(3,3), intent(out), optional :: transformationmatrix
    end subroutine
    module subroutine readisotopefromfile(p)
        class(lo_crystalstructure), intent(inout) :: p
    end subroutine
    module subroutine write_bz_to_hdf5(p,filename,input_id)
        class(lo_crystalstructure), intent(in) :: p
        character(len=*), intent(in) :: filename
        integer(HID_T), intent(in), optional :: input_id
    end subroutine
    module subroutine write_structure_to_hdf5(p,filename,input_id)
        class(lo_crystalstructure), intent(in) :: p
        character(len=*), intent(in), optional :: filename
        integer(HID_T), intent(in), optional :: input_id
    end subroutine
end interface
! Interfaces to things in type_crystalstructure_symmetry
interface
    module function coordinate_from_high_symmetry_point_label(p,label,previous) result(qpoint)
        class(lo_crystalstructure), intent(in) :: p
        character(len=*), intent(in) :: label
        real(flyt), dimension(3), intent(in), optional :: previous
        real(flyt), dimension(3) :: qpoint
    end function
    module subroutine classify(p,how,uc,tolerance,refine,timereversal)
        class(lo_crystalstructure), intent(inout) :: p
        character(len=*), intent(in) :: how
        type(lo_crystalstructure), intent(in), optional :: uc
        real(flyt), intent(in), optional :: tolerance
        logical, intent(in), optional :: refine
        logical, intent(in), optional :: timereversal
    end subroutine
end interface
! Interfaces to utility alloy routines
interface
    module subroutine alloy_site_permutation(ss,sqs,permutation)
        class(lo_crystalstructure), intent(in) :: ss
        type(lo_crystalstructure), intent(in) :: sqs
        integer, dimension(:), intent(out) :: permutation
    end subroutine
    module subroutine permute_positions(p,permutation,forward)
        class(lo_crystalstructure), intent(inout) :: p
        integer, dimension(:), intent(in) :: permutation
        logical, intent(in) :: forward
    end subroutine
end interface

contains

!> Transform lists of stuff to a proper object. How a crystal structure should be initialized.
subroutine generate(p,latticevectors,positions,atomic_numbers,enhet,velocities,header,&
                    collmag,noncollmag,cmatom,collmagmom,noncollmagmom,verbosity,alloy,&
                    alloy_componentcounter,alloy_components,alloy_concentrations)
    !> crystal structure
    class(lo_crystalstructure), intent(out) :: p
    !> basis
    real(flyt), dimension(3,3), intent(in) :: latticevectors
    !> positions of atoms, in fractional coordinates
    real(flyt), dimension(:,:), intent(in) :: positions
    !> the atomic numbers of each atom
    integer, dimension(:), intent(in) :: atomic_numbers
    !> in what unit is the lattice vectors?
    integer, intent(in) :: enhet
    !> specify velocities?
    real(flyt), dimension(:,:), intent(in), optional :: velocities
    !> give it a name?
    character(len=*), intent(in), optional :: header
    !> collinear magnetism?
    logical, intent(in), optional :: collmag
    !> non-collinear magnetism?
    logical, intent(in), optional :: noncollmag
    !> if so, which atoms are to be randomized?
    logical, dimension(:), intent(in), optional :: cmatom
    !> some magnetic moments
    real(flyt), dimension(:), intent(in), optional :: collmagmom
    !> and some non-collinear magnetic moments
    real(flyt), dimension(:,:), intent(in), optional :: noncollmagmom
    !> is it an alloy
    logical, intent(in), optional :: alloy
    !> number of alloy components per site
    integer, dimension(:), intent(in), optional :: alloy_componentcounter
    !> what species are there per site
    integer, dimension(:,:), intent(in), optional :: alloy_components
    !> concentration of these species
    real(flyt), dimension(:,:), intent(in), optional :: alloy_concentrations
    !> talk a lot?
    integer, intent(in), optional :: verbosity

    ! Since it is intent(out), nothing is allocated. Set everything to nothing, to be on the
    ! safe side, or something if it is obvious. Also, have some heuristics to figure out if the
    ! input is stupid.
    if ( size(positions,2) .eq. size(atomic_numbers,1) ) then
        p%na=size(positions,2)
    else
        call lo_stop_gracefully(['You need the same number of positions as atomic numbers when creating a structure.'],&
                                lo_exitcode_param,__FILE__,__LINE__)
    endif

    ! Right now, I have done nothing, and know nothing.
    if ( present(verbosity) ) then
        p%info%verbosity=verbosity
    else
        p%info%verbosity=0
    endif
    if ( present(header) ) then
        p%info%title=trim(adjustl(header))
    else
        p%info%title='cell' ! This has to be something.
    endif
    if ( present(alloy) ) then
        p%info%alloy=alloy
    else
        p%info%alloy=.false.
    endif
    if ( present(collmag) ) then
        p%info%collmag=collmag
    else
        p%info%collmag=.false.
    endif
    if ( present(noncollmag) ) then
        p%info%noncollmag=noncollmag
    else
        p%info%noncollmag=.false.
    endif
    p%info%havespacegroup=.false.
    p%info%havebz=.false.
    p%info%havewedge=.false.
    p%info%havebravais=.false.
    p%info%supercell=.false.
    p%info%pointslabelled=.false.
    p%info%decidedtimereversal=.false.
    p%info%supercellmatrix=-1
    p%info%unique_primitive_basis=-1
    p%info%unique_conventional_basis=-1
    p%info%permutation_to_unique=-1
    p%info%transformation_to_unique=-1
    p%info%bravaislattice='dunno'
    p%info%bravaislatticelongname='dunno'
    p%info%unitcell_lattice_parameter=-1.0_flyt

    ! Get the lattice vectors every which way, with some safety checks, as well as setting the positions
    ! the volume and stuff like that.
    fixlattice: block
        real(flyt) :: posunitfactor,velunitfactor
        integer :: i
        ! Make sure we have nonzero volume
        if ( abs(lo_determ(latticevectors))/lo_frobnorm(latticevectors) .lt. lo_tol ) then
            call lo_stop_gracefully(['The determinant of the basis is really small: '//tochar(lo_determ(latticevectors))],&
                                    lo_exitcode_param,__FILE__,__LINE__)
        endif
        ! What unit are they in?
        select case(enhet)
        case(1)
            ! This is vasp-ish units, positions in A and velocitied in A/fs
            ! convert that to Hartree atomic units
            posunitfactor=lo_A_to_bohr
            velunitfactor=lo_velocity_Afs_to_au
        case(2)
            ! Input is in atomic units. Distances in Bohr, velocities in something weird.
            posunitfactor=1.0_flyt
            velunitfactor=1.0_flyt
        case default
            call lo_stop_gracefully(['Unknown unit for positions and velocities'],&
                                    lo_exitcode_param,__FILE__,__LINE__)
        end select
        ! Store them
        p%latticevectors=latticevectors*posunitfactor
        p%inv_latticevectors=lo_invert3x3matrix( p%latticevectors )
        p%reciprocal_latticevectors=lo_reciprocal_basis(p%latticevectors)
        p%inv_reciprocal_latticevectors=lo_invert3x3matrix( p%reciprocal_latticevectors )
        ! Clean a little, for cosmetic reasons
        p%latticevectors=lo_chop(p%latticevectors,lo_sqtol)
        p%inv_latticevectors=lo_chop(p%inv_latticevectors,lo_sqtol)
        p%reciprocal_latticevectors=lo_chop(p%reciprocal_latticevectors,lo_sqtol)
        p%inv_reciprocal_latticevectors=lo_chop(p%inv_reciprocal_latticevectors,lo_sqtol)
        ! Get the volume
        p%volume=abs(lo_determ(p%latticevectors))
        ! store the positions. Do the cleaning thing just to be sure.
        lo_allocate(p%r(3,p%na))
        lo_allocate(p%rcart(3,p%na))
        lo_allocate(p%u(3,p%na))
        lo_allocate(p%v(3,p%na))
        lo_allocate(p%f(3,p%na))
        p%r=0.0_flyt
        p%rcart=0.0_flyt
        p%u=0.0_flyt
        p%v=0.0_flyt
        p%f=0.0_flyt
        do i=1,p%na
            p%r(:,i)=lo_chop( lo_clean_fractional_coordinates(positions(:,i)) , lo_sqtol )
            p%rcart(:,i)=lo_chop( p%fractional_to_cartesian(p%r(:,i)) , lo_sqtol )
        enddo
        ! also, velocities and stuff
        if ( present(velocities) ) then
            p%v=velocities*velunitfactor
        endif
    end block fixlattice

    if ( p%info%alloy ) then
    ! If I have an alloy this gets rather tedious.
    fixalloyatoms: block
        type(lo_isotope_distribution), dimension(:), allocatable :: dumiso
        integer :: i,j,k,l,nc

        ! Make some space for the things that need to be specified
        lo_allocate(p%species(p%na))
        lo_allocate(p%isotope(p%na))
        lo_allocate(p%mass(p%na))
        lo_allocate(p%invsqrtmass(p%na))
        lo_allocate(p%atomic_number(p%na))
        lo_allocate(p%inelastic_neutron_cross_section(p%na))

        ! First step is to figure out how many distinct species I have
        p%species=0
        p%nelements=0
        do i=1,p%na
            if ( p%species(i) .ne. 0 ) cycle
            ! new species!
            p%nelements=p%nelements+1
            p%species(i)=p%nelements
            do j=i+1,p%na
                if ( alloy_componentcounter(i) .eq. alloy_componentcounter(j) ) then
                if ( sum(abs(alloy_components(:,i)-alloy_components(:,j))) .eq. 0 ) then
                if ( sum(abs(alloy_concentrations(:,i)-alloy_concentrations(:,j))) .lt. lo_tol ) then
                    p%species(j)=p%species(i)
                endif
                endif
                endif
            enddo
        enddo
        ! Get a symbol for each species
        lo_allocate(p%atomic_symbol(p%nelements))
        do i=1,p%nelements
            p%atomic_symbol(i)='dunno'
        enddo

        ! Not sure if atomic number make sense here
        do i=1,p%na
            p%atomic_number(i)=-1
        enddo

        ! Count the number of components per species
        lo_allocate(p%element_counter(p%nelements))
        p%element_counter=0
        do i=1,p%na
            p%element_counter(p%species(i))=p%element_counter(p%species(i))+1
        enddo
        ! Get the masses. This is far too careful to be sane, but whatever.
        do i=1,p%na
            nc=alloy_componentcounter(i)
            lo_allocate(dumiso(nc))
            ! Get the natural isotope distribution for all components of the alloy
            do j=1,nc
                call dumiso(j)%naturaldistribution( trim(z_to_symbol(alloy_components(j,i))) )
            enddo
            ! Count total number of isotopes
            l=0
            do j=1,nc
                l=l+dumiso(j)%n
            enddo
            ! Make space
            p%isotope(i)%n=l
            lo_allocate(p%isotope(i)%conc( p%isotope(i)%n ))
            lo_allocate(p%isotope(i)%mass( p%isotope(i)%n ))
            ! Store concentrations and masses
            l=0
            do j=1,nc
            do k=1,dumiso(j)%n
                l=l+1
                p%isotope(i)%conc(l)=dumiso(j)%conc(k)*alloy_concentrations(j,i)
                p%isotope(i)%mass(l)=dumiso(j)%mass(k)
            enddo
            enddo
            ! ensure it adds up to 1
            p%isotope(i)%conc=p%isotope(i)%conc/sum(p%isotope(i)%conc)
            p%isotope(i)%mean_mass=sum(p%isotope(i)%conc*p%isotope(i)%mass)
            p%isotope(i)%disorderparameter=p%isotope(i)%mass_disorder_parameter()
            lo_deallocate(dumiso)
            ! Store mean masses
            p%mass(i)=p%isotope(i)%mean_mass
            p%invsqrtmass(i)=1.0_flyt/sqrt(p%isotope(i)%mean_mass)
            p%inelastic_neutron_cross_section(i)=0.0_flyt !neutron_cross_section(p%atomic_number(i))
        enddo

        ! Perhaps store all the alloy information somehow.
        lo_allocate(p%alloyspecies(p%nelements))
        lo_allocate(dumiso(1))
        do i=1,p%nelements
            ! locate an atom of this species
            k=-1
            do j=1,p%na
                if ( p%species(j) .eq. i ) then
                    k=j
                    exit
                endif
            enddo
            p%alloyspecies(i)%n=alloy_componentcounter(k)
            lo_allocate(p%alloyspecies(i)%mass( p%alloyspecies(i)%n ))
            lo_allocate(p%alloyspecies(i)%concentration( p%alloyspecies(i)%n ))
            lo_allocate(p%alloyspecies(i)%atomic_symbol( p%alloyspecies(i)%n ))
            lo_allocate(p%alloyspecies(i)%atomic_number( p%alloyspecies(i)%n ))
            do j=1,p%alloyspecies(i)%n
                p%alloyspecies(i)%atomic_number(j)=alloy_components(j,k)
                p%alloyspecies(i)%atomic_symbol(j)=trim(z_to_symbol( alloy_components(j,k) ))
                p%alloyspecies(i)%concentration(j)=alloy_concentrations(j,k)
                call dumiso(1)%naturaldistribution( trim(z_to_symbol(alloy_components(j,k))) )
                p%alloyspecies(i)%mass(j)=dumiso(1)%mean_mass
            enddo
        enddo
        lo_deallocate(dumiso)

        ! There might be non-alloy sublattices, fix the species and atomic number of those
        do i=1,p%nelements
            if ( p%alloyspecies(i)%n .ne. 1 ) cycle
            j=p%alloyspecies(i)%atomic_number(1)
            do k=1,p%na
                if ( p%species(k) .eq. i ) p%atomic_number(k)=j
            enddo
            p%atomic_symbol(i)=trim(z_to_symbol( j ))
        enddo

    end block fixalloyatoms
    else
    ! Figure out how classify the species in a neat way. Once I fix alloys, I have to do something else.
    ! This gets the symbols, species in the right way, masses and so on.
    fixatoms: block
        integer, dimension(:), allocatable :: dum
        integer :: i,j
        character(len=1000) :: trams
        ! start by getting the union of the atomic numbers
        call lo_return_unique(atomic_numbers,dum)
        ! number of elements in the system
        p%nelements=size(dum,1)
        lo_allocate(p%isotope(p%na))
        lo_allocate(p%atomic_symbol(p%nelements))
        lo_allocate(p%element_counter(p%nelements))
        lo_allocate(p%species(p%na))
        lo_allocate(p%mass(p%na))
        lo_allocate(p%invsqrtmass(p%na))
        lo_allocate(p%atomic_number(p%na))
        lo_allocate(p%inelastic_neutron_cross_section(p%na))
        ! Set the species
        p%element_counter=0
        do i=1,p%na
            do j=1,p%nelements
                if ( dum(j) .eq. atomic_numbers(i) ) then
                    p%species(i)=j
                    p%element_counter(j)=p%element_counter(j)+1
                endif
            enddo
        enddo
        ! Get the symbols
        do i=1,p%nelements
            p%atomic_symbol(i)=trim(adjustl(z_to_symbol( dum(i) )))
        enddo
        ! And the natural isotope distributions, and masses, and atomic numbers. Also neutron cross-section.
        do i=1,p%na
            call p%isotope(i)%naturaldistribution( p%atomic_symbol(p%species(i)) )
            p%mass(i)=p%isotope(i)%mean_mass
            p%invsqrtmass(i)=1.0_flyt/sqrt(p%mass(i))
            p%atomic_number(i)=atomic_numbers(i)
            p%inelastic_neutron_cross_section(i)=neutron_cross_section(p%atomic_number(i))
        enddo
        ! And maybe talk a little
        if ( p%info%verbosity .gt. 0 ) then
            write(*,*) 'Generating structure, lattice vectors:'
            write(*,"(1X,A,3(2X,F18.12))") '    a1:',p%latticevectors(:,1)
            write(*,"(1X,A,3(2X,F18.12))") '    a2:',p%latticevectors(:,2)
            write(*,"(1X,A,3(2X,F18.12))") '    a3:',p%latticevectors(:,3)
            write(*,*) '... reciprocal lattice vectors:'
            write(*,"(1X,A,3(2X,F18.12))") '    b1:',p%reciprocal_latticevectors(:,1)
            write(*,"(1X,A,3(2X,F18.12))") '    b2:',p%reciprocal_latticevectors(:,2)
            write(*,"(1X,A,3(2X,F18.12))") '    b3:',p%reciprocal_latticevectors(:,3)
            trams=''
            do i=1,p%nelements
                trams=trim(trams)//' '//trim(p%atomic_symbol(i))//': '//tochar(p%element_counter(i))
            enddo
            write(*,*) '... with composition '//trim(trams)
        endif
    end block fixatoms
    ! And for now, I am satisfied. I don't classify unless I really have to, I suppose.
    ! maybe add an option to agressively classify later on.
    endif

    ! If I specify magnetic moments, that is useful to take care of.
    if ( p%info%collmag .or. p%info%noncollmag ) then
    fixmag: block
        real(flyt), dimension(:,:), allocatable :: unmom,dum
        integer, dimension(:), allocatable :: dumi1,dumi2
        integer :: i,j,ii
        logical :: coll
        ! Quick check of one or three magnetic moments
        if ( p%info%collmag ) then
            coll=.true.
        else
            coll=.false.
        endif
        ! Small sanity check
        if ( p%info%collmag .and. p%info%noncollmag ) then
            call lo_stop_gracefully(['Please specifify either collinear or noncollinear magnetism, not both.'],&
                                    lo_exitcode_param,__FILE__,__LINE__)
        endif
        ! Count the number of magnetic atoms, and get their moments
        lo_allocate(dum(3,p%na))
        dum=0.0_flyt
        ii=0
        do i=1,p%na
            if ( cmatom(i) ) then
                ii=ii+1
                if ( coll ) then
                    dum(1,ii)=collmagmom(i)
                else
                    dum(:,ii)=noncollmagmom(:,i)
                endif
            endif
        enddo
        ! Get the unique moments
        call lo_return_unique(dum(:,1:ii),unmom,lo_tol)
        ! Build a new kind of species
        lo_allocate(dumi1(p%na))
        do i=1,p%na
            if ( cmatom(i) ) then
                do j=1,size(unmom,2)
                    if ( lo_sqnorm(unmom(:,j)-dum(:,i)) .lt. lo_sqtol ) then
                        dumi1(i)=atomic_numbers(i)*1000+j
                    endif
                enddo
            else
                dumi1(i)=atomic_numbers(i)
            endif
        enddo
        ! The new, unique species
        call lo_return_unique(dumi1,dumi2)
        ! Start storing some things
        p%mag%nmagspecies=size(dumi2,1)
        lo_allocate(p%mag%magspeciescounter( p%mag%nmagspecies ))
        lo_allocate(p%mag%magspecies( p%na ))
        p%mag%magspeciescounter=0
        do i=1,p%na
            do j=1,p%mag%nmagspecies
                if ( dumi1(i) .eq. dumi2(j) ) then
                    p%mag%magspeciescounter(j)=p%mag%magspeciescounter(j)+1
                    p%mag%magspecies(i)=j
                    exit
                endif
            enddo
        enddo
        lo_allocate(p%mag%atom_has_moment(p%na))
        if ( coll ) then
            lo_allocate(p%mag%collinear_moment(p%na))
            p%mag%collinear_moment=collmagmom
        else
            lo_allocate(p%mag%noncollinear_moment(3,p%na))
            p%mag%noncollinear_moment=noncollmagmom
        endif
        p%mag%atom_has_moment=cmatom
        !
        lo_deallocate(unmom)
        lo_deallocate(dum)
        lo_deallocate(dumi1)
        lo_deallocate(dumi2)
    end block fixmag
    endif
end subroutine

!> Build a supercell. The cell will be returned in ss.
subroutine build_supercell(p,ss,dimensions,nondiagdimensions)
    !> unitcell
    class(lo_crystalstructure), intent(inout) :: p
    !> supercell
    type(lo_crystalstructure), intent(out) :: ss
    !> how many times to repeat the unit cell along a1,a2,a3, the basis vectors.
    integer, dimension(3), intent(in), optional :: dimensions
    !> non-diagonal supercell. Fancy.
    integer, dimension(3,3), intent(in), optional :: nondiagdimensions

    real(flyt), dimension(:,:), allocatable :: r,noncollmagmom,alloy_concentrations
    real(flyt), dimension(:), allocatable :: collmagmom
    real(flyt), dimension(3,3) :: basis
    integer, dimension(:,:), allocatable :: alloy_components
    integer, dimension(:), allocatable :: atomic_numbers,alloy_componentcounter
    logical, dimension(:), allocatable :: magatom
    logical :: diagonal
    character(len=1000) :: hdr

    ! Basic heuristics
    if ( present(dimensions) ) then
        diagonal=.true.
    elseif ( present(nondiagdimensions) ) then
        diagonal=.false.
    else
        call lo_stop_gracefully(['You need to specify dimensions when building a supercell somehow.'],lo_exitcode_param,__FILE__,__LINE__)
    endif

    ! I want to know what Bravais lattice it is
    if ( p%info%havebravais .eqv. .false. ) call p%classify('bravais')

    if ( diagonal ) then
        call diagonal_cell(p,dimensions,basis,atomic_numbers,r,magatom,collmagmom,noncollmagmom,&
                          alloy_componentcounter,alloy_components,alloy_concentrations)
        hdr=tochar(dimensions(1))//'x'//tochar(dimensions(2))//'x'//tochar(dimensions(3))//' supercell'
    else
        call nondiagonal_cell(p,nondiagdimensions,basis,atomic_numbers,r,magatom,collmagmom,noncollmagmom,&
                              alloy_componentcounter,alloy_components,alloy_concentrations)
        hdr='non-diagonal supercell'
    endif

    ! Create a nice object from this
    if ( p%info%collmag ) then
        call ss%generate(basis,r,atomic_numbers,enhet=2,header=trim(hdr),verbosity=p%info%verbosity,collmag=.true.,&
                         cmatom=magatom,collmagmom=collmagmom)
    elseif ( p%info%noncollmag ) then
        call ss%generate(basis,r,atomic_numbers,enhet=2,header=trim(hdr),verbosity=p%info%verbosity,&
                         noncollmag=.true.,cmatom=magatom,noncollmagmom=noncollmagmom)
    elseif ( p%info%alloy ) then
        call ss%generate(basis,r,atomic_numbers,enhet=2,header=trim(hdr),verbosity=p%info%verbosity,&
                        alloy=.true.,alloy_componentcounter=alloy_componentcounter,&
                        alloy_components=alloy_components,alloy_concentrations=alloy_concentrations)
    else
        call ss%generate(basis,r,atomic_numbers,enhet=2,header=trim(hdr),verbosity=p%info%verbosity)
    endif

    ! Store the latticeparameter, perhaps
    ss%info%unitcell_lattice_parameter=p%info%unitcell_lattice_parameter
    contains

    !> build normal diagonal supercell
    subroutine diagonal_cell(p,dimensions,basis,atomic_numbers,r,magatom,collmagmom,noncollmagmom,alloy_componentcounter,alloy_components,alloy_concentrations)
        !> unit cell
        type(lo_crystalstructure), intent(in) :: p
        !> supercell dimensions
        integer, dimension(3), intent(in) :: dimensions
        !> supercell positions
        real(flyt), dimension(:,:), allocatable, intent(out) :: r
        !> supercell species
        integer, dimension(:), allocatable, intent(out) :: atomic_numbers
        !> supercell lattice vectors
        real(flyt), dimension(3,3), intent(out) :: basis
        !> which atoms are magnetic
        logical, dimension(:), allocatable, intent(out) :: magatom
        !> what their collinear magnetic moments
        real(flyt), dimension(:), allocatable, intent(out) :: collmagmom
        !> or the non-collinear magnetic moments
        real(flyt), dimension(:,:), allocatable, intent(out) :: noncollmagmom
        !> alloy component counter
        integer, dimension(:), allocatable, intent(out) :: alloy_componentcounter
        !> what are the components
        integer, dimension(:,:), allocatable, intent(out) :: alloy_components
        !> concentration of said components
        real(flyt), dimension(:,:), allocatable, intent(out) :: alloy_concentrations

        integer :: i,j,k,l,a1,na,ii,jj

        ! build the cell
        basis=p%latticevectors
        basis(:,1)=basis(:,1)*dimensions(1)
        basis(:,2)=basis(:,2)*dimensions(2)
        basis(:,3)=basis(:,3)*dimensions(3)
        na=product(dimensions)*p%na
        lo_allocate(r(3,na))
        lo_allocate(atomic_numbers(na))
        lo_allocate(magatom(na))
        lo_allocate(collmagmom(na))
        lo_allocate(noncollmagmom(3,na))
        lo_allocate(alloy_componentcounter(na))
        lo_allocate(alloy_components(max_n_components,na))
        lo_allocate(alloy_concentrations(max_n_components,na))
        r=0.0_flyt
        atomic_numbers=0
        magatom=.false.
        collmagmom=0.0_flyt
        noncollmagmom=0.0_flyt
        alloy_componentcounter=0
        alloy_components=0
        alloy_concentrations=0.0_flyt
        ! loop and build everything
        l=0
        do a1=1,p%na
            do i=1,dimensions(1)
            do j=1,dimensions(2)
            do k=1,dimensions(3)
                l=l+1
                r(:,l)=[i-1,j-1,k-1]*1.0_flyt+p%r(:,a1)
                do ii=1,3
                    r(ii,l)=r(ii,l)/(1.0_flyt*dimensions(ii))
                enddo
                atomic_numbers(l)=p%atomic_number(a1)
                if ( p%info%collmag ) then
                    magatom(l)=p%mag%atom_has_moment(a1)
                    collmagmom(l)=p%mag%collinear_moment(a1)
                elseif ( p%info%noncollmag ) then
                    magatom(l)=p%mag%atom_has_moment(a1)
                    noncollmagmom(:,l)=p%mag%noncollinear_moment(:,a1)
                elseif ( p%info%alloy ) then
                    ii=p%species(a1)
                    jj=p%alloyspecies(ii)%n
                    alloy_componentcounter(l)=p%alloyspecies(ii)%n
                    alloy_concentrations(1:jj,l)=p%alloyspecies(ii)%concentration(1:jj)
                    alloy_components(1:jj,l)=p%alloyspecies(ii)%atomic_number(1:jj)
                else
                    magatom(l)=.false.
                    collmagmom(l)=0.0_flyt
                    noncollmagmom(:,l)=0.0_flyt
                endif
            enddo
            enddo
            enddo
        enddo
        ! Clean the positions a little:
        r=lo_chop(lo_clean_fractional_coordinates(r),lo_sqtol)
        basis=lo_chop(basis,lo_sqtol)
    end subroutine

    !> build fancy non-diagonal supercell
    subroutine nondiagonal_cell(p,dimensions,basis,atomic_numbers,r,magatom,collmagmom,noncollmagmom,alloy_componentcounter,alloy_components,alloy_concentrations)
        !> unit cell
        type(lo_crystalstructure), intent(in) :: p
        !> supercell dimensions
        integer, dimension(3,3), intent(in) :: dimensions
        !> supercell lattice vectors
        real(flyt), dimension(3,3), intent(out) :: basis
        !> supercell species
        integer, dimension(:), allocatable, intent(out) :: atomic_numbers
        !> supercell positions
        real(flyt), dimension(:,:), allocatable, intent(out) :: r
        !> which atoms are magnetic
        logical, dimension(:), allocatable, intent(out) :: magatom
        !> what are their magnetic moments
        real(flyt), dimension(:), allocatable, intent(out) :: collmagmom
        !> or the non-collinear magnetic moments
        real(flyt), dimension(:,:), allocatable, intent(out) :: noncollmagmom
        !> alloy component counter
        integer, dimension(:), allocatable, intent(out) :: alloy_componentcounter
        !> what are the components
        integer, dimension(:,:), allocatable, intent(out) :: alloy_components
        !> concentration of said components
        real(flyt), dimension(:,:), allocatable, intent(out) :: alloy_concentrations

        real(flyt), dimension(:,:), allocatable :: dumr1,dumr2
        real(flyt), dimension(3,3) :: invbasis,m0,m1
        real(flyt), dimension(3) :: v0
        real(flyt) :: f0
        integer, dimension(:), allocatable :: ind
        integer :: i,j,k,l,a1,na,nrep,ctr,ii,jj

        ! Get the new basis
        i=lo_determ(dimensions)
        if ( i .eq. 0 ) then
            call lo_stop_gracefully(['Determinant of supercell matrix is zero'],lo_exitcode_param,__FILE__,__LINE__)
        elseif ( i .lt. 0 ) then
            write(*,*) 'NOTE: I flipped sign of the supercell matrix to not invert the supercell'
            basis=-lo_chop(matmul(p%latticevectors,dimensions),lo_sqtol)
        else
            basis=lo_chop(matmul(p%latticevectors,dimensions),lo_sqtol)
        endif
        ! Possibly, this could be cleaned in a smart way, to get good precision for
        ! the more obvious cases, such as cubic/tetragonal/hexagonal
        invbasis=lo_invert3x3matrix( basis )

        ! Figure out how many times the unitcell should be repeated
        f0=lo_bounding_sphere_of_box(basis)*2
        nrep=1
        do
            m0=p%latticevectors*(2*nrep+1)
            if ( lo_inscribed_sphere_in_box(m0) .gt. f0 ) then
                ! got it, was enough
                exit
            else
                ! check larger
                nrep=nrep+1
            endif
        enddo
        nrep=nrep+1 ! just to be on the safe side.
        ! Expected number of atoms
        na=abs(int(anint(lo_determ(dimensions*1.0_flyt)))*p%na)
        m1=lo_chop( matmul(invbasis,p%latticevectors) ,lo_sqtol)

        ! Count atoms
        ctr=0
        do a1=1,p%na
            do i=-nrep,nrep
            do j=-nrep,nrep
            do k=-nrep,nrep
                v0=matmul(m1,[i,j,k]+p%r(:,a1)) ! fractional coordinates in new system
                if ( minval(v0) .gt. -lo_sqtol .and. maxval(v0) .lt. 1.0_flyt+lo_sqtol ) ctr=ctr+1
            enddo
            enddo
            enddo
        enddo

        ! Store these
        lo_allocate(dumr1(6,ctr))
        ctr=0
        do a1=1,p%na
            do i=-nrep,nrep
            do j=-nrep,nrep
            do k=-nrep,nrep
                v0=matmul(m1,[i,j,k]+p%r(:,a1)) ! fractional coordinates in new system
                if ( minval(v0) .gt. -lo_sqtol .and. maxval(v0) .lt. 1.0_flyt+lo_sqtol ) then
                    ctr=ctr+1
                    dumr1(1:3,ctr)=lo_chop(lo_clean_fractional_coordinates(v0),lo_sqtol)
                    dumr1(4,ctr)=p%atomic_number(a1)
                    dumr1(5,ctr)=p%species(a1)
                    if ( p%info%collmag .or. p%info%noncollmag ) then
                        dumr1(6,ctr)=p%mag%magspecies(a1)
                    else
                        dumr1(6,ctr)=0
                    endif
                endif
            enddo
            enddo
            enddo
        enddo

        ! get the unique
        call lo_return_unique(dumr1,dumr2)

        ! Simple sanity check
        if ( size(dumr2,2) .ne. na ) then
            call lo_stop_gracefully(['Failed building supercell. Probably a pathological supercell matrix'],lo_exitcode_symmetry,__FILE__,__LINE__)
        endif

        ! A bit more serious sanity check: this is how many of each
        ! atom I should find.
        l=abs(lo_determ(dimensions))
        do i=1,p%nelements
            k=0
            do j=1,size(dumr2,2)
                if ( nint(dumr2(5,j)) .eq. i ) k=k+1
            enddo
            if ( k .ne. l*p%element_counter(i) ) then
                call lo_stop_gracefully(['Got the correct number of atoms, but the composition is off. Not good.'],lo_exitcode_symmetry,__FILE__,__LINE__)
            endif
        enddo

        ! Now I suppose things are ok. Prepare thing for output:
        basis=lo_choplarge(basis/p%info%unitcell_lattice_parameter,lo_sqtol)*p%info%unitcell_lattice_parameter
        lo_allocate(ind(na))
        lo_allocate(r(3,na))
        lo_allocate(atomic_numbers(na))
        lo_allocate(magatom(na))
        lo_allocate(collmagmom(na))
        lo_allocate(noncollmagmom(3,na))
        lo_allocate(alloy_componentcounter(na))
        lo_allocate(alloy_components(max_n_components,na))
        lo_allocate(alloy_concentrations(max_n_components,na))
        ind=0
        r=0.0_flyt
        atomic_numbers=0
        magatom=.false.
        collmagmom=0.0_flyt
        noncollmagmom=0.0_flyt
        alloy_componentcounter=0
        alloy_components=0
        alloy_concentrations=0.0_flyt
        call qsort(dumr2(5,:),ind)
        do i=1,na
            j=ind(i)
            r(:,i)=dumr2(1:3,j)
            atomic_numbers(i)=nint(dumr2(4,j))
            if ( p%info%collmag ) then
                do k=1,p%na
                    if ( nint(dumr2(6,j)) .eq. p%mag%magspecies(k) ) then
                        magatom(i)=p%mag%atom_has_moment(k)
                        collmagmom(i)=p%mag%collinear_moment(k)
                    endif
                enddo
            elseif ( p%info%noncollmag ) then
                do k=1,p%na
                    if ( nint(dumr2(6,j)) .eq. p%mag%magspecies(k) ) then
                        magatom(i)=p%mag%atom_has_moment(k)
                        noncollmagmom(:,i)=p%mag%noncollinear_moment(:,k)
                    endif
                enddo
            elseif ( p%info%alloy ) then
                ii=nint(dumr2(5,i)) ! note i here, not j
                jj=p%alloyspecies(ii)%n
                alloy_componentcounter(i)=p%alloyspecies(ii)%n
                alloy_concentrations(1:jj,i)=p%alloyspecies(ii)%concentration(1:jj)
                alloy_components(1:jj,i)=p%alloyspecies(ii)%atomic_number(1:jj)
            endif
        enddo
        ! and things are done!
        lo_deallocate(ind)
        lo_deallocate(dumr1)
        lo_deallocate(dumr2)
    end subroutine
end subroutine

!> Get the kinetic energy of the atoms in the cell. Will return it in Hartree/cell, not per atom.
real(flyt) function kinetic_energy(p)
    !> the crystal structure
    class(lo_crystalstructure), intent(in) :: p

    integer :: i
    real(flyt) :: mass,normv,ek

    ek=0.0_flyt
    do i=1,p%na
        mass=p%mass(i)
        normv=lo_sqnorm(p%v(:,i))
        ek=ek+mass*normv*0.5_flyt
    enddo
    kinetic_energy=ek
end function

!> The shortest distance in a cell, should be the shortest neighbour distance.
#ifdef AGRESSIVE_SANITY
function mincutoff(p) result(r)
#else
pure function mincutoff(p) result(r)
#endif
    !> the crystal structure
    class(lo_crystalstructure), intent(in) :: p
    !> the cutoff
    real(flyt) :: r

    integer :: i,j
    real(flyt), dimension(3) :: v
    real(flyt) :: f0,f1

    f0=-lo_huge
    do i=1,p%na
        f1=lo_huge
        do j=1,p%na
            v=p%r(:,j)-p%r(:,i)
            v=p%displacement_fractional_to_cartesian(v)
            if ( norm2(v) .gt. lo_tol ) then
                f1=min(f1,norm2(v))
            endif
        enddo
        f0=max(f0,f1)
    enddo
    ! Catch pathological case of feeding a unitcell to this routine
    if ( f0 .lt. lo_tol .or. p%na .eq. 1 ) then
        f0=lo_bounding_sphere_of_box(p%latticevectors)
    endif
    ! And a little bit more for good measure.
    r=f0*1.001_flyt
end function

!> calculate the nearest neighbour distance. Not fast.
#ifdef AGRESSIVE_SANITY
function nearest_neighbour_distance(p) result(r)
#else
pure function nearest_neighbour_distance(p) result(r)
#endif
    !> structure
    class(lo_crystalstructure), intent(in) :: p
    !> nearest neighbour distance
    real(flyt) :: r

    real(flyt), dimension(3,3) :: m0
    real(flyt), dimension(3) :: v0
    real(flyt) :: f0
    integer :: nrep,a1,a2,i,j,k

    f0=lo_bounding_sphere_of_box(p%latticevectors)
    do nrep=1,100
        m0=p%latticevectors*(2*nrep+1)
        if ( lo_inscribed_sphere_in_box(m0) .gt. f0 ) exit
    enddo

    r=lo_huge
    do a1=1,p%na
    do a2=1,p%na
        do i=-nrep,nrep
        do j=-nrep,nrep
        do k=-nrep,nrep
            v0=[i,j,k]+p%r(:,a2)-p%r(:,a1)
            f0=norm2(matmul(p%latticevectors,v0))
            if ( f0 .gt. lo_tol ) r=min(r,f0)
        enddo
        enddo
        enddo
    enddo
    enddo
end function

!> The longest distance in a cell. Will return the radius of the largest sphere that can fit safely inside the cell
#ifdef AGRESSIVE_SANITY
function maxcutoff(p) result(r)
#else
pure function maxcutoff(p) result(r)
#endif
    !> the crystal structure
    class(lo_crystalstructure), intent(in) :: p
    !> the cutoff
    real(flyt) :: r
    ! just get the radius of the inscribed sphere, with some tolerance
    r=lo_inscribed_sphere_in_box(p%latticevectors)-10*lo_tol
end function

!> Converts a vector from cartesian to fractional coordinates. Not fast at all, don't use for something speed-sensitive.
pure function cartesian_to_fractional(p,v,reciprocal,pbc) result(r)
    !> the crystal structure
    class(lo_crystalstructure), intent(in) :: p
    !> the vector in cartesian coordinates
    real(flyt), dimension(3), intent(in) :: v
    !> in reciprocal space?
    logical, intent(in), optional :: reciprocal
    !> should the pbc checks be done
    logical, intent(in), optional :: pbc
    !> the vector in fractional coordinates
    real(flyt), dimension(3) :: r
    !
    logical :: check,reclat

    ! Choice wether to care about periodic boundary conditions
    if ( present(pbc) ) then
        check=pbc
    else
        check=.true.
    endif
    ! Real or reciprocal space?
    if ( present(reciprocal) ) then
        reclat=reciprocal
    else
        reclat=.false.
    endif

    if ( reclat ) then
        r=matmul(p%inv_reciprocal_latticevectors,v)
    else
        r=matmul(p%inv_latticevectors,v)
    endif
    if ( check ) r=lo_clean_fractional_coordinates(r)
end function

!> Convert a displacement vector from fractional to cartesian coordinates. This means that a displacement in fractional coordinates larger than 0.5 will be shifted, as will those smaller than -0.5
pure function displacement_fractional_to_cartesian(p,v,reciprocal) result(r)
    !> the crystal structure
    class(lo_crystalstructure), intent(in) :: p
    !> displacement vector to be converted
    real(flyt), dimension(3), intent(in) :: v
    !> in reciprocal space?
    logical, intent(in), optional :: reciprocal
    !> cartesian displacement vector
    real(flyt), dimension(3) :: r

    ! local
    real(flyt), dimension(3) :: w
    logical :: reclat

    if ( present(reciprocal) ) then
        reclat=reciprocal
    else
        reclat=.false.
    endif

    w=lo_clean_fractional_coordinates(v+0.5_flyt)-0.5_flyt
    if ( reclat ) then
        r=matmul(p%reciprocal_latticevectors,w)
    else
        r=matmul(p%latticevectors,w)
    endif
end function

!> Converts a vector from fractional to cartesian coordinates.
pure function fractional_to_cartesian(p,v,reciprocal) result(r)
    !> the crystal structure
    class(lo_crystalstructure), intent(in) :: p
    !> vector in fractional coordinates
    real(flyt), dimension(3), intent(in) :: v
    !> in reciprocal space?
    logical, intent(in), optional :: reciprocal
    !> vector in cartesian coordinates
    real(flyt), dimension(3) :: r

    logical :: reclat
    if ( present(reciprocal) ) then
        reclat=reciprocal
    else
        reclat=.false.
    endif
    if ( reclat ) then
        r=matmul(p%reciprocal_latticevectors,v)
    else
        r=matmul(p%latticevectors,v)
    endif
end function

!> Gives the reciprocal lattice vector that moves q into the first BZ. Uses the information from the Brilloin zone polyhedron and moves around through edges in a pretty neat way. It is pretty fast.
pure function gshift(bz,point) result(g)
    !> brillouin zone
    class(lo_brillouin_zone), intent(in) :: bz
    !> point
    real(flyt), dimension(3), intent(in) :: point
    !> shift
    real(flyt), dimension(3) :: g

    integer :: i,facectr
    real(flyt), dimension(3) :: sh,pt,nv

    ! First a very cheap check, if the point is inside the inscribed sphere of the BZ, no need
    ! to do expensive checks, the shift is zero.
    if ( lo_sqnorm(point) .lt. bz%rmin**2 ) then
        g=0.0_flyt
        return
    endif

    pt=point
    sh=0.0_flyt
    shiftloop: do
        facectr=0
        faceloop: do i=1,bz%nfaces
            ! Check the signed distance between the point and all the planes constituting
            ! the Brillouin zone.
            if ( bz%face(i)%plane%distance_to_point(pt) .gt. lo_tol ) then
                ! it means the point is outside this face, get the vector pointing to the
                ! Gamma-point on the other side of this face
                nv=bz%face(i)%neighbourvector
                ! add this G-vector to the shift
                sh=sh+nv
                ! move the point with this vector
                pt=pt-nv
                exit faceloop
            else
                ! the point is on the correct side of this face
                facectr=facectr+1
            endif
        enddo faceloop
        ! We might be back in the first BZ if we are on the negative side of all faces.
        if ( facectr .eq. bz%nfaces ) exit shiftloop
    enddo shiftloop
    ! And we are done!
    g=sh
end function

!> Calculates distance to the Brillouin zone edge from a given wave vector. Will return a value between 0 and 1, where 1 is the zone center and 0 is exactly on the edge.
pure function distance_to_zone_edge(bz,q) result(d)
    !> the Brillouin zone
    class(lo_brillouin_zone), intent(in) :: bz
    !> the q-point
    real(flyt), dimension(3), intent(in) :: q
    !> the distance
    real(flyt) :: d
    !
    real(flyt), dimension(3) :: v0
    real(flyt) :: f0,f1,f2,qnorm
    integer :: i

    ! first check if we are at gamma, then we don't have to care at all.
    qnorm=norm2(q)
    if ( qnorm .lt. lo_tol ) then
        ! this is gamma
        d=1.0_flyt
        return
    endif

    ! unit vector for the line from Gamma to q.
    v0=q/qnorm
    ! Check the intersection between this line and all possible planes
    f0=lo_huge*0.5_flyt
    d=lo_huge
    do i=1,bz%nfaces
        f1=dot_product(v0,bz%face(i)%plane%normal)
        if ( abs(f1) .gt. lo_tol ) then
            ! not perpendicular, intersection at finite distance
            f2=bz%face(i)%plane%p/f1
        else
            ! intersection infinitely far away
            f2=lo_huge
        endif
        !
        if ( f2 .lt. f0 .and. f2 .gt. 0.0_flyt ) then
            f0=f2
            d=min(d,f2)
        endif
    enddo
    ! Straight ratio, 0 at zone center, 1 at edge
    d=lo_chop(1.0_flyt-qnorm/d,lo_sqtol)
end function

!> return a string containing the names of unique atoms separated by spaces
subroutine unique_atom_label(p,labelstring)
    !> structure
    class(lo_crystalstructure), intent(in) :: p
    !> labels
    character(len=2000), intent(out) :: labelstring

    character(len=100), dimension(:), allocatable :: unique_atom_names
    integer :: a1,a2,s1,s2,ctr

    ! Get names for the unique atoms?
    allocate(unique_atom_names(p%sym%n_irreducible_atom))
    unique_atom_names='dunno'
    a1l: do a1=1,p%sym%n_irreducible_atom
        ! If already decided, don't bother.
        if ( trim(unique_atom_names(a1)) .ne. 'dunno' ) cycle
        s1=p%species(p%sym%irr_to_all(a1))
        ctr=0
        do a2=1,p%sym%n_irreducible_atom
            s2=p%species(p%sym%irr_to_all(a2))
            if ( s1 .eq. s2 ) ctr=ctr+1
        enddo

        if ( ctr .eq. 1 ) then
            ! If only one, just give it the normal name
            unique_atom_names(a1)=trim(adjustl(p%atomic_symbol(s1)))
            ! and move on to the next atom
            cycle a1l
        endif

        ! More than one unique atom of this kind, decorate the names
        ctr=0
        do a2=1,p%sym%n_irreducible_atom
            s2=p%species(p%sym%irr_to_all(a2))
            if ( s1 .eq. s2 ) then
                ctr=ctr+1
                unique_atom_names(a2)=trim(adjustl(p%atomic_symbol(s2)))//"_"//tochar(ctr)
            endif
        enddo
    enddo a1l

    ! Stuff the names into a string
    labelstring=""
    do a1=1,p%sym%n_irreducible_atom
        labelstring=trim(adjustl(labelstring))//" "//trim(adjustl(unique_atom_names(a1)))
    enddo
    deallocate(unique_atom_names)
end subroutine

!> measure size in memory, in bytes
function structure_size_in_mem(p) result(mem)
    !> dispersions
    class(lo_crystalstructure), intent(in) :: p
    !> memory in bytes
    integer(i8) :: mem

    integer :: i

    mem=0
    mem=mem+storage_size(p)

    if ( allocated(p%flavor                         ) ) mem=mem+storage_size(p%flavor                         )*size(p%flavor                         )
    if ( allocated(p%inelastic_neutron_cross_section) ) mem=mem+storage_size(p%inelastic_neutron_cross_section)*size(p%inelastic_neutron_cross_section)
    if ( allocated(p%invsqrtmass                    ) ) mem=mem+storage_size(p%invsqrtmass                    )*size(p%invsqrtmass                    )
    if ( allocated(p%mass                           ) ) mem=mem+storage_size(p%mass                           )*size(p%mass                           )
    if ( allocated(p%u                              ) ) mem=mem+storage_size(p%u                              )*size(p%u                              )
    if ( allocated(p%f                              ) ) mem=mem+storage_size(p%f                              )*size(p%f                              )
    if ( allocated(p%v                              ) ) mem=mem+storage_size(p%v                              )*size(p%v                              )
    if ( allocated(p%rcart                          ) ) mem=mem+storage_size(p%rcart                          )*size(p%rcart                          )
    if ( allocated(p%r                              ) ) mem=mem+storage_size(p%r                              )*size(p%r                              )
    if ( allocated(p%species                        ) ) mem=mem+storage_size(p%species                        )*size(p%species                        )
    if ( allocated(p%atomic_number                  ) ) mem=mem+storage_size(p%atomic_number                  )*size(p%atomic_number                  )
    if ( allocated(p%atomic_symbol                  ) ) mem=mem+storage_size(p%atomic_symbol                  )*size(p%atomic_symbol                  )
    if ( allocated(p%element_counter                ) ) mem=mem+storage_size(p%element_counter                )*size(p%element_counter                )
    if ( allocated(p%info%cellindex                 ) ) mem=mem+storage_size(p%info%cellindex                 )*size(p%info%cellindex                 )
    if ( allocated(p%info%index_in_unitcell         ) ) mem=mem+storage_size(p%info%index_in_unitcell         )*size(p%info%index_in_unitcell         )
    if ( allocated(p%mag%atom_has_moment            ) ) mem=mem+storage_size(p%mag%atom_has_moment            )*size(p%mag%atom_has_moment            )
    if ( allocated(p%mag%collinear_moment           ) ) mem=mem+storage_size(p%mag%collinear_moment           )*size(p%mag%collinear_moment           )
    if ( allocated(p%mag%noncollinear_moment        ) ) mem=mem+storage_size(p%mag%noncollinear_moment        )*size(p%mag%noncollinear_moment        )
    if ( allocated(p%mag%magspecies                 ) ) mem=mem+storage_size(p%mag%magspecies                 )*size(p%mag%magspecies                 )
    if ( allocated(p%mag%magspeciescounter          ) ) mem=mem+storage_size(p%mag%magspeciescounter          )*size(p%mag%magspeciescounter          )

    if ( allocated(p%isotope) ) then
        do i=1,size(p%isotope)
            mem=mem+storage_size(p%isotope(i))
            if ( allocated(p%isotope(i)%conc ) ) mem=mem+storage_size(p%isotope(i)%conc)*size(p%isotope(i)%conc)
            if ( allocated(p%isotope(i)%mass ) ) mem=mem+storage_size(p%isotope(i)%mass)*size(p%isotope(i)%mass)
        enddo
    endif
    if ( allocated(p%alloyspecies) ) then
        do i=1,size(p%alloyspecies)
            mem=mem+storage_size(p%alloyspecies(i))
            if ( allocated(p%alloyspecies(i)%atomic_symbol) ) mem=mem+storage_size(p%alloyspecies(i)%atomic_symbol)*size(p%alloyspecies(i)%atomic_symbol)
            if ( allocated(p%alloyspecies(i)%concentration) ) mem=mem+storage_size(p%alloyspecies(i)%concentration)*size(p%alloyspecies(i)%concentration)
            if ( allocated(p%alloyspecies(i)%atomic_number) ) mem=mem+storage_size(p%alloyspecies(i)%atomic_number)*size(p%alloyspecies(i)%atomic_number)
            if ( allocated(p%alloyspecies(i)%mass         ) ) mem=mem+storage_size(p%alloyspecies(i)%mass         )*size(p%alloyspecies(i)%mass         )
        enddo
    endif
    ! take care of the bz
    if ( allocated(p%bz%highsymmetrypoints) ) mem=mem+storage_size(p%bz%highsymmetrypoints)*size(p%bz%highsymmetrypoints)
    if ( allocated(p%bz%label             ) ) mem=mem+storage_size(p%bz%label             )*size(p%bz%label             )
    if ( allocated(p%bz%r                 ) ) mem=mem+storage_size(p%bz%r                 )*size(p%bz%r                 )
    if ( allocated(p%bz%face) ) then
        do i=1,size(p%bz%face)
            mem=mem+storage_size(p%bz%face(i))
            if ( allocated(p%bz%face(i)%ind) ) mem=mem+storage_size(p%bz%face(i)%ind)*size(p%bz%face(i)%ind)
        enddo
    endif
    if ( allocated(p%bz%node) ) then
        do i=1,size(p%bz%node)
            mem=mem+storage_size(p%bz%node(i))
            if ( allocated(p%bz%node(i)%faceind        ) ) mem=mem+storage_size(p%bz%node(i)%faceind        )*size(p%bz%node(i)%faceind        )
            if ( allocated(p%bz%node(i)%neighbourind   ) ) mem=mem+storage_size(p%bz%node(i)%neighbourind   )*size(p%bz%node(i)%neighbourind   )
            if ( allocated(p%bz%node(i)%neighbourvector) ) mem=mem+storage_size(p%bz%node(i)%neighbourvector)*size(p%bz%node(i)%neighbourvector)
        enddo
    endif
    if ( allocated(p%bz%edge) ) then
        mem=mem+storage_size(p%bz%edge)*size(p%bz%edge)
    endif
    ! irreducible wedge?
    if ( allocated(p%irrw%label) ) mem=mem+storage_size(p%irrw%label)*size(p%irrw%label)
    if ( allocated(p%irrw%r    ) ) mem=mem+storage_size(p%irrw%r    )*size(p%irrw%r    )
    if ( allocated(p%irrw%face  ) ) then
        do i=1,size(p%irrw%face)
            mem=mem+storage_size(p%irrw%face(i))
            if ( allocated(p%irrw%face(i)%ind) ) mem=mem+storage_size(p%irrw%face(i)%ind)*size(p%irrw%face(i)%ind)
        enddo
    endif

    ! get it to bytes
    mem=mem/8
    ! Finally add the spacegroup
    mem=mem+p%sym%size_in_mem()
end function

!> Mass disorder parameter. This function takes a lo_isotope_distribution and calculates the mass disorder parameter.
pure function mass_disorder_parameter(self) result(g)
    !> the isotope distribution
    class(lo_isotope_distribution), intent(in) :: self
    !> the mass disorder parameter
    real(flyt) :: g
    !
    integer :: i
    real(flyt) :: f

    f=0.0_flyt
    do i=1,self%n
        f=f+self%conc(i)*((self%mass(i)-self%mean_mass)**2)
    enddo
    g=f/(self%mean_mass**2)
end function

!> Returns the natural isotope distribution for the specified element in atomic mass units:
subroutine naturaldistribution(self,symbol)
    !> the isotope distribution
    class(lo_isotope_distribution), intent(out) :: self
    !> Atomic symbol, e.g. "Si"
    character(len=*), intent(in) :: symbol

    ! Make sure they are not already allocated
    if ( allocated(self%conc) ) lo_deallocate(self%conc)
    if ( allocated(self%mass) ) lo_deallocate(self%mass)
    !
    select case(trim(symbol))
        case("H")
            self%n=2
           allocate(self%conc(self%n),self%mass(self%n))
            self%conc(1)=0.9998850000_flyt
            self%conc(2)=0.0001150000_flyt
            self%mass(1)=    1.007825032070_flyt
            self%mass(2)=    2.014101777800_flyt
        case("He")
            self%n=2
           allocate(self%conc(self%n),self%mass(self%n))
            self%conc(1)=0.0000013400_flyt
            self%conc(2)=0.9999986600_flyt
            self%mass(1)=    3.016029319100_flyt
            self%mass(2)=    4.002603254150_flyt
        case("Li")
            self%n=2
           allocate(self%conc(self%n),self%mass(self%n))
            self%conc(1)=0.0759000000_flyt
            self%conc(2)=0.9241000000_flyt
            self%mass(1)=    6.015122795000_flyt
            self%mass(2)=    7.016004550000_flyt
        case("Be")
            self%n=1
           allocate(self%conc(self%n),self%mass(self%n))
            self%conc(1)=1.0000000000_flyt
            self%mass(1)=    9.012182200000_flyt
        case("B")
            self%n=2
           allocate(self%conc(self%n),self%mass(self%n))
            self%conc(1)=0.1990000000_flyt
            self%conc(2)=0.8010000000_flyt
            self%mass(1)=   10.012937000000_flyt
            self%mass(2)=   11.009305400000_flyt
        case("C")
            self%n=2
           allocate(self%conc(self%n),self%mass(self%n))
            self%conc(1)=0.9893000000_flyt
            self%conc(2)=0.0107000000_flyt
            self%mass(1)=   12.000000000000_flyt
            self%mass(2)=   13.003354837800_flyt
        case("N")
            self%n=2
           allocate(self%conc(self%n),self%mass(self%n))
            self%conc(1)=0.9963600000_flyt
            self%conc(2)=0.0036400000_flyt
            self%mass(1)=   14.003074004800_flyt
            self%mass(2)=   15.000108898200_flyt
        case("O")
            self%n=3
           allocate(self%conc(self%n),self%mass(self%n))
            self%conc(1)=0.9975700000_flyt
            self%conc(2)=0.0003800000_flyt
            self%conc(3)=0.0020500000_flyt
            self%mass(1)=   15.994914619560_flyt
            self%mass(2)=   16.999131700000_flyt
            self%mass(3)=   17.999161000000_flyt
        case("F")
            self%n=1
           allocate(self%conc(self%n),self%mass(self%n))
            self%conc(1)=1.0000000000_flyt
            self%mass(1)=   18.998403220000_flyt
        case("Ne")
            self%n=3
           allocate(self%conc(self%n),self%mass(self%n))
            self%conc(1)=0.9048000000_flyt
            self%conc(2)=0.0027000000_flyt
            self%conc(3)=0.0925000000_flyt
            self%mass(1)=   19.992440175400_flyt
            self%mass(2)=   20.993846680000_flyt
            self%mass(3)=   21.991385114000_flyt
        case("Na")
            self%n=1
           allocate(self%conc(self%n),self%mass(self%n))
            self%conc(1)=1.0000000000_flyt
            self%mass(1)=   22.989769280900_flyt
        case("Mg")
            self%n=3
           allocate(self%conc(self%n),self%mass(self%n))
            self%conc(1)=0.7899000000_flyt
            self%conc(2)=0.1000000000_flyt
            self%conc(3)=0.1101000000_flyt
            self%mass(1)=   23.985041700000_flyt
            self%mass(2)=   24.985836920000_flyt
            self%mass(3)=   25.982592929000_flyt
        case("Al")
            self%n=1
           allocate(self%conc(self%n),self%mass(self%n))
            self%conc(1)=1.0000000000_flyt
            self%mass(1)=   26.981538630000_flyt
        case("Si")
            self%n=3
           allocate(self%conc(self%n),self%mass(self%n))
            self%conc(1)=0.9222300000_flyt
            self%conc(2)=0.0468500000_flyt
            self%conc(3)=0.0309200000_flyt
            self%mass(1)=   27.976926532500_flyt
            self%mass(2)=   28.976494700000_flyt
            self%mass(3)=   29.973770170000_flyt
        case("P")
            self%n=1
           allocate(self%conc(self%n),self%mass(self%n))
            self%conc(1)=1.0000000000_flyt
            self%mass(1)=   30.973761630000_flyt
        case("S")
            self%n=4
           allocate(self%conc(self%n),self%mass(self%n))
            self%conc(1)=0.9499000000_flyt
            self%conc(2)=0.0075000000_flyt
            self%conc(3)=0.0425000000_flyt
            self%conc(4)=0.0001000000_flyt
            self%mass(1)=   31.972071000000_flyt
            self%mass(2)=   32.971458760000_flyt
            self%mass(3)=   33.967866900000_flyt
            self%mass(4)=   35.967080760000_flyt
        case("Cl")
            self%n=2
           allocate(self%conc(self%n),self%mass(self%n))
            self%conc(1)=0.7576000000_flyt
            self%conc(2)=0.2424000000_flyt
            self%mass(1)=   34.968852680000_flyt
            self%mass(2)=   36.965902590000_flyt
        case("Ar")
            self%n=3
           allocate(self%conc(self%n),self%mass(self%n))
            self%conc(1)=0.0033650000_flyt
            self%conc(2)=0.0006320000_flyt
            self%conc(3)=0.9960030000_flyt
            self%mass(1)=   35.967545106000_flyt
            self%mass(2)=   37.962732400000_flyt
            self%mass(3)=   39.962383122500_flyt
        case("K")
            self%n=3
           allocate(self%conc(self%n),self%mass(self%n))
            self%conc(1)=0.9325810000_flyt
            self%conc(2)=0.0001170000_flyt
            self%conc(3)=0.0673020000_flyt
            self%mass(1)=   38.963706680000_flyt
            self%mass(2)=   39.963998480000_flyt
            self%mass(3)=   40.961825760000_flyt
        case("Ca")
            self%n=6
           allocate(self%conc(self%n),self%mass(self%n))
            self%conc(1)=0.9694100000_flyt
            self%conc(2)=0.0064700000_flyt
            self%conc(3)=0.0013500000_flyt
            self%conc(4)=0.0208600000_flyt
            self%conc(5)=0.0000400000_flyt
            self%conc(6)=0.0018700000_flyt
            self%mass(1)=   39.962590980000_flyt
            self%mass(2)=   41.958618010000_flyt
            self%mass(3)=   42.958766600000_flyt
            self%mass(4)=   43.955481800000_flyt
            self%mass(5)=   45.953692600000_flyt
            self%mass(6)=   47.952534000000_flyt
        case("Sc")
            self%n=1
           allocate(self%conc(self%n),self%mass(self%n))
            self%conc(1)=1.0000000000_flyt
            self%mass(1)=   44.955911900000_flyt
        case("Ti")
            self%n=5
           allocate(self%conc(self%n),self%mass(self%n))
            self%conc(1)=0.0825000000_flyt
            self%conc(2)=0.0744000000_flyt
            self%conc(3)=0.7372000000_flyt
            self%conc(4)=0.0541000000_flyt
            self%conc(5)=0.0518000000_flyt
            self%mass(1)=   45.952631600000_flyt
            self%mass(2)=   46.951763100000_flyt
            self%mass(3)=   47.947946300000_flyt
            self%mass(4)=   48.947870000000_flyt
            self%mass(5)=   49.944791200000_flyt
        case("V")
            self%n=2
           allocate(self%conc(self%n),self%mass(self%n))
            self%conc(1)=0.0025000000_flyt
            self%conc(2)=0.9975000000_flyt
            self%mass(1)=   49.947158500000_flyt
            self%mass(2)=   50.943959500000_flyt
        case("Cr")
            self%n=4
           allocate(self%conc(self%n),self%mass(self%n))
            self%conc(1)=0.0434500000_flyt
            self%conc(2)=0.8378900000_flyt
            self%conc(3)=0.0950100000_flyt
            self%conc(4)=0.0236500000_flyt
            self%mass(1)=   49.946044200000_flyt
            self%mass(2)=   51.940507500000_flyt
            self%mass(3)=   52.940649400000_flyt
            self%mass(4)=   53.938880400000_flyt
        case("Mn")
            self%n=1
           allocate(self%conc(self%n),self%mass(self%n))
            self%conc(1)=1.0000000000_flyt
            self%mass(1)=   54.938045100000_flyt
        case("Fe")
            self%n=4
           allocate(self%conc(self%n),self%mass(self%n))
            self%conc(1)=0.0584500000_flyt
            self%conc(2)=0.9175400000_flyt
            self%conc(3)=0.0211900000_flyt
            self%conc(4)=0.0028200000_flyt
            self%mass(1)=   53.939610500000_flyt
            self%mass(2)=   55.934937500000_flyt
            self%mass(3)=   56.935394000000_flyt
            self%mass(4)=   57.933275600000_flyt
        case("Co")
            self%n=1
           allocate(self%conc(self%n),self%mass(self%n))
            self%conc(1)=1.0000000000_flyt
            self%mass(1)=   58.933195000000_flyt
        case("Ni")
            self%n=5
           allocate(self%conc(self%n),self%mass(self%n))
            self%conc(1)=0.6807690000_flyt
            self%conc(2)=0.2622310000_flyt
            self%conc(3)=0.0113990000_flyt
            self%conc(4)=0.0363450000_flyt
            self%conc(5)=0.0092560000_flyt
            self%mass(1)=   57.935342900000_flyt
            self%mass(2)=   59.930786400000_flyt
            self%mass(3)=   60.931056000000_flyt
            self%mass(4)=   61.928345100000_flyt
            self%mass(5)=   63.927966000000_flyt
        case("Cu")
            self%n=2
           allocate(self%conc(self%n),self%mass(self%n))
            self%conc(1)=0.6915000000_flyt
            self%conc(2)=0.3085000000_flyt
            self%mass(1)=   62.929597500000_flyt
            self%mass(2)=   64.927789500000_flyt
        case("Zn")
            self%n=5
           allocate(self%conc(self%n),self%mass(self%n))
            self%conc(1)=0.4826800000_flyt
            self%conc(2)=0.2797500000_flyt
            self%conc(3)=0.0410200000_flyt
            self%conc(4)=0.1902400000_flyt
            self%conc(5)=0.0063100000_flyt
            self%mass(1)=   63.929142200000_flyt
            self%mass(2)=   65.926033400000_flyt
            self%mass(3)=   66.927127300000_flyt
            self%mass(4)=   67.924844200000_flyt
            self%mass(5)=   69.925319300000_flyt
        case("Ga")
            self%n=2
           allocate(self%conc(self%n),self%mass(self%n))
            self%conc(1)=0.6010800000_flyt
            self%conc(2)=0.3989200000_flyt
            self%mass(1)=   68.925573600000_flyt
            self%mass(2)=   70.924701300000_flyt
        case("Ge")
            self%n=5
           allocate(self%conc(self%n),self%mass(self%n))
            self%conc(1)=0.2038000000_flyt
            self%conc(2)=0.2731000000_flyt
            self%conc(3)=0.0776000000_flyt
            self%conc(4)=0.3672000000_flyt
            self%conc(5)=0.0783000000_flyt
            self%mass(1)=   69.924247400000_flyt
            self%mass(2)=   71.922075800000_flyt
            self%mass(3)=   72.923458900000_flyt
            self%mass(4)=   73.921177800000_flyt
            self%mass(5)=   75.921402600000_flyt
        case("As")
            self%n=1
           allocate(self%conc(self%n),self%mass(self%n))
            self%conc(1)=1.0000000000_flyt
            self%mass(1)=   74.921596500000_flyt
        case("Se")
            self%n=6
           allocate(self%conc(self%n),self%mass(self%n))
            self%conc(1)=0.0089000000_flyt
            self%conc(2)=0.0937000000_flyt
            self%conc(3)=0.0763000000_flyt
            self%conc(4)=0.2377000000_flyt
            self%conc(5)=0.4961000000_flyt
            self%conc(6)=0.0873000000_flyt
            self%mass(1)=   73.922476400000_flyt
            self%mass(2)=   75.919213600000_flyt
            self%mass(3)=   76.919914000000_flyt
            self%mass(4)=   77.917309100000_flyt
            self%mass(5)=   79.916521300000_flyt
            self%mass(6)=   81.916699400000_flyt
        case("Br")
            self%n=2
           allocate(self%conc(self%n),self%mass(self%n))
            self%conc(1)=0.5069000000_flyt
            self%conc(2)=0.4931000000_flyt
            self%mass(1)=   78.918337100000_flyt
            self%mass(2)=   80.916290600000_flyt
        case("Kr")
            self%n=6
           allocate(self%conc(self%n),self%mass(self%n))
            self%conc(1)=0.0035500000_flyt
            self%conc(2)=0.0228600000_flyt
            self%conc(3)=0.1159300000_flyt
            self%conc(4)=0.1150000000_flyt
            self%conc(5)=0.5698700000_flyt
            self%conc(6)=0.1727900000_flyt
            self%mass(1)=   77.920364800000_flyt
            self%mass(2)=   79.916379000000_flyt
            self%mass(3)=   81.913483600000_flyt
            self%mass(4)=   82.914136000000_flyt
            self%mass(5)=   83.911507000000_flyt
            self%mass(6)=   85.910610730000_flyt
        case("Rb")
            self%n=2
           allocate(self%conc(self%n),self%mass(self%n))
            self%conc(1)=0.7217000000_flyt
            self%conc(2)=0.2783000000_flyt
            self%mass(1)=   84.911789738000_flyt
            self%mass(2)=   86.909180527000_flyt
        case("Sr")
            self%n=4
           allocate(self%conc(self%n),self%mass(self%n))
            self%conc(1)=0.0056000000_flyt
            self%conc(2)=0.0986000000_flyt
            self%conc(3)=0.0700000000_flyt
            self%conc(4)=0.8258000000_flyt
            self%mass(1)=   83.913425000000_flyt
            self%mass(2)=   85.909260200000_flyt
            self%mass(3)=   86.908877100000_flyt
            self%mass(4)=   87.905612100000_flyt
        case("Y")
            self%n=1
           allocate(self%conc(self%n),self%mass(self%n))
            self%conc(1)=1.0000000000_flyt
            self%mass(1)=   88.905848300000_flyt
        case("Zr")
            self%n=5
           allocate(self%conc(self%n),self%mass(self%n))
            self%conc(1)=0.5145000000_flyt
            self%conc(2)=0.1122000000_flyt
            self%conc(3)=0.1715000000_flyt
            self%conc(4)=0.1738000000_flyt
            self%conc(5)=0.0280000000_flyt
            self%mass(1)=   89.904704400000_flyt
            self%mass(2)=   90.905645800000_flyt
            self%mass(3)=   91.905040800000_flyt
            self%mass(4)=   93.906315200000_flyt
            self%mass(5)=   95.908273400000_flyt
        case("Nb")
            self%n=1
           allocate(self%conc(self%n),self%mass(self%n))
            self%conc(1)=1.0000000000_flyt
            self%mass(1)=   92.906378100000_flyt
        case("Mo")
            self%n=7
           allocate(self%conc(self%n),self%mass(self%n))
            self%conc(1)=0.1477000000_flyt
            self%conc(2)=0.0923000000_flyt
            self%conc(3)=0.1590000000_flyt
            self%conc(4)=0.1668000000_flyt
            self%conc(5)=0.0956000000_flyt
            self%conc(6)=0.2419000000_flyt
            self%conc(7)=0.0967000000_flyt
            self%mass(1)=   91.906811000000_flyt
            self%mass(2)=   93.905088300000_flyt
            self%mass(3)=   94.905842100000_flyt
            self%mass(4)=   95.904679500000_flyt
            self%mass(5)=   96.906021500000_flyt
            self%mass(6)=   97.905408200000_flyt
            self%mass(7)=   99.907477000000_flyt
        case("Tc")
            self%n=1
            allocate(self%conc(self%n),self%mass(self%n))
            self%conc(1)=1.0_flyt
            self%mass(1)=96.906_flyt
        case("Ru")
            self%n=7
           allocate(self%conc(self%n),self%mass(self%n))
            self%conc(1)=0.0554000000_flyt
            self%conc(2)=0.0187000000_flyt
            self%conc(3)=0.1276000000_flyt
            self%conc(4)=0.1260000000_flyt
            self%conc(5)=0.1706000000_flyt
            self%conc(6)=0.3155000000_flyt
            self%conc(7)=0.1862000000_flyt
            self%mass(1)=   95.907598000000_flyt
            self%mass(2)=   97.905287000000_flyt
            self%mass(3)=   98.905939300000_flyt
            self%mass(4)=   99.904219500000_flyt
            self%mass(5)=  100.905582100000_flyt
            self%mass(6)=  101.904349300000_flyt
            self%mass(7)=  103.905433000000_flyt
        case("Rh")
            self%n=1
           allocate(self%conc(self%n),self%mass(self%n))
            self%conc(1)=1.0000000000_flyt
            self%mass(1)=  102.905504000000_flyt
        case("Pd")
            self%n=6
           allocate(self%conc(self%n),self%mass(self%n))
            self%conc(1)=0.0102000000_flyt
            self%conc(2)=0.1114000000_flyt
            self%conc(3)=0.2233000000_flyt
            self%conc(4)=0.2733000000_flyt
            self%conc(5)=0.2646000000_flyt
            self%conc(6)=0.1172000000_flyt
            self%mass(1)=  101.905609000000_flyt
            self%mass(2)=  103.904036000000_flyt
            self%mass(3)=  104.905085000000_flyt
            self%mass(4)=  105.903486000000_flyt
            self%mass(5)=  107.903892000000_flyt
            self%mass(6)=  109.905153000000_flyt
        case("Ag")
            self%n=2
           allocate(self%conc(self%n),self%mass(self%n))
            self%conc(1)=0.5183900000_flyt
            self%conc(2)=0.4816100000_flyt
            self%mass(1)=  106.905097000000_flyt
            self%mass(2)=  108.904752000000_flyt
        case("Cd")
            self%n=8
           allocate(self%conc(self%n),self%mass(self%n))
            self%conc(1)=0.0125000000_flyt
            self%conc(2)=0.0089000000_flyt
            self%conc(3)=0.1249000000_flyt
            self%conc(4)=0.1280000000_flyt
            self%conc(5)=0.2413000000_flyt
            self%conc(6)=0.1222000000_flyt
            self%conc(7)=0.2873000000_flyt
            self%conc(8)=0.0749000000_flyt
            self%mass(1)=  105.906459000000_flyt
            self%mass(2)=  107.904184000000_flyt
            self%mass(3)=  109.903002100000_flyt
            self%mass(4)=  110.904178100000_flyt
            self%mass(5)=  111.902757800000_flyt
            self%mass(6)=  112.904401700000_flyt
            self%mass(7)=  113.903358500000_flyt
            self%mass(8)=  115.904756000000_flyt
        case("In")
            self%n=2
           allocate(self%conc(self%n),self%mass(self%n))
            self%conc(1)=0.0429000000_flyt
            self%conc(2)=0.9571000000_flyt
            self%mass(1)=  112.904058000000_flyt
            self%mass(2)=  114.903878000000_flyt
        case("Sn")
            self%n=10
           allocate(self%conc(self%n),self%mass(self%n))
            self%conc(1)=0.0097000000_flyt
            self%conc(2)=0.0066000000_flyt
            self%conc(3)=0.0034000000_flyt
            self%conc(4)=0.1454000000_flyt
            self%conc(5)=0.0768000000_flyt
            self%conc(6)=0.2422000000_flyt
            self%conc(7)=0.0859000000_flyt
            self%conc(8)=0.3258000000_flyt
            self%conc(9)=0.0463000000_flyt
            self%conc(10)=0.0579000000_flyt
            self%mass(1)=  111.904818000000_flyt
            self%mass(2)=  113.902779000000_flyt
            self%mass(3)=  114.903342000000_flyt
            self%mass(4)=  115.901741000000_flyt
            self%mass(5)=  116.902952000000_flyt
            self%mass(6)=  117.901603000000_flyt
            self%mass(7)=  118.903308000000_flyt
            self%mass(8)=  119.902194700000_flyt
            self%mass(9)=  121.903439000000_flyt
            self%mass(10)=  123.905273900000_flyt
        case("Sb")
            self%n=2
           allocate(self%conc(self%n),self%mass(self%n))
            self%conc(1)=0.5721000000_flyt
            self%conc(2)=0.4279000000_flyt
            self%mass(1)=  120.903815700000_flyt
            self%mass(2)=  122.904214000000_flyt
        case("Te")
            self%n=8
           allocate(self%conc(self%n),self%mass(self%n))
            self%conc(1)=0.0009000000_flyt
            self%conc(2)=0.0255000000_flyt
            self%conc(3)=0.0089000000_flyt
            self%conc(4)=0.0474000000_flyt
            self%conc(5)=0.0707000000_flyt
            self%conc(6)=0.1884000000_flyt
            self%conc(7)=0.3174000000_flyt
            self%conc(8)=0.3408000000_flyt
            self%mass(1)=  119.904020000000_flyt
            self%mass(2)=  121.903043900000_flyt
            self%mass(3)=  122.904270000000_flyt
            self%mass(4)=  123.902817900000_flyt
            self%mass(5)=  124.904430700000_flyt
            self%mass(6)=  125.903311700000_flyt
            self%mass(7)=  127.904463100000_flyt
            self%mass(8)=  129.906224400000_flyt
        case("I")
            self%n=1
           allocate(self%conc(self%n),self%mass(self%n))
            self%conc(1)=1.0000000000_flyt
            self%mass(1)=  126.904473000000_flyt
        case("Xe")
            self%n=9
           allocate(self%conc(self%n),self%mass(self%n))
            self%conc(1)=0.0009520000_flyt
            self%conc(2)=0.0008900000_flyt
            self%conc(3)=0.0191020000_flyt
            self%conc(4)=0.2640060000_flyt
            self%conc(5)=0.0407100000_flyt
            self%conc(6)=0.2123240000_flyt
            self%conc(7)=0.2690860000_flyt
            self%conc(8)=0.1043570000_flyt
            self%conc(9)=0.0885730000_flyt
            self%mass(1)=  123.905893000000_flyt
            self%mass(2)=  125.904274000000_flyt
            self%mass(3)=  127.903531300000_flyt
            self%mass(4)=  128.904779400000_flyt
            self%mass(5)=  129.903508000000_flyt
            self%mass(6)=  130.905082400000_flyt
            self%mass(7)=  131.904153500000_flyt
            self%mass(8)=  133.905394500000_flyt
            self%mass(9)=  135.907219000000_flyt
        case("Cs")
            self%n=1
           allocate(self%conc(self%n),self%mass(self%n))
            self%conc(1)=1.0000000000_flyt
            self%mass(1)=  132.905451933000_flyt
        case("Ba")
            self%n=7
           allocate(self%conc(self%n),self%mass(self%n))
            self%conc(1)=0.0010600000_flyt
            self%conc(2)=0.0010100000_flyt
            self%conc(3)=0.0241700000_flyt
            self%conc(4)=0.0659200000_flyt
            self%conc(5)=0.0785400000_flyt
            self%conc(6)=0.1123200000_flyt
            self%conc(7)=0.7169800000_flyt
            self%mass(1)=  129.906320800000_flyt
            self%mass(2)=  131.905061300000_flyt
            self%mass(3)=  133.904508400000_flyt
            self%mass(4)=  134.905688600000_flyt
            self%mass(5)=  135.904575900000_flyt
            self%mass(6)=  136.905827400000_flyt
            self%mass(7)=  137.905247200000_flyt
        case("La")
            self%n=2
           allocate(self%conc(self%n),self%mass(self%n))
            self%conc(1)=0.0009000000_flyt
            self%conc(2)=0.9991000000_flyt
            self%mass(1)=  137.907112000000_flyt
            self%mass(2)=  138.906353300000_flyt
        case("Ce")
            self%n=4
           allocate(self%conc(self%n),self%mass(self%n))
            self%conc(1)=0.0018500000_flyt
            self%conc(2)=0.0025100000_flyt
            self%conc(3)=0.8845000000_flyt
            self%conc(4)=0.1111400000_flyt
            self%mass(1)=  135.907172000000_flyt
            self%mass(2)=  137.905991000000_flyt
            self%mass(3)=  139.905438700000_flyt
            self%mass(4)=  141.909244000000_flyt
        case("Pr")
            self%n=1
           allocate(self%conc(self%n),self%mass(self%n))
            self%conc(1)=1.0000000000_flyt
            self%mass(1)=  140.907652800000_flyt
        case("Nd")
            self%n=7
           allocate(self%conc(self%n),self%mass(self%n))
            self%conc(1)=0.2720000000_flyt
            self%conc(2)=0.1220000000_flyt
            self%conc(3)=0.2380000000_flyt
            self%conc(4)=0.0830000000_flyt
            self%conc(5)=0.1720000000_flyt
            self%conc(6)=0.0570000000_flyt
            self%conc(7)=0.0560000000_flyt
            self%mass(1)=  141.907723300000_flyt
            self%mass(2)=  142.909814300000_flyt
            self%mass(3)=  143.910087300000_flyt
            self%mass(4)=  144.912573600000_flyt
            self%mass(5)=  145.913116900000_flyt
            self%mass(6)=  147.916893000000_flyt
            self%mass(7)=  149.920891000000_flyt
        case("Pm")
            self%n=1
            allocate(self%conc(self%n),self%mass(self%n))
            self%conc(1)=1.0_flyt
            self%mass(1)=144.91_flyt
        case("Sm")
            self%n=7
           allocate(self%conc(self%n),self%mass(self%n))
            self%conc(1)=0.0307000000_flyt
            self%conc(2)=0.1499000000_flyt
            self%conc(3)=0.1124000000_flyt
            self%conc(4)=0.1382000000_flyt
            self%conc(5)=0.0738000000_flyt
            self%conc(6)=0.2675000000_flyt
            self%conc(7)=0.2275000000_flyt
            self%mass(1)=  143.911999000000_flyt
            self%mass(2)=  146.914897900000_flyt
            self%mass(3)=  147.914822700000_flyt
            self%mass(4)=  148.917184700000_flyt
            self%mass(5)=  149.917275500000_flyt
            self%mass(6)=  151.919732400000_flyt
            self%mass(7)=  153.922209300000_flyt
        case("Eu")
            self%n=2
           allocate(self%conc(self%n),self%mass(self%n))
            self%conc(1)=0.4781000000_flyt
            self%conc(2)=0.5219000000_flyt
            self%mass(1)=  150.919850200000_flyt
            self%mass(2)=  152.921230300000_flyt
        case("Gd")
            self%n=7
           allocate(self%conc(self%n),self%mass(self%n))
            self%conc(1)=0.0020000000_flyt
            self%conc(2)=0.0218000000_flyt
            self%conc(3)=0.1480000000_flyt
            self%conc(4)=0.2047000000_flyt
            self%conc(5)=0.1565000000_flyt
            self%conc(6)=0.2484000000_flyt
            self%conc(7)=0.2186000000_flyt
            self%mass(1)=  151.919791000000_flyt
            self%mass(2)=  153.920865600000_flyt
            self%mass(3)=  154.922622000000_flyt
            self%mass(4)=  155.922122700000_flyt
            self%mass(5)=  156.923960100000_flyt
            self%mass(6)=  157.924103900000_flyt
            self%mass(7)=  159.927054100000_flyt
        case("Tb")
            self%n=1
           allocate(self%conc(self%n),self%mass(self%n))
            self%conc(1)=1.0000000000_flyt
            self%mass(1)=  158.925346800000_flyt
        case("Dy")
            self%n=7
           allocate(self%conc(self%n),self%mass(self%n))
            self%conc(1)=0.0005600000_flyt
            self%conc(2)=0.0009500000_flyt
            self%conc(3)=0.0232900000_flyt
            self%conc(4)=0.1888900000_flyt
            self%conc(5)=0.2547500000_flyt
            self%conc(6)=0.2489600000_flyt
            self%conc(7)=0.2826000000_flyt
            self%mass(1)=  155.924283000000_flyt
            self%mass(2)=  157.924409000000_flyt
            self%mass(3)=  159.925197500000_flyt
            self%mass(4)=  160.926933400000_flyt
            self%mass(5)=  161.926798400000_flyt
            self%mass(6)=  162.928731200000_flyt
            self%mass(7)=  163.929174800000_flyt
        case("Ho")
            self%n=1
           allocate(self%conc(self%n),self%mass(self%n))
            self%conc(1)=1.0000000000_flyt
            self%mass(1)=  164.930322100000_flyt
        case("Er")
            self%n=6
           allocate(self%conc(self%n),self%mass(self%n))
            self%conc(1)=0.0013900000_flyt
            self%conc(2)=0.0160100000_flyt
            self%conc(3)=0.3350300000_flyt
            self%conc(4)=0.2286900000_flyt
            self%conc(5)=0.2697800000_flyt
            self%conc(6)=0.1491000000_flyt
            self%mass(1)=  161.928778000000_flyt
            self%mass(2)=  163.929200000000_flyt
            self%mass(3)=  165.930293100000_flyt
            self%mass(4)=  166.932048200000_flyt
            self%mass(5)=  167.932370200000_flyt
            self%mass(6)=  169.935464300000_flyt
        case("Tm")
            self%n=1
           allocate(self%conc(self%n),self%mass(self%n))
            self%conc(1)=1.0000000000_flyt
            self%mass(1)=  168.934213300000_flyt
        case("Yb")
            self%n=7
           allocate(self%conc(self%n),self%mass(self%n))
            self%conc(1)=0.0013000000_flyt
            self%conc(2)=0.0304000000_flyt
            self%conc(3)=0.1428000000_flyt
            self%conc(4)=0.2183000000_flyt
            self%conc(5)=0.1613000000_flyt
            self%conc(6)=0.3183000000_flyt
            self%conc(7)=0.1276000000_flyt
            self%mass(1)=  167.933897000000_flyt
            self%mass(2)=  169.934761800000_flyt
            self%mass(3)=  170.936325800000_flyt
            self%mass(4)=  171.936381500000_flyt
            self%mass(5)=  172.938210800000_flyt
            self%mass(6)=  173.938862100000_flyt
            self%mass(7)=  175.942571700000_flyt
        case("Lu")
            self%n=2
           allocate(self%conc(self%n),self%mass(self%n))
            self%conc(1)=0.9741000000_flyt
            self%conc(2)=0.0259000000_flyt
            self%mass(1)=  174.940771800000_flyt
            self%mass(2)=  175.942686300000_flyt
        case("Hf")
            self%n=6
           allocate(self%conc(self%n),self%mass(self%n))
            self%conc(1)=0.0016000000_flyt
            self%conc(2)=0.0526000000_flyt
            self%conc(3)=0.1860000000_flyt
            self%conc(4)=0.2728000000_flyt
            self%conc(5)=0.1362000000_flyt
            self%conc(6)=0.3508000000_flyt
            self%mass(1)=  173.940046000000_flyt
            self%mass(2)=  175.941408600000_flyt
            self%mass(3)=  176.943220700000_flyt
            self%mass(4)=  177.943698800000_flyt
            self%mass(5)=  178.945816100000_flyt
            self%mass(6)=  179.946550000000_flyt
        case("Ta")
            self%n=2
           allocate(self%conc(self%n),self%mass(self%n))
            self%conc(1)=0.0001200000_flyt
            self%conc(2)=0.9998800000_flyt
            self%mass(1)=  179.947464800000_flyt
            self%mass(2)=  180.947995800000_flyt
        case("W")
            self%n=5
           allocate(self%conc(self%n),self%mass(self%n))
            self%conc(1)=0.0012000000_flyt
            self%conc(2)=0.2650000000_flyt
            self%conc(3)=0.1431000000_flyt
            self%conc(4)=0.3064000000_flyt
            self%conc(5)=0.2843000000_flyt
            self%mass(1)=  179.946704000000_flyt
            self%mass(2)=  181.948204200000_flyt
            self%mass(3)=  182.950223000000_flyt
            self%mass(4)=  183.950931200000_flyt
            self%mass(5)=  185.954364100000_flyt
        case("Re")
            self%n=2
           allocate(self%conc(self%n),self%mass(self%n))
            self%conc(1)=0.3740000000_flyt
            self%conc(2)=0.6260000000_flyt
            self%mass(1)=  184.952955000000_flyt
            self%mass(2)=  186.955753100000_flyt
        case("Os")
            self%n=7
           allocate(self%conc(self%n),self%mass(self%n))
            self%conc(1)=0.0002000000_flyt
            self%conc(2)=0.0159000000_flyt
            self%conc(3)=0.0196000000_flyt
            self%conc(4)=0.1324000000_flyt
            self%conc(5)=0.1615000000_flyt
            self%conc(6)=0.2626000000_flyt
            self%conc(7)=0.4078000000_flyt
            self%mass(1)=  183.952489100000_flyt
            self%mass(2)=  185.953838200000_flyt
            self%mass(3)=  186.955750500000_flyt
            self%mass(4)=  187.955838200000_flyt
            self%mass(5)=  188.958147500000_flyt
            self%mass(6)=  189.958447000000_flyt
            self%mass(7)=  191.961480700000_flyt
        case("Ir")
            self%n=2
           allocate(self%conc(self%n),self%mass(self%n))
            self%conc(1)=0.3730000000_flyt
            self%conc(2)=0.6270000000_flyt
            self%mass(1)=  190.960594000000_flyt
            self%mass(2)=  192.962926400000_flyt
        case("Pt")
            self%n=6
           allocate(self%conc(self%n),self%mass(self%n))
            self%conc(1)=0.0001400000_flyt
            self%conc(2)=0.0078200000_flyt
            self%conc(3)=0.3296700000_flyt
            self%conc(4)=0.3383200000_flyt
            self%conc(5)=0.2524200000_flyt
            self%conc(6)=0.0716300000_flyt
            self%mass(1)=  189.959932000000_flyt
            self%mass(2)=  191.961038000000_flyt
            self%mass(3)=  193.962680300000_flyt
            self%mass(4)=  194.964791100000_flyt
            self%mass(5)=  195.964951500000_flyt
            self%mass(6)=  197.967893000000_flyt
        case("Au")
            self%n=1
           allocate(self%conc(self%n),self%mass(self%n))
            self%conc(1)=1.0000000000_flyt
            self%mass(1)=  196.966568700000_flyt
        case("Hg")
            self%n=7
           allocate(self%conc(self%n),self%mass(self%n))
            self%conc(1)=0.0015000000_flyt
            self%conc(2)=0.0997000000_flyt
            self%conc(3)=0.1687000000_flyt
            self%conc(4)=0.2310000000_flyt
            self%conc(5)=0.1318000000_flyt
            self%conc(6)=0.2986000000_flyt
            self%conc(7)=0.0687000000_flyt
            self%mass(1)=  195.965833000000_flyt
            self%mass(2)=  197.966769000000_flyt
            self%mass(3)=  198.968279900000_flyt
            self%mass(4)=  199.968326000000_flyt
            self%mass(5)=  200.970302300000_flyt
            self%mass(6)=  201.970643000000_flyt
            self%mass(7)=  203.973493900000_flyt
        case("Tl")
            self%n=2
           allocate(self%conc(self%n),self%mass(self%n))
            self%conc(1)=0.2952000000_flyt
            self%conc(2)=0.7048000000_flyt
            self%mass(1)=  202.972344200000_flyt
            self%mass(2)=  204.974427500000_flyt
        case("Pb")
            self%n=4
           allocate(self%conc(self%n),self%mass(self%n))
            self%conc(1)=0.0140000000_flyt
            self%conc(2)=0.2410000000_flyt
            self%conc(3)=0.2210000000_flyt
            self%conc(4)=0.5240000000_flyt
            self%mass(1)=  203.973043600000_flyt
            self%mass(2)=  205.974465300000_flyt
            self%mass(3)=  206.975896900000_flyt
            self%mass(4)=  207.976652100000_flyt
        case("Bi")
            self%n=1
           allocate(self%conc(self%n),self%mass(self%n))
            self%conc(1)=1.0000000000_flyt
            self%mass(1)=  208.980398700000_flyt
        case("Po")
            self%n=1
           allocate(self%conc(self%n),self%mass(self%n))
            self%conc(1)=1.0000000000_flyt
            self%mass(1)=208.98_flyt
        case("At")
            self%n=1
           allocate(self%conc(self%n),self%mass(self%n))
            self%conc(1)=1.0000000000_flyt
            self%mass(1)=209.99_flyt
        case("Rn")
            self%n=1
           allocate(self%conc(self%n),self%mass(self%n))
            self%conc(1)=1.0000000000_flyt
            self%mass(1)=222.02_flyt
        case("Fr")
            self%n=1
           allocate(self%conc(self%n),self%mass(self%n))
            self%conc(1)=1.0000000000_flyt
            self%mass(1)=223.02_flyt
        case("Ra")
            self%n=1
           allocate(self%conc(self%n),self%mass(self%n))
            self%conc(1)=1.0000000000_flyt
            self%mass(1)=226.03_flyt
        case("Ac")
            self%n=1
           allocate(self%conc(self%n),self%mass(self%n))
            self%conc(1)=1.0000000000_flyt
            self%mass(1)=227.03_flyt
        case("Th")
            self%n=1
           allocate(self%conc(self%n),self%mass(self%n))
            self%conc(1)=1.0000000000_flyt
            self%mass(1)=  232.038055300000_flyt
        case("Pa")
            self%n=1
           allocate(self%conc(self%n),self%mass(self%n))
            self%conc(1)=1.0000000000_flyt
            self%mass(1)=  231.035884000000_flyt
        case("U")
            self%n=3
           allocate(self%conc(self%n),self%mass(self%n))
            self%conc(1)=0.0000540000_flyt
            self%conc(2)=0.0072040000_flyt
            self%conc(3)=0.9927420000_flyt
            self%mass(1)=  234.040952100000_flyt
            self%mass(2)=  235.043929900000_flyt
            self%mass(3)=  238.050788200000_flyt
        case("Np")
            self%n=1
           allocate(self%conc(self%n),self%mass(self%n))
            self%conc(1)=1.0000000000_flyt
            self%mass(1)=237.05_flyt
        case("Pu")
            self%n=1
           allocate(self%conc(self%n),self%mass(self%n))
            self%conc(1)=1.0000000000_flyt
            self%mass(1)=244.06_flyt
        case("Am")
            self%n=1
           allocate(self%conc(self%n),self%mass(self%n))
            self%conc(1)=1.0000000000_flyt
            self%mass(1)=243.06_flyt
        case("Cm")
            self%n=1
           allocate(self%conc(self%n),self%mass(self%n))
            self%conc(1)=1.0000000000_flyt
            self%mass(1)=247.07_flyt
        case("Bk")
            self%n=1
           allocate(self%conc(self%n),self%mass(self%n))
            self%conc(1)=1.0000000000_flyt
            self%mass(1)=247.07_flyt
        case("Cf")
            self%n=1
           allocate(self%conc(self%n),self%mass(self%n))
            self%conc(1)=1.0000000000_flyt
            self%mass(1)=251.08_flyt
        case("Es")
            self%n=1
           allocate(self%conc(self%n),self%mass(self%n))
            self%conc(1)=1.0000000000_flyt
            self%mass(1)=252.08_flyt
        case("Fm")
            self%n=1
           allocate(self%conc(self%n),self%mass(self%n))
            self%conc(1)=1.0000000000_flyt
            self%mass(1)=257.1_flyt
        case("Md")
            self%n=1
           allocate(self%conc(self%n),self%mass(self%n))
            self%conc(1)=1.0000000000_flyt
            self%mass(1)=258.1_flyt
        case("No")
            self%n=1
           allocate(self%conc(self%n),self%mass(self%n))
            self%conc(1)=1.0000000000_flyt
            self%mass(1)=259.1_flyt
        case("Lr")
            self%n=1
           allocate(self%conc(self%n),self%mass(self%n))
            self%conc(1)=1.0000000000_flyt
            self%mass(1)=262.11_flyt
        case default
            call lo_stop_gracefully(['no isotope distributions available for '//trim(symbol)//', are you sure it is a stable element?'],lo_exitcode_param,__FILE__,__LINE__)
    end select

    ! Convert to atomic units
    self%mass=self%mass*lo_amu_to_emu
    ! Get the actual mass
    self%conc=self%conc/sum(self%conc)
    self%mean_mass=sum(self%conc*self%mass)
    self%disorderparameter=self%mass_disorder_parameter()
end subroutine

end module
#include "precompilerdefinitions"
submodule (type_crystalstructure) type_crystalstructure_alloy
implicit none
contains

!> creates the permutation that orders the SQS the same way as the reference
module subroutine alloy_site_permutation(ss,sqs,permutation)
    !> reference ideal structure
    class(lo_crystalstructure), intent(in) :: ss
    !> ideal SQS
    type(lo_crystalstructure), intent(in) :: sqs
    !> permutation, how to order the atoms in the sqs as the atoms in the supercell
    integer, dimension(:), intent(out) :: permutation


    init: block
        ! First some sanity checks
        if ( ss%info%alloy .eqv. .false. ) then
            call lo_stop_gracefully(['The reference structure needs to be an alloy.'],&
                                    lo_exitcode_param,__FILE__,__LINE__)
        endif
        if ( ss%na .ne. sqs%na ) then
            call lo_stop_gracefully(['Need the same number of atoms to generate permutation.'],&
                                    lo_exitcode_param,__FILE__,__LINE__)
        endif
        if ( norm2(ss%latticevectors-sqs%latticevectors) .gt. lo_tol ) then
            call lo_stop_gracefully(['The lattices need to match when generating permutations.'],&
                                    lo_exitcode_param,__FILE__,__LINE__)
        endif
        ! Add more if needed eventually.

    end block init

    ! Naive and slow way of matching. Plenty of ways to make it faster if I have to at some point.
    ! Ideas include, but are not limited to
    ! a) Sort both, then compare quickly. Annoying to sort vectors with PBC though.
    ! b) Verlet box for lookups, O(N) and sensibly fast I think. Maybe doublebox.
    ! Whatever. Just double-loop for now. Should not be time-sensitive in any way, but you never know.
    stupidandslow: block
        real(flyt), dimension(3) :: v0
        real(flyt) :: f0
        integer :: i,j,k
        logical, dimension(:), allocatable :: assigned

        allocate(assigned(ss%na))
        assigned=.false.
        permutation=0

        do i=1,ss%na
            k=0
            do j=1,sqs%na
                if ( assigned(j) ) cycle
                v0=ss%r(:,i)-sqs%r(:,j)
                v0=lo_clean_fractional_coordinates(v0+0.5_flyt)-0.5_flyt
                f0=lo_sqnorm(v0)
                if ( f0 .lt. lo_sqtol ) then
                    k=j
                    exit
                endif
            enddo

            if ( k .eq. 0 ) then
                call lo_stop_gracefully(['Could not find permutation.'],&
                                        lo_exitcode_param,__FILE__,__LINE__)
            else
                permutation(i)=k
            endif
        enddo

        ! Check that nothing is unassigned
        if ( count(permutation==0) .ne. 0 ) then
            call lo_stop_gracefully(['Could not find permutation.'],&
                                    lo_exitcode_param,__FILE__,__LINE__)
        endif
    end block stupidandslow

end subroutine

!> Change the order of the atoms according to some permutation
module subroutine permute_positions(p,permutation,forward)
    !> structure
    class(lo_crystalstructure), intent(inout) :: p
    !> how to permute them
    integer, dimension(:), intent(in) :: permutation
    !> forward or backwards permutation
    logical, intent(in) :: forward

    if ( forward ) then
        p%atomic_number                     = p%atomic_number                   (permutation)
        p%species                           = p%species                         (permutation)
        p%r                                 = p%r                               (:,permutation)
        p%rcart                             = p%rcart                           (:,permutation)
        p%v                                 = p%v                               (:,permutation)
        p%f                                 = p%f                               (:,permutation)
        p%u                                 = p%u                               (:,permutation)
        p%mass                              = p%mass                            (permutation)
        p%invsqrtmass                       = p%invsqrtmass                     (permutation)
        p%inelastic_neutron_cross_section   = p%inelastic_neutron_cross_section (permutation)
        !p%flavor                            = p%flavor                          (permutation)
        if ( allocated(p%info%index_in_unitcell) ) then
            p%info%index_in_unitcell = p%info%index_in_unitcell (permutation)
        endif
        if ( allocated(p%info%cellindex) ) then
            p%info%cellindex = p%info%cellindex(:,permutation)
        endif
    else
        p%atomic_number                   (permutation) = p%atomic_number
        p%species                         (permutation) = p%species
        p%r                               (:,permutation) = p%r
        p%rcart                           (:,permutation) = p%rcart
        p%v                               (:,permutation) = p%v
        p%f                               (:,permutation) = p%f
        p%u                               (:,permutation) = p%u
        p%mass                            (permutation) = p%mass
        p%invsqrtmass                     (permutation) = p%invsqrtmass
        p%inelastic_neutron_cross_section (permutation) = p%inelastic_neutron_cross_section
        !p%flavor                          (permutation) = p%flavor
        if ( allocated(p%info%index_in_unitcell) ) then
            p%info%index_in_unitcell(permutation) = p%info%index_in_unitcell
        endif
        if ( allocated(p%info%cellindex) ) then
            p%info%cellindex(:,permutation) = p%info%cellindex
        endif
    endif

end subroutine

end submodule
#include "precompilerdefinitions"
submodule (type_crystalstructure) type_crystalstructure_atomdata
implicit none
contains

!> Returns the atomic symbol from atomic number
module function z_to_symbol(z_nucleus) result(symbol)
    !> Nuclear charge, e.g. 14
    integer, intent(in) :: z_nucleus
    !> Atomic symbol, e.g. "Si"
    character(len=2) :: symbol

    character(len=2), parameter :: symbols_of_z(103) = [&
        &'H ', 'He', 'Li', 'Be', 'B ', 'C ', 'N ', 'O ', 'F ', 'Ne', &
        &'Na', 'Mg', 'Al', 'Si', 'P ', 'S ', 'Cl', 'Ar', 'K ', 'Ca', &
        &'Sc', 'Ti', 'V ', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', &
        &'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y ', 'Zr', &
        &'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', &
        &'Sb', 'Te', 'I ', 'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd', &
        &'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', &
        &'Lu', 'Hf', 'Ta', 'W ', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', &
        &'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr', 'Ra', 'Ac', 'Th', &
        &'Pa', 'U ', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm', &
        &'Md', 'No', 'Lr']
    symbol = symbols_of_z(z_nucleus)
end function

!> Translates an atomic symbol into atomic number
module function symbol_to_z(symbol) result(z_nucleus)
    !> Atomic symbol, e.g. "Si"
    character(len=*), intent(in) :: symbol
    !> Nuclear charge, e.g. 14
    integer :: z_nucleus
    !
    select case(trim(adjustl(symbol)))
        case('H' )  ;  z_nucleus=1
        case('He')  ;  z_nucleus=2
        case('Li')  ;  z_nucleus=3
        case('Be')  ;  z_nucleus=4
        case('B' )  ;  z_nucleus=5
        case('C' )  ;  z_nucleus=6
        case('N' )  ;  z_nucleus=7
        case('O' )  ;  z_nucleus=8
        case('F' )  ;  z_nucleus=9

        case('Ne')  ;  z_nucleus=10
        case('Na')  ;  z_nucleus=11
        case('Mg')  ;  z_nucleus=12
        case('Al')  ;  z_nucleus=13
        case('Si')  ;  z_nucleus=14
        case('P' )  ;  z_nucleus=15
        case('S' )  ;  z_nucleus=16
        case('Cl')  ;  z_nucleus=17
        case('Ar')  ;  z_nucleus=18
        case('K' )  ;  z_nucleus=19

        case('Ca')  ;  z_nucleus=20
        case('Sc')  ;  z_nucleus=21
        case('Ti')  ;  z_nucleus=22
        case('V' )  ;  z_nucleus=23
        case('Cr')  ;  z_nucleus=24
        case('Mn')  ;  z_nucleus=25
        case('Fe')  ;  z_nucleus=26
        case('Co')  ;  z_nucleus=27
        case('Ni')  ;  z_nucleus=28
        case('Cu')  ;  z_nucleus=29

        case('Zn')  ;  z_nucleus=30
        case('Ga')  ;  z_nucleus=31
        case('Ge')  ;  z_nucleus=32
        case('As')  ;  z_nucleus=33
        case('Se')  ;  z_nucleus=34
        case('Br')  ;  z_nucleus=35
        case('Kr')  ;  z_nucleus=36
        case('Rb')  ;  z_nucleus=37
        case('Sr')  ;  z_nucleus=38
        case('Y' )  ;  z_nucleus=39

        case('Zr')  ;  z_nucleus=40
        case('Nb')  ;  z_nucleus=41
        case('Mo')  ;  z_nucleus=42
        case('Tc')  ;  z_nucleus=43
        case('Ru')  ;  z_nucleus=44
        case('Rh')  ;  z_nucleus=45
        case('Pd')  ;  z_nucleus=46
        case('Ag')  ;  z_nucleus=47
        case('Cd')  ;  z_nucleus=48
        case('In')  ;  z_nucleus=49

        case('Sn')  ;  z_nucleus=50
        case('Sb')  ;  z_nucleus=51
        case('Te')  ;  z_nucleus=52
        case('I' )  ;  z_nucleus=53
        case('Xe')  ;  z_nucleus=54
        case('Cs')  ;  z_nucleus=55
        case('Ba')  ;  z_nucleus=56
        case('La')  ;  z_nucleus=57
        case('Ce')  ;  z_nucleus=58
        case('Pr')  ;  z_nucleus=59

        case('Nd')  ;  z_nucleus=60
        case('Pm')  ;  z_nucleus=61
        case('Sm')  ;  z_nucleus=62
        case('Eu')  ;  z_nucleus=63
        case('Gd')  ;  z_nucleus=64
        case('Tb')  ;  z_nucleus=65
        case('Dy')  ;  z_nucleus=66
        case('Ho')  ;  z_nucleus=67
        case('Er')  ;  z_nucleus=68
        case('Tm')  ;  z_nucleus=69

        case('Yb')  ;  z_nucleus=70
        case('Lu')  ;  z_nucleus=71
        case('Hf')  ;  z_nucleus=72
        case('Ta')  ;  z_nucleus=73
        case('W' )  ;  z_nucleus=74
        case('Re')  ;  z_nucleus=75
        case('Os')  ;  z_nucleus=76
        case('Ir')  ;  z_nucleus=77
        case('Pt')  ;  z_nucleus=78
        case('Au')  ;  z_nucleus=79

        case('Hg')  ;  z_nucleus=80
        case('Tl')  ;  z_nucleus=81
        case('Pb')  ;  z_nucleus=82
        case('Bi')  ;  z_nucleus=83
        case('Po')  ;  z_nucleus=84
        case('At')  ;  z_nucleus=85
        case('Rn')  ;  z_nucleus=86
        case('Fr')  ;  z_nucleus=87
        case('Ra')  ;  z_nucleus=88
        case('Ac')  ;  z_nucleus=89

        case('Th')  ;  z_nucleus=90
        case('Pa')  ;  z_nucleus=91
        case('U' )  ;  z_nucleus=92
        case('Np')  ;  z_nucleus=93
        case('Pu')  ;  z_nucleus=94
        case('Am')  ;  z_nucleus=95
        case('Cm')  ;  z_nucleus=96
        case('Bk')  ;  z_nucleus=97
        case('Cf')  ;  z_nucleus=98
        case('Es')  ;  z_nucleus=89

        case('Fm')  ;  z_nucleus=100
        case('Md')  ;  z_nucleus=101
        case('No')  ;  z_nucleus=102
        case('Lr')  ;  z_nucleus=103
        case default
            write(*,*) 'TROUBLE!!!, no atomic number available for '//trim(adjustl(symbol))
            write(*,*) 'Are you really sure that is an element?'
            stop
    end select
    !
end function symbol_to_z

! !> inelastic neutron cross-section
module function neutron_cross_section(atomic_number) result(xs)
    !> Nuclear charge, e.g. 14
    integer, intent(in) :: atomic_number
    !> inelastic cross section
    real(flyt) :: xs

    select case(trim(z_to_symbol(atomic_number)))
        case("H")
            xs=1.756800_flyt
        case("He")
            xs=1.340000_flyt
        case("Li")
            xs=0.454000_flyt
        case("Be")
            xs=7.630000_flyt
        case("B")
            xs=3.540000_flyt
        case("C")
            xs=5.551000_flyt
        case("N")
            xs=11.010000_flyt
        case("O")
            xs=4.232000_flyt
        case("F")
            xs=4.017000_flyt
        case("Ne")
            xs=2.620000_flyt
        case("Na")
            xs=1.660000_flyt
        case("Mg")
            xs=3.631000_flyt
        case("Al")
            xs=1.495000_flyt
        case("Si")
            xs=2.163000_flyt
        case("P")
            xs=3.307000_flyt
        case("S")
            xs=1.018600_flyt
        case("Cl")
            xs=11.525700_flyt
        case("Ar")
            xs=0.458000_flyt
        case("K")
            xs=1.690000_flyt
        case("Ca")
            xs=2.780000_flyt
        case("Sc")
            xs=19.000000_flyt
        case("Ti")
            xs=1.485000_flyt
        case("V")
            xs=0.018400_flyt
        case("Cr")
            xs=1.660000_flyt
        case("Mn")
            xs=1.750000_flyt
        case("Fe")
            xs=11.220000_flyt
        case("Co")
            xs=0.779000_flyt
        case("Ni")
            xs=13.300000_flyt
        case("Cu")
            xs=7.485000_flyt
        case("Zn")
            xs=4.054000_flyt
        case("Ga")
            xs=6.675000_flyt
        case("Ge")
            xs=8.420000_flyt
        case("As")
            xs=5.440000_flyt
        case("Se")
            xs=7.980000_flyt
        case("Br")
            xs=5.800000_flyt
        case("Kr")
            xs=7.670000_flyt
        case("Rb")
            xs=6.320000_flyt
        case("Sr")
            xs=6.190000_flyt
        case("Y")
            xs=7.550000_flyt
        case("Zr")
            xs=6.440000_flyt
        case("Nb")
            xs=6.253000_flyt
        case("Mo")
            xs=5.670000_flyt
        case("Tc")
            xs=5.800000_flyt
        case("Ru")
            xs=6.210000_flyt
        case("Rh")
            xs=4.340000_flyt
        case("Pd")
            xs=4.390000_flyt
        case("Ag")
            xs=4.407000_flyt
        case("Cd")
            xs=3.040000_flyt
        case("In")
            xs=2.080000_flyt
        case("Sn")
            xs=4.871000_flyt
        case("Sb")
            xs=3.900000_flyt
        case("Te")
            xs=4.230000_flyt
        case("I")
            xs=3.500000_flyt
        case("Xe")
            xs=2.960000_flyt
        case("Cs")
            xs=3.690000_flyt
        case("Ba")
            xs=3.230000_flyt
        case("La")
            xs=8.530000_flyt
        case("Ce")
            xs=2.940000_flyt
        case("Pr")
            xs=2.640000_flyt
        case("Nd")
            xs=7.430000_flyt
        case("Pm")
            xs=1.0_flyt
        case("Sm")
            xs=0.422000_flyt
        case("Eu")
            xs=6.570000_flyt
        case("Gd")
            xs=29.300000_flyt
        case("Tb")
            xs=6.840000_flyt
        case("Dy")
            xs=35.900000_flyt
        case("Ho")
            xs=8.060000_flyt
        case("Er")
            xs=7.630000_flyt
        case("Tm")
            xs=6.280000_flyt
        case("Yb")
            xs=19.420000_flyt
        case("Lu")
            xs=6.530000_flyt
        case("Hf")
            xs=7.600000_flyt
        case("Ta")
            xs=6.000000_flyt
        case("W")
            xs=2.970000_flyt
        case("Re")
            xs=10.600000_flyt
        case("Os")
            xs=14.400000_flyt
        case("Ir")
            xs=14.100000_flyt
        case("Pt")
            xs=11.580000_flyt
        case("Au")
            xs=7.320000_flyt
        case("Hg")
            xs=20.240000_flyt
        case("Tl")
            xs=9.678000_flyt
        case("Pb")
            xs=11.115000_flyt
        case("Bi")
            xs=9.148000_flyt
        case("Po")
            xs=0.000000_flyt
        case("At")
            xs=0.000000_flyt
        case("Rn")
            xs=0.000000_flyt
        case("Fr")
            xs=0.000000_flyt
        case("Ra")
            xs=1.0_flyt
        case("Ac")
            xs=0.000000_flyt
        case("Th")
            xs=13.360000_flyt
        case("Pa")
            xs=10.400000_flyt
        case("U")
            xs=8.903000_flyt
        case("Np")
            xs=14.000000_flyt
        case("Pu")
            xs=1.0_flyt
        case("Am")
            xs=8.700000_flyt
        case("Cm")
            xs=0.000000_flyt
    end select
end function

end submodule
#include "precompilerdefinitions"
submodule (type_crystalstructure) type_crystalstructure_io
use hdf5_wrappers, only: lo_hdf5_helper,lo_h5_store_data,lo_h5_store_attribute
use konstanter, only: lo_volume_bohr_to_A
implicit none
contains

!> Reads a vasp poscar from file
module subroutine readfromfile(p,filename,verbosity)
    !> the crystal structure
    class(lo_crystalstructure), intent(out) :: p
    !> filename to be read
    character(len=*), intent(in) :: filename
    !> verbosity
    integer, intent(in), optional :: verbosity

    real(flyt), dimension(:,:), allocatable :: r, noncollmagmom, alloy_concentrations
    real(flyt), dimension(:), allocatable :: collmagmom
    real(flyt), dimension(3,3) :: m
    integer, dimension(:,:), allocatable :: alloy_components
    integer, dimension(:), allocatable :: atomic_number,alloy_componentcounter
    integer :: n_elem,verb
    character(len=1000), dimension(:), allocatable :: symbols
    character(len=1000) :: header
    logical, dimension(:), allocatable :: cmatom
    logical :: seldyn,alloy,collmag,noncollmag

    if ( present(verbosity) ) then
        verb=verbosity
    else
        verb=0
    endif

    ! Since this is fortran, I have to dummy open the file and count some stuff first.
    ! the only thing this should do is get the number of different elements, and check
    ! if selective dynamics are used.
    dummyopen: block
        integer, parameter :: max_n_elements=20000
        integer :: u,i,j,k
        character(len=10*max_n_elements) :: trams
        !
        u=open_file('in',trim(filename))
            read(u,*) trams
            read(u,*) trams
            read(u,*) trams
            read(u,*) trams
            read(u,*) trams
            read(u,*) trams
            read(u,'(A)') trams
            j=0
            i=0
            ! figure out how many specified elements there are, this is surprisingly robust.
            do
                if ( j > max_n_elements ) exit
                j=j+1
                if ( trams(j:j) .ne. ' ' ) then
                   i=i+1
                   do k=1,max_n_elements
                        if ( trams(j+1:j+1) .eq. ' ' ) exit
                        if ( trams(j+1:j+1) .ne. ' ' ) j=j+1
                        if ( k .eq. max_n_elements ) then
                            write(*,*) 'more than '//tochar(max_n_elements)//' different elements in the structure? really?'
                            stop
                        endif
                   enddo
                endif
            enddo
            n_elem=i
            ! perhaps it is a POSCAR with selective dynamics
            ! I will ignore that, but it's nice if it does not crash
            read(u,*) trams
            seldyn=.false.
            if ( trams(1:1) .eq. 's' .or. trams(1:1) .eq. 'S' ) seldyn=.true.
        close(u)
    end block dummyopen

    ! Then we open the file for reals this time, and get stuff
    readstuff: block
        real(flyt) :: latpar_or_volume,f0
        real(flyt), dimension(3,3) :: im
        real(flyt), dimension(3) :: v
        integer, dimension(max_n_components) :: di
        integer, dimension(n_elem) :: elemcount
        integer :: i,j,k,l,u,na,ii,jj
        character(len=1) :: trams
        character(len=10), dimension(max_n_components) :: dvsp
        character(len=10) :: dumsp
        character(len=2000) :: dumstr
        logical :: cartesian

        u=open_file('in',trim(filename))
            ! the header
            read(u,*) header
            read(u,*) latpar_or_volume
            read(u,*) m(:,1)
            read(u,*) m(:,2)
            read(u,*) m(:,3)
            lo_allocate(symbols(n_elem))
            read(u,*) symbols
            read(u,*) elemcount
            read(u,*) trams
            ! Figure out if input is in fractional or cartesian coordinates.
            if ( trams .eq. 'C' .or. trams .eq. 'c' ) then
                cartesian=.true.
            elseif ( trams .eq. 'D' .or. trams .eq. 'd' ) then
                cartesian=.false.
            else
                call lo_stop_gracefully(['Specify either Cartesian or direct coordinates in '//trim(filename)],lo_exitcode_io)
            endif
            ! Maybe skip a line
            if ( seldyn ) read(u,*) trams

            ! Now we parsed the header, some small figuring out to do, such as how the input was
            ! specified and that stuff.
            alloy=.false.
            do i=1,n_elem
                if ( symbols(i)(1:5) .eq. 'alloy' .or. symbols(i)(1:5) .eq. 'Alloy' .or. symbols(i)(1:5) .eq. 'ALLOY' ) then
                    symbols(i)='ALLOY'
                    alloy=.true.
                endif
            enddo
            collmag=.false.
            do i=1,n_elem
                if ( symbols(i)(1:2) .eq. 'CM' ) then
                    symbols(i)='CM'
                    collmag=.true.
                endif
            enddo
            noncollmag=.false.
            do i=1,n_elem
                if ( symbols(i)(1:3) .eq. 'NCM' ) then
                    symbols(i)='NCM'
                    noncollmag=.true.
                endif
            enddo

            if ( alloy .and. collmag ) then
                call lo_stop_gracefully(['I will bother with combined alloy and magnetic stuff on a rainy day. Very rainy day.'],lo_exitcode_io)
            endif

            na=sum(elemcount)
            ! figure out if I have lattice parameter or volume
            if ( latpar_or_volume .gt. 0.0_flyt ) then
                m=m*latpar_or_volume
            else
                f0=( abs(latpar_or_volume)/abs(lo_determ(m)) )**(1.0_flyt/3.0_flyt)
                m=m*f0
            endif
            im=lo_invert3x3matrix( m )

            ! Some temporary space for all possible variants of things I want to read in
            lo_allocate(r(3,na))
            lo_allocate(atomic_number(na))
            lo_allocate(collmagmom(na))
            lo_allocate(cmatom(na))
            lo_allocate(noncollmagmom(3,na))
            lo_allocate(alloy_concentrations(max_n_components,na))
            lo_allocate(alloy_components(max_n_components,na))
            lo_allocate(alloy_componentcounter(na))

            r=0.0_flyt
            atomic_number=0
            collmagmom=0.0_flyt
            noncollmagmom=0.0_flyt
            cmatom=.false.
            alloy_concentrations=0.0_flyt
            alloy_components=0
            alloy_componentcounter=0

            ! The basic stuff has been figured out now, read the actual positions and stuff
            if ( alloy ) then
                l=0
                do i=1,n_elem
                do j=1,elemcount(i)
                    l=l+1
                    if ( trim(symbols(i)) .eq. 'ALLOY' ) then
                        read(u,'(A)') dumstr ! fetch this line to a buffer
                        read(dumstr,*) v,k
                        ! Quick sanity test
                        if ( k .le. 1 ) then
                            call lo_stop_gracefully(['It makes no sense to have an alloy with one component.'],lo_exitcode_io)
                        endif

                        ! Not the most elegant solution, but it works.
                        select case(k)
                            case(2)
                                read(dumstr,*) v,alloy_componentcounter(l),&
                                               dvsp(1),alloy_concentrations(1,l),&
                                               dvsp(2),alloy_concentrations(2,l)
                            case(3)
                                read(dumstr,*) v,alloy_componentcounter(l),&
                                               dvsp(1),alloy_concentrations(1,l),&
                                               dvsp(2),alloy_concentrations(2,l),&
                                               dvsp(3),alloy_concentrations(3,l)
                            case(4)
                                read(dumstr,*) v,alloy_componentcounter(l),&
                                               dvsp(1),alloy_concentrations(1,l),&
                                               dvsp(2),alloy_concentrations(2,l),&
                                               dvsp(3),alloy_concentrations(3,l),&
                                               dvsp(4),alloy_concentrations(4,l)
                            case(5)
                                read(dumstr,*) v,alloy_componentcounter(l),&
                                               dvsp(1),alloy_concentrations(1,l),&
                                               dvsp(2),alloy_concentrations(2,l),&
                                               dvsp(3),alloy_concentrations(3,l),&
                                               dvsp(4),alloy_concentrations(4,l),&
                                               dvsp(5),alloy_concentrations(5,l)
                            case(6)
                                read(dumstr,*) v,alloy_componentcounter(l),&
                                               dvsp(1),alloy_concentrations(1,l),&
                                               dvsp(2),alloy_concentrations(2,l),&
                                               dvsp(3),alloy_concentrations(3,l),&
                                               dvsp(4),alloy_concentrations(4,l),&
                                               dvsp(5),alloy_concentrations(5,l),&
                                               dvsp(6),alloy_concentrations(6,l)
                            case default
                                call lo_stop_gracefully(['I have not bothered with input for more than six components.'],lo_exitcode_io,__FILE__,__LINE__)
                        end select

                        ! transform component labels to atomic numbers
                        do k=1,alloy_componentcounter(l)
                            alloy_components(k,l)=symbol_to_z( trim(adjustl(dvsp(k))) )
                        enddo
                        ! Some sanity checks right away, the components have to be unique
                        do ii=1,alloy_componentcounter(l)
                            do jj=ii+1,alloy_componentcounter(l)
                                if ( alloy_components(ii,l) .eq. alloy_components(jj,l) ) then
                                    call lo_stop_gracefully(['You can not alloy something with itself.'],lo_exitcode_io,__FILE__,__LINE__)
                                endif
                            enddo
                        enddo
                        ! concentrations have to add up to 1
                        if ( abs(sum(alloy_concentrations(1:alloy_componentcounter(l),l))-1.0_flyt) .gt. lo_tol ) then
                            call lo_stop_gracefully(['Alloy concentrations have to add up to 1.'],lo_exitcode_io,__FILE__,__LINE__)
                        endif
                        ! Sort components by atomic number
                        k=alloy_componentcounter(l)
                        di=0
                        call qsort(alloy_components(1:k,l),di(1:k))
                        alloy_concentrations(1:k,l)=alloy_concentrations(di(1:k),l)
                        ! Make sure the concentrations add up to 1
                        alloy_concentrations(1:k,l)=lo_chop( alloy_concentrations(1:k,l)/sum(alloy_concentrations(1:k,l)) ,lo_tol )
                    else
                        ! no need to bother, just read normally
                        read(u,*) v
                        atomic_number(l)=symbol_to_z( trim(adjustl(symbols(i))) )
                        alloy_componentcounter(l)=1
                        alloy_concentrations(1,l)=1.0_flyt
                        alloy_components(1,l)=symbol_to_z( trim(adjustl(symbols(i))) )
                    endif
                    ! fix cartesian stuff right away
                    if ( cartesian ) then
                        v=matmul(im,v)
                    endif
                    ! make sure the fractional coordinates really are that.
                    r(:,l)=lo_clean_fractional_coordinates(v)
                enddo
                enddo
            elseif ( collmag ) then
                ! Now there are some sort of magnetic moments specified.
                cmatom=.false.
                collmagmom=0.0_flyt
                l=0
                do i=1,n_elem
                do j=1,elemcount(i)
                    l=l+1
                    if ( trim(symbols(i)) .eq. 'CM' ) then
                        ! reading a DLM atom
                        read(u,*) v,dumsp,collmagmom(l)
                        cmatom(l)=.true.
                        atomic_number(l)=symbol_to_z(trim(dumsp))
                    else
                        read(u,*) v
                        atomic_number(l)=symbol_to_z( trim(adjustl(symbols(i))) )
                    endif
                    ! fix cartesian stuff right away
                    if ( cartesian ) then
                        v=matmul(im,v)
                    endif
                    ! make sure the fractional coordinates really are that.
                    r(:,l)=lo_clean_fractional_coordinates(v)
                enddo
                enddo
            elseif ( noncollmag ) then
                ! Non-collinear magnetic moments specified!
                cmatom=.false.
                noncollmagmom=0.0_flyt
                l=0
                do i=1,n_elem
                do j=1,elemcount(i)
                    l=l+1
                    if ( trim(symbols(i)) .eq. 'NCM' ) then
                        ! reading an atom with non-collinear moment
                        read(u,*) v,dumsp,noncollmagmom(:,l)
                        cmatom(l)=.true.
                        atomic_number(l)=symbol_to_z(trim(dumsp))
                    else
                        read(u,*) v
                        atomic_number(l)=symbol_to_z( trim(adjustl(symbols(i))) )
                    endif
                    ! fix cartesian stuff right away
                    if ( cartesian ) then
                        v=matmul(im,v)
                    endif
                    ! make sure the fractional coordinates really are that.
                    r(:,l)=lo_clean_fractional_coordinates(v)
                enddo
                enddo
            else
                ! Not an alloy or dlm, easy.
                ! get the atomic number
                l=0
                do i=1,n_elem
                do j=1,elemcount(i)
                    l=l+1
                    atomic_number(l)=symbol_to_z( trim(adjustl(symbols(i))) )
                enddo
                enddo
                ! now read the positions
                do i=1,na
                    read(u,*) v
                    ! fix cartesian stuff right away
                    if ( cartesian ) then
                        v=matmul(im,v)
                    endif
                    ! make sure the fractional coordinates really are that.
                    r(:,i)=lo_clean_fractional_coordinates(v)
                enddo
            endif
        close(u)
        ! Maybe say that this was mildly successful
        if ( verb .gt. 0 ) then
            write(*,*) 'Parsed POSCAR header, found '//tochar(sum(elemcount))//' atoms.'
        endif
    end block readstuff

    ! Do the real parsing, with all the classification and stuff.
    call p%generate(m,r,atomic_number,enhet=1,verbosity=verb,collmag=collmag,cmatom=cmatom,collmagmom=collmagmom,&
                    noncollmag=noncollmag,noncollmagmom=noncollmagmom,alloy=alloy,&
                    alloy_componentcounter=alloy_componentcounter,alloy_components=alloy_components,&
                    alloy_concentrations=alloy_concentrations)

    ! And some cleanup
    lo_deallocate(symbols)
    lo_deallocate(r)
    lo_deallocate(atomic_number)
    lo_deallocate(cmatom)
    lo_deallocate(collmagmom)
    lo_deallocate(noncollmagmom)
end subroutine

!> Writes a structure to file or stdout. Use 'stdout' as the filename if you want to write it to screen.
module subroutine writetofile(p,filename,output_format,write_velocities,transformationmatrix)
    !> crystal structure
    class(lo_crystalstructure), intent(in) :: p
    !> the filename
    character(len=*), intent(in) :: filename
    !> what format to write in
    integer, intent(in) :: output_format
    !> if velocities should be written. Default false.
    logical, intent(in), optional :: write_velocities
    !> the structure might have the need to get transformed. If so, how was it transformed?
    real(flyt), dimension(3,3), intent(out), optional :: transformationmatrix

    ! Just pass it on the the appropriate routine
    select case(output_format)
        case(1) ! VASP
            call writetofile_vasp(p,filename,write_velocities)
        case(2) ! Abinit
            call writetofile_abinit(p,filename,write_velocities)
        case(3) ! LAMMPS
            if(present(transformationmatrix)) then
              call writetofile_lammps(p,filename,write_velocities,transformationmatrix)
            else
              call writetofile_lammps(p,filename,write_velocities)
            end if
        case(4) ! FHI Aims
            call writetofile_aims(p,filename,write_velocities)
        case(5) ! Siesta 
            call writetofile_siesta(p,filename,write_velocities)
        case default
            call lo_stop_gracefully(['Unknown output format: '//tochar(output_format)],lo_exitcode_io,__FILE__,__LINE__)
    end select

    contains

    ! The code-specific versions
    subroutine writetofile_vasp(p,filename,write_velocities)
        !> crystal structure
        class(lo_crystalstructure), intent(in) :: p
        !> the filename
        character(len=*), intent(in) :: filename
        !> if velocities should be written. Default false.
        logical, intent(in), optional :: write_velocities

        ! local
        integer :: u, i,j,jj
        real(flyt) :: latpar
        character(len=1000) :: opf
        character(len=1000) :: dum
        character(len=4000) :: dum1

        !@TODO MAKE SURE SPECIES ARE IN THE CORRECT ORDER

        ! Write to a file, or stdout.
        if ( filename .eq. 'stdout' ) then
            u=0
        else
            u=open_file('out',trim(filename))
        endif

        write(u,*) trim(p%info%title)
        ! Print header with lattice parameter and lattice vectors
        opf="(1X,3(F20.14,' '))"
        if ( p%info%unitcell_lattice_parameter .gt. 0.0_flyt ) then
            latpar=p%info%unitcell_lattice_parameter
            write(u,"(2X,F20.12)") latpar*lo_bohr_to_A
            write(u,opf) lo_chop(p%latticevectors(:,1)/latpar,lo_sqtol)
            write(u,opf) lo_chop(p%latticevectors(:,2)/latpar,lo_sqtol)
            write(u,opf) lo_chop(p%latticevectors(:,3)/latpar,lo_sqtol)
        else
            latpar=1.0_flyt
            write(u,"(2X,F20.12)") latpar
            write(u,opf) lo_chop(p%latticevectors(:,1)*lo_bohr_to_A,lo_sqtol)
            write(u,opf) lo_chop(p%latticevectors(:,2)*lo_bohr_to_A,lo_sqtol)
            write(u,opf) lo_chop(p%latticevectors(:,3)*lo_bohr_to_A,lo_sqtol)
        endif

        dum=" "
        if ( p%info%alloy ) then
            do i=1,p%nelements
                if ( p%alloyspecies(i)%n .gt. 1 ) then
                    dum=trim(dum)//" ALLOY"
                else
                    dum=trim(dum)//" "//trim(p%atomic_symbol(i))
                endif
            enddo
        else
            do i=1,p%nelements
                dum=trim(dum)//" "//trim(p%atomic_symbol(i))
            enddo
        endif
        write(u,*) trim(dum)
        write(u,*) tochar(p%element_counter)
        write(u,*) 'Direct coordinates'
        if ( p%info%alloy ) then
            !! a bit fiddly to write nice alloy output, but ok
            do i=1,p%na
                jj=p%species(i)

                if ( p%alloyspecies(jj)%n .eq. 1 ) then
                    j=ceiling(log10(p%na*1.0_flyt+0.1_flyt))+1
                    opf="(1X,3(F18.14,' '),' site'I"//tochar(j)//",' species',I2,': ',A2)"
                    write(u,opf) p%r(:,i),i,p%species(i),p%atomic_symbol(p%species(i))
                else
                    dum1=""
                    do j=1,3
                        dum1=trim(dum1)//"   "//tochar(p%r(j,i),14)
                    enddo
                    dum1(1:499)=dum1(2:500)
                    dum1=trim(dum1)//" "//tochar(p%alloyspecies(jj)%n)
                    do j=1,p%alloyspecies(jj)%n
                        dum1=trim(dum1)//" "//trim(p%alloyspecies(jj)%atomic_symbol(j))
                        dum1=trim(dum1)//" "//tochar(p%alloyspecies(jj)%concentration(j),12)
                    enddo
                    write(u,"(1X,A)") trim(dum1)
                endif
            enddo
        elseif ( p%info%collmag ) then
            j=ceiling(log10(p%na*1.0_flyt+0.1_flyt))+1
            opf="(1X,3(F18.14,' '),' site'I"//tochar(j)//",' magmom: ',F8.4)"
            do i=1,p%na
                write(u,opf) p%r(:,i),i,p%mag%collinear_moment(i)
            enddo
        elseif ( p%info%noncollmag ) then
            j=ceiling(log10(p%na*1.0_flyt+0.1_flyt))+1
            opf="(1X,3(F18.14,' '),' site'I"//tochar(j)//",' magmom: ',3(1X,F8.4))"
            do i=1,p%na
                write(u,opf) p%r(:,i),i,p%mag%noncollinear_moment(:,i)
            enddo
        else
            j=ceiling(log10(p%na*1.0_flyt+0.1_flyt))+1
            opf="(1X,3(F18.14,' '),' site'I"//tochar(j)//",' species',I2,': ',A2)"
            do i=1,p%na
                write(u,opf) p%r(:,i),i,p%species(i),p%atomic_symbol(p%species(i))
            enddo
        endif

        ! Maybe write the velocities
        if ( present(write_velocities) ) then
        if ( write_velocities ) then
            write(u,*) ' '
            do i=1,p%na
                write(u,*) p%v(:,i)*lo_velocity_au_to_Afs
            enddo
        endif
        endif

        if ( filename .ne. 'stdout' ) close(u)
    end subroutine

    subroutine writetofile_lammps(p,filename,write_velocities,transformationmatrix)
        !> crystal structure
        class(lo_crystalstructure), intent(in) :: p
        !> filename
        character(len=*), intent(in) :: filename
        !> should be written. Default false.
        logical, intent(in), optional :: write_velocities
        !> optionally, return the magical transformation matrix
        real(flyt), dimension(3,3), intent(out), optional :: transformationmatrix
        !
        real(flyt), dimension(3,3) :: tm,basis
        real(flyt) :: a,b,c,al,be,gm
        real(flyt) :: ax,bx,cx,by,cy,cz
        logical :: writevel
        integer :: i,u,it

        call lo_stop_gracefully(['Native LAMMPS IO was removed, please use external converters.'], 8)

        ! if ( present(write_velocities) ) then
        !     writevel=write_velocities
        ! else
        !     writevel=.false.
        ! endif

        ! ! Get the lattice parameter thingy
        ! call lo_get_axis_angles(p%latticevectors,a,b,c,al,be,gm)

        ! ax=a
        ! bx=lo_chop(b*cos(gm),lo_sqtol)
        ! by=sqrt(b**2-bx**2)
        ! cx=lo_chop( c*cos(be) ,lo_sqtol)
        ! cy=lo_chop( ( dot_product(p%latticevectors(2,:),p%latticevectors(3,:))-bx*cx )/by ,lo_sqtol)
        ! cz=lo_chop( sqrt(c**2-cx**2-cy**2) , lo_sqtol )

        ! ! Convert to Angstrom
        ! ax=ax*lo_bohr_to_A
        ! bx=bx*lo_bohr_to_A
        ! by=by*lo_bohr_to_A
        ! cx=cx*lo_bohr_to_A
        ! cy=cy*lo_bohr_to_A
        ! cz=cz*lo_bohr_to_A

        ! basis=0.0_flyt
        ! basis(:,1)=[ ax , 0.0_flyt, 0.0_flyt]
        ! basis(:,2)=[ bx, by, 0.0_flyt]
        ! basis(:,3)=[ cx, cy, cz ]
        ! basis=lo_chop(basis,lo_sqtol)

        ! ! Return the coordinate transformation to what LAMMPS likes.
        ! if ( present(transformationmatrix) ) then
        !     transformationmatrix=matmul(basis,p%inv_latticevectors)
        ! endif

        ! ! figure out how to convert stuff
        ! tm=matmul(basis,p%inv_latticevectors)

        ! ! and print it
        ! u=open_file('out',trim(filename))
        !     write(u,*) '# Something'
        !     write(u,*) tochar(p%na),' atoms'
        !     write(u,*) tochar(maxval(p%species)),' atom types'
        !     write(u,*) '0',ax,'xlo xhi '
        !     write(u,*) '0',by,'ylo yhi '
        !     write(u,*) '0',cz,'zlo zhi '
        !     write(u,'(a,3f20.8,a)') '#', bx,cx,cy,'xy xz yz'
        !     write(u,*) 'Masses'
        !     do it=1,p%nelements
        !       i = 1
        !       do i=1,p%na
        !          if (p%species(i) /= it) exit
        !       end do
        !       write(u,'(i10,e20.8,2a)') it, lo_chop(p%mass(i),lo_sqtol), ' # ', trim(p%atomic_symbol(it))
        !     end do
        !     write(u,*) 'Atoms'
        !     write(u,*) ''
        !     do i=1,p%na
        !         write(u,'(2i10,3e20.8)') i,p%species(i),lo_chop(matmul(basis,p%r(:,i)),lo_sqtol)
        !     enddo
        ! close(u)

        ! if ( writevel ) then
        !     write(*,*) 'Have not fixed velocity output for LAMMMPS yet'
        ! endif
    end subroutine

    !> Writes a structure to abinit output file or stdout
    subroutine writetofile_abinit(p,filename,write_velocities)
        !> crystal structure
        class(lo_crystalstructure), intent(in) :: p
        !> filename
        character(len=*), intent(in) :: filename
        !> should be written. Default false.
        logical, intent(in), optional :: write_velocities
        ! local
        integer :: i
        !
        real(flyt) :: latpar
        character(len=1000) :: opf
        character(len=10000) :: dum
        integer :: u,l

        ! file or stdout
        if ( filename .eq. 'stdout' ) then
            u=0
        else
            u=open_file('out',filename)
        endif
        ! It is quite neat to figure a 'lattice parameter' of sorts.
        if ( p%info%unitcell_lattice_parameter .gt. 0.0_flyt ) then
            latpar=p%info%unitcell_lattice_parameter
        else
            latpar=1.0_flyt
        endif

        write(u,*) "#"
        write(u,*) "# you still need to fill in converged values for ecut, tolerance, k-points, MD mode, etc..."
        write(u,*) "# header string ", trim(adjustl(p%info%title))
        write(u,*) "# e.g. :"
        write(u,*) "# ngkpt 1 1 1    shiftk 0 0 0"
        write(u,*) "# ecut 20"
        write(u,*) "# toldfe 1.e-12"
        write(u,*) "acell 3*",tochar(latpar,12)
        write(u,*) "rprim "
        write(u,*) lo_choplarge(p%latticevectors(:,1)/latpar,lo_sqtol)
        write(u,*) lo_choplarge(p%latticevectors(:,2)/latpar,lo_sqtol)
        write(u,*) lo_choplarge(p%latticevectors(:,3)/latpar,lo_sqtol)
        dum="znucl "
        do i=1,p%nelements
            dum=trim(dum)//' '//tochar(symbol_to_z(p%atomic_symbol(i)))
        enddo
        write(u,*) trim(dum)
        write(u,*) "natom ",tochar(p%na)
        write(u,*) "ntypat ",tochar(p%nelements)
        dum="typat "
        l=0 ! this is a counter so that the lines don't get silly long
        do i=1,p%na
            l=l+1
            dum=trim(dum)//' '//tochar(p%species(i))
            if ( l .gt. 45 ) then
                dum=trim(dum)//""
                write(u,*) trim(dum)
                dum=" "
                l=0
            endif
        enddo
        write(u,*) trim(dum)
        write(u,*) 'xred'
        if ( p%info%alloy ) then
            call lo_stop_gracefully(['Alloys in abinit output format are not coded yet'],lo_exitcode_param)
        else
            opf="(1X,3(F18.12,' '),'# site:',I5,' species:',I2,' ',A2)"
            do i=1,p%na
                write(u,opf) lo_chop(p%r(:,i),lo_sqtol),i,p%species(i),p%atomic_symbol(p%species(i))
            enddo
        endif

        ! Maybe write the velocities
        if ( present(write_velocities) ) then
        if ( write_velocities ) then
            write(u,*) 'vel'
            do i=1,p%na
                write(u,*) p%v(:,i)
            enddo
        endif
        endif

        if ( filename .ne. 'stdout' ) close(u)
    end subroutine

    !> Writes a structure in AIMS format
    subroutine writetofile_aims(p,filename,write_velocities)
        !> crystal structure
        class(lo_crystalstructure), intent(in) :: p
        !> filename
        character(len=*), intent(in) :: filename
        !> should velocities be written. Default false.
        logical, intent(in), optional :: write_velocities
        ! local
        integer :: i
        !
        logical :: vel
        integer :: u

        if ( present(write_velocities) ) then
            vel=write_velocities
        else
            vel=.false.
        endif

        ! file or stdout
        if ( filename .eq. 'stdout' ) then
            u=0
        else
            u=open_file('out',filename)
        endif

        ! Write velocities in A/ps = A/fc * 1e3
        write(u,*) "# kommentar"
        write(u,"(1X,'lattice_vector',3(1X,E19.12))") p%latticevectors(:,1)*lo_bohr_to_A
        write(u,"(1X,'lattice_vector',3(1X,E19.12))") p%latticevectors(:,2)*lo_bohr_to_A
        write(u,"(1X,'lattice_vector',3(1X,E19.12))") p%latticevectors(:,3)*lo_bohr_to_A
        do i=1,p%na
            write(u,"(1X,'atom_frac',3(1X,E19.12),1X,A)") p%r(:,i),trim(p%atomic_symbol( p%species(i) ))
            if ( vel ) write(u,"(1X,'velocity',3(1X,E19.12))") p%v(:,i)*lo_velocity_au_to_Afs*1E3_flyt
        enddo

        if ( filename .ne. 'stdout' ) close(u)
    end subroutine

    !> Writes a structure to Siest output file or stdout
    subroutine writetofile_siesta(p,filename,write_velocities)
        !> crystal structure
        class(lo_crystalstructure), intent(in) :: p
        !> filename
        character(len=*), intent(in) :: filename
        !> should be written. Default false.
        logical, intent(in), optional :: write_velocities
        ! local
        integer :: i

        real(flyt) :: latpar
        character(len=1000) :: opf
        integer :: u

        ! I don't know how to write velocities in Siesta. Someone should tell me that.
        ! MJV: fixed at least partly, by putting them in XV files (pos, veloc).
        ! not sure you can add them to the main fdf input files as well...
        ! NB: there is a flag to use the XV file MD.UseSaveXV .true. but you have to fill 
        !   the coordinates block anyway, so might as well fill both here.
        !   Main advantage is that you can add the velocities in XV, but for configuration
        !   generation it's not important, you just get forces out.

        ! file or stdout
        if ( filename .eq. 'stdout' ) then
            u=0
        else
            u=open_file('out',filename)
        endif
        ! It is quite neat to figure a 'lattice parameter' of sorts.
        if ( p%info%unitcell_lattice_parameter .gt. 0.0_flyt ) then
            latpar=p%info%unitcell_lattice_parameter
        else
            latpar=1.0_flyt
        endif

        write(u,*) "# General system descriptors"
        write(u,*) "#"
        write(u,*) "SystemName  TDEPConfiguration"
        write(u,*) "SystemLabel ", filename
        write(u,*) "NumberOfAtoms ", tochar(p%na)
        write(u,*) "NumberOfSpecies ", tochar(p%nelements)
        write(u,*) "%block ChemicalSpeciesLabel"
        do i=1,p%nelements
            write(u,*) tochar(i),'  ',tochar(symbol_to_z(p%atomic_symbol(i))),'  ',p%atomic_symbol(i)
        enddo
        write(u,*) "%endblock ChemicalSpeciesLabel"
        write(u,*) "# Lattice, coordinates"

        write(u,*) "LatticeConstant    1.0  Ang"
        write(u,*) "%block LatticeVectors "
        write(u,*) lo_choplarge(p%latticevectors(:,1)*lo_bohr_to_A,lo_sqtol)
        write(u,*) lo_choplarge(p%latticevectors(:,2)*lo_bohr_to_A,lo_sqtol)
        write(u,*) lo_choplarge(p%latticevectors(:,3)*lo_bohr_to_A,lo_sqtol)
        write(u,*) "%endblock LatticeVectors"
        write(u,*) ""
        write(u,*) "# The positions and velocities are in the corresponding XV file"
        write(u,*) "MD.UseSaveXV true"
        write(u,*) ""
        write(u,*) "# SIESTA requires the following block to be present, so we give it 0s"
        write(u,*) ""
        write(u,*) "AtomicCoordinatesFormat   Cartesian"
        write(u,*) "%block AtomicCoordinatesAndAtomicSpecies"
        if ( p%info%alloy ) then
            call lo_stop_gracefully(['Alloys in siesta output format are not coded yet'],lo_exitcode_param)
        else
            opf="(1X,3(F18.12,' '),I2,'    # site:',I5,' species:',' ',A2)"
            do i=1,p%na
                write(u,opf) lo_chop(p%rcart(:,i),lo_sqtol), p%species(i), i, p%atomic_symbol(p%species(i))
                !write(u,opf) lo_chop(p%r(:,i),lo_sqtol),p%species(i),i,p%atomic_symbol(p%species(i))
            enddo
        endif
        write(u,*) "%endblock AtomicCoordinatesAndAtomicSpecies"

        write(u,*) "XC.functional GGA"
        write(u,*) "XC.authors    PBE"
        write(u,*)
 
        write(u,*) "%block kgrid_Monkhorst_Pack"
        write(u,*) "   1   0   0   0."
        write(u,*) "   0   1   0   0."
        write(u,*) "   0   0   1   0."
        write(u,*) "%endblock kgrid_Monkhorst_Pack"
        write(u,*)
        write(u,*) "SCF.Mix density"
        write(u,*) "SCF.DM.Converge  true"
        write(u,*) "SCF.DM.Tolerance 10e-4"
        write(u,*) "SCF.Mixer.Method Pulay"
        write(u,*) "SCF.Mixer.Weight  0.1"
        write(u,*) "SCF.Mixer.History  6"
        write(u,*) "DM.UseSaveDM True"
        write(u,*) "MeshCutoff 600.0 Ry"

        if ( filename .ne. 'stdout' ) close(u)

! add a XV format file as well: Uses bohr for the cartesian positions.
!    what are the units for the velocity?
!
! latice as 3x3 block, with additional 3x3 block of 0s after
! natom
! itype zatom xcoord(3) veloc(3)
!
! NB: I _think_ they use Angstr for the positions, needs to be checked. 
!    for the velocities, no idea what they expect...
        u=open_file('out',filename//".XV")
        write(u,'(3E20.10,2x,3E20.10)') lo_choplarge(p%latticevectors(:,1),lo_sqtol), 0.000000000, 0.000000000, 0.000000000
        write(u,'(3E20.10,2x,3E20.10)') lo_choplarge(p%latticevectors(:,2),lo_sqtol), 0.000000000, 0.000000000, 0.000000000
        write(u,'(3E20.10,2x,3E20.10)') lo_choplarge(p%latticevectors(:,3),lo_sqtol), 0.000000000, 0.000000000, 0.000000000
        write(u,*) tochar(p%na) 
        opf="(1X,2(I5),2X,2(3F18.9,'  '))"
        do i=1,p%na
            write(u,opf) p%species(i), symbol_to_z(p%atomic_symbol(p%species(i))), lo_chop(p%rcart(:,i),lo_sqtol), lo_chop(p%v(:,i),lo_sqtol)
        enddo
    end subroutine
end subroutine

!> Read isotope distribution from file
module subroutine readisotopefromfile(p)
    !> crystal structure
    class(lo_crystalstructure), intent(inout) :: p

    integer :: i,j,u

    ! Destroy the default isotope distribution
    if ( allocated(p%isotope) ) then
        do i=1,p%na
            if ( allocated(p%isotope(i)%conc) ) lo_deallocate(p%isotope(i)%conc)
            if ( allocated(p%isotope(i)%mass) ) lo_deallocate(p%isotope(i)%mass)
        enddo
        lo_deallocate(p%isotope)
    endif

    ! Get the new, desired isotope distribution
    lo_allocate(p%isotope(p%na))
    u=open_file('in','infile.isotopes')
        do i=1,p%na
            read(u,*) p%isotope(i)%n
            lo_allocate(p%isotope(i)%conc(p%isotope(i)%n))
            lo_allocate(p%isotope(i)%mass(p%isotope(i)%n))
            do j=1,p%isotope(i)%n
                read(u,*) p%isotope(i)%conc(j),p%isotope(i)%mass(j)
            enddo
            ! convert to atomic units
            p%isotope(i)%mass=p%isotope(i)%mass*lo_amu_to_emu
        enddo
    close(u)

    ! Update all the other things that need updating.
    do i=1,p%na
        p%isotope(i)%conc=p%isotope(i)%conc/sum(p%isotope(i)%conc)
        p%isotope(i)%mean_mass=sum(p%isotope(i)%conc*p%isotope(i)%mass)
        p%mass(i)=p%isotope(i)%mean_mass
        p%invsqrtmass(i)=1.0_flyt/sqrt(p%mass(i))
        p%isotope(i)%disorderparameter=p%isotope(i)%mass_disorder_parameter()
    enddo
end subroutine

!> write the Brillouin zone to HDF5
module subroutine write_bz_to_hdf5(p,filename,input_id)
    !> the brillouin zone
    class(lo_crystalstructure), intent(in) :: p
    !> the filename
    character(len=*), intent(in) :: filename
    !> in case I want to write it as a part of another file
    integer(HID_T), intent(in), optional :: input_id

    integer :: i
    character(len=1000) :: dum
    type(lo_hdf5_helper) :: h5

    ! Write in some other file, or create a new one
    if ( present(input_id) ) then
        ! now I assume hdf is open, and a file is open
        !@todo insert sanity check that this is really the case
        h5%file_id=input_id
    else
        ! open a new file
        call h5%init(__FILE__,__LINE__)
        call h5%open_file('write',trim(filename))
    endif

    ! store the reciprocal lattice vectors in case I want to plot them too?
    call h5%store_data(p%reciprocal_latticevectors/lo_Bohr_to_A,h5%file_id,'reciprocal_latticevectors',enhet='1/A')

    ! some meta
    call lo_h5_store_attribute(p%bz%nnodes,h5%file_id,'number_of_zone_nodes')
    call lo_h5_store_attribute(p%bz%nfaces,h5%file_id,'number_of_zone_faces')
    ! store the nodes of the full zone
    call lo_h5_store_data(p%bz%r/lo_bohr_to_A,h5%file_id,'zone_nodes',enhet='1/A')
    ! and the faces. Not the prettiest way, but it does not really matter.
    do i=1,p%bz%nfaces
        dum='zone_face_'//tochar(i)
        call lo_h5_store_data(p%bz%face(i)%ind,h5%file_id,trim(dum))
    enddo
    ! and store the wedge
    call lo_h5_store_attribute(p%irrw%nnodes,h5%file_id,'number_of_wedge_nodes')
    call lo_h5_store_attribute(p%irrw%nfaces,h5%file_id,'number_of_wedge_faces')
    ! store the nodes of the full zone
    call lo_h5_store_data(p%irrw%r/lo_bohr_to_A,h5%file_id,'wedge_nodes',enhet='1/A')
    ! and the faces. Not the prettiest way, but it does not really matter.
    do i=1,p%irrw%nfaces
        dum='wedge_face_'//tochar(i)
        call lo_h5_store_data(p%irrw%face(i)%ind,h5%file_id,trim(dum))
    enddo
    ! also, the labels
    dum=''
    do i=1,p%irrw%nnodes
        dum=trim(dum)//' '//trim(adjustl(p%irrw%label(i)))
    enddo
    call lo_h5_store_attribute(trim(adjustl(dum)),h5%file_id,'wedge_node_labels')

    ! maybe close
    if ( present(input_id) .eqv. .false. ) then
        call h5%close_file()
        call h5%destroy()
    endif
end subroutine

!> write structure information to hdf5.
module subroutine write_structure_to_hdf5(p,filename,input_id)
    !> structure
    class(lo_crystalstructure), intent(in) :: p
    !> filename
    character(len=*), intent(in), optional :: filename
    !> in case I want to write it as a part of another file
    integer(HID_T), intent(in), optional :: input_id

    integer :: i
    character(len=2000) :: atomnames
    type(lo_hdf5_helper) :: h5

    ! Write in some other file, or create a new one
    if ( present(input_id) ) then
        ! now I assume hdf is open, and a file is open
        !@todo insert sanity check that this is really the case
        h5%file_id=input_id
    elseif ( present(filename) ) then
        ! open a new file
        call h5%init(__FILE__,__LINE__)
        call h5%open_file('write',trim(filename))
    else
        call lo_stop_gracefully(['Provide filename or input id, but only ony of them.'],lo_exitcode_param,__FILE__,__LINE__)
    endif

    ! Start with basic things.
    call h5%store_data(p%atomic_number,h5%file_id,'atomic_numbers')
    call h5%store_data(p%r,h5%file_id,'fractional_coordinates',enhet='dimensionless',dimensions='atom,xyz')
    call h5%store_data(p%rcart*lo_Bohr_to_A,h5%file_id,'cartesian_coordinates',enhet='A',dimensions='atom,xyz')
    call h5%store_attribute(p%na,h5%file_id,'number_of_atoms')
    call h5%store_attribute(p%volume*lo_volume_bohr_to_A,h5%file_id,'volume_of_cell')
    call h5%store_data(p%reciprocal_latticevectors/lo_Bohr_to_A,h5%file_id,'reciprocal_latticevectors',enhet='1/A')
    call h5%store_data(p%latticevectors/lo_Bohr_to_A,h5%file_id,'latticevectors',enhet='1/A')
    call p%unique_atom_label(atomnames)
    call h5%store_attribute(trim(adjustl(atomnames)),h5%file_id,'unique_atom_labels')

    ! Can add more things later if needed for some reason.

    ! maybe close
    if ( present(input_id) .eqv. .false. ) then
        call h5%close_file()
        call h5%destroy()
    endif
end subroutine

end submodule
#include "precompilerdefinitions"
submodule (type_crystalstructure) type_crystalstructure_symmetry
implicit none
contains

!> Classify the structure in different ways. Will erase any classification that might have existed earlier.
module subroutine classify(p,how,uc,tolerance,refine,timereversal)
    !> crystal structure
    class(lo_crystalstructure), intent(inout) :: p
    !> classify how?
    character(len=*), intent(in) :: how
    !> perhaps a unitcell
    type(lo_crystalstructure), intent(in), optional :: uc
    !> non-default tolerance?
    real(flyt), intent(in), optional :: tolerance
    !> while classifying, refine the cell?
    logical, intent(in), optional :: refine
    !> with the symmetry stuff, consider time-reveral symmetry as well?
    logical, intent(in), optional :: timereversal

    real(flyt) :: tol
    logical :: refinecell

    if ( present(tolerance) ) then
        tol=tolerance
    else
        tol=lo_tol
    endif

    if ( present(refine) ) then
        refinecell=refine
    else
        refinecell=.false.
    endif

    select case(trim(how))
    case('bravais')
        ! figure out the name of the Bravais lattice
        call find_bravais_lattice(p,tolerance=tol,refine=refinecell)
    case('bz')
        ! calculate the Brillouin zone
        if ( p%info%havebravais .eqv. .false. ) call find_bravais_lattice(p)
        call get_brillouin_zone(p)
    case('spacegroup')
        ! Get the spacegroup
        if ( present(timereversal) .eqv. .false. ) then
            call lo_stop_gracefully(['Time reversal symmetry needs to be specified for finding the space group!'],lo_exitcode_param,__FILE__,__LINE__)
        endif
        call p%sym%generate(p%latticevectors,timereversal,p%r,p%species,verbosity=p%info%verbosity,tolerance=tol)
        p%info%decidedtimereversal=.true.
        p%info%havespacegroup=.true.
    case('wedge')
        if ( present(timereversal) .eqv. .false. ) then
            call lo_stop_gracefully(['Time reversal symmetry needs to be specified defining the irreducible wedge.'],lo_exitcode_param,__FILE__,__LINE__)
        endif
        ! Get the irreducible wedge. Need quite some stuff for this.
        if ( p%info%havebravais .eqv. .false. ) call find_bravais_lattice(p)
        if ( p%info%havebz .eqv. .false. ) call get_brillouin_zone(p)
        if ( p%info%havespacegroup .eqv. .false. ) then
            call p%sym%generate(p%latticevectors,timereversal,p%r,p%species,verbosity=p%info%verbosity)
            p%info%decidedtimereversal=.true.
            p%info%havespacegroup=.true.
        endif
        if ( p%info%havewedge .eqv. .false. ) then
            call get_irreducible_wedge(p,timereversal)
            call label_highsymmetry_points(p,timereversal)
        endif
    case('supercell')
        ! Figure out how a supercell and unitcell might be related.
        call classify_unitcell_supercell_pair(uc,p)
    case('irrep')
        ! Get the character table for the spacegroup
        if ( present(timereversal) .eqv. .false. ) then
            call lo_stop_gracefully(['Time reversal symmetry needs to be specified to get the irrep'],lo_exitcode_param,__FILE__,__LINE__)
        endif
        if ( p%info%havebravais .eqv. .false. ) call find_bravais_lattice(p)
        if ( p%info%havebz .eqv. .false. ) call get_brillouin_zone(p)
        if ( p%info%havespacegroup .eqv. .false. ) then
            call p%sym%generate(p%latticevectors,timereversal,p%r,p%species,verbosity=p%info%verbosity)
            p%info%decidedtimereversal=.true.
            p%info%havespacegroup=.true.
        endif
        if ( p%info%havewedge .eqv. .false. ) then
            call get_irreducible_wedge(p,timereversal)
            call label_highsymmetry_points(p,timereversal)
        endif
        call p%sym%get_character_table(p%info%verbosity)
    case default
        call lo_stop_gracefully(['Not a known way to classify things.'],lo_exitcode_param,__FILE__,__LINE__)
    end select
end subroutine

!> Returns the Cartesian coordinates for a high symmetry point from its label
module function coordinate_from_high_symmetry_point_label(p,label,previous) result(qpoint)
    !> the crystal structure
    class(lo_crystalstructure), intent(in) :: p
    !> the label
    character(len=*), intent(in) :: label
    !> The previous point
    real(flyt), dimension(3), intent(in), optional :: previous
    !> the q-point
    real(flyt), dimension(3) :: qpoint
    !
    integer :: i,j
    real(flyt) :: f0,f1

    ! First make sure that all things are initialized properly.
    if ( p%info%havewedge .eqv. .false. ) then
        call lo_stop_gracefully(['Need the wedge to fetch coordinates from labels.'],lo_exitcode_param,__FILE__,__LINE__)
    endif
    if ( p%info%pointslabelled .eqv. .false. ) then
        call lo_stop_gracefully(['Points need to be labelled before they can be fetched'],lo_exitcode_param,__FILE__,__LINE__)
    endif

    if ( present(previous) ) then
        ! find the closest point that matches
        j=0
        f0=lo_huge
        do i=1,p%bz%nhighsymmetrypoints
            if ( trim(label) .eq. trim(p%bz%label(i)) ) then
                ! found a point matching.
                f1=norm2(previous-p%bz%highsymmetrypoints(:,i))
                if ( f1 .lt. f0 ) then
                    ! I want the closest matching point.
                    j=j+1
                    qpoint=p%bz%highsymmetrypoints(:,i)
                    f0=f1
                endif
            endif
        enddo
        !
        if ( j .eq. 0 ) then
            ! I failed
            call lo_stop_gracefully(['Could not find a tabulated value for the point '//trim(label)],lo_exitcode_param,__FILE__,__LINE__)
        else
            ! I'm done
            return
        endif
    else
        ! just get the first point that matches
        do i=1,p% irrw%nnodes
            if ( trim(label) .eq. trim(p%irrw%label(i)) ) then
                ! found the point
                qpoint=p%irrw%r(:,i)
                return
            endif
        enddo
        do i=1,p%bz%nhighsymmetrypoints
            if ( trim(label) .eq. trim(p%bz%label(i)) ) then
                qpoint=p%bz%highsymmetrypoints(:,i)
                return
            endif
        enddo
    endif
    ! if I made it here, I failed
    call lo_stop_gracefully(['Could not find a tabulated value for the point '//trim(label)],lo_exitcode_param,__FILE__,__LINE__)
end function

! Below are not exposed outside this submodule.

!> figure out how a unitcell and a supercell is related
subroutine classify_unitcell_supercell_pair(uc,ss)
    !> unitcell
    type(lo_crystalstructure), intent(in) :: uc
    !> supercell
    type(lo_crystalstructure), intent(inout) :: ss

    real(flyt), dimension(3,3) :: ssm
    integer :: multiple

    ! Maybe I already did this? Does not matter, do it again, it's fast anyway.
    if ( allocated(ss%info%index_in_unitcell) ) deallocate(ss%info%index_in_unitcell)
    if ( allocated(ss%info%cellindex) ) deallocate(ss%info%cellindex)

    ! Start by assuming it is not a supercell
    ss%info%supercell=.false.

    if ( ss%info%verbosity .gt. 0 ) then
        write(*,*) ' '
        write(*,*) 'Classifying crystal structure: matching unit and supercell '
    endif

    ! A series of sanity checks, to make really sure they actually are a unit-supercell pair
    checklatticevectors: block
        real(flyt) :: f0
        integer :: i,j
        logical :: dl

        ! First that the volume of the supercell is an integer multiple of the unitcell
        f0=ss%volume/uc%volume
        if ( abs(f0-anint(f0)) .gt. lo_tol ) then
            call lo_stop_gracefully(['Supercell/unitcell volume ratio is non-integer.'],lo_exitcode_symmetry,__FILE__,__LINE__)
        endif

        ! The number of atoms must also have the same multiple as the volume
        multiple=int(anint(ss%volume/uc%volume))
        if ( uc%na*multiple .ne. ss%na ) then
            call lo_stop_gracefully(['Inconsistent ratio of atoms between supercell and unitcell.'],lo_exitcode_symmetry,__FILE__,__LINE__)
        endif

        ! Get the supercell latticevectors in terms of the unitcell latticevectors, these
        ! must also have integer values
        ssm=matmul(uc%inv_latticevectors,ss%latticevectors)
        dl=.false.
        do j=1,3
        do i=1,3
            if ( abs( ssm(i,j)-anint(ssm(i,j)) ) .gt. lo_tol ) dl=.true.
        enddo
        enddo
        if ( dl ) then
            call lo_stop_gracefully(['Supercell matrix with non-integer values.'],lo_exitcode_symmetry,__FILE__,__LINE__)
        endif
        ! Clean it a little
        ssm=lo_chop(anint(ssm),lo_sqtol)
    end block checklatticevectors

    if ( ss%info%verbosity .gt. 0 ) then
        write(*,*) '... passed the basic sanity checks'
    endif

    ! So at this point I am pretty sure it's a valid supercell.
    ! Now try to match it per atom
    matchatoms: block
        real(flyt), dimension(:,:), allocatable :: ucv,ssv
        integer, dimension(:,:), allocatable :: ssi
        integer, dimension(:), allocatable :: ucctr,ssctr,ssuci
        integer :: i,j
        logical :: dl
        !
        lo_allocate(ucv(5,uc%na))
        lo_allocate(ssv(5,ss%na))
        lo_allocate(ssi(3,ss%na))
        lo_allocate(ucctr(uc%na))
        lo_allocate(ssctr(ss%na))
        lo_allocate(ssuci(ss%na))

        ! Build some vectors that are good for comparisons!
        do i=1,uc%na
            ucv(1:3,i)=lo_chop(lo_clean_fractional_coordinates(uc%r(:,i),lo_sqtol),lo_sqtol)
            ucv(4,i)=uc%atomic_number(i)
            if ( allocated(uc%flavor) ) then
                ucv(5,i)=uc%flavor(i)
            else
                ucv(5,i)=0
            endif
        enddo

        ! Same for the supercell, but with some more stuff
        do i=1,ss%na
            ssi(:,i)=floor(matmul(ssm,ss%r(:,i))+lo_sqtol)
            ssv(1:3,i)=lo_chop( lo_clean_fractional_coordinates( matmul(ssm,ss%r(:,i))-ssi(:,i)*1.0_flyt, lo_sqtol) ,lo_sqtol)
            ssv(4,i)=ss%atomic_number(i)
            if ( allocated(ss%flavor) ) then
                ssv(5,i)=uc%flavor(i)
            else
                ssv(5,i)=0
            endif
        enddo

        ! Now match it
        ucctr=0
        ssctr=0
        ssuci=0
        do i=1,ss%na
            do j=1,uc%na
                if ( sum(abs(ucv(:,j)-ssv(:,i))) .lt. lo_tol ) then
                    ucctr(j)=ucctr(j)+1
                    ssuci(i)=j
                    ssctr(i)=ssctr(i)+1
                endif
            enddo
        enddo

        ! Plenty of sanity checks
        dl=.false.
        do i=1,ss%na
            if ( ssctr(i) .eq. 0 ) then
                write(*,*) 'Atom ',tochar(i),' in the supercell was not matched'
                dl=.true.
            elseif ( ssctr(i) .gt. 1 ) then
                write(*,*) 'Atom ',tochar(i),' in the supercell was matched ',tochar(ssctr(i)),' times'
                dl=.true.
            endif
        enddo
        do i=1,uc%na
            if ( ucctr(i) .ne. multiple ) then
                write(*,*) 'Atom ',tochar(i),' in the unitcell was matched ',tochar(ucctr(i)),' times, I expected ',tochar(multiple)
                dl=.true.
            endif
        enddo
        if ( dl ) stop
        if ( ss%info%verbosity .gt. 0 ) then
            write(*,*) '... matched all the atoms'
        endif

        ! If I made it here, everything is ok I hope.
        ss%info%supercell=.true.
        ss%info%supercellmatrix=int(anint(ssm))
        lo_allocate(ss%info%cellindex(3,ss%na))
        lo_allocate(ss%info%index_in_unitcell(ss%na))
        ss%info%cellindex=ssi+1
        ss%info%index_in_unitcell=ssuci

        ! And some cleanup
        lo_deallocate(ucv)
        lo_deallocate(ssv)
        lo_deallocate(ssi)
        lo_deallocate(ucctr)
        lo_deallocate(ssctr)
        lo_deallocate(ssuci)
    end block matchatoms

    if ( ss%info%verbosity .gt. 0 ) then
        write(*,*) 'Successfully matched unit and supercell'
    endif
end subroutine

!> Generates the Brillouin zone. I just build the Voronoi diagram of the reciprocal lattice.
subroutine get_brillouin_zone(p)
    !> The crystal structure
    type(lo_crystalstructure), intent(inout) :: p
    !
    type(lo_voronoi_diagram) :: voro
    real(flyt), dimension(3,1) :: dumpts
    real(flyt), dimension(3) :: v0
    real(flyt) :: rc
    integer :: i,j,l

    if ( p%info%verbosity .gt. 0 ) then
        write(*,*) ' '
        write(*,*) 'Generating Brillouin zone.'
    endif

    ! get the Voronoi diagram for the reciprocal lattice
    rc=lo_bounding_sphere_of_box(p%reciprocal_latticevectors)*3.0_flyt ! just a little safety margin
    dumpts=0.0_flyt

    call voro%generate(dumpts,p%reciprocal_latticevectors,cutoff=rc,verbosity=p%info%verbosity)
    ! stuff this into the BZ structure
    p%bz%nnodes=voro%cell(1)%nnodes
    p%bz%nfaces=voro%cell(1)%nfaces
    p%bz%nedges=voro%cell(1)%nedges
    lo_allocate(p%bz%node( p%bz%nnodes ))
    lo_allocate(p%bz%face( p%bz%nfaces ))
    lo_allocate(p%bz%edge( p%bz%nedges ))
    lo_allocate(p%bz%r( 3,p%bz%nnodes ))

    do i=1,p%bz%nfaces
        p%bz%face(i)%n=voro%cell(1)%face(i)%n
        lo_allocate(p%bz%face(i)%ind( p%bz%face(i)%n ))
    enddo
    p%bz%r=voro%cell(1)%r
    p%bz%node=voro%cell(1)%node
    p%bz%face=voro%cell(1)%face
    p%bz%edge=voro%cell(1)%edge
    p%bz%rmin=voro%cell(1)%rmin
    p%bz%rmax=voro%cell(1)%rmax
    ! I'm not 100% sure that this assignment is valid, but it seems to work.
    ! It should according to the f2008 standard, if I got it right.
    p%bz%polyhedron=voro%cell(1)%polyhedron

    ! How many high-symmetry points are there in total? (+1 is gamma)
    p%bz%nhighsymmetrypoints=p%bz%nnodes+p%bz%nedges+p%bz%nfaces+1
    ! Get a list of all high-symmetry points
    lo_allocate(p%bz%highsymmetrypoints(3,p%bz%nhighsymmetrypoints))
    ! add all nodes
    l=0
    do i=1,p%bz%nnodes
        l=l+1
        p%bz%highsymmetrypoints(:,l)=p%bz%r(:,i)
    enddo
    ! add centers of all edges
    do i=1,p%bz%nedges
        l=l+1
        p%bz%highsymmetrypoints(:,l)=(p%bz%r( :, p%bz%edge(i)%i1)+p%bz%r( :, p%bz%edge(i)%i2 ))*0.5_flyt
    enddo
    ! and centers of all faces
    do i=1,p%bz%nfaces
        l=l+1
        v0=0.0_flyt
        do j=1,p%bz%face(i)%n
            v0=v0+p%bz%r(:, p%bz%face(i)%ind(j) )
        enddo
        v0=v0/(p%bz%face(i)%n*1.0_flyt)
        p%bz%highsymmetrypoints(:,l)=v0
    enddo
    ! and gamma
    p%bz%highsymmetrypoints(:,l+1)=0.0_flyt

    ! and tell the world we did it
    p%info%havebz=.true.

    if ( p%info%verbosity .gt. 0 ) then
        write(*,*) 'Got Brillouin zone with '//tochar(p%bz%nnodes)//' nodes, '//&
                   tochar(p%bz%nfaces)//' faces, '//tochar(p%bz%nedges)//' edges and '//&
                   tochar(p%bz%nhighsymmetrypoints)//' high symmetry points.'
    endif
end subroutine

!> Will return the irreducible wedge of the BZ. This might be unnecesarily complicated, but it seems to work quite reliably.
subroutine get_irreducible_wedge(p,timereversal)
    !> the crystal structure
    type(lo_crystalstructure), intent(inout) :: p
    !> consider timereversal symmetry
    logical, intent(in) :: timereversal
    !
    real(flyt), dimension(:,:), allocatable :: tetr,untetr
    real(flyt) :: t0
    integer, dimension(:,:), allocatable :: teti,unteti
    integer :: ntet,nuntet

    if ( p%info%verbosity .gt. 0 ) then
        t0=walltime()
        write(*,*) ''
        write(*,*) 'Generating the irreducible wedge'
    endif

    ! Check the prereqs
    if ( p%info%havebravais .eqv. .false. ) then
        call lo_stop_gracefully(['Need the Bravais lattice to get the wedge'],lo_exitcode_param,__FILE__,__LINE__)
    endif
    if ( p%info%havespacegroup .eqv. .false. ) then
        call lo_stop_gracefully(['Need the spacegroup to get the wedge'],lo_exitcode_param,__FILE__,__LINE__)
    endif
    if ( p%info%havebz .eqv. .false. ) then
        call lo_stop_gracefully(['Need the BZ to get the wedge'],lo_exitcode_param,__FILE__,__LINE__)
    endif
    if ( p%sym%timereversal .neqv. timereversal ) then
        call lo_stop_gracefully(['Conflicting information regarding time-reversal symmetry'],lo_exitcode_param,__FILE__,__LINE__)
    endif

    ! First thing to do is to chop the BZ into tetrahedrons
    buildtetrahedrons: block
        real(flyt), dimension(:,:), allocatable :: points
        real(flyt), dimension(3) :: facecentroid
        real(flyt) :: f0
        integer, dimension(:), allocatable :: di
        real(flyt), dimension(3,4) :: t1,t2
        real(flyt), dimension(3) :: v0
        integer :: np,i_gamma,tetcount,i_1,i_2,i_centroid
        integer :: i,j,k,l,ii,jj,kk

        ! fetch the points
        np=p%bz%nhighsymmetrypoints
        lo_allocate(tetr(3,np))
        tetr=p%bz%highsymmetrypoints
        if ( p%info%verbosity .gt. 0 ) write(*,*) '... found '//tochar(np)//' high symmetry points.'

        ! count number of tetrahedrons
        ntet=0
        do i=1,p%bz%nfaces
            ntet=ntet+p%bz%face(i)%n*2
        enddo
        lo_allocate(teti(4,ntet))

        ! need to keep track of where gamma is
        i_gamma=np

        ! Construct the tetrahedrons
        tetcount=0
        do i=1,p%bz%nfaces
            ! Get the points on this face that are along the edge: the nodes and midpoints of edges.
            lo_allocate(points( 3, p%bz%face(i)%n*2 ))
            l=0
            facecentroid=0.0_flyt
            do j=1,p%bz%face(i)%n
                ! add the point
                l=l+1
                points(:,l)=p%bz%r( : , p%bz%face(i)%ind(j) )
                ! and the one on the middle of the path
                l=l+1
                k=lo_index_in_periodic_array(j+1,p%bz%face(i)%n)
                points(:,l)=( p%bz%r( : , p%bz%face(i)%ind(j) ) + p%bz%r( : , p%bz%face(i)%ind(k) ) )*0.5_flyt
                ! and add to the centroid of the face
                facecentroid=facecentroid+p%bz%r(:, p%bz%face(i)%ind(j) )/(p%bz%face(i)%n*1.0_flyt)
            enddo
            ! Anglesort the points. Should not be necessary, but you never know.
            call p%bz%face(i)%plane%anglesort(points)

            ! Find these points in the list of tetrahedron nodes.
            lo_allocate(di( p%bz%face(i)%n*2 ))
            i_centroid=0
            do j=1,np
                ! find the centroid
                v0=facecentroid-tetr(:,j)
                if ( lo_sqnorm(v0) .lt. lo_sqtol ) i_centroid=j
                ! find the others
                do k=1,size(points,2)
                    v0=tetr(:,j)-points(:,k)
                    if ( lo_sqnorm(v0) .lt. lo_sqtol ) di(k)=j
                enddo
            enddo

            ! Build some tetrahedrons, walk along the perimeter of the face
            do j=1,size(di,1)
                tetcount=tetcount+1
                ! get the two points not gamma and the centroid
                i_1=di( j )
                i_2=di( lo_index_in_periodic_array(j+1, size(di,1) ) )
                ! store the tetrahedron
                teti(:,tetcount)=[i_gamma,i_centroid,i_1,i_2]
            enddo

            lo_deallocate(points)
            lo_deallocate(di)
        enddo
        ! sanity checks, first check the number of tetrahedrons
        if ( tetcount .ne. ntet ) then
            call lo_stop_gracefully(['Bad number of tetrahedrons, expected '//tochar(ntet)//' found '//tochar(tetcount)],lo_exitcode_symmetry,__FILE__,__LINE__)
        endif
        ! check that they make up the entire BZ
        f0=0.0_flyt
        do i=1,ntet
            t1=tetr(:,teti(:,i))
            f0=f0+lo_unsigned_tetrahedron_volume(t1)
        enddo
        if ( abs(f0-1.0_flyt/p%volume) .gt. lo_tol ) then
            call lo_stop_gracefully(['Volume of tetrahedrons does not add up to the volume of the BZ'],lo_exitcode_symmetry,__FILE__,__LINE__)
        endif
        ! They can not be degenerate
        do i=1,ntet
        do j=i+1,ntet
            t1=tetr(:,teti(:,i))
            t2=tetr(:,teti(:,j))
            ! To avoid annoying things with the Cray compiler I had to inline this function.
            kk=0
            do ii=1,4
                f0=lo_huge
                do jj=1,4
                    v0=t1(:,ii)-t2(:,jj)
                    f0=min(f0,lo_sqnorm(v0))
                enddo
                if ( f0 .lt. lo_sqtol ) kk=kk+1
            enddo
            if ( kk .eq. 4 ) then
                call lo_stop_gracefully(['Tetrahedron '//tochar(i)//' and '//tochar(j)//' are degenerate'],lo_exitcode_symmetry,__FILE__,__LINE__)
            endif
        enddo
        enddo
        ! The only thing left to check if if there is overlap, but that has never happened yet, so I don't bother.
    end block buildtetrahedrons
    if ( p%info%verbosity .gt. 0 ) write(*,*) '... chopped the BZ into '//tochar(ntet)//' tetrahedrons'

    ! Reduce the tetrahedrons by symmetry
    symmtet: block
        real(flyt), dimension(:,:), allocatable :: tetcentroids,dr1
        real(flyt), dimension(3,6) :: strangepts
        real(flyt), dimension(3,4) :: t1,t2,t3
        real(flyt), dimension(3) :: v0
        real(flyt) :: f0,f1,bsrad
        integer, dimension(:), allocatable :: di1,di2,tetclass,untet
        integer :: i,j,l,o,pt,ii,jj,kk

        lo_allocate(di1(ntet))
        di1=1
        ! Start by findin the unique tetrahedrons
        tetloop1: do i=1,ntet
            t1=tetr(:,teti(:,i))
            do j=i+1,ntet
                do o=1,p%sym%n
                    t2=tetr(:,teti(:,j))
                    ! Rotate the tetrahedron, this used to be a function but the Cray compilers got angry
                    do ii=1,4
                        t2(:,ii)=lo_operate_on_vector(p%sym%op(o),t2(:,ii),reciprocal=.true.)
                    enddo
                    ! Inlined comparison function
                    kk=0
                    do ii=1,4
                        f0=lo_huge
                        do jj=1,4
                            v0=t1(:,ii)-t2(:,jj)
                            f0=min(f0,lo_sqnorm(v0))
                        enddo
                        if ( f0 .lt. lo_sqtol ) kk=kk+1
                    enddo
                    if ( kk .eq. 4 ) then
                        di1(i)=0
                        cycle tetloop1
                    endif
                    if ( timereversal ) then
                        ! Inlined comparison function
                        t3=-t2
                        kk=0
                        do ii=1,4
                            f0=lo_huge
                            do jj=1,4
                                v0=t1(:,ii)-t3(:,jj)
                                f0=min(f0,lo_sqnorm(v0))
                            enddo
                            if ( f0 .lt. lo_sqtol ) kk=kk+1
                        enddo
                        if ( kk .eq. 4 ) then
                            di1(i)=0
                            cycle tetloop1
                        endif
                    endif
                enddo
            enddo
        enddo tetloop1
        ! Store the unique tetrahedrons
        nuntet=sum(di1)
        lo_allocate(untet(nuntet))
        l=0
        do i=1,ntet
            if ( di1(i) .eq. 1 ) then
                l=l+1
                untet(l)=i
            endif
        enddo
        lo_deallocate(di1)
        ! Figure out what kind each tetrahedron is
        lo_allocate(tetclass(ntet))
        tetclass=0
        tetloop2: do i=1,ntet
            t1=tetr(:,teti(:,i))
            do j=1,nuntet
                do o=1,p%sym%n
                    t2=tetr(:,teti(:,untet(j)))
                    do ii=1,4
                        t2(:,ii)=lo_operate_on_vector(p%sym%op(o),t2(:,ii),reciprocal=.true.)
                    enddo
                    t3=-t2

                    ! First comparison, inlined because Cray
                    kk=0
                    do ii=1,4
                        f0=lo_huge
                        do jj=1,4
                            v0=t1(:,ii)-t2(:,jj)
                            f0=min(f0,lo_sqnorm(v0))
                        enddo
                        if ( f0 .lt. lo_sqtol ) kk=kk+1
                    enddo
                    if ( kk .eq. 4 ) then
                        tetclass(i)=j
                        cycle tetloop2
                    endif
                    if ( timereversal ) then
                        ! Second comparison, inlined because Cray
                        kk=0
                        do ii=1,4
                            f0=lo_huge
                            do jj=1,4
                                v0=t1(:,ii)-t3(:,jj)
                                f0=min(f0,lo_sqnorm(v0))
                            enddo
                            if ( f0 .lt. lo_sqtol ) kk=kk+1
                        enddo
                        if ( kk .eq. 4 ) then
                            tetclass(i)=j
                            cycle tetloop2
                        endif
                    endif
                enddo
            enddo
            if ( tetclass(i) .eq. 0 ) then
                write(*,*) 'Failed classifying tetrahedron',i
                stop
            endif
        enddo tetloop2
        if ( p%info%verbosity .gt. 0 ) write(*,*) '... reduced to '//tochar(nuntet)//' tetrahedrons'

        ! The wedge is built by picking one of each of the irreducible tetrahedrons. There are many many ways to
        ! do this, but I only need one, and one that generates a polyhedron that is as convex as possible.
        ! My algorithm seems dumber than it is, but it actually has not failed (so far).

        ! The centroids of each tetrahedron
        lo_allocate(tetcentroids(3,ntet))
        do i=1,ntet
            do j=1,3
                tetcentroids(j,i)=lo_mean( tetr(j,teti(:,i)) )
            enddo
        enddo
        ! some random points
        strangepts(:,1)=[1010, 123  ,34  ]*1.0_flyt
        strangepts(:,2)=[10,   1023 ,34  ]*1.0_flyt
        strangepts(:,3)=[10,   123  ,1034]*1.0_flyt
        strangepts(:,4)=[1010, 1023 ,34  ]*1.0_flyt
        strangepts(:,5)=[1010, 123  ,1034]*1.0_flyt
        strangepts(:,6)=[1110, 1023 ,1036]*1.0_flyt

        ! So, the algorithm is that I pick one random point, and then the irreducible tetrahedrons closest
        ! to this point. Then I measure the bounding sphere of the centroids of the irreducible. The smallest
        ! bounding sphere is the most compact representation of the wedge, at last I hope so.
        lo_allocate(di1(nuntet))
        lo_allocate(di2(nuntet))
        di1=0
        di2=0
        bsrad=lo_huge
        do pt=1,6
            ! For each unique tetrahedron, pick the one closest to the strange point
            do i=1,nuntet
                f0=lo_huge
                do j=1,ntet
                    if ( tetclass(j) .ne. i ) cycle
                    v0=tetcentroids(:,j)-strangepts(:,pt)
                    f1=lo_sqnorm(v0)
                    if ( f1 .lt. f0 ) then
                        di1(i)=j
                        f0=f1
                    endif
                enddo
            enddo
            ! Now calculate the bounding sphere
            f0=0.0_flyt
            do i=1,nuntet
            do j=i+1,nuntet
                v0=tetcentroids(:,di1(i))-tetcentroids(:,di1(j))
                f0=max(f0,lo_sqnorm( v0 ))
            enddo
            enddo
            if ( f0+1E-8_flyt .lt. bsrad ) then
                bsrad=f0
                di2=di1
            endif
        enddo

        ! There is no need to keep all the points anymore, the irreducible tetrahedra are all we need
        ! First reduce the points to the irreducible
        lo_allocate(dr1(3,nuntet*4))
        l=0
        do i=1,nuntet
        do j=1,4
            l=l+1
            dr1(:,l)=tetr(:,teti(j,di2(i)))
        enddo
        enddo
        call lo_return_unique(dr1,untetr,lo_tol)
        ! and fix the tetrahedron indices
        lo_allocate(unteti(4,nuntet))
        unteti=0
        do i=1,nuntet
        do j=1,4
            do l=1,size(untetr,2)
                v0=untetr(:,l)-tetr(:, teti(j,di2(i)) )
                if ( lo_sqnorm(v0) .lt. lo_sqtol ) then
                    unteti(j,i)=l
                    cycle
                endif
            enddo
            if ( unteti(j,i) .eq. 0 ) then
                call lo_stop_gracefully(['Failed finding irreducible tetrahedron index'],lo_exitcode_symmetry,__FILE__,__LINE__)
            endif
        enddo
        enddo
        ! And now I should have the wedge as a reduced list of tetrahedrons
    end block symmtet

    ! With the short list of tetrahedrons, build the actual wedge!
    buildwedge: block
        type(lo_plane) :: plane
        real(flyt), dimension(:,:), allocatable :: dum,dumun,points
        real(flyt), dimension(3,3) :: planepts
        real(flyt), dimension(3) :: v0,v1,v2
        integer, dimension(:), allocatable :: di
        integer :: i,j,k,l,ii,jj

        ! First, store the points
        p%irrw%nnodes=size(untetr,2)
        lo_allocate(p%irrw%r(3,p%irrw%nnodes))
        p%irrw%r=untetr

        ! Now I want to find the planes that enclose these points. This is faster than I thought it would be.
        ! I just construct all possible planes from the unique points, and keep those planes where all the
        ! unique points are on the negative side, or on the plane.
        lo_allocate(dum(4,p%irrw%nnodes**3))
        l=0
        do i=1,p%irrw%nnodes
        do j=1,p%irrw%nnodes
        do k=1,p%irrw%nnodes
            if ( i .eq. j ) cycle
            if ( j .eq. k ) cycle
            if ( k .eq. i ) cycle
            planepts(:,1)=p%irrw%r(:,i)
            planepts(:,2)=p%irrw%r(:,j)
            planepts(:,3)=p%irrw%r(:,k)
            ! first a dummy check that the points are not coplanaer:
            v0=planepts(:,1)-planepts(:,3)
            v1=planepts(:,2)-planepts(:,3)
            v2=lo_cross(v0,v1)
            if ( lo_sqnorm( v2 ) .lt. lo_sqtol ) cycle
            ! construct this plane
            call plane%generate(planepts)
            ! check that all points are on the negative side, or on the plane
            jj=0
            do ii=1,p%irrw%nnodes
                if ( plane%distance_to_point( p%irrw%r(:,ii) ) .lt. lo_tol ) jj=jj+1
            enddo
            ! if all points are below, it is a keeper!
            if ( jj .eq. p%irrw%nnodes ) then
                l=l+1
                dum(1:3,l)=plane%normal
                dum(4,l)=plane%p
            endif
        enddo
        enddo
        enddo

        ! it's enough to have the unique planes
        call lo_return_unique(dum(:,1:l),dumun)
        lo_deallocate(dum)

        ! Now I know the number of faces!
        p%irrw%nfaces=size(dumun,2)
        lo_allocate(p%irrw%face( p%irrw%nfaces ))

        lo_allocate(points(3,p%irrw%nnodes))
        ! Just have to reconstruct the planes again, from the Hessian normal form. With that
        ! I build the polygons that enclose the irreducible wedge.
        do i=1,p%irrw%nfaces
            ! Define the enclosing planes
            v0=dumun(1:3,i)
            v1=dumun(1:3,i)*(-dumun(4,i))
            call p%irrw%face(i)%plane%generate( normal=v0, point=v1 )
            ! count points on this plane
            points=0.0_flyt
            l=0
            do j=1,p%irrw%nnodes
                if ( abs(p%irrw%face(i)%plane%distance_to_point( p%irrw%r(:,j) )) .lt. lo_tol ) then
                    l=l+1
                    points(:,l)=p%irrw%r(:,j)
                endif
            enddo
            ! find the centroid
            v0=0.0_flyt
            do j=1,l
                v0=v0+points(:,j)/l
            enddo
            ! See if any of the points is the centroid
            k=0
            do j=1,l
                if ( lo_sqnorm(v0-points(:,j)) .lt. lo_sqtol ) k=k+1
            enddo
            ! the centroid should not be added.
            l=l-k

            ! make some space
            if ( l .eq. 0 ) then
                write(*,*) 'Warning, bad plane in the irreducible wedge. Something is strange.'
                ! I don't have to kill it since this really only matters for plots and stuff.
                ! but maybe I should. Dunno.
                p%irrw%face(i)%n=0
                cycle
            endif

            ! Store them
            p%irrw%face(i)%n=l
            lo_allocate(p%irrw%face(i)%ind(l))
            lo_allocate(di(l))
            lo_allocate(dum(3,l))
            l=0
            do j=1,p%irrw%nnodes
                if ( abs(p%irrw%face(i)%plane%distance_to_point( p%irrw%r(:,j) )) .lt. lo_tol ) then
                if ( lo_sqnorm(v0-p%irrw%r(:,j)) .gt. lo_sqtol ) then
                    l=l+1
                    p%irrw%face(i)%ind(l)=j
                    dum(:,l)=p%irrw%r(:,j)
                endif
                endif
            enddo
            ! Anglesort to make plots nice
            call p%irrw%face(i)%plane%anglesort(dum,di)
            p%irrw%face(i)%ind=p%irrw%face(i)%ind( di )
            lo_deallocate(di)
            lo_deallocate(dum)
        enddo
    end block buildwedge
    ! And I have the irreducible wedge!
    p%info%havewedge=.true.
    if ( p%info%verbosity .gt. 0 ) then
        write(*,*) 'Generated the irreducible wedge in '//tochar(walltime()-t0)//'s'
    endif
end subroutine

!> Figures out the Bravais family of the lattice and specific member of that family. Converts the basis to an a,b,c,alpha,beta,gamma representation. From this the Bravais family is deduced. Using that information, a standard orientation of that lattice is generated to figure out the space group and the nomenclature of the high symmetry points. Not stable at all.
subroutine find_bravais_lattice(p,tolerance,refine)
    !> The crystal structure
    class(lo_crystalstructure), intent(inout) :: p
    !> Tolerance, if not default
    real(flyt), intent(in), optional :: tolerance
    !> Should the lattice vectors be refined to match their Bravais family?
    logical, intent(in), optional :: refine
    !
    real(flyt) :: tol,breakpoint
    integer :: i
    !
    real(flyt), dimension(3,3) :: pri,con
    real(flyt) :: f0
    character(len=10) :: family,familymember,oldfamily,newfamily
    character(len=50) :: familymemberlongname
    logical :: refinevectors

    ! default tolerance or something else?
    if ( present(tolerance) ) then
        tol=tolerance
    else
        tol=lo_tol
    endif

    ! Should I adjust the lattice vectors to match the Bravais lattice, to high precision?
    if ( present(refine) ) then
        refinevectors=refine
    else
        refinevectors=.false.
    endif

    ! Get the family and stuff
    call find_bravais_family(p%latticevectors,tol,family,familymember,familymemberlongname)


    ! Perhaps refine it properly?
    if ( refinevectors ) then

        ! First check where the break point is, as in what tolerance changes the determined Bravais lattice
        f0=1E-4_flyt ! Very rough tolerance to start with
        call find_bravais_family(p%latticevectors,f0,oldfamily,familymember,familymemberlongname)
        breakpoint=0.0_flyt
        do i=1,50
            f0=f0*0.5_flyt
            call find_bravais_family(p%latticevectors,f0,newfamily,familymember,familymemberlongname)
            if ( oldfamily .ne. newfamily ) then
                breakpoint=f0
                exit
            else
                oldfamily=newfamily
            endif
        enddo

        ! Set the tolerance above the breakpoint, I guess
        ! Not sure what to do about the break point. Could be useful for some heuristics in the future.
        if ( breakpoint .gt. lo_tiny ) then
            if ( p%info%verbosity .gt. 0 ) write(*,*) '... Bravais lattice changes from ',trim(oldfamily),' to ',trim(newfamily), ' at tolerance:',breakpoint
            f0=breakpoint*2
        else
            f0=lo_tol
        endif

        ! Temporarily get the canonical basis, and the transformation to it
        call get_conventional_abc_and_standard_basis(p%latticevectors,family,pri,con,p%info%unitcell_lattice_parameter)
        p%info%unique_primitive_basis=lo_chop(pri,lo_sqtol)
        call find_transformation_between_current_and_conventional(p)
        pri=matmul(pri,p%info%permutation_to_unique)

        ! Write the old lattice vectors
        if ( p%info%verbosity .gt. 0 ) then
             write(*,*) '... input lattice vectors (scaled with lattice parameter):'
             do i=1,3
                write(*,*) p%latticevectors(:,i)/p%info%unitcell_lattice_parameter
             enddo
        endif

        ! Now adjust the latticevectors to match this exactly
        call adjust_latticevectors(p%latticevectors,pri)

        ! Reclassify with the adjusted lattice vectors
        call find_bravais_family(p%latticevectors,f0,newfamily,familymember,familymemberlongname)
        call get_conventional_abc_and_standard_basis(p%latticevectors,family,pri,con,p%info%unitcell_lattice_parameter)
        if ( p%info%verbosity .gt. 0 ) then
            write(*,*) '... slightly refined lattice vectors:'
            do i=1,3
                write(*,*) p%latticevectors(:,i)/p%info%unitcell_lattice_parameter
            enddo
        endif
        ! Add a final touch, to get e.g. sqrt(3)/2 to all digits
        p%latticevectors=lo_chop(p%latticevectors/p%info%unitcell_lattice_parameter,1E-8_flyt)*p%info%unitcell_lattice_parameter

        if ( p%info%verbosity .gt. 0 ) then
            write(*,*) '... completely refined lattice vectors:'
            do i=1,3
                write(*,*) lo_chop(p%latticevectors(:,i)/p%info%unitcell_lattice_parameter,1E-13_flyt)
            enddo
        endif

        ! Now update all the coordinates and things that depend on the latticevectors
        p%inv_latticevectors=lo_chop(lo_invert3x3matrix(p%latticevectors),1E-13_flyt)
        p%reciprocal_latticevectors=lo_chop(lo_reciprocal_basis(p%latticevectors),1E-13_flyt)
        p%inv_reciprocal_latticevectors=lo_chop(lo_invert3x3matrix(p%reciprocal_latticevectors),1E-13_flyt)
    endif

    ! It's neat to store the standardized conventional and primitive basis
    call get_conventional_abc_and_standard_basis(p%latticevectors,family,pri,con,p%info%unitcell_lattice_parameter)
    p%info%bravaislattice=trim(adjustl(familymember))
    p%info%bravaislatticelongname=trim(adjustl(familymemberlongname))
    p%info%unique_primitive_basis=lo_chop(pri,lo_sqtol)
    p%info%unique_conventional_basis=lo_chop(con,lo_sqtol)
    ! It's also convenient to have the transformation from the current format to the
    ! standardized one. It is an annoying transformation+permutation.
    call find_transformation_between_current_and_conventional(p)
    ! Now I have the Bravais lattice
    p%info%havebravais=.true.
    if ( p%info%verbosity .gt. 0 ) then
        write(*,*) 'Found Bravais lattice: ',trim(p%info%bravaislatticelongname),' (',trim(p%info%bravaislattice),')'
    endif
end subroutine

!> Make sure the distances and angles make sense to really high precision
subroutine adjust_latticevectors(basis,reference)
    real(flyt), dimension(3,3), intent(inout) :: basis
    real(flyt), dimension(3,3), intent(in) :: reference
    !
    real(flyt), dimension(3,3) :: grad,m1
    real(flyt) :: step,f0,f1
    integer :: iter,swctr

    ! Set some starting parameters
    m1=basis
    f0=pardist(basis,reference)
    grad=gradient(basis,reference,basis)
    step=1E-12_flyt/lo_frobnorm(grad)
    swctr=1
    ! Actual minimization
    do iter=1,10000
        if ( f0 .lt. 1E-14_flyt ) exit
        ! update the basis
        m1=lo_chop(m1-step*grad,1E-14_flyt)
        ! new distance
        f1=pardist(m1,reference)
        ! Might be satisified here
        if ( f1 .lt. f0 ) then
            step=step*1.1_flyt
        else
            m1=lo_chop(m1+step*grad,1E-14_flyt)
            swctr=swctr+1
            step=1E-12_flyt/lo_frobnorm(grad)
            step=step/sqrt(swctr*1.0_flyt)
            grad=gradient(m1,reference,basis)
        endif
        f0=f1
    enddo
    ! And make sure the volume does not change
    f0=lo_determ(basis)
    f1=lo_determ(m1)
    basis=m1*f0/f1
    contains

    ! Calculate the gradient of the distance difference with a 9-point stencil
    function gradient(m0,m1,m2) result(grad)
        real(flyt), dimension(3,3) :: m0,m1,m2,grad
        !
        integer :: i,j,k
        real(flyt) :: delta
        real(flyt), dimension(3,3) :: dm
        real(flyt), dimension(8) :: sc_wt,sc_dl,sc_val
        !
        delta=1E-8_flyt
        !
        sc_wt=[3,-32,168,-672,672,-168,32,-3]/(840.0_flyt*delta)
        sc_dl=[-4,-3,-2,-1,1,2,3,4]*delta
        !
        do i=1,3
        do j=1,3
            ! Only provide a gradient for non-zero entries in the original
            if ( abs(m2(j,i)) .gt. 1E-10_flyt ) then
                do k=1,8
                    dm=0.0_flyt
                    dm(j,i)=sc_dl(k)
                    sc_val(k)=pardist(m0+dm,m1)
                enddo
                grad(j,i)=sum(sc_wt*sc_val)
            else
                grad(j,i)=0.0_flyt
            endif
        enddo
        enddo
    end function

    !> distance in parameters such as a,b,c,alpha,beta,gamma
    function pardist(m0,m1) result(d)
        real(flyt), dimension(3,3), intent(in) :: m0,m1
        real(flyt) :: d
        !
        real(flyt), dimension(6) :: p0,p1

        ! Calculate the two a,b,c,alpha,beta,gamma
        call lo_get_axis_angles(m0,p0(1),p0(2),p0(3),p0(4),p0(5),p0(6))
        call lo_get_axis_angles(m1,p1(1),p1(2),p1(3),p1(4),p1(5),p1(6))
        ! Squared difference
        p0=(p0-p1)**2
        d=sqrt(sum(p0))
    end function
end subroutine

!> Find the transformation that takes the current definition of the lattice to a standardized one.
subroutine find_transformation_between_current_and_conventional(p)
    type(lo_crystalstructure), intent(inout) :: p

    integer :: i
    real(flyt), dimension(3,3) :: m0,m1
    real(flyt), dimension(3,3,6) :: permutationmatrices
    real(flyt), dimension(6) :: y0,y1
    real(flyt) :: a,b,c,al,be,gm,f0,f1

    ! Premutation matrices
    permutationmatrices(1,:,1)=[1,0,0]*1.0_flyt
    permutationmatrices(2,:,1)=[0,1,0]*1.0_flyt
    permutationmatrices(3,:,1)=[0,0,1]*1.0_flyt

    permutationmatrices(1,:,2)=[1,0,0]*1.0_flyt
    permutationmatrices(2,:,2)=[0,0,1]*1.0_flyt
    permutationmatrices(3,:,2)=[0,1,0]*1.0_flyt

    permutationmatrices(1,:,3)=[0,1,0]*1.0_flyt
    permutationmatrices(2,:,3)=[1,0,0]*1.0_flyt
    permutationmatrices(3,:,3)=[0,0,1]*1.0_flyt

    permutationmatrices(1,:,4)=[0,0,1]*1.0_flyt
    permutationmatrices(2,:,4)=[1,0,0]*1.0_flyt
    permutationmatrices(3,:,4)=[0,1,0]*1.0_flyt

    permutationmatrices(1,:,5)=[0,1,0]*1.0_flyt
    permutationmatrices(2,:,5)=[0,0,1]*1.0_flyt
    permutationmatrices(3,:,5)=[1,0,0]*1.0_flyt

    permutationmatrices(1,:,6)=[0,0,1]*1.0_flyt
    permutationmatrices(2,:,6)=[0,1,0]*1.0_flyt
    permutationmatrices(3,:,6)=[1,0,0]*1.0_flyt

    ! Find the correct permutation.
    m0=p%info%unique_primitive_basis
    call lo_get_axis_angles(m0,a,b,c,al,be,gm)
    y0=[a,b,c,al,be,gm]
    f0=lo_huge
    do i=1,6
        m1=matmul(p%latticevectors,permutationmatrices(:,:,i))
        call lo_get_axis_angles(m1,a,b,c,al,be,gm)
        y1=[a,b,c,al,be,gm]
        f1=sum(abs(y0-y1))
        if ( f1 .lt. f0 ) then
            p%info%permutation_to_unique=permutationmatrices(:,:,i)
            f0=f1
        endif
    enddo
    ! And the transformation
    p%info%transformation_to_unique=matmul(p%inv_latticevectors,p%info%unique_primitive_basis)
end subroutine

!> Given three latticevectors and a tolerance, figure out the Bravais family.
subroutine find_bravais_family(latticevectors,tolerance,family,member,memberlongname)
    !> latticevectors
    real(flyt), dimension(3,3), intent(in) :: latticevectors
    !> the tolerance, in A
    real(flyt), intent(in) :: tolerance
    !> Bravais family
    character(len=10) :: family
    !> Specific member in that family
    character(len=10) :: member
    !> Long name for this familymemeber
    character(len=50) :: memberlongname
    !
    real(flyt), dimension(3,3) :: m0
    real(flyt), dimension(3,3) :: pri,con
    real(flyt) :: a,b,c,al,be,gm
    real(flyt) :: ka,kb,kc,kal,kbe,kgm
    real(flyt) :: ap,bp,cp
    real(flyt) :: a0,b0,c0,al0,be0,gm0
    real(flyt) :: d90,d60,d120,dbcc
    real(flyt) :: f0,f1
    real(flyt) :: tol_dist,tol_angle

    ! some angles to test against
    d90=1.570796326794897_flyt ! 90 degrees in radians
    d60=1.047197551196598_flyt ! 60 degrees in radians
    d120=2.094395102393195_flyt ! 120 degrees in radians
    dbcc=1.910633236249019_flyt ! the bcc angle

    ! Scale the distance and angle tolerances together somehow
    tol_dist=tolerance
    tol_angle=tolerance*180.0_flyt/lo_pi

    ! Get a,b,c,alpha,beta,gamma
    call get_a_lt_b_lt_c_order(latticevectors,m0)
    call lo_get_axis_angles(m0,a0,b0,c0,al0,be0,gm0)

    ! Guess the family first. The trick was to move these around to catch things in the correct order.
    if ( lo_three_equal(a0,b0,c0,tol_dist) .and. lo_three_equal_to_num(al0,be0,gm0,d90,tol_angle) ) then
        ! 1 CUB Simple cubic 1.000000 1.000000 1.000000 ang 90.000000 90.000000 90.000000
        family='CUB'
        member='CUB'
        memberlongname='simple cubic'
    elseif ( lo_three_equal(a0,b0,c0,tol_dist) .and. lo_three_equal_to_num(al0,be0,gm0,d60,tol_angle) ) then
        ! 2 FCC Face centered cubic 0.707107 0.707107 0.707107 ang 60.000000 60.000000 60.000000
        family='FCC'
        member='FCC'
        memberlongname='face centered cubic'
    elseif ( lo_three_equal(a0,b0,c0,tol_dist) .and. lo_three_equal_to_num(al0,be0,gm0,dbcc,tol_angle) ) then
        ! 3 BCC Body centered cubic 0.866025 0.866025 0.866025 ang 109.471221 109.471221 109.471221
        family='BCC'
        member='BCC'
        memberlongname=''
    elseif ( lo_two_equal(a0,b0,c0,tol_dist) .and. lo_three_equal_to_num(al0,be0,gm0,d90,tol_angle) ) then
        ! 4 TET Tetragonal 1.000000 1.000000 1.100000 ang 90.000000 90.000000 90.000000
        family='TET'
        member='TET'
        memberlongname='tetragonal'
    elseif ( lo_three_equal(a0,b0,c0,tol_dist) .and. lo_two_equal(al0,be0,gm0,tol_angle)) then ! .and. lo_three_larger_than_num(a0,be0,gm0,d90,tol_angle) ) then
        ! 5 BCT1 Body-centered tetragonal 0.838153 0.838153 0.838153 ang 106.753588 106.753588 115.054967
        family='BCT'
        call get_conventional_abc_and_standard_basis(latticevectors,family,pri,con)
        call lo_get_axis_angles(con,a,b,c,al,be,gm)
        if ( a .gt. c ) then
            member='BCT1'
        else
            member='BCT2'
        endif
        memberlongname='body-centered tetragonal'
    elseif ( lo_none_equal(a0,b0,c0,tol_dist) .and. lo_three_equal_to_num(al0,be0,gm0,d90,tol_angle) ) then
        ! 7 ORC Orthorhombic 1.000000 1.100000 1.200000 ang 90.000000 90.000000 90.000000
        family='ORC'
        member='ORC'
        memberlongname='orthorhombic'
    elseif ( lo_one_equal_to_num(al0,be0,gm0,d120,tol_dist) .and. lo_two_equal_to_num(al0,be0,gm0,d90,tol_angle) .and. lo_two_equal(a0,b0,c0,tol_dist) ) then
        ! 13 HEX Hexagonal 1.000000 1.000000 1.200000 ang 90.000000 90.000000 120.000000
        family='HEX'
        member='HEX'
        memberlongname='hexagonal'
    elseif ( lo_one_equal_to_num(al0,be0,gm0,d60,tol_dist) .and. lo_two_equal_to_num(al0,be0,gm0,d90,tol_angle) .and. lo_two_equal(a0,b0,c0,tol_dist) ) then
        ! 13 HEX Hexagonal 1.000000 1.000000 1.200000 ang 90.000000 90.000000 60.000000
        family='HEX2'
        member='HEX2'
        memberlongname='hexagonal'
    elseif ( lo_two_equal(a0,b0,c0,tol_dist) .and. lo_two_equal_to_num(al0,be0,gm0,d90,tol_angle) ) then
        ! 12 ORCC C-centered orthorhombic 0.743303 0.743303 1.300000 ang 90.000000 90.000000 95.452622
        family='ORCC'
        call get_conventional_abc_and_standard_basis(latticevectors,family,pri,con)
        call lo_get_axis_angles(pri,a,b,c,al,be,gm)
        if ( gm .lt. lo_pi*0.5_flyt ) then
            member='ORCC1'
        else
            member='ORCC2'
        endif
        memberlongname='c-centered orthorhombic'
    elseif ( lo_none_equal(a0,b0,c0,tol_dist) .and. lo_two_equal_to_num(al0,be0,gm0,d90,tol_angle) ) then
        ! 16 MCL Monoclinic 1.000000 1.100000 1.200000 ang 80.000000 90.000000 90.000000
        family='MCL'
        member='MCL'
        memberlongname='monoclinic'
    elseif ( lo_three_equal(a0,b0,c0,tol_dist) .and. lo_none_equal(al0,be0,gm0,tol_angle) ) then
        ! 11 ORCI Body-centered orthorhombic 0.955249 0.955249 0.955249 ang 116.875594 109.693369 102.178552
        family='ORCI'
        member='ORCI'
        memberlongname='body-centered orthorhombic'
    elseif ( lo_three_equal(a0,b0,c0,tol_dist) .and. lo_three_equal(al0,be0,gm0,tol_angle) ) then
        ! 14 RHL1 Rhombohedral 1.000000 1.000000 1.000000 ang 80.000000 80.000000 80.000000
        family='RHL'
        call get_conventional_abc_and_standard_basis(latticevectors,family,pri,con)
        call lo_get_axis_angles(con,a,b,c,al,be,gm)
        if ( al .lt. d90 ) then
            member='RHL1'
        else
            member='RHL2'
        endif
        memberlongname='rhombohedral'
    elseif ( lo_two_equal(a0,b0,c0,tol_dist) .and. lo_two_equal(al0,be0,gm0,tol_angle) ) then
        ! 17 MCLC1 C-centered monoclinic 0.743303 0.743303 1.200000 ang 82.617700 82.617700 84.547378
        family='MCLC'
        ! this one is annoying to differentiate
        call get_conventional_abc_and_standard_basis(latticevectors,family,pri,con)
        call lo_get_axis_angles(con,a,b,c,al,be,gm)
        ! get the reciprocal basis as well
        m0=lo_reciprocal_basis(pri)
        call lo_get_axis_angles(m0,ka,kb,kc,kal,kbe,kgm)
        ! sort stuff out
        if ( abs(kgm-d90) .lt. tol_angle ) then
            member='MCLC2'
        elseif ( kgm .gt. d90 ) then
            member='MCLC1'
        else
            f0=b*cos(al)/c+(b*sin(al)/a)**2
            if ( abs(f0-1.0_flyt) .lt. tol_dist ) then
                member='MCLC4'
            elseif ( f0 .lt. 1.0_flyt ) then
                member='MCLC3'
            else
                member='MCLC5'
            endif
        endif
        memberlongname='c-centered monoclinic'
    else
        ! final annoying one, differentiate ORCF and TRI.
        call get_c_lt_b_lt_a_order(latticevectors,m0)
        call lo_get_axis_angles(m0,ap,bp,cp,al,be,gm)
        ! get the ORCF conventional lattice parameters
        a=2*Sqrt(bp*cp*Cos(al))
        b=2*Sqrt(ap*cp*Cos(be))
        c=2*Sqrt(ap*bp*Cos(gm))
        ! build the primitive again
        m0(:,1)=(/0.0_flyt, b/2, c/2/)
        m0(:,2)=(/a/2, 0.0_flyt, c/2/)
        m0(:,3)=(/a/2, b/2, 0.0_flyt/)
        ! this should not match for a trigonal lattice
        if ( abs(ap-norm2(m0(:,1))) .lt. tol_dist .and. &
             abs(bp-norm2(m0(:,2))) .lt. tol_dist .and. &
             abs(cp-norm2(m0(:,3))) .lt. tol_dist ) then
            ! 8 ORCF1 Face-centered orthorhombic 1.250000 1.118034 0.901388 ang 75.636697 60.050918 44.312385
            family='ORCF'
            ! check stuff
            call get_conventional_abc_and_standard_basis(latticevectors,family,pri,con)
            call lo_get_axis_angles(con,a,b,c,al,be,gm)
            f0=1.0_flyt/(a**2)
            f1=1.0_flyt/(b**2)+1.0_flyt/(c**2)
            if ( abs(f0-f1) .lt. tol_dist ) then
                member='ORCF3'
            elseif ( f0 .gt. f1 ) then
                member='ORCF1'
            else
                member='ORCF2'
            endif
            memberlongname='face-centered orthorhombic'
        else
            ! 22 TRI1a Triclinic 1.000000 1.100000 1.200000 ang 60.000000 70.000000 75.000000
            family='TRI'
            call get_conventional_abc_and_standard_basis(latticevectors,family,pri,con)
            call lo_get_axis_angles(con,a,b,c,al,be,gm)
            m0=lo_reciprocal_basis(con)
            call lo_get_axis_angles(m0,ka,kb,kc,kal,kbe,kgm)
            if ( lo_three_smaller_than_num(kal,kbe,kgm,d90,tol_angle) ) then
                member='TRI1b'
            elseif ( lo_three_larger_than_num(kal,kbe,kgm,d90,tol_angle) ) then
                member='TRI1a'
            elseif ( lo_two_smaller_than_num(kal,kbe,kgm,d90,tol_angle) ) then
                member='TRI2b'
            else
                member='TRI2a'
            endif
            memberlongname='triclinic'
        endif
    endif
end subroutine

!> Given the Bravais family, I can figure out the lattice parameter of the conventional cell, as well as a standardized primitive basis.
subroutine get_conventional_abc_and_standard_basis(latticevectors,family,pri,con,alat)
    !> the lattice vectors
    real(flyt), dimension(3,3), intent(in) :: latticevectors
    !> Bravais family name
    character(len=*), intent(in) :: family
    !> standardized primitive basis
    real(flyt), dimension(3,3), intent(out) :: pri
    !> standardized conventional basis
    real(flyt), dimension(3,3), intent(out) :: con
    !> the "a" lattice parameter, in some standard sense
    real(flyt), intent(out), optional :: alat
    !
    real(flyt), dimension(3,3) :: m0
    real(flyt) :: a,b,c,ap,bp,cp,ra,rb,rg

    select case(trim(family))
    case('CUB')
    ! CUB  a=b=c, al=be=gm=90
        call get_a_lt_b_lt_c_order(latticevectors,m0)
        call lo_get_axis_angles(m0,ap,bp,cp,ra,rb,rg)
        a = ap
        b = bp
        c = cp
        pri(:,1)=[a,0.0_flyt,0.0_flyt]
        pri(:,2)=[0.0_flyt,a,0.0_flyt]
        pri(:,3)=[0.0_flyt,0.0_flyt,a]
        con=pri
    case('FCC')
    ! FCC  a=b=c, al=be=gm
        call get_a_lt_b_lt_c_order(latticevectors,m0)
        call lo_get_axis_angles(m0,ap,bp,cp,ra,rb,rg)
        a = sqrt(2.0_flyt)*ap
        b = sqrt(2.0_flyt)*bp
        c = sqrt(2.0_flyt)*cp
        pri(:,1)=[0.0_flyt,a/2,a/2]
        pri(:,2)=[a/2,0.0_flyt,a/2]
        pri(:,3)=[a/2,a/2,0.0_flyt]

        con(:,1)=[a,0.0_flyt,0.0_flyt]
        con(:,2)=[0.0_flyt,a,0.0_flyt]
        con(:,3)=[0.0_flyt,0.0_flyt,a]
    case('BCC')
    ! BCC  a=b=c, al=be=gm
        call get_a_lt_b_lt_c_order(latticevectors,m0)
        call lo_get_axis_angles(m0,ap,bp,cp,ra,rb,rg)
        a = (2*ap)/sqrt(3.0_flyt)
        b = (2*bp)/sqrt(3.0_flyt)
        c = (2*cp)/sqrt(3.0_flyt)
        pri(:,1)=[-a/2,a/2,a/2]
        pri(:,2)=[a/2,-a/2,a/2]
        pri(:,3)=[a/2,a/2,-a/2]

        con(:,1)=[a,0.0_flyt,0.0_flyt]
        con(:,2)=[0.0_flyt,a,0.0_flyt]
        con(:,3)=[0.0_flyt,0.0_flyt,a]
    case('TET')
    ! TET  a=b, al=be=gm=90
        call get_a_eq_b_order(latticevectors,m0)
        call lo_get_axis_angles(m0,ap,bp,cp,ra,rb,rg)
        a = ap
        b = bp
        c = cp
        pri(:,1)=[a,0.0_flyt,0.0_flyt]
        pri(:,2)=[0.0_flyt,a,0.0_flyt]
        pri(:,3)=[0.0_flyt,0.0_flyt,c]
        con=pri
    case('BCT')
    ! BCT  a=b=c, al=be
        call get_al_eq_be_order(latticevectors,m0)
        call lo_get_axis_angles(m0,ap,bp,cp,ra,rb,rg)
        a = sqrt(2.0_flyt)*ap*Sqrt(-Cos(ra) - Cos(rg))
        b = sqrt(2.0_flyt)*bp*Sqrt(-Cos(ra) - Cos(rg))
        c = 2*cp*Sqrt(-Cos(rb))
        pri(:,1)=[-a/2,a/2,c/2]
        pri(:,2)=[a/2,-a/2,c/2]
        pri(:,3)=[a/2,a/2,-c/2]

        con(:,1)=[a,0.0_flyt,0.0_flyt]
        con(:,2)=[0.0_flyt,a,0.0_flyt]
        con(:,3)=[0.0_flyt,0.0_flyt,c]
    case('ORC')
    ! ORC  a<b<c al=be=gm=90
        call get_a_lt_b_lt_c_order(latticevectors,m0)
        call lo_get_axis_angles(m0,ap,bp,cp,ra,rb,rg)
        a = ap
        b = bp
        c = cp
        pri(:,1)=[a,0.0_flyt,0.0_flyt]
        pri(:,2)=[0.0_flyt,b,0.0_flyt]
        pri(:,3)=[0.0_flyt,0.0_flyt,c]
        con=pri
    case('ORCF')
    ! ORCF  c<b<a al<be<gm
        call get_c_lt_b_lt_a_order(latticevectors,m0)
        call lo_get_axis_angles(m0,ap,bp,cp,ra,rb,rg)
        a = 2*Sqrt(bp*cp*Cos(ra))
        b = 2*Sqrt(ap*cp*Cos(rb))
        c = 2*Sqrt(ap*bp*Cos(rg))
        pri(:,1)=[0.0_flyt,b/2,c/2]
        pri(:,2)=[a/2,0.0_flyt,c/2]
        pri(:,3)=[a/2,b/2,0.0_flyt]
        !
        con(:,1)=[a,0.0_flyt,0.0_flyt]
        con(:,2)=[0.0_flyt,b,0.0_flyt]
        con(:,3)=[0.0_flyt,0.0_flyt,c]
    case('ORCI')
    ! ORCI  a=b=c gm<be<al
        call get_gm_lt_be_lt_al_order(latticevectors,m0)
        call lo_get_axis_angles(m0,ap,bp,cp,ra,rb,rg)
        a = sqrt(2.0_flyt)*Sqrt(ap**2*(-Cos(rb) - Cos(rg)))
        b = sqrt(2.0_flyt)*Sqrt(ap**2*(-Cos(ra) - Cos(rg)))
        c = sqrt(2.0_flyt)*Sqrt(ap**2*(-Cos(ra) - Cos(rb)))
        pri(:,1)=[-a/2,b/2,c/2]
        pri(:,2)=[a/2,-b/2,c/2]
        pri(:,3)=[a/2,b/2,-c/2]

        con(:,1)=[a,0.0_flyt,0.0_flyt]
        con(:,2)=[0.0_flyt,b,0.0_flyt]
        con(:,3)=[0.0_flyt,0.0_flyt,c]
    case('ORCC')
        call get_al_eq_be_order(latticevectors,m0)
        call lo_get_axis_angles(m0,ap,bp,cp,ra,rb,rg)
        a = sqrt(2.0_flyt)*Sqrt(dot_product(m0(1,:),m0(1,:)) + dot_product(m0(1,:),m0(2,:)))
        b = sqrt(2.0_flyt)*Sqrt(-dot_product(m0(1,:),m0(2,:)) + dot_product(m0(2,:),m0(2,:)))
        c = cp

        pri(:,1)=[a/2,b/2,0.0_flyt]
        pri(:,2)=[a/2,-b/2,0.0_flyt]
        pri(:,3)=[0.0_flyt,0.0_flyt,c]

        con(:,1)=[a,0.0_flyt,0.0_flyt]
        con(:,2)=[0.0_flyt,b,0.0_flyt]
        con(:,3)=[0.0_flyt,0.0_flyt,c]
    case('HEX')
    ! HEX  a=b al=be
        call get_a_eq_b_order(latticevectors,m0)
        call lo_get_axis_angles(m0,ap,bp,cp,ra,rb,rg)
        a = ap
        b = bp
        c = cp
        pri(:,1)=[ a/2,-(sqrt(3.0_flyt)*a)/2 ,0.0_flyt]
        pri(:,2)=[ a/2, (sqrt(3.0_flyt)*a)/2 ,0.0_flyt]
        pri(:,3)=[0.0_flyt,0.0_flyt,c]
        con=pri
    case('HEX2')
    ! HEX  a=b al=be
        call get_a_eq_b_order(latticevectors,m0)
        call lo_get_axis_angles(m0,ap,bp,cp,ra,rb,rg)
        a = ap
        b = bp
        c = cp
        pri(:,1)=[ a, 0.0_flyt ,0.0_flyt]
        pri(:,2)=[ a/2, (sqrt(3.0_flyt)*a)/2 ,0.0_flyt]
        pri(:,3)=[0.0_flyt,0.0_flyt,c]
        con=pri
    case('RHL')
    ! RHL  a=b=c al=be=gm
        call get_a_lt_b_lt_c_order(latticevectors,m0)
        call lo_get_axis_angles(m0,ap,bp,cp,ra,rb,rg)
        a = ap
        b = bp
        c = cp
        pri(:,1)=[a*Cos(ra/2), -(a*Sin(ra/2)), 0.0_flyt ]
        pri(:,2)=[a*Cos(ra/2),  a*Sin(ra/2),   0.0_flyt ]
        pri(:,3)=[a*Cos(ra)*(1/Cos(ra/2)), 0.0_flyt, a*Sqrt(1 - Cos(ra)**2*(1/Cos(ra/2))**2) ]
        con=pri
    case('MCL')
    ! MCL  a<b<c be=gm
        call get_mcl_order(latticevectors,m0)
        call lo_get_axis_angles(m0,ap,bp,cp,ra,rb,rg)
        a = ap
        b = bp
        c = cp
        pri(:,1)=[a, 0.0_flyt, 0.0_flyt]
        pri(:,2)=[0.0_flyt, b, 0.0_flyt]
        pri(:,3)=[0.0_flyt, c*Cos(ra), c*Sin(ra)]
        con=pri
    case('MCLC')
    ! MCLC  a=b al=be
        call get_a_eq_b_order(latticevectors,m0)
        call lo_get_axis_angles(m0,ap,bp,cp,ra,rb,rg)
        a = norm2(m0(:,1)-m0(:,2))
        b = norm2(m0(:,1)+m0(:,2))
        c = cp
        ra = acos(dot_product(m0(:,1)+m0(:,2),m0(:,3))/(b*c))
        pri(:,1)=[a/2, b/2, 0.0_flyt]
        pri(:,2)=[-a/2, b/2, 0.0_flyt]
        pri(:,3)=[0.0_flyt, c*Cos(ra), c*Sin(ra)]
        !
        con(:,1)=[a, 0.0_flyt, 0.0_flyt]
        con(:,2)=[0.0_flyt, b, 0.0_flyt]
        con(:,3)=[0.0_flyt, c*Cos(ra), c*Sin(ra)]
    case('TRI')
    ! TRI  a<b<c, al/=be/=gm
        call get_a_lt_b_lt_c_order(latticevectors,m0)
        call lo_get_axis_angles(m0,ap,bp,cp,ra,rb,rg)
        a = ap
        b = bp
        c = cp

        pri(:,1)=[a, 0.0_flyt, 0.0_flyt]
        pri(:,2)=[b*Cos(rg), b*Sin(rg), 0.0_flyt]
        pri(:,3)=[c*Cos(rb), c*(Cos(ra) - Cos(rb)*Cos(rg))*(1/Sin(rg)),&
                 c*(1/Sin(rg))*Sqrt(-Cos(ra)**2 - Cos(rb)**2 + 2*Cos(ra)*Cos(rb)*Cos(rg) + Sin(rg)**2)]
        con=pri
    end select
    ! keep the "a" lattice parameter
    if ( present(alat) ) then
        alat=a
    endif
end subroutine

!> Sort the lattice vectors so that a<b<c
subroutine get_a_lt_b_lt_c_order(original_basis,reordered_basis)
    !> the input basis
    real(flyt), dimension(3,3), intent(in) :: original_basis
    !> the sorted basis
    real(flyt), dimension(3,3), intent(out) :: reordered_basis
    !
    real(flyt) :: a,b,c,al,be,gm
    real(flyt), dimension(3) :: abc
    integer, dimension(3) :: ind
    !
    call lo_get_axis_angles(original_basis(:,:),a,b,c,al,be,gm)
    abc(1)=a
    abc(2)=b
    abc(3)=c
    call qsort(abc,ind)
    reordered_basis=original_basis(:,ind)
end subroutine

!> sort the lattice vectors in MCL order
subroutine get_mcl_order(original_basis,reordered_basis)
    !> the input basis
    real(flyt), dimension(3,3), intent(in) :: original_basis
    !> the sorted basis
    real(flyt), dimension(3,3), intent(out) :: reordered_basis
    !
    real(flyt) :: a,b,c,al,be,gm,f0,f1
    integer, dimension(:,:), allocatable :: perms
    integer :: i
    !
    reordered_basis=0.0_flyt
    call lo_permutations(perms,3)
    f0=lo_huge
    do i=1,size(perms,2)
        call lo_get_axis_angles(original_basis(:,perms(:,i)),a,b,c,al,be,gm)
        f1=abs(a-c)+abs(b-c)+abs(be-0.5_flyt*lo_pi)+abs(gm-0.5_flyt*lo_pi)
        if ( f1 .lt. f0 ) then
            reordered_basis=original_basis(:,perms(:,i))
            f0=f1
        endif
    enddo
end subroutine

!> sort the lattice vectors so that c<b<a
subroutine get_c_lt_b_lt_a_order(original_basis,reordered_basis)
    !> the input basis
    real(flyt), dimension(3,3), intent(in) :: original_basis
    !> the sorted basis
    real(flyt), dimension(3,3), intent(out) :: reordered_basis
    !
    real(flyt) :: a,b,c,al,be,gm
    real(flyt), dimension(3) :: abc
    integer, dimension(3) :: ind
    !
    call lo_get_axis_angles(original_basis,a,b,c,al,be,gm)
    abc(1)=a
    abc(2)=b
    abc(3)=c
    call qsort(abc,ind)
    reordered_basis(:,1)=original_basis(:,ind(3))
    reordered_basis(:,2)=original_basis(:,ind(2))
    reordered_basis(:,3)=original_basis(:,ind(1))
end subroutine

!> sort the lattice vectors so that gamma<beta<allpha
subroutine get_gm_lt_be_lt_al_order(original_basis,reordered_basis)
    !> the input basis
    real(flyt), dimension(3,3), intent(in) :: original_basis
    !> the sorted basis
    real(flyt), dimension(3,3), intent(out) :: reordered_basis
    !
    real(flyt) :: a,b,c,al,be,gm
    real(flyt), dimension(3) :: abc
    integer, dimension(3) :: ind
    !
    call lo_get_axis_angles(original_basis(:,:),a,b,c,al,be,gm)
    abc(1)=al
    abc(2)=be
    abc(3)=gm
    call qsort(abc,ind)
    reordered_basis(:,1)=original_basis(:,ind(3))
    reordered_basis(:,2)=original_basis(:,ind(2))
    reordered_basis(:,3)=original_basis(:,ind(1))
end subroutine

!> sort the lattice vectors so that a=b
subroutine get_a_eq_b_order(original_basis,reordered_basis)
    !> the input basis
    real(flyt), dimension(3,3), intent(in) :: original_basis
    !> the sorted basis
    real(flyt), dimension(3,3), intent(out) :: reordered_basis
    !
    real(flyt) :: a,b,c,al,be,gm
    integer, dimension(:,:), allocatable :: perms
    integer :: i
    !
    reordered_basis=0.0_flyt
    call lo_permutations(perms,3)
    do i=1,size(perms,2)
        call lo_get_axis_angles(original_basis(:,perms(:,i)),a,b,c,al,be,gm)
        if ( abs(a-b) .lt. lo_tol ) then
            reordered_basis=original_basis(:,perms(:,i))
            exit
        endif
    enddo
end subroutine

!> sort the lattice vectors so that alpha=beta
subroutine get_al_eq_be_order(original_basis,reordered_basis)
    !> the input basis
    real(flyt), dimension(3,3), intent(in) :: original_basis
    !> the sorted basis
    real(flyt), dimension(3,3), intent(out) :: reordered_basis
    !
    real(flyt) :: a,b,c,al,be,gm
    integer, dimension(:,:), allocatable :: perms
    integer :: i
    !
    reordered_basis=0.0_flyt
    call lo_permutations(perms,3)
    do i=1,size(perms,2)
        call lo_get_axis_angles(original_basis(:,perms(:,i)),a,b,c,al,be,gm)
        if ( abs(al-be) .lt. lo_radiantol ) then
            reordered_basis=original_basis(:,perms(:,i))
            exit
        endif
    enddo
end subroutine

!> Tedious routine to label the high symmetry points of the BZ. There is no logic to it whatsoever. I want to retire it and just use Miller indices instead, at least they are logical.
subroutine label_highsymmetry_points(p,timereversal)
    !> the crystal structure
    type(lo_crystalstructure), intent(inout) :: p
    !> consider time reversal symmetry in the labelling?
    logical, intent(in) :: timereversal

    type(lo_symset) :: sym
    real(flyt), dimension(:,:), allocatable :: pts
    real(flyt), dimension(:,:), allocatable :: r,rw
    real(flyt), dimension(3,3) :: bas,ibas
    real(flyt), dimension(3) :: v0
    real(flyt) :: zero,one,half,fourth,fiveovereight,threeoverfour,threeovereight,twooverthree,oneoverthree
    real(flyt) :: m,g,s,f,d,x,w,l,q
    real(flyt) :: a,b,c,al,be,gm
    integer :: i,j,o,np,nhit,nmiss
    character(len=10), dimension(:), allocatable :: lbl

    if ( p%info%verbosity .gt. 0 ) then
        write(*,*) ''
        write(*,*) 'Trying to label high-symmetry points'
    endif

    ! First block is to fetch the tabulated values
    zero=0.0_flyt
    one=1.0_flyt
    half=0.5_flyt
    fourth=0.25_flyt
    fiveovereight=5.0_flyt/8.0_flyt
    threeoverfour=3.0_flyt/4.0_flyt
    threeovereight=3.0_flyt/8.0_flyt
    twooverthree=2.0_flyt/3.0_flyt
    oneoverthree=1.0_flyt/3.0_flyt

    ! Structural parameters from the conventional cell in the correct order
    call lo_get_axis_angles(p%info%unique_conventional_basis,a,b,c,al,be,gm)
    select case(trim(p%info%bravaislattice))
        case('CUB ')
            ! Table 2 3 CUB
            lo_allocate(pts(3,3))
            lo_allocate(lbl(3))
            pts(:, 1)=(/   half,   half,   half/)
            pts(:, 2)=(/   half,   half,   zero/)
            pts(:, 3)=(/   zero,   half,   zero/)
            lbl( 1)='R'
            lbl( 2)='M'
            lbl( 3)='X'
        case('FCC ')
            ! Table 3 5 FCC
            lo_allocate(pts(3,5))
            lo_allocate(lbl(5))
            pts(:, 1)=(/fiveovereight, fourth,fiveovereight/)
            pts(:, 2)=(/threeovereight,threeovereight,threeoverfour/)
            pts(:, 3)=(/   half, fourth,threeoverfour/)
            pts(:, 4)=(/   half,   half,   half/)
            pts(:, 5)=(/   half,   zero,   half/)
            lbl( 1)='U'
            lbl( 2)='K'
            lbl( 3)='W'
            lbl( 4)='L'
            lbl( 5)='X'
        case('BCC ')
        ! Table 4 3 BCC
            lo_allocate(pts(3,3))
            lo_allocate(lbl(3))
            pts(:, 1)=(/ fourth, fourth, fourth/)
            pts(:, 2)=(/   half,  -half,   half/)
            pts(:, 3)=(/   zero,   zero,   half/)
            lbl( 1)='P'
            lbl( 2)='H'
            lbl( 3)='N'
        case('TET ')
        ! Table 5 5 TET
            lo_allocate(pts(3,5))
            lo_allocate(lbl(5))
            pts(:, 1)=(/   zero,   half,   half/)
            pts(:, 2)=(/   half,   half,   half/)
            pts(:, 3)=(/   zero,   half,   zero/)
            pts(:, 4)=(/   half,   half,   zero/)
            pts(:, 5)=(/   zero,   zero,   half/)
            lbl( 1)='R'
            lbl( 2)='A'
            lbl( 3)='X'
            lbl( 4)='M'
            lbl( 5)='Z'
        case('BCT1 ')
        ! Table 6 6 BCT1
            lo_allocate(pts(3,6))
            lo_allocate(lbl(6))
            !#g=(1+c^2/a^2)/4
            g=(one+(c/a)**2)*0.25_flyt
            pts(:, 1)=(/   zero,   zero,   half/)
            pts(:, 2)=(/  -half,   half,   half/)
            pts(:, 3)=(/      g,      g,     -g/)
            pts(:, 4)=(/   zero,   half,   zero/)
            pts(:, 5)=(/     -g,  one-g,      g/)
            pts(:, 6)=(/ fourth, fourth, fourth/)
            lbl( 1)='X'
            lbl( 2)='M'
            lbl( 3)='Z'
            lbl( 4)='N'
            lbl( 5)='Z1'
            lbl( 6)='P'
        case('BCT2 ')
        ! Table 7 8 BCT2
            lo_allocate(pts(3,8))
            lo_allocate(lbl(8))
            !#g=(1+a^2/c^2)/4
            !#f=a^2/(2*c^2)
            g=(one+(a/c)**2)*0.25_flyt
            f=(a**2)/(2.0_flyt*c**2)
            pts(:, 1)=(/   zero,   zero,   half/)
            pts(:, 2)=(/   zero,   half,   zero/)
            pts(:, 3)=(/     -f,      f,   half/)
            pts(:, 4)=(/ fourth, fourth, fourth/)
            pts(:, 5)=(/   half,   half,      f/)
            pts(:, 6)=(/     -g,      g,      g/)
            pts(:, 7)=(/   half,   half,  -half/)
            pts(:, 8)=(/      g,  one-g,      g/)
            lbl( 1)='X'
            lbl( 2)='N'
            lbl( 3)='Y'
            lbl( 4)='P'
            lbl( 5)='Y1'
            lbl( 6)='R'
            lbl( 7)='Z'
            lbl( 8)='R1'
        case('ORC ')
        ! Table 8 7 ORC
            lo_allocate(pts(3,7))
            lo_allocate(lbl(7))
            pts(:, 1)=(/   half,   zero,   half/)
            pts(:, 2)=(/   half,   half,   half/)
            pts(:, 3)=(/   half,   zero,   zero/)
            pts(:, 4)=(/   half,   half,   zero/)
            pts(:, 5)=(/   zero,   half,   zero/)
            pts(:, 6)=(/   zero,   half,   half/)
            pts(:, 7)=(/   zero,   zero,   half/)
            lbl( 1)='U'
            lbl( 2)='R'
            lbl( 3)='X'
            lbl( 4)='S'
            lbl( 5)='Y'
            lbl( 6)='T'
            lbl( 7)='Z'
        case('ORCF1')
        ! Table 9 8 ORCF1 ORCF3
            lo_allocate(pts(3,8))
            lo_allocate(lbl(8))
            !#f=(1+a^2/b^2-a^2/c^2)/4
            !#g=(1+a^2/b^2+a^2/c^2)/4
            f=(one+a**2/b**2-a**2/c**2)*0.25_flyt
            g=(one+a**2/b**2+a**2/c**2)*0.25_flyt
            pts(:, 1)=(/   zero,      g,      g/)
            pts(:, 2)=(/   half, half+f,      f/)
            pts(:, 3)=(/    one,  one-g,  one-g/)
            pts(:, 4)=(/   half, half-f,  one-f/)
            pts(:, 5)=(/   half,   zero,   half/)
            pts(:, 6)=(/   half,   half,   half/)
            pts(:, 7)=(/   half,   half,   zero/)
            pts(:, 8)=(/    one,   half,   half/)
            lbl( 1)='X'
            lbl( 2)='A'
            lbl( 3)='X1'
            lbl( 4)='A1'
            lbl( 5)='Y'
            lbl( 6)='L'
            lbl( 7)='Z'
            lbl( 8)='T'
        case('ORCF3')
        ! Table 9 8 ORCF1 ORCF3
            lo_allocate(pts(3,8))
            lo_allocate(lbl(8))
            !#f=(1+a^2/b^2-a^2/c^2)/4
            !#g=(1+a^2/b^2+a^2/c^2)/4
            f=(one+a**2/b**2-a**2/c**2)*0.25_flyt
            g=(one+a**2/b**2+a**2/c**2)*0.25_flyt
            pts(:, 1)=(/   zero,      g,      g/)
            pts(:, 2)=(/   half, half+f,      f/)
            pts(:, 3)=(/    one,  one-g,  one-g/)
            pts(:, 4)=(/   half, half-f,  one-f/)
            pts(:, 5)=(/   half,   zero,   half/)
            pts(:, 6)=(/   half,   half,   half/)
            pts(:, 7)=(/   half,   half,   zero/)
            pts(:, 8)=(/    one,   half,   half/)
            lbl( 1)='X'
            lbl( 2)='A'
            lbl( 3)='X1'
            lbl( 4)='A1'
            lbl( 5)='Y'
            lbl( 6)='L'
            lbl( 7)='Z'
            lbl( 8)='T'
        case('ORCF2 ')
        ! Table 10 10 ORCF2
            lo_allocate(pts(3,10))
            lo_allocate(lbl(10))
            !#g=(1+a^2/b^2-a^2/c^2)/4
            !#d=(1+b^2/a^2-b^2/c^2)/4
            !#s=(1+c^2/b^2-c^2/a^2)/4
            g=(one+a**2/b**2-a**2/c**2)*0.25_flyt
            d=(one+b**2/a**2-b**2/c**2)*0.25_flyt
            s=(one+c**2/b**2-c**2/a**2)*0.25_flyt
            pts(:, 1)=(/  one-s, half-s,   half/)
            pts(:, 2)=(/   half, half-g,  one-g/)
            pts(:, 3)=(/      s, half+s,   half/)
            pts(:, 4)=(/   half, half+g,      g/)
            pts(:, 5)=(/   zero,   half,   half/)
            pts(:, 6)=(/ half-d,   half,  one-d/)
            pts(:, 7)=(/   half,   zero,   half/)
            pts(:, 8)=(/ half+d,   half,      d/)
            pts(:, 9)=(/   half,   half,   zero/)
            pts(:,10)=(/   half,   half,   half/)
            lbl( 1)='H'
            lbl( 2)='C'
            lbl( 3)='H1'
            lbl( 4)='C1'
            lbl( 5)='X'
            lbl( 6)='D'
            lbl( 7)='Y'
            lbl( 8)='D1'
            lbl( 9)='Z'
            lbl(10)='L'
        case('ORCI ')
        ! Table 11 12 ORCI
            lo_allocate(pts(3,12))
            lo_allocate(lbl(12))
            !#f=(1+a^2/c^2)/4
            !#g=(1+b^2/c^2)/4
            !#d=(b^2-a^2)/(4*c^2)
            !#l=(a^2+b^2)/(4*c^2)
            f=(one+a**2/c**2)*0.25_flyt
            g=(one+b**2/c**2)*0.25_flyt
            d=(b**2-a**2)/(4*c**2)
            l=(a**2+b**2)/(4*c**2)
            pts(:, 1)=(/ fourth, fourth, fourth/)
            pts(:, 2)=(/     -l,      l, half-d/)
            pts(:, 3)=(/     -f,      f,      f/)
            pts(:, 4)=(/      l,     -l, half+d/)
            pts(:, 5)=(/      f,  one-f,      f/)
            pts(:, 6)=(/ half-d, half+d,     -l/)
            pts(:, 7)=(/      g,     -g,      g/)
            pts(:, 8)=(/   zero,   half,   zero/)
            pts(:, 9)=(/  one-g,      g,     -g/)
            pts(:,10)=(/   half,   zero,   zero/)
            pts(:,11)=(/   half,   half,  -half/)
            pts(:,12)=(/   zero,   zero,   half/)
            lbl( 1)='W'
            lbl( 2)='L'
            lbl( 3)='X'
            lbl( 4)='L1'
            lbl( 5)='X1'
            lbl( 6)='L2'
            lbl( 7)='Y'
            lbl( 8)='R'
            lbl( 9)='Y1'
            lbl(10)='S'
            lbl(11)='Z'
            lbl(12)='T'
        case('ORCC1')
        ! Table 12 9 ORCC
            lo_allocate(pts(3,10))
            lo_allocate(lbl(10))
            f=(one+b**2/a**2)*0.25_flyt
            pts(:, 1)=(/   half,   half,   half/)
            pts(:, 2)=(/      f,  one-f,   half/)
            pts(:, 3)=(/     -f,      f,   zero/)
            pts(:, 4)=(/     -f,      f,   half/)
            pts(:, 5)=(/      f,  one-f,   zero/)
            pts(:, 6)=(/   zero,   half,   half/)
            pts(:, 7)=(/  -half,   half,   zero/)
            pts(:, 8)=(/   zero,   half,   zero/)
            pts(:, 9)=(/   zero,   zero,   half/)
            pts(:,10)=(/   half,   half,   zero/)
            lbl( 1)='T'
            lbl( 2)='A'
            lbl( 3)='X'
            lbl( 4)='A1'
            lbl( 5)='X1'
            lbl( 6)='R'
            lbl( 7)='Y'
            lbl( 8)='S'
            lbl( 9)='Z'
            lbl(10)='Y'
        case('ORCC2')
            ! Table 12 9 ORCC
            lo_allocate(pts(3,10))
            lo_allocate(lbl(10))
            f=(one+a**2/b**2)*0.25_flyt
            pts(:, 1)=(/  -half,   half,   half/)
            pts(:, 2)=(/      f,      f,   half/)
            pts(:, 3)=(/      f,      f,   zero/)
            pts(:, 4)=(/     -f,  one-f,   half/)
            pts(:, 5)=(/     -f,  one-f,   zero/)
            pts(:, 6)=(/   zero,   half,   half/)
            pts(:, 7)=(/  -half,   half,   zero/)
            pts(:, 8)=(/   zero,   half,   zero/)
            pts(:, 9)=(/   zero,   zero,   half/)
            pts(:,10)=(/  -half,   half,   zero/)
            lbl( 1)='T'
            lbl( 2)='A'
            lbl( 3)='X'
            lbl( 4)='A1'
            lbl( 5)='X1'
            lbl( 6)='R'
            lbl( 7)='Y'
            lbl( 8)='S'
            lbl( 9)='Z'
            lbl(10)='Y'
        case('HEX ')
        ! Table 13 5 HEX
            lo_allocate(pts(3,5))
            lo_allocate(lbl(5))
            pts(:, 1)=(/  oneoverthree,  oneoverthree,   zero/)
            pts(:, 2)=(/   zero,   zero,   half/)
            pts(:, 3)=(/   half,   zero,   half/)
            pts(:, 4)=(/ oneoverthree,  oneoverthree,   half/)
            pts(:, 5)=(/   half,   zero,   zero/)
            lbl( 1)='K'
            lbl( 2)='A'
            lbl( 3)='L'
            lbl( 4)='H'
            lbl( 5)='M'
        case('HEX2 ')
        ! Table 13 5 HEX
            lo_allocate(pts(3,5))
            lo_allocate(lbl(5))
            pts(:, 1)=(/   zero,   zero,   half/)
            pts(:, 2)=(/  oneoverthree,  twooverthree,   half/)
            pts(:, 3)=(/  half,  half,  half/)
            pts(:, 4)=(/   half,   half,   zero/)
            pts(:, 5)=(/  oneoverthree,   twooverthree,   zero/)

            lbl( 1)='A'
            lbl( 2)='L'
            lbl( 3)='H'
            lbl( 4)='K'
            lbl( 5)='M'
        case('RHL1 ')
        ! Table 14 11 RHL1
            lo_allocate(pts(3,10))
            lo_allocate(lbl(10))
            !#g=(1+4*cos(alpha))/(2+4*cos(alpha))
            !#m=3/4-g/2
            g=(1.0_flyt+4.0_flyt*cos(al))/(2.0_flyt+4.0_flyt*cos(al))
            m=3.0_flyt/4.0_flyt-g/2.0_flyt
            pts(:, 1)=(/      g,      m,      m/)
            pts(:, 2)=(/      g,   half,  one-g/)
            pts(:, 3)=(/  one-m,  one-m,  one-g/)
            pts(:, 4)=(/   half,  one-g,  g-one/)
            pts(:, 5)=(/      m,      m,  g-one/)
            pts(:, 6)=(/   half,   half,   zero/)
            pts(:, 7)=(/  one-m,      m,   zero/)
            pts(:, 8)=(/   half,   zero,   zero/)
            pts(:, 9)=(/      m,   zero,     -m/)
            pts(:,10)=(/   half,   half,   half/)
            !pts(:,11)=(/   zero,   zero,  -half/)
            lbl( 1)='P'
            lbl( 2)='B'
            lbl( 3)='P1'
            lbl( 4)='B1'
            lbl( 5)='P2'
            lbl( 6)='F'
            lbl( 7)='Q'
            lbl( 8)='L'
            lbl( 9)='X'
            lbl(10)='Z'
            !lbl(11)='L1'
        case('RHL2 ')
        ! Table 15 7 RHL2
            lo_allocate(pts(3,7))
            lo_allocate(lbl(7))
            !#g=1/( 2*tan(alpha/2)^2 )
            !#m=3/4-g/2
            g=1.0_flyt/( 2*tan(al/2)**2 )
            m=3.0_flyt/4.0_flyt-g/2.0_flyt
            pts(:, 1)=(/      m,  m-one,  m-one/)
            pts(:, 2)=(/   half,  -half,   zero/)
            pts(:, 3)=(/      g,      g,      g/)
            pts(:, 4)=(/   half,   zero,   zero/)
            pts(:, 5)=(/  one-g,     -g,     -g/)
            pts(:, 6)=(/  one-m,     -m,  one-m/)
            pts(:, 7)=(/   half,  -half,   half/)
            lbl( 1)='P1'
            lbl( 2)='F'
            lbl( 3)='Q'
            lbl( 4)='L'
            lbl( 5)='Q1'
            lbl( 6)='P'
            lbl( 7)='Z'
        case('MCL ')
        ! Table 16 15 MCL
            lo_allocate(pts(3,15))
            lo_allocate(lbl(15))
            !#g=(1-b*cos(alpha)/c)/(2*sin(alpha)^2)
            !#m=1/2-g*c*cos(alpha)/b
            ! I think the curtarolo paper is a little wrong, had to modify this
            if ( al .lt. lo_pi*0.5_flyt ) then
                g=(1-b*cos(al)/c)/(2*sin(al)**2)
                m=0.5_flyt-g*c*cos(al)/b
               ! write(*,*) ' eta',g
               ! write(*,*) '  nu',m
            else
                !g=(1+c*cos(alpha)/b)/(2*sin(alpha)^2)
                !m=1/2+g*b*cos(alpha)/c
                m=(1.0_flyt+c*cos(al)/b)/(2*sin(al)*sin(al))
                g=0.5_flyt+m*b*cos(al)/c
                m=1.0_flyt-m
                !g=1.0_flyt-g
                !write(*,*) ' eta',g
                !write(*,*) '  nu',m
            endif
            pts(:, 1)=(/   zero,      g,     -m/)
            pts(:, 2)=(/   half,   half,   zero/)
            pts(:, 3)=(/   half,      g,  one-m/)
            pts(:, 4)=(/   zero,   half,   half/)
            pts(:, 5)=(/   half,  one-g,      m/)
            pts(:, 6)=(/   half,   zero,   half/)
            pts(:, 7)=(/   half,      g,     -m/)
            pts(:, 8)=(/   half,   zero,  -half/)
            pts(:, 9)=(/   zero,   half,   zero/)
            pts(:,10)=(/   half,   half,   half/)
            pts(:,11)=(/   zero,   zero,   half/)
            pts(:,12)=(/   zero,      g,  one-m/)
            pts(:,13)=(/   zero,   zero,  -half/)
            pts(:,14)=(/   zero,  one-g,      m/)
            pts(:,15)=(/   half,   zero,   zero/)
            lbl( 1)='H2'
            lbl( 2)='A'
            lbl( 3)='M'
            lbl( 4)='C'
            lbl( 5)='M1'
            lbl( 6)='D'
            lbl( 7)='M2'
            lbl( 8)='D1'
            lbl( 9)='X'
            lbl(10)='E'
            lbl(11)='Y'
            lbl(12)='H'
            lbl(13)='Y1:q'
            lbl(14)='H1'
            lbl(15)='Z'
        case('MCLC1')
        ! Table 17 16 MCLC1  MCLC2
            lo_allocate(pts(3,16))
            lo_allocate(lbl(16))
            !#f=(2-b*cos(alpha)/c)/(4*sin(alpha)^2)
            !#g=1/2+2*f*c*cos(alpha)/b
            !#w=3/4-(a/(2*b*sin(alpha)))^2
            !#s=w+(3/4-w)*b*cos(alpha)/c
            f=(2.0_flyt-b*cos(al)/c)/(4*sin(al)**2)
            g=0.5_flyt+2*f*c*cos(al)/b
            w=0.75_flyt-(a/(2*b*sin(al)))**2
            s=w+(0.75_flyt-w)*b*cos(al)/c
            pts(:, 1)=(/   half,   half,   half/)
            pts(:, 2)=(/   half,   zero,   zero/)
            pts(:, 3)=(/   half,   zero,   half/)
            pts(:, 4)=(/   zero,  -half,   zero/)
            pts(:, 5)=(/  one-w,  w-one,   zero/)
            pts(:, 6)=(/  one-f,  one-f,  one-g/)
            pts(:, 7)=(/      w,  one-w,   zero/)
            pts(:, 8)=(/      f,      f,      g/)
            pts(:, 9)=(/  w-one,     -w,   zero/)
            pts(:,10)=(/     -f,     -f,  one-g/)
            pts(:,11)=(/   half,   half,   zero/)
            pts(:,12)=(/  one-f,     -f,  one-g/)
            pts(:,13)=(/  -half,  -half,   zero/)
            pts(:,14)=(/      s,  one-s,   half/)
            pts(:,15)=(/   zero,   zero,   half/)
            pts(:,16)=(/  one-s,  s-one,   half/)
            lbl( 1)='L'
            lbl( 2)='N'
            lbl( 3)='M'
            lbl( 4)='N1'
            lbl( 5)='X'
            lbl( 6)='F'
            lbl( 7)='X1'
            lbl( 8)='F1'
            lbl( 9)='X2'
            lbl(10)='F2'
            lbl(11)='Y'
            lbl(12)='F3'
            lbl(13)='Y1'
            lbl(14)='I'
            lbl(15)='Z'
            lbl(16)='I1'
        case('MCLC2')
        ! Table 17 16 MCLC1  MCLC2
            lo_allocate(pts(3,16))
            lo_allocate(lbl(16))
            !#f=(2-b*cos(alpha)/c)/(4*sin(alpha)^2)
            !#g=1/2+2*f*c*cos(alpha)/b
            !#w=3/4-(a/(2*b*sin(alpha)))^2
            !#s=w+(3/4-w)*b*cos(alpha)/c
            f=(2.0_flyt-b*cos(al)/c)/(4*sin(al)**2)
            g=0.5_flyt+2*f*c*cos(al)/b
            w=0.75_flyt-(a/(2*b*sin(al)))**2
            s=w+(0.75_flyt-w)*b*cos(al)/c
            pts(:, 1)=(/   half,   half,   half/)
            pts(:, 2)=(/   half,   zero,   zero/)
            pts(:, 3)=(/   half,   zero,   half/)
            pts(:, 4)=(/   zero,  -half,   zero/)
            pts(:, 5)=(/  one-w,  w-one,   zero/)
            pts(:, 6)=(/  one-f,  one-f,  one-g/)
            pts(:, 7)=(/      w,  one-w,   zero/)
            pts(:, 8)=(/      f,      f,      g/)
            pts(:, 9)=(/  w-one,     -w,   zero/)
            pts(:,10)=(/     -f,     -f,  one-g/)
            pts(:,11)=(/   half,   half,   zero/)
            pts(:,12)=(/  one-f,     -f,  one-g/)
            pts(:,13)=(/  -half,  -half,   zero/)
            pts(:,14)=(/      s,  one-s,   half/)
            pts(:,15)=(/   zero,   zero,   half/)
            pts(:,16)=(/  one-s,  s-one,   half/)
            lbl( 1)='L'
            lbl( 2)='N'
            lbl( 3)='M'
            lbl( 4)='N1'
            lbl( 5)='X'
            lbl( 6)='F'
            lbl( 7)='X1'
            lbl( 8)='F1'
            lbl( 9)='X2'
            lbl(10)='F2'
            lbl(11)='Y'
            lbl(12)='F3'
            lbl(13)='Y1'
            lbl(14)='I'
            lbl(15)='Z'
            lbl(16)='I1'
        case('MCLC3')
        ! Table 18 16 MCLC3 MCLC4
            lo_allocate(pts(3,16))
            lo_allocate(lbl(16))
            !#l  =  (1 + b^2/a^2)/4
            !#d  =  b*c*cos(alpha)/(2*a^2)
            !#f  =  l - 1/4 + (1-b*cos(alpha)/c)/(4*sin(alpha)^2)
            !#g  =  1/2 + 2*f*c*cos(alpha)/b
            !#s  =  1 + f - 2*l
            !#w  =  g - 2d
            l  =  (1.0_flyt + b**2/a**2)*0.25_flyt
            d  =  b*c*cos(al)/(2*a**2)
            f  =  l - 0.25_flyt + (1.0_flyt-b*cos(al)/c)/(4*sin(al)**2)
            g  =  0.5_flyt + 2*f*c*cos(al)/b
            s  =  1.0_flyt + f - 2*l
            w  =  g - 2*d
            pts(:, 1)=(/   half,   zero,   zero/)
            pts(:, 2)=(/  one-s,  one-s,  one-w/)
            pts(:, 3)=(/   zero,   half,   zero/)
            pts(:, 4)=(/      s,      s,  one-w/)
            pts(:, 5)=(/   half,   half,   zero/)
            pts(:, 6)=(/  one-s,      s,  one-w/)
            pts(:, 7)=(/      l,      l,      d/)
            pts(:, 8)=(/      f,      f,      g/)
            pts(:, 9)=(/  one-l,     -l,      d/)
            pts(:,10)=(/  one-f,     -f,  one-g/)
            pts(:,11)=(/      l,      l,      d/)
            pts(:,12)=(/      f,      f,  one-g/)
            pts(:,13)=(/      l,  l-one,      d/)
            pts(:,14)=(/   half,   half,   half/)
            pts(:,15)=(/   zero,   zero,   half/)
            pts(:,16)=(/   half,   zero,   half/)
            lbl( 1)='N'
            lbl( 2)='F'
            lbl( 3)='N1'
            lbl( 4)='F1'
            lbl( 5)='X'
            lbl( 6)='F2'
            lbl( 7)='Y'
            lbl( 8)='H'
            lbl( 9)='Y1'
            lbl(10)='H1'
            lbl(11)='Y2'
            lbl(12)='H2'
            lbl(13)='Y3'
            lbl(14)='I'
            lbl(15)='Z'
            lbl(16)='M'
        case('MCLC4')
        ! Table 18 16 MCLC3 MCLC4
            lo_allocate(pts(3,16))
            lo_allocate(lbl(16))
            !#l  =  (1 + b^2/a^2)/4
            !#d  =  b*c*cos(alpha)/(2*a^2)
            !#f  =  l - 1/4 + (1-b*cos(alpha)/c)/(4*sin(alpha)^2)
            !#g  =  1/2 + 2*f*c*cos(alpha)/b
            !#s  =  1 + f - 2*l
            !#w  =  g - 2d
            l  =  (1.0_flyt + b**2/a**2)*0.25_flyt
            d  =  b*c*cos(al)/(2*a**2)
            f  =  l - 0.25_flyt + (1.0_flyt-b*cos(al)/c)/(4*sin(al)**2)
            g  =  0.5_flyt + 2*f*c*cos(al)/b
            s  =  1.0_flyt + f - 2*l
            w  =  g - 2*d
            pts(:, 1)=(/   half,   zero,   zero/)
            pts(:, 2)=(/  one-s,  one-s,  one-w/)
            pts(:, 3)=(/   zero,   half,   zero/)
            pts(:, 4)=(/      s,      s,  one-w/)
            pts(:, 5)=(/   half,   half,   zero/)
            pts(:, 6)=(/  one-s,      s,  one-w/)
            pts(:, 7)=(/      l,      l,      d/)
            pts(:, 8)=(/      f,      f,      g/)
            pts(:, 9)=(/  one-l,     -l,      d/)
            pts(:,10)=(/  one-f,     -f,  one-g/)
            pts(:,11)=(/      l,      l,      d/)
            pts(:,12)=(/      f,      f,  one-g/)
            pts(:,13)=(/      l,  l-one,      d/)
            pts(:,14)=(/   half,   half,   half/)
            pts(:,15)=(/   zero,   zero,   half/)
            pts(:,16)=(/   half,   zero,   half/)
            lbl( 1)='N'
            lbl( 2)='F'
            lbl( 3)='N1'
            lbl( 4)='F1'
            lbl( 5)='X'
            lbl( 6)='F2'
            lbl( 7)='Y'
            lbl( 8)='H'
            lbl( 9)='Y1'
            lbl(10)='H1'
            lbl(11)='Y2'
            lbl(12)='H2'
            lbl(13)='Y3'
            lbl(14)='I'
            lbl(15)='Z'
            lbl(16)='M'
        case('MCLC5 ')
        ! Table 19 16 MCLC5
            lo_allocate(pts(3,16))
            lo_allocate(lbl(16))
            !#f  = ( (b^2/a^2)  + (1-b*cos(alpha)/c)/(sin(alpha)^2 )/4
            !#l  =  g/2 + (b/(2*a))^2 - b*c*cos(alpha)/(2*a^2)
            !#x  =  (4*m - 1 - (b*sin(alpha)/a)^2) * c/(2*b*cos(alpha))
            !#g  =  1/2 + 2*f*c*cos(a)/b
            !#m  =  2*l - f
            !#d  =  f*c*cos(alpha)/b + x/2 - 1/4
            !#q  =  1 - f*(a/b)^2
            f  = ( (b**2/a**2)  + (1.0_flyt-b*cos(al)/c)/(sin(al)**2) )*0.25_flyt
            g  =  0.5_flyt + 2*f*c*cos(al)/b
            l  =  g*0.5_flyt + (b/(2*a))**2 - b*c*cos(al)/(2*a**2)
            m  =  2*l - f
            x  =  (4*m - 1.0_flyt - (b*sin(al)/a)**2)*c/(2*b*cos(al))
            d  =  f*c*cos(al)/b + x/2 - 0.25_flyt
            q  =  1.0_flyt - f*(a/b)**2
            !
            pts(:, 1)=(/   half,   zero,   half/)
            pts(:, 2)=(/      m,      m,      x/)
            pts(:, 3)=(/   half,   zero,   zero/)
            pts(:, 4)=(/  one-m,  one-m,  one-x/)
            pts(:, 5)=(/   zero,  -half,   zero/)
            pts(:, 6)=(/      m,  m-one,      x/)
            pts(:, 7)=(/   half,  -half,   zero/)
            pts(:, 8)=(/      f,      f,      g/)
            pts(:, 9)=(/      l,      l,      d/)
            pts(:,10)=(/  one-f,     -f,  one-g/)
            pts(:,11)=(/  one-l,     -l,     -d/)
            pts(:,12)=(/     -f,     -f,  one-g/)
            pts(:,13)=(/     -l,     -l,     -d/)
            pts(:,14)=(/      q,  one-q,   half/)
            pts(:,15)=(/      l,  l-one,      d/)
            pts(:,16)=(/  one-q,  q-one,   half/)
            lbl( 1)='M'
            lbl( 2)='F'
            lbl( 3)='N'
            lbl( 4)='F1'
            lbl( 5)='N1'
            lbl( 6)='F2'
            lbl( 7)='X'
            lbl( 8)='H'
            lbl( 9)='Y'
            lbl(10)='H1'
            lbl(11)='Y1'
            lbl(12)='H2'
            lbl(13)='Y2'
            lbl(14)='I'
            lbl(15)='Y3'
            lbl(16)='I1'
        case('TRI1a')
        ! Table 20 7 TRI1a TRI2a
            lo_allocate(pts(3,7))
            lo_allocate(lbl(7))
            pts(:, 1)=(/   half,   half,   zero/)
            pts(:, 2)=(/   zero,   half,   half/)
            pts(:, 3)=(/   half,   zero,   half/)
            pts(:, 4)=(/   half,   half,   half/)
            pts(:, 5)=(/   half,   zero,   zero/)
            pts(:, 6)=(/   zero,   half,   zero/)
            pts(:, 7)=(/   zero,   zero,   half/)
            lbl( 1)='L'
            lbl( 2)='M'
            lbl( 3)='N'
            lbl( 4)='R'
            lbl( 5)='X'
            lbl( 6)='Y'
            lbl( 7)='Z'
        case('TRI2a')
        ! Table 20 7 TRI1a TRI2a
            lo_allocate(pts(3,7))
            lo_allocate(lbl(7))
            pts(:, 1)=(/   half,   half,   zero/)
            pts(:, 2)=(/   zero,   half,   half/)
            pts(:, 3)=(/   half,   zero,   half/)
            pts(:, 4)=(/   half,   half,   half/)
            pts(:, 5)=(/   half,   zero,   zero/)
            pts(:, 6)=(/   zero,   half,   zero/)
            pts(:, 7)=(/   zero,   zero,   half/)
            lbl( 1)='L'
            lbl( 2)='M'
            lbl( 3)='N'
            lbl( 4)='R'
            lbl( 5)='X'
            lbl( 6)='Y'
            lbl( 7)='Z'
        case('TRI1b')
        ! Table 21 7 TRI1b TRI2b
            lo_allocate(pts(3,7))
            lo_allocate(lbl(7))
            pts(:, 1)=(/   half,  -half,   zero/)
            pts(:, 2)=(/   zero,   zero,   half/)
            pts(:, 3)=(/  -half,  -half,   half/)
            pts(:, 4)=(/   zero,  -half,   half/)
            pts(:, 5)=(/   zero,  -half,   zero/)
            pts(:, 6)=(/   half,   zero,   zero/)
            pts(:, 7)=(/  -half,   zero,   half/)
            lbl( 1)='L'
            lbl( 2)='M'
            lbl( 3)='N'
            lbl( 4)='R'
            lbl( 5)='X'
            lbl( 6)='Y'
            lbl( 7)='Z'
        case('TRI2b')
        ! Table 21 7 TRI1b TRI2b
            lo_allocate(pts(3,7))
            lo_allocate(lbl(7))
            pts(:, 1)=(/   half,  -half,   zero/)
            pts(:, 2)=(/   zero,   zero,   half/)
            pts(:, 3)=(/  -half,  -half,   half/)
            pts(:, 4)=(/   zero,  -half,   half/)
            pts(:, 5)=(/   zero,  -half,   zero/)
            pts(:, 6)=(/   half,   zero,   zero/)
            pts(:, 7)=(/  -half,   zero,   half/)
            lbl( 1)='L'
            lbl( 2)='M'
            lbl( 3)='N'
            lbl( 4)='R'
            lbl( 5)='X'
            lbl( 6)='Y'
            lbl( 7)='Z'
        case default
            write(*,*) 'Unknown Bravais lattice type'
            lo_allocate(pts(3,1))
            lo_allocate(lbl(1))
            pts=0.0_flyt
            lbl(1)='GM'
    end select

    ! First try to label the full BZ, set the labels to nothing:
    if ( allocated(p%bz%label) ) deallocate(p%bz%label)
    lo_allocate(p%bz%label( p%bz%nhighsymmetrypoints ))
    do i=1,p%bz%nhighsymmetrypoints
        p%bz%label(i)='NP'
    enddo

    ! get the tabulated version of the reciprocal lattice
    bas=lo_reciprocal_basis(p%info%unique_primitive_basis)
    ibas=lo_invert3x3matrix( bas )

    ! fetch the coordinates for the bz nodes
    lo_allocate(r(3,p%bz%nhighsymmetrypoints))
    do i=1,p%bz%nhighsymmetrypoints
        ! to fractional
        r(:,i)=p%cartesian_to_fractional( p%bz%highsymmetrypoints(:,i), reciprocal=.true., pbc=.false. )
        ! then strange permutation
        r(:,i)=matmul(r(:,i),p%info%permutation_to_unique)
        ! to Cartesian again, in the new coordinate system
        r(:,i)=matmul(bas,r(:,i))
    enddo

    ! Also get the tabulated points to cartesian
    np=size(pts,2)
    do i=1,np
        pts(:,i)=matmul(bas,pts(:,i))
    enddo
    ! Get all possible symmetry operations in this represenation
    call sym%generate(bas,timereversal)
    if ( p%info%verbosity .gt. 0 ) write(*,*) '... generated symmorphic group'

    ! Try to label things:
    nhit=1
    nmiss=0
    ! minus 1 because I know that the last point is gamma
    bzptloop: do i=1,p%bz%nhighsymmetrypoints-1
        do o=1,sym%n
            ! rotate the BZ points
            v0=lo_operate_on_vector(sym%op(o),r(:,i),reciprocal=.true.)
            do j=1,np
                if ( lo_sqnorm(v0-pts(:,j)) .lt. lo_sqtol ) then
                    ! I found a point, I guess
                    nhit=nhit+1
                    p%bz%label(i)=trim(lbl(j))
                    cycle bzptloop
                endif
                if ( sym%timereversal ) then
                    if ( lo_sqnorm(v0+pts(:,j)) .lt. lo_sqtol ) then
                        ! I found a point, I guess
                        nhit=nhit+1
                        p%bz%label(i)=trim(lbl(j))
                        cycle bzptloop
                    endif
                endif
            enddo
        enddo
        ! if I made it here, I guess I did not find a label.
        nmiss=nmiss+1
    enddo bzptloop
    p%bz%label(p%bz%nhighsymmetrypoints)='GM'

    if ( p%info%verbosity .gt. 0 ) then
        !write(*,*) '... tried labelling: '//trim(int2char(nhit))//' hits, '//trim(int2char(nmiss))//' misses'
        write(*,*) '... tried labelling: ',nhit,' hits, ',nmiss,' misses'
        if ( nmiss .eq. 0 ) then
            write(*,*) '... highly successful, all points found'
        elseif ( nhit .eq. 0 ) then
            write(*,*) '... did not go very well, no points found'
        else
            write(*,*) '... moderate success, some points found'
        endif
    endif

    ! Test backward, see if all points got assigned
    nhit=0
    nmiss=0
    bwloop: do i=1,np
        do j=1,p%bz%nhighsymmetrypoints
            if ( trim(lbl(i)) .eq. trim(p%bz%label(j)) ) then
                nhit=nhit+1
                cycle bwloop
            endif
        enddo
        nmiss=nmiss+1
    enddo bwloop

    if ( p%info%verbosity .gt. 0 ) then
        write(*,*) '... testing assignment: '//tochar(nhit)//' hits, '//tochar(nmiss)//' misses'
        if ( nmiss .eq. 0 ) then
            write(*,*) '... highly successful, all points got assigned'
        elseif ( nhit .eq. 0 ) then
            write(*,*) '... did not go very well, nothing assigned'
        else
            write(*,*) '... moderate success, some assigned'
        endif
    endif

    ! Now that the full BZ is somewhat labelled, try to label the wedge as well.
    ! fetch the coordinates for the wedge
    if ( allocated(p%irrw%label) ) deallocate(p%irrw%label)
    lo_allocate(p%irrw%label(p%irrw%nnodes))
    lo_allocate(rw(3,p%irrw%nnodes))
    do i=1,p%irrw%nnodes
        p%irrw%label(i)='NP'
    enddo
    do i=1,p%irrw%nnodes
        ! to fractional
        rw(:,i)=p%cartesian_to_fractional( p%irrw%r(:,i), reciprocal=.true., pbc=.false. )
        ! then strange permutation
        rw(:,i)=matmul(rw(:,i),p%info%permutation_to_unique)
        ! to Cartesian again, in the new coordinate system
        rw(:,i)=matmul(bas,rw(:,i))
    enddo

    ! Figure out labels for the wedge nodes
    nhit=0
    nmiss=0
    ! first try without any operations, might have success with that.
    wedgeptloop1: do i=1,p%irrw%nnodes
        do j=1,p%bz%nhighsymmetrypoints
            if ( lo_sqnorm(rw(:,i)-r(:,j)) .lt. lo_sqtol ) then
                ! I found a point, if it is labelled
                if ( trim(p%bz%label(j)) .ne. 'NP' ) then
                    nhit=nhit+1
                    p%irrw%label(i)=trim(p%bz%label(j))
                    cycle wedgeptloop1
                endif
            endif
        enddo
        ! if I made it here, I guess I did not find a label.
        nmiss=nmiss+1
    enddo wedgeptloop1

    ! If not successful, try again with all operations. It might work.
    if ( nmiss .ne. 0 ) then
        !
        if ( p%info%verbosity .gt. 0 ) then
            write(*,*) '... not complete success first pass, trying again'
        endif
        ! New attempt, with all operations
        wedgeptloop2: do i=1,p%irrw%nnodes
            ! skip if already fixed
            if ( trim(p%irrw%label(i)) .ne. 'NP' ) cycle
            ! test vs all BZ points
            do j=1,p%bz%nhighsymmetrypoints
                ! skip if this one also failed
                if ( trim(p%bz%label(j)) .eq. 'NP' ) cycle
                ! test vs all operations
                do o=1,sym%n
                    ! rotate the wedge point
                    v0=lo_operate_on_vector(sym%op(o),rw(:,i),reciprocal=.true.)
                    if ( lo_sqnorm(v0-r(:,j)) .lt. lo_sqtol ) then
                        ! I found a point, I guess
                        nhit=nhit+1
                        p%irrw%label(i)=trim(p%bz%label(j))
                        cycle wedgeptloop2
                    endif
                    if ( sym%timereversal ) then
                        if ( lo_sqnorm(v0+r(:,j)) .lt. lo_sqtol ) then
                            nhit=nhit+1
                            p%irrw%label(i)=trim(p%bz%label(j))
                            cycle wedgeptloop2
                        endif
                    endif
                enddo
            enddo
            ! if I made it here, I guess I did not find a label.
            nmiss=nmiss+1
        enddo wedgeptloop2
    endif

    ! Report on success
    if ( p%info%verbosity .gt. 0 ) then
        write(*,*) '... labelled wedge: '//tochar(nhit)//' hits, '//tochar(nmiss)//' misses'
    endif

    ! Name the new points in the wedge, maybe
    if ( nmiss .ne. 0 ) then
        ! name the points in the wedge that don't have real names
        j=0
        do i=1,p%irrw%nnodes
            if ( trim(p%irrw%label(i)) .eq. 'NP' ) then
                j=j+1
                p%irrw%label(i)='NP'//tochar(j)
            endif
        enddo

        ! Now go over the wedge and spread these names:
        bzmissloop: do i=1,p%bz%nhighsymmetrypoints
            if ( trim(p%bz%label(i)) .ne. 'NP' ) cycle
            do j=1,p%irrw%nnodes
                do o=1,sym%n
                    v0=lo_operate_on_vector(sym%op(o),rw(:,j),reciprocal=.true.)
                    if ( lo_sqnorm(v0-r(:,i)) .lt. lo_sqtol ) then
                        p%bz%label(i)=trim(p%irrw%label(j))
                        cycle bzmissloop
                    endif
                    if ( sym%timereversal ) then
                        if ( lo_sqnorm(v0+r(:,i)) .lt. lo_sqtol ) then
                            p%bz%label(i)=trim(p%irrw%label(j))
                            cycle bzmissloop
                        endif
                    endif
                enddo
            enddo
            ! Hopefully, I will never make it here
            write(*,*) 'Completely failed labelling point ',tochar(i)
            stop
        enddo bzmissloop
    endif
    p%info%pointslabelled=.true.
    if ( p%info%verbosity .gt. 0 ) write(*,*) 'Did my best at labelling points'
    !
end subroutine


!> Check if three numbers are equal within some tolerance. If no tolerance is specified, I guess one.
logical pure function lo_three_equal(a,b,c,tolerance)
    !> first number
    real(flyt), intent(in) :: a
    !> second number
    real(flyt), intent(in) :: b
    !> third number
    real(flyt), intent(in) :: c
    !> the tolerance
    real(flyt), intent(in), optional :: tolerance
    !
    integer :: i
    real(flyt) :: thres
    ! maybe guess a tolerance
    if ( present(tolerance) ) then
        thres=tolerance
    else
        thres=(abs(a)+abs(b)+abs(c))/1E8_flyt
    endif
    ! check
    i=0
    if ( abs(a-b) .lt. thres ) i=i+1
    if ( abs(a-c) .lt. thres ) i=i+1
    if ( abs(b-c) .lt. thres ) i=i+1
    !
    if ( i .eq. 3 ) then
        lo_three_equal=.true.
    else
        lo_three_equal=.false.
    endif
end function

!> Check if exactly one number is equal to some number within some tolerance. If no tolerance is specified, I guess one.
logical pure function lo_one_equal_to_num(a,b,c,num,tolerance)
    !> first number
    real(flyt), intent(in) :: a
    !> second number
    real(flyt), intent(in) :: b
    !> third number
    real(flyt), intent(in) :: c
    !> reference number
    real(flyt), intent(in) :: num
    !> the tolerance
    real(flyt), intent(in), optional :: tolerance
    !
    integer :: i
    real(flyt) :: thres
    ! maybe guess a tolerance
    if ( present(tolerance) ) then
        thres=tolerance
    else
        thres=(abs(a)+abs(b)+abs(c))/1E8_flyt
    endif
    ! check
    i=0
    if ( abs(a-num) .lt. thres ) i=i+1
    if ( abs(b-num) .lt. thres ) i=i+1
    if ( abs(c-num) .lt. thres ) i=i+1
    !
    if ( i .eq. 1 ) then
        lo_one_equal_to_num=.true.
    else
        lo_one_equal_to_num=.false.
    endif
end function

!> Check if exactly two numbers are equal to some number within some tolerance. If no tolerance is specified, I guess one.
logical pure function lo_two_equal_to_num(a,b,c,num,tolerance)
    !> first number
    real(flyt), intent(in) :: a
    !> second number
    real(flyt), intent(in) :: b
    !> third number
    real(flyt), intent(in) :: c
    !> reference number
    real(flyt), intent(in) :: num
    !> the tolerance
    real(flyt), intent(in), optional :: tolerance
    !
    integer :: i
    real(flyt) :: thres
    ! maybe guess a tolerance
    if ( present(tolerance) ) then
        thres=tolerance
    else
        thres=(abs(a)+abs(b)+abs(c))/1E8_flyt
    endif
    ! check
    i=0
    if ( abs(a-num) .lt. thres ) i=i+1
    if ( abs(b-num) .lt. thres ) i=i+1
    if ( abs(c-num) .lt. thres ) i=i+1
    !
    if ( i .eq. 2 ) then
        lo_two_equal_to_num=.true.
    else
        lo_two_equal_to_num=.false.
    endif
end function

!> Check if three numbers are equal to some number within some tolerance. If no tolerance is specified, I guess one.
logical pure function lo_three_equal_to_num(a,b,c,num,tolerance)
    !> first number
    real(flyt), intent(in) :: a
    !> second number
    real(flyt), intent(in) :: b
    !> third number
    real(flyt), intent(in) :: c
    !> reference number
    real(flyt), intent(in) :: num
    !> the tolerance
    real(flyt), intent(in), optional :: tolerance
    !
    integer :: i
    real(flyt) :: thres
    ! maybe guess a tolerance
    if ( present(tolerance) ) then
        thres=tolerance
    else
        thres=(abs(a)+abs(b)+abs(c))/1E8_flyt
    endif
    ! check
    i=0
    if ( abs(a-num) .lt. thres ) i=i+1
    if ( abs(b-num) .lt. thres ) i=i+1
    if ( abs(c-num) .lt. thres ) i=i+1
    !
    if ( i .eq. 3 ) then
        lo_three_equal_to_num=.true.
    else
        lo_three_equal_to_num=.false.
    endif
end function

!> Check if three numbers are less than some number within some tolerance. If no tolerance is specified, I guess one.
logical pure function lo_three_smaller_than_num(a,b,c,num,tolerance)
    !> first number
    real(flyt), intent(in) :: a
    !> second number
    real(flyt), intent(in) :: b
    !> third number
    real(flyt), intent(in) :: c
    !> reference number
    real(flyt), intent(in) :: num
    !> the tolerance
    real(flyt), intent(in), optional :: tolerance
    !
    integer :: i
    real(flyt) :: thres
    ! maybe guess a tolerance
    if ( present(tolerance) ) then
        thres=tolerance
    else
        thres=(abs(a)+abs(b)+abs(c))/1E8_flyt
    endif
    ! check how many equal
    i=0
    if ( a .lt. num-thres ) i=i+1
    if ( b .lt. num-thres ) i=i+1
    if ( c .lt. num-thres ) i=i+1
    !
    if ( i .eq. 3 ) then
        lo_three_smaller_than_num=.true.
    else
        lo_three_smaller_than_num=.false.
    endif
end function

!> Check if three numbers are greater than some number within some tolerance. If no tolerance is specified, I guess one.
logical pure function lo_three_larger_than_num(a,b,c,num,tolerance)
    !> first number
    real(flyt), intent(in) :: a
    !> second number
    real(flyt), intent(in) :: b
    !> third number
    real(flyt), intent(in) :: c
    !> reference number
    real(flyt), intent(in) :: num
    !> the tolerance
    real(flyt), intent(in), optional :: tolerance
    !
    integer :: i
    real(flyt) :: thres
    ! maybe guess a tolerance
    if ( present(tolerance) ) then
        thres=tolerance
    else
        thres=(abs(a)+abs(b)+abs(c))/1E8_flyt
    endif
    ! check how many equal
    i=0
    if ( a .gt. num-thres ) i=i+1
    if ( b .gt. num-thres ) i=i+1
    if ( c .gt. num-thres ) i=i+1
    !
    if ( i .eq. 3 ) then
        lo_three_larger_than_num=.true.
    else
        lo_three_larger_than_num=.false.
    endif
end function

!> Check if three numbers are less than some number within some tolerance. If no tolerance is specified, I guess one.
logical pure function lo_two_smaller_than_num(a,b,c,num,tolerance)
    !> first number
    real(flyt), intent(in) :: a
    !> second number
    real(flyt), intent(in) :: b
    !> third number
    real(flyt), intent(in) :: c
    !> reference number
    real(flyt), intent(in) :: num
    !> the tolerance
    real(flyt), intent(in), optional :: tolerance
    !
    integer :: i
    real(flyt) :: thres
    ! maybe guess a tolerance
    if ( present(tolerance) ) then
        thres=tolerance
    else
        thres=(abs(a)+abs(b)+abs(c))/1E8_flyt
    endif
    ! check how many equal
    i=0
    if ( a .lt. num-thres ) i=i+1
    if ( b .lt. num-thres ) i=i+1
    if ( c .lt. num-thres ) i=i+1
    !
    if ( i .eq. 2 ) then
        lo_two_smaller_than_num=.true.
    else
        lo_two_smaller_than_num=.false.
    endif
end function

!> Check if two out of three numbers are equal within some tolerance. If no tolerance is specified, I guess one. Will return false for three equal numbers.
logical pure function lo_two_equal(a,b,c,tolerance)
    !> first number
    real(flyt), intent(in) :: a
    !> second number
    real(flyt), intent(in) :: b
    !> third number
    real(flyt), intent(in) :: c
    !> the tolerance
    real(flyt), intent(in), optional :: tolerance
    !
    integer :: i
    real(flyt) :: thres
    ! maybe guess a tolerance
    if ( present(tolerance) ) then
        thres=tolerance
    else
        thres=(abs(a)+abs(b)+abs(c))/1E8_flyt
    endif
    !
    i=0
    if ( abs(a-b) .lt. thres ) i=i+1
    if ( abs(a-c) .lt. thres ) i=i+1
    if ( abs(b-c) .lt. thres ) i=i+1
    !
    if ( i .eq. 1 ) then
        lo_two_equal=.true.
    else
        lo_two_equal=.false.
    endif
end function

!> Check that none of three numbers are equal within some tolerance. If no tolerance is specified, I guess one.
logical pure function lo_none_equal(a,b,c,tolerance)
    !> first number
    real(flyt), intent(in) :: a
    !> second number
    real(flyt), intent(in) :: b
    !> third number
    real(flyt), intent(in) :: c
    !> the tolerance
    real(flyt), intent(in), optional :: tolerance
    !
    integer :: i
    real(flyt) :: thres
    ! maybe guess a tolerance
    if ( present(tolerance) ) then
        thres=tolerance
    else
        thres=(abs(a)+abs(b)+abs(c))/1E8_flyt
    endif
    !
    i=0
    if ( abs(a-b) .lt. thres ) i=i+1
    if ( abs(a-c) .lt. thres ) i=i+1
    if ( abs(b-c) .lt. thres ) i=i+1
    !
    if ( i .eq. 0 ) then
        lo_none_equal=.true.
    else
        lo_none_equal=.false.
    endif
end function

end submodule
