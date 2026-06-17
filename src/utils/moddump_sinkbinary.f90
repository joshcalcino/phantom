!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2026 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module moddump
!
! Add a sink particle binary to the dump
!
! :References: None
!
! :Owner: Mike Lau
!
! :Runtime parameters:
!   - a    : *semi-major axis [code units]*
!   - e    : *eccentricity*
!   - m1   : *mass of point mass [Msun]*
!   - racc : *accretion radius of point mass [code units]*
!
! :Dependencies: centreofmass, extern_gwinspiral, externalforces,
!   infile_utils, io, options, part, physcon, prompting, setbinary, timestep,
!   units
!
 implicit none
 character(len=*), parameter, public :: moddump_flags = ''

 real :: m1   = 1.4     ! mass of point mass in Msun
 real :: a    = 2500.   ! semi-major axis in code units
 real :: racc = 12.5    ! accretion radius of point mass in code units
 real :: e    = 0.      ! eccentricity

contains

subroutine modify_dump(npart,npartoftype,massoftype,xyzh,vxyzu)
 use part,      only:xyzmh_ptmass,vxyz_ptmass,nptmass,igas
 use setbinary, only:set_binary
 use units,     only:umass,udist
 use physcon,   only:solarm,solarr,pi
 use io,        only:id,master,fileprefix
 use centreofmass, only:reset_centreofmass
 use timestep, only:dtmax,tmax
 use options,   only:iexternalforce
 use externalforces, only:iext_gwinspiral
 use extern_gwinspiral, only:Nstar
 use infile_utils, only:get_options
 integer, intent(inout) :: npart
 integer, intent(inout) :: npartoftype(:)
 real,    intent(inout) :: massoftype(:)
 real,    intent(inout) :: xyzh(:,:),vxyzu(:,:)
 real :: m2,period
 integer :: i,ierr

 ! find current mass from existing particles
 m2 = npart*massoftype(igas)
 print*,' Mass of star from existing file in Msun = ',m2*umass/solarm
 print*,' code unit of distance in Rsun = ',udist/solarr
 call reset_centreofmass(npart,xyzh,vxyzu)
 !
 ! read the binary parameters (or write a template and stop)
 !
 call get_options(trim(fileprefix)//'.moddump',id==master,ierr,&
                  read_moddumpfile,write_moddumpfile,read_interactive_moddumpfile)
 if (ierr /= 0) stop 'rerun phantommoddump with the new .moddump file'

 print*,' accretion radius in solar radii = ',racc*udist/solarr
 !
 ! add sink particle binary
 !
 nptmass = 0
 call set_binary(m1,m2,a,e,racc,racc,xyzmh_ptmass,vxyz_ptmass,nptmass,ierr)
 !
 ! delete one of the sink particles and replace it with our polytrope
 !
 nptmass = nptmass - 1
 do i=1,npart
    xyzh(1:3,i) = xyzh(1:3,i) + xyzmh_ptmass(1:3,2)
    vxyzu(1:3,i) = vxyzu(1:3,i) + vxyz_ptmass(1:3,2)
 enddo
 call reset_centreofmass(npart,xyzh,vxyzu,nptmass,xyzmh_ptmass,vxyz_ptmass)

 if (iexternalforce==iext_gwinspiral) then
    Nstar(1) = npart
 endif

 period = 2.*pi*sqrt(a**3/(m1 + m2))
 print*,' orbital period = ',period
 tmax = 1000.*period
 dtmax = 0.1*period

end subroutine modify_dump

!
!---Interactively set the moddump parameters--------------------------------
!
subroutine read_interactive_moddumpfile()
 use prompting, only:prompt

 call prompt('enter mass of point mass in Msun ',m1,0.)
 call prompt('enter semi-major axis in code units',a,0.)
 racc = a/200.
 call prompt('enter accretion radius of point mass',racc,0.01)
 call prompt('enter eccentricity ',e,0.)

end subroutine read_interactive_moddumpfile

!
!---Write the moddump parameter file----------------------------------------
!
subroutine write_moddumpfile(filename)
 use infile_utils, only:write_inopt,write_moddump_header
 character(len=*), intent(in) :: filename
 integer, parameter :: iunit = 20

 print "(a)",' writing moddump params file '//trim(filename)
 open(unit=iunit,file=filename,status='replace',form='formatted')
 call write_moddump_header(iunit)
 call write_inopt(m1,'m1','mass of point mass [Msun]',iunit)
 call write_inopt(a,'a','semi-major axis [code units]',iunit)
 call write_inopt(racc,'racc','accretion radius of point mass [code units]',iunit)
 call write_inopt(e,'e','eccentricity',iunit)
 close(iunit)

end subroutine write_moddumpfile

!
!---Read the moddump parameter file-----------------------------------------
!
subroutine read_moddumpfile(filename,ierr)
 use infile_utils, only:open_db_from_file,inopts,read_inopt,close_db
 character(len=*), intent(in)  :: filename
 integer,          intent(out) :: ierr
 integer, parameter :: iunit = 21
 integer :: nerr
 type(inopts), allocatable :: db(:)

 print "(a)",' reading moddump options from '//trim(filename)
 nerr = 0
 call open_db_from_file(db,filename,iunit,ierr)
 if (ierr /= 0) return
 call read_inopt(m1,'m1',db,min=0.,errcount=nerr)
 call read_inopt(a,'a',db,min=0.,errcount=nerr)
 call read_inopt(racc,'racc',db,min=0.01,errcount=nerr)
 call read_inopt(e,'e',db,min=0.,errcount=nerr)
 call close_db(db)
 if (nerr > 0) ierr = nerr

end subroutine read_moddumpfile

end module moddump
