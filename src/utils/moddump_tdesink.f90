!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2026 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module moddump
!
! Add a sink particle to the dump
!
! :References: None
!
! :Owner: Megha Sharma
!
! :Runtime parameters:
!   - sink_radius : *sink radius [Rsun]*
!
! :Dependencies: centreofmass, infile_utils, io, part, physcon, prompting,
!   sortutils, units
!
 implicit none
 character(len=*), parameter, public :: moddump_flags = ''

 real :: sink_radius = 1.0   ! sink radius in Rsun

contains

subroutine modify_dump(npart,npartoftype,massoftype,xyzh,vxyzu)
 use part, only:xyzmh_ptmass,vxyz_ptmass,nptmass,igas,ihacc,ihsoft,rhoh
 use part, only:delete_particles_inside_radius
 use io,             only:fatal,id,master,fileprefix
 use centreofmass,   only:reset_centreofmass
 use units,          only:umass,udist
 use sortutils,      only:set_r2func_origin,indexxfunc,r2func_origin
 use physcon,        only:solarr,solarm
 use infile_utils,   only:get_options

 integer, intent(inout) :: npart
 integer, intent(inout) :: npartoftype(:)
 real,    intent(inout) :: massoftype(:)
 real,    intent(inout) :: xyzh(:,:),vxyzu(:,:)

 real :: mcore,rcore,xpos(3),vpos(3)
 real :: den_all(npart),pmass,r
 integer :: j,n,location,iorder(npart),i,ierr
 !
 ! mass of gas particle
 !
 pmass = massoftype(1)
 mcore = 0.0
 !
 ! check gas particles exist
 !
 if (npart <= 0) then
    call fatal('moddump','no gas particles present in file')
 endif
 do j = 1,npart
    den_all(j) = rhoh(xyzh(4,j),pmass)
 enddo

 location = maxloc(den_all,dim=1)
 !
 ! read the sink radius (or write a template and stop)
 !
 call get_options(trim(fileprefix)//'.moddump',id==master,ierr,&
                  read_moddumpfile,write_moddumpfile,read_interactive_moddumpfile)
 if (ierr /= 0) stop 'rerun phantommoddump with the new .moddump file'

 if (sink_radius <= 0.0) then
    call fatal('moddump','Invalid sink radius entered')
 endif
 rcore = sink_radius * solarr / udist

 !
 ! sort particles by radius from dense core
 !
 xpos(:) = xyzh(1:3,location)
 vpos(:) = vxyzu(1:3,location)
 !
 ! sorting particles by radius
 !
 call set_r2func_origin(xpos(1),xpos(2),xpos(3))
 call indexxfunc(npart,r2func_origin,xyzh,iorder)
 !
 ! find particles within the rcore and determine mcore
 !
 do i = 1, npart
    r = sqrt(dot_product(xyzh(1:3,i)-xpos,xyzh(1:3,i)-xpos))
    if (r < rcore) then
       mcore = mcore + pmass
    endif
 enddo
 print*,'Mass of sink: ', mcore*umass/solarm, 'Msun , Radius of sink: ', rcore*udist/solarr, 'Rsun'
 !
 ! set sink particle
 !
 nptmass = 1
 n = nptmass
 xyzmh_ptmass(:,n)   = 0.
 xyzmh_ptmass(1:3,n) = xpos(:)
 xyzmh_ptmass(4,n)   = mcore
 xyzmh_ptmass(ihsoft,n) = rcore
 vxyz_ptmass(1:3, n) = vpos(:)
 !
 ! delete the gas particles from the region <= rcore
 !
 call delete_particles_inside_radius(xpos,rcore,npart,npartoftype)
 print*, 'Number of gas particles:',npart, ', Number of sink:',nptmass

end subroutine modify_dump

!
!---Interactively set the moddump parameters--------------------------------
!
subroutine read_interactive_moddumpfile()
 use prompting, only:prompt

 call prompt('Enter sink radius in Rsun', sink_radius, 0.)

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
 call write_inopt(sink_radius,'sink_radius','sink radius [Rsun]',iunit)
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
 call read_inopt(sink_radius,'sink_radius',db,min=0.,errcount=nerr)
 call close_db(db)
 if (nerr > 0) ierr = nerr

end subroutine read_moddumpfile

end module moddump
