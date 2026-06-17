!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2026 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module moddump
!
! Remove particles outside a cylinder
!
! :References: None
!
! :Owner: Daniel Price
!
! :Runtime parameters:
!   - radius : *cylinder radius [code units]*
!   - xcen   : *cylinder centre x [code units]*
!   - ycen   : *cylinder centre y [code units]*
!   - zcen   : *cylinder centre z [code units]*
!   - zmax   : *cylinder half-height [code units]*
!
! :Dependencies: infile_utils, io, part, prompting
!
 implicit none
 character(len=*), parameter, public :: moddump_flags = ''

 ! runtime parameters (written to / read from the prefix.moddump file)
 real :: xcen   = 0.      ! cylinder centre x [code units]
 real :: ycen   = 0.      ! cylinder centre y [code units]
 real :: zcen   = 0.      ! cylinder centre z [code units]
 real :: radius = 1500.   ! cylinder radius [code units]
 real :: zmax   = 1500.   ! cylinder half-height [code units]

contains

subroutine modify_dump(npart,npartoftype,massoftype,xyzh,vxyzu)
 use part,          only:delete_particles_outside_cylinder
 use io,            only:id,master,fileprefix
 use infile_utils,  only:get_options
 integer, intent(inout) :: npart
 integer, intent(inout) :: npartoftype(:)
 real,    intent(inout) :: massoftype(:)
 real,    intent(inout) :: xyzh(:,:),vxyzu(:,:)
 real :: center(3)
 integer :: ierr

 print*,' Phantommoddump: Remove particles outside a cylinder'

 call get_options(trim(fileprefix)//'.moddump',id==master,ierr,&
                  read_moddumpfile,write_moddumpfile,read_interactive_moddumpfile)
 if (ierr /= 0) stop 'rerun phantommoddump with the new .moddump file'

 center = (/xcen,ycen,zcen/)
 !
 !--removing particles
 !
 print*,'Removing particles outside the cylinder centered in ( ', center,' ), with radius ',radius,' and zmax ',zmax,' : '
 call delete_particles_outside_cylinder(center, radius, zmax, npartoftype)

end subroutine modify_dump

subroutine read_interactive_moddumpfile()
 use prompting, only:prompt

 call prompt('Enter cylinder centre x (code units)',xcen)
 call prompt('Enter cylinder centre y (code units)',ycen)
 call prompt('Enter cylinder centre z (code units)',zcen)
 call prompt('Enter cylinder radius (code units)',radius,0.)
 call prompt('Enter cylinder half-height zmax (code units)',zmax,0.)

end subroutine read_interactive_moddumpfile

subroutine read_moddumpfile(filename,ierr)
 use infile_utils, only:open_db_from_file,inopts,read_inopt,close_db
 character(len=*), intent(in)  :: filename
 integer,          intent(out) :: ierr
 integer, parameter :: iunit = 23
 type(inopts), allocatable :: db(:)
 integer :: nerr

 nerr = 0
 call open_db_from_file(db,filename,iunit,ierr)
 if (ierr /= 0) return
 call read_inopt(xcen,'xcen',db,errcount=nerr)
 call read_inopt(ycen,'ycen',db,errcount=nerr)
 call read_inopt(zcen,'zcen',db,errcount=nerr)
 call read_inopt(radius,'radius',db,errcount=nerr,min=0.)
 call read_inopt(zmax,'zmax',db,errcount=nerr,min=0.)
 call close_db(db)
 if (nerr > 0) ierr = nerr

end subroutine read_moddumpfile

subroutine write_moddumpfile(filename)
 use infile_utils, only:write_inopt,write_moddump_header
 character(len=*), intent(in) :: filename
 integer, parameter :: iunit = 23

 open(unit=iunit,file=filename,status='replace',form='formatted')
 call write_moddump_header(iunit)
 write(iunit,"(/,a)") '# remove-outside-cylinder parameters'
 call write_inopt(xcen,'xcen','cylinder centre x [code units]',iunit)
 call write_inopt(ycen,'ycen','cylinder centre y [code units]',iunit)
 call write_inopt(zcen,'zcen','cylinder centre z [code units]',iunit)
 call write_inopt(radius,'radius','cylinder radius [code units]',iunit)
 call write_inopt(zmax,'zmax','cylinder half-height [code units]',iunit)
 close(iunit)

end subroutine write_moddumpfile

end module moddump
