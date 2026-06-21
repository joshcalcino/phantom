!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2026 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module moddump
!
! Changes particle mass
!
! :References: None
!
! :Owner: Daniel Price
!
! :Runtime parameters:
!   - disc_mass : *desired total disc mass [code units]*
!
! :Dependencies: infile_utils, io, part, prompting, units
!
 implicit none
 character(len=*), parameter, public :: moddump_flags = ''

 ! runtime parameter (written to / read from the prefix.moddump file)
 real :: disc_mass = 0.05   ! desired total disc mass [code units]

contains

subroutine modify_dump(npart,npartoftype,massoftype,xyzh,vxyzu)
 use part,          only:igas,isdead_or_accreted,kill_particle,shuffle_part
 use units,         only:umass
 use io,            only:id,master,fileprefix
 use infile_utils,  only:get_options
 integer, intent(inout) :: npart
 integer, intent(inout) :: npartoftype(:)
 real,    intent(inout) :: massoftype(:)
 real,    intent(inout) :: xyzh(:,:),vxyzu(:,:)
 real     :: current_disc_mass, mass_factor
 integer  :: i,ierr

 ! Remove particles that are dead or accreted
 do i=1,npart
    if (isdead_or_accreted(xyzh(4,i))) then
       call kill_particle(i)
    endif
 enddo

 call shuffle_part(npart)
 npartoftype(igas) = npart

 call get_options(trim(fileprefix)//'.moddump',id==master,ierr,&
                  read_moddumpfile,write_moddumpfile,read_interactive_moddumpfile)
 if (ierr /= 0) stop 'rerun phantommoddump with the new .moddump file'

 current_disc_mass = npartoftype(igas)*massoftype(igas)
 mass_factor = disc_mass/current_disc_mass

 massoftype(igas) = mass_factor*massoftype(igas)
 print*,'Particle mass is now ', massoftype(igas)*umass, ' g'
 print*,'Total disc mass is now ', npartoftype(igas)*massoftype(igas)*umass, ' g'
 print*,'Total disc mass is now ', npartoftype(igas)*massoftype(igas), 'Msun'

end subroutine modify_dump

subroutine read_interactive_moddumpfile()
 use prompting, only:prompt

 call prompt('Enter desired total disc mass in code units',disc_mass,0.)

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
 call read_inopt(disc_mass,'disc_mass',db,errcount=nerr,min=0.)
 call close_db(db)
 if (nerr > 0) ierr = nerr

end subroutine read_moddumpfile

subroutine write_moddumpfile(filename)
 use infile_utils, only:write_inopt,write_moddump_header
 character(len=*), intent(in) :: filename
 integer, parameter :: iunit = 23

 open(unit=iunit,file=filename,status='replace',form='formatted')
 call write_moddump_header(iunit)
 write(iunit,"(/,a)") '# changemass parameters'
 call write_inopt(disc_mass,'disc_mass','desired total disc mass [code units]',iunit)
 close(iunit)

end subroutine write_moddumpfile

end module moddump
