!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2026 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module moddump
!
! merges particles; input simulation with npart, get npart/nchild back
!
! :References: Vacondio et al. 2013
!
! :Owner: Rebecca Nealon
!
! :Runtime parameters:
!   - nchild : *number of children per merged particle (>= 2)*
!
! :Dependencies: infile_utils, io, part, prompting, splitpart
!
 implicit none
 character(len=*), parameter, public :: moddump_flags = ''

 ! runtime parameter
 integer :: nchild = 2   ! number of children merged into each particle (>= 2)

contains

subroutine modify_dump(npart,npartoftype,massoftype,xyzh,vxyzu)
 use splitpart,   only:merge_all_particles
 use part,        only:igas,kill_particle,delete_dead_or_accreted_particles
 use part,        only:isdead_or_accreted,copy_particle
 use io,          only:fatal,error,id,master,fileprefix
 use infile_utils,only:get_options
 integer, intent(inout) :: npart
 integer, intent(inout) :: npartoftype(:)
 real,    intent(inout) :: massoftype(:)
 real,    intent(inout) :: xyzh(:,:),vxyzu(:,:)
 integer :: i,nactive,ierr

 call get_options(trim(fileprefix)//'.moddump',id==master,ierr,&
                  read_moddumpfile,write_moddumpfile,read_interactive_moddumpfile)
 if (ierr /= 0) stop 'rerun phantommoddump with the new .moddump file'
 if (nchild < 2) stop 'error nchild cannot be < 2'

 !-- how many active particles
 nactive = 0
 do i = 1,npart
    if (.not.isdead_or_accreted(xyzh(4,i))) then
       nactive = nactive + 1
    else
       call kill_particle(i,npartoftype)
    endif
 enddo

 if (nactive < npart) then
    call delete_dead_or_accreted_particles(npart,npartoftype)
    print*,' discarding inactive particles'
 endif

 ! Merge 'em!
 call merge_all_particles(npart,npartoftype,massoftype,xyzh,vxyzu, &
                          nchild,nactive)

 print*,' new npart = ',npart

end subroutine modify_dump

subroutine read_interactive_moddumpfile()
 use prompting, only:prompt

 call prompt('Enter number of children per merged particle',nchild,2)

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
 call read_inopt(nchild,'nchild',db,errcount=nerr,min=2)
 call close_db(db)
 if (nerr > 0) ierr = nerr

end subroutine read_moddumpfile

subroutine write_moddumpfile(filename)
 use infile_utils, only:write_inopt,write_moddump_header
 character(len=*), intent(in) :: filename
 integer, parameter :: iunit = 23

 open(unit=iunit,file=filename,status='replace',form='formatted')
 call write_moddump_header(iunit)
 write(iunit,"(/,a)") '# merge parameters'
 call write_inopt(nchild,'nchild','number of children per merged particle (>= 2)',iunit)
 close(iunit)

end subroutine write_moddumpfile

end module moddump
