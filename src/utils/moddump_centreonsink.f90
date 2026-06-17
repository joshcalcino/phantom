!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2026 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module moddump
!
! Center positions and velocities of particles around one given sink
!
! :References: None
!
! :Owner: Antoine Alaguero
!
! :Runtime parameters:
!   - sink_ind : *index of the sink to centre on*
!
! :Dependencies: infile_utils, io, part, prompting
!
 implicit none
 character(len=*), parameter, public :: moddump_flags = ''

 ! runtime parameter
 integer :: sink_ind = 1   ! index of the sink to centre positions/velocities on

contains

subroutine modify_dump(npart,npartoftype,massoftype,xyzh,vxyzu)
 use part,         only:xyzmh_ptmass,vxyz_ptmass,nptmass
 use io,           only:id,master,fileprefix
 use infile_utils, only:get_options
 integer, intent(inout) :: npart
 integer, intent(inout) :: npartoftype(:)
 real,    intent(inout) :: massoftype(:)
 real,    intent(inout) :: xyzh(:,:),vxyzu(:,:)
 integer                :: i,ierr

 sink_ind = 1
 call get_options(trim(fileprefix)//'.moddump',id==master,ierr,&
                  read_moddumpfile,write_moddumpfile,read_interactive_moddumpfile)
 if (ierr /= 0) stop 'rerun phantommoddump with the new .moddump file'

 if (nptmass < sink_ind) then
    print*,'Selected sink index larger than number of sinks'
    return
 endif

 do i=1,npart
    xyzh(1:3,i) = xyzh(1:3,i) - xyzmh_ptmass(1:3,sink_ind)
    vxyzu(1:3,i) = vxyzu(1:3,i) - vxyz_ptmass(1:3,sink_ind)
 enddo

 do i=1,nptmass
 !skip sink_ind, because useless if put to 0 first
     if (i==sink_ind) then
         cycle
     else
         xyzmh_ptmass(1:3,i) = xyzmh_ptmass(1:3,i) - xyzmh_ptmass(1:3,sink_ind)
         vxyz_ptmass(1:3,i) = vxyz_ptmass(1:3,i) - vxyz_ptmass(1:3,sink_ind)
     endif
 enddo

 xyzmh_ptmass(1:3,sink_ind) = xyzmh_ptmass(1:3,sink_ind) - xyzmh_ptmass(1:3,sink_ind)
 vxyz_ptmass(1:3,sink_ind) = vxyz_ptmass(1:3,sink_ind) - vxyz_ptmass(1:3,sink_ind)

 return
end subroutine modify_dump

!----------------------------------------------------------------
!+
!  set parameters interactively (when no .moddump file is found)
!+
!----------------------------------------------------------------
subroutine read_interactive_moddumpfile()
 use prompting, only:prompt

 call prompt('Enter index of the sink to centre on',sink_ind,1)

end subroutine read_interactive_moddumpfile

!----------------------------------------------------------------
!+
!  write options to .moddump file
!+
!----------------------------------------------------------------
subroutine write_moddumpfile(filename)
 use infile_utils, only:write_inopt,write_moddump_header
 character(len=*), intent(in) :: filename
 integer, parameter :: iunit = 23

 open(unit=iunit,file=filename,status='replace',form='formatted')
 call write_moddump_header(iunit)
 write(iunit,"(/,a)") '# centre-on-sink parameters'
 call write_inopt(sink_ind,'sink_ind','index of the sink to centre on',iunit)
 close(iunit)

end subroutine write_moddumpfile

!----------------------------------------------------------------
!+
!  read options from .moddump file
!+
!----------------------------------------------------------------
subroutine read_moddumpfile(filename,ierr)
 use infile_utils, only:open_db_from_file,inopts,read_inopt,close_db
 character(len=*), intent(in)  :: filename
 integer,          intent(out) :: ierr
 integer, parameter :: iunit = 21
 type(inopts), allocatable :: db(:)
 integer :: nerr

 nerr = 0
 call open_db_from_file(db,filename,iunit,ierr)
 if (ierr /= 0) return
 call read_inopt(sink_ind,'sink_ind',db,errcount=nerr,min=1)
 call close_db(db)
 if (nerr > 0) ierr = nerr

end subroutine read_moddumpfile

end module moddump
