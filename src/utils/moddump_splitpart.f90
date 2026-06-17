!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2026 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module moddump
!
! split every particle in a dump into nchild children
!
! :References: None
!
! :Owner: Daniel Price
!
! :Runtime parameters:
!   - lattice_type : *child arrangement (0=regular lattice, 1=random)*
!   - nchild       : *number of children per particle (>= 2; forced to 13 for lattice)*
!
! :Dependencies: infile_utils, io, part, prompting, splitpart
!
 implicit none
 character(len=*), parameter, public :: moddump_flags = ''

 ! runtime parameters
 integer :: nchild       = 13   ! number of children per particle
 integer :: lattice_type = 0    ! 0 for regular lattice, 1 for random

contains

subroutine modify_dump(npart,npartoftype,massoftype,xyzh,vxyzu)
 use splitpart,    only:split_all_particles
 use io,           only:fatal,error,id,master,fileprefix
 use part,         only:delete_dead_or_accreted_particles
 use infile_utils, only:get_options
 integer, intent(inout) :: npart
 integer, intent(inout) :: npartoftype(:)
 real,    intent(inout) :: massoftype(:)
 real,    intent(inout) :: xyzh(:,:),vxyzu(:,:)
 integer :: ierr

 call get_options(trim(fileprefix)//'.moddump',id==master,ierr,&
                  read_moddumpfile,write_moddumpfile,read_interactive_moddumpfile)
 if (ierr /= 0) stop 'rerun phantommoddump with the new .moddump file'

 !-- the regular lattice requires a specific number of children
 if (lattice_type == 0) nchild = 13
 if (nchild < 2) stop 'error nchild cannot be < 2'

 !-- don't split accreted particles
 call delete_dead_or_accreted_particles(npart,npartoftype)

 ! Split 'em!
 print "(/,a,i0,a)", ' >>> splitting all particles into ',nchild,' children <<<'

 if (lattice_type==0) then
    print "(a,/)", ' >>> placing children on regular lattice <<<'
 else
    print "(a,/)", ' >>> placing children using random arrangement <<<'
 endif
 call split_all_particles(npart,npartoftype,massoftype,xyzh,vxyzu, &
                          nchild,lattice_type,1)

 print "(a,i0,/)",' new npart = ',npart

end subroutine modify_dump

subroutine read_interactive_moddumpfile()
 use prompting, only:prompt

 call prompt('Enter child arrangement (0=regular lattice, 1=random)',lattice_type,0,1)
 if (lattice_type /= 0) call prompt('Enter number of children per particle',nchild,2)

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
 call read_inopt(lattice_type,'lattice_type',db,errcount=nerr,min=0,max=1)
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
 write(iunit,"(/,a)") '# split parameters'
 call write_inopt(lattice_type,'lattice_type','child arrangement (0=regular lattice, 1=random)',iunit)
 call write_inopt(nchild,'nchild','number of children per particle (>= 2; forced to 13 for lattice)',iunit)
 close(iunit)

end subroutine write_moddumpfile

end module moddump
