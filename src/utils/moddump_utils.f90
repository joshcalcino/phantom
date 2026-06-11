!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2026 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module moddump_utils
!
! Shared helpers for moddumps: handles the prefix.mod parameter file in the
! same way get_options handles .setup files, and carries the flag that tells
! a moddump whether to prompt interactively
!
! :References: None
!
! :Owner: Josh Calcino
!
! :Runtime parameters: None
!
! :Dependencies: None
!
 implicit none

 ! set by the driver via get_moddump_options, read by modify_dump:
 ! .true.  -> no parameter file was found, prompt the user interactively
 ! .false. -> a complete parameter file was read, do not prompt
 logical, public :: prompt_for_params = .true.

 ! set by the driver, written into the .mod file as a record (comments)
 character(len=120), public :: moddump_dumpfile_in = ''  ! input dump being modified
 real,               public :: moddump_time = 0.         ! output time of the new dump

 public :: get_moddump_options

 abstract interface
  subroutine read_func(filename,ierr)
   character(len=*), intent(in)  :: filename
   integer,          intent(out) :: ierr
  end subroutine read_func
 end interface

 abstract interface
  subroutine write_func(filename)
   character(len=*), intent(in) :: filename
  end subroutine write_func
 end interface

 private

contains

!----------------------------------------------------------------
!+
!  high level read/write of the prefix.mod parameter file, modelled
!  on get_options for .setup files:
!   - file complete   -> read it, do not prompt (prompt_for_params=.false.)
!   - file incomplete -> add the missing options and stop so the user
!                        can edit it (no modified dump is produced)
!   - file absent     -> prompt interactively (prompt_for_params=.true.)
!+
!----------------------------------------------------------------
subroutine get_moddump_options(modfile,master,read_pars,write_pars)
 character(len=*), intent(in) :: modfile
 logical,          intent(in) :: master
 procedure(read_func)  :: read_pars
 procedure(write_func) :: write_pars
 integer :: ierr
 logical :: iexist

 inquire(file=modfile,exist=iexist)
 if (iexist) then
    call read_pars(modfile,ierr)
    if (ierr /= 0) then
       ! missing or invalid options: top up the file and stop
       if (master) then
          call write_pars(modfile)
          print "(/,2a)", ' >> added missing parameters to ', trim(modfile)
          print "(2a,/)", ' >> edit it and rerun phantommoddump'
       endif
       stop
    endif
    prompt_for_params = .false.
    if (master) print "(/,2a,/)", ' Reading moddump parameters from ', trim(modfile)
 else
    ! no parameter file: modify_dump will prompt the user interactively
    prompt_for_params = .true.
 endif

end subroutine get_moddump_options

end module moddump_utils
