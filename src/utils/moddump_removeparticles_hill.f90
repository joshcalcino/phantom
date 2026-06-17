!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2026 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module moddump
!
! Delete particles within a fraction of the Hill radius of a sink (planet)
!
! :References: None
!
! :Owner: Josh Calcino
!
! :Runtime parameters:
!   - delete_inside : *delete particles within the Hill radius fraction*
!   - hill_fraction : *act on particles within this fraction of the Hill radius*
!   - isink         : *index of the sink treated as the planet*
!   - istar         : *index of the sink treated as the central star*
!   - zero_vz       : *set vz=0 for particles within the Hill radius fraction (kept)*
!
! :Dependencies: infile_utils, io, part
!
 implicit none
 character(len=*), parameter, public :: moddump_flags = ''

 integer :: isink         = 2       ! index of the planet sink
 integer :: istar         = 1       ! index of the central star sink
 real    :: hill_fraction = 1.0     ! fraction of the Hill radius to act within
 logical :: delete_inside = .true.  ! delete particles within the cut radius
 logical :: zero_vz       = .false. ! set vz=0 for particles within the cut radius (kept)

contains

subroutine modify_dump(npart,npartoftype,massoftype,xyzh,vxyzu)
 use part,         only:xyzmh_ptmass,nptmass,delete_particles_inside_radius
 use io,           only:id,master,fileprefix
 use infile_utils, only:get_options
 integer, intent(inout) :: npart
 integer, intent(inout) :: npartoftype(:)
 real,    intent(inout) :: massoftype(:)
 real,    intent(inout) :: xyzh(:,:),vxyzu(:,:)
 integer :: ierr,i,nzero
 real    :: mplanet,mstar,sep,r_hill,r_cut,dx(3),r

 isink         = 2
 istar         = 1
 hill_fraction = 1.0
 delete_inside = .true.
 zero_vz       = .false.
 call get_options(trim(fileprefix)//'.moddump',id==master,ierr,&
                  read_moddumpfile,write_moddumpfile,read_interactive_moddumpfile)
 if (ierr /= 0) stop 'rerun phantommoddump with the new .moddump file'

 if (nptmass < max(isink,istar)) then
    print*,'ERROR: requested sink index larger than number of sinks (nptmass=',nptmass,')'
    return
 endif

 mplanet = xyzmh_ptmass(4,isink)
 mstar   = xyzmh_ptmass(4,istar)
 if (mstar <= 0.) then
    print*,'ERROR: central star (sink ',istar,') has non-positive mass'
    return
 endif

 !--orbital separation of the planet from the star
 dx  = xyzmh_ptmass(1:3,isink) - xyzmh_ptmass(1:3,istar)
 sep = sqrt(dot_product(dx,dx))

 !--Hill radius: R_H = a (m_planet / (3 M_star))^(1/3)
 r_hill = sep*(mplanet/(3.*mstar))**(1./3.)
 r_cut  = hill_fraction*r_hill

 print "(a,i3)",   ' planet sink index       = ',isink
 print "(a,es10.3)",' planet mass             = ',mplanet
 print "(a,es10.3)",' star mass               = ',mstar
 print "(a,es10.3)",' orbital separation      = ',sep
 print "(a,es10.3)",' Hill radius             = ',r_hill
 print "(a,es10.3)",' deletion radius (R_cut) = ',r_cut

 !--optionally zero the vz of particles inside R_cut, keeping them
 !  (done before any deletion so it acts on the original particle set)
 if (zero_vz) then
    nzero = 0
    do i=1,npart
       dx = xyzh(1:3,i) - xyzmh_ptmass(1:3,isink)
       r  = sqrt(dot_product(dx,dx))
       if (r < r_cut) then
          vxyzu(3,i) = 0.
          nzero = nzero + 1
       endif
    enddo
    print "(a,i10,a)",' set vz=0 for ',nzero,' particles within R_cut of the planet'
 endif

 if (delete_inside) then
    print "(a)",' removing particles within R_cut of the planet'
    call delete_particles_inside_radius(xyzmh_ptmass(1:3,isink),r_cut,npart,npartoftype)
 endif

end subroutine modify_dump

!----------------------------------------------------------------
!+
!  set parameters interactively (when no .moddump file is found)
!+
!----------------------------------------------------------------
subroutine read_interactive_moddumpfile()
 use prompting, only:prompt

 call prompt('Enter index of the sink treated as the planet',isink,1)
 call prompt('Enter index of the sink treated as the central star',istar,1)
 call prompt('Enter fraction of the Hill radius to act within',hill_fraction,0.)
 call prompt('Delete particles within the Hill radius fraction ?',delete_inside)
 call prompt('Set vz=0 for particles within the Hill radius fraction (kept) ?',zero_vz)

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
 write(iunit,"(/,a)") '# remove-particles-within-Hill-radius parameters'
 call write_inopt(isink,'isink','index of the sink treated as the planet',iunit)
 call write_inopt(istar,'istar','index of the sink treated as the central star',iunit)
 call write_inopt(hill_fraction,'hill_fraction','act on particles within this fraction of the Hill radius',iunit)
 call write_inopt(delete_inside,'delete_inside','delete particles within the Hill radius fraction',iunit)
 call write_inopt(zero_vz,'zero_vz','set vz=0 for particles within the Hill radius fraction (kept)',iunit)
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
 call read_inopt(isink,'isink',db,errcount=nerr,min=1)
 call read_inopt(istar,'istar',db,errcount=nerr,min=1)
 call read_inopt(hill_fraction,'hill_fraction',db,errcount=nerr,min=0.)
 call read_inopt(delete_inside,'delete_inside',db,errcount=nerr)
 call read_inopt(zero_vz,'zero_vz',db,errcount=nerr)
 call close_db(db)
 if (nerr > 0) ierr = nerr

end subroutine read_moddumpfile

end module moddump
