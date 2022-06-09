!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2022 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
module moddump
!
! Changes particle mass
!
! :References: None
!
! :Owner: Daniel Price
!
! :Runtime parameters: None
!
! :Dependencies: part, units
!
 implicit none

contains

subroutine modify_dump(npart,npartoftype,massoftype,xyzh,vxyzu)
 use part, only : igas,isdead_or_accreted,kill_particle,shuffle_part
 use units, only : umass
 integer, intent(inout) :: npart
 integer, intent(inout) :: npartoftype(:)
 real,    intent(inout) :: massoftype(:)
 real,    intent(inout) :: xyzh(:,:),vxyzu(:,:)
 real     :: disc_mass, current_disc_mass, mass_factor
 integer  :: i

  ! Remove particles that are dead or accreted
  do i=1,npartoftype(igas)
    if (isdead_or_accreted(xyzh(4,i))) then
      call kill_particle(i)
    endif
  enddo

  call shuffle_part(npart)
  npartoftype(igas) = npart

 ! Specify the disc mass, or explicitly set the mass factor to scale gas mass
 disc_mass = 0.05
 current_disc_mass = npartoftype(igas)*massoftype(igas)
 mass_factor = disc_mass/current_disc_mass
 ! mass_factor = 5.

 massoftype(igas) = mass_factor*massoftype(igas)
 print*,'Particle mass is now ', massoftype(igas)*umass, ' g'
 print*,'Total disc mass is now ', npartoftype(igas)*massoftype(igas)*umass, ' g'
 print*,'Total disc mass is now ', npartoftype(igas)*massoftype(igas), 'Msun'



 return
end subroutine modify_dump

end module moddump
