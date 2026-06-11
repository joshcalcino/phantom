!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2026 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module moddump
!
! Adds either a sphere or cylinder of material to infall
!
! :References: None
!
! :Owner: Josh Calcino
!
! :Runtime parameters:
!   - add_turbulence : *add turbulence (0=no, 1=yes)*
!   - b              : *impact parameter*
!   - b_frac         : *impact parameter b as fraction of b_crit*
!   - eccentricity   : *eccentricity*
!   - in_mass        : *infall mass*
!   - in_orbit       : *orbit type (0=parabolic, 1=hyperbolic)*
!   - in_shape       : *infall material shape (0=sphere, 1=ellipse)*
!   - incx           : *rotation on x axis (deg)*
!   - incy           : *rotation on y axis (deg)*
!   - incz           : *rotation on z axis (deg)*
!   - r_a            : *semi-major axis of ellipse*
!   - r_close        : *closest approach*
!   - r_in           : *radius of shape (or semi-minor axis)*
!   - r_init         : *initial radial distance*
!   - r_slope        : *density power law index*
!   - r_soft         : *softening radius*
!   - rho_mode       : *density mode (0=current, 1=Dullemond Eq4/Eq5)*
!   - rms_mach       : *rms Mach number*
!   - tfact          : *tfact*
!   - v_inf          : *velocity at infinity (code units)*
!
! :Dependencies: centreofmass, datafiles, dim, eos, infile_utils, io,
!   kernel, moddump_utils, options, part, partinject, physcon, prompting,
!   setvfield, spherical, stretchmap, units, vectorutils, velfield
!
 implicit none
 character(len=*), parameter, public :: moddump_flags = ''

 integer, parameter :: nr = 200
 ! runtime parameters, written to/read from the prefix.mod file.
 ! defaults are set here (module-level), so they are in place before the
 ! driver reads the .mod file via get_moddump_options - do NOT re-set them
 ! at the top of modify_dump or the values read from file would be clobbered
 integer :: in_shape = 1, in_orbit = 1, add_turbulence = 0
 integer :: rho_mode = 0, cloud_control_mode = 0, n_add = 0
 real    :: in_mass = 0.001, r_in = 250.0, r_a = 3500.0, r_init = 4000.0, r_close = 100.0
 real    :: v_inf = 1.0, b = 0.0, b_frac = 1.0, ecc = 0.0
 real    :: incx = 0.0, incy = 0.0, incz = 0.0, rms_mach = 1.0, tfact = 0.0
 real    :: r_slope = 0.0
 real    :: r_soft = 100.0

contains

subroutine modify_dump(npart,npartoftype,massoftype,xyzh,vxyzu)
 use dim,            only:use_dust,maxdusttypes,maxdustlarge,maxdustsmall,use_dustgrowth
 use partinject,     only:add_or_update_particle
 use options,        only:use_dustfrac
 use part,           only:igas,isdead_or_accreted,xyzmh_ptmass,nptmass,ihacc,ihsoft,gravity,&
                          dustfrac
 use units,          only:udist,utime,get_G_code,umass
 use io,             only:id,master,fatal
 use spherical,      only:set_sphere,set_ellipse
 use stretchmap,     only:rho_func
 use kernel,         only:hfact_default
 use physcon,        only:pi,mass_proton_cgs,gg,au
 use vectorutils,    only:rotatevec
 use prompting,      only:prompt
 use centreofmass,   only:reset_centreofmass,get_total_angular_momentum
 use eos,            only:ieos,isink,get_spsound
 use velfield,       only:set_velfield_from_cubes
 use datafiles,      only:find_phantom_datafile
 use setvfield,      only:normalise_vfield
 use moddump_utils,  only:prompt_for_params
 integer, intent(inout) :: npart
 integer, intent(inout) :: npartoftype(:)
 real,    intent(inout) :: massoftype(:)
 real,    intent(inout) :: xyzh(:,:),vxyzu(:,:)
 real, allocatable :: xyzh_add(:,:),vxyzu_add(:,:)
 integer :: ipart,i,np,ierr
 integer(kind=8) :: nptot
 real    :: pmass,delta,r_init_min,big_omega,b_crit
 real    :: vp(3), xp(3), rot_axis(3), rellipsoid(3)
 real    :: dma,n0,pf,m0,x0,y0,z0,r0,vx0,vy0,vz0,mtot,tiny_number,n1
 real    :: y1,x1,dx,x_prime,y_prime
 real    :: unit_velocity,G,rms_in,vol_obj,rhoi,spsound,factor,my_vrms,vxi,vyi,vzi
 real    :: v_inf_cgs,b_crit_cgs
 real    :: rho_cloud_cgs, rho_cloud, mu_cloud, r_equiv
 real    :: dustfrac_tmp
 logical :: lrhofunc,empty_sim,using_modin
 character(len=20), parameter :: filevx = 'cube_v1.dat'
 character(len=20), parameter :: filevy = 'cube_v2.dat'
 character(len=20), parameter :: filevz = 'cube_v3.dat'
 character(len=120)           :: filex,filey,filez
 procedure(rho_func), pointer :: prhofunc

 ! note: the .mod parameters (in_orbit, in_shape, in_mass, ...) are NOT
 ! re-initialised here - their defaults live at module level so they survive
 ! the read performed by the driver before this routine runs. Only local
 ! working variables are initialised below.
 big_omega = 0.
 tiny_number = 1e-4
 lrhofunc = .false.
 empty_sim = .false.
 ierr = 0
 my_vrms = 0.
 pf = 0.
 rho_cloud_cgs = 0.
 rho_cloud = 0.
 mu_cloud = 2.3
 r_equiv = 0.

 ! Gas particle properties
 pmass = massoftype(igas)

 x0 = 0.
 y0 = 0.
 z0 = 0.

 vol_obj = 0.

 ! the prefix.mod parameter file is read by the driver (phantom_moddump) via
 ! get_moddump_options before this routine is called; prompt_for_params tells
 ! us whether to prompt interactively (no file found) or use the values read in
 using_modin = .not. prompt_for_params

 if (npartoftype(igas) <= 0) then
    empty_sim = .true.
    pmass = 0.0
    if (prompt_for_params) then
       write(*,*) "No gas particles detected"
       call prompt('Enter number of particles to add:', n_add, 0)
    endif
 endif

 ! udist default is cm
 unit_velocity = udist/utime ! cm/s
 G = get_G_code()

 if (gravity) then
    write(*,*) "Disc self-gravity is on. Including disc mass in cloud orbit calculation."
    mtot=sum(xyzmh_ptmass(4,:)) + npartoftype(igas)*massoftype(igas)
 else
    mtot=sum(xyzmh_ptmass(4,:))
 endif

 if (prompt_for_params) then
    ! Prompt user for infall material shape
    call prompt('Enter the infall material shape (0=sphere, 1=ellipse)',in_shape,0,1)

    call prompt('Enter cloud control mode (0=manual mass+size (rho not fixed), 1=mass/n_add set radius (fixed rho), '&
               //'2=radius sets mass/n_add (fixed rho))',cloud_control_mode,0,2)

    if (cloud_control_mode == 0 .or. cloud_control_mode == 2) then
       if (in_shape == 0) then
          call prompt('Enter radius of shape:', r_in, 0.1)
       elseif (in_shape == 1) then
          call prompt('Enter semi-minor axis of ellipse:', r_in, 0.1)
          call prompt('Enter semi-major axis of ellipse:', r_a, 0.1)
       endif
    endif
 endif

 if ((cloud_control_mode == 1 .or. cloud_control_mode == 2) .and. empty_sim) then
    write(*,*) "WARNING: Simulation has no mass, massoftype(igas) is not set."
    if (prompt_for_params) then
       write(*,*) "You must set massoftype(igas) to a non-zero value now."
       call prompt('Enter gas particle mass in Msun:', pmass, 0.0)
    endif
    if (pmass <= 0.) then
       call fatal('moddump','pmass must be > 0')
    else
       massoftype(igas) = pmass
    endif
 elseif (cloud_control_mode == 0 .and. empty_sim) then
    if (prompt_for_params) then
       call prompt('Enter infall mass in Msun:', in_mass, 0.0)
    endif
    if (n_add <= 0) call fatal('moddump_infall','n_add must be > 0 for empty simulations')
    pmass = in_mass/real(n_add)
 endif

 if (cloud_control_mode == 0 .or. (cloud_control_mode == 1 .and. .not. using_modin)) then
    if (prompt_for_params) then
       call prompt('Enter infall mass in Msun:', in_mass, 0.0)
    endif
    n_add = int(in_mass/pmass)
 endif

 if (cloud_control_mode == 1 .and. .not. empty_sim) then
    ! Ask how many particles user wants to add
    ! if npartoftype(igas) is 0 and pmass is not set, we set it later
    if (prompt_for_params) then
       call prompt('Enter number of particles to add:', n_add, 0)
    elseif (n_add <= 0 .and. pmass > 0.) then
       n_add = int(in_mass/pmass)
    endif
 endif

!'Enter cloud control mode (0=manual mass+size, 1=N sets size (const rho), '&
!             //'2=size sets mass (const rho))'

 if (cloud_control_mode == 1) then
! n_add sets size (basically in_mass), this sets r_equic.
    in_mass = real(n_add)*pmass
    r_equiv = (in_mass/0.01)**(1./2.3) * 5000.
 elseif (in_shape == 0) then
    r_equiv = r_in
 else
    r_equiv = (r_in*r_in*r_a)**(1./3.)
 endif

 if (cloud_control_mode == 2) then
    ! size sets mass, we already set r_in, r_a, so get the mass
    in_mass = 0.01*(r_equiv/5000.)**2.3
    if (pmass > 0.) n_add = int(in_mass/pmass)
 endif

 if (n_add <= 0) call fatal('moddump_infall','number of particles to add must be > 0')

 vol_obj = (4.0/3.0)*pi*r_equiv**3

 if (cloud_control_mode == 0) then
    if (prompt_for_params) then
       call prompt('Enter value of power-law density along radius:', r_slope, 0.0)
    endif
 else
    r_slope = 0.
    lrhofunc = .false.
 endif
 if (cloud_control_mode == 0 .and. r_slope > tiny_number) then
    prhofunc => rhofunc
    lrhofunc = .true.
    if (prompt_for_params) then
       call prompt('Enter softening radius:', r_soft, 0.1)
    endif
 endif

 ! Prompt user for the infall material orbit
 if (prompt_for_params) then
    call prompt('Enter orbit type (0=parabolic, 1=hyperbolic)', in_orbit,0)
 endif

 if (prompt_for_params) then
    if (in_orbit == 0) then
       print*, "Parabolic orbit"
       call prompt('Enter closest approach in au:', r_close, 0.)
    endif
 endif

 if (in_orbit == 1) then
    write(*,*) "Hyperbolic orbit, see Dullemond+2019 for parameter definitions."
    if (prompt_for_params) then
       call prompt('Enter cloud velocity at infinity, v_inf, in km/s:', v_inf, 0.0)
    endif

    v_inf_cgs = v_inf * 1.e5
    b_crit_cgs = gg * (mtot*umass) / v_inf_cgs**2
    b_crit = b_crit_cgs / au
    write(*,*) "Critical impact parameter, b_crit, is ", b_crit, " au"

    if (prompt_for_params) then
       call prompt('Enter impact parameter b as a ratio of b_crit:', b_frac, 0.0)
    endif
    b = b_frac * b_crit
    ecc = sqrt(1 + b**2/b_crit**2)
    r_close = b * sqrt((ecc-1)/(ecc+1))
    write(*,*) "Eccentricity of the cloud is ", ecc
    write(*,*) "Closest approach of cloud center will be ", r_close, " au."
 endif

 if (prompt_for_params) then
    write(*,*) "Initial radial distance is centre of star/sphere or ellipse."
    if (in_shape == 0) then
       r_init_min = r_in + r_close
    else
       r_init_min = r_a + r_close
    endif
    if (in_orbit == 1) then
       ! Dullemond+2019: L = 4.6 Rcloud, initial cloud centre at x = -1.7 L
       ! => r_init_default = 1.7 * 4.6 * Rcloud
       r_init = 1.7 * 4.6 * r_in
    else
       r_init = 2*r_init_min
    endif
    if (r_init < r_init_min) r_init = r_init_min
    call prompt('Enter initial radial distance in au:', r_init, r_init_min)
 endif

 select case (in_orbit)
 case (0)
    ! Parabolic orbit, taken from set_flyby
    dma = r_close
    ! if (in_shape==1) then
    !   n0  = (r_init+r_a)/r_close
    ! else
    !   n0  = r_init/r_close
    ! endif
    n0  = r_init/r_close
    !--focal parameter dma = pf/2
    pf = 2*dma

    !--define m0 = -x0/dma such that r0 = n0*dma
    !  companion starts at negative x and y
    !  positive root of 1/8*m**4 + m**2 + 2(1-n0**2) = 0
    !  for n0 > 1
    m0 = 2*sqrt(n0-1.0)

    !--perturber initial position
    x0 = -m0*dma
    y0 = dma*(1.0-(x0/pf)**2)
    z0 = 0.0
    xp = (/x0,y0,z0/)
 case (1)
    ! Dullemond+2019
    ! Initial position is x=r_init and y=b (impact parameter)
    if (in_shape==1) then
       x0 = -(r_init+r_in)
    else
       x0 = -r_init
    endif
    y0 = -b
    z0 = 0.0
    xp = (/x0, y0, z0/)
 end select
 write(*,*) "Initial centre is: ", xp

 allocate(xyzh_add(4,n_add+int(0.1*n_add)),vxyzu_add(4,n_add+int(0.1*n_add)))
 delta = 1.0 ! no idea what this is
 nptot = n_add + npartoftype(igas)
 np = 0

 if (in_shape == 0) then
    if (lrhofunc) then
       call set_sphere('random',id,master,0.,r_in,delta,hfact_default,np,xyzh_add,xyz_origin=xp,&
                     np_requested=n_add, nptot=nptot, rhofunc=prhofunc)
    else
       call set_sphere('random',id,master,0.,r_in,delta,hfact_default,np,xyzh_add,xyz_origin=xp,&
                     np_requested=n_add, nptot=nptot)
    endif
    write(*,*) "The sphere has been succesfully initialised."
 elseif (in_shape == 1) then
    rellipsoid(1) = r_in
    rellipsoid(2) = r_a
    rellipsoid(3) = r_in
    call set_ellipse('random',id,master,rellipsoid,delta,hfact_default,xyzh_add,np,&
                    np_requested=n_add, nptot=nptot)

    ! Need to correct the ellipse
    do i = 1,n_add
       xyzh_add(1, i) = xyzh_add(1, i) + xp(1)
       xyzh_add(2, i) = xyzh_add(2, i) + xp(2)
       xyzh_add(3, i) = xyzh_add(3, i) + xp(3)
    enddo
    if (in_orbit == 0) then
       do i = 1,n_add
          x1 = xyzh_add(1, i)
          y1 = xyzh_add(2, i)
          dma = r_close
          n0  = (sqrt(xp(1)**2 + xp(2)**2))/dma
          pf = 2*dma
          n1 = (xp(2)-y1)/dma
          m0 = 2*sqrt(n0-n1-1.0)
          x0 = -m0*dma
          y0 = dma*(1.0-(x0/pf)**2)
          dx = xyzh_add(1, i) - xp(1)
          y_prime = 4*dma/x0 *dx
          x_prime = dx

          xyzh_add(1, i) = x0 + x_prime
          xyzh_add(2, i) = y0 + y_prime
       enddo
    endif
    write(*,*) "The ellipse has been succesfully initialised."
 endif

 !--Set velocities (from pre-made velocity cubes)
 if (prompt_for_params) then
    call prompt('Add turbulence to the gas?:', add_turbulence, 0, 1)
    if (add_turbulence == 1) then
       call prompt('Enter rms Mach number:', rms_mach, 0., 20.)
       call prompt('Enter tfact:', tfact, 0.0)
    endif
 endif
 vxyzu_add(:,:) = 0.

 if (add_turbulence==1) then

    write(*,"(1x,a)") 'Setting up velocity field on the particles...'

    filex = find_phantom_datafile(filevx,'velfield')
    filey = find_phantom_datafile(filevy,'velfield')
    filez = find_phantom_datafile(filevz,'velfield')

    call set_velfield_from_cubes(xyzh_add,vxyzu_add,n_add,filex,filey,filez,1.,tfact*r_in,.false.,ierr)

    if (ierr /= 0) call fatal('setup','error setting up velocity field')

    if (in_shape == 0) then
       vol_obj = (4.0/3.0)*pi*r_in**3
    elseif (in_shape == 1) then
       vol_obj = (4.0/3.0)*pi*r_in*r_in*r_a
    endif

    rhoi = in_mass/vol_obj
    spsound = get_spsound(ieos,xp,rhoi,vxyzu_add(:,1)) ! eos_type,xyzi,rhoi,vxyzui
    rms_in = spsound*rms_mach

    do i=1,n_add
       vxi  = vxyzu_add(1,i)
       vyi  = vxyzu_add(2,i)
       vzi  = vxyzu_add(3,i)
       my_vrms = my_vrms + vxi*vxi + vyi*vyi + vzi*vzi
    enddo

    ! Normalise velocity field
    my_vrms = sqrt(1/real(n_add) * my_vrms)
    factor = rms_in/my_vrms
    do i=1,n_add
       vxyzu_add(1:3,i) = vxyzu_add(1:3,i)*factor
    enddo
 endif

 ! Set up velocities
 if (in_orbit == 0) then

    !--perturber initial velocity
    r0  = sqrt(x0**2+y0**2+z0**2)
    vx0 = (1. + (y0/r0))*sqrt(mtot/pf)
    vy0 = -(x0/r0)*sqrt(mtot/pf)
    vz0 = 0.0
    vp  = (/vx0,vy0,vz0/)
    if (in_shape == 0) then
       ! Initiate initial velocity of the particles in the shape
       vxyzu_add(1, :) = vxyzu_add(1, :) + vx0
       vxyzu_add(2, :) = vxyzu_add(2, :) + vy0
       vxyzu_add(3, :) = vxyzu_add(3, :) + vz0

    elseif (in_shape == 1) then
       do i=1,n_add
          x0 = xyzh_add(1, i)
          y0 = xyzh_add(2, i)
          z0 = xyzh_add(3, i)

          r0  = sqrt(x0**2+y0**2+z0**2)
          vx0 = (1. + (y0/r0))*sqrt(mtot/pf)
          vy0 = -(x0/r0)*sqrt(mtot/pf)
          vz0 = 0.0

          vxyzu_add(1, i) = vxyzu_add(1, i) + vx0
          vxyzu_add(2, i) = vxyzu_add(2, i) + vy0
          vxyzu_add(3, i) = vxyzu_add(3, i) + vz0

       enddo
       vxyzu_add(4, :) = vxyzu(4, 1)
    endif
 elseif (in_orbit == 1) then
    ! Dullemond+2019
    ! Initial velocity, all initially in x direction
    vx0 = v_inf
    vy0 = 0.0
    vz0 = 0.0
    vp = (/vx0, vy0, vz0/)
    vxyzu_add(1, :) = vxyzu_add(1, :) + vx0
    vxyzu_add(2, :) = vxyzu_add(2, :) + vy0
    vxyzu_add(3, :) = vxyzu_add(3, :) + vz0
 endif

 write(*,*) "Initial velocity of object centre is ", vp

 if (use_dust) then
    if (use_dustfrac) then
       write(*,*) "Detected one-fluid dust in the simulation, adding smallest dust to infall."
       ! Set the dustfrac to the global dust to gas ratio

       dustfrac_tmp = sum(dustfrac)/npartoftype(igas)

       write(*,*) "The total dustfrac is ", dustfrac_tmp

       ! Set the dustfrac to the dustfrac of the smallest bin
       dustfrac_tmp = sum(dustfrac(1,:))/npartoftype(igas)
       write(*,*) "The single bin dustfrac is ", dustfrac_tmp

    endif
 endif

 ! Incline the infall
 if (prompt_for_params) then
    write(*,*) "Rotating the infalling gas."
    write(*,*) "Convention: clock-wise rotation in the xy-plane."
    call prompt('Enter rotation on x axis:', incx, -360., 360.)
    call prompt('Enter rotation on y axis:', incy, -360., 360.)
    call prompt('Enter rotation on z axis:', incz, -360., 360.)
    ! call prompt('Enter position angle of ascending node:', big_omega, 0., 360.)
 endif

 ! Now rotate and add those new particles to existing disc
 ipart = npart ! The initial particle number (post shuffle)
 incx = incx*pi/180.
 incy = incy*pi/180.
 incz = incz*pi/180.
 rot_axis = (/1.,1.,0./)
 do i = 1,n_add
    ! Rotate particle to correct position and velocity
    ! First rotate to get the right initial position
    ! Need to do this due to the parabolic orbit notation
    ! xyzh_add(4,i) = 1.0
    if (in_orbit == 0) then
       call rotatevec(xyzh_add(1:3,i),(/0.,-1.,0./),pi)
       call rotatevec(vxyzu_add(1:3,i),(/0.,-1.,0./),pi)
    endif

    ! Now rotate around x axis
    call rotatevec(xyzh_add(1:3,i),(/1.,0.,0./),incx)
    call rotatevec(vxyzu_add(1:3,i),(/1.,0.,0./),incx)

    call rotatevec(xyzh_add(1:3,i),(/0.,1.,0./),incy)
    call rotatevec(vxyzu_add(1:3,i),(/0.,1.,0./),incy)

    call rotatevec(xyzh_add(1:3,i),(/0.,0.,1./),incz)
    call rotatevec(vxyzu_add(1:3,i),(/0.,0.,1./),incz)

    ! Add the particle
    ipart = ipart + 1
    call  add_or_update_particle(igas, xyzh_add(1:3,i), vxyzu_add(1:3,i), xyzh_add(4,i), &
                                vxyzu_add(4,i), ipart, npart, npartoftype, xyzh, vxyzu)
    if (use_dust) then
       if (use_dustfrac) then
          dustfrac(1, ipart) = dustfrac_tmp
       endif
    endif

 enddo

 ! Update if ieos=3 since this will no longer make sense
 if (ieos==3) then
    ! centred at 0,0,0, change to centred on isink=1 if nptmass == 1
    if (nptmass==1) then
       write(*,*) "WARNING: Changing ieos from 3 to 6."
       ieos = 6
       isink = 1
    elseif (nptmass==2) then
       write(*,*) "WARNING: Changing ieos from 3 to 14."
       ieos = 14
    endif
 endif
 write(*,*)  " ###### Added infall successfully ###### "
 deallocate(xyzh_add,vxyzu_add)

end subroutine modify_dump

real function rhofunc(r)
 real, intent(in) :: r

 rhofunc = 1./(abs(r) + r_soft)**(r_slope)

end function rhofunc

!----------------------------------------------------------------
!+
!  read the moddump parameters from the prefix.mod file; nerr
!  counts missing or invalid options (so the caller can top up
!  the file and ask the user to edit it)
!+
!----------------------------------------------------------------
subroutine read_moddump(filename,nerr)
 use infile_utils, only:open_db_from_file,inopts,read_inopt,close_db
 character(len=*), intent(in) :: filename
 integer, intent(out) :: nerr
 type(inopts), allocatable :: db(:)
 integer :: iunit,ierr

 nerr = 0
 iunit = 23
 call open_db_from_file(db,filename,iunit,ierr)
 if (ierr /= 0) then
    nerr = nerr + 1
    return
 endif

 ! control options: read these first as they determine which
 ! dependent options below are required
 call read_inopt(in_shape,'in_shape',db,errcount=nerr,min=0,max=1)
 call read_inopt(in_orbit,'in_orbit',db,errcount=nerr,min=0,max=1)
 call read_inopt(rho_mode,'rho_mode',db,errcount=nerr,min=0)
 call read_inopt(cloud_control_mode,'cloud_control_mode',db,errcount=nerr,min=0,max=2)
 call read_inopt(add_turbulence,'add_turbulence',db,errcount=nerr,min=0,max=1)

 call read_inopt(n_add,'n_add',db,errcount=nerr,min=0)
 call read_inopt(in_mass,'in_mass',db,errcount=nerr,min=0.)
 call read_inopt(r_in,'r_in',db,errcount=nerr,min=0.)
 if (in_shape==1) call read_inopt(r_a,'r_a',db,errcount=nerr,min=0.)

 call read_inopt(r_slope,'r_slope',db,errcount=nerr,min=0.)
 if (abs(r_slope) > 0.) call read_inopt(r_soft,'r_soft',db,errcount=nerr,min=0.)

 call read_inopt(r_init,'r_init',db,errcount=nerr,min=0.)

 if (in_orbit==0) call read_inopt(r_close,'r_close',db,errcount=nerr,min=0.)
 if (in_orbit==1) then
    call read_inopt(v_inf,'v_inf',db,errcount=nerr,min=0.)
    call read_inopt(b_frac,'b_frac',db,errcount=nerr,min=0.)
    call read_inopt(b,'b',db,errcount=nerr,min=0.)
    call read_inopt(ecc,'ecc',db,errcount=nerr,min=0.)
    call read_inopt(r_close,'r_close',db,errcount=nerr,min=0.)
 endif

 call read_inopt(incx,'incx',db,errcount=nerr)
 call read_inopt(incy,'incy',db,errcount=nerr)
 call read_inopt(incz,'incz',db,errcount=nerr)

 if (add_turbulence==1) then
    call read_inopt(rms_mach,'rms_mach',db,errcount=nerr,min=0.)
    call read_inopt(tfact,'tfact',db,errcount=nerr,min=0.)
 endif

 call close_db(db)

end subroutine read_moddump

!----------------------------------------------------------------
!+
!  write the moddump parameters to the prefix.mod file
!+
!----------------------------------------------------------------
subroutine write_moddump(filename)
 use infile_utils,  only:write_inopt
 use physcon,       only:pi
 use moddump_utils, only:moddump_dumpfile_in,moddump_time
 character(len=*), intent(in) :: filename
 integer, parameter :: iunit = 23
 real :: rad_to_deg

 rad_to_deg = 180./pi

 open(unit=iunit,file=filename,status='replace',form='formatted')
 write(iunit,"(a)") '# moddump parameters file for moddump_infall'
 write(iunit,"(2a)") '# input dump file modified: ',trim(moddump_dumpfile_in)
 write(iunit,"(a,g0)") '# output time of modified dump: ',moddump_time
 write(iunit,"(/,a)") '# Infall parameters'
 call write_inopt(in_shape,'in_shape','infall material shape (0=sphere, 1=ellipse)',iunit)
 call write_inopt(in_orbit,'in_orbit','orbit type (0=parabolic, 1=hyperbolic)',iunit)
 call write_inopt(rho_mode,'rho_mode','density mode (0=current, 1=Dullemond Eq4/Eq5)',iunit)
 call write_inopt(cloud_control_mode,'cloud_control_mode',&
                  'cloud control mode (0=manual mass+size, 1=N sets size, 2=size sets mass)',iunit)
 call write_inopt(n_add,'n_add','number of particles added',iunit)
 call write_inopt(in_mass,'in_mass','infall mass',iunit)
 call write_inopt(r_in,'r_in','radius of shape (or semi-minor axis)',iunit)
 if (in_shape==1) call write_inopt(r_a,'r_a','semi-major axis of ellipse',iunit)

 call write_inopt(r_slope,'r_slope','density power law index',iunit)
 if (abs(r_slope) > 0.) call write_inopt(r_soft,'r_soft','softening radius',iunit)

 call write_inopt(r_init,'r_init','initial radial distance',iunit)

 if (in_orbit==0) call write_inopt(r_close,'r_close','closest approach',iunit)
 if (in_orbit==1) then
    call write_inopt(v_inf,'v_inf','velocity at infinity (code units)',iunit)
    call write_inopt(b_frac,'b_frac','impact parameter b as fraction of b_crit',iunit)
    call write_inopt(b,'b','impact parameter',iunit)
    call write_inopt(ecc,'ecc','eccentricity',iunit)
    call write_inopt(r_close,'r_close','closest approach',iunit)
 endif

 call write_inopt(incx*rad_to_deg,'incx','rotation on x axis (deg)',iunit)
 call write_inopt(incy*rad_to_deg,'incy','rotation on y axis (deg)',iunit)
 call write_inopt(incz*rad_to_deg,'incz','rotation on z axis (deg)',iunit)

 call write_inopt(add_turbulence,'add_turbulence','add turbulence (0=no, 1=yes)',iunit)
 if (add_turbulence==1) then
    call write_inopt(rms_mach,'rms_mach','rms Mach number',iunit)
    call write_inopt(tfact,'tfact','tfact',iunit)
 endif

 close(iunit)

end subroutine write_moddump

end module moddump
