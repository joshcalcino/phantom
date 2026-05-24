!--------------------------------------------------------------------------!
! Rectangular BHL wind-tunnel injection for a Phantom disc setup            !
!--------------------------------------------------------------------------!
module inject
!
! Pair this with setup_disc.f90.  The normal Phantom disc setup creates the
! sink(s) and disc.  This module injects/removes the ambient BHL wind.
!
 use bhldisc_options, only:bhl_vinf,bhl_rhoinf,bhl_wind_dir,bhl_box_x,bhl_box_y, &
                           bhl_z_upstream,bhl_z_downstream,bhl_initial_layers, &
                           write_options_bhldisc_inject,read_options_bhldisc_inject
 implicit none
 character(len=*), parameter, public :: inject_type = 'BHLdisc'

 public :: init_inject, inject_particles, write_options_inject, &
           read_options_inject, set_default_options_inject, &
           update_injected_par

 private
 ! With PERIODIC=yes Phantom has some unconditional minimum-image sections.
 ! Make the formal z box huge so those sections never wrap z separations.
 real, parameter :: zbox_factor = 1000.0
 real :: dx_layer, dy_layer, dz_layer, dt_layer
 real :: h_inf, u_inf, cs_inf, v_inf_code, vwind
 real :: mach_inf
 real :: rho_inf_code
 real :: r_bhl, box_x_code, box_y_code, z_in_code, z_out_code
 real :: x_ref, y_ref, z_ref
 real :: xmin_w, xmax_w, ymin_w, ymax_w
 integer :: nx_layer, ny_layer
 integer(kind=8) :: last_layer = -huge(1_8)
 logical :: first_run = .true.
 logical :: have_layer_state = .false.
 logical, parameter :: verbose = .false.

contains

!-----------------------------------------------------------------------
subroutine init_inject(ierr)
 use boundary,  only:set_boundary, print_boundaries
 use dim,       only:maxp, periodic
 use eos,       only:ieos, get_spsound
 use io,        only:fatal
 use mpidomain, only:isperiodic
 use part,      only:hfact, massoftype, igas, nptmass, xyzmh_ptmass
 use physcon,   only:km
 use units,     only:get_G_code, unit_density, unit_velocity
 implicit none
 integer, intent(out) :: ierr
 real :: pmass, element_volume, target_sep, area, mgrav, gcode
 real :: r_bhl_old, xyzi(3), vxyzui(4)
 real :: zmid, zhalf_formal
 integer :: i, iter

 ierr = 0

 if (bhl_vinf <= 0.) call fatal('inject_BHLdisc','bhl_vinf must be positive')
 if (bhl_rhoinf <= 0.) call fatal('inject_BHLdisc','bhl_rhoinf must be positive')
 if (bhl_box_x <= 0.) call fatal('inject_BHLdisc','bhl_box_x must be positive')
 if (bhl_box_y <= 0.) call fatal('inject_BHLdisc','bhl_box_y must be positive')
 if (bhl_z_upstream <= 0.) call fatal('inject_BHLdisc','bhl_z_upstream must be positive')
 if (bhl_z_downstream <= 0.) call fatal('inject_BHLdisc','bhl_z_downstream must be positive')
 if (bhl_wind_dir /= 1 .and. bhl_wind_dir /= -1) then
    call fatal('inject_BHLdisc','bhl_wind_dir must be +1 or -1')
 endif

 if (unit_velocity <= 0.) call fatal('inject_BHLdisc','unit_velocity must be positive')
 v_inf_code = bhl_vinf*km/unit_velocity
 vwind = real(bhl_wind_dir)*v_inf_code
 pmass = massoftype(igas)
 if (pmass <= 0.) then
    call fatal('inject_BHLdisc','massoftype(igas) must be set before injection')
 endif

 u_inf = 0.0

 mgrav = 0.0
 x_ref = 0.0
 y_ref = 0.0
 z_ref = 0.0
 do i = 1, nptmass
    if (xyzmh_ptmass(4,i) > 0.) then
       mgrav = mgrav + xyzmh_ptmass(4,i)
       x_ref = x_ref + xyzmh_ptmass(4,i)*xyzmh_ptmass(1,i)
       y_ref = y_ref + xyzmh_ptmass(4,i)*xyzmh_ptmass(2,i)
       z_ref = z_ref + xyzmh_ptmass(4,i)*xyzmh_ptmass(3,i)
    endif
 enddo
 if (mgrav <= 0.) then
    call fatal('inject_BHLdisc','need at least one positive-mass sink to define R_BHL')
 endif
 x_ref = x_ref/mgrav
 y_ref = y_ref/mgrav
 z_ref = z_ref/mgrav
 gcode = get_G_code()
 if (unit_density <= 0.) call fatal('inject_BHLdisc','unit_density must be positive')
 rho_inf_code = bhl_rhoinf/unit_density

 r_bhl = 2.0*gcode*mgrav/v_inf_code**2
 vxyzui = 0.0
 vxyzui(3) = vwind
 do iter = 1, 20
    r_bhl_old = r_bhl
    if (bhl_wind_dir == 1) then
       xyzi = (/ x_ref, y_ref, z_ref - bhl_z_upstream*r_bhl /)
    else
       xyzi = (/ x_ref, y_ref, z_ref + bhl_z_upstream*r_bhl /)
    endif
    cs_inf = get_spsound(ieos,xyzi,rho_inf_code,vxyzui)
    if (cs_inf < 0.) call fatal('inject_BHLdisc','computed EOS sound speed is negative')
    r_bhl = 2.0*gcode*mgrav/(v_inf_code**2 + cs_inf**2)
    if (abs(r_bhl - r_bhl_old) <= 1.e-10*max(r_bhl,1.0)) exit
 enddo
 if (r_bhl <= 0.) call fatal('inject_BHLdisc','computed R_BHL is not positive')
 mach_inf = v_inf_code/max(cs_inf,tiny(cs_inf))

 box_x_code = bhl_box_x*r_bhl
 box_y_code = bhl_box_y*r_bhl
 if (bhl_wind_dir == 1) then
    z_in_code  = z_ref - bhl_z_upstream*r_bhl
    z_out_code = z_ref + bhl_z_downstream*r_bhl
 else
    z_in_code  = z_ref + bhl_z_upstream*r_bhl
    z_out_code = z_ref - bhl_z_downstream*r_bhl
 endif

 element_volume = pmass/rho_inf_code
 target_sep = element_volume**(1.0/3.0)
 h_inf = hfact*target_sep

 ! Fill the x/y periodic rectangle with an integer lattice.  Then choose the
 ! z spacing so that dx*dy*dz = pmass/rho_inf exactly.
 nx_layer = max(1, nint(box_x_code/target_sep))
 ny_layer = max(1, nint(box_y_code/target_sep))
 dx_layer = box_x_code/real(nx_layer)
 dy_layer = box_y_code/real(ny_layer)
 dz_layer = element_volume/(dx_layer*dy_layer)
 dt_layer = dz_layer/v_inf_code

 xmin_w = x_ref - 0.5*box_x_code
 xmax_w = x_ref + 0.5*box_x_code
 ymin_w = y_ref - 0.5*box_y_code
 ymax_w = y_ref + 0.5*box_y_code

 if (periodic) then
    isperiodic = (/ .true., .true., .false. /)
    zmid = 0.5*(z_in_code + z_out_code)
    zhalf_formal = 0.5*zbox_factor*abs(z_out_code - z_in_code)
    call set_boundary(xmin_w, xmax_w, ymin_w, ymax_w, &
                      zmid - zhalf_formal, zmid + zhalf_formal)
    call print_boundaries(6,.true.)
 endif

 if (int(max(1,bhl_initial_layers),kind=8)*int(nx_layer,kind=8) &
       *int(ny_layer,kind=8) > int(maxp,kind=8)) then
    call fatal('inject_BHLdisc','maxp too small for the initial injected layers')
 endif

 area = box_x_code*box_y_code
 print*, 'BHLdisc rectangular wind injection'
 print*, '  v_inf km/s, code; mach          = ', bhl_vinf, v_inf_code, mach_inf
 print*, '  EOS c_s km/s, code at inlet     = ', cs_inf*unit_velocity/km, cs_inf
 print*, '  rho_inf cgs, code               = ', bhl_rhoinf, rho_inf_code
 print*, '  R_BHL=2GM/(v^2+c_s^2), ref xyz  = ', r_bhl, x_ref, y_ref, z_ref
 print*, '  u_inf, h_inf                    = ', u_inf, h_inf
 print*, '  x/y layer particles             = ', nx_layer, ny_layer
 print*, '  dx, dy, dz, dt_layer            = ', dx_layer, dy_layer, dz_layer, dt_layer
 print*, '  z_in, z_out, wind_dir           = ', z_in_code, z_out_code, bhl_wind_dir
 print*, '  wind mass flux                  = ', rho_inf_code*v_inf_code*area

end subroutine init_inject

!-----------------------------------------------------------------------
subroutine inject_particles(time, dtlast, xyzh, vxyzu, xyzmh_ptmass, &
                            vxyz_ptmass, npart, npart_old, &
                            npartoftype, dtinject)
 real,    intent(in)    :: time, dtlast
 real,    intent(inout) :: xyzh(:,:), vxyzu(:,:)
 real,    intent(inout) :: xyzmh_ptmass(:,:), vxyz_ptmass(:,:)
 integer, intent(inout) :: npart, npart_old
 integer, intent(inout) :: npartoftype(:)
 real,    intent(out)   :: dtinject
 integer :: ierr
 integer(kind=8) :: target_layer, layer

 if (first_run) then
    call init_inject(ierr)
    first_run = .false.
 endif

 call delete_outflow_particles(npart, npartoftype, xyzh)

 ! Deletion compacts the particle array.  Tell Phantom that only particles
 ! appended after this point are new injected particles.
 npart_old = npart

 if (.not.have_layer_state) then
    if (dtlast <= 0. .and. abs(time) <= 10.0*tiny(time)) then
       last_layer = -int(max(0,bhl_initial_layers),kind=8)
    else
       ! Restart: do not duplicate the inlet reservoir from the dump.
       last_layer = int(floor(time/dt_layer),kind=8)
    endif
    have_layer_state = .true.
 endif

 target_layer = int(floor(time/dt_layer),kind=8)
 do layer = last_layer + 1_8, target_layer
    call append_rectangular_layer(layer, time, npart, npartoftype, xyzh, vxyzu)
 enddo
 last_layer = max(last_layer, target_layer)

 dtinject = 0.5*dt_layer

end subroutine inject_particles

!-----------------------------------------------------------------------
subroutine append_rectangular_layer(layer, time, npart, npartoftype, &
                                    xyzh, vxyzu)
 use part,       only:igas
 use partinject, only:add_or_update_particle
 implicit none
 integer(kind=8), intent(in)    :: layer
 real,            intent(in)    :: time
 integer,         intent(inout) :: npart, npartoftype(:)
 real,            intent(inout) :: xyzh(:,:), vxyzu(:,:)
 integer :: ix, iy
 real :: xyzi(3), vxyzi(3), zlayer

 zlayer = z_in_code + vwind*(time - real(layer)*dt_layer)
 if (.not.inside_open_z_domain(zlayer)) return

 vxyzi = (/0.0, 0.0, vwind/)
 xyzi(3) = zlayer
 do iy = 1, ny_layer
    xyzi(2) = ymin_w + (real(iy)-0.5)*dy_layer
    do ix = 1, nx_layer
       xyzi(1) = xmin_w + (real(ix)-0.5)*dx_layer
       call add_or_update_particle(igas, xyzi, vxyzi, h_inf, u_inf, &
                                   npart+1, npart, npartoftype, &
                                   xyzh, vxyzu)
    enddo
 enddo

 if (verbose) print*, 'BHLdisc injected layer ', layer, ' z=', zlayer

end subroutine append_rectangular_layer

!-----------------------------------------------------------------------
logical function inside_open_z_domain(z)
 implicit none
 real, intent(in) :: z
 if (bhl_wind_dir > 0) then
    inside_open_z_domain = (z >= z_in_code .and. z <= z_out_code)
 else
    inside_open_z_domain = (z <= z_in_code .and. z >= z_out_code)
 endif
end function inside_open_z_domain

!-----------------------------------------------------------------------
subroutine delete_outflow_particles(npart, npartoftype, xyzh)
 use part, only:igas, iamtype, iphase, kill_particle, shuffle_part
 implicit none
 integer, intent(inout) :: npart, npartoftype(:)
 real,    intent(inout) :: xyzh(:,:)
 integer :: i, nkill
 logical :: killme

 nkill = 0
 do i = 1, npart
    if (iamtype(iphase(i)) /= igas) cycle
    if (xyzh(4,i) <= 0.) cycle
    if (bhl_wind_dir > 0) then
       killme = (xyzh(3,i) > z_out_code)
    else
       killme = (xyzh(3,i) < z_out_code)
    endif
    if (killme) then
       call kill_particle(i,npartoftype)
       nkill = nkill + 1
    endif
 enddo
 if (nkill > 0) call shuffle_part(npart)
 if (verbose .and. nkill > 0) print*, 'BHLdisc deleted ', nkill, ' particles'

end subroutine delete_outflow_particles

!-----------------------------------------------------------------------
subroutine write_options_inject(iunit)
 implicit none
 integer, intent(in) :: iunit

 call write_options_bhldisc_inject(iunit)

end subroutine write_options_inject

!-----------------------------------------------------------------------
subroutine read_options_inject(db,nerr)
 use infile_utils, only:inopts
 implicit none
 type(inopts), intent(inout) :: db(:)
 integer,      intent(inout) :: nerr

 call read_options_bhldisc_inject(db,nerr)

end subroutine read_options_inject

!-----------------------------------------------------------------------
subroutine set_default_options_inject(flag)
 implicit none
 integer, intent(in), optional :: flag
end subroutine set_default_options_inject

!-----------------------------------------------------------------------
subroutine update_injected_par
 implicit none
end subroutine update_injected_par

end module inject
