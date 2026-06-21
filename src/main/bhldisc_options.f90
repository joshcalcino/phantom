!--------------------------------------------------------------------------!
! Shared options for the BHLdisc setup/injector pair                       !
!--------------------------------------------------------------------------!
module bhldisc_options
 implicit none

 public :: bhl_vinf, bhl_rhoinf, bhl_wind_dir, bhl_box_x, bhl_box_y
 public :: bhl_z_upstream, bhl_z_downstream, bhl_initial_layers
 public :: gas_mass_scale, bhl_dust_to_gas
 public :: write_options_bhldisc_setup, read_options_bhldisc_setup
 public :: write_options_bhldisc_inject, read_options_bhldisc_inject

 ! Wind speed is in km/s; density is in cgs.
 real    :: bhl_vinf     = 1.0
 real    :: bhl_rhoinf   = 3.0e-20
 integer :: bhl_wind_dir = 1

 ! Rectangular wind tunnel, in Bondi-Hoyle radii.
 real :: bhl_box_x = 8.0
 real :: bhl_box_y = 8.0
 real :: bhl_z_upstream   = 4.0
 real :: bhl_z_downstream = 4.0

 ! Number of layers to prefill at t=0.
 integer :: bhl_initial_layers = 4

 ! BHLdisc-only global gas particle mass scale.
 real :: gas_mass_scale = 1.0

 ! Dust seeded into the wind: smallest bin only, one-fluid (dustfrac) runs only.
 real :: bhl_dust_to_gas = 0.01

contains

!-----------------------------------------------------------------------
subroutine write_options_bhldisc_setup(iunit)
 use infile_utils, only:write_inopt
 implicit none
 integer, intent(in) :: iunit

 write(iunit,"(/,a)") '# BHL wind tunnel setup'
 call write_options_bhldisc_common(iunit)
 call write_inopt(gas_mass_scale,'gas_mass_scale','multiply all gas particle masses by this factor',iunit)

end subroutine write_options_bhldisc_setup

!-----------------------------------------------------------------------
subroutine read_options_bhldisc_setup(db,nerr)
 use infile_utils, only:inopts, read_inopt
 use dim,          only:use_dust
 implicit none
 type(inopts), intent(inout) :: db(:)
 integer,      intent(inout) :: nerr

 call read_inopt(bhl_vinf,'bhl_vinf',db,errcount=nerr,min=epsilon(0.),default=bhl_vinf)
 call read_inopt(bhl_rhoinf,'bhl_rhoinf',db,errcount=nerr,min=tiny(0.),default=bhl_rhoinf)
 call read_inopt(bhl_wind_dir,'bhl_wind_dir',db,errcount=nerr,min=-1,max=1,default=bhl_wind_dir)
 call read_inopt(bhl_box_x,'bhl_box_x',db,errcount=nerr,min=epsilon(0.),default=bhl_box_x)
 call read_inopt(bhl_box_y,'bhl_box_y',db,errcount=nerr,min=epsilon(0.),default=bhl_box_y)
 call read_inopt(bhl_z_upstream,'bhl_z_upstream',db,errcount=nerr,min=epsilon(0.),default=bhl_z_upstream)
 call read_inopt(bhl_z_downstream,'bhl_z_downstream',db,errcount=nerr,min=epsilon(0.),default=bhl_z_downstream)
 call read_inopt(bhl_initial_layers,'bhl_initial_layers',db,errcount=nerr,min=0,default=bhl_initial_layers)
 call read_inopt(gas_mass_scale,'gas_mass_scale',db,min=tiny(0.),errcount=nerr,default=gas_mass_scale)
 if (use_dust) call read_inopt(bhl_dust_to_gas,'bhl_dust_to_gas',db,errcount=nerr,min=0.,default=bhl_dust_to_gas)

end subroutine read_options_bhldisc_setup

!-----------------------------------------------------------------------
subroutine write_options_bhldisc_inject(iunit)
 use infile_utils, only:write_inopt
 implicit none
 integer, intent(in) :: iunit

 call write_options_bhldisc_common(iunit)
 call write_inopt(gas_mass_scale,'gas_mass_scale','BHL gas mass scale from setup',iunit)

end subroutine write_options_bhldisc_inject

!-----------------------------------------------------------------------
subroutine read_options_bhldisc_inject(db,nerr)
 use infile_utils, only:inopts, read_inopt
 use dim,          only:use_dust
 implicit none
 type(inopts), intent(inout) :: db(:)
 integer,      intent(inout) :: nerr

 call read_inopt(bhl_vinf,'bhl_vinf',db,errcount=nerr,min=epsilon(0.))
 call read_inopt(bhl_rhoinf,'bhl_rhoinf',db,errcount=nerr,min=tiny(0.))
 call read_inopt(bhl_wind_dir,'bhl_wind_dir',db,errcount=nerr,min=-1,max=1)
 call read_inopt(bhl_box_x,'bhl_box_x',db,errcount=nerr,min=epsilon(0.))
 call read_inopt(bhl_box_y,'bhl_box_y',db,errcount=nerr,min=epsilon(0.))
 call read_inopt(bhl_z_upstream,'bhl_z_upstream',db,errcount=nerr,min=epsilon(0.))
 call read_inopt(bhl_z_downstream,'bhl_z_downstream',db,errcount=nerr,min=epsilon(0.))
 call read_inopt(bhl_initial_layers,'bhl_initial_layers',db,errcount=nerr,min=0)
 call read_inopt(gas_mass_scale,'gas_mass_scale',db,min=tiny(0.),errcount=nerr,default=gas_mass_scale)
 if (use_dust) call read_inopt(bhl_dust_to_gas,'bhl_dust_to_gas',db,errcount=nerr,min=0.,default=bhl_dust_to_gas)

end subroutine read_options_bhldisc_inject

!-----------------------------------------------------------------------
subroutine write_options_bhldisc_common(iunit)
 use infile_utils, only:write_inopt
 use dim,          only:use_dust
 implicit none
 integer, intent(in) :: iunit

 call write_inopt(bhl_vinf,'bhl_vinf','BHL wind speed in km/s',iunit)
 call write_inopt(bhl_rhoinf,'bhl_rhoinf','BHL ambient density in g cm^-3',iunit)
 call write_inopt(bhl_wind_dir,'bhl_wind_dir','+1 or -1 flow direction in z',iunit)
 call write_inopt(bhl_box_x,'bhl_box_x','periodic/injection box size in x [R_BHL]',iunit)
 call write_inopt(bhl_box_y,'bhl_box_y','periodic/injection box size in y [R_BHL]',iunit)
 call write_inopt(bhl_z_upstream,'bhl_z_upstream','inlet distance upstream of sink COM [R_BHL]',iunit)
 call write_inopt(bhl_z_downstream,'bhl_z_downstream','outlet distance downstream of sink COM [R_BHL]',iunit)
 call write_inopt(bhl_initial_layers,'bhl_initial_layers','wind layers prefilled at t=0',iunit)
 if (use_dust) call write_inopt(bhl_dust_to_gas,'bhl_dust_to_gas',&
    'dust-to-gas ratio in the wind (one-fluid dust only)',iunit)

end subroutine write_options_bhldisc_common

end module bhldisc_options
