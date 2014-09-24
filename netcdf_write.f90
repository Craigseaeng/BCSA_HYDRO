SUBROUTINE netcdf_write

! Writes EFDC output to NetCDF file

! Set for 3D variables only

! 7/3/2014 Chris Flanary

USE GLOBAL
USE netcdf
IMPLICIT NONE

LOGICAL,SAVE::FIRST_NETCDF=.FALSE.
CHARACTER (len = *), PARAMETER :: FILE_NAME = "efdc_his.nc"
INTEGER, SAVE :: nc_step, ncid
INTEGER :: I, J, L, S, status
INTEGER, SAVE :: I_dimid,J_dimid,k_dimid,time_dimid,mask_varid
INTEGER, SAVE :: ts_varid,time_varid,X_varid,Y_varid
INTEGER, SAVE :: surfel_varid,u_varid,v_varid,sal_varid, dye_varid
INTEGER, SAVE :: tss_varid, tau_varid,d50_varid,thick_varid

INTEGER :: nstep, start(1), start_3d(3), start_4d(4)
REAL, DIMENSION(1) :: deltat, time_efdc
REAL, DIMENSION(LCM) :: zeta,wet_dry_mask

! Allocate 2D array variables
REAL, ALLOCATABLE, DIMENSION(:,:) :: lat,lon,mask,wl,u_2d,v_2d,salt,dye_2d
REAL, ALLOCATABLE, DIMENSION(:,:) :: bed_shear,grain_size,sed_thick
REAL, ALLOCATABLE, DIMENSION(:,:,:) :: tss
ALLOCATE(mask(JC-2,IC-2))
ALLOCATE(lat(JC-2,IC-2))
ALLOCATE(lon(JC-2,IC-2))
ALLOCATE(wl(JC-2,IC-2))
ALLOCATE(u_2d(JC-2,IC-2))
ALLOCATE(v_2d(JC-2,IC-2))
ALLOCATE(salt(JC-2,IC-2))
ALLOCATE(dye_2d(JC-2,IC-2))
ALLOCATE(bed_shear(JC-2,IC-2))
ALLOCATE(grain_size(JC-2,IC-2))
ALLOCATE(sed_thick(JC-2,IC-2))
ALLOCATE(tss(JC-2,IC-2,NSCM))

! EFDC time parameter
time_efdc=DT*FLOAT(N)+TCON*TBEGIN
time_efdc=time_efdc/86400.

! Timing parameters
deltat=tidalp/float(ntsptc)
nc_step=nc_step+1

! Current write step
start=(/nc_step/)
start_3d=(/1,1,nc_step/)
start_4d=(/1,1,1,nc_step/)

DO L=2,LA
    ! Calculate surface elevations for all active cells
    zeta(L)=(HP(L)+BELV(L))

    ! Create wet dry mask
    IF (HP(L)<=HDRY) THEN
        wet_dry_mask(L)=-99.
    ELSE
        wet_dry_mask(L)=1.
    ENDIF

ENDDO

! Create 2D arrays for all time variable parameters
DO I=3,IC-2
   DO J=3,JC-2
      IF(LIJ(I,J)>0) THEN
		 mask(J,I)=wet_dry_mask(LIJ(I,J))
         wl(J,I)=zeta(LIJ(I,J))
         u_2d(J,I)=U(LIJ(I,J),1)
         v_2d(J,I)=V(LIJ(I,J),1)
         salt(J,I)=SAL(LIJ(I,J),1)
         dye_2d(J,I)=DYE(LIJ(I,J),1)
         bed_shear(J,I)=TAU(LIJ(I,J))
         grain_size(J,I)=D50AVG(LIJ(I,J))
         sed_thick(J,I)=THCK(LIJ(I,J))
         DO S=1,NSCM
            tss(J,I,S)=SED(LIJ(I,J),1,S)
         ENDDO
      ELSE
         ! Flag for inactive cells
		 mask(J,I)=-7999.
         wl(J,I)=-7999.
         u_2d(J,I)=-7999.
         v_2d(J,I)=-7999.
         salt(J,I)=-7999.
         dye_2d(J,I)=-7999.
         bed_shear(J,I)=-7999.
         grain_size(J,I)=-7999.
         sed_thick(J,I)=-7999.
      ENDIF
   ENDDO
ENDDO

IF(.NOT.FIRST_NETCDF)THEN

   ! Create 2D lat lon array for all cells
   DO I=3,IC-2
      DO J=3,JC-2
         IF(LIJ(I,J)>0) THEN
            lat(J,I)=DLAT(LIJ(I,J))
            lon(J,I)=DLON(LIJ(I,J))
         ELSE
            lat(J,I)=-7999.
            lon(J,I)=-7999.
         ENDIF
      ENDDO
   ENDDO

    ! Create NetCDF file and overwrite if exists
    status=nf90_create(FILE_NAME, NF90_CLOBBER, ncid)
    if(status /= nf90_NoErr) call handle_err(status)

    !*************************************************************************************
    ! Variable definitions

    ! Define global attributes
    status=nf90_put_att(ncid, NF90_GLOBAL, 'title', 'test file')
    status=nf90_put_att(ncid, NF90_GLOBAL, 'file', FILE_NAME)
    status=nf90_put_att(ncid, NF90_GLOBAL, 'format', 'netCDF-3 64bit offset file')
    status=nf90_put_att(ncid, NF90_GLOBAL, 'os', 'Linux')
    status=nf90_put_att(ncid, NF90_GLOBAL, 'arch', 'x86_64')

    ! Define deltat
    status=nf90_def_var(ncid,'timestep',nf90_real,ts_varid)
    status=nf90_put_att(ncid, ts_varid, 'long_name', 'Numerical Timestep')
    status=nf90_put_att(ncid, ts_varid, 'units', 's')

    ! Define dimensions
    status=nf90_def_dim(ncid,'I',IC-2,I_dimid)
    status=nf90_def_dim(ncid,'J',JC-2,J_dimid)
    status=nf90_def_dim(ncid,'NSCM',NSCM,k_dimid)
    status=nf90_def_dim(ncid,'efdc_time',nf90_unlimited,time_dimid)

    ! Define model time
    status=nf90_def_var(ncid,'efdc_time',nf90_real,time_dimid,time_varid)
    status=nf90_put_att(ncid, time_varid, 'long_name', 'Time')
    status=nf90_put_att(ncid, time_varid, 'units', 'Days')
	status=nf90_put_att(ncid, time_varid, 'ref_time', '20090516')
    status=nf90_put_att(ncid, time_varid, 'calendar', 'gregorian')

    ! Define coordinates variables
    status=nf90_def_var(ncid,'X',nf90_float,(/J_dimid,I_dimid/),X_varid)
    status=nf90_put_att(ncid, X_varid, 'long_name', 'X coordinate')
    status=nf90_put_att(ncid, X_varid, 'units', 'm')
	status=nf90_put_att(ncid, X_varid, 'coord_sys', 'UTM')
    status=nf90_put_att(ncid, X_varid, 'fill_value', -7999)

    status=nf90_def_var(ncid,'Y',nf90_float,(/J_dimid,I_dimid/),Y_varid)
    status=nf90_put_att(ncid, Y_varid, 'long_name', 'Y coordinate')
    status=nf90_put_att(ncid, Y_varid, 'units', 'm')
	status=nf90_put_att(ncid, Y_varid, 'coord_sys', 'UTM')
    status=nf90_put_att(ncid, Y_varid, 'fill_value', -7999)
	
	! Define wet dry mask
    status=nf90_def_var(ncid,'wet_dry_mask',nf90_real,(/J_dimid,I_dimid, time_dimid/),mask_varid)
    status=nf90_put_att(ncid, mask_varid, 'long_name', 'Wet Dry Mask')
	status=nf90_put_att(ncid, mask_varid, 'caxis_lablel', 'Wet Dry Mask')
    status=nf90_put_att(ncid, mask_varid, 'dry_flag', '-99')
	status=nf90_put_att(ncid, mask_varid, 'wet_flag', '1')

    ! Define water surface elevation
    status=nf90_def_var(ncid,'zeta',nf90_real,(/J_dimid,I_dimid, time_dimid/),surfel_varid)
    status=nf90_put_att(ncid, surfel_varid, 'long_name', 'Surface Elevation')
	status=nf90_put_att(ncid, surfel_varid, 'caxis_label', 'Surface Elevation (m)')
    status=nf90_put_att(ncid, surfel_varid, 'units', 'm')

    ! Define u component velocity
    status=nf90_def_var(ncid,'u',nf90_float,(/J_dimid,I_dimid, time_dimid/),u_varid)
    status=nf90_put_att(ncid, u_varid, 'long_name', 'U Component Velocity')
    status=nf90_put_att(ncid, u_varid, 'units', 'm/s')

    ! Define v component velocity
    status=nf90_def_var(ncid,'v',nf90_float,(/J_dimid,I_dimid, time_dimid/),v_varid)
    status=nf90_put_att(ncid, v_varid, 'long_name', 'V Component Velocity')
    status=nf90_put_att(ncid, v_varid, 'units', 'm/s')

    ! Define salinity
    status=nf90_def_var(ncid,'salt',nf90_float,(/J_dimid,I_dimid, time_dimid/),sal_varid)
    status=nf90_put_att(ncid, sal_varid, 'long_name', 'Salinity')
	status=nf90_put_att(ncid, sal_varid, 'caxis_label', 'Salinity (psu)')
    status=nf90_put_att(ncid, sal_varid, 'units', 'psu')

    ! Define dye concentration
    status=nf90_def_var(ncid,'dye',nf90_float,(/J_dimid,I_dimid, time_dimid/),dye_varid)
    status=nf90_put_att(ncid, dye_varid, 'long_name', 'Dye Concentration')
	status=nf90_put_att(ncid, dye_varid, 'caxis_label', 'Dye Concentration (mg/L)')
    status=nf90_put_att(ncid, dye_varid, 'units', 'mg/L')

    ! Define water column TSS
    status=nf90_def_var(ncid,'tss',nf90_float,(/J_dimid,I_dimid,k_dimid,time_dimid/),tss_varid)
    status=nf90_put_att(ncid, tss_varid, 'long_name', 'TSS')
	status=nf90_put_att(ncid, tss_varid, 'caxis_label', 'TSS (mg/L)')
    status=nf90_put_att(ncid, tss_varid, 'units', 'mg/L')

    ! Define bed shear
    status=nf90_def_var(ncid,'tau',nf90_float,(/J_dimid,I_dimid, time_dimid/),tau_varid)
    status=nf90_put_att(ncid, tau_varid, 'long_name', 'Shear Stress')
	status=nf90_put_att(ncid, tau_varid, 'caxis_label', 'Shear Stress (dynes/cm^2)')
    status=nf90_put_att(ncid, tau_varid, 'units', 'dynes/cm^2')

    ! Define D50 grain size
    status=nf90_def_var(ncid,'D50',nf90_float,(/J_dimid,I_dimid, time_dimid/),d50_varid)
    status=nf90_put_att(ncid, d50_varid, 'long_name', 'Median Grain Size')
	status=nf90_put_att(ncid, d50_varid, 'long_name', 'Median Grain Size (um)')
    status=nf90_put_att(ncid, d50_varid, 'units', 'um')

    ! Define sediment thickness
    status=nf90_def_var(ncid,'thickness',nf90_float,(/J_dimid,I_dimid, time_dimid/),thick_varid)
    status=nf90_put_att(ncid, thick_varid, 'long_name', 'Bed Thickness')
	status=nf90_put_att(ncid, thick_varid, 'caxis_label', 'Bed Thickness (m)')
    status=nf90_put_att(ncid, thick_varid, 'units', 'm')

    ! End define mode
    status=nf90_enddef(ncid)

    !*************************************************************************************
    ! Put first step variables into file

    ! Put EFDC time step into file once
    status=nf90_put_var(ncid, ts_varid, deltat)
    if(status /= nf90_NoErr) call handle_err(status)

    ! Put first EFDC time in days into file
    status=nf90_put_var(ncid, time_varid, time_efdc, start=start)
    if(status /= nf90_NoErr) call handle_err(status)

    ! Put coordinates into file
    status=nf90_put_var(ncid,X_varid,lon)
    if(status /= nf90_NoErr) call handle_err(status)
    status=nf90_put_var(ncid,Y_varid,lat)
    if(status /= nf90_NoErr) call handle_err(status)
	
	! Put first wet dry mask into file
    status=nf90_put_var(ncid,mask_varid,mask,start=start_3d)
    if(status /= nf90_NoErr) call handle_err(status)

    ! Put first surface elevation into file
    status=nf90_put_var(ncid,surfel_varid,wl,start=start_3d)
    if(status /= nf90_NoErr) call handle_err(status)

    ! Put first u component into file
    status=nf90_put_var(ncid,u_varid,u_2d,start=start_3d)
    if(status /= nf90_NoErr) call handle_err(status)
    
    ! Put first v component into file
    status=nf90_put_var(ncid,v_varid,v_2d,start=start_3d)
    if(status /= nf90_NoErr) call handle_err(status)

    ! Put first salinity into file
    status=nf90_put_var(ncid,sal_varid,salt,start=start_3d)
    if(status /= nf90_NoErr) call handle_err(status)

    ! Put first dye into file
    status=nf90_put_var(ncid,dye_varid,dye_2d,start=start_3d)
    if(status /= nf90_NoErr) call handle_err(status)

    ! Put first TSS into file
    status=nf90_put_var(ncid,tss_varid,tss,start=start_4d)
    if(status /= nf90_NoErr) call handle_err(status)

    ! Put first tau into file
    status=nf90_put_var(ncid,tau_varid,bed_shear,start=start_3d)
    if(status /= nf90_NoErr) call handle_err(status)

    ! Put first D50 into file
    status=nf90_put_var(ncid,d50_varid,grain_size,start=start_3d)
    if(status /= nf90_NoErr) call handle_err(status)

    ! Put first sediment thickness into file
    status=nf90_put_var(ncid,thick_varid,sed_thick,start=start_3d)
    if(status /= nf90_NoErr) call handle_err(status)

    FIRST_NETCDF=.TRUE.
    
    ! Close NetCDF file after first write
    status=nf90_close(ncid)
    if(status /= nf90_NoErr) call handle_err(status)

ELSE

   !*************************************************************************************
    ! Put subsequent step variables into file

    ! Open NetCDF file
    status=nf90_open(FILE_NAME, nf90_write, ncid)
    if(status /= nf90_NoErr) call handle_err(status)

    ! Put EFDC times into file
    status=nf90_put_var(ncid, time_varid, time_efdc, start=start)
    if(status /= nf90_NoErr) call handle_err(status)
	
	! Put wet dry mask into file
    status=nf90_put_var(ncid,mask_varid,mask,start=start_3d)
    if(status /= nf90_NoErr) call handle_err(status)

    ! Put surface elevations into file
    status=nf90_put_var(ncid,surfel_varid,wl,start=start_3d)
    if(status /= nf90_NoErr) call handle_err(status)

    ! Put u components into file
    status=nf90_put_var(ncid,u_varid,u_2d,start=start_3d)
    if(status /= nf90_NoErr) call handle_err(status)
    
    ! Put v components into file
    status=nf90_put_var(ncid,v_varid,v_2d,start=start_3d)
    if(status /= nf90_NoErr) call handle_err(status)

    ! Put salinities into file
    status=nf90_put_var(ncid,sal_varid,salt,start=start_3d)
    if(status /= nf90_NoErr) call handle_err(status)

    ! Put dye into file
    status=nf90_put_var(ncid,dye_varid,dye_2d,start=start_3d)
    if(status /= nf90_NoErr) call handle_err(status)

    ! Put TSS into file
    status=nf90_put_var(ncid,tss_varid,tss,start=start_4d)
    if(status /= nf90_NoErr) call handle_err(status)

    ! Put tau into file
    status=nf90_put_var(ncid,tau_varid,bed_shear,start=start_3d)
    if(status /= nf90_NoErr) call handle_err(status)

    ! Put D50 into file
    status=nf90_put_var(ncid,d50_varid,grain_size,start=start_3d)
    if(status /= nf90_NoErr) call handle_err(status)

    ! Put sediment thickness into file
    status=nf90_put_var(ncid,thick_varid,sed_thick,start=start_3d)
    if(status /= nf90_NoErr) call handle_err(status)

    ! Close NetCDF file each time
    status=nf90_close(ncid)
    if(status /= nf90_NoErr) call handle_err(status)

ENDIF

CONTAINS
    SUBROUTINE handle_err(status)
        INTEGER, INTENT ( in) :: status

        IF(status /= nf90_noerr) THEN
            PRINT *, TRIM(nf90_strerror(status))
            STOP 'Stopped'
        END IF
    END SUBROUTINE handle_err
END SUBROUTINE netcdf_write
