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
INTEGER :: I, J, K, L, S, status
INTEGER, SAVE :: I_dimid,J_dimid,k_dimid,time_dimid
INTEGER, SAVE :: ts_varid,time_varid,X_varid,Y_varid,depth_varid,mask_varid
INTEGER, SAVE :: surfel_varid,u_varid,v_varid,sal_varid,dye_varid
INTEGER, SAVE :: tss_varid,tau_varid,taumax_varid,d50_varid,thick_varid
INTEGER, SAVE :: vmax_varid

INTEGER :: start(1), start_3d(3), start_4d(4)
REAL, DIMENSION(1) :: deltat, time_efdc
REAL, DIMENSION(LCM) :: zeta,wet_dry_mask
REAL, DIMENSION(LCM) :: utmps,vtmps
REAL :: utmpa,vtmpa
REAL, DIMENSION(LCM) :: vel_max,vel_magc

! Allocate 2D array variables
REAL, ALLOCATABLE, DIMENSION(:,:) :: lat,lon,mask,hz,wl,u_2d,v_2d,salt,dye_2d
REAL, ALLOCATABLE, DIMENSION(:,:) :: shear,maxshear,grain_size,sed_thick,mag
REAL, ALLOCATABLE, DIMENSION(:,:,:) :: tss
ALLOCATE(mask(JC-2,IC-2))
ALLOCATE(hz(JC-2,IC-2))
ALLOCATE(lat(JC-2,IC-2))
ALLOCATE(lon(JC-2,IC-2))
ALLOCATE(wl(JC-2,IC-2))
ALLOCATE(u_2d(JC-2,IC-2))
ALLOCATE(v_2d(JC-2,IC-2))
ALLOCATE(mag(JC-2,IC-2))
ALLOCATE(salt(JC-2,IC-2))
ALLOCATE(dye_2d(JC-2,IC-2))
ALLOCATE(shear(JC-2,IC-2))
ALLOCATE(maxshear(JC-2,IC-2))
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

! Create NetCDF file and define attributes
IF(.NOT.FIRST_NETCDF)THEN

   ! Create 2D lat, lon, and depth array for all cells
   DO I=3,IC-2
      DO J=3,JC-2
         IF(LIJ(I,J)>0) THEN
            lat(J,I)=DLAT(LIJ(I,J))
            lon(J,I)=DLON(LIJ(I,J))
            hz(J,I)=BELV(LIJ(I,J))
         ELSE
            lat(J,I)=-7999.
            lon(J,I)=-7999.
			hz(J,I)=-7999.
         ENDIF
      ENDDO
   ENDDO

    ! Create NetCDF file and overwrite if exists
    status=nf90_create(FILE_NAME, NF90_CLOBBER, ncid)
    if(status /= nf90_NoErr) call handle_err(status)

    !*************************************************************************************
    ! Variable definitions

    ! Define global attributes
    status=nf90_put_att(ncid, NF90_GLOBAL, 'title', 'EFDC Output')
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
    !status=nf90_def_var(ncid,'X',nf90_float,(/J_dimid,I_dimid/),X_varid)
    !status=nf90_put_att(ncid, X_varid, 'long_name', 'X coordinate')
    !status=nf90_put_att(ncid, X_varid, 'units', 'm')
    !status=nf90_put_att(ncid, X_varid, 'coord_sys', 'UTM')
    !status=nf90_put_att(ncid, X_varid, 'fill_value', -7999)

    !status=nf90_def_var(ncid,'Y',nf90_float,(/J_dimid,I_dimid/),Y_varid)
    !status=nf90_put_att(ncid, Y_varid, 'long_name', 'Y coordinate')
    !status=nf90_put_att(ncid, Y_varid, 'units', 'm')
    !status=nf90_put_att(ncid, Y_varid, 'coord_sys', 'UTM')
    !status=nf90_put_att(ncid, Y_varid, 'fill_value', -7999)
	
    ! Define depth
    status=nf90_def_var(ncid,'depth',nf90_real,(/J_dimid,I_dimid/),depth_varid)
    status=nf90_put_att(ncid, depth_varid, 'long_name', 'Bottom Elevation')
    status=nf90_put_att(ncid, depth_varid, 'caxis_label', 'Bottom Elevation (m)')
    status=nf90_put_att(ncid, depth_varid, 'units', 'm')
	
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
	
	! Define max. velocity
    status=nf90_def_var(ncid,'vel_max',nf90_float,(/J_dimid,I_dimid/),vmax_varid)
    status=nf90_put_att(ncid, vmax_varid, 'long_name', 'Max. Velocity')
	status=nf90_put_att(ncid, vmax_varid, 'caxis_label', 'Max. Velocity (m/s)')
    status=nf90_put_att(ncid, vmax_varid, 'units', 'm/s')

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

    ! Define shear stress
    status=nf90_def_var(ncid,'tau',nf90_float,(/J_dimid,I_dimid, time_dimid/),tau_varid)
    status=nf90_put_att(ncid, tau_varid, 'long_name', 'Instant Shear Stress')
    status=nf90_put_att(ncid, tau_varid, 'caxis_label', 'Shear Stress (dynes/cm^2)')
    status=nf90_put_att(ncid, tau_varid, 'units', 'dynes/cm^2')
	
	! Define maximum shear stress
    status=nf90_def_var(ncid,'taumax',nf90_float,(/J_dimid,I_dimid/),taumax_varid)
    status=nf90_put_att(ncid, taumax_varid, 'long_name', 'Maximum Shear Stress')
    status=nf90_put_att(ncid, taumax_varid, 'caxis_label', 'Max. Shear Stress (dynes/cm^2)')
    status=nf90_put_att(ncid, taumax_varid, 'units', 'dynes/cm^2')

    ! Define D50 grain size
    status=nf90_def_var(ncid,'D50',nf90_float,(/J_dimid,I_dimid, time_dimid/),d50_varid)
    status=nf90_put_att(ncid, d50_varid, 'long_name', 'Median Grain Size')
    status=nf90_put_att(ncid, d50_varid, 'caxis_label', 'Median Grain Size (um)')
    status=nf90_put_att(ncid, d50_varid, 'units', 'um')

    ! Define sediment thickness
    status=nf90_def_var(ncid,'thickness',nf90_float,(/J_dimid,I_dimid, time_dimid/),thick_varid)
    status=nf90_put_att(ncid, thick_varid, 'long_name', 'Bed Thickness')
    status=nf90_put_att(ncid, thick_varid, 'caxis_label', 'Bed Thickness (cm)')
    status=nf90_put_att(ncid, thick_varid, 'units', 'cm')

    ! End define mode
    status=nf90_enddef(ncid)

    !*************************************************************************************
    ! Put first step variables into file

    ! Put EFDC time step into file once
    status=nf90_put_var(ncid, ts_varid, deltat)
    if(status /= nf90_NoErr) call handle_err(status)

    ! Put coordinates into file
    !status=nf90_put_var(ncid,X_varid,lon)
    !if(status /= nf90_NoErr) call handle_err(status)
    !status=nf90_put_var(ncid,Y_varid,lat)
    !if(status /= nf90_NoErr) call handle_err(status)
	
    ! Put depth
    status=nf90_put_var(ncid,depth_varid,hz)
    if(status /= nf90_NoErr) call handle_err(status)

    FIRST_NETCDF=.TRUE.
	
	! Initialize variables once
	vel_magc=0.0
	vel_max=0.0
	TAUMAX=0.0

ENDIF

DO L=2,LA
    ! Calculate surface elevations for all active cells
    zeta(L)=(HP(L)+BELV(L))

    ! Create wet dry mask
    IF (HP(L)<=HDRY) THEN
        wet_dry_mask(L)=-99.
    ELSE
        wet_dry_mask(L)=1.
    ENDIF
	
	! Calculate sediment thickness
	TSET0T(L)=SUM(TSED0(1:KB,L)/BULKDENS(1:KB,L))
	TSEDT(L)=SUM(TSED(1:KB,L)/BULKDENS(1:KB,L))
	THCK(L)=TSEDT(L)-TSET0T(L)
ENDDO

! Create 2D arrays for all time variable parameters
DO J=3,JC-2
   DO I=3,IC-2
      IF(LIJ(I,J)>0) THEN
		L=LIJ(I,J)
		
		! Initialize temp. velocity variables
		utmpa=0.0
		vtmpa=0.0
		
		IF(LMASKDRY(L).AND.HP(L).GT.0.3) THEN
		
			! Calculate temp. velocities for each water layer
			DO K=1,KC
				utmps=U(L,K) ! m/s
				vtmps=V(L,K) 
				utmpa=utmps(L)*DZC(K)+utmpa
				vtmpa=vtmps(L)*DZC(K)+vtmpa
			ENDDO
		 
			! Cell max. velocity
			vel_magc(L)=SQRT(utmpa**2+vtmpa**2)
		 
			IF(vel_magc(L).GT.vel_max(L)) THEN
				vel_max(L)=vel_magc(L)
			ENDIF

			IF(TAU(L).GT.TAUMAX(L)) THEN
				TAUMAX(L)=TAU(L)
			ENDIF
		ENDIF
		
        mask(J,I)=wet_dry_mask(LIJ(I,J))
        wl(J,I)=zeta(LIJ(I,J))
        u_2d(J,I)=U(LIJ(I,J),1)
        v_2d(J,I)=V(LIJ(I,J),1)
		mag(J,I)=vel_max(LIJ(I,J))
        salt(J,I)=SAL(LIJ(I,J),1)
        dye_2d(J,I)=DYE(LIJ(I,J),1)
        shear(J,I)=TAU(LIJ(I,J))
		maxshear(J,I)=TAUMAX(LIJ(I,J))
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
		 mag(J,I)=-7999.
         salt(J,I)=-7999.
         dye_2d(J,I)=-7999.
         shear(J,I)=-7999.
		 maxshear(J,I)=-7999.
         grain_size(J,I)=-7999.
         sed_thick(J,I)=-7999.
         DO S=1,NSCM
            tss(J,I,S)=-7999.
         ENDDO
      ENDIF
   ENDDO
ENDDO

! Open NetCDF file
status=nf90_open(FILE_NAME, nf90_write, ncid)
if(status /= nf90_NoErr) call handle_err(status)

!*************************************************************************************
! Put time stepped variables into file

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

!Put max. velocity into file
status=nf90_put_var(ncid,vmax_varid,mag,start=(/1,1,1/))
if(status /= nf90_NoErr) call handle_err(status)

! Put salinity into file
status=nf90_put_var(ncid,sal_varid,salt,start=start_3d)
if(status /= nf90_NoErr) call handle_err(status)

! Put dye into file
status=nf90_put_var(ncid,dye_varid,dye_2d,start=start_3d)
if(status /= nf90_NoErr) call handle_err(status)

! Put TSS into file
status=nf90_put_var(ncid,tss_varid,tss,start=start_4d)
if(status /= nf90_NoErr) call handle_err(status)

! Put shear stress into file
status=nf90_put_var(ncid,tau_varid,shear,start=start_3d)
if(status /= nf90_NoErr) call handle_err(status)

! Put max. shear stress into file
status=nf90_put_var(ncid,taumax_varid,maxshear,start=(/1,1,1/))
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


CONTAINS
    SUBROUTINE handle_err(status)
        INTEGER, INTENT ( in) :: status

        IF(status /= nf90_noerr) THEN
            PRINT *, TRIM(nf90_strerror(status))
            STOP 'Stopped'
        END IF
    END SUBROUTINE handle_err
END SUBROUTINE netcdf_write
