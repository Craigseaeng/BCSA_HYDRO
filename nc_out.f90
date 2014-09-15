SUBROUTINE nc_out

! Writes EFDC output to NetCDF file


! 7/3/2014 Chris Flanary

USE GLOBAL
USE netcdf
IMPLICIT NONE

LOGICAL,SAVE::FIRST_NETCDF=.FALSE.
CHARACTER (len = *), PARAMETER :: FILE_NAME = "efdc_his.nc"
INTEGER, SAVE :: nc_step, ncid
INTEGER :: L, status
INTEGER, SAVE :: time_dimid,lcm_dimid,lc_dimid,k_dimid
INTEGER, SAVE :: ts_varid,time_varid
INTEGER, SAVE :: zeta_varid,u_varid,v_varid,sal_varid, dye_varid
INTEGER, SAVE :: mask_varid, tau_varid,d50_varid,thick_varid,tss_varid

INTEGER :: start(1),start_2d(2),start_3d(3),count_2d_lcm(2),count_2d_lc(2),count_3d_lcm(3)
REAL, DIMENSION(1) :: deltat,time_efdc
REAL, DIMENSION(LCM) :: zeta, wet_dry_mask

! EFDC time parameter
time_efdc=DT*FLOAT(N)+TCON*TBEGIN 
time_efdc=time_efdc/86400.

! Timing parameters
deltat=tidalp/float(ntsptc)
nc_step=nc_step+1

! Current write step
start=(/nc_step/)
start_2d=(/1,nc_step/)
start_3d=(/1,1,nc_step/)
count_2d_lcm=(/LCM,1/)
count_2d_lc=(/LC,1/)
count_3d_lcm=(/LCM,NSCM,1/)

IF(.NOT.FIRST_NETCDF)THEN

    ! Create NetCDF file and overwrite if exists
    status=nf90_create(FILE_NAME, NF90_CLOBBER, ncid)
    if(status /= nf90_NoErr) call handle_err(status)

    !*************************************************************************************
    ! Variable definitions

    ! Define global attributes
    status=nf90_put_att(ncid, NF90_GLOBAL, 'title', 'BCSA EFDC Simulation')
    status=nf90_put_att(ncid, NF90_GLOBAL, 'file', FILE_NAME)
    status=nf90_put_att(ncid, NF90_GLOBAL, 'format', 'netCDF-3 64bit offset file')
    status=nf90_put_att(ncid, NF90_GLOBAL, 'os', 'Linux')
    status=nf90_put_att(ncid, NF90_GLOBAL, 'arch', 'x86_64')

    ! Define deltat
    status=nf90_def_var(ncid,'timestep',nf90_real,ts_varid)
    status=nf90_put_att(ncid, ts_varid, 'long_name', 'EFDC numerical timestep')
    status=nf90_put_att(ncid, ts_varid, 'units', 'seconds')

    ! Define dimensions
    status=nf90_def_dim(ncid,'LCM',LCM,lcm_dimid)
    status=nf90_def_dim(ncid,'LC',LC,lc_dimid)
    status=nf90_def_dim(ncid,'NSCM',NSCM,k_dimid)
    status=nf90_def_dim(ncid,'efdc_time',nf90_unlimited,time_dimid)

    ! Define model time
    status=nf90_def_var(ncid,'efdc_time',nf90_real,time_dimid,time_varid)
    status=nf90_put_att(ncid, time_varid, 'long_name', 'EFDC time since initialization')
    status=nf90_put_att(ncid, time_varid, 'units', 'seconds since ???')
    status=nf90_put_att(ncid, time_varid, 'calendar', 'gregorian')

    ! Define water surface elevation
    status=nf90_def_var(ncid,'zeta',nf90_real,(/lcm_dimid,time_dimid/),zeta_varid)
    status=nf90_put_att(ncid, zeta_varid, 'long_name', 'EFDC surface elevation')
    status=nf90_put_att(ncid, zeta_varid, 'units', 'meters')

    ! Define u component velocity
    status=nf90_def_var(ncid,'u',nf90_float,(/lcm_dimid,time_dimid/),u_varid)
    status=nf90_put_att(ncid, u_varid, 'long_name', 'EFDC u component velocity')
    status=nf90_put_att(ncid, u_varid, 'units', 'meters second-1')

    ! Define v component velocity
    status=nf90_def_var(ncid,'v',nf90_float,(/lcm_dimid,time_dimid/),v_varid)
    status=nf90_put_att(ncid, v_varid, 'long_name', 'EFDC v component velocity')
    status=nf90_put_att(ncid, v_varid, 'units', 'meters second-1')

    ! Define salinity
    status=nf90_def_var(ncid,'salt',nf90_float,(/lcm_dimid,time_dimid/),sal_varid)
    status=nf90_put_att(ncid, sal_varid, 'long_name', 'EFDC salinity')
    status=nf90_put_att(ncid, sal_varid, 'units', 'psu')

    ! Define dye concentration
    status=nf90_def_var(ncid,'dye',nf90_float,(/lcm_dimid,time_dimid/),dye_varid)
    status=nf90_put_att(ncid, dye_varid, 'long_name', 'EFDC dye concentration')
    status=nf90_put_att(ncid, dye_varid, 'units', 'mg/L')

    ! Define water column TSS
    status=nf90_def_var(ncid,'tss',nf90_float,(/lcm_dimid,k_dimid,time_dimid/),tss_varid)
    status=nf90_put_att(ncid, tss_varid, 'long_name', 'EFDC water column TSS')
    status=nf90_put_att(ncid, tss_varid, 'units', 'mg/L') 

    ! Define wet dry mask
    status=nf90_def_var(ncid,'wet_dry',nf90_int,(/lcm_dimid, time_dimid/),mask_varid)
    status=nf90_put_att(ncid, mask_varid, 'long_name', 'EFDC wet dry mask')
    status=nf90_put_att(ncid, mask_varid, 'units', '0 = dry, 1 = wet')

    ! Define bed shear
    status=nf90_def_var(ncid,'tau',nf90_float,(/lc_dimid,time_dimid/),tau_varid)
    status=nf90_put_att(ncid, tau_varid, 'long_name', 'EFDC bed shear')
    status=nf90_put_att(ncid, tau_varid, 'units', 'dynes/cm^2')

    ! Define D50 grain size
    status=nf90_def_var(ncid,'D50',nf90_float,(/lc_dimid,time_dimid/),d50_varid)
    status=nf90_put_att(ncid, d50_varid, 'long_name', 'EFDC D50 grain size')
    status=nf90_put_att(ncid, d50_varid, 'units', 'um')

    ! Define sediment thickness
    status=nf90_def_var(ncid,'thickness',nf90_float,(/lc_dimid,time_dimid/),thick_varid)
    status=nf90_put_att(ncid, thick_varid, 'long_name', 'EFDC sediment thickness')
    status=nf90_put_att(ncid, thick_varid, 'units', 'meters')   
    
    ! End define mode
    status=nf90_enddef(ncid)

    !*************************************************************************************
    ! Put first step variables into file

    ! Wet dry mask
    DO L=2,LA
       ! Calculate surface elevations for all active cells
       zeta(L)=(HP(L)+BELV(L))

       IF (HP(L)<=HDRY) THEN
          wet_dry_mask(L)=0
       ELSE
          wet_dry_mask(L)=1
       ENDIF
    ENDDO

    ! Put EFDC time step into file once
    status=nf90_put_var(ncid, ts_varid, deltat)
    if(status /= nf90_NoErr) call handle_err(status) 

    ! Put first mask into file
    status=nf90_put_var(ncid, mask_varid, wet_dry_mask, start=start_2d, count=count_2d_lcm)
    if(status /= nf90_NoErr) call handle_err(status)

    ! Put first EFDC time in days into file
    status=nf90_put_var(ncid, time_varid, time_efdc, start=start)
    if(status /= nf90_NoErr) call handle_err(status)

    ! Put first surface elevation into file
    status=nf90_put_var(ncid,zeta_varid,zeta,start=start_2d, count=count_2d_lcm)
    if(status /= nf90_NoErr) call handle_err(status)

    ! Put first u component into file
    status=nf90_put_var(ncid,u_varid,U,start=start_2d, count=count_2d_lcm)
    if(status /= nf90_NoErr) call handle_err(status)
    
    ! Put first v component into file
    status=nf90_put_var(ncid,v_varid,V,start=start_2d, count=count_2d_lcm)
    if(status /= nf90_NoErr) call handle_err(status)

    ! Put first salinity into file
    status=nf90_put_var(ncid,sal_varid,SAL,start=start_2d, count=count_2d_lcm)
    if(status /= nf90_NoErr) call handle_err(status)

    ! Put first dye into file
    status=nf90_put_var(ncid,dye_varid,DYE,start=start_2d, count=count_2d_lcm)
    if(status /= nf90_NoErr) call handle_err(status)

    ! Put first TSS into file
    status=nf90_put_var(ncid,tss_varid,SED(:,1,:),start=start_3d, count=count_3d_lcm)
    if(status /= nf90_NoErr) call handle_err(status)

    ! Put first tau into file
    status=nf90_put_var(ncid,tau_varid,TAU,start=start_2d, count=count_2d_lc)
    if(status /= nf90_NoErr) call handle_err(status)

    ! Put first D50 into file
    status=nf90_put_var(ncid,d50_varid,D50AVG,start=start_2d, count=count_2d_lc)
    if(status /= nf90_NoErr) call handle_err(status)

    ! Put first sediment thickness into file
    status=nf90_put_var(ncid,thick_varid,THCK,start=start_2d, count=count_2d_lc)
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

    ! Wet dry mask
    DO L=2,LA
       ! Calculate surface elevations for all active cells
       zeta(L)=(HP(L)+BELV(L))

       IF (HP(L)<=HDRY) THEN
          wet_dry_mask(L)=0
       ELSE
          wet_dry_mask(L)=1
       ENDIF
    ENDDO

    ! Put EFDC times into file
    status=nf90_put_var(ncid, time_varid, time_efdc, start=start)
    if(status /= nf90_NoErr) call handle_err(status)

    ! Put first mask into file
    status=nf90_put_var(ncid, mask_varid, wet_dry_mask, start=start_2d, count=count_2d_lcm)
    if(status /= nf90_NoErr) call handle_err(status)

    ! Put surface elevations into file
    status=nf90_put_var(ncid,zeta_varid,zeta,start=start_2d, count=count_2d_lcm)
    if(status /= nf90_NoErr) call handle_err(status)

    ! Put u components into file
    status=nf90_put_var(ncid,u_varid,U,start=start_2d, count=count_2d_lcm)
    if(status /= nf90_NoErr) call handle_err(status)
    
    ! Put v components into file
    status=nf90_put_var(ncid,v_varid,V,start=start_2d, count=count_2d_lcm)
    if(status /= nf90_NoErr) call handle_err(status)

    ! Put salinities into file
    status=nf90_put_var(ncid,sal_varid,SAL,start=start_2d, count=count_2d_lcm)
    if(status /= nf90_NoErr) call handle_err(status)

    ! Put dye into file
    status=nf90_put_var(ncid,dye_varid,DYE,start=start_2d, count=count_2d_lcm)
    if(status /= nf90_NoErr) call handle_err(status)

    ! Put TSS into file
    status=nf90_put_var(ncid,tss_varid,SED(:,1,:),start=start_3d, count=count_3d_lcm)
    if(status /= nf90_NoErr) call handle_err(status)

    ! Put tau into file
    status=nf90_put_var(ncid,tau_varid,TAU,start=start_2d, count=count_2d_lc)
    if(status /= nf90_NoErr) call handle_err(status)

    ! Put D50 into file
    status=nf90_put_var(ncid,d50_varid,D50AVG,start=start_2d, count=count_2d_lc)
    if(status /= nf90_NoErr) call handle_err(status)

    ! Put sediment thickness into file
    status=nf90_put_var(ncid,thick_varid,THCK,start=start_2d, count=count_2d_lc)
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
END SUBROUTINE nc_out


