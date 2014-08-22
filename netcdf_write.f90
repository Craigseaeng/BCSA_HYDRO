SUBROUTINE netcdf_write(nstep,deltat,time,surfel)

! Writes EFDC output to NetCDF file

! Set for 3D variables only

! 7/3/2014 Chris Flanary

USE GLOBAL
USE netcdf
IMPLICIT NONE

LOGICAL,SAVE::FIRST_NETCDF=.FALSE.
CHARACTER (len = *), PARAMETER :: FILE_NAME = "efdc_his.nc"
INTEGER, SAVE :: ncid
INTEGER :: I, J, status
INTEGER, SAVE :: I_dimid,J_dimid,time_dimid
INTEGER, SAVE :: ts_varid,time_varid,X_varid,Y_varid
INTEGER, SAVE :: surfel_varid,u_varid,v_varid,sal_varid, dye_varid
INTEGER, SAVE :: tau_varid,d50_varid,thick_varid

INTEGER :: nstep, start(1), start_3d(3), count_3d(3)
REAL, DIMENSION(1) :: deltat, time
REAL, DIMENSION(LA) :: SURFEL

! Allocate 2D array variables
REAL, ALLOCATABLE, DIMENSION(:,:) :: lat,lon,wl,u_2d,v_2d,salt,dye_2d,bed_shear
REAL, ALLOCATABLE, DIMENSION(:,:) :: grain_size,sed_thick
ALLOCATE(lat(IC-2,JC-2))
ALLOCATE(lon(IC-2,JC-2))
ALLOCATE(wl(IC-2,JC-2))
ALLOCATE(u_2d(IC-2,JC-2))
ALLOCATE(v_2d(IC-2,JC-2))
ALLOCATE(salt(IC-2,JC-2))
ALLOCATE(dye_2d(IC-2,JC-2))
ALLOCATE(bed_shear(IC-2,JC-2))
ALLOCATE(grain_size(IC-2,JC-2))
ALLOCATE(sed_thick(IC-2,JC-2))
!ALLOCATE(dye_2d(IC-2,JC-2))

! Current write step
start=(/nstep/)
start_3d=(/1,1,nstep/)
count_3d=(/IC-4,JC-4,1/)

! Create 2D arrays for all time variable parameters
DO J=3,JC-2
   DO I=3,IC-2
      IF(LIJ(I,J)>0) THEN
         wl(I,J)=surfel(LIJ(I,J))
         u_2d(I,J)=U(LIJ(I,J),1)
         v_2d(I,J)=V(LIJ(I,J),1)
         salt(I,J)=SAL(LIJ(I,J),1)
         dye_2d(I,J)=DYE(LIJ(I,J),1)
         bed_shear(I,J)=TAU(LIJ(I,J))
         grain_size(I,J)=D50AVG(LIJ(I,J))
         sed_thick(I,J)=THCK(LIJ(I,J))
      ELSE
         ! Flag for inactive cells
         wl(I,J)=-7999
         u_2d(I,J)=-7999
         v_2d(I,J)=-7999
         salt(I,J)=-7999
         dye_2d(I,J)=-7999
         bed_shear(I,J)=-7999
         grain_size(I,J)=-7999
         sed_thick(I,J)=-7999
      ENDIF
   ENDDO
ENDDO

IF(.NOT.FIRST_NETCDF)THEN

   ! Create 2D lat lon array for all cells
   DO J=3,JC-2
      DO I=3,IC-2
         IF(LIJ(I,J)>0) THEN
            lat(I,J)=DLAT(LIJ(I,J))
            lon(I,J)=DLON(LIJ(I,J))
         ELSE
            lat(I,J)=-7999
            lon(I,J)=-7999
         ENDIF
      ENDDO
   ENDDO

    ! Create NetCDF file and overwrite if exists
    status=nf90_create(FILE_NAME, NF90_CLOBBER, ncid)

    !*************************************************************************************
    ! Variable definitions

    ! Define global attributes
    status=nf90_put_att(ncid, NF90_GLOBAL, 'title', PROJTITLE)
    status=nf90_put_att(ncid, NF90_GLOBAL, 'file', FILE_NAME)
    status=nf90_put_att(ncid, NF90_GLOBAL, 'format', 'netCDF-3 64bit offset file')
    status=nf90_put_att(ncid, NF90_GLOBAL, 'os', 'Linux')
    status=nf90_put_att(ncid, NF90_GLOBAL, 'arch', 'x86_64')

    ! Define deltat
    status=nf90_def_var(ncid,'timestep',nf90_real,ts_varid)
    status=nf90_put_att(ncid, ts_varid, 'long_name', 'EFDC numerical timestep')
    status=nf90_put_att(ncid, ts_varid, 'units', 'seconds')

    ! Define dimensions
    status=nf90_def_dim(ncid,'I',IC-2,I_dimid)
    status=nf90_def_dim(ncid,'J',JC-2,J_dimid)
    status=nf90_def_dim(ncid,'efdc_time',nf90_unlimited,time_dimid)

    ! Define model time
    status=nf90_def_var(ncid,'efdc_time',nf90_real,time_dimid,time_varid)
    status=nf90_put_att(ncid, time_varid, 'long_name', 'EFDC time since initialization')
    status=nf90_put_att(ncid, time_varid, 'units', 'seconds since ???')
    status=nf90_put_att(ncid, time_varid, 'calendar', 'gregorian')

    ! Define coordinates variables
    status=nf90_def_var(ncid,'X',nf90_float,(/I_dimid,J_dimid/),X_varid)
    status=nf90_put_att(ncid, X_varid, 'long_name', 'X coordinate')
    status=nf90_put_att(ncid, X_varid, 'units', 'meters')
    status=nf90_put_att(ncid, X_varid, 'fill_value', -7999)

    status=nf90_def_var(ncid,'Y',nf90_float,(/I_dimid,J_dimid/),Y_varid)
    status=nf90_put_att(ncid, Y_varid, 'long_name', 'Y coordinate')
    status=nf90_put_att(ncid, Y_varid, 'units', 'meters')
    status=nf90_put_att(ncid, Y_varid, 'fill_value', -7999)

    ! Define water surface elevation
    status=nf90_def_var(ncid,'surfel',nf90_real,(/I_dimid, J_dimid, time_dimid/),surfel_varid)
    status=nf90_put_att(ncid, surfel_varid, 'long_name', 'EFDC surface elevation')
    status=nf90_put_att(ncid, surfel_varid, 'units', 'meters')

    ! Define u component velocity
    status=nf90_def_var(ncid,'u',nf90_float,(/I_dimid, J_dimid, time_dimid/),u_varid)
    status=nf90_put_att(ncid, u_varid, 'long_name', 'EFDC u component velocity')
    status=nf90_put_att(ncid, u_varid, 'units', 'meters second-1')

    ! Define v component velocity
    status=nf90_def_var(ncid,'v',nf90_float,(/I_dimid, J_dimid, time_dimid/),v_varid)
    status=nf90_put_att(ncid, v_varid, 'long_name', 'EFDC v component velocity')
    status=nf90_put_att(ncid, v_varid, 'units', 'meters second-1')

    ! Define salinity
    status=nf90_def_var(ncid,'salt',nf90_float,(/I_dimid, J_dimid, time_dimid/),sal_varid)
    status=nf90_put_att(ncid, sal_varid, 'long_name', 'EFDC salinity')
    status=nf90_put_att(ncid, sal_varid, 'units', 'psu')

    ! Define dye concentration
    status=nf90_def_var(ncid,'dye',nf90_float,(/I_dimid, J_dimid, time_dimid/),dye_varid)
    status=nf90_put_att(ncid, dye_varid, 'long_name', 'EFDC dye concentration')
    status=nf90_put_att(ncid, dye_varid, 'units', 'mg/L')

    ! Define bed shear
    status=nf90_def_var(ncid,'tau',nf90_float,(/I_dimid, J_dimid, time_dimid/),tau_varid)
    status=nf90_put_att(ncid, tau_varid, 'long_name', 'EFDC bed shear')
    status=nf90_put_att(ncid, tau_varid, 'units', 'dynes/cm^2')

    ! Define D50 grain size
    status=nf90_def_var(ncid,'D50',nf90_float,(/I_dimid, J_dimid, time_dimid/),d50_varid)
    status=nf90_put_att(ncid, d50_varid, 'long_name', 'EFDC D50 grain size')
    status=nf90_put_att(ncid, d50_varid, 'units', 'um')

    ! Define sediment thickness
    status=nf90_def_var(ncid,'thickness',nf90_float,(/I_dimid, J_dimid, time_dimid/),thick_varid)
    status=nf90_put_att(ncid, thick_varid, 'long_name', 'EFDC sediment thickness')
    status=nf90_put_att(ncid, thick_varid, 'units', 'meters')
    
    ! End define mode
    status=nf90_enddef(ncid)

    !*************************************************************************************
    ! Put first step variables into file

    ! Put EFDC time step into file once
    status=nf90_put_var(ncid, ts_varid, deltat)

    ! Put first EFDC time in days into file
    status=nf90_put_var(ncid, time_varid, time, start=start)

    ! Put coordinates into file
    status=nf90_put_var(ncid,X_varid,lon)
    status=nf90_put_var(ncid,Y_varid,lat)

    ! Put first surface elevation into file
    status=nf90_put_var(ncid,surfel_varid,wl,start=start_3d, count=count_3d)

    ! Put first u component into file
    status=nf90_put_var(ncid,u_varid,u_2d,start=start_3d, count=count_3d)
    
    ! Put first v component into file
    status=nf90_put_var(ncid,v_varid,v_2d,start=start_3d, count=count_3d)

    ! Put first salinity into file
    status=nf90_put_var(ncid,sal_varid,salt,start=start_3d, count=count_3d)

    ! Put first dye into file
    status=nf90_put_var(ncid,dye_varid,dye_2d,start=start_3d, count=count_3d)

    ! Put first tau into file
    status=nf90_put_var(ncid,tau_varid,bed_shear,start=start_3d, count=count_3d)

    ! Put first D50 into file
    status=nf90_put_var(ncid,d50_varid,grain_size,start=start_3d, count=count_3d)

    ! Put first sediment thickness into file
    status=nf90_put_var(ncid,thick_varid,sed_thick,start=start_3d, count=count_3d)

    FIRST_NETCDF=.TRUE.
    
    ! Close NetCDF file after first write
    status=nf90_close(ncid)

ELSE

   !*************************************************************************************
    ! Put subsequent step variables into file

    ! Open NetCDF file
    status=nf90_open(FILE_NAME, nf90_write, ncid)

    ! Put EFDC times into file
    status=nf90_put_var(ncid, time_varid, time, start=start)

    ! Put surface elevations into file
    status=nf90_put_var(ncid,surfel_varid,wl,start=start_3d, count=count_3d)

    ! Put u components into file
    status=nf90_put_var(ncid,u_varid,u_2d,start=start_3d, count=count_3d)
    
    ! Put v components into file
    status=nf90_put_var(ncid,v_varid,v_2d,start=start_3d, count=count_3d)

    ! Put salinities into file
    status=nf90_put_var(ncid,sal_varid,salt,start=start_3d, count=count_3d)

    ! Put dye into file
    status=nf90_put_var(ncid,dye_varid,dye_2d,start=start_3d, count=count_3d)

    ! Put tau into file
    status=nf90_put_var(ncid,tau_varid,bed_shear,start=start_3d, count=count_3d)

    ! Put D50 into file
    status=nf90_put_var(ncid,d50_varid,grain_size,start=start_3d, count=count_3d)

    ! Put sediment thickness into file
    status=nf90_put_var(ncid,thick_varid,sed_thick,start=start_3d, count=count_3d)

    ! Close NetCDF file each time
    status=nf90_close(ncid)

ENDIF

!DEALLOCATE(lat(IC-2,JC-2))
!DEALLOCATE(lon(IC-2,JC-2))
!DEALLOCATE(wl(IC-2,JC-2))
!DEALLOCATE(u_2d(IC-2,JC-2))
!DEALLOCATE(v_2d(IC-2,JC-2))
!DEALLOCATE(salt(IC-2,JC-2))
!DEALLOCATE(dye_2d(IC-2,JC-2))
!DEALLOCATE(bed_shear(IC-2,JC-2))
!DEALLOCATE(grain_size(IC-2,JC-2))
!DEALLOCATE(sed_thick(IC-2,JC-2))

RETURN
END SUBROUTINE netcdf_write
