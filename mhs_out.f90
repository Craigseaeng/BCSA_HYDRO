SUBROUTINE MHS_OUT

! Output MHS station data in calibration files
!***************************************************************
!***************************************************************

USE GLOBAL
IMPLICIT NONE

INTEGER::I,J,L,K,LOC
INTEGER::ITEMP1,ITEMP2,JTEMP1,JTEMP2
REAL,DIMENSION(LCM)::zeta,vel_maxc,vel_max,tau_max
REAL::time_efdc, tss_flux_u_tmp, tss_flux_v_tmp
REAL,DIMENSION(8)::tss_flux_u, tss_flux_v
LOGICAL,SAVE::FIRSTTIME=.FALSE.	

! Create files if first call
IF(.NOT.FIRSTTIME)THEN

    ! Velocity calibration file with surface elevation and velocity components at all MHS stations
    OPEN (UNIT=112,FILE='vel_cal.dat', STATUS='REPLACE')

    ! Tracer calibration file with salinity and dye
    IF (ISTRAN(1).EQ.1.OR.ISTRAN(3).EQ.1) OPEN (UNIT=115,FILE='tracer_cal.dat', STATUS='REPLACE')

    IF (ISTRAN(6).EQ.1) THEN
        ! TSS calibration file
        OPEN (UNIT=105, FILE='tss_cal.dat', STATUS='REPLACE')
    
        ! Bed thickness calibration file
        OPEN (UNIT=106, FILE='thick_cal.dat', STATUS='REPLACE')
    
        ! TSS flux file
        OPEN (UNIT=108, FILE='tss_flux.dat', STATUS='REPLACE')
    ENDIF

    ! Shear stress calibration file
    OPEN (UNIT=107, FILE='shear_cal.dat', STATUS='REPLACE')

    ! Initialize variable for first time step
    vel_maxc=0.0
    vel_max=0.0
    tau_max=0.0

    FIRSTTIME=.TRUE.
ENDIF

! EFDC time parameter
time_efdc=DT*FLOAT(N)+TCON*TBEGIN 
time_efdc=time_efdc/86400.

IF (ISTRAN(6).EQ.1) THEN
! Calculate TSS flux at MHS01, MHS06, and MHS07
DO LOC=1,4
    tss_flux_u(LOC)=0.0
    tss_flux_v(LOC)=0.0

    SELECT CASE (LOC)
    CASE(1)
    ITEMP1=116;ITEMP2=129;JTEMP1=114;JTEMP2=114

    CASE(2)
    ITEMP1=35;ITEMP2=41;JTEMP1=203;JTEMP2=203

    CASE(3)
    ITEMP1=119;ITEMP2=119;JTEMP1=311;JTEMP2=314

    CASE(4)
    ITEMP1=128;ITEMP2=131;JTEMP1=349;JTEMP2=349
    END SELECT

    DO I=ITEMP1,ITEMP2
        DO J=JTEMP1,JTEMP2
            IF(LMASKDRY(LIJ(I,J))) THEN
                DO K=1,NSCM
                    tss_flux_u_tmp=U(LIJ(I,J),1)*HP(LIJ(I,J))*DXU(LIJ(I,J))*SED(LIJ(I,J),1,K)
                    tss_flux_u(LOC)=tss_flux_u(LOC)+tss_flux_u_tmp
                    tss_flux_v_tmp=V(LIJ(I,J),1)*HP(LIJ(I,J))*DYV(LIJ(I,J))*SED(LIJ(I,J),1,K)
                    tss_flux_v(LOC)=tss_flux_v(LOC)+tss_flux_v_tmp                    
                ENDDO
            ENDIF
        ENDDO
    ENDDO
ENDDO
ENDIF

DO  L=2,LA
    ! Sediment thickness
    TSET0T(L)=SUM(TSED0(1:KB,L)/BULKDENS(1:KB,L))
    TSEDT(L)=SUM(TSED(1:KB,L)/BULKDENS(1:KB,L))
    THCK(L)=TSEDT(L)-TSET0T(L)

    ! Calculate surface elevations and velocities for all active cells
    zeta(L)=(HP(L)+BELV(L))
ENDDO 

! Write velocity calibration data each call
WRITE(112,'(16F7.3)') time_efdc,zeta(LIJ(122,114)), &
    U(LIJ(122,114),1),V(LIJ(122,114),1),zeta(LIJ(45,29)),U(LIJ(45,29),1),V(LIJ(45,29),1), &
    zeta(LIJ(39,202)),U(LIJ(39,202),1),V(LIJ(39,202),1),zeta(LIJ(119,312)),U(LIJ(119,312),1), &
    V(LIJ(119,312),1),zeta(LIJ(130,349)),U(LIJ(130,349),1),V(LIJ(130,349),1)

IF (ISTRAN(1).EQ.1.OR.ISTRAN(3).EQ.1) THEN
    ! Write tracer calibration data each call
    WRITE(115,'(11F7.3)') time_efdc,SAL(LIJ(122,114),1), &
    DYE(LIJ(122,114),1),SAL(LIJ(45,29),1),DYE(LIJ(45,29),1),SAL(LIJ(39,202),1),DYE(LIJ(39,202),1), &
    SAL(LIJ(119,312),1),DYE(LIJ(119,312),1),SAL(LIJ(130,349),1),DYE(LIJ(130,349),1)
ENDIF

! If sediment is activated
IF(ISTRAN(6).EQ.1) THEN
    ! TSS calibration file with SEDZLJ (hard coded for 1 water layer (KC))
    DO K=1,NSCM
        WRITE(105,299)  time_efdc, K, SED(LIJ(122,114),1,K), &
        SED(LIJ(45,29),1,K), SED(LIJ(39,202),1,K), SED(LIJ(119,312),1,K), SED(LIJ(130,349),1,K)
    ENDDO

    WRITE(106,'(6F7.3)')  time_efdc, THCK(LIJ(122,114)), THCK(LIJ(45,29)), &
    THCK(LIJ(39,202)), THCK(LIJ(119,312)), THCK(LIJ(130,349))

    WRITE(108,'(9F12.3)')  time_efdc, tss_flux_u(1),tss_flux_v(1), tss_flux_u(2),tss_flux_v(2), &
    tss_flux_u(3),tss_flux_v(3), tss_flux_u(4),tss_flux_v(4)
ENDIF

DO J=3,JC-2
    DO I=3,IC-2
        IF(LIJ(I,J)>0) THEN
            L=LIJ(I,J)
            IF(LMASKDRY(L).AND.HP(L).GT.0.3) THEN

                IF(vel_maxc(L).GT.vel_max(L)) THEN
                    vel_max(L)=vel_maxc(L)
                ENDIF

                IF(TAU(L).GT.tau_max(L)) THEN
                    tau_max(L)=TAU(L)
                ENDIF

            ENDIF
        ENDIF
    ENDDO
ENDDO

! Shear stress calibration file
WRITE(107,'(11F10.3)') time_efdc,TAU(LIJ(122,114)),tau_max(LIJ(122,114)),TAU(LIJ(45,29)),tau_max(LIJ(45,29)), &
    TAU(LIJ(39,202)),tau_max(LIJ(39,202)),TAU(LIJ(119,312)),tau_max(LIJ(119,312)), &
    TAU(LIJ(130,349)),tau_max(LIJ(130,349))

! Format for tss_cal.dat file
299 FORMAT(F7.3,2X,I1,F12.3,F12.3,F12.3,F12.3,F12.3)

FLUSH(112)
IF (ISTRAN(1).EQ.1.OR.ISTRAN(3).EQ.1) FLUSH(115)

IF (ISTRAN(6).EQ.1) THEN
    FLUSH(105)
    FLUSH(106)
    FLUSH(108)
ENDIF

FLUSH(107)

RETURN
END SUBROUTINE MHS_OUT
