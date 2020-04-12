module ALL_MEGAN

      USE MEG_EXT      
      IMPLICIT NONE
      
CONTAINS    !!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

 ! FSB These subroutines are needed in MEGCAN
        
! ______________ SUBROUTINE get_CBETA_____________
! This subroutine calculates the solar zenith angle , its sine and the eccentricity
!    INPUT variables:
!    JD       current Julian Day.
!    LAT      current latitude [deg]
!    HOUR     current hour [ hr ] 

!    OUTPUT variables:
!    BETA            curent solar zenith angle
!    SINBETA         Sine of zenith angle
!    ECCEBTRICITY    orbital eccentricity (solar distance) [ AU ] 
!  coded by Dr. Francis S.Binkowski on April 2, 2019, modified Abril 4, 2019
!  based upon earlier code and updated to latest algorithms from USNO

!   Reference:
!   The algorithm for the solar position  are from 
!   https://aa.usno.navy.mil/faq/docs/SunApprox.php
!   Us NAVAL OBSERVATORY The acuracy is good for two centuries 
!   before and after 2000 CE. 


      SUBROUTINE get_BETA( JD, LAT, HOUR, BETA, SINBETA, ECCENTRICITY )
      IMPLICIT NONE
!*--CALCBETA575
! INPUTS: 
      REAL, INTENT(IN) :: JD   ! True Julian day of interest
      REAL, INTENT(IN) :: LAT  ! Latitude [ degrees ]
      REAL, INTENT(IN) :: HOUR ! local standard time
      
! OUTPUT:
      REAL, INTENT(OUT) :: BETA ! the zenith angle of the sun at lat, hour in degrees
      REAL, INTENT(OUT) :: SINBETA ! sine of BETA
      REAL, INTENT(OUT) :: ECCENTRICITY ! orbital eccentricity 
                                        ! earth distance from sun {AU}
! LOCAL VARIBLES            
      REAL :: sindelta , cosdelta , a , b ,  D, R, num, den, decl 
      REAL :: sinepsilon, cosepsilon, sinlamda, coslamda, RA, EQT, hangle
      REAL :: g, q, L, lamda, epsilon, e, sing, singsq, sin2g, cosg, cos2g, latrad      
      REAL, PARAMETER :: PI = 4.0 * ATAN(1.0)
      REAL, PARAMETER :: RAD2DEG = 180.0 / PI,  DEG2RAD = 1.0 / RAD2DEG
      REAL, PARAMETER :: JD0 = 2451545.5 ! True Julian day for January 1, 2000 at midnight
      REAL, PARAMETER :: ONE15 = 1.0  / 15.0        
      REAL, PARAMETER :: RADCONV = RAD2DEG * ONE15       

       
!--------------------------------------------------------------------

      D = JD  - JD0 
      
      g = 357.529 + 0.98560028 * D  ! Mean anomaly of the Sun

      g = MODULO(g,360.0)
      IF ( g.LT.0.0 ) g = g + 360.0
 
      g = deg2RAD * g   ! g in radians now

!     calculate trig functions of g using identities
!     this speeds up the calculations
 
      sing = SIN(g)
      singsq = sing * sing
      cosg  = sqrt( 1 - singsq)
      sin2g = 2.0 * sing * cosg
      cos2g = cosg*cosg - sing*sing

      q = 280.459 + 0.98564736 *  D  ! Mean longitude of the Sun:

! *** now force L to be betweeon 0.0 and 360. degrees
      q = MODULO( q,360.0)
      IF ( q.LT.0.0 ) q = q + 360.0

      lamda = q + 1.915 * sing + 0.020 * sin2g ! apparent longitude of the sun
      lamda = MODULO(lamda,360.0)
      IF ( lamda.LT.0.0 ) lamda = lamda + 360.0     
      
      epsilon = 23.439 - 3.6e-7 * D  ! obliquity of ecliptic

! FSB convert to radian
       epsilon   = DEG2RAD * epsilon
       lamda     = DEG2RAD * lamda
      
      sinepsilon = SIN(epsilon)  
      sinlamda   = SIN(lamda)
      coslamda   = SQRT( 1.0 - sinlamda *sinlamda )
            
      sindelta = sinepsilon * sinlamda           ! sine of solar declination
      cosdelta = SQRT(1.0 - sindelta * sindelta) ! cosine of solar declination

      cosepsilon = SQRT( 1.0 - sinepsilon * sinepsilon )      

! FSB calculate Right Ascension of the sun

       num = cosepsilon * sinlamda
       den = coslamda
       RA  = atan2(num,den) 

! FSB Ignore EQT  adds at most plus or minus 15 minutes over the year
!       EQT = ONE15 * ( q - RA)   ! Equation  of Time        

!      hangle = ( 15.0*(HOUR-12.0) + EQT ) * DEG2RAD ! hour angle in radians
      hangle = ( 15.0*(HOUR-12.0) ) * DEG2RAD ! hour angle in radians
           
      latrad = DEG2RAD * LAT
      
 
      a = SIN(latrad) * sindelta
      b = COS(Latrad) * cosdelta
      SINBETA = a + b * COS( hangle) 
      
      BETA = ASIN(SINBETA) * RAD2DEG ! [ degrees]
     
! FSB calculate solar distance [ Astronomical Units ] This does
!     change over the seasons and more importantly over
!     annual and longer time scales because it is a function of
!     the Mean anomaly of the Sun.
 
      R            = 1.00014 - 0.01671*cosg - 0.00014*cos2g 
      ECCENTRICITY = R    
      RETURN
      END SUBROUTINE  get_BETA
 


     REAL FUNCTION getJD (YEAR,MONTH,DAY)
!
!---COMPUTES THE JULIAN DATE (JD) GIVEN A GREGORIAN CALENDAR
!   DATE (YEAR,MONTH,DAY).

!   REFERENCE:
!   Reda,Ibrahim, and Andreas Afshin, 2008, Solar position algorithm for solar
!    radiation applications, NREL/TP-550-34302, Revised January 2008l, National 
!    Renewable Energy Laboratory, Golden CO. 
!    Coded April 10, 2019 by Dr. Francis S. Binkowski 
 
    REAL YEAR,MONTH,DAY,Y, M, D, A, B,  JD
!
    Y  = YEAR
    M  = MONTH
    D  = DAY

! FSB The following is from Equation (4) on Page 3 of the reference.
  
    A  = AINT( Y / 100.0)
    B = 2.0 -A +AINT(A/4) 
    
  getJD = AINT( 365.25 *( Y + 4716.0 ) ) + AINT( 30.6001*( M + 1.0) ) + D    &
            + B - 1524.5
!
    RETURN
    END FUNCTION getJD  





!///////////////_______----------------------- 
! FSB            This SUBROUTINE converts a date on the Gregorian
!                 calendar to a day of the year.
! REFERENCE:

! Original  Programmer:   David G. Simpson
!                NASA Goddard Space Flight Center
!                Greenbelt, Maryland  2077  Date:         November 20, 2001


!  Modified April 13, 2019 by Dr Francis S. Binkowski to do only Gregorian years


      SUBROUTINE get_DOY( YEAR, MONTH, DAY, DOY)

      IMPLICIT NONE
! FSB INPUT:      
      REAL, INTENT(IN)  :: YEAR, MONTH, DAY   ! GREGORIAN DATE                                                    

! FSB OUTPUT:
      INTEGER, INTENT(OUT)  :: DOY                ! DAY OF the YEAR
      
! FSB LOCAL:
      
      INTEGER :: Y                               ! year
      INTEGER :: M                               ! month (1-12)
      INTEGER :: D                               ! day of month 
      INTEGER :: K
   
      
      LOGICAL :: LEAP

! FSB BEGIN CODE:
     
       Y = YEAR
       M = MONTH
       D = DAY
       
      LEAP = .FALSE.
      
! FSB TEST FOR LEAP YEARS
      
      IF ( MOD(Y,4)   .EQ. 0) LEAP = .TRUE.
      IF ( MOD(Y,100) .EQ. 0) LEAP = .FALSE.
      IF  (MOD(Y,400) .EQ. 0) LEAP = .TRUE.

      IF (LEAP) THEN
         K = 1
      ELSE
         K = 2
      END IF

! FSB CALCULATE DAY OF THE YEAR using INTEGER arithmetic

      DOY = ( ( 275 * M) / 9 ) - K * ( ( M + 9) / 12 ) + D - 30

      RETURN
      
      END SUBROUTINE get_DOY
      
     SUBROUTINE get_date(YEAR,DOY,MM,DD)
     
!=============WHEN GIVEN A VALID YEAR, YYYY, AND DAY OF THE
!             YEAR, DDD, RETURNS THE MONTH, MM, AND DAY OF THE
!             MONTH, DD.

! REference:
!             SEE ACM ALGORITHM 398, TABLELESS DATE CONVERSION, BY
!            {*filter*} STONE, CACM 13(10):621.ACM 1970; DOI:10.1145/355598.362779. 7. 
!            Summary We have introduced a formalism which allows us to 
!            explicate certain rather gross properties of language ...

! FSB Modified to f90 by Dr. Francis S. Binkowski on April 13, 2019

! FSB    INPUT:
      REAL, INTENT(in)  ::  YEAR , DOY  ! year and day of year 
      
! FSB    OUTPUT:      
      INTEGER, INTENT(out) :: MM , DD

! FSB   LOCAL:       
      INTEGER              :: T
      INTEGER              ::  YYYY , DDD

! FSB Start code 

      YYYY = YEAR
      DDD  = DOY 
      T    = 0
      
      IF( MOD(YYYY,4) .EQ. 0)  T = 1 ! test for leap year
      
!-----------THE FOLLOWING STATEMENT IS NECESSARY IF YYYY IS LESS TNAN
!           1900 OR GREATER THAN 2100.

      IF( MOD(YYYY,400) .NE.0 .AND. MOD(YYYY,100) .EQ. 0 ) T = 0
      
      DD = DDD
      
      IF( DDD.GT. 59 + T ) DD = DD + 2 - T
      
      MM =  (  (DD + 91 ) * 100 ) / 3055
      
      DD = ( DD + 91 ) - ( MM * 3055 ) / 100
      
      MM = MM - 2
!----------MM WILL BE CORRECT IFF DDD IS CORRECT FOR YYYY.
      
      RETURN
      
     
      END  SUBROUTINE get_date
      
     
      
      
! FSB These subroutines are from MEGCAN
 
!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
!
!   SUBROUTINE GaussianDist
!
!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
 
      SUBROUTINE GAUSSIANDIST(Distgauss,Layers)
 
      IMPLICIT NONE
!*--GAUSSIANDIST91
 
      INTEGER , INTENT(IN) :: Layers
 
      REAL , DIMENSION(Layers) , INTENT(OUT) :: Distgauss
 
! local variables
      INTEGER :: i
!----------------------------------------------------------------
 
      IF ( Layers.EQ.1 ) THEN
         Distgauss(1) = 0.5
      ELSEIF ( Layers.EQ.3 ) THEN
         Distgauss(1) = 0.112702
         Distgauss(2) = 0.5
         Distgauss(3) = 0.887298
      ELSEIF ( Layers.EQ.5 ) THEN
         Distgauss(1) = 0.0469101
         Distgauss(2) = 0.2307534
         Distgauss(3) = 0.5
         Distgauss(4) = 0.7692465
         Distgauss(5) = 0.9530899
      ELSE
         DO i = 1 , Layers
            Distgauss(i) = (i-0.5)/Layers
         ENDDO
      ENDIF
      
      RETURN
      
      END  SUBROUTINE GAUSSIANDIST
 
 
!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
!
!   SUBROUTINE SolarFractions
!   Based on actual and potential max solar radiation:
!   Determine the fraction of solar radiation that is
!   diffuse PPFD, direct PPFD, diffuse near IR, direct near IR
!
!   Originally developed by Alex Guenther in 1990s
!   Modified in 2010
!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
 
      SUBROUTINE SOLARFRACTIONS(Solar,Maxsolar,Qdiffv,Qbeamv,Qdiffn,    &
                              & Qbeamn)
 
      IMPLICIT NONE
!*--SOLARFRACTIONS137
 
      REAL , INTENT(IN) :: Solar , Maxsolar
      REAL , INTENT(OUT) :: Qdiffv , Qbeamv , Qdiffn , Qbeamn
      REAL :: fracdiff , ppfdfrac , ppfddiffrac , qv , qn
! internal variables
      REAL :: transmis
!-----------------------------------------------------
      IF ( Maxsolar<=0 ) THEN
         transmis = 0.5
      ELSEIF ( Maxsolar<Solar ) THEN
         transmis = 1.0
      ELSE
         transmis = Solar/Maxsolar
      ENDIF
 
!FracDiff is based on Lizaso 2005
      fracdiff = 0.156 + 0.86/(1+EXP(11.1*(transmis-0.53)))
 
!PPFDfrac is based on Goudrian and Laar 1994
      ppfdfrac = 0.55 - transmis*0.12
 
!PPFDdifFrac is based on data in Jacovides 2007
      ppfddiffrac = fracdiff*(1.06+transmis*0.4)
 
! Calculate  Qdiffv,Qbeamv, Qdiffn, Qbeamn in the subroutine
 
      IF ( ppfddiffrac>1.0 ) ppfddiffrac = 1.0
 
      qv = ppfdfrac*Solar
      Qdiffv = qv*ppfddiffrac
      Qbeamv = qv - Qdiffv
      qn = Solar - qv
      Qdiffn = qn*fracdiff
      Qbeamn = qn - Qdiffn
 
      RETURN
 
      END  SUBROUTINE SOLARFRACTIONS
!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo!
!   Subroutine CanopyRad
!
!   Canopy light environment model
!   Code originally developed by Alex Guenther in 1990s
!   Coded into FORTRAN by Xuemei Wang
!   based on Spitters et al. (1986),
!   Goudrian and van Laar (1994), Leuning (1997)
!   Initial code 8-99,
!   modified 7-2000, 12-2001, 1-2017
!
!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
 
      SUBROUTINE CANOPYRAD(Distgauss,Layers,Lai,Sinbeta,Qbeamv,Qdiffv,  &
                          Qbeamn,Qdiffn,Cantype,Sunfrac,    &
                          Qbabsv,Qdabsv,Qsabsv,Qbabsn,Qdabsn,Qsabsn,   &
                          Sunqv,Shadeqv,Sunqn,Shadeqn,Sunppfd,         &
                          Shadeppfd,Nrcha,Nrtyp)
 
      IMPLICIT NONE
!*--CANOPYRAD194
 
! input
      INTEGER , INTENT(IN) :: Layers , Nrcha , Nrtyp , Cantype
      REAL , INTENT(IN) :: Qbeamv , Qdiffv , Sinbeta , Lai , Qbeamn ,   &
                          Qdiffn
      REAL , DIMENSION(Layers) , INTENT(IN) :: Distgauss
 
! output
      REAL , INTENT(OUT) :: Qbabsv , Qbabsn
 
      REAL , DIMENSION(Layers) , INTENT(OUT) :: Shadeppfd , Sunppfd ,   &
                                       Qdabsv , Qsabsv , Qsabsn ,      &
                                       Shadeqv , Sunqn , Qdabsn ,      &
                                       Sunqv , Shadeqn , Sunfrac
 
!     REAL , DIMENSION(Nrcha,Nrtyp) , INTENT(OUT) :: Canopychar
 
! internal variables
      INTEGER :: i
      REAL :: scatv , scatn , refldv , refldn , reflbv , reflbn , kb ,  &
             kd , kbpv , kbpn , kdpv , kdpn , laidepth , cluster ,     &
             qdabsvl , qsabsvl , qdabsnl , qsabsnl , cantran , laiadj
 
! Stefan-boltzman constant  W m-2 K-4
      REAL , PARAMETER :: CONVERTSHADEPPFD = 4.6
      REAL , PARAMETER :: CONVERTSUNPPFD = 4.0
!---------------------------------------------------------------------
 
! adjust LAI for canopy transparency
      cantran = Canopychar(17,Cantype)
      laiadj = Lai/(1-cantran)
 
      IF ( ((Qbeamv+Qdiffv)>0.001) .AND. (Sinbeta>0.002) .AND.          &
          (laiadj>0.001) ) THEN         ! Daytime
 
! Scattering coefficients (scatV,scatN), diffuse and beam reflection
! coefficients (ref..) for visible or near IR
         scatv = Canopychar(5,Cantype)
         scatn = Canopychar(6,Cantype)
         refldv = Canopychar(7,Cantype)
         refldn = Canopychar(8,Cantype)
         cluster = Canopychar(9,Cantype)
!        print*,'cluster',  Cluster
! Extinction coefficients for black leaves for beam (kb) or diffuse (kd)
         kb = cluster*0.5/Sinbeta
! (0.5 assumes a spherical leaf angle distribution (0.5 = cos (60 deg))
         kd = 0.8*cluster
! (0.8 assumes a spherical leaf angle distribution)
 
         CALL CALCEXTCOEFF(Qbeamv,scatv,kb,kd,reflbv,kbpv,kdpv,Qbabsv)
         CALL CALCEXTCOEFF(Qbeamn,scatn,kb,kd,reflbn,kbpn,kdpn,Qbabsn)
 
         DO i = 1 , Layers
! LAI depth at this layer
            laidepth = laiadj*Distgauss(i)
!fraction of leaves that are sunlit
            Sunfrac(i) = EXP(-kb*laidepth)
 
            CALL CALCRADCOMPONENTS(Qdiffv,Qbeamv,kdpv,kbpv,kb,scatv,    &
                                  refldv,reflbv,laidepth,qdabsvl,       &
                                  qsabsvl)
 
            CALL CALCRADCOMPONENTS(Qdiffn,Qbeamn,kdpn,kbpn,kb,scatn,    &
                                 & refldn,reflbn,laidepth,qdabsnl,      &
                                 & qsabsnl)
 
            Shadeppfd(i) = (qdabsvl+qsabsvl)*CONVERTSHADEPPFD/(1-scatv)
            Sunppfd(i) = Shadeppfd(i)                                   &
                       & + (Qbabsv*CONVERTSUNPPFD/(1-scatv))
            Qdabsv(i) = qdabsvl
            Qsabsv(i) = qsabsvl
            Qdabsn(i) = qdabsnl
            Qsabsn(i) = qsabsnl
            Shadeqv(i) = qdabsvl + qsabsvl
            Sunqv(i) = Shadeqv(i) + Qbabsv
            Shadeqn(i) = qdabsnl + qsabsnl
            Sunqn(i) = Shadeqn(i) + Qbabsn
         ENDDO
 
      ELSE                           ! Night time
 
         Qbabsv = 0
         Qbabsn = 0
 
         Sunfrac(:) = 0.2
         Sunqn(:) = 0
         Shadeqn(:) = 0
         Sunqv(:) = 0
         Shadeqv(:) = 0
         Sunppfd(:) = 0
         Shadeppfd(:) = 0
         Qdabsv(:) = 0
         Qsabsv(:) = 0
         Qdabsn(:) = 0
         Qsabsn(:) = 0
 
      ENDIF
      
      Return
      
      END  SUBROUTINE CANOPYRAD
 
!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
!
!   Subroutine CalcExtCoeff
!   Calculate light extinction coefficients
!   Code originally developed by Alex Guenther in 1990s
!   Coded into FORTRAN by Xuemei Wang
!
!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
 
      SUBROUTINE CALCEXTCOEFF(Qbeam,Scat,Kb,Kd,Reflb,Kbp,Kdp,           &
                            & Qbeamabsorb)
 
      IMPLICIT NONE
!*--CALCEXTCOEFF308
 
      REAL , INTENT(IN) :: Qbeam , Scat , Kb , Kd
      REAL , INTENT(OUT) :: Reflb , Kbp , Kdp , Qbeamabsorb
 
! local variables
      REAL :: p
!-------------------------------------------------------------------
 
      p = (1-Scat)**0.5
      Reflb = 1 - EXP((-2*((1-p)/(1+p))*Kb)/(1+Kb))
 
! Extinction coefficients
      Kbp = Kb*p
      Kdp = Kd*p
! Absorbed beam radiation
      Qbeamabsorb = Kb*Qbeam*(1-Scat)

      RETURN
 
      END  SUBROUTINE CALCEXTCOEFF
 
!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
!
!   Subroutine CalcRadComponents
!   Code originally developed by Alex Guenther in 1990s
!   Coded into FORTRAN by Xuemei Wang
!
!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
 
      SUBROUTINE CALCRADCOMPONENTS(Qdiff,Qbeam,Kdp,Kbp,Kb,Scat,Refld,   &
                                  Reflb,Laidepth,Qdabs,Qsabs)
 
      IMPLICIT NONE
!*--CALCRADCOMPONENTS340
 
      REAL , INTENT(IN) :: Qdiff , Qbeam , Kdp , Kbp , Kb , Scat ,      &
                           Refld , Reflb , Laidepth
      REAL , INTENT(OUT) :: Qdabs , Qsabs
!-------------------------------------------------------------------
 
      Qdabs = Qdiff*Kdp*(1-Refld)*EXP(-Kdp*Laidepth)
      Qsabs = Qbeam*((Kbp*(1-Reflb)*EXP(-Kbp*Laidepth))                 &
              -(Kb*(1-Scat)*EXP(-Kb*Laidepth)))
      
      RETURN
      
      END  SUBROUTINE CALCRADCOMPONENTS
 
!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
!
!   Subroutine CanopyEB
!
!   Canopy energy balance model for estimating leaf temperature
!   Coded into FORTRAN by Xuemei Wang
!   Code developed by Alex Guenther in 1990s
!   based on Goudrian and Laar (1994) and Leuning (1997)
!   Initial code 8-99, modified 7-2000 and 12-2001
!   Modified in 1-2017 by Alex Guenther and Ling Huang
!   to correct IR balance and atmos. emissivity
!   Note: i denotes an array containing a vertical profile
!         through the canopy with 0 (above canopy conditions)
!         plus 1 to number of canopy layers
!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
 
      SUBROUTINE CANOPYEB(Trate,Layers,Distgauss,Canopychar,Cantype,   &
                         Tairk,Humidairpa,Ws,Sunppfd,Shadeppfd,Sunqv,  &
                         Shadeqv,Sunqn,Shadeqn,Sunleaftk,Sunleafsh,    &
                         Sunleaflh,Sunleafir,Shadeleaftk,Shadeleafsh,  &
                         Shadeleaflh,Shadeleafir,Nrcha,Nrtyp,Ws0,      &
                         Tairk0,Humidairpa0)
 
      IMPLICIT NONE
!*--CANOPYEB377
 
! inputs
      INTEGER , INTENT(IN) :: Nrcha , Nrtyp , Layers , Cantype
      REAL , INTENT(IN) :: Trate , Tairk0 , Humidairpa0 , Ws0
      REAL , DIMENSION(Layers) , INTENT(IN) :: Distgauss , Sunqv ,      &
                                       Shadeqv , Sunqn , Shadeqn ,     &
                                       Sunppfd , Shadeppfd
      REAL , DIMENSION(Nrcha,Nrtyp) , INTENT(IN) :: Canopychar
 
! outputs
      REAL , DIMENSION(Layers) , INTENT(OUT) :: Humidairpa , Ws ,      &
                                       Sunleaftk , Sunleafsh ,         &
                                       Sunleaflh , Sunleafir , Tairk , &
                                       Shadeleaftk , Shadeleafsh ,     &
                                       Shadeleaflh , Shadeleafir
 
! local variables
      INTEGER :: i
!     &         Deltah, UnexposedLeafIRin, ExposedLeafIRin, IRin,IRout
      REAL :: cdepth , lwidth , llength , cheight , eps ,               &
               transpiretype , deltah , emissatm , irin , irout
      REAL , DIMENSION(Layers) :: ldepth , wsh
!-----------------------------------------------------------------------
 
      cdepth = Canopychar(1,Cantype)
      lwidth = Canopychar(2,Cantype)
      llength = Canopychar(3,Cantype)
      cheight = Canopychar(4,Cantype)
      eps = Canopychar(10,Cantype)
      transpiretype = Canopychar(11,Cantype)
 
      IF ( Tairk0>288 ) THEN
! Pa m-1  (humidity profile for T < 288)
         deltah = Canopychar(14,Cantype)/cheight
      ELSEIF ( Tairk0>278 ) THEN
         deltah = (Canopychar(14,Cantype)-((288-Tairk0)/10)             &
                  *(Canopychar(14,Cantype)-Canopychar(15,Cantype)))     &
                  /cheight
      ELSE
! Pa m-1  (humidity profile for T <278)
         deltah = Canopychar(15,Cantype)/cheight
      ENDIF
 
      ldepth(:) = cdepth*Distgauss(:)
      Tairk(:) = Tairk0 + (Trate*ldepth(:))               ! check this
      Humidairpa(:) = Humidairpa0 + (deltah*ldepth(:))
 
      wsh(:) = (cheight-ldepth(:)) - (Canopychar(16,Cantype)*cheight)
      Ws(:) = (Ws0*LOG(wsh(:))/LOG(cheight-Canopychar(16,Cantype)*      &
            & cheight))
      WHERE(wsh(:)<0.001)Ws(:) = 0.05
 
      DO i = 1 , Layers
 
! REVISE - Replace UnexposedLeafIR with LeafIR
 
!        IRin     = UnexposedLeafIRin(TairK(i), Eps)
!        ShadeleafIR(i) = 2 * IRin
!        SunleafIR(i) = 0.5*ExposedLeafIRin(HumidairPa0,TairK0)+1.5*IRin
 
! Apparent atmospheric emissivity for clear skies:
! function of water vapor pressure (Pa)
! and ambient Temperature (K) based on Brutsaert(1975)
! referenced in Leuning (1997)
         emissatm = 0.642*(Humidairpa(i)/Tairk(i))**(1./7.)
         irin = LEAFIR(Tairk(i),emissatm)
         Shadeleafir(i) = irin
         Sunleafir(i) = irin
 
      ! Sun
         CALL LEAFEB(Sunppfd(i),Sunqv(i)+Sunqn(i),Sunleafir(i),eps,     &
                     transpiretype,lwidth,llength,Tairk(i),Humidairpa(i)&
                     ,Ws(i),Sunleaftk(i),Sunleafsh(i),Sunleaflh(i),     &
                     irout)
 
         Sunleafir(i) = Sunleafir(i) - irout
 
      ! Shade
         CALL LEAFEB(Shadeppfd(i),Shadeqv(i)+Shadeqn(i),Shadeleafir(i), &
                     eps,transpiretype,lwidth,llength,Tairk(i),         &
                     Humidairpa(i),Ws(i),Shadeleaftk(i),Shadeleafsh(i), &
                     Shadeleaflh(i),irout)
 
         Shadeleafir(i) = Shadeleafir(i) - irout
      ENDDO
      
      RETURN
      
      END  SUBROUTINE CANOPYEB
 
!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
!
!   Subroutine LeafEB
!
!   Leaf energy balance
!   Code originally developed by Alex Guenther in 1990s
!   Coded into FORTRAN by Xuemei Wang
!
!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
 
      SUBROUTINE LEAFEB(Ppfd,Q,Irin,Eps,Transpiretype,Lwidth,Llength,   &
                        Tairk,Humidairpa,Ws,Tleaf,Sh,Lh,Irout)
 
      IMPLICIT NONE
!*--LEAFEB480
 
      REAL , INTENT(IN) :: Eps , Transpiretype , Lwidth , Llength ,     &
                           Ppfd , Q , Irin , Tairk , Humidairpa , Ws
      REAL , INTENT(OUT) :: Irout , Tleaf , Sh , Lh
 
! local variables
 
      INTEGER :: i
!     &        LHairT,Tdelt,Balance,LeafBLC,LeafH,LeafLE,LeafIRout,
      REAL :: humidairkgm3 , ghforced , stomres , iroutairt ,     &
              lathv , ws1 , lhairt , tdelt , balance , gh1 , sh1 , lh1 , e1 ,   &
              irout1 , gh
!----------------------------------------------------
 
      IF ( Ws<=0 ) THEN
         ws1 = 0.001
      ELSE
         ws1 = Ws
      ENDIF
 
      ! Air vapor density kg m-3
      humidairkgm3 = CONVERTHUMIDITYPA2KGM3(Humidairpa,Tairk)
 
      ! Heat convection coefficient (W m-2 K-1) for forced convection.
      ! Nobel page 366
      ghforced = 0.0259/(0.004*((Llength/Ws)**0.5))
 
      ! Stomatal resistence s m-1
      stomres = RESSC(Ppfd)
 
! REVISE - Replace LeafIRout with LeafIR
!      IRoutairT = LeafIROut(tairK, eps)
      iroutairt = LEAFIR(Tairk+tdelt,Eps)
 
      ! Latent heat of vaporization (J Kg-1)
      lathv = LHV(Tairk)
 
      ! Latent heat flux
      lhairt = LEAFLE(Tairk,humidairkgm3,lathv,ghforced,stomres,        &
               Transpiretype)
 
      e1 = (Q+Irin-iroutairt-lhairt)
      IF ( e1.EQ.0. ) e1 = -1.
 
      tdelt = 1.
      balance = 10.
      DO i = 1 , 10
         IF ( ABS(balance)>2 ) THEN
          ! Boundary layer conductance
            gh1 = LEAFBLC(ghforced,tdelt,Llength)
          ! Convective heat flux
            sh1 = LEAFH(tdelt,gh1)
          ! Latent heat of vaporization (J Kg-1)
            lathv = LHV(Tairk+tdelt)
            Lh = LEAFLE(Tairk+tdelt,humidairkgm3,lathv,gh1,stomres,     &
                 Transpiretype)
            lh1 = Lh - lhairt
 
! REVISE - Replace LeafIROut with LeafIR
!          IRout  = LeafIROut(TairK + Tdelt, Eps)
            Irout = LEAFIR(Tairk+tdelt,Eps)
            irout1 = Irout - iroutairt
            tdelt = e1/((sh1+lh1+irout1)/tdelt)
            balance = Q + Irin - Irout - sh1 - Lh
         ENDIF
      ENDDO
 
      IF ( tdelt>10 ) tdelt = 10
      IF ( tdelt<-10 ) tdelt = -10
 
      Tleaf = Tairk + tdelt
      gh = LEAFBLC(ghforced,Tleaf-Tairk,Llength)
      Sh = LEAFH(Tleaf-Tairk,gh)
      Lh = LEAFLE(Tleaf,humidairkgm3,lathv,gh,stomres,Transpiretype)
 
! REVISE - Replace LeafIROut with LeafIR
!      IRout = LeafIROut(Tleaf, Eps)
      Irout = LEAFIR(Tleaf,Eps)
      
      RETURN
      
      END  SUBROUTINE LEAFEB

! FSB the following is the original MEGAN3 code for zenith angle Beta
 
!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
!
!   FUNCTION Calcbeta
!   Calculates the solar zenith angle
!   Code originally developed by Alex Guenther in 1990s
!   Coded into FORTRAN by Xuemei Wang
!
!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
 
      FUNCTION CALCBETA(Day,Lat,Hour)
 
      IMPLICIT NONE
!*--CALCBETA575
 
      INTEGER :: Day
 
      REAL :: Hour , Lat , sindelta , cosdelta , a , b , sinbeta ,      &
              CALCBETA
      REAL , PARAMETER :: PI = 3.14159 , RPI180 = 57.29578
!--------------------------------------------------------------------
      sindelta = -SIN(0.40907)*COS(6.28*(Day+10)/(365))
      cosdelta = (1-sindelta**2.)**0.5
 
      a = SIN(Lat/RPI180)*sindelta
      b = COS(Lat/RPI180)*cosdelta
      sinbeta = a + b*COS(2*PI*(Hour-12)/24)
      CALCBETA = ASIN(sinbeta)*57.29578
 
      END FUNCTION CALCBETA
 
! The following is the original code for eccentricity.
 
!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
!
!   FUNCTION CalcEccentricity
!
!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
 
      FUNCTION CALCECCENTRICITY(Day)
 
      IMPLICIT NONE
!*--CALCECCENTRICITY605
      INTEGER :: Day
      REAL :: CALCECCENTRICITY
!----------------------------------------------------------------
 
      CALCECCENTRICITY = 1 + 0.033*COS(2*3.14*(Day-10)/365)
 
      END  FUNCTION CALCECCENTRICITY
 
!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
!
!   FUNCTION WaterVapPres
!
!   Convert water mixing ratio (kg/kg) to water vapor pressure
!   (Pa or Kpa depending on units of input )
!   Mixing ratio (kg/kg), temp (C), pressure (KPa)
!
!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
 
      FUNCTION WATERVAPPRES(Dens,Pres,Waterairratio)
 
      IMPLICIT NONE
!*--WATERVAPPRES627
      REAL :: Dens , Pres , WATERVAPPRES , Waterairratio
!----------------------------------------------------------------
 
      WATERVAPPRES = (Dens/(Dens+Waterairratio))*Pres
 
      END FUNCTION WATERVAPPRES
 
!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
!
!   FUNCTION Stability
!
!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
 
      FUNCTION STABILITY(Canopychar,Cantype,Solar,Nrcha,Nrtyp)
 
      IMPLICIT NONE
!*--STABILITY644
      INTEGER :: Cantype , Nrcha , Nrtyp
      REAL :: Solar , trateboundary , STABILITY
      REAL , DIMENSION(Nrcha,Nrtyp) :: Canopychar
!----------------------------------------------------------------
 
      trateboundary = 500
 
      IF ( Solar>trateboundary ) THEN
            ! Daytime temperature lapse rate
         STABILITY = Canopychar(12,Cantype)
      ELSEIF ( Solar>0 ) THEN
         STABILITY = Canopychar(12,Cantype)                             &
                     - ((trateboundary-Solar)/trateboundary)            &
                     *(Canopychar(12,Cantype)-Canopychar(13,Cantype))
      ELSE
            ! Nightime temperature lapse rate
         STABILITY = Canopychar(13,Cantype)
      ENDIF
 
      END  FUNCTION STABILITY
 
!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
!
!   FUNCTION ConvertHumidityPa2kgm3
!
!   Saturation vapor density  (kg/m3)
!
!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
 
      FUNCTION CONVERTHUMIDITYPA2KGM3(Pa,Tk)
 
      IMPLICIT NONE
!*--CONVERTHUMIDITYPA2KGM3677
      REAL :: CONVERTHUMIDITYPA2KGM3 , Pa , Tk
!--------------------------------------------------------------------
 
      CONVERTHUMIDITYPA2KGM3 = 0.002165*Pa/Tk
 
      END FUNCTION CONVERTHUMIDITYPA2KGM3
 
!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
!
!   FUNCTION ResSC
!
!   Leaf stomatal cond. resistance s m-1
!
!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
 
      FUNCTION RESSC(Par)
 
      IMPLICIT NONE
!*--RESSC696
      REAL :: Par , scadj , RESSC
!----------------------------------------------------------------
 
      scadj = ((0.0027*1.066*Par)/((1+0.0027*0.0027*Par**2.)**0.5))
 
      IF ( scadj<0.1 ) THEN
         RESSC = 2000.0
      ELSE
         RESSC = 200.0/scadj
      ENDIF
 
      END  FUNCTION RESSC
 
 
!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
!
!   FUNCTION LeafIR
!
!   Calculate IR transfer between leaf and air
!   Added by Alex Guenther and Ling Huang to replace previous
!   MEGAN2.1 IR balance functions
!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
 
      FUNCTION LEAFIR(Tk,Eps)
 
      IMPLICIT NONE
!*--LEAFIR723
      REAL :: Eps , Tk , LEAFIR
! Stefan-boltzman constant  W m-2 K-4
      REAL , PARAMETER :: SB = 0.0000000567
!----------------------------------------------------------------
 
      LEAFIR = Eps*SB*(2*(Tk**4.))
 
      END FUNCTION LEAFIR
 
 
!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
!
!   FUNCTION LHV
!
!   Latent Heat of vaporization(J Kg-1) from Stull p641
!
!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
 
      FUNCTION LHV(Tk)
 
      IMPLICIT NONE
!*--LHV745
      REAL :: Tk , LHV
!----------------------------------------------------------------
 
! REVISE - Replace 273 with 273.15
!      LHV = 2501000 - (2370 * (Tk - 273))
      LHV = 2501000 - (2370*(Tk-273.15))
 
      END  FUNCTION LHV
 
!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
!
!   FUNCTION LeafLE
!
!   Latent energy term in Energy balance
!
!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
 
      FUNCTION LEAFLE(Tleaf,Ambvap,Lathv,Gh,Stomres,Transpiretype)
 
      IMPLICIT NONE
!*--LEAFLE766
      REAL :: Tleaf , Ambvap , Lathv , Gh , Stomres , Transpiretype ,   &
              leafres , vapdeficit , LEAFLE , le
!----------------------------------------------------------------
 
      leafres = (1/(1.075*(Gh/1231))) + Stomres
      vapdeficit = (SVDTK(Tleaf)-Ambvap)
 
! Latent heat of vap (J Kg-1) * vap deficit(Kg m-3) /
!                 leaf resistence (s m-1)
      le = Transpiretype*(1/leafres)*Lathv*vapdeficit
      IF ( le<0 ) THEN
         LEAFLE = 0
      ELSE
         LEAFLE = le
      ENDIF
 
      END FUNCTION LEAFLE
 
!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
!
!   FUNCTION LeafBLC
!
!   Boundary layer conductance
!
!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
 
      FUNCTION LEAFBLC(Ghforced,Tdelta,Llength)
 
      IMPLICIT NONE
!*--LEAFBLC796
      REAL :: Ghforced , Tdelta , Llength , ghfree , LEAFBLC
!----------------------------------------------------------------
 
! This is based on Leuning 1995 p.1198 except using molecular
! conductivity (.00253 W m-1 K-1 Stull p 640) instead of molecular
! diffusivity so that you end up with a heat convection coefficient
! (W m-2 K-1) instead of a conductance for free convection
 
      IF ( Tdelta>=0 ) THEN
         ghfree = 0.5*0.00253*((160000000*Tdelta/(Llength**3.))**0.25)  &
                  /Llength
      ELSE
         ghfree = 0
      ENDIF
      LEAFBLC = Ghforced + ghfree
 
      END  FUNCTION LEAFBLC
 
!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
!
!   FUNCTION LeafH
!
!   Convective energy term in Energy balance (W m-2 heat flux
!      from both sides of leaf)
!
!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
 
      FUNCTION LEAFH(Tdelta,Gh)
 
      IMPLICIT NONE
!*--LEAFH827
      REAL :: Tdelta , Gh , LEAFH
!----------------------------------------------------------------
 
! 2 sides X conductance X Temperature gradient
      LEAFH = 2.0 * Gh * Tdelta
 
      END FUNCTION LEAFH
 
!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
!
!   FUNCTION SvdTk
!
!   Saturation vapor density  (kg/m3)
!
!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
 
      FUNCTION SVDTK(Tk)
 
      IMPLICIT NONE
!*--SVDTK847
      REAL :: Tk , svp , SVDTK
!----------------------------------------------------------------
 
! Saturation vapor pressure (millibars)
      svp = 10**((-2937.4/Tk)-(4.9283*LOG10(Tk))+23.5518)
      SVDTK = 0.2165*svp/Tk
 
      END  FUNCTION SVDTK


    ! FSB These subroutines were in file megvea.f

    !----------------------------------------------------------------
    !
    !   SUBROUTINE GAMMA_CD
    !       Emission response to canopy depath
    !----------------------------------------------------------------
    SUBROUTINE GAMMA_CD(NCOLS,NROWS,Layers,LAI,CDEA)

        IMPLICIT NONE
        ! input
        INTEGER,INTENT(IN)                             :: NCOLS,NROWS,Layers
        REAL,DIMENSION(NCOLS,NROWS),INTENT(IN)         :: LAI
        ! output
        REAL,DIMENSION(NCOLS,NROWS,Layers),INTENT(OUT) :: CDEA

        ! local
        REAL,DIMENSION(Layers) :: Cdepth
        REAL                   :: LAIdepth
        INTEGER                 :: I,J,K

        IF ( Layers .EQ. 5 ) THEN
            Cdepth (1)   = 0.0469101
            Cdepth (2)   = 0.2307534
            Cdepth (3)   = 0.5
            Cdepth (4)   = 0.7692465
            Cdepth (5)   = 0.9530899
        ELSE
            DO K = 1,Layers
                Cdepth(K) =(K - 0.5) /Layers
            END DO
        ENDIF

        DO K = 1, Layers
        DO J = 1, NROWS
        DO I = 1, NCOLS
            LAIdepth = MIN( LAI(I,J) * Cdepth(K), 3.0 )
            CDEA(I,J,K) = CCD1 * LAIdepth + CCD2
        END DO
        END DO
        END DO

        RETURN

    END SUBROUTINE GAMMA_CD
    !----------------------------------------------------------------

    !----------------------------------------------------------------
    !
    !   FUNCTION GAMTLD
    !       EA Temperature response (light dependent emission)
    !----------------------------------------------------------------
    FUNCTION GAMTLD(T1,T24,S)

        IMPLICIT NONE
        REAL,PARAMETER :: Ct2 = 230
        INTEGER        :: S
        REAL           :: T1,T24,T240,Topt, X, Eopt, GAMTLD

        T240 = T24

        IF (T1 < 260.0) THEN
            GAMTLD = 0.0
        ELSE
            ! Temperature at which maximum emission occurs
            Topt = 312.5 + 0.6 * (T240 - 297.0)
            X    = ((1.0 / Topt) - (1.0 / T1)) / 0.00831
            ! Maximum emission (relative to emission at 30 C)
            Eopt = Cleo(S) * EXP(0.05 * (T24 - 297.0)) *          &
                  Exp(0.05*(T240-297.0))

            GAMTLD= Eopt * Ct2 * Exp(Ct1(S) * X) /                &
                  (Ct2 - Ct1(S) * (1.0 - EXP(Ct2 * X)))
        ENDIF

    END FUNCTION GAMTLD
    !----------------------------------------------------------------

    !----------------------------------------------------------------
    !
    !   FUNCTION GAMTLI
    !       EA Temperature response (light independent emission)
    !----------------------------------------------------------------


    FUNCTION GAMTLI(temp,S)

        IMPLICIT NONE

        REAL           :: temp, GAMTLI
        REAL,PARAMETER :: Ts = 303.15
        INTEGER        :: S

        GAMTLI = exp( beta(S)*(temp-Ts) )

    END FUNCTION GAMTLI
    !----------------------------------------------------------------


    !----------------------------------------------------------------
    !
    !   FUNCTION GAMP
    !       EA Light response
    !----------------------------------------------------------------

    FUNCTION GAMP(PPFD1,PPFD24)

        IMPLICIT NONE
        REAL            :: PPFD1, PPFD24, Alpha, C1, GAMP

        IF (PPFD24 < 0.01) THEN
            GAMP= 0.0
        ELSE
            Alpha  = 0.004
            !        C1     = 0.0468 * EXP(0.0005 * (PPFD24 - PSTD))
            !     &          * (PPFD24 ** 0.6)
            C1 = 1.03
            !        GAMP= (Alpha * C1 * PPFD1) / ((1 + Alpha**2. * PPFD1**2.)**0.5)
            
! FSB use SQRT her for clarity and efficiency
            
            GAMP= (Alpha * C1 * PPFD1) / SQRT(1.0 + Alpha**2 * PPFD1**2)
        ENDIF

    END FUNCTION GAMP

    !----------------------------------------------------------------
    !
    !   SUBROUTINE GAMMA_HT
    !   EA response to high temperature
    !
    !----------------------------------------------------------------

    SUBROUTINE GAMMA_HT(NCOLS, NROWS, S, MaxT, GAMHT)

        IMPLICIT NONE
        ! input
        INTEGER,INTENT(IN)                            :: NCOLS, NROWS, S
        REAL,DIMENSION(NCOLS,NROWS),INTENT(IN)        :: MaxT
        ! output
        REAL,DIMENSION(NCOLS,NROWS),INTENT(OUT)       :: GAMHT
        ! local
        INTEGER     :: I,J
        REAL        :: THTK, t1

        DO J = 1,NROWS
        DO I = 1,NCOLS
            THTK = 273.15 + THT(S)
            t1 = THTK + DTHT(S)
            IF (MaxT(I,J) <= THTK) THEN
                GAMHT(I,J) = 1.0
            ELSE IF ( MaxT(I,J) < t1) THEN
                GAMHT(I,J) = 1.0 + (CHT(S) - 1.0)* (MaxT(I,J) -  THTK)/DTHT(S)
            ELSE
                GAMHT(I,J) = CHT(S)
            ENDIF
        END DO
        END DO

        RETURN
    END SUBROUTINE GAMMA_HT
    !----------------------------------------------------------------

    !----------------------------------------------------------------
    !
    !   SUBROUTINE GAMMA_LT
    !   EA response to low temperature
    !
    !----------------------------------------------------------------

    SUBROUTINE GAMMA_LT(NCOLS, NROWS, S, MinT, GAMLT)

        IMPLICIT NONE
        ! input
        INTEGER,INTENT(IN)                       :: NCOLS, NROWS, S
        REAL,DIMENSION(NCOLS,NROWS),INTENT(IN)   :: MinT
        ! output
        REAL,DIMENSION(NCOLS,NROWS),INTENT(OUT)  :: GAMLT
        ! local
        INTEGER      :: I,J
        REAL         :: TLTK, t1

        DO J = 1,NROWS
        DO I = 1,NCOLS
            TLTK = 273.15 + TLT(S)
            t1 = TLTK - DTLT(S)
            IF (MinT(I,J) >= TLTK) THEN
                GAMLT(I,J) = 1.0
            ELSE IF ( MinT(I,J) > t1) THEN
                GAMLT(I,J) = 1.0 + (CLT(S) - 1.0)* (TLTK - MinT(I,J))/DTLT(S)
            ELSE
                GAMLT(I,J) = CLT(S)
            ENDIF
        END DO
        END DO

        RETURN
    END SUBROUTINE GAMMA_LT
    !----------------------------------------------------------------


    !----------------------------------------------------------------
    !
    !   SUBROUTINE GAMMA_HW
    !   EA response to high wind speed
    !
    !----------------------------------------------------------------

    SUBROUTINE GAMMA_HW(NCOLS, NROWS, S, MaxWS, GAMHW)

        IMPLICIT NONE
        ! input
        INTEGER,INTENT(IN)                        :: NCOLS, NROWS, S
        REAL,DIMENSION(NCOLS,NROWS),INTENT(IN)    :: MaxWS
        ! output
        REAL,DIMENSION(NCOLS,NROWS),INTENT(OUT)   :: GAMHW
        ! local
        INTEGER     :: I,J
        REAL        :: t1

        DO J = 1,NROWS
        DO I = 1,NCOLS
            t1 = THW(S) + DTHW(S)
            IF (MaxWS(I,J) <= THW(S)) THEN
                GAMHW(I,J) = 1.0
            ELSE IF ( MaxWS(I,J) < t1) THEN
                GAMHW(I,J) = 1.0 + (CHW(S) - 1.0)* (MaxWs(I,J) - THW(S))/ DTHW(S)
            ELSE
                GAMHW(I,J) = CHW(S)
            ENDIF
        END DO
        END DO

        RETURN
    END SUBROUTINE GAMMA_HW
    !----------------------------------------------------------------



    !----------------------------------------------------------------
    !
    !   SUBROUTINE GAMMA_AQ
    !   EA response to air quality
    !
    !----------------------------------------------------------------

    SUBROUTINE GAMMA_AQ(NCOLS, NROWS, S, AQI, GAMAQ)

        IMPLICIT NONE
        ! input
        INTEGER, INTENT(IN)                       :: NCOLS, NROWS, S
        REAL, DIMENSION(NCOLS,NROWS),INTENT(IN)   :: AQI
        ! output
        REAL, DIMENSION(NCOLS,NROWS),INTENT(OUT)   :: GAMAQ
        ! local
        INTEGER    :: I,J
        REAL       :: t1

        DO J = 1, NROWS
        DO I = 1, NCOLS
            t1 = TAQ(S) + DTAQ(S)
            IF (AQI(I,J) <= TAQ(S)) THEN
                GAMAQ(I,J) = 1.0
            ELSE IF ( AQI(I,J) < t1) THEN
                GAMAQ(I,J) = 1.0 + (CAQ(S) - 1.0)* (AQI(I,J) - TAQ(S))/DTAQ(S)
            ELSE
                GAMAQ(I,J) = CAQ(S)
            ENDIF
        END DO
        END DO

        RETURN
    END SUBROUTINE GAMMA_AQ
    !-----------------------------------------------------------------------

    !-----------------------------------------------------------------------
    !
    ! Subroutine GAMMA_CO2
    !-----------------------------------------------------------------------
    !From Alex Guenther 2017-03-11
    SUBROUTINE GAMMA_CO2( NCOLS, NROWS, GAMCO2 )

        IMPLICIT NONE

        INTEGER, INTENT(IN)                         :: NCOLS, NROWS
        REAL,DIMENSION(NCOLS,NROWS),INTENT(OUT)     :: GAMCO2

        ! local
        REAL    :: Ci, CO2temp, cxxx, cyyy
        INTEGER :: C, R

        CO2temp = CO2
        Ci      = 0.7 * CO2

        IF ( CO2 .EQ. 400.0 ) THEN
            GAMCO2 = 1.0
        ELSE
            DO R = 1, NROWS
            DO C = 1, NCOLS
! FSB set common factors for pipeline                 
            
                cxxx =  Ci**CO2h
                cyyy =  Cstar**CO2h
                
     !      GAMCO2 = ISmax- ((ISmax * Ci**CO2h ) / (Cstar**CO2h + Ci **CO2h))
                GAMCO2(C,R) = ISmax- ((ISmax * cxxx) / (cyyy + cxxx))
                
            END DO
            END DO
        END IF

        RETURN

    END SUBROUTINE GAMMA_CO2

    !-----------------------------------------------------------------------


    !-----------------------------------------------------------------------
    !
    ! Subroutine GAMMA_LAIbidir(gam_LAIbidir,LAI)
    !-----------------------------------------------------------------------
    !From Alex Guenther 2010-01-26
    !If lai < 2 Then
    !gammaLAIbidir= 0.5 * lai
    !ElseIf lai <= 6 Then
    !gammaLAIbidir= 1 - 0.0625 * (lai - 2)
    !Else
    !gammaLAIbidir= 0.75
    !End If
    !
    !     SUBROUTINE GAMMA_LAIbidir returns the gam_LAIbidir values
    !    Xuemei Wang-2010-01-28
    !
    !-----------------------------------------------------------------------
    SUBROUTINE GAMMA_LAIbidir(NCOLS, NROWS,LAI,GAMBD)

        IMPLICIT NONE
        ! input
        INTEGER,INTENT(IN)                          :: NCOLS, NROWS
        REAL,DIMENSION(NCOLS, NROWS),INTENT(IN)     ::  LAI

        ! output
        REAL,DIMENSION(NCOLS, NROWS),INTENT(OUT)    :: GAMBD

        ! local
        INTEGER                                     :: I,J

        DO J = 1, NROWS
        DO I = 1, NCOLS

            IF(LAI(I,J) < 2) THEN
                GAMBD(I,J) =  0.5 * LAI(I,J)
            ELSEIF (LAI(I,J) .LE. 6 ) THEN
                GAMBD(I,J) = 1 - 0.0625 * (LAI(I,J) - 2)
            ELSE
                GAMBD(I,J) = 0.75
            ENDIF

        END DO
        END DO

        RETURN
    END SUBROUTINE GAMMA_LAIbidir
    !----------------------------------------------------------------


    !----------------------------------------------------------------
    !
    !   SUBROUTINE GAMLA
    !
    !     EA leaf age response
    !----------------------------------------------------------------
    !
    !       GAMLA = Fnew*Anew + Fgro*Agro + Fmat*Amat + Fold*Aold
    !       where Fnew = new foliage fraction
    !             Fgro = growing foliage fraction
    !             Fmat = mature foliage fraction
    !             Fold = old foliage fraction
    !             Anew = emission activity for new foliage
    !             Agro = emission activity for growing foliage
    !             Amat = emission activity for mature foliage
    !             Aold = emission activity for old foliage
    !           Age class fractions are determined from LAI changes
    !             LAIc = current LAI
    !             LAIp = past LAI
    !             t  = length of the time step (days)
    !             ti = days between budbreak and emission induction
    !             tm = days between budbreak and peak emission
    !             Tt = average above canopy temperature (K)
    !
    !----------------------------------------------------------------

    SUBROUTINE GAMMA_A( NCOLS, NROWS, S,                      &
          LAIp, LAIc, D_TEMP, GAMLA )

        IMPLICIT NONE
 ! input
        INTEGER,INTENT(IN)                       :: NCOLS,NROWS, S
        REAL,DIMENSION(NCOLS,NROWS),INTENT(IN)   :: D_TEMP, LAIp, LAIc
 ! output
        REAL,DIMENSION(NCOLS,NROWS),INTENT(OUT)  :: GAMLA

        INTEGER :: C, R

        REAL :: Fnew, Fgro, Fmat, Fold
        REAL :: ti,tm
        REAL :: Tt

        !---------------------------------------------------
        ! local parameter arrays
        
        DO R = 1, NROWS
        DO C = 1, NCOLS

            Tt = D_TEMP(C,R)

            !... Calculate foliage fraction

            IF (LAIp(C,R) .LT. LAIc(C,R)) THEN

                !        Calculate ti and tm
                IF (Tt .LE. 303.0) THEN
                    ti = 5.0 + 0.7*(300-Tt)
                ELSE
                    ti = 2.9
                END IF
                tm = 2.3*ti

                !       Calculate Fnew and Fmat, then Fgro and Fold
                !       Fnew
                IF (ti .GE. TSTLEN) THEN
                    Fnew = 1.0 - (LAIp(C,R)/LAIc(C,R))
                ELSE
                    Fnew = (ti/TSTLEN) * ( 1-(LAIp(C,R)/LAIc(C,R)) )
                END IF

                !       Fmat
                IF (tm .GE. TSTLEN) THEN
                    Fmat = LAIp(C,R)/LAIc(C,R)
                ELSE
                    Fmat = (LAIp(C,R)/LAIc(C,R)) +                             &
                          ( (TSTLEN-tm)/TSTLEN ) * ( 1-(LAIp(C,R)/LAIc(C,R)) )
                END IF

                Fgro = 1.0 - Fnew - Fmat
                Fold = 0.0

            ELSE IF (LAIp(C,R) .EQ. LAIc(C,R)) THEN

                Fnew = 0.0
                Fgro = 0.1
                Fmat = 0.8
                Fold = 0.1

            ELSE IF (LAIp(C,R) .GT. LAIc(C,R)) THEN

                Fnew = 0.0
                Fgro = 0.0
                Fold = ( LAIp(C,R)-LAIc(C,R) ) / LAIp(C,R)
                Fmat = 1-Fold

            END IF

            !...  Calculate GAMLA
            GAMLA(C,R) = Fnew*Anew(S) + Fgro*Agro(S) + Fmat*Amat(S) + Fold*Aold(S)
        
        END DO
        END DO

        RETURN
    END SUBROUTINE GAMMA_A

! FSB subroutines used by MEGSEA

!=======================================================================
!=======================================================================
      REAL FUNCTION FERTLZ_ADJ( DATE, LAT )

!***********************************************************************
!  DESCRIPTION:
!    This internal function computes a fertilizer adjustment factor
!    for the given date in yyyyddd format. If it is not growing 
!    season, the adjustment factor is 0; otherwise, it ranges from
!    0.0 to 1.0.
!
!  CALL:
!    GROWSEASON
!
!  HISTORY:
!    07/21/11 : Imported from SMOKE-BEIS v3.14 and modified  (Tan)
!***********************************************************************

      IMPLICIT NONE
            
!.... Function arguments
      INTEGER, INTENT(IN) :: DATE
      REAL,    INTENT(IN) :: LAT

!.... Local variables
      INTEGER  GDAY, GLEN

      CHARACTER(256)  MESG         ! message buffer
!-----------------------------------------------------------------------------

      CALL GROWSEASON( DATE, LAT, GDAY, GLEN )

      IF( GDAY == 0 ) THEN
          FERTLZ_ADJ = 0.
      ELSE IF( GDAY >= 1 .AND. GDAY < 30 ) THEN
          ! first month of growing season
          FERTLZ_ADJ = 1.
      ELSE IF( GDAY >= 30 .AND. GDAY <= 366) THEN
          ! later month of growing season
          FERTLZ_ADJ = 1. + 30. / FLOAT(GLEN) - FLOAT(GDAY) / FLOAT(GLEN)
      ELSE
          WRITE( MESG,94010 ) 'Invalid date specified; date = ',  &
                              DATE, 'growing season day = ',      &
                              GDAY
          CALL M3EXIT( 'FERTLZ_ADJ', 0, 0, MESG, 2 )
      ENDIF

!******************  FORMAT  STATEMENTS   ******************************
94010 FORMAT( A, F10.2, 1X, A, I3, ',', I3 )


      RETURN

      END FUNCTION FERTLZ_ADJ
!=======================================================================
!=======================================================================


!=======================================================================
!=======================================================================
      REAL FUNCTION VEG_ADJ( LAI )

!***********************************************************************
!  DESCRIPTION
!    This internal function computes a vegetation adjustment factor
!    based on LAIv.  See Yienger and Levy 1995
!    VEG_ADJ = (EXP(-0.24*LAIv)+EXP(-0.0525*LAIv))*0.5 
!
!  CALL
!    NONE
!
!  HISTORY:
!***********************************************************************

      IMPLICIT NONE
      
!...  Function arguments
      REAL,    INTENT(IN) :: LAI

!-----------------------------------------------------------------------------

      VEG_ADJ = (EXP(-0.24*LAI)+EXP(-0.0525*LAI))*0.5 

!******************  FORMAT  STATEMENTS   ******************************

      RETURN
      END FUNCTION VEG_ADJ
!=======================================================================
!=======================================================================
            


!=======================================================================
!=======================================================================

!=======================================================================
!=======================================================================
      SUBROUTINE GROWSEASON ( DATE, LAT, GDAY, GLEN )

!***********************************************************************
!  DESCRIPTION
!    This internal function computes the day of the growing season
!    corresponding to the given date in yyyyddd format.
!
!  CALL
!    G2J
!
!  HISTORY:
!    07/21/11 : Imported from SMOKE-BEIS v3.14 and modified  (Tan)
!               Variation of growing season depends on latitude
!               (Guenther)
!    04/22/2019 Converted to 90 format and redid error repotring
!               modified G2j to be completely internal. 
!               DR. FRANCIS S. BINKOWSKI, IE, UNC-CHAPEL HILL
!***********************************************************************

      IMPLICIT NONE

!.......  Function arguments
      INTEGER, INTENT(IN)  :: DATE
      REAL,    INTENT(IN)  :: LAT
      INTEGER, INTENT(OUT) :: GDAY, GLEN 
!.......  External functions

!.......  Local parameters
      INTEGER            :: GSEASON_START
      INTEGER            :: GSEASON_END

!.......  Local variables
      INTEGER  YEAR, MONTH, DAY
      INTEGER  JDAY
      INTEGER  GSJULIAN_START
      INTEGER  GSJULIAN_END
      OPEN(UNIT = 10, FILE = 'growseason_error_outoputs.txt')
!-----------------------------------------------------------------------------
! FSB NOTE: The use of "julian Day" to describe the day of tHE year is
!     technically incorrect. 

 ! The Julian Day Number (JDN) is the integer assigned to a whole solar 
 ! day in the Julian day count starting from noon Universal time, with 
 ! Julian day number 0 assigned to the day starting at noon on Monday, 
 ! January 1, 4713 BCE, proleptic Julian calendar (November 24, 4714 BCE, 
 ! in the proleptic Gregorian calendar), a date at which three 
 ! multi-year cycles started (which are: Indiction, Solar, and Lunar cycles)
 !  and which preceded any dates in recorded history. 
 !
 !  For example for January 1st, 2000 CE  at 00:00:00.0 UT1 the Julian Day 
 !  is 2451544.500000 according to  the U.S. Naval Observatory.



      YEAR = INT( FLOAT( DATE ) / 1000. )
      JDAY = DATE - YEAR * 1000

      IF( JDAY .LT. 1 .OR. JDAY .GT. 366 ) THEN
         WRITE( 10,* ) 'Invalid date specified; date = ',        &
                             DATE, 'jday = ', JDAY
      ENDIF

      IF ( LAT .LE. 23.0 .AND. LAT .GE. -23.0 ) THEN
      ! tropical regions, year round
         GSEASON_START = 0101
         GSEASON_END   = 1231

         GSJULIAN_START = G2J(YEAR, GSEASON_START)
         GSJULIAN_END   = G2J(YEAR, GSEASON_END)
         GDAY = JDAY - GSJULIAN_START + 1
         GLEN = GSJULIAN_END - GSJULIAN_START + 1
      ELSE IF ( LAT .LT. -23.0 ) THEN
      ! southern hemisphere
         IF ( LAT .LT. -60.0 ) THEN
         ! antarctic start = 0 end = 0, no growing
            GDAY = 0
            GLEN = 0
         ELSE
         ! southern hemisphere temperate, NOV, DEC, JAN-MAY
            IF (JDAY .GE. 1101 .AND. JDAY .LE. 1231 ) THEN
              GSEASON_START = 1101
              GSEASON_END   = 1231

              GSJULIAN_START = G2J(YEAR, GSEASON_START)
              GSJULIAN_END   = G2J(YEAR, GSEASON_END)
              GDAY = JDAY - GSJULIAN_START + 1
            ELSE IF (JDAY .GE. 0101 .AND. JDAY .LE. 0531) THEN
              GSEASON_START = 0101
              GSEASON_END   = 0531

              GSJULIAN_START = G2J(YEAR, GSEASON_START)
              GSJULIAN_END   = G2J(YEAR, GSEASON_END)
              GDAY = JDAY - GSJULIAN_START + 1 + 61
            ELSE
              GDAY = 0
            ENDIF
            GLEN = 30 + 31 + G2J(YEAR,0531) - G2J(YEAR,0101) + 1

         ENDIF
      ELSE IF ( LAT .GT. 23.0 ) THEN
      ! northern hemisphere
         IF ( LAT .GT. 65.0 ) THEN
         ! arctic start = 0 end = 0, no growing season
            GDAY = 0
            GLEN = 0
         ELSE
         ! northern hemisphere temperate
         ! start= (lat-23)*4.5            189
         ! end = 365 -((lat-23)*3.3)      226
            GSEASON_START = 0
            GSEASON_END   = 1231

            GSJULIAN_START = 0
            GSJULIAN_END   = G2J(YEAR, GSEASON_END)

            GSJULIAN_START = INT( (LAT-23.0) * 4.5 )
            GSJULIAN_END   = GSJULIAN_END - INT( (LAT-23.0) * 3.3 )
            
            ! UNC added to avoid GDAY excede 366
            IF ( JDAY == 366 .AND. GSJULIAN_START==0 ) GSJULIAN_END = GSJULIAN_END - 1
                        
            IF (JDAY .GE. GSJULIAN_START .AND. JDAY .LE. GSJULIAN_END) THEN
               GDAY = JDAY - GSJULIAN_START + 1
            ELSE
               GDAY = 0
            ENDIF
            GLEN = GSJULIAN_END - GSJULIAN_START + 1
         ENDIF
      ELSE
         WRITE(10,*) 'Invalid LAT = ', LAT
      ENDIF
  
     RETURN

      END SUBROUTINE GROWSEASON

!=======================================================================
!=======================================================================
      
! FSB This is a modified version of function G2J.

      
      INTEGER FUNCTION G2J( YYYY, MMDD )
      IMPLICIT NONE

!.......  Function arguments
      INTEGER, INTENT(IN) :: YYYY
      INTEGER, INTENT(IN) :: MMDD


!.......  Local parameters
      INTEGER :: MM
      INTEGER :: DD
      INTEGER :: K
      INTEGER :: DOY
      LOGICAL :: LEAP


      MM = INT( FLOAT( MMDD ) / 100.0 )
      DD = MMDD - MM * 100

! FSB use internal code to get G2j

 !     G2J = JULIAN( YYYY, MM , DD )

! FSB The following code is taken from NASA subroutine get_DOY 
! Original  Programmer:   David G. Simpson
!            NASA Goddard Space Flight Center
!            Greenbelt, Maryland  2077  Date:         November 20, 2001
!  Modified April 13, 2019 by Dr Francis S. Binkowski to do only Gregorian years


          LEAP = .FALSE.
      
! FSB TEST FOR LEAP YEARS
      
      IF ( MOD(YYYY,4)   .EQ. 0) LEAP = .TRUE.
      IF ( MOD(YYYY,100) .EQ. 0) LEAP = .FALSE.
      IF  (MOD(YYYY,400) .EQ. 0) LEAP = .TRUE.

      IF (LEAP) THEN
         K = 1
      ELSE
         K = 2
      END IF

! FSB CALCULATE DAY OF THE YEAR

      DOY = ( ( 275 * MM) / 9 ) - K * ( ( MM + 9) / 12 ) + DD - 30
      
      G2J = DOY  

      END FUNCTION G2J

!=======================================================================



      SUBROUTINE SOILNOX( JDATE, JTIME, NX, NY,            &
                      TA, LSOIL, ISLTYP, SOILM, SOILT,     &
                      LAIc, LAT,                           &
                      PRECADJ,                             &
                      CFNO, CFNOG )

!***********************************************************************
!  DESCRIPTION:
!  
!     Uses new NO algorithm NO = Normalized*Tadj*Padj*Fadj*Cadj
!     to estimate NO emissions 
!     Information needed to estimate NO emissions
!     Julian Day          (integer)    JDATE
!     Surface Temperature (MCIP field) TA    (K)
!     Soil Moisture       (MCIP field) SOILM (M**3/M**3) (LSOIL)
!          (ratio of volume of water per volume of soil)
!     Soil Temperature    (MCIP field) SOILT (K)         (LSOIL)
!     Soil Type           (MCIP field) ISLTYP            (LSOIL)
!
!     saturation values for soil types (constants)       (LSOIL)
!     FOR PX Version, the Temperature adjustment factor accounts for wet and dry soils
!                and  the precipitation adjustment factor accounts for saturated soils
!     FOR the non-PX version, the basic algorithm remains with a temperature adjustment factor (dry soil)
!                     and no adjustment for saturated soils
!
!
!     The following arrays are updated after each call to SOILNOX
!     PULTYPE   type of NO emission pulse 
!     PULSEDATE julian date for the beginning of an NO pulse 
!     PULSETIME        time for the beginning of an NO pulse
!  
!     The calculation are based on the following paper by J.J. Yienger and H. Levy II
!     J.J. Yienger and H. Levy II, Journal of Geophysical Research, vol 100,11447-11464,1995
!
!     The Temperature Adjustment Factor is based on section 4.2 for wet and dry soils with
!       the following modification (PX version):
!       Instead of classifying soils as either 'wet' or 'dry', the wet and dry adjustment is 
!       calculated at each grid cell.  A linear interpolation between the wet and dry adjustment
!       factor is made using the relative amount of soil moisture in the top layer (1cm)
!       as the interpolating factor.  The relative amount of soil moisture is determined by
!       taking the MCIP soil moisture field and dividing by the saturation value defined for each
!       soil type in the PX version of MCIP
!       the soil temperature is used in PX version
!
!     The Precipation Adjustment factor is based on section 4.1 with the following modifications.
!       The rainrate is computed from the MCIP directly using a 24 hr daily total. 
!       THe types of Pulses as described in YL95 were used to estimate the NO emission
!       rate.  
!
!    Also see the following paper for more information:
!    Proceedings of the Air and Waste Management Association/U.S. Environmental Protection
!    Agency EMission Inventory Conference, Raleigh October 26-28, 1999 Raleigh NC
!    by Tom Pierce and Lucille Bender       
!
!    REFERENCES
!
!    JACQUEMIN B. AND NOILHAN J. (1990), BOUND.-LAYER METEOROL., 52, 93-134.
!    J.J. Yienger and H. Levy II, Journal of Geophysical Research, vol 100,11447-11464,1995
!    T. Pierce and L. Bender, Examining the Temporal Variability of Ammonia and Nitric Oxide Emissions from Agricultural Processes
!       Proceedings of the Air and Waste Management Association/U.S. Environmental Protection
!        Agency EMission Inventory Conference, Raleigh October 26-28, 1999 Raleigh NC
!
!  PRECONDITIONS REQUIRED:
!     Normalized NO emissions, Surface Temperature, Soil Moisture, Soil type,
!     NO emission pulse type, soil moisture from previous time step, julian date
!     of NO emission pulse start, time of NO emission pulse start,
!     soil type, SOIL TYPES, Land use data
!
!  SUBROUTINES AND FUNCTIONS CALLED (directly or indirectly):
!     FERTILIZER_ADJ computes fertlizer adjustment factor
!     VEG_ADJ        computes vegatation adjustment factor
!     GROWSEASON     computes day of growing season
!     
!  REVISION  HISTORY:
!    10/01 : Prototype by GAP
!    10/03 : modified transition to non growing season for jul-oct of the year
!    08/04 : Converted to SMOKE code style by C. Seppanen
!    07/21/11 : Imported form SMOKE-BEIS v3.14 for MEGAN v2.10
!    MAY 13, 2019 made inot f90 format and  improved efficiency -FSB
! 
!***********************************************************************

!        USE SOILNOX_FX

        IMPLICIT NONE
        

!.........  ARGUMENTS and their descriptions
        INTEGER, INTENT (IN)  :: JDATE   !  current simulation date (YYYYDDD)
        INTEGER, INTENT (IN)  :: JTIME   !  current simulation time (HHMMSS)
        INTEGER, INTENT (IN)  :: NX      !  no. columns
        INTEGER, INTENT (IN)  :: NY      !  no. rows

        REAL, INTENT (IN)  ::  TA      ( NX, NY )    !  air temperature (K)
        REAL, INTENT (IN)  ::  SOILM   ( NX, NY )    !  soil moisture (m3/m3)
        REAL, INTENT (IN)  ::  SOILT   ( NX, NY )    !  soil temperature (K)
        REAL, INTENT (IN)  ::  PRECADJ ( NX, NY )    !  precip adjustment
        REAL, INTENT (IN)  ::  LAIc    ( NX, NY )    !  soil temperature (K)
        REAL, INTENT (IN)  ::  LAT     ( NX, NY )    !  Latitude
        REAL, INTENT (IN OUT)  ::  CFNO    ( NX, NY )    !  NO correction factor
        REAL, INTENT (IN OUT)  ::  CFNOG   ( NX, NY )    !  NO correction factor for grass
        
        INTEGER, INTENT (IN)  ::  ISLTYP  ( NX, NY )    !  soil type

        LOGICAL, INTENT (IN) :: LSOIL              ! true: using PX version of MCIP
        
!.........  Local ARRAYS
! Saturation values for 11 soil types from pxpbl.F  (MCIP PX version)
!       PLEIM-XIU LAND-SURFACE AND PBL MODEL (PX-LSM)
! See JACQUEMIN B. AND NOILHAN J. (1990), BOUND.-LAYER METEOROL., 52, 93-134.
        INTEGER, PARAMETER :: MAXSTYPES = 11
        REAL, PARAMETER    :: SATURATION( MAXSTYPES )     =  (/   &       
                              0.395, 0.410, 0.435, 0.485,         &
                              0.451, 0.420, 0.477, 0.476,         &
                              0.426, 0.482, 0.482            /)       

!.........  SCRATCH LOCAL VARIABLES and their descriptions:
        INTEGER       ::   R, C, L      ! counters
        INTEGER       ::   SOILCAT      ! soil category
        
        REAL          ::   CF           ! NO correction factor
        REAL          ::   CFG          ! NO correction factor for grasslands
        REAL          ::  TAIR         ! surface temperature
        REAL          ::   TSOI         ! soil temperature
        REAL          ::   CFNOWET, CFNODRY, RATIO, FAC1, FAC2 
        REAL, PARAMETER ::  const1 = (1. / 3.0 )  * (1.0 / 30.0)
        REAL, PARAMETER ::  const2 =EXP(-0.103 * 30.0)
        CHARACTER(256)  MESG         ! message buffer
        
        CHARACTER(16) :: PROGNAME = 'SOILNOX'   !  program name

!***********************************************************************

 
!.....  Loop through cells
        DO R = 1, NY
        DO C = 1, NX

          TAIR = TA( C, R )         ! unit in degree K

!.......  Check max and min bounds for temperature
          IF (TAIR < 200.0) THEN
             WRITE( 10,* ) 'TAIR=', TAIR,                &
                   'out of range at (C,R)=', C, R
          END IF

          IF (TAIR > 315.0 ) THEN
              WRITE( 10,* ) 'TAIR=', TAIR,               &
                    'out of range at (C,R)=', C, R,             &
                    ' resetting to 315K'
              TAIR = 315.0
          END IF

!.......  CFNOG
          IF( TAIR > 303.00 ) TAIR = 303.00

          IF ( TAIR > 268.8690 ) THEN  
              CFG = EXP( 0.04686 * TAIR - 14.30579 ) ! grass (from BEIS2)
          ELSE
              CFG = 0.0
          END IF

          CFNOG(C,R) = CFG
          
! FSB pre calculate common factors
              FAC1 = (TSOI- 273.16)
              FAC2 = const2


!.......  CFNO
          IF( .NOT. LSOIL ) THEN
          ! no soil

             TSOI = 0.72 * TAIR + 82.28
             IF (TSOI <= 273.16) TSOI = 273.16
             IF (TSOI >= 303.16) TSOI = 303.16
             
             
              
!             CFNODRY = (1./3.) * (1./30.) * (TSOI-273.16)  ! see YL 1995 Equa 9a p. 11452             
              CFNODRY = const1 * FAC1  ! see YL 1995 Equa 9a p. 11452
            
             IF (TSOI <= 283.16) THEN         ! linear cold case
!                 CFNOWET = (TSOI-273.16)*EXP(-0.103*30.0)*0.28 ! see YL 1995 Equ 7b
                 CFNOWET =  FAC1 * FAC2 * 0.28 ! see YL 1995 Equ 7b
                 
             ELSE                             ! exponential case
!                 CFNOWET = EXP(0.103 * (TSOI-273.16)) *  EXP(-0.103 * 30.0)
                 CFNOWET = EXP(0.103 * FAC1) *  FAC2
             END IF
             CF = 0.5 * CFNOWET + 0.5 * CFNODRY

          ELSE
          ! soil

             TSOI = SOILT( C,R )
             IF (TSOI <= 273.16) TSOI = 273.16
             IF (TSOI >= 303.16) TSOI = 303.16

!             CFNODRY = (1./3.)*(1./30.)*(TSOI-273.16)  ! see YL 1995 Equa 9a p. 11452
             CFNODRY = const1 * FAC1  ! see YL 1995 Equa 9a p. 11452
                          
             IF (TSOI <= 283.16) THEN         ! linear cold case
!                CFNOWET = (TSOI-273.16)*EXP(-0.103*30.0)*0.28 ! see YL 1995 Equ 7b
                CFNOWET = FAC1 * FAC2 * 0.28 ! see YL 1995 Equ 7b
                
             ELSE                             ! exponential case
!                CFNOWET = EXP(0.103 * (TSOI-273.16)) * EXP(-0.103 * 30.0)
                 CFNOWET = EXP(0.103 * FAC1 ) * FAC2
               
             END IF

             SOILCAT = INT( ISLTYP( C,R ) )
             IF( SOILCAT > 0 .AND. SOILCAT <= MAXSTYPES ) THEN
             
                 RATIO = SOILM( C,R ) / SATURATION( SOILCAT )
                 CF = RATIO*CFNOWET + (1.0 - RATIO ) * CFNODRY
                 
             ELSE
             
                 CF = 0.0
                 
             END IF

          END IF  ! Endif LSOIL

          CFNO(C,R) = CF *                                      &
                     FERTLZ_ADJ( JDATE, LAT(C,R) ) *           &
                     VEG_ADJ( LAIc(C,R) ) * PRECADJ(C,R)

        END DO  ! loop over columns
        END DO  ! loop over rows

        RETURN

        END SUBROUTINE SOILNOX
        
!//////////////////////////////////////////////////////////////////  
    

       SUBROUTINE CHECKMEM( MSTATUS, ONVAR, CALLER )

       IMPLICIT NONE

!...........   ARGUMENTS and their descriptions:

       INTEGER       MSTATUS !  ALLOCATE function exit status
       CHARACTER*(*) ONVAR   !  Variable name of previous ALLOCATE statement
       CHARACTER*(*) CALLER  !  Name of calling program

!...........   ARGUMENTS and their descriptions:
       INTEGER      TRIMLEN
       EXTERNAL     TRIMLEN

!...........   Local variables

       INTEGER         L1
       INTEGER         L2
       CHARACTER*256   MESG

       CHARACTER*16 :: PROGNAME = 'CHECKMEM' ! program name


!***********************************************************************
!   begin body of function CHECKMEM

!.........  Get lengths of input character strings
        L1 = TRIMLEN( ONVAR )
        L2 = TRIMLEN( CALLER )

!.........  Abort if memory status is non-zero

        IF( MSTATUS .GT. 0 ) THEN
            MESG = 'Failure allocating memory for "' // ONVAR( 1:L1 ) // &
                   '" variable'
            CALL M3EXIT( CALLER( 1:L2 ), 0, 0, MESG, 2 )
        ENDIF

        RETURN

        END

!-----------------------------------------------------------------------
!   Created by Tan 07/28/11
!   Updated by Ling Huang 02/18/17 for MEGAN3: LAI data is saved as
!   LAI1, LAI2, ... LAIS92, instead of one variable with multiple time
!   step.
!-----------------------------------------------------------------------
      SUBROUTINE FINDLAI( IDATE, MXLAI, LAIp_I, LAIc_I)

      IMPLICIT NONE

! input
      INTEGER,INTENT(IN) :: IDATE  ! YYYYJJJ
      INTEGER,INTENT(IN) :: MXLAI
! output
      INTEGER,INTENT(OUT) :: LAIp_I, LAIc_I
! Local
      INTEGER :: JJJ
      REAL    :: XXX

! Calculation

      JJJ = MOD(IDATE,1000)
      XXX = JJJ/8.0
      LAIc_I = CEILING(XXX)

      IF (LAIc_I .EQ. 1) THEN
        LAIp_I = MXLAI
      ELSE
        LAIp_I = LAIc_I - 1
      ENDIF

      RETURN
      END SUBROUTINE FINDLAI
!-----------------------------------------------------------------------`


end module ALL_MEGAN
