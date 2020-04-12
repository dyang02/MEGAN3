       PROGRAM test_MEG
  
! FSB This is the MEGAN program driver, including: 
! DAYMET,  MEGSEA,  MEGVEA, MEGCAN, MGN2MECH
! FSB March, 2020 UNC Institute of the Environment..

! prerequisite: IOAPI2UAM, MET2MGN

! Based on code initiated by Alex Guenther in 1990s
! Coded in FORTRAN by Xuemei Wang--Nov. 2007
! Revised by Alex Guenther and Ling Huang in Feb 2017
! to correct, modify, and update the code and make it
! a stand-alone program
!
!*****************************************************************
!
!   Select Input/output files before running the program
!*****************************************************************
!   Input varibles
!
!   Day                  Julian day
!   Lat                  Latitude
!   Long                 Longitude
!   Hour                 Hour of the day
!   TEMP                 Temperature [K]
!   PPFD           Incoming photosynthetic active radiation [umol/m2/s1]
!   Wind                 Wind speed [m s-1]
!   Humidity             Relative humidity [%]
!   Cantype             Defines set of canopy characteristics
!   LAI                  Leaf area index [m2 per m2 ground area]
!   Pres                 Pressure [Pa]
!
!*****************************************************************
! Variables used
!
!   PPFDfrac             Fraction of solar radiation that is PPFD
!   Solar                Solar radiation [W/m2]
!   Maxsolar             Maximum of solar radiation
!   Sinbeta              Sin of solar angle above horizon
!   Beta                 Solar angle above horizon
!   TairK0               Above canopy air temperature [K]
!   TairK                Array of canopy air temperature [K]
!   Ws0                  Above canopy wind speed [m/s]
!   Ws                   Array of canopy wind speed [m/s]
!   HumidairPA0          Above canopy ambient humidity [Pa]
!   HumidairPa           Array of canopy ambient humidity in [Pa]
!   Transmis             Transmission of PPFD that is diffuse
!   Difffrac             Fraction of PPFD that is diffuse
!   PPFDfrac             Fraction of solar rad that is PPFD
!   Trate			   temperature vertical profile
!   QbAbsV, QbAbsN       Absorbed direct beam visible/near IR
!   QdAbsV, QdAbsN       Absorbed diffuse visible/near IR
!   QsAbsV, QsAbsN       Absorbed scattered visible//near IR
!   QBeamV, QBeamN       Above canopy direct beam visible/near IR
!   QDiffV, QDiffN       Above canopy diffuse visible/near IR
!
! Arrays with values for each canopy layer (vertical profile)
!   SunleafSH            sensible heat flux for sun leaves [W/m2]
!   SunleafLH            latent heat flux for sun leaves [W/m2]
!   SunleafIR            infrared flux for sun leaves [W/m2]
!   ShadeleafSH          sensible heat for shade leaves [W/m2]
!   ShadeleafLH          latent heat flux for shade leaves [W/m2]
!   ShadeleafIR          infrared flux for shade leaves [W/m2]
!   VPgausDis            gaussian weighting factors for distance
!   SunQv                visible radiation on sun leaves
!   ShadeQv              visible radiation on shade leaves
!   SunQn                near IR radiation on sun leaves
!   ShadeQn              near IR radiation on shade leaves
!   sun_ppfd             Array of incoming (NOT absorbed) PPFD on a sun leaf [umol/m2/s]
!   shade_ppfd           Array of incoming (NOT absorbed) PPFD on a shade leaf [umol/m2/s]
!   sun_tk               Array of leaf temperature for sun leaves [K]
!   shade_tk             Array of leaf temperature for shade leaves [K]
!   sun_frac             Array of the fraction of sun leaves. i = 1 is the top canopy layer, 2 is the next layer, etc.

!*****************************************************************
! OUTPUT
! For each time step and location
! Each variable is an array with a value for each canopy layer
!			       (vertical profile)
! i = 1 is the top canopy layer, 2 is the next layer, etc.
!   ShadeleafTK          leaf temperature for shade leaves [K] (weighted by canopy type)
!   SunleafTK            leaf temperature for sun leaves [K] (weighted by canopy type)
!   SunFrac              fraction of sun leaves (weighted by canopy type)
!   SunPPFD              PPFD on a sun leaf [umol/m2/s] (weighted by canopy type)
!   ShadePPFD            PPFD on a shade leaf [umol/m2/s] (weighted by canopy type)
!
!*****************************************************************
! FUNCTIONS
!   Calcbeta             Calculation of solar zenith angle
!   WaterVapPres         Convert water mixing ratio (kg/kg) to water vapor
!   pressure
!   Stability            Temperature lapse rate in canopy
!   CalcEccentricity     Eccentricity of earth's orbit
!

!...  INCLUDES:
           USE ALL_MEGAN
           USE MEG_EXT
           USE M3UTILIO

      IMPLICIT NONE

!...  EXTERNAL FUNCTIONS and their descriptions:
!     INTEGER, EXTERNAL   ::   ENVINT
!     LOGICAL      DSCGRID
!     EXTERNAL     DSCGRID

!...  Program I/O files: From run script
! Program name
      CHARACTER*16  :: PROGNAME = 'TEST_MEGCAN'
! Netcdf file
      CHARACTER*16  :: CANTYP = 'CANTYP'     ! canopy type file logical name
      CHARACTER*16  :: LAIS46 = 'LAIS46'     ! LAI file logical name
      CHARACTER*16  :: LDFILE   = 'LDFILE'   ! Light dependent Cfraction file
                                             ! for certain MT/SQT
! Met files
      CHARACTER*16  :: MGNMET = 'MGNMET'     ! Met file logical name

! Output file (daily meteorology)
      CHARACTER*16  :: DailyMET = 'DailyMET'   ! Output file logical Cname
!     CHARACTER*16  :: DMET_SPC(5) =  &
!      ( / 'D_TEMP','D_PPFD','MaxT','MinT','MaxWS' /)

! Output file (canopy meteorology)
      CHARACTER*16  :: CANMET = 'CANMET'     ! Output file logical name
!     CHARACTER*16 ::  CAN_SPC(2) =  &
!      (/'SunleafTK', 'ShadeleafTK','SunPPFD','ShadePPFD','SunFrac'/)

! Output file (soil nox)
      CHARACTER*16  :: MGNSEA = 'MGNSEA'     ! Output file logical name
!     CHARACTER*16  :: SEA_SPC(5) =  &
!      (/'GAMNO', 'GAMSM'/)

! output file, MGNVEA
      CHARACTER*16  :: MGNERS = 'MGNERS'       ! Emission activity

!...  Parameters for file units
      INTEGER  LOGDEV                      ! Logfile unit number

! ... External parameters
! From run script
      INTEGER       SDATE          ! Start date YYYYDDD
      INTEGER       STIME          ! Start time HHMMSS
      INTEGER       RLENG          ! Run length HHMMSS

! I/O API file parameters
      INTEGER, SAVE :: JDATE        ! Date YYYYDDD from inpname
      INTEGER, SAVE :: JTIME        ! Time HHMMSS from inpname
      INTEGER, SAVE :: NCOLS        ! Number of columns
      INTEGER, SAVE :: NROWS        ! Number of rows
      INTEGER :: MXREC        ! Total number of timesteps
      INTEGER :: TSTEP        ! Time step

!... Internal parameters
! Internal parameters (status and buffer)
      INTEGER       IOS                    ! i/o status
      CHARACTER*256 MESG                   ! message buffer

! Parameters for output species
      INTEGER, PARAMETER :: NOUT_DAYMET = 5
      INTEGER, PARAMETER :: NOUT_MEGCAN = 5
      INTEGER, PARAMETER :: LAYERS      = 5
                          ! number of output variables: sunleaftk, shadeleaftk,
                          ! sunPPFD, shadePPFD, sunfrac, minT, maxT, maxWS
!
      CHARACTER*16  :: GDNAM
      CHARACTER*16  :: CNAME        ! Coord name

! variables from input/output file

!     variable from LAI file
      REAL, ALLOCATABLE :: LAT( :,: )     ! Latitude of grid cell
      REAL, ALLOCATABLE :: LONG( :,: )    ! Longitude of grid cell
      REAL, ALLOCATABLE :: LAIc( :,: )    ! Current step LAI

!     variable from MGNMET
      REAL, ALLOCATABLE :: TEMP( :,:,: )    ! Temperature (K)
      REAL, ALLOCATABLE :: PPFD( :,:,: )    ! Calculated PAR (umol/m2.s)
      REAL, ALLOCATABLE :: WIND( :,:,: )
      REAL, ALLOCATABLE :: PRES( :,: )
      REAL, ALLOCATABLE :: QV( :,: )
      REAL, ALLOCATABLE :: CTF( :, :, : ) ! Canopy type factor array

      REAL, ALLOCATABLE :: TEMP_C(:,: )    ! Temperature (K)
      REAL, ALLOCATABLE :: PPFD_C(:,: )    ! Calculated PAR (umol/m2.s)
      REAL, ALLOCATABLE :: WIND_C(:,: )

      REAL, ALLOCATABLE :: PRECADJ (:,:)
      INTEGER, ALLOCATABLE :: SLTYP ( :,: ) ! soil type
      REAL, ALLOCATABLE :: SOILM ( :,: ) ! soil moisture
      REAL, ALLOCATABLE :: SOILT ( :,: ) ! soil temperature
      REAL, ALLOCATABLE    :: RSTYP( :,: )

      REAL :: TotalCT

!     variable from MGNVEA
      REAL, ALLOCATABLE :: ER( :,: )      ! Output emission buffer
      REAL, ALLOCATABLE :: NON_DIMGARMA (:,:,:)

!     REAL, ALLOCATABLE :: GAMLA( :,: )      ! EA leaf age response 
      REAL, ALLOCATABLE :: GAMSM( :,: )      ! emission activity response to soil moisture
!     REAL, ALLOCATABLE :: GAMAQ( :,: )      ! EA response to air pollution
      REAL, ALLOCATABLE :: AQI( :,: )        ! air quality index (i.e.W126)
!     REAL, ALLOCATABLE :: GAMBD( :,: )      ! EA bidirectional exchange LAI response
!     REAL, ALLOCATABLE :: GAMHT( :,: )      ! EA response to high temperature
!     REAL, ALLOCATABLE :: GAMLT( :,: )      ! EA response to low temperature
!     REAL, ALLOCATABLE :: GAMHW( :,: )      ! EA response to high wind speed
!     REAL, ALLOCATABLE :: GAMCO2( :,: )     ! EA response to CO2
      REAL, ALLOCATABLE :: LDFMAP( :, :,: )     ! light depenedent fraction map
      REAL, ALLOCATABLE :: LAIp( :,: )    ! Previous time step LAI


      INTEGER :: IDATE, ITIME      ! Looping
      INTEGER :: I_CT, N, T, I, J, S, K
      INTEGER :: LAIp_I, LAIc_I
      INTEGER :: MXCT, MXLAI
      INTEGER :: YR
!     REAL, ALLOCATABLE :: VPGWT(:), Ea1L(:), Ea2L(:)

      LOGICAL :: LSOIL = .TRUE.
        
! daily met output
      REAL, ALLOCATABLE :: D_TEMP  ( :,: )   ! Daily average temperature [K]
      REAL, ALLOCATABLE :: D_PPFD  ( :,: )   ! Daily average solar radiation [umol/m2.s]
      REAL, ALLOCATABLE :: MaxT    ( :,: )   ! Daily maximum temperature [K]
      REAL, ALLOCATABLE :: MinT    ( :,: )   ! Daily minimum temperature [K]
      REAL, ALLOCATABLE :: MaxWS   ( :,: )   ! Daily maximum wind speed [m/s]

! megan canopy output variable
      REAL, ALLOCATABLE :: SunleafTK  ( :,:,: )  ! Sun leaf temperature [K]
      REAL, ALLOCATABLE :: ShadeleafTK( :,:,: )  ! Shade leaf temperature [K]
      REAL, ALLOCATABLE :: SunPPFD    ( :,:,: )  ! PPFD on a sun leaf [umol/m2.s]
      REAL, ALLOCATABLE :: ShadePPFD  ( :,:,: )  ! PPFD on a shade leaf [(umol/m2.s]
      REAL, ALLOCATABLE :: SunFrac    ( :,:,: )   ! fraction of sun leaves
      REAL, ALLOCATABLE :: CDEA (:,:,:)   ! Emission response to canopy depth

! megsea output variable
      REAL, ALLOCATABLE :: CFNO ( :,: )      ! Emission activity for crop
      REAL, ALLOCATABLE :: CFNOG ( :,: )     ! Emission activity for grass
      REAL, ALLOCATABLE :: GAMNO ( :,: )     ! Final NO emission activity


! local variables and their descriptions:

      INTEGER :: Day
      REAL    :: DOY,  ZTIME, YEAR, Hour           
      INTEGER :: AVEBY           ! Divider for daily average

!**************************************************************************
      INTERFACE
!     SUBROUTINE MEGCANOPY( NCOLS, NROWS, YEAR, LAYERS, DOY, ZTIME, &
!                LAT, LONG, LAIc, TEMP, PPFD, WIND, PRES, QV, CTF, &
!                ShadeleafTK, SunleafTK, SunFrac, SunPPFD, ShadePPFD  )

!      INTEGER, INTENT(IN)   :: NCOLS
!      INTEGER, INTENT(IN)   :: NROWS
!      INTEGER, INTENT(IN)   :: LAYERS
!      REAL, INTENT(IN)      :: YEAR
!      REAl, INTENT(IN)      :: DOY
!      REAL, INTENT(IN)      :: ZTIME
!      REAL, INTENT(IN)      :: LAT ( NCOLS, NROWS )
!      REAL, INTENT(IN)      :: LONG( NCOLS, NROWS )
!      REAL, INTENT(IN)      :: LAIc( NCOLS, NROWS )
!      REAL, INTENT(IN)      :: TEMP( NCOLS, NROWS )
!      REAL, INTENT(IN)      :: PPFD( NCOLS, NROWS )
!      REAL, INTENT(IN)      :: WIND( NCOLS, NROWS )
!      REAL, INTENT(IN)      :: PRES( NCOLS, NROWS )
!      REAL, INTENT(IN)      :: QV  ( NCOLS, NROWS )
!      REAL, INTENT(IN)      :: CTF(:, :, : ) ! Canopy type factor array

!      REAL, INTENT(OUT)     ::  ShadeleafTK ( NCOLS, NROWS, LAYERS )
!      REAL, INTENT(OUT)     ::  SunleafTK   ( NCOLS, NROWS, LAYERS )
!      REAL, INTENT(OUT)     ::  SunFrac     ( NCOLS, NROWS, LAYERS )
!      REAL, INTENT(OUT)     ::  SunPPFD     ( NCOLS, NROWS, LAYERS )
!      REAL, INTENT(OUT)     ::  ShadePPFD   ( NCOLS, NROWS, LAYERS )

!      END SUBROUTINE MEGCANOPY

      END INTERFACE

!--========================================================================
!... Begin program
!--========================================================================

!--------------------------------------------------------------------------
!.....1) File set up and assign I/O parameters
!--------------------------------------------------------------------------
!...  Initialize log file unit
      LOGDEV = INIT3()
           !  Now I/O API is set up, and LOGUNIT is the unit number
           !  for the log file (or it 6 for st'd output).

!...  Get input parameters from run script
      CALL ENVSTR( 'GDNAM3D', MESG, 'ASACA36km', GDNAM, IOS )
      IF( .NOT. DSCGRID( GDNAM, CNAME, GDTYP3D,                     &
                   P_ALP3D, P_BET3D, P_GAM3D, XCENT3D, YCENT3D,     &
                   XORIG3D, YORIG3D, XCELL3D, YCELL3D,              &      
                   NCOLS3D, NROWS3D, NTHIK3D ) ) THEN
         MESG = 'Could not get grid description.'
         CALL M3EXIT ( PROGNAME, 0, 0, MESG, 2 )
      ENDIF

!...  Open files
      WRITE(MESG,1030) 'Checking up files',0,0,0
      CALL M3MESG( MESG )

! Canopy type file
      IF ( .NOT. OPEN3( CANTYP, FSREAD3, PROGNAME ) ) THEN
         CALL NAMEVAL (CANTYP, MESG)  ! get input file name and path
         MESG = 'Could not open file '//TRIM(MESG)
         CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
      ENDIF
      ! Check grid
      IF ( .NOT. FILCHK3 ( CANTYP,                                 &
                   GRDDED3, NCOLS3D, NROWS3D, 1, NTHIK3D))  THEN
         MESG = 'CANTYP has differenet grid definition'
         CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
      ENDIF
      IF ( .NOT. DESC3( CANTYP ) ) THEN
         CALL NAMEVAL (CANTYP, MESG)  ! get input file name and path
         MESG = 'Could not get description of '//TRIM(MESG)
         CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
      ENDIF
      MXCT = MXREC3D

! LAI file
      IF ( .NOT. OPEN3( LAIS46, FSREAD3, PROGNAME ) ) THEN
         CALL NAMEVAL (LAIS46, MESG)  ! get input file name and path
         MESG = 'Could not open file '//TRIM(MESG)
         CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
      ENDIF
      ! Check grid
      IF ( .NOT. FILCHK3 ( LAIS46,                                 &
                   GRDDED3, NCOLS3D, NROWS3D, 1, NTHIK3D))  THEN
         MESG = 'LAIS46 has differenet grid definition'
         CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
      ENDIF
      IF ( .NOT. DESC3( LAIS46 ) ) THEN
         CALL NAMEVAL (LAIS46, MESG)  ! get input file name and path
         MESG = 'Could not get description of '//TRIM(MESG)
         CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
      ENDIF
      MXLAI = MXREC3D

! LDF file
      WRITE(MESG,1030) 'Checking up LDF file',0,0,0
      CALL M3MESG( MESG )
      IF ( .NOT. OPEN3( LDFILE, FSREAD3, PROGNAME ) ) THEN
         CALL NAMEVAL (LDFILE, MESG)  ! get input file name and path
         MESG = 'Could not open file '//TRIM(MESG)
         CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
      ENDIF
      ! Check grid
      IF ( .NOT. FILCHK3 ( LDFILE,                              &
                  GRDDED3, NCOLS3D, NROWS3D, 1, NTHIK3D))  THEN
         MESG = 'LDFILE has differenet grid definition'
         CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
      ENDIF
      IF ( .NOT. DESC3( LDFILE ) ) THEN
         CALL NAMEVAL (LDFILE, MESG)  ! get input file name and path
         MESG = 'Could not get description of '//TRIM(MESG)
         CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
      ENDIF

! Met file
      IF ( .NOT. OPEN3( MGNMET, FSREAD3, PROGNAME ) ) THEN
         CALL NAMEVAL (MGNMET, MESG)  ! get input file name and path
         MESG = 'Could not open file '//TRIM(MESG)
         CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
      ENDIF
      ! Check grid
      IF ( .NOT. FILCHK3 ( MGNMET, GRDDED3, NCOLS3D, NROWS3D, 1, NTHIK3D))  THEN
         MESG = 'MGNMET has differenet grid definition'
         CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
      ENDIF
      IF ( .NOT. DESC3( MGNMET ) ) THEN
         CALL NAMEVAL (MGNMET, MESG)  ! get input file name and path
         MESG = 'Could not get description of '//TRIM(MESG)
         CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
      ENDIF
      MXREC = MXREC3D
      NCOLS = NCOLS3D
      NROWS = NROWS3D
      TSTEP = TSTEP3D

!...  Get input parameters from run script
      MESG = 'Model start date (YYYYDDD)'
      SDATE = ENVINT( 'SDATE', MESG, JDATE, IOS )

      MESG = 'Model start time (HHMMSS)'
      STIME = ENVINT( 'STIME', MESG, JTIME, IOS )

      RLENG = MXREC*10000

!...  Check start date, start time, end date, end time in MGNMET
      WRITE(MESG,1030) 'Checking up MGNMET',0,0,0
      CALL M3MESG( MESG )

      IDATE = SDATE; ITIME = STIME
      IF ( .NOT. CHECK3( MGNMET, 'TEMP2', IDATE, ITIME ) ) THEN
         MESG = 'Starting time not on met file'
         CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
      ENDIF
      CALL NEXTIME ( IDATE, ITIME, RLENG-10000 )
      IF ( .NOT. CHECK3( MGNMET, 'TEMP2', IDATE, ITIME ) ) THEN
         MESG = 'Ending time not on met file'
         CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
      ENDIF

!... Set outputs parameters 
      SDATE3D = SDATE                ! From run-script
      STIME3D = STIME                ! From run-script
!     MXREC3D = RLENG / 10000        ! From run-script

!... Set outputs parameters for daily met
      NVARS3D = NOUT_DAYMET
      NLAYS3D = 1 
      NCOLS = NCOLS3D
      NROWS = NROWS3D
      TSTEP = TSTEP3D
      TSTEP3D = 240000
      FDESC3D(:) = ' '
      FDESC3D(1) = 'Daily meteorology with daily average temperature/PPFD, &
                   daily max/min temperature/wind speed'
   
! Define daymet output variables
      VNAME3D(1) = 'D_TEMP'
      UNITS3D(1) = 'K'
      VTYPE3D(1) = M3REAL
      VDESC3D(1) = 'Daily average temperature (K)'

      VNAME3D(2) = 'D_PPFD'
      UNITS3D(2) = 'umol/m2.s'
      VTYPE3D(2) = M3REAL
      VDESC3D(2) = 'Daily average PPFD (umol/m2.s)'

      VNAME3D(3) = 'MaxT'
      UNITS3D(3) = 'K'
      VTYPE3D(3) = M3REAL
      VDESC3D(3) = 'Daily maximum temperature (K)'

      VNAME3D(4) = 'MinT'
      UNITS3D(4) = 'K'
      VTYPE3D(4) = M3REAL
      VDESC3D(4) = 'Daily minimum temperature (K)'

      VNAME3D(5) = 'MaxWS'
      UNITS3D(5) = 'K'
      VTYPE3D(5) = M3REAL
      VDESC3D(5) = 'Daily maximum wind speed (m/s)'

      IF ( .NOT. OPEN3( DailyMET, FSCREA3, PROGNAME ) ) THEN
         CALL NAMEVAL (DailyMET, MESG)  ! get input file name and path
         MESG = 'Could not open file '//TRIM(MESG)
         CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
      ENDIF

! Allocate memory for DailyMET
      ALLOCATE ( TEMP(24,NCOLS,NROWS), STAT = IOS)
      CALL CHECKMEM    ( IOS, 'TEMP',     PROGNAME )
      ALLOCATE ( PPFD(24,NCOLS,NROWS), STAT = IOS)
      CALL CHECKMEM    ( IOS, 'PPFD',     PROGNAME )
      ALLOCATE ( WIND(24,NCOLS,NROWS), STAT = IOS)
      CALL CHECKMEM    ( IOS, 'WIND',     PROGNAME )
      ALLOCATE ( D_TEMP  ( NCOLS, NROWS ), STAT = IOS )
      CALL CHECKMEM    ( IOS, 'D_TEMP',     PROGNAME )
      ALLOCATE ( D_PPFD ( NCOLS, NROWS ), STAT = IOS )
      CALL CHECKMEM    ( IOS, 'D_PPFD',     PROGNAME )
      ALLOCATE ( MaxT  ( NCOLS, NROWS ), STAT = IOS )
      CALL CHECKMEM    ( IOS, 'MaxT',     PROGNAME )
      ALLOCATE ( MinT  ( NCOLS, NROWS ), STAT = IOS )
      CALL CHECKMEM    ( IOS, 'MinT',     PROGNAME )
      ALLOCATE ( MaxWS   ( NCOLS, NROWS ), STAT = IOS )
      CALL CHECKMEM    ( IOS, 'MaxWS',     PROGNAME )



!...  Read in hourly data
      AVEBY = MIN(24,MXREC)
      IDATE = SDATE; ITIME = STIME
!     CALL NEXTIME ( IDATE, ITIME, (I-1)*240000 )
      TEMP = 0.0
      WIND = 0.0
      PPFD = 0.0
      ! Start the loop over the hours
      DO T = 1, AVEBY
          IF ( .NOT. READ3(MGNMET,'TEMP2',1,IDATE,ITIME,TEMP(T,:,:))) THEN
            MESG = 'Error reading TEMP2'
            CALL M3EXIT(PROGNAME,IDATE,ITIME,MESG,2)
          ENDIF

          IF ( .NOT. READ3(MGNMET,'PAR',1,IDATE,ITIME,PPFD(T,:,:))) THEN
             MESG = 'Error reading PAR'
             CALL M3EXIT(PROGNAME,IDATE,ITIME,MESG,2)
          ENDIF

          IF ( .NOT. READ3(MGNMET,'WINDSPD',1,IDATE,ITIME,WIND(T,:,:))) THEN
             MESG = 'Error reading WIND'
             CALL M3EXIT(PROGNAME,IDATE,ITIME,MESG,2)
          ENDIF

          CALL NEXTIME( IDATE, ITIME, TSTEP )
      ENDDO ! End loop over hours

      MaxT = MAXVAL(TEMP,1)
      MinT = MINVAL(TEMP,1)
      MaxWS = MAXVAL(WIND,1)
      D_TEMP = (SUM(TEMP,1))/AVEBY
      ! Convert incoming PAR in W/m2 to umol/m2.s by multiplying 4.5
      D_PPFD = (SUM(PPFD,1))*4.5/AVEBY

!-----------------------------------------------------------------------
!... Write met data to file
        IDATE = SDATE;ITIME = STIME
!       WRITE(MESG,1000) 'Writing met data at ',IDATE,ITIME
!       CALL M3MESG( MESG )
! #1
        IF ( .NOT. WRITE3(DailyMET,'D_TEMP',IDATE,ITIME, D_TEMP)) THEN
          CALL NAMEVAL (DailyMET, MESG)  ! get input file name and path
          MESG = 'Error writing to file: '//TRIM(MESG)
          CALL M3EXIT(PROGNAME,IDATE,ITIME,MESG,2)
        ENDIF
! #2
        IF ( .NOT. WRITE3(DailyMET,'D_PPFD',IDATE,ITIME, D_PPFD)) THEN
          CALL NAMEVAL (DailyMET, MESG)  ! get input file name and path
          MESG = 'Error writing to file: '//TRIM(MESG)
          CALL M3EXIT(PROGNAME,IDATE,ITIME,MESG,2)
        ENDIF
! #3
        IF ( .NOT. WRITE3(DailyMET,'MaxT',IDATE,ITIME, MaxT)) THEN
          CALL NAMEVAL (DailyMET, MESG)  ! get input file name and path
          MESG = 'Error writing to file: '//TRIM(MESG)
          CALL M3EXIT(PROGNAME,IDATE,ITIME,MESG,2)
        ENDIF
! #4
        IF ( .NOT. WRITE3(DailyMET,'MinT',IDATE,ITIME, MinT)) THEN
          CALL NAMEVAL (DailyMET, MESG)  ! get input file name and path
          MESG = 'Error writing to file: '//TRIM(MESG)
          CALL M3EXIT(PROGNAME,IDATE,ITIME,MESG,2)
        ENDIF
! #5
        IF ( .NOT. WRITE3(DailyMET,'MaxWS',IDATE,ITIME, MaxWS)) THEN
          CALL NAMEVAL (DailyMET, MESG)  ! get input file name and path
          MESG = 'Error writing to file: '//TRIM(MESG)
          CALL M3EXIT(PROGNAME,IDATE,ITIME,MESG,2)
        ENDIF

!     DEALLOCATE ( D_TEMP )
!     DEALLOCATE ( D_PPFD )
!     DEALLOCATE ( MaxT )
!     DEALLOCATE ( MinT )
!     DEALLOCATE ( MaxWS )
!     DEALLOCATE ( TEMP )
!     DEALLOCATE ( WIND )
!     DEALLOCATE ( PPFD )
! ... Exit and close file
!     CALL M3EXIT(PROGNAME,0,0,'Finished DAYMET',0)

!--=====================================================================
!...  FORMAT
!--=====================================================================
1000  FORMAT (A20,I8,X,I8)

! METCAN
!--=====================================================================
!...  Set output parameters that are different from met file and open file
      NVARS3D = NOUT_MEGCAN
      NLAYS3D = LAYERS
      MXREC3D = MXREC
      TSTEP3D = 10000

      VNAME3D(1) = 'SunleafTK'
      UNITS3D(1) = 'K'
      VTYPE3D(1) = M3REAL
      VDESC3D(1) = 'Sun leaf temperature (K)'

      VNAME3D(2) = 'ShadeleafTK'
      UNITS3D(2) = 'K'
      VTYPE3D(2) = M3REAL
      VDESC3D(2) = 'Shade leaf temperature (K)'

      VNAME3D(3) = 'SunPPFD'
      UNITS3D(3) = 'umol/m2.s'
      VTYPE3D(3) = M3REAL
      VDESC3D(3) = 'Sun leaf PPFD (umol/m2.s)'

      VNAME3D(4) = 'ShadePPFD'
      UNITS3D(4) = 'umol/m2.s'
      VTYPE3D(4) = M3REAL
      VDESC3D(4) = 'Shade leaf PPFD (umol/m2.s)'

      VNAME3D(5) = 'SunFrac'
      UNITS3D(5) = 'fraction'
      VTYPE3D(5) = M3REAL
      VDESC3D(5) = 'Sun leaf fraction'

      CALL NAMEVAL (CANMET, MESG)  ! get output file name and path
      FDESC3D(:) = ' '
      FDESC3D(1) = 'Output CANMET file: '//TRIM(MESG)

      IF ( .NOT. OPEN3( CANMET, FSCREA3, PROGNAME ) ) THEN
         CALL NAMEVAL (CANMET, MESG)  ! get output file name and path
         MESG = 'Could not open file '//TRIM(MESG)
         CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
      ENDIF

!...  Set output parameters for soil
      NLAYS3D = 1
      NVARS3D = 2

      VNAME3D(1) = 'GAMNO'
      UNITS3D(1) = ' '
      VTYPE3D(1) = M3REAL
      VDESC3D(1) = ' '

      VNAME3D(2) = 'GAMSM'
      UNITS3D(2) = ' '
      VTYPE3D(2) = M3REAL
      VDESC3D(2) = ' '

      IF ( .NOT. OPEN3( MGNSEA, FSCREA3, PROGNAME ) ) THEN
         CALL NAMEVAL (MGNSEA, MESG)  ! get input file name and path
         MESG = 'Could not open file '//TRIM(MESG)
         CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
      ENDIF


!...  Set output parameters that are different from met file and open
!file
      SDATE3D = SDATE                ! From run-script
      STIME3D = STIME                ! From run-script
!     MXREC3D = RLENG / 10000
      NLAYS3D = 1
      NVARS3D = 20

      DO S = 1, NCLASS
         VNAME3D(S) = TRIM( MGN_SPC( S ) )
         VDESC3D(s) = 'Environmental activity factor for '// TRIM( MGN_SPC(S) )
         UNITS3D(s) = 'Non-Dimension '
         VTYPE3D(s) = M3REAL
!         print*,'VNAME=', vname3d(s),VDESC3d(s),UNITS3d(s)
      ENDDO

      IF ( .NOT. OPEN3( MGNERS, FSCREA3, PROGNAME ) ) THEN
         CALL NAMEVAL (MGNERS, MESG)  ! get input file name and path
         MESG = 'Could not open file '//TRIM(MESG)
         CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
      ENDIF

!-----------------------------------------------------------------------
!.....2) Process canopy model
!-----------------------------------------------------------------------
!...  Allocate memory
      ALLOCATE ( SunleafTK ( NCOLS, NROWS, Layers ), STAT = IOS )
      CALL CHECKMEM    ( IOS, 'SunleafTK',    PROGNAME )
      ALLOCATE ( ShadeleafTK ( NCOLS, NROWS, Layers ), STAT = IOS )
      CALL CHECKMEM    ( IOS, 'ShadeleafTK',    PROGNAME )
      ALLOCATE ( SunPPFD ( NCOLS, NROWS, Layers ), STAT = IOS )
      CALL CHECKMEM    ( IOS, 'SunPPFD',    PROGNAME )
      ALLOCATE ( ShadePPFD ( NCOLS, NROWS, Layers ), STAT = IOS )
      CALL CHECKMEM    ( IOS, 'ShadePPFD',    PROGNAME )
      ALLOCATE ( SunFrac ( NCOLS, NROWS, Layers ), STAT = IOS )
      CALL CHECKMEM    ( IOS, 'SunFrac',    PROGNAME )
      ALLOCATE ( LAT   ( NCOLS, NROWS ), STAT = IOS )
      CALL CHECKMEM    ( IOS, 'LAT',    PROGNAME )
      ALLOCATE ( LONG  ( NCOLS, NROWS ), STAT = IOS )
      CALL CHECKMEM    ( IOS, 'LONG',   PROGNAME )
      ALLOCATE ( LAIc  ( NCOLS, NROWS ), STAT = IOS )
      CALL CHECKMEM    ( IOS, 'LAIc',   PROGNAME )
      ALLOCATE ( PPFD_C( NCOLS, NROWS ), STAT = IOS )
      CALL CHECKMEM    ( IOS, 'PPFD_C',   PROGNAME )
      ALLOCATE ( TEMP_C( NCOLS, NROWS ), STAT = IOS )
      CALL CHECKMEM    ( IOS, 'TEMP_C',   PROGNAME )
      ALLOCATE ( WIND_C( NCOLS, NROWS ), STAT = IOS )
      CALL CHECKMEM    ( IOS, 'WIND_C',   PROGNAME )
      ALLOCATE ( PRES  ( NCOLS, NROWS ), STAT = IOS )
      CALL CHECKMEM    ( IOS, 'PRES',   PROGNAME )
      ALLOCATE ( QV    ( NCOLS, NROWS ), STAT = IOS )
      CALL CHECKMEM    ( IOS, 'QV',   PROGNAME )
      ALLOCATE ( CTF( NRTYP, NCOLS, NROWS ), STAT = IOS )
      CALL CHECKMEM ( IOS, 'CTF', PROGNAME )

!.....2) Calculate emission activity
!----------------------------------------------------------------
!...  Allocate memory
      ALLOCATE ( PRECADJ    ( NCOLS, NROWS ), STAT = IOS )
      CALL CHECKMEM    ( IOS, 'PRECADJ',        PROGNAME )
      ALLOCATE ( SOILT    ( NCOLS, NROWS ), STAT = IOS )
      CALL CHECKMEM    ( IOS, 'SOILT',  PROGNAME )
      ALLOCATE ( SOILM ( NCOLS, NROWS ), STAT = IOS )
      CALL CHECKMEM    ( IOS, 'SOILM',  PROGNAME )
      ALLOCATE ( SLTYP ( NCOLS, NROWS ), STAT = IOS )
      CALL CHECKMEM    ( IOS, 'SLTYP',  PROGNAME )
      ALLOCATE ( RSTYP ( NCOLS, NROWS ), STAT = IOS )
      CALL CHECKMEM    ( IOS, 'RSTYP',  PROGNAME )
      ALLOCATE ( CFNO  ( NCOLS, NROWS ), STAT = IOS )
      CALL CHECKMEM    ( IOS, 'CFNO',   PROGNAME )
      ALLOCATE ( CFNOG ( NCOLS, NROWS ), STAT = IOS )
      CALL CHECKMEM    ( IOS, 'CFNOG',  PROGNAME )
      ALLOCATE ( GAMNO ( NCOLS, NROWS ), STAT = IOS )
      CALL CHECKMEM    ( IOS, 'GAMNO',  PROGNAME )
      ALLOCATE ( GAMSM ( NCOLS, NROWS ), STAT = IOS )
      CALL CHECKMEM    ( IOS, 'GAMSM',  PROGNAME )

!----------------------------------------------------------------
!.....2) Calculate emission activity
!----------------------------------------------------------------
!...  Allocate memory
      ALLOCATE ( NON_DIMGARMA (NCLASS, NCOLS, NROWS ), STAT = IOS )
      CALL CHECKMEM ( IOS, 'NON_DIMGARMA', PROGNAME )
      ALLOCATE ( ER ( NCOLS, NROWS ), STAT = IOS )
      CALL CHECKMEM ( IOS, 'ER', PROGNAME )
      ALLOCATE ( LAIp    ( NCOLS, NROWS ), STAT = IOS )
      CALL CHECKMEM    ( IOS, 'LAIp',        PROGNAME )
!     ALLOCATE ( GAMLA   ( NCOLS, NROWS ), STAT = IOS )
!     CALL CHECKMEM    ( IOS, 'GAMLA',       PROGNAME )
      ALLOCATE ( LDFMAP  ( NCLASS, NCOLS, NROWS ), STAT = IOS )
      CALL CHECKMEM    ( IOS, 'LDFMAP',       PROGNAME )
      ALLOCATE ( AQI   ( NCOLS, NROWS ), STAT = IOS )
      CALL CHECKMEM    ( IOS, 'AQI',       PROGNAME )
!     ALLOCATE ( GAMAQ   ( NCOLS, NROWS ), STAT = IOS )
!     CALL CHECKMEM    ( IOS, 'GAMAQ',       PROGNAME )
!     ALLOCATE ( GAMBD   ( NCOLS, NROWS ), STAT = IOS )
!     CALL CHECKMEM    ( IOS, 'GAMBD',       PROGNAME )
!     ALLOCATE ( GAMHT   ( NCOLS, NROWS ), STAT = IOS )
!     CALL CHECKMEM    ( IOS, 'GAMHT',       PROGNAME )
!     ALLOCATE ( GAMLT   ( NCOLS, NROWS ), STAT = IOS )
!     CALL CHECKMEM    ( IOS, 'GAMLT',       PROGNAME )
!     ALLOCATE ( GAMHW   ( NCOLS, NROWS ), STAT = IOS )
!     CALL CHECKMEM    ( IOS, 'GAMHW',       PROGNAME )
!     ALLOCATE ( GAMCO2  ( NCOLS, NROWS ), STAT = IOS )
!     CALL CHECKMEM    ( IOS, 'GAMCO2',      PROGNAME )
!     ALLOCATE ( CDEA  ( NCOLS, NROWS, Layers ),      STAT = IOS )
!     CALL CHECKMEM    ( IOS, 'CDEA',  PROGNAME )
!     ALLOCATE ( VPGWT  ( Layers ), STAT = IOS )
!     CALL CHECKMEM    ( IOS, 'VPGWT',   PROGNAME )

!...  Read LAIS46

      IF ( .NOT. READ3(LAIS46,'LAT',1,0,0,LAT)) THEN
         MESG = 'Error reading LAT'
         CALL M3EXIT(PROGNAME,IDATE,ITIME,MESG,2)
      ENDIF

      IF ( .NOT. READ3(LAIS46,'LONG',1,0,0,LONG)) THEN
         MESG = 'Error reading LONG'
         CALL M3EXIT(PROGNAME,IDATE,ITIME,MESG,2)
      ENDIF

!...  Read CANTYP
      DO N = 1, MXCT
         IF ( .NOT. READ3(CANTYP,'CTS',1,0,(N-1)*10000,CTF(N,:,:))) THEN
            MESG = 'Error reading CTS'
            CALL M3EXIT(PROGNAME,IDATE,ITIME,MESG,2)
         ENDIF
      ENDDO

!...  Start the loop over the time period
      IDATE = SDATE
      ITIME = STIME
      MXREC3D = MXREC
      DO T = 1, MXREC
        WRITE(MESG,1030) 'Processing: ',T,IDATE,ITIME
        CALL M3MESG( MESG )

!...  Initialize hourly variables
        TEMP_C = 0
        PPFD_C = 0
        LAIc = 0

        IF ( .NOT. READ3(MGNMET,'TEMP2',  ALLAYS3,IDATE,ITIME,TEMP_C)) THEN
          MESG = 'Error reading temperature'
          CALL M3EXIT(PROGNAME,IDATE,ITIME,MESG,2)
        ENDIF

        IF ( .NOT. READ3(MGNMET,'PAR',   ALLAYS3,IDATE,ITIME,PPFD_C)) THEN
          MESG = 'Error reading PAR'
          CALL M3EXIT(PROGNAME,IDATE,ITIME,MESG,2)
        ENDIF
        !PPFD = PPFD * 4.766
        PPFD_C = PPFD_C * 4.5

        IF( .NOT. READ3(MGNMET,'WINDSPD',ALLAYS3,IDATE,ITIME,WIND_C)) THEN
          MESG = 'Error reading wind speed'
          CALL M3EXIT(PROGNAME,IDATE,ITIME,MESG,2)
        ENDIF

        IF ( .NOT. READ3(MGNMET,'PRES',  ALLAYS3,IDATE,ITIME,PRES)) THEN
          MESG = 'Error reading pressure'
          CALL M3EXIT(PROGNAME,IDATE,ITIME,MESG,2)
        ENDIF

        IF ( .NOT. READ3(MGNMET,'QV',    ALLAYS3,IDATE,ITIME,QV)) THEN
          MESG = 'Error reading QV'
          CALL M3EXIT(PROGNAME,IDATE,ITIME,MESG,2)
        ENDIF

        IF ( .NOT. READ3(MGNMET,'PREC_ADJ',ALLAYS3,IDATE,ITIME,PRECADJ)) THEN
          MESG = 'Error reading precipitation adjustment'
          CALL M3EXIT(PROGNAME,IDATE,ITIME,MESG,2)
        ENDIF

! ... Calculate emission activity factor for NO
        WRITE(MESG,1030) 'Estimating soil NOx adj: ',T,IDATE,ITIME
        CALL M3MESG( MESG )
        IF ( READ3(MGNMET,'SOIM1', ALLAYS3,IDATE,ITIME,SOILM) .AND. &
             READ3(MGNMET,'SOIT1', ALLAYS3,IDATE,ITIME,SOILT) .AND. &
             READ3(MGNMET,'SLTYP', ALLAYS3,IDATE,ITIME,RSTYP) ) THEN

          MESG = 'Using SOIL parameters in NOx adjustment'
          CALL M3MESG( MESG )
          LSOIL = .TRUE.
          SLTYP = INT(RSTYP)
        ELSE
          MESG = 'SOIL parameters are not available'
          CALL M3MESG( MESG )
          LSOIL = .FALSE.
        ENDIF


        ! Find LAIc from date
        CALL FINDLAI(IDATE,MXLAI,LAIp_I,LAIc_I)
        WRITE(MESG,1020) 'Found LAI current period for YYYYJJJ:', IDATE,LAIc_I
        CALL M3MESG( MESG )

        WRITE(MESG,'(I0.2)') LAIc_I
        IF ( .NOT. READ3(LAIS46,'LAI'//TRIM(MESG),ALLAYS3,0,0,LAIc)) THEN
          MESG = 'Error reading LAI'
          CALL M3EXIT(PROGNAME,IDATE,ITIME,MESG,2)
        ENDIF

        WRITE(MESG,'(I0.2)')  LAIp_I
        IF ( .NOT. READ3(LAIS46,'LAI'//TRIM(MESG),ALLAYS3,0,0,LAIp)) THEN
          MESG = 'Error reading LAI at previous time step'
          CALL M3EXIT(PROGNAME,IDATE,ITIME,MESG,2)
        ENDIF

! FSB The loops over NCOLS, and NROWS are now inside the subroutine
! FSB capture the year and day of the year from IDATE

        GAMSM = 1.0
        DAY  = MOD(IDATE,1000)
        YR   = (IDATE - DAY) / 1000
        DOY  = float(DAY)
        YEAR = float(YR)
                
!       Convert from XXXXXX format to XX.XX (solar hour)
!       HOUR = 0 -> 23.xx
!       Solar hour
        Hour  = float(ITIME)/10000.0 
        ZTIME = Hour
                
        CALL MEGCANOPY( NCOLS, NROWS, YEAR, LAYERS, DOY, ZTIME,  &
                 LAT,LONG,LAIc,TEMP_C,PPFD_C,WIND_C,PRES,QV,CTF, &
                 ShadeleafTK,SunleafTK,SunFrac,SunPPFD,ShadePPFD  ) 
           
        CALL M3MESG( " "  )
!       WRITE(MESG,1030) 'Finished CANOPY run: ',IDATE,ITIME
!       CALL M3MESG( MESG )

        CALL MEGSEA ( IDATE,ITIME,NCOLS,NROWS,SLTYP, CTF,LAIc,LAT, &
                    TEMP_C, SOILM, SOILT, PRECADJ,                 &
                    CFNO, CFNOG, GAMSM, GAMNO )
!       WRITE(MESG,1030) 'Finished MEGSEA run: ',IDATE,ITIME
!       CALL M3MESG( MESG )

!       ! Go over all the emission classes
        DO S = 1, NCLASS 
          IF ( S .EQ. 3 .OR. S .EQ. 4 .OR. S .EQ. 5 .OR. S .EQ. 6 .OR. &
                        S .EQ. 9 .OR. S .EQ. 10 ) THEN
! Read LDF from LDFILE for MT_PINE, MT_ACYC, MT_CAMP, MT_SABI,
! SQT_HR,SQT_LR
            WRITE(MESG,*)  'Reading LDF from LDFILE for '//VNAME3D(S)
            CALL M3MESG(MESG)
            WRITE(MESG,'(A)') VNAME3D(S)
            WRITE(MESG,'(I0.2)') S
            IF ( .NOT. READ3(LDFILE,'LDF'//TRIM(MESG),                   &
                                        ALLAYS3,0,0,LDFMAP(S,:,:) )) THEN
              MESG = 'Error reading LDF for'//VNAME3D(S)
              CALL M3EXIT(PROGNAME,IDATE,ITIME,MESG,2)
            ENDIF
          ELSE
! Read LDF from MEGVEA.EXT
            LDFMAP(S,:,:) = LDF(S)
          ENDIF
       ENDDO 

       CALL MEGVEA(NCOLS,NROWS, LAYERS, YEAR, DOY, ZTIME,           &
                    LAIp, LAIc, LDFMAP,GAMSM, MaxT, MinT, MaxWS,    &
                    AQI, D_TEMP, D_PPFD, SunleafTK, ShadeleafTK,    &
                    SunFrac, SunPPFD,ShadePPFD,ER, NON_DIMGARMA  )
 
   
!-----------------------------------------------------------------------
!.....3) Write out the calculated canmet data
!-----------------------------------------------------------------------
!... Write met data to file
!       WRITE(MESG,1030) 'Writing met data at ',T,IDATE,ITIME
!       CALL M3MESG( MESG )

! #1
        IF ( .NOT. WRITE3(CANMET,'SunleafTK',IDATE,ITIME,         &
                 SunleafTK)) THEN
          CALL NAMEVAL (CANMET, MESG)  ! get input file name and path
          MESG = 'Error writing to file: '//TRIM(MESG)
          CALL M3EXIT(PROGNAME,IDATE,ITIME,MESG,2)
        ENDIF
! #2
        IF ( .NOT. WRITE3(CANMET,'ShadeleafTK',IDATE,ITIME,        &
                 ShadeleafTK)) THEN
          CALL NAMEVAL (CANMET, MESG)  ! get input file name and path
          MESG = 'Error writing to file: '//TRIM(MESG)
          CALL M3EXIT(PROGNAME,IDATE,ITIME,MESG,2)
        ENDIF
! #3
        IF ( .NOT. WRITE3(CANMET,'SunPPFD',IDATE,ITIME,        &
                 SunPPFD)) THEN
          CALL NAMEVAL (CANMET, MESG)  ! get input file name and path
          MESG = 'Error writing to file: '//TRIM(MESG)
          CALL M3EXIT(PROGNAME,IDATE,ITIME,MESG,2)
        ENDIF
!#4
        IF ( .NOT. WRITE3(CANMET,'ShadePPFD',IDATE,ITIME,         &
                 ShadePPFD)) THEN
          CALL NAMEVAL (CANMET, MESG)  ! get input file name and path
          MESG = 'Error writing to file: '//TRIM(MESG)
          CALL M3EXIT(PROGNAME,IDATE,ITIME,MESG,2)
        ENDIF
! #5
        IF ( .NOT. WRITE3(CANMET,'SunFrac',IDATE,ITIME,          &
                    SunFrac)) THEN
          CALL NAMEVAL (CANMET, MESG)  ! get input file name and path
          MESG = 'Error writing to file: '//TRIM(MESG)
          CALL M3EXIT(PROGNAME,IDATE,ITIME,MESG,2)
        ENDIF


!...  Allocate memory
!----------------------------------------------------------------
!.....3) Write out the calculated EA
!----------------------------------------------------------------
!... Write emission to file
!       WRITE(MESG,1030) 'Writing emission at ',T,IDATE,ITIME
!       CALL M3MESG( MESG )
! #1
        IF ( .NOT. WRITE3(MGNSEA,'GAMNO',IDATE,ITIME,GAMNO(:,:))) THEN
          CALL NAMEVAL (MGNSEA, MESG)  ! get input file name and path
          MESG = 'Error writing to file: '//TRIM(MESG)
          CALL M3EXIT(PROGNAME,IDATE,ITIME,MESG,2)
        ENDIF
! #2
        IF ( .NOT. WRITE3(MGNSEA,'GAMSM',IDATE,ITIME, GAMSM(:,:))) THEN
          CALL NAMEVAL (MGNSEA, MESG)  ! get input file name and path
          MESG = 'Error writing to file: '//TRIM(MESG)
          CALL M3EXIT(PROGNAME,IDATE,ITIME,MESG,2)
        ENDIF
!----------------------------------------------------------------
!.....3) Write out the calculated EA, MEGVEA
!----------------------------------------------------------------
!... Write emission to file
!       WRITE(MESG,1030) 'Writing emission at ',T,IDATE,ITIME
!       CALL M3MESG( MESG )
        DO S = 1,NCLASS
          IF (.NOT. WRITE3(MGNERS, VNAME3D(S),IDATE,ITIME,         &
                   NON_DIMGARMA (S,:,:))) THEN
            CALL NAMEVAL (MGNERS, MESG)  ! get input file name and path
            MESG = 'Error writing to file: '//TRIM(MESG)
            CALL M3EXIT(PROGNAME,IDATE,ITIME,MESG,2)
          ENDIF

        ENDDO ! End loop for emission species (S)


       CALL NEXTIME( IDATE, ITIME, TSTEP )
      ENDDO ! End loop for time step (T)
! ... Exit and close file
      CALL M3EXIT(PROGNAME,0,0,' ',0)

!
! ... Exit and close file
!     CALL M3EXIT(PROGNAME,0,0,' ',0)

!     DEALLOCATE( SunleafTK     )
!     DEALLOCATE( ShadeleafTK   )
!     DEALLOCATE( SunPPFD       )
!     DEALLOCATE( ShadePPFD     )
!     DEALLOCATE( SunFrac       )
!     DEALLOCATE ( LAT     )   ! input latitude of grid cell
!     DEALLOCATE ( LONG    )   ! input longitude of grid cell
!     DEALLOCATE ( LAIc    )   ! current monthly LAI
!     DEALLOCATE ( TEMP_C  )   ! input hourly temperature (K)
!     DEALLOCATE ( PPFD_C  )   ! calculated PAR (umol/m2.s)
!     DEALLOCATE ( WIND_C  )
!     DEALLOCATE ( PRES    )
!     DEALLOCATE ( QV      )
!     DEALLOCATE ( CTF    )

!--=====================================================================
!...  FORMAT
!--=====================================================================
1001  FORMAT( A )
1010  FORMAT( 43( A, :, I8, :, 1X ) )
1020  FORMAT (A40,I8,X,I8,X,I8)
1030  FORMAT (A20,I8,X,I8,X,I8)

!--=====================================================================
!...  End program MEGCAN
!--=====================================================================      

      END PROGRAM test_MEG
