!------------------------------------------------------------------------!
!  The Community Multiscale Air Quality (CMAQ) system software is in     !
!  continuous development by various groups and is based on information  !
!  from these groups: Federal Government employees, contractors working  !
!  within a United States Government contract, and non-Federal sources   !
!  including research institutions.  These groups give the Government    !
!  permission to use, prepare derivative works of, and distribute copies !
!  of their work in the CMAQ system to the public and to permit others   !
!  to do so.  The United States Environmental Protection Agency          !
!  therefore grants similar permission to use the CMAQ system software,  !
!  but users are requested to provide copies of derivative works or      !
!  products designed to operate in the CMAQ system to the United States  !
!  Government without restrictions as to use by others.  Software        !
!  that is used with the CMAQ system but distributed under the GNU       !
!  General Public License or the GNU Lesser General Public License is    !
!  subject to their copyright restrictions.                              !
!------------------------------------------------------------------------!


! ///////////////////////////////////////////////////////////////////////////

      SUBROUTINE MEGCANOPY( NCOLS, NROWS, YEAR, LAYERS, DOY, ZTIME, &
                 LAT, LONG, LAIc, TEMP, PAR, WIND, PRES, QV, CTF,   &
                 ShadeleafTK, SunleafTK, SunFrac, SunPPFD, ShadePPFD  )

! ----------------------------------------------------------------------------
! FSB START subroutine MEGCANOPY
!     Coded March 19, 2019 by Dr. Francis S. Binkowski
!     modifies 05/07/2019 by FSB to correct loop order
!     Institue for the Environment UNC, Chapel Hill

! FSB Based upon        PROGRAM MEGCAN

!   Based on code initiated by Alex Guenther in 1990s
!   Coded in FORTRAN by Xuemei Wang--Nov. 2007
!   Revised by Alex Guenther and Ling Huang in Feb 2017
!   to correct, modify, and update the code and make it
!   a stand-alone program
!
!*****************************************************************

!
!*****************************************************************
!   Input varibles
!   
!   NROWs, NCOLS         grid location
!   LAYERS               number of vertical layers in canopy
!   DOY                  day of the year
!   Lat                  Latitude
!   Long                 Longitude
!   ZTIME                Hour of the day (UTC/GMT)
!   TEMP                 Temperature [K]
!   PAR                  Photosynthetically active radiation [ W/m**2]
!   Wind                 Wind speed [m s-1]
!   Humidity             Relative humidity [%]
!   Cantype              Defines set of canopy characteristics
!   LAI                  Leaf area index [m2 per m2 ground area]
!   Pres                 Pressure [Pa]
!
!*****************************************************************
! Variables used
!   PPFD           Incoming photosynthetic active radiation [umol/m2/s1]
!   PPFDfrac             Fraction of solar radiation that is PPFD
!   Solar                Solar radiation [W/m2]
!   Maxsolar             Maximum of solar radiation
!   Sinbeta              Sine of solar angle above horizon
!   Betav                Solar angle above horizon
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
! i = 1 is the top canopy layer, 2 is the next layer down , etc.
!   ShadeleafTK          leaf temperature for shade leaves [K] (weighted by canopy type)
!   SunleafTK            leaf temperature for sun leaves [K] (weighted by canopy type)
!   SunFrac              fraction of sun leaves (weighted by canopy type)
!   SunPPFD              PPFD on a sun leaf [umol/m2/s] (weighted by canopy type)
!   ShadePPFD            PPFD on a shade leaf [umol/m2/s] (weighted by canopy type)
!
!*****************************************************************
! FUNCTION S
!   Calcbeta             Calculation of solar zenith angle
!   WaterVapPres         Convert water mixing ratio (kg/kg) to water vapor
!   pressure
!   Stability            Temperature lapse rate in canopy
!   get_BETA             solar position 

       USE MEG_EXT 
       USE ALL_MEGAN
       IMPLICIT NONE
     
! FSB INPUT VARIABLES

       INTEGER, INTENT(IN)   :: NCOLS
       INTEGER, INTENT(IN)   :: NROWS
       INTEGER, INTENT(IN)   :: LAYERS
       REAL, INTENT(IN)      :: YEAR
       REAl, INTENT(IN)      :: DOY
       REAL, INTENT(IN)      :: ZTIME
       REAL, INTENT(IN)      :: LAT ( NCOLS, NROWS )     
       REAL, INTENT(IN)      :: LONG( NCOLS, NROWS )   
       REAL, INTENT(IN)      :: LAIc( NCOLS, NROWS )
       REAL, INTENT(IN)      :: TEMP( NCOLS, NROWS )    
       REAL, INTENT(IN)      :: PAR ( NCOLS, NROWS )     
       REAL, INTENT(IN)      :: WIND( NCOLS, NROWS )
       REAL, INTENT(IN)      :: PRES( NCOLS, NROWS )
       REAL, INTENT(IN)      :: QV  ( NCOLS, NROWS )
       REAL, INTENT(IN)      :: CTF(NRTYP, NCOLS, NROWS ) ! Canopy type factor array

! FSB OUTPUT VARIABLES 

       REAL, INTENT(OUT)     ::  ShadeleafTK ( NCOLS, NROWS, LAYERS ) 
       REAL, INTENT(OUT)     ::  SunleafTK   ( NCOLS, NROWS, LAYERS ) 
       REAL, INTENT(OUT)     ::  SunFrac     ( NCOLS, NROWS, LAYERS ) 
       REAL, INTENT(OUT)     ::  SunPPFD     ( NCOLS, NROWS, LAYERS ) 
       REAL, INTENT(OUT)     ::  ShadePPFD   ( NCOLS, NROWS, LAYERS ) 

! FSB LOCAL VARIABLES
      INTEGER                :: I, I_CT, J, MM, DD
      INTEGER                :: IDAY         ! For using original solar method
      REAL                   :: TotalCT
      REAL                   :: month, Date, JDAY
      REAL                   :: Sinbeta, Betav, HOUR, DAY
      REAL,DIMENSION(LAYERS) ::  VPgausWt, VPgausDis2,             &
       VPgausDis, VPslwWT,  QdAbsV, QsAbsV, QdAbsn,                &
       QsAbsn, SunQv, ShadeQv, SunQn, ShadeQn,                     &
       TairK, HumidairPa, Ws, SunleafSH, sun_ppfd,shade_ppfd,      &           
       SunleafLH,SunleafIR, ShadeleafSH, sun_tk,shade_tk,sun_frac, &
       ShadeleafLH,ShadeleafIR, sun_ppfd_total, shade_ppfd_total,  &
      sun_tk_total, shade_tk_total, sun_frac_total
  
      REAL :: Solar, Maxsolar, Eccentricity,                       &
              Difffrac, PPFDfrac, QbAbsn,                          &                 
               Trate, Qbeamv,Qdiffv, Qbeamn, Qdiffn,               &
               QbAbsV,Ea1tCanopy, Ea1pCanopy,                      &                 
               TairK0, HumidairPa0,Ws0, SH                         
                       
!      REAL ::  CalcEccentricity,WaterVapPres,                      &      
!               Stability, Calcbeta



! FSB Start code        

! FSB calculate the date from Year and day of the year 

      call get_date(YEAR,DOY, MM, DD)
      MONTH = MM ! conver to REAL
      DATE  = DD ! conver to REAL
! FSB Get authentic  Julian Day Number

      JDAY =  getJD (YEAR,MONTH,Date)
    
       DAY  = DOY 
 
       
! FSB revised  loop order 05/07/2019 to have stride 1. 
  
          DO J=1, NROWS
          DO I=1, NCOLS
         
 ! FSB calculate beta as originally done:
 
!                Beta       = Calcbeta(IDAY , Lat(I,J) , Hour )
!                Sinbeta    = SIN(Beta  / 57.29578)
!                Maxsolar   = Sinbeta  * SolarConstant * CalcEccentricity(IDAY )

! FSB calculate the Local Standard Time at the particular gridcell from ZTIME.
            Hour  = ZTIME + LONG(I,J) / 15.0
            IF ( Hour .LT. 0.0 ) THEN
                Hour  = Hour  + 24.0
                DAY  = DAY  - 1.0
            ELSEIF ( Hour  .GT. 24.0 ) THEN
                print*,'Invalid hour: HOUR  is ', Hour ! should write to logfile
            ENDIF

! FSB use better solar position codes here than used in MEGCAN.
            Call  get_BETA( JDAY, Lat(I,J), Hour, Betav, Sinbeta, Eccentricity ) 
            Maxsolar   = Sinbeta  * SolarConstant * Eccentricity
!           write(*,*)  "test Maxsolar: ", Betav, Maxsolar
          
           
            SunleafTK(I,J,:)   = TEMP(I,J)
            ShadeleafTK(I,J,:) = TEMP(I,J)
            SunPPFD(I,J,:)     = PAR (I,J)
            ShadePPFD(I,J,:)   = PAR (I,J)
            SunFrac(I,J,:)     = 1.0
            TotalCT            = 0.0
 
          
            DO I_CT = 1,NRTYP   !canopy types
              TotalCT = TotalCT + CTF(I_CT,I,J) * 0.01
            ENDDO   ! ENDDO I_CT
 

!           Only invoke canopy model when both CT and LAI are valid

            IF (TotalCT .GT. 0.0 .AND. LAIc(I,J) .GT. 0.0) THEN
!           only invoke canopy model when both CT and LAI are valid

              sun_ppfd_total     = 0.0
              shade_ppfd_total   = 0.0
              sun_tk_total       = 0.0
              shade_tk_total     = 0.0
              sun_frac_total     = 0.0
 
              DO I_CT = 1,NRTYP   ! loop overcanopy types
              IF (CTF(I_CT,I,J) .NE. 0.0) THEN
                sun_ppfd           = 0.0
                shade_ppfd         = 0.0
                sun_tk             = 0.0
                shade_tk           = 0.0
                sun_frac           = 0.0
!            Solar angle
                TairK0     = TEMP(I,J)
                Ws0        = WIND(I,J)
                Solar      = PAR(I,J)/2.25  ! Follows original MEGAN3 code

                Call GaussianDist(VPgausDis, Layers)
                Call SolarFractions(Solar, Maxsolar, Qdiffv,Qbeamv,Qdiffn,Qbeamn)
                Call CanopyRad(VPgausDis, Layers, LAIc(I,J), Sinbeta, Qbeamv,    &
                     Qdiffv, Qbeamn, Qdiffn,I_CT , sun_frac,                     &
                     QbAbsV, QdAbsV, QsAbsV, QbAbsn, QdAbsn, QsAbsn, SunQv,      &
                     ShadeQv, SunQn, ShadeQn, sun_ppfd, shade_ppfd,              &
                     NrCha,NrTyp)

                 HumidairPa0  =  WaterVapPres(QV(I,J), PRES(I,J), WaterAirRatio)
                 Trate    =  Stability(Canopychar, I_CT, Solar , NrCha, NrTyp)
                 Call CanopyEB(Trate, Layers, VPgausDis, Canopychar, I_CT,       &
                               TairK, HumidairPa, Ws, sun_ppfd,                  &
                               shade_ppfd, SunQv, ShadeQv, SunQn, ShadeQn,       &
                               sun_tk, SunleafSH, SunleafLH, SunleafIR,          &
                               shade_tk,ShadeleafSH,ShadeleafLH,ShadeleafIR,     &
                               NrCha, NrTyp, Ws0, TairK0, HumidairPa0)
     
                     sun_ppfd_total(:)   = sun_ppfd_total(:) +                  &
                                      0.01*CTF(I_CT,I,J)*sun_ppfd(:)            
     
                     shade_ppfd_total(:) = shade_ppfd_total(:) +                &
                                      0.01*CTF(I_CT,I,J)*shade_ppfd(:)
     
                     sun_tk_total(:)     = sun_tk_total(:) +                    &
                                      0.01*CTF(I_CT,I,J)*sun_tk(:)
                                      
                     shade_tk_total(:)   = shade_tk_total(:) +                  &
                                      0.01*CTF(I_CT,I,J)*shade_tk(:)
                                      
                     sun_frac_total(:)   = sun_frac_total(:) +                  &
                                      0.01*CTF(I_CT,I,J)*sun_frac(:)
              ENDIF
              ENDDO  ! ENDDO I_CT loop


! Calculate outputs
              
              SunleafTK(I,J,:)   = sun_tk_total(:)/TotalCT
              ShadeleafTK(I,J,:) = shade_tk_total(:)/TotalCT
              SunPPFD(I,J,:)     = sun_ppfd_total(:)/TotalCT
              ShadePPFD(I,J,:)   = shade_ppfd_total(:)/TotalCT
              SunFrac(I,J,:)     = sun_frac_total(:)/TotalCT

            ELSE
            ! total CT is zero
            SunleafTK(I,J,:)   = TEMP(I,J)
            ShadeleafTK(I,J,:) = TEMP(I,J)
            SunPPFD(I,J,:)     = PAR (I,J)
            ShadePPFD(I,J,:)   = PAR (I,J)
            SunFrac(I,J,:)     = 1

            ENDIF

           ENDDO   ! ENDDO I
          ENDDO   ! ENDDO J
          
      RETURN 
      END SUBROUTINE MEGCANOPY
      
!//////////////////////////////////////////////////////////////////  
    
