        
!//////////////////////////////////////////////////////////////////  
    
SUBROUTINE MEGVEA(NCOLS,NROWS, LAYERS, YEAR, DOY, ZTIME,            &
                    LAIp, LAIc, LDFMAP, OGAMSM, MaxT, MinT,          &
                    MaxWS, AQI, D_TEMP, D_PPFD, SUNT, SHAT,         &
                    SUNF, SUNP, SHAP,ER, NON_DIMGARMA  )
! PURPOSE: Calculate Vegetation emission activity (EA) for each emission
!		class as the product of EA for individual drivers
!		calculated by the following functions
!
! Vegetation Emission Activity (EA) algorithm FUNCTIONS
!
!   GAMTLD: EA Temperature response (light dependent emission) 
!   GAMTLI: EA Temperature response (light independent emission)
!   GAMP: EA Light response
!   GAMTP: combines GAMLD, GAMLI, GAMP to get canopy average
!
!   GAMLA: EA leaf age response 
!   GAMBD: EA bidirectional exchange LAI response
!
!   CDEA: Canopy depth emission response
!
!   GAMHW: EA response to high wind storms
!   GAMHT: EA resposne to high temperature
!   GAMLT: EA response to low temperature
!
!   GAMAQ: EA response to air pollution
!
!   GAMCO2: EA CO2 response (only applied to isoprene)
!   GAMSM: EA response to soil moisture (multiplied with LDF)
!
! INCLUDE FILES
!     'PARMS3.EXT'   ! I/O API parameters
!     'IODECL3.EXT'  ! I/O API function declarations
!     'FDESC3.EXT'   ! I/O API file desc. data structures
!     'MEGVEA.EXT'    ! coefficients
!
!  INPUT Files
!	Single value for each location
!		LDF: Light dependent fraction (for categories other than
!		monoterpene, use constant values from MEGVEA.EXT)
!		AQ:  Air Quality indicator
!	Time series for each location
!		LAI: Leaf Area Index 
!		GAMSM: emission activity response to soil moisture
!		MaxT: Daily maximum temperature (K)
!		MinT: Daily minimum temperature (K)
!		MaxWS: Daily mximum wind speed (m/s)
!               D_TEMP: Daily average temperature (K)
!               D_PPFD: Daily averaged PPFD (umol/m2.s)
!	Hourly time series for each canopy layer in each location
!		sunfrac: fraction of leaves that are sunlit
!		SUNT: leaf Temperature of sunlit leaves (K)
!		SUNP: sun leaf PPFD visible light (micromol/m2/s)
!		SHAT: leaf Temperature of shade leaves (K)
!		SHAP: shade leaf PPFD visible light (micromol/m2/s)
!	
!  OUTPUT Files
!	Hourly time series for each location 
!		Emission activity for each of the 20 emission types
!
!
! HISTORY:
!   Based on code initiated by Alex Guenther in 1990s
!   Coded in FORTRAN as
!   MEGAN: Jack Chen 11/04
!   MEGANv2.04: Tan 11/21/06
!   MEGANv2.1: X. Wang 11/04/2007
!       Modified by Julia Lee-Taylor 03/18/2008
!       Modified by Xuemei Wang 09/30/2008
!       Modified by Tan 07/28/2011
!   MEGAN3.0:
!   Alex Guenther and Ling Huang Feb 2017
!   FSB converted program to subroutine

    USE MEG_EXT
    USE ALL_MEGAN
    IMPLICIT NONE
! INPUT VARIABLES
    INTEGER, INTENT(IN)  ::  NCOLS
    INTEGER, INTENT(IN)  ::  NROWS
    INTEGER, INTENT(IN)  ::  LAYERS
    REAL, INTENT(IN)     ::  YEAR
    REAl, INTENT(IN)     ::  DOY
    REAL, INTENT(IN)     ::  ZTIME
    REAL, INTENT(IN)     ::  LAIp       ( NCOLS, NROWS )
    REAL, INTENT(IN)     ::  LAIc       ( NCOLS, NROWS )
!   REAL, INTENT(IN)     ::  LDF_in      ( NCOLS, NROWS )
    REAL, INTENT(IN)     ::  OGAMSM      ( NCOLS, NROWS )
    REAL, INTENT(IN)     ::  MaxT        ( NCOLS, NROWS )
    REAL, INTENT(IN)     ::  MinT        ( NCOLS, NROWS )
    REAL, INTENT(IN)     ::  MaxWS       ( NCOLS, NROWS )
    REAL, INTENT(IN)     ::  AQI         ( NCOLS, NROWS )
    REAL, INTENT(IN)     ::  D_TEMP      ( NCOLS, NROWS )
    REAL, INTENT(IN)     ::  D_PPFD      ( NCOLS, NROWS )
    REAL, INTENT(IN)     ::  SUNT        ( NCOLS, NROWS, LAYERS )
    REAL, INTENT(IN)     ::  SHAT        ( NCOLS, NROWS, LAYERS )
    REAL, INTENT(IN)     ::  SUNF        ( NCOLS, NROWS, LAYERS )
    REAL, INTENT(IN)     ::  SUNP        ( NCOLS, NROWS, LAYERS )
    REAL, INTENT(IN)     ::  SHAP       ( NCOLS, NROWS, LAYERS )
    REAL, INTENT(IN)     ::  LDFMAP     ( NCLASS, NCOLS, NROWS )    ! light depenedent fraction map

! OUTPUT VARIABLES
    REAL, INTENT(OUT)     :: ER(NCOLS, NROWS )
    REAL, INTENT(OUT)     :: NON_DIMGARMA ( NCLASS, NCOLS, NROWS )

 !LOCAL VARIABLES

    LOGICAL, PARAMETER    :: GAMBD_YN  = .true.
    LOGICAL, PARAMETER    :: GAMAQ_YN  = .true.
    LOGICAL, PARAMETER    :: GAMHT_YN  = .true.
    LOGICAL, PARAMETER    :: GAMLT_YN  = .true.
    LOGICAL, PARAMETER    :: GAMHW_YN  = .true.
    LOGICAL, PARAMETER    :: GAMCO2_YN = .true.

    REAL                  :: VPGWT(LAYERS), Ea1L, Ea2L
    !REAL                  :: VPGWT(:), Ea1L(:), Ea2L(:)

    REAL  :: CDEA   ( NCOLS, NROWS, LAYERS ) ! Emission response to canopy depth


    REAL  :: GAMLA  ( NCOLS, NROWS )     ! EA leaf age response
    REAL  :: GAMSM  ( NCOLS, NROWS )
    REAL  :: GAMAQ  ( NCOLS, NROWS )     ! EA response to air pollution
    REAL  :: GAMBD  ( NCOLS, NROWS )     ! EA bidirectional exchange LAI response
    REAL  :: GAMHT  ( NCOLS, NROWS )     ! EA response to high temperature
    REAL  :: GAMLT  ( NCOLS, NROWS )     ! EA response to low temperature
    REAL  :: GAMHW  ( NCOLS, NROWS )     ! EA response to high wind speed
    REAL  :: GAMCO2 ( NCOLS, NROWS )     ! EA response to CO2
    REAL  :: GAMTP  ( NCOLS, NROWS )     ! combines GAMLD, GAMLI, GAMP to get canopy average
    REAL ::  SUM1, SUM2


    ! loop indices
    INTEGER :: IDATE, ITIME
    INTEGER :: S, T, I, J, K

    !ALLOCATE ( VPGWT  ( Layers ), STAT = IOS )
    !CALL CHECKMEM    ( IOS, 'VPGWT',   PROGNAME )
    !ALLOCATE ( Ea1L  ( Layers ), STAT = IOS )
    !CALL CHECKMEM    ( IOS, 'Ea1L',   PROGNAME )
    !ALLOCATE ( Ea2L  ( Layers ), STAT = IOS )
    !CALL CHECKMEM    ( IOS, 'Ea2L',   PROGNAME )


    ! EA response to canopy temperature/light

    IF ( LAYERS .EQ. 5 ) THEN
        VPGWT(1) = 0.1184635
        VPGWT(2) = 0.2393144
        VPGWT(3) = 0.284444444
        VPGWT(4) = 0.2393144
        VPGWT(5) = 0.1184635
    ELSE
        DO K = 1,LAYERS
          VPGWT(K) = 1.0 / FLOAT( LAYERS )
        END DO
    ENDIF
    write(*,*) "test :",  VPGWT

! First process Factors independent of species emission classes S :

    
! Emission response to canopy depth
    CALL GAMMA_CD( NCOLS, NROWS, LAYERS, LAIc, CDEA )

! EA bidirectional exchange LAI response
!   IF ( GAMBD_YN ) THEN
!       CALL GAMMA_LAIbidir(NCOLS, NROWS, LAIc, GAMBD)
!   ELSE
        GAMBD = 1.0
!   ENDIF

!   IF ( GAMCO2_YN ) THEN
!       CALL GAMMA_CO2(NCOLS, NROWS, GAMCO2)
!   ELSE
        GAMCO2 = 1.0
!   ENDIF

!  Now process all factors dependent on S:

    DO S = 1, NCLASS  ! Loop over all the emission classes

        ! leaf age activity factor:  dependent upon S
        CALL GAMMA_A( NCOLS, NROWS, S, LAIp, LAIc, D_TEMP, GAMLA )

        ! emission activity response to air quality
!       IF ( GAMAQ_YN ) THEN
!           CALL GAMMA_AQ(NCOLS, NROWS, S, AQI, GAMAQ)
!       ELSE
            GAMAQ = 1.0
!       ENDIF

        ! EA response to high temperature
!       IF ( GAMHT_YN ) THEN
            CALL GAMMA_HT(NCOLS, NROWS, S, MaxT, GAMHT)
!       ELSE
            GAMHT = 1.0
!       ENDIF

        ! EA response to low temperature
!       IF ( GAMLT_YN ) THEN
!           CALL GAMMA_LT(NCOLS, NROWS, S, MinT, GAMLT)
!       ELSE
            GAMLT = 1.0
!       ENDIF

        ! EA response to high wind speed
!       IF ( GAMHW_YN ) THEN
!           CALL GAMMA_HW(NCOLS, NROWS, S, MaxWS, GAMHW)
!       ELSE
            GAMHW = 1.0
!       ENDIF

! soil moisture activity factor

!       IF ( GAMSM_YN ) THEN
!          GAMSM = OGAMSM 
!       ELSE
           GAMSM = 1.0
!       ENDIF

  
        DO J = 1, NROWS
        DO I = 1, NCOLS ! Preserve stride 1 for output arrays

            SUM1 = 0.0
            SUM2 = 0.0

            DO K = 1, LAYERS
                Ea1L = CDEA(I,J,K) *                                  &
                      GAMTLD(SunT(I,J,K),D_TEMP(I,J),S) *             &
                      GAMP(SunP(I,J,K),D_PPFD(I,J)) *  SunF(I,J,K) +  &
                      GAMTLD(ShaT(I,J,K),D_TEMP(I,J),S) *             &
                      GAMP(ShaP(I,J,K),D_PPFD(I,J))                   &
                      * (1.0-SunF(I,J,K))
                SUM1 = SUM1 + Ea1L* VPGWT(K)

                Ea2L = GAMTLI(SunT(I,J,K),S)* SunF(I,J,K) +           &
                      GAMTLI(ShaT(I,J,K),S)*(1.0-SunF(I,J,K))
                SUM2 = SUM2 + Ea2L*VPGWT(K)

            END DO   ! END DO canopy layers

            GAMTP(I,J) = SUM1*LDFMAP(S,I,J) + SUM2*( 1.0-LDFMAP(S,I,J) )

 ! ... Calculate emission activity factors
 
            IF ( S .EQ. 1 ) THEN
            
!    GAMCO2 only applied to isoprene
            
                ER(:,:) = LAIc(I,J) * GAMTP(I,J) * GAMCO2(I,J) * GAMLA(I,J) *       &
                          GAMHW(I,J) * GAMAQ(I,J) * GAMHT(I,J) * GAMLT(I,J) *  &
                          GAMSM(I,J) * LDFMAP(S,I,J)

            ELSE IF ( S .EQ. 13 ) THEN
            
 !   GAMBD only applied to ethanol and acetaldehyde
            
                ER(I,J) = LAIc(I,J) * GAMTP(I,J) * GAMBD(I,J) * GAMLA(I,J) *        &
                      GAMHW(I,J) * GAMAQ(I,J) * GAMHT(I,J) * GAMLT(I,J) *      &
                      GAMSM(I,J) * LDFMAP(S,I,J)

 !       For CO (S=19), LDFMAP not applied
 
            ELSE IF ( S .EQ. 19 ) THEN
 !       For CO (S=19), LDFMAP not applied
            
                ER(I,J) = LAIc(I,J) * GAMTP(I,J) * GAMLA(I,J) *                     &
                      GAMHW(I,J) * GAMAQ(I,J) * GAMHT(I,J) * GAMLT(I,J) * GAMSM(I,J)
                    
            ELSE
!  Process remaining species            
            
                ER(I,J) = LAIc(I,J) * GAMTP(I,J) * GAMLA(I,J) *                     &
                      GAMHW(I,J) * GAMAQ(I,J) * GAMHT(I,J) * GAMLT(I,J) *      &
                       GAMSM(I,J) * LDFMAP(S,I,J)

            END IF

            IF ( ER(I,J) .GT. 0.0 ) THEN
                NON_DIMGARMA (S,I,J) = ER(I,J)
            ELSE
                NON_DIMGARMA (S,I,J) = 0.0
            END IF
    
            IF ( I == 118 .AND. J == 143 ) THEN
!           GAMTP(I,J) = SUM1*LDFMAP(S,I,J) + SUM2*( 1.0-LDFMAP(S,I,J) )
            write(*,*) "test, LDFMAP: ", layers,S,LDFMAP(S,I,J), VPGWT(S)
            write(*,*) "test, LDFMAP: ",Ea1L
            write(*,*) "test, LDFMAP: ",Ea2L
            write(*,*) "test vga: ", S,I,J, NON_DIMGARMA(S,I,J)
            write(*,*) "test vga: ",LAIc(I,J),GAMTP(I,J),GAMLA(I,J), GAMHW(I,J),&
                   GAMAQ(I,J),GAMHT(I,J),GAMLT(I,J),GAMSM(I,J),LDFMAP(S,I,J)
            endif

        END DO   ! NCOLS
        END DO ! NROWS

    END DO  ! End loop of species (S)
 
    RETURN
    
   END SUBROUTINE MEGVEA

! ///////////////////////////////////////////////////////////////////////////

