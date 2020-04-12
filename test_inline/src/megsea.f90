
! ///////////////////////////////////////////////////////////////////////////
      SUBROUTINE MEGSEA (IDATE,ITIME,NCOLS,NROWS,SLTYP, CTF,LAIc, LAT,       &
                    TEMP, SOILM, SOILT, PRECADJ,                             &
                    CFNO, CFNOG, GAMSM, GAMNO )


!***********************************************************************
!   This subroutine computes soil NO emission activity factor and isoprene
!   soil moisture activity using MCIP output variables.
!
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
!     FOR PX Version, the Temperature adjustment factor accounts for wet
!     and dry soils
!                and  the precipitation adjustment factor accounts for
!                saturated soils
!     FOR the non-PX version, the basic algorithm remains with a
!     temperature adjustment factor (dry soil)
!                     and no adjustment for saturated soils
!
!     The following arrays are updated after each call to SOILNOX
!     PULTYPE   type of NO emission pulse
!     PULSEDATE julian date for the beginning of an NO pulse
!     PULSETIME        time for the beginning of an NO pulse
!
!     The calculation are based on the following paper by J.J. Yienger
!     and H. Levy II
!     J.J. Yienger and H. Levy II, Journal of Geophysical Research, vol
!     100,11447-11464,1995
!
!     The Temperature Adjustment Factor is based on section 4.2 for wet
!     and dry soils with the following modification (PX version):
!       Instead of classifying soils as either 'wet' or 'dry', the wet
!       and dry adjustment is calculated at each grid cell.  A linear 
!       interpolation between the wet and dry adjustment factor is made 
!       using the relative amount of soil moisture in the top layer (1cm)
!       as the interpolating factor.  The relative amount of soil moisture 
!       is determined by taking the MCIP soil moisture field and dividing by the
!       saturation value defined for each soil type in the PX version of MCIP
!       the soil temperature is used in PX version
!
!     The Precipation Adjustment factor is based on section 4.1 with the
!     following modifications.
!       The rainrate is computed from the MCIP directly using a 24 hr daily total.
!       THe types of Pulses as described in YL95 were used to estimate
!       the NO emission rate.
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
!    T. Pierce and L. Bender, Examining the Temporal Variability of Ammonia and 
!      Nitric Oxide Emissions from Agricultural Proc Proceedings of the Air and Waste 
!      Management Association/U.S. Environmental Protection Agency EMission Inventory 
!      Conference, Raleigh October 26-28, 1999 Raleigh NC
!  PRECONDITIONS REQUIRED:
!     Normalized NO emissions, Surface Temperature, Soil Moisture, Soil type,
!     NO emission pulse type, soil moisture from previous time step, julian date
!     of NO emission pulse start, time of NO emission pulse start, soil type, 
!     SOIL TYPES, Land use data
!
!  SUBROUTINES AND FUNCTIONS CALLED (directly or indirectly):
!     FERTILIZER_ADJ computes fertlizer adjustment factor
!     VEG_ADJ        computes vegatation adjustment factor
!     GROWSEASON     computes day of growing season
!
! HISTORY:
!   07/21/11: Imported from SMOKE-BEIS v3.14 for MEGEAN v2.10 (Tan)
!   03/19/17: Make as an indpendent program (MEGSEA) (Ling Huang)
!   03/31/17: Add calculation for soil moisture activity (Ling Huang)
!*********************************************************************

      USE MEG_EXT
      USE ALL_MEGAN
      IMPLICIT NONE


!    input variables
     
     INTEGER, INTENT(IN) :: IDATE, ITIME  
     INTEGER, INTENT(IN) :: NCOLS, NROWS
     
      INTEGER, INTENT(IN) :: SLTYP  (NCOLS, NROWS)  ! soil type
      REAL, INTENT(IN)    :: CTF( NrTyp, NCOLS, NROWS ) ! Canopy type factor arra
      REAL, INTENT(IN)    :: LAIc( NCOLS, NROWS )    ! Current time step LAI
      REAL, INTENT(IN)    :: LAT (NCOLS, NROWS )    ! Latitude
      REAL, INTENT(IN)    :: TEMP (NCOLS, NROWS)   ! Temperautre (K)
      REAL, INTENT(IN)    :: SOILM  (NCOLS, NROWS)  ! soil moisture
      REAL, INTENT(IN)    :: SOILT  (NCOLS, NROWS)  ! soil temperature
      REAL, INTENT(IN)    :: PRECADJ (NCOLS, NROWS)   



!     output variable
      REAL, INTENT(OUT)   :: CFNO  (NCOLS, NROWS)       ! Emission activity for crop
      REAL, INTENT(OUT)   :: CFNOG  (NCOLS, NROWS)      ! Emission activity for grass
      REAL, INTENT(OUT)   :: GAMSM  (NCOLS, NROWS)      ! Soil moisture activity for isoprene
      REAL, INTENT(OUT)   :: GAMNO  (NCOLS, NROWS)      ! Final NO emission activity

! Local variables and their descriptions:
      CHARACTER*16  :: GDNAM
      CHARACTER*16  :: CNAME        ! Coord name



      INTEGER :: GDAY, GLEN
      INTEGER :: MXLAI,MXCT
      REAL :: T1,WILT,TMO1,TMO2

      LOGICAL :: LSOIL = .TRUE.

      INTEGER :: T,I,J,I_CT


        CALL SOILNOX(IDATE,ITIME,NCOLS,NROWS,            &
                     TEMP,LSOIL,SLTYP, SOILM, SOILT,     &
                     LAIc, LAT, PRECADJ,                 &
                     CFNO, CFNOG )

          DO J = 1,NROWS
          DO I = 1,NCOLS
            CALL GROWSEASON(IDATE,LAT(I,J),GDAY,GLEN)
            IF (GDAY .EQ. 0) THEN
             ! non growing season
             ! CFNOG for everywhere
               GAMNO(I,J) = CFNOG(I,J)

             ELSE IF (GDAY .GT. 0 .AND. GDAY .LE. 366) THEN
             ! growing season
             ! CFNOG for everywhere except crops
             TMO1 = 0.
             TMO2 = 0.
             DO I_CT = 1,5
               TMO1 = TMO1 + CTF(I_CT,I,J)
               TMO2 = TMO2 + CTF(I_CT,I,J) * CFNOG(I,J)
             ENDDO
             ! CFNO for crops
             TMO1 = TMO1 + CTF(6,I,J)
             TMO2 = TMO2 + CTF(6,I,J) * CFNO(I,J)
             IF (TMO1 .EQ. 0.0) THEN
                GAMNO(I,J) = 0.0
             ELSE
                GAMNO(I,J) = TMO2 / TMO1
             ENDIF
             ENDIF
 
           ENDDO  !NCOLS
           ENDDO  !NNOWS


           DO J = 1, NROWS
           DO I = 1, NCOLS
               WILT = WWLT(SLTYP(I,J))
               T1 = WILT + D1
               IF ( SOILM(I,J) < WILT ) THEN
                   GAMSM(I,J) = 0
               ELSE IF ( SOILM(I,J) >= WILT .AND. SOILM(I,J) < T1 ) THEN
                   GAMSM(I,J) = (SOILM(I,J) - WILT)/D1
               ELSE
                   GAMSM(I,J) = 1
               END IF
!          if ( I == 118 .AND. J == 143 ) then
!             write(*,*) "test data: ", IDATE, ITIME
!             write(*,*) "test oMEGSEA: ", WILT,D1,SOILM(I,J),GAMSM(I,J)
!          ENDIF

           END DO ! NCOLS
           END DO ! NROWS
         
END SUBROUTINE MEGSEA


