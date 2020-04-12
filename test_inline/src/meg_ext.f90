module MEG_EXT

      IMPLICIT NONE

!=======================================================================
!  MEGCAN.EXT
!  This include file contains information required
!  for running MEGCAN canopy environment model
!
!  Who             When       What
!  ---------------------------------------------------------------------
!  Xuemei Wang     06/16/2009 - Created inputs for CANOPY.EXT
!                               Some of these are used in this file
!  Alex Guenther   01/28/2017 - Created this file
!  UNC             01/28/2020 - Modified CONVERTWM2TOUMOLM2S to 4.5
!=======================================================================

      REAL , PARAMETER :: CONVERTWM2TOUMOLM2S = 4.5 ,                   &
                          SOLARCONSTANT = 1367 ,                        &
                          WATERAIRRATIO = 18.016/28.97
                                              ! Ratio between water and air
                                              ! molecules
!     INTEGER,PARAMETER :: LAYERS =5          ! Number of layers in canopy model

! Canopy characteristics for MEGCAN canopy types
      INTEGER, PARAMETER ::  NRTYP = 6          ! Number of canopy types
      INTEGER , PARAMETER :: NRCHA = 17         ! Number of canopy characteristics

! 16 variables are assigned for each canopy type
! 1  = canopy depth
! 2  = leaf width
! 3  = leaf length
! 4  = canopy height
! 5  = scattering coefficient for PPFD
! 6  = scattering coefficient for near IR
! 7  = reflection coefficient for diffuse PPFD
! 8  = reflection coefficient for diffuse near IR
! 9  = clustering coefficient (accounts for leaf clumping influence on mean
!    projected leaf area in the direction of the suns beam)
! 10 = leaf IR emissivity
! 11 = leaf stomata and cuticle factor: 1=hypostomatous, 2=amphistomatous,
!     1.25=hypostomatous but with some transpiration through cuticle
! 12 = daytime temperature lapse rate (K m-1)
! 13 = nighttime temperature lapse rate (K m-1)
! 14 = warm (>283K) canopy total humidity change (Pa)
! 15 = cool (>= 283K) canopy total humidity change (Pa)
! 16 = normalized canopy depth where wind is negligible
! 17 = canopy transparency
!
! Six canopy types currently used in MEGCAN:
! 1  = Needleleaf trees
! 2  = Tropical forest trees,
! 3  = Temperate broadleaf trees
! 4  = shrubs
! 5  = herbaceous
! 6  = crops

      REAL , DIMENSION(NRCHA,NRTYP) :: Canopychar = reshape((/16.,16., &
                                      16.,1.,0.5,1.,0.005,0.05,0.05,   &
                                      0.015,0.01,0.02,0.1,0.1,0.1,0.1, &
                                      0.15,0.15,24.,24.,24.,2.,0.5,1.0,&
                                      0.2,0.2,0.2,0.2,0.2,0.2,0.8,0.8, &
                                      0.8,0.8,0.8,0.8,0.057,0.057,     &
                                      0.057,0.057,0.057,0.057,0.389,   &
                                      0.389,0.389,0.389,0.389,0.389,   &
                                      0.85,1.1,0.9,0.85,0.7,0.65,0.95, &
                                      0.95,0.95,0.95,0.95,0.95,1.25,   &
                                      1.25,1.25,1.,1.25,1.25,0.06,0.06,&
                                      0.06,0.06,0.06,0.06,-0.06,-0.06, &
                                      -0.06,-0.06,-0.06,-0.06,700.,    &
                                      700.,700.,700.,700.,700.,150.,   &
                                      150.,150.,150.,150.,150.,0.7,0.7,&
                                      0.7,0.7,0.7,0.7,0.2,0.2,0.2,0.2, &
                                      0.2,0.2/),SHAPE=(/NRCHA,NRTYP/), &
                                      ORDER=(/2,1/))


!=======================================================================

! FSB for Soil NOx 
! =======================================================================
!  MEGSEA.EXT
!  This include file contains wilting point information
!  for calculating soil moisture activity factor
!
! Created by Alex Guenther and Ling Huang in March 2017
!=======================================================================

      REAL,   PARAMETER  :: d1 = 0.04

!-- WWLT is wilting point (M^3/M^3) (JN90)

      REAL, PARAMETER  ::    WWLT(16) = (/               &        
                            0.068, 0.075, 0.114, 0.179,  &
                            0.155, 0.175, 0.218, 0.250,  &
                            0.219, 0.283, 0.286, 0.286,  &
                            0.286, 0.286, 0.286, 0.286   /)           

    ! FSB Based upon  MEGVEA.EXT includes suggestions from CJC

!=======================================================================
!  MEGVEA.EXT
!  This include file contains information required
!  for running MEGEAV module for calculating emission activity responses

!  Created by Alex Guenther and Ling Huang in Feb 2017
!
! 
!=======================================================================

    !Time step of LAI data
    REAL, PARAMETER :: TSTLEN  = 8.0

    !Number of emission classes
    INTEGER, PARAMETER :: NCLASS = 20
!   INTEGER, PARAMETER :: NCLASS  = NCLASS
    ! number of emission classes

    ! CO2 related emission activity factor parameters
    REAL,PARAMETER :: CO2   = 400.0
    REAL,PARAMETER :: ISmax =   1.344
    REAL,PARAMETER :: CO2h  =   1.4614
    REAL,PARAMETER :: Cstar = 585.0

    ! PSTD
    REAL,PARAMETER :: PSTD = 200

    ! canopy depth emission response
    REAL,PARAMETER :: CCD1 = -0.2
    REAL,PARAMETER :: CCD2 =  1.3

    !Light and temperature emission activity response coefficients for each emission class
    !LDF: light dependent fraction
    REAL           LDF(NCLASS)
    !CT1: temperature coefficient (emission type 1: light dependent)
    REAL           CT1(NCLASS)
    !Cleo: temperature coefficient (emission type 1: light dependent)
    REAL           Cleo(NCLASS)
    !beta: temperature coefficient (emission type 2: light independent)
    REAL           beta(NCLASS)

    DATA    beta(1),LDF(1),CT1(1),Cleo(1)        /  0.13,1.0,95,2  /
    DATA    beta(2),LDF(2),CT1(2),Cleo(2)        /  0.13,1.0,95,2  /
    DATA    beta(3),LDF(3),CT1(3),Cleo(3)        /  0.10,0.6,80,1.83  /
    DATA    beta(4),LDF(4),CT1(4),Cleo(4)        /  0.10,0.9,80,1.83  /
    DATA    beta(5),LDF(5),CT1(5),Cleo(5)        /  0.10,0.2,80,1.83  /
    DATA    beta(6),LDF(6),CT1(6),Cleo(6)        /  0.10,0.4,80,1.83  /
    DATA    beta(7),LDF(7),CT1(7),Cleo(7)        /  0.10,0.6,80,1.83  /
    DATA    beta(8),LDF(8),CT1(8),Cleo(8)        /  0.10,0.2,80,1.83  /
    DATA    beta(9),LDF(9),CT1(9),Cleo(9)        /  0.17,0.6,130,2.37  /
    DATA    beta(10),LDF(10),CT1(10),Cleo(10)    /  0.17,0.6,130,2.37  /
    DATA    beta(11),LDF(11),CT1(11),Cleo(11)    /  0.08,1.0,60,1.6  /
    DATA    beta(12),LDF(12),CT1(12),Cleo(12)    /  0.10,0.2,80,1.83  /
    DATA    beta(13),LDF(13),CT1(13),Cleo(13)    /  0.13,0.8,95,2  /
    DATA    beta(14),LDF(14),CT1(14),Cleo(14)    /  0.13,0.8,95,2  /
    DATA    beta(15),LDF(15),CT1(15),Cleo(15)    /  0.10,0.2,80,1.83  /
    DATA    beta(16),LDF(16),CT1(16),Cleo(16)    /  0.10,0.2,80,1.83  /
    DATA    beta(17),LDF(17),CT1(17),Cleo(17)    /  0.10,0.8,80,1.83  /
    DATA    beta(18),LDF(18),CT1(18),Cleo(18)    /  0.10,0.2,80,1.83  /
    DATA    beta(19),LDF(19),CT1(19),Cleo(19)    /  0.08,1.0,60,1.6  /
    DATA    beta(20),LDF(20),CT1(20),Cleo(20)    /  0.10,0,80,1.83  /


    ! Parameters for leaf age algorithm for each emission activity classes
      REAL           Anew(NCLASS)
      REAL           Agro(NCLASS)
      REAL           Amat(NCLASS)
      REAL           Aold(NCLASS)

      DATA    Anew(1),Agro(1),Amat(1),Aold(1)  / 0.05,0.6,1.0,0.9 /
      DATA    Anew(2),Agro(2),Amat(2),Aold(2)  / 0.05,0.6,1.0,0.9 /
      DATA    Anew(3),Agro(3),Amat(3),Aold(3)  / 2.0,1.8,1.0,1.05 /
      DATA    Anew(4),Agro(4),Amat(4),Aold(4)  / 2.0,1.8,1.0,1.05 /
      DATA    Anew(5),Agro(5),Amat(5),Aold(5)  / 2.0,1.8,1.0,1.05 /
      DATA    Anew(6),Agro(6),Amat(6),Aold(6)  / 2.0,1.8,1.0,1.05 /
      DATA    Anew(7),Agro(7),Amat(7),Aold(7)  / 2.0,1.8,1.0,1.05 /
      DATA    Anew(8),Agro(8),Amat(8),Aold(8)  / 2.0,1.8,1.0,1.05 /
      DATA    Anew(9),Agro(9),Amat(9),Aold(9)  / 0.4,0.6,1.0,0.95 /
      DATA    Anew(10),Agro(10),Amat(10),Aold(10) /0.4,0.6,1.0,0.95 /
      DATA    Anew(11),Agro(11),Amat(11),Aold(11) /3.5,3.0,1.0,1.2  /
      DATA    Anew(12),Agro(12),Amat(12),Aold(12) /1.0,1.0,1.0,1.0  /
      DATA    Anew(13),Agro(13),Amat(13),Aold(13) /1.0,1.0,1.0,1.0  /
      DATA    Anew(14),Agro(14),Amat(14),Aold(14) /1.0,1.0,1.0,1.0  /
      DATA    Anew(15),Agro(15),Amat(15),Aold(15) /1.0,1.0,1.0,1.0  /
      DATA    Anew(16),Agro(16),Amat(16),Aold(16) /1.0,1.0,1.0,1.0  /
      DATA    Anew(17),Agro(17),Amat(17),Aold(17) /1.0,1.0,1.0,1.0  /
      DATA    Anew(18),Agro(18),Amat(18),Aold(18) /1.0,1.0,1.0,1.0  /
      DATA    Anew(19),Agro(19),Amat(19),Aold(19) /1.0,1.0,1.0,1.0  /
      DATA    Anew(20),Agro(20),Amat(20),Aold(20) /1.0,1.0,1.0,1.0  /


    !stress emission activity response coefficients for each emission class
    !CAQ: coefficient for poor Air Quality stress
    REAL           CAQ(NCLASS)
    !CHW: coefficient for high wind speed stress
    REAL           CHW(NCLASS)
    !CHT: coefficient for high temperature stress
    REAL           CHT(NCLASS)
    !CLT: coefficient for high temperature stress
    REAL           CLT(NCLASS)

    DATA    CAQ(1),CHW(1),CHT(1),CLT(1)           /  1,1,1,1  /
    DATA    CAQ(2),CHW(2),CHT(2),CLT(2)           /  1,1,1,1  /
    DATA    CAQ(3),CHW(3),CHT(3),CLT(3)           /  1,5,1,1  /
    DATA    CAQ(4),CHW(4),CHT(4),CLT(4)           /  5,5,5,5  /
    DATA    CAQ(5),CHW(5),CHT(5),CLT(5)           /  1,5,1,1  /
    DATA    CAQ(6),CHW(6),CHT(6),CLT(6)           /  1,5,1,1  /
    DATA    CAQ(7),CHW(7),CHT(7),CLT(7)           /  1,5,1,1  /
    DATA    CAQ(8),CHW(8),CHT(8),CLT(8)           /  1,5,1,1  /
    DATA    CAQ(9),CHW(9),CHT(9),CLT(9)           /  5,5,5,5  /
    DATA    CAQ(10),CHW(10),CHT(10),CLT(10)       /  5,5,5,5  /
    DATA    CAQ(11),CHW(11),CHT(11),CLT(11)       /  1,1,1,1  /
    DATA    CAQ(12),CHW(12),CHT(12),CLT(12)       /  1,1,1,1  /
    DATA    CAQ(13),CHW(13),CHT(13),CLT(13)       /  1,1,1,1  /
    DATA    CAQ(14),CHW(14),CHT(14),CLT(14)       /  1,1,1,1  /
    DATA    CAQ(15),CHW(15),CHT(15),CLT(15)       /  1,1,1,1  /
    DATA    CAQ(16),CHW(16),CHT(16),CLT(16)       /  1,1,1,1  /
    DATA    CAQ(17),CHW(17),CHT(17),CLT(17)       /  5,5,5,5  /
    DATA    CAQ(18),CHW(18),CHT(18),CLT(18)       /  1,1,1,1  /
    DATA    CAQ(19),CHW(19),CHT(19),CLT(19)       /  1,1,1,1  /
    DATA    CAQ(20),CHW(20),CHT(20),CLT(20)       /  1,1,1,1  /


    !stress emission activity thresholds for each emission class
    !TAQ: threshold for poor Air Quality stress (ppm-hours)
    REAL           TAQ(NCLASS)
    !THW: threshold for high wind speed stress (m/s)
    REAL           THW(NCLASS)
    !THT: threshold for high temperature stress (Celsius degree)
    REAL           THT(NCLASS)
    !TLT: threshold for high temperature stress (Celsius degree)
    REAL           TLT(NCLASS)

    DATA    TAQ(1),THW(1),THT(1),TLT(1)           /  20,12,40,10  /
    DATA    TAQ(2),THW(2),THT(2),TLT(2)           /  20,12,40,10  /
    DATA    TAQ(3),THW(3),THT(3),TLT(3)           /  20,12,40,10  /
    DATA    TAQ(4),THW(4),THT(4),TLT(4)           /  20,12,40,10  /
    DATA    TAQ(5),THW(5),THT(5),TLT(5)           /  20,12,40,10  /
    DATA    TAQ(6),THW(6),THT(6),TLT(6)           /  20,12,40,10  /
    DATA    TAQ(7),THW(7),THT(7),TLT(7)           /  20,12,40,10  /
    DATA    TAQ(8),THW(8),THT(8),TLT(8)           /  20,12,40,10  /
    DATA    TAQ(9),THW(9),THT(9),TLT(9)           /  20,12,40,10  /
    DATA    TAQ(10),THW(10),THT(10),TLT(10)       /  20,12,40,10  /
    DATA    TAQ(11),THW(11),THT(11),TLT(11)       /  20,12,40,10  /
    DATA    TAQ(12),THW(12),THT(12),TLT(12)       /  20,12,40,10  /
    DATA    TAQ(13),THW(13),THT(13),TLT(13)       /  20,12,40,10  /
    DATA    TAQ(14),THW(14),THT(14),TLT(14)       /  20,12,40,10  /
    DATA    TAQ(15),THW(15),THT(15),TLT(15)       /  20,12,40,10  /
    DATA    TAQ(16),THW(16),THT(16),TLT(16)       /  20,12,40,10  /
    DATA    TAQ(17),THW(17),THT(17),TLT(17)       /  20,12,40,10  /
    DATA    TAQ(18),THW(18),THT(18),TLT(18)       /  20,12,40,10  /
    DATA    TAQ(19),THW(19),THT(19),TLT(19)       /  20,12,40,10  /
    DATA    TAQ(20),THW(20),THT(20),TLT(20)       /  20,12,40,10  /

    !stress emission activity delta thresholds for each emission class
    !DTAQ: delta threshold for poor Air Quality stress (ppm-hours)
    REAL           DTAQ(NCLASS)
    !DTHW: delta threshold for high wind speed stress (m/s)
    REAL           DTHW(NCLASS)
    !DTHT: delta threshold for high temperature stress (Celsius degree)
    REAL           DTHT(NCLASS)
    !DTLT: delta threshold for low temperature stress (Celsius degree)
    REAL           DTLT(NCLASS)

    DATA    DTAQ(1),DTHW(1),DTHT(1),DTLT(1)               /  30,8,8,8  /
    DATA    DTAQ(2),DTHW(2),DTHT(2),DTLT(2)               /  30,8,8,8  /
    DATA    DTAQ(3),DTHW(3),DTHT(3),DTLT(3)               /  30,8,8,8  /
    DATA    DTAQ(4),DTHW(4),DTHT(4),DTLT(4)               /  30,8,8,8  /
    DATA    DTAQ(5),DTHW(5),DTHT(5),DTLT(5)               /  30,8,8,8  /
    DATA    DTAQ(6),DTHW(6),DTHT(6),DTLT(6)               /  30,8,8,8  /
    DATA    DTAQ(7),DTHW(7),DTHT(7),DTLT(7)               /  30,8,8,8  /
    DATA    DTAQ(8),DTHW(8),DTHT(8),DTLT(8)               /  30,8,8,8  /
    DATA    DTAQ(9),DTHW(9),DTHT(9),DTLT(9)               /  30,8,8,8  /
    DATA    DTAQ(10),DTHW(10),DTHT(10),DTLT(10)           /  30,8,8,8  /
    DATA    DTAQ(11),DTHW(11),DTHT(11),DTLT(11)           /  30,8,8,8  /
    DATA    DTAQ(12),DTHW(12),DTHT(12),DTLT(12)           /  30,8,8,8  /
    DATA    DTAQ(13),DTHW(13),DTHT(13),DTLT(13)           /  30,8,8,8  /
    DATA    DTAQ(14),DTHW(14),DTHT(14),DTLT(14)           /  30,8,8,8  /
    DATA    DTAQ(15),DTHW(15),DTHT(15),DTLT(15)           /  30,8,8,8  /
    DATA    DTAQ(16),DTHW(16),DTHT(16),DTLT(16)           /  30,8,8,8  /
    DATA    DTAQ(17),DTHW(17),DTHT(17),DTLT(17)           /  30,8,8,8  /
    DATA    DTAQ(18),DTHW(18),DTHT(18),DTLT(18)           /  30,8,8,8  /
    DATA    DTAQ(19),DTHW(19),DTHT(19),DTLT(19)           /  30,8,8,8  /
    DATA    DTAQ(20),DTHW(20),DTHT(20),DTLT(20)           /  30,8,8,8  /

    ! MEGAN species
    ! Based on Alex Guenther's "Description Class.xlsx" for MEGANv3

    CHARACTER*16, PARAMETER :: MGN_SPC(NCLASS) =    &
    (/  'ISOP            ',                         &      !  1 isoprene
        'MBO             ',                         &      !  2 MBO
        'MT_PINE         ',                         &      !  3 monoterpenes: pines (alpha and beta)
        'MT_ACYC         ',                         &      !  4 monoterpenes: acyclic, 3 = (e.g., myrcene, ocimenes)
        'MT_CAMP         ',                         &      !  5 monoterpenes: carene, camphene, others
        'MT_SABI         ',                         &      !  6 monoterpenes: sabinene, limonene, terpinenes, others
        'MT_AROM         ',                         &      !  7  aromatic: cymenes, cymenenes
        'MT_OXY          ',                         &      !  8 C8-C13 oxygenated (e.g., camphor)
        'SQT_HR          ',                         &      !  9 Highly reactive SQT (e.g., caryophyllene)
        'SQT_LR          ',                         &      ! 10 less reactive SQT  (e.g., longifolene, copaene) and salates
        'MEOH            ',                         &      ! 11 methanol
        'ACTO            ',                         &      ! 12 acetone
        'ETOH            ',                         &      ! 13 acetaldehyde and ethanol
        'ACID            ',                         &      ! 14 organic acids: formic acid, acetic acid, pyruvic acid
        'LVOC            ',                         &      ! 15 C2 to C4 HC (e.g., ethene, ethane)
        'OXPROD          ',                         &      ! 16 oxidation products: aldehydes
        'STRESS          ',                         &      ! 17 Stress compounds (e.g., linalool)
        'OTHER           ',                         &      ! 18 other VOC (e.g., indole, pentane, methyl bromide)
        'CO              ',                         &      ! 19 carbon monoxide
        'NO              '  /)                             ! 20 Nitric oxide

end module MEG_EXT