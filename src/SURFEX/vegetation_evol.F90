!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!     #########
    SUBROUTINE VEGETATION_EVOL(IO, DTI, PK, PEK, DEK, OAGRIP, PTSTEP, KMONTH, KDAY,    &
                               PTIME, PLAT, PCO2, ISSK, PRESP_BIOMASS_INST, OBIOM_REAP,&
                               OIRRIG_ONLY, OTEB_HVEG_ONLY                             )  
!   ###############################################################
!!****  *VEGETATION EVOL*
!!
!!    PURPOSE
!!    -------
!
!     performs the time evolution of vegetation parameters
!     at solar midnight in the case of interactive vegetation (ISBA-Ags)
!              
!!**  METHOD
!!    ------
!!
!!    EXTERNAL
!!    --------
!!    none
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      
!!    none
!!
!!    REFERENCE
!!    ---------
!!
!!      
!!    AUTHOR
!!    ------
!!
!!      V. Masson          * Meteo-France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    01/03/03 
!!      P. Le Moigne 12/2004 : NIT version 
!!      P Le Moigne  09/2005 : AGS modifs of L. Jarlan
!!      A.L. Gibelin 04/2009 : BIOMASS and RESP_BIOMASS arrays
!!      A.L. Gibelin 04/2009 : Add NCB option 
!!      D. Carrer    01/2012 : representation of nitrogen dilution fct of CO2 (from Calvet et al. 2008)
!!      B. Decharme  05/2012 : Optimization and ISBA-DIF coupling
!!      C. Delire    01/2014 : IBIS respiration for tropical evergreen
!!      R. Seferian  05/2015 : expanding of Nitrogen dilution option to the complete formulation proposed by Yin et al. GCB 2002 
!!Seferian & Delire  06/2015 : accouting for living woody biomass respiration (expanding work of E Joetzjer to all woody PFTs) 
!!      B. Decharme  01/2016 : Bug when vegetation veg, z0 and emis are imposed whith interactive vegetation
!!      A. Druel     02/2019 : Streamlines the code and adapt it to be compatible with new irrigation and multi-season
!!      R. Seferian  08/2016 : Remove dependency to SW in Yin ; Crop harvesting
!!      B. Decharme  01/2018 : Compute H_VEG
!!
!-------------------------------------------------------------------------------
!
!*       0.     DECLARATIONS
!               ------------
!
USE MODD_ISBA_OPTIONS_n,   ONLY : ISBA_OPTIONS_t
USE MODD_DATA_ISBA_n,      ONLY : DATA_ISBA_t
USE MODD_ISBA_n,           ONLY : ISBA_P_t, ISBA_PE_t
USE MODD_DIAG_EVAP_ISBA_n, ONLY : DIAG_EVAP_ISBA_t
!
USE MODD_SSO_n,          ONLY : SSO_t
!
USE MODD_CO2V_PAR,       ONLY : XMC, XMCO2, XPCCO2, XRESPFACTOR_NIT,       &
                                XCOEFF_MAINT_RESP_ZERO, XSLOPE_MAINT_RESP, &
                                XDILUDEC, ITRANSFERT_ESG, XGTOKG
USE MODD_CSTS,           ONLY : XDAY, XTT, XMD
!
USE MODI_ALBEDO
USE MODI_LAIGAIN
USE MODI_LAILOSS
USE MODI_NITRO_DECLINE
USE MODI_EMIS_FROM_VEG
USE MODI_VEG_FROM_LAI
USE MODI_Z0V_FROM_LAI
USE MODI_H_VEG_FROM_LAI
USE MODI_SUBSCALE_Z0EFF
USE MODD_TYPE_DATE_SURF
USE MODD_DATA_COVER_PAR, ONLY : NVEGTYPE_ECOSG, NVEGTYPE, &
                                NVT_TEBD, NVT_TRBE, NVT_BONE,   &
                                NVT_TRBD, NVT_TEBE, NVT_TENE,   &
                                NVT_BOBD, NVT_BOND, NVT_SHRB,   &
                                NVT_TRBE, NVT_C3,   NVT_C4,     &
                                NVT_IRR
USE MODD_AGRI,           ONLY : NVEG_IRR, LIRRIGMODE, LMULTI_SEASON
USE MODD_SURF_PAR
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
TYPE(ISBA_OPTIONS_t),   INTENT(INOUT) :: IO
TYPE(DATA_ISBA_t),      INTENT(INOUT) :: DTI
TYPE(ISBA_P_t),         INTENT(INOUT) :: PK
TYPE(ISBA_PE_t),        INTENT(INOUT) :: PEK
TYPE(DIAG_EVAP_ISBA_t), INTENT(INOUT) :: DEK
TYPE(SSO_t),            INTENT(INOUT) :: ISSK
!
LOGICAL,              INTENT(IN)    :: OAGRIP  ! agricultural practices
!
REAL,                 INTENT(IN)    :: PTSTEP  ! time step
INTEGER,              INTENT(IN)    :: KMONTH  ! current month
INTEGER,              INTENT(IN)    :: KDAY    ! current day
REAL,                 INTENT(IN)    :: PTIME   ! current time since midnight
REAL,   DIMENSION(:), INTENT(IN)    :: PLAT    ! latitude of each grid point
REAL,   DIMENSION(:), INTENT(IN)    :: PCO2    ! CO2 concentration [ppmm]
!
REAL, DIMENSION(:,:), INTENT(INOUT) :: PRESP_BIOMASS_INST   ! instantaneous respiration of biomass (kgCO2/m2/s)
!
LOGICAL,              INTENT(IN)    :: OBIOM_REAP ! .TRUE. + OAGRIP = Biomass harvested
!
LOGICAL,              INTENT(IN),   OPTIONAL :: OIRRIG_ONLY    ! To only done the multi-irrigation computation
LOGICAL,                            OPTIONAL :: OTEB_HVEG_ONLY ! flag for TEB high vegetation computations only
!
!*      0.2    declarations of local parameter
!
REAL, PARAMETER                   :: ZCOEF1 = 10.0
REAL, PARAMETER                   :: ZCOEF2 = 25.0
REAL, PARAMETER                   :: ZDEPTH = 1.0   !Temp depth m
!
REAL, PARAMETER                   :: ZWOOD_IBIS=0.0125
REAL, PARAMETER                   :: ZROOT_IBIS=1.25 
REAL, PARAMETER                   :: ZCIBIS1   =3500.
REAL, PARAMETER                   :: ZCIBIS2   =1./288.
REAL, PARAMETER                   :: ZNDAY     =365.
!
REAL, PARAMETER                   :: ZCDILU2 = 6.3
REAL, PARAMETER                   :: ZCDILU3 = 288.
!
REAL, PARAMETER                   :: ZDEPTH_VEG = 0.40    !Depth in meters for daily temperature
REAL, PARAMETER                   :: ZTEMP_VEG  = 23.     !Average temperature of the vegetation
REAL, PARAMETER                   :: ZTMAX_VEG  = 33.     !Maximum temperature of the vegetation used for the Yin parameterization
REAL, PARAMETER                   :: ZDECIDUS   = 0.75    !Coef for decidus trees
REAL, PARAMETER                   :: ZCOEF      = 0.33    !Coef for Yin et al., nitrogen dilu param
!
!*      0.3    declarations of local variables
!
REAL, DIMENSION(SIZE(PEK%XRESP_BIOMASS,1),SIZE(PEK%XRESP_BIOMASS,2)) :: ZRESP_BIOMASS_LAST ! biomass at t-1 (kg_DM/m2/day)
!
REAL, DIMENSION(SIZE(PEK%XTG,1)) :: ZBIOMASS_LEAF ! temporary leaf biomass 
REAL, DIMENSION(SIZE(PEK%XTG,1)) :: ZBSLAI_NITRO  ! (Calvet et al. 2008) ratio of biomass to LAI with representation of nitrogen dilution
REAL, DIMENSION(SIZE(PEK%XTG,1)) :: ZCO2          ! CO2 concentration in the air (ppm)
REAL, DIMENSION(SIZE(PEK%XTG,1)) :: ZCNA_NITRO    ! fct of CO2        
REAL, DIMENSION(SIZE(PEK%XTG,1)) :: ZPARAM
REAL, DIMENSION(SIZE(PEK%XTG,1)) :: ZHTREE        ! tree height used for estimation of sapwood fraction
REAL, DIMENSION(SIZE(PEK%XTG,1)) :: ZSAPFRAC      ! tree sap fraction used for estimation of sapwood fraction
REAL, DIMENSION(SIZE(PEK%XTG,1)) :: ZTG_VEG       ! surface temperature   (C)
REAL, DIMENSION(SIZE(PEK%XTG,1)) :: ZTG_SOIL      ! soil temperature   (C)
REAL, DIMENSION(SIZE(PEK%XTG,1)) :: ZDG_SOIL      ! soil depth for DIF (m)
REAL, DIMENSION(SIZE(PEK%XTG,1)) :: ZVEG_TREES
REAL, DIMENSION(SIZE(PEK%XTG,1)) :: ZVEG_CROP
REAL, DIMENSION(SIZE(PEK%XTG,1)) :: ZFERT         ! Azote
!
LOGICAL, DIMENSION(SIZE(PEK%XTG,1)) :: GWOOD,GHERB,GCROP
!
LOGICAL :: GMASK, GTEB_HVEG_ONLY, GIRRIG_ONLY
!
INTEGER :: INI, INL, JI, JL, JTYPE, JTYPE2
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-----------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('VEGETATION_EVOL',0,ZHOOK_HANDLE)
!
!*      1.     Preliminaries
!              -------------
!
INI=SIZE(PEK%XTG,1)
INL=SIZE(PEK%XTG,2)
!
GIRRIG_ONLY=.FALSE.
IF(PRESENT(OIRRIG_ONLY))THEN
  GIRRIG_ONLY=OIRRIG_ONLY
ENDIF
!
GTEB_HVEG_ONLY=.FALSE.
IF(PRESENT(OTEB_HVEG_ONLY))THEN
  GTEB_HVEG_ONLY = OTEB_HVEG_ONLY
ENDIF
!
! Mask where vegetation evolution is performed (just before solar midnight)
GMASK = ( PTIME - PTSTEP < 0. ) .AND. ( PTIME >= 0. )
!
!-----------------------------------------------------------------
!
!*      2.     Interactive vegetation
!              ----------------------
!
IF (.NOT.GIRRIG_ONLY) THEN
   !
   ! Define herbaceous and woody patches, coniferous and deciduous woody species (BOND excluded)
   !
   ZVEG_TREES(:)=0.
   ZVEG_CROP (:)=0.
   DO JTYPE = 1, SIZE(PK%XVEGTYPE_PATCH,2)
      !
      JTYPE2 = JTYPE
      !
      IF (JTYPE > NVEGTYPE) THEN
         JTYPE2 = DTI%NPAR_VEG_IRR_USE(JTYPE-NVEGTYPE)
      ENDIF
      !
      DO JI=1,INI
         !
         IF (JTYPE2==NVT_TEBD .OR. JTYPE2==NVT_TRBE .OR. JTYPE2==NVT_BONE .OR. &
             JTYPE2==NVT_TRBD .OR. JTYPE2==NVT_TEBE .OR. JTYPE2==NVT_TENE .OR. &
             JTYPE2==NVT_BOBD .OR. JTYPE2==NVT_BOND .OR. JTYPE2==NVT_SHRB      ) THEN
            ZVEG_TREES(JI) = ZVEG_TREES(JI) + PK%XVEGTYPE_PATCH(JI,JTYPE)
         ENDIF
         !
         IF (JTYPE2==NVT_C3 .OR. JTYPE2==NVT_C4 .OR. JTYPE2==NVT_IRR) THEN
            ZVEG_CROP(JI) = ZVEG_CROP(JI) + PK%XVEGTYPE_PATCH(JI,JTYPE)
         ENDIF
         !
      ENDDO
   ENDDO
   GHERB(:) = (ZVEG_TREES(:)<0.5)
   GWOOD(:) = (.NOT.GHERB (:))
   GCROP(:) = (ZVEG_CROP(:)>=0.5)
   !
   ! * Save RESP_BIOMASS at t-1
   !
   IF (GMASK) THEN
     PEK%XRESP_BIOMASS (:,1) = 0.0
     ZRESP_BIOMASS_LAST(:,:) = 0.0
   ELSE
     PEK%XRESP_BIOMASS (:,1) = PEK%XRESP_BIOMASS(:,1) + PRESP_BIOMASS_INST(:,1) * (PTSTEP*XMC/(XPCCO2*XMCO2))
     ZRESP_BIOMASS_LAST(:,:) = PEK%XRESP_BIOMASS(:,:)
   ENDIF
   !
   ! * Comput interactive vegetation
   !
   CALL INTERACTIVE_VEGETATION
   !
   ! * Instantaneous respiration (kgCO2/m2/s)
   !
   DO JL=2,SIZE(PEK%XRESP_BIOMASS(:,:),2)
      DO JI=1,INI
         PRESP_BIOMASS_INST(JI,JL) = (PEK%XRESP_BIOMASS(JI,JL)-ZRESP_BIOMASS_LAST(JI,JL)) * (XPCCO2*XMCO2/(PTSTEP*XMC))
      ENDDO
   ENDDO
   !
   ! * Simple harvesting for ISBA-CC (kgC/m2/s)
   !
   IF (GMASK.AND.IO%CPHOTO=='NCB'.AND.IO%LLULCC) THEN
      DEK%XFHARVEST(:) = 0.
      WHERE(GCROP(:)) 
            DEK%XFHARVEST(:) = (PK%XTURNOVER(:,1) + PK%XTURNOVER(:,2) + PK%XTURNOVER(:,3) + PK%XTURNOVER(:,5)) * XGTOKG
            ! reset turnover since AGB has been harvested
            PK%XTURNOVER(:,1) = 0.
            PK%XTURNOVER(:,2) = 0.
            PK%XTURNOVER(:,3) = 0.
            PK%XTURNOVER(:,5) = 0.
      ENDWHERE
   ENDIF
   !
ENDIF
!
!-----------------------------------------------------------------
!
!*      3.     Agricultural or irrigation practices
!              ------------------------------------
!
IF (OAGRIP .OR. (LIRRIGMODE .AND. LMULTI_SEASON) ) THEN
   CALL AGRI_OR_IRRIG
ENDIF
!
!-----------------------------------------------------------------
!
!*      4.     Physical parameters depending on vegetation
!              -------------------------------------------
!
IF (GMASK .AND. (.NOT.GIRRIG_ONLY) ) THEN
  !
  ! Evolution of vegetation fraction and roughness length due to LAI change
  !
  IF(.NOT.DTI%LIMP_Z0) THEN
    WHERE( PEK%XVEG(:) > 0. ) 
      PEK%XZ0 (:) = Z0V_FROM_LAI(PEK%XLAI(:),PK%XH_TREE(:),PK%XVEGTYPE_PATCH(:,:),DTI%NPAR_VEG_IRR_USE)
    ENDWHERE
  ENDIF
  IF(.NOT.DTI%LIMP_VEG .AND. (.NOT.GTEB_HVEG_ONLY)) THEN
    WHERE( PEK%XVEG(:) > 0. )
      PEK%XVEG(:) = VEG_FROM_LAI(PEK%XLAI(:),PK%XVEGTYPE_PATCH(:,:),DTI%NPAR_VEG_IRR_USE)
    ENDWHERE
  ENDIF
  IF(.NOT.DTI%LIMP_H_VEG .AND. (.NOT.GTEB_HVEG_ONLY)) THEN
    WHERE( PEK%XVEG(:) > 0. )
      PEK%XH_VEG(:) = H_VEG_FROM_LAI(PEK%XLAI(:),PK%XH_TREE(:),PK%XVEGTYPE_PATCH(:,:),DTI%NPAR_VEG_IRR_USE)
    ENDWHERE
  ENDIF
  !
  ! Evolution of radiative parameters due to vegetation fraction change
  !
  IF(.NOT.DTI%LIMP_EMIS) THEN
    WHERE( PEK%XVEG(:) > 0. ) PEK%XEMIS(:)= EMIS_FROM_VEG(PEK%XVEG(:),PK%XVEGTYPE_PATCH(:,:),DTI%NPAR_VEG_IRR_USE)
  ENDIF
  !
  IF (.NOT.GTEB_HVEG_ONLY) THEN
     CALL ALBEDO(IO%CALBEDO,PEK)
  ENDIF
  !
  ! Evolution of effective roughness length due to new surface roughness length
  !
  IF (ASSOCIATED(ISSK%XAOSIP)) THEN
    IF (SIZE(ISSK%XAOSIP)>0) THEN
      CALL SUBSCALE_Z0EFF(ISSK,PEK%XZ0(:),.FALSE. )
    ENDIF
  ENDIF
  !
ENDIF
!
!
!-----------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('VEGETATION_EVOL',1,ZHOOK_HANDLE)
!
!
!-------------------------------------------------------------------------------
CONTAINS
!-------------------------------------------------------------------------------
!
SUBROUTINE INTERACTIVE_VEGETATION
!
IMPLICIT NONE
!
REAL :: ZWGHT_SOIL, ZPARAM_TYPE, &
        ZLOG2, ZDILUDEC 
!
INTEGER :: IDEPTH

REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!
IF (LHOOK) CALL DR_HOOK('VEGETATION_EVOL:INTERACTIVE_VEGETATION',0,ZHOOK_HANDLE)
!
ZLOG2 = LOG(2.0)
!
ZTG_SOIL(:) = 0.0
ZTG_VEG (:) = 0.0
!
!  LAI daily mortality and assimilation
!
ZBIOMASS_LEAF(:) = PEK%XBIOMASS(:,1)
!
IF (GMASK) THEN
  !        
  PK%XINCREASE(:,:) = 0.0
  PK%XTURNOVER(:,:) = 0.0
  ZBSLAI_NITRO(:  ) = PK%XBSLAI_NITRO(:) 
  !
  IF(IO%LNITRO_DILU)THEN
    !
    !* Compute Vegetation temperature
    !  We use the temperature of the second layer of the soil (<40cm)
    !  since the parametrization employs a daily temperature
    !
    IF(IO%CISBA/='DIF')THEN        
      ZTG_VEG(:) = PEK%XTG(:,2)
    ELSE 
      DO JI=1,INI
         IDEPTH=PK%NWG_LAYER(JI)
         ZDG_SOIL(JI)=MIN(ZDEPTH_VEG,PK%XDG(JI,IDEPTH))
      ENDDO  
      DO JL=1,INL
        DO JI=1,INI     
          ZWGHT_SOIL=MIN(PK%XDZG(JI,JL),MAX(0.0,ZDG_SOIL(JI)-PK%XDG(JI,JL)+PK%XDZG(JI,JL)))        
          ZTG_VEG(JI)=ZTG_VEG(JI)+PEK%XTG(JI,JL)*ZWGHT_SOIL/ZDG_SOIL(JI)
        ENDDO
      ENDDO 
    ENDIF
    !
    !* set carbon-nitrogen limitation
    !
    ZPARAM(:) = 0.0
    ZFERT (:) = 0.0   ! #strange: Oo All time = 0 ? never implemented ?
    !
    DO JTYPE=1,SIZE(PK%XVEGTYPE_PATCH,2)
       !
       JTYPE2 = JTYPE
       !
       IF (JTYPE > NVEGTYPE) THEN
          JTYPE2 = DTI%NPAR_VEG_IRR_USE( JTYPE - NVEGTYPE )
       ENDIF
       !
       IF (NVEGTYPE==NVEGTYPE_ECOSG) THEN
         ZDILUDEC = XDILUDEC(ITRANSFERT_ESG(JTYPE2))
       ELSE
         ZDILUDEC = XDILUDEC(JTYPE2)
       ENDIF
       !
       DO JI = 1,INI
          ZPARAM_TYPE =        ZDILUDEC  * (ZDECIDUS +  (MIN(ZTG_VEG(JI)-XTT,ZTMAX_VEG))/ZTEMP_VEG - ZCOEF * ZFERT(JI)) &
                      + (1.0 - ZDILUDEC) * (            (MIN(ZTG_VEG(JI)-XTT,ZTMAX_VEG))/ZTEMP_VEG - ZCOEF * ZFERT(JI))
          ZPARAM(JI) = ZPARAM(JI) + ZPARAM_TYPE * PK%XVEGTYPE_PATCH(JI,JTYPE)
       ENDDO
       !
    ENDDO  
    !
    WHERE((PEK%XCE_NITRO(:)*PEK%XCNA_NITRO(:)+PEK%XCF_NITRO(:))/=0.0.AND.PEK%XCNA_NITRO(:)/=0.0)
      ZCO2        (:) = PCO2(:)*(XMD/(1.E-6*XMCO2))  ! (kg/kg ->  ppm)
      ZCNA_NITRO  (:) = PEK%XCNA_NITRO(:) * &
                   EXP(IO%XCNLIM*EXP(ZPARAM(:)-PEK%XCNA_NITRO(:)/ZCDILU2) * ALOG(MAX(1.,ZCO2(:)/ZCDILU3)))
      ZBSLAI_NITRO(:) = 1. / (PEK%XCE_NITRO(:)*ZCNA_NITRO(:)+PEK%XCF_NITRO(:))
    ENDWHERE
    !
  ENDIF
  !
  IF(ANY(PEK%XLAI(:)/=XUNDEF))THEN
    CALL NITRO_DECLINE(IO, PK, PEK, GWOOD, ZBSLAI_NITRO, PLAT, ZBIOMASS_LEAF)
    CALL LAIGAIN(ZBSLAI_NITRO, PEK, ZBIOMASS_LEAF)
  ENDIF
  !  
  ! reinitialise  PEK%XANDAY(:) and PEK%XANFM(:) 
  PEK%XANDAY(:)=0.0
  PEK%XANFM(:) =0.0
  !
ENDIF
!
!
! * soil temperature in K (over 1m depth for DIF)
!
ZTG_VEG(:) = PEK%XTG(:,1)
!
IF(IO%CISBA/='DIF')THEN        
  ZTG_SOIL(:) = PEK%XTG(:,2)
ELSE       
  DO JI=1,INI
    IDEPTH=PK%NWG_LAYER(JI)
    ZDG_SOIL(JI)=MIN(ZDEPTH,PK%XDG(JI,IDEPTH))
  ENDDO  
  DO JL=1,INL
    DO JI=1,INI     
      ZWGHT_SOIL=MIN(PK%XDZG(JI,JL),MAX(0.0,ZDG_SOIL(JI)-PK%XDG(JI,JL)+PK%XDZG(JI,JL)))        
      ZTG_SOIL(JI)=ZTG_SOIL(JI)+PEK%XTG(JI,JL)*ZWGHT_SOIL/ZDG_SOIL(JI)
    ENDDO
  ENDDO 
ENDIF
!
!
! * Respiration of structural biomass pools
!
WHERE(GWOOD(:))
  ! IBIS respiration with either respiration factor rwood=0.0125 - otherwise rroot=1.25 
  ! (Kucharik et al, 2000, eq 6-8) Soil temp in K         
  PEK%XRESP_BIOMASS(:,2) = PEK%XRESP_BIOMASS(:,2) + PEK%XBIOMASS(:,2) * PTSTEP &
                              * MAX(0.,ZROOT_IBIS*EXP(ZCIBIS1*(ZCIBIS2-1./ZTG_VEG(:)))/(ZNDAY*XDAY)) 
ELSEWHERE 
  PEK%XRESP_BIOMASS(:,2) = PEK%XRESP_BIOMASS(:,2) + PEK%XBIOMASS(:,2) * XRESPFACTOR_NIT    &
                              * EXP((ZLOG2/ZCOEF1)*(ZTG_VEG(:)-XTT-ZCOEF2)) * PTSTEP  
  ! before optimization                   * 2.0**((PEK%XTG(:,2)-XTT-ZCOEF2)/ZCOEF1) * PTSTEP               
ENDWHERE
!
IF (IO%CPHOTO == 'NIT') THEN
  !
  PEK%XRESP_BIOMASS(:,3) = PEK%XRESP_BIOMASS(:,3) + PEK%XBIOMASS(:,3) * XRESPFACTOR_NIT &
                            * EXP((ZLOG2/ZCOEF1)*(ZTG_SOIL(:)-XTT-ZCOEF2)) * PTSTEP  
  ! before optimization                   * 2.0**((PEK%XTG(:,2)-XTT-ZCOEF2)/ZCOEF1) * PTSTEP               
  !
ELSEIF (IO%CPHOTO == 'NCB') THEN
  !
  PEK%XRESP_BIOMASS(:,2) = MIN(PEK%XRESP_BIOMASS(:,2), PEK%XBIOMASS(:,2))
  ! 
  PEK%XRESP_BIOMASS(:,3) = PEK%XRESP_BIOMASS(:,3) + PEK%XBIOMASS(:,3) * &
            MAX( 0., XCOEFF_MAINT_RESP_ZERO * (1. + XSLOPE_MAINT_RESP*(ZTG_VEG(:)-XTT))) * PTSTEP  
  PEK%XRESP_BIOMASS(:,3) = MIN(PEK%XRESP_BIOMASS(:,3), PEK%XBIOMASS(:,3))
  ! 
  WHERE(GWOOD(:))
    ! Resp IBIS (Soil temp in K)
    PEK%XRESP_BIOMASS(:,4) = PEK%XRESP_BIOMASS(:,4) + PEK%XBIOMASS(:,4) * PTSTEP &
                        * MAX(0.,ZROOT_IBIS * EXP(ZCIBIS1*(ZCIBIS2-1./ZTG_SOIL(:)))/(ZNDAY*XDAY))
  ELSEWHERE 
    PEK%XRESP_BIOMASS(:,4) = PEK%XRESP_BIOMASS(:,4) + PEK%XBIOMASS(:,4) * &
             MAX( 0., XCOEFF_MAINT_RESP_ZERO * (1. + XSLOPE_MAINT_RESP*(ZTG_SOIL(:)-XTT))) * PTSTEP  
  ENDWHERE
  !
  PEK%XRESP_BIOMASS(:,4) = MIN(PEK%XRESP_BIOMASS(:,4), PEK%XBIOMASS(:,4))
  !
  WHERE( (GWOOD(:)).AND.(PEK%XBIOMASS(:,5)>0.) )
    ! IBIS estimation of sapwood fraction based on the height of tree, sapspeed and 
    ! max transpiration rates. Conversion from DM to C. To be changed with DGVM.  (Soil temp in K)        
    ZHTREE(:) = 2.5*0.75*(PEK%XBIOMASS(:,1)+PEK%XBIOMASS(:,2)+PEK%XBIOMASS(:,3)+&
                          PEK%XBIOMASS(:,4)+PEK%XBIOMASS(:,5)+PEK%XBIOMASS(:,6))*0.4
    ZSAPFRAC(:) = MIN(0.5, MAX(0.05,0.0025/25.*ZHTREE(:)*0.75*400/(PEK%XBIOMASS(:,5)*0.4)))
    !ZSAPFRAC(:) = 0.5
    
    PEK%XRESP_BIOMASS(:,5) = PEK%XRESP_BIOMASS(:,5) + PEK%XBIOMASS(:,5) * ZSAPFRAC(:) * PTSTEP &
                               * MAX(0.,ZWOOD_IBIS*EXP(ZCIBIS1*(ZCIBIS2-1./ZTG_VEG(:)))/(ZNDAY*XDAY))
    PEK%XRESP_BIOMASS(:,5) = MIN(PEK%XRESP_BIOMASS(:,5), PEK%XBIOMASS(:,5))
  ELSEWHERE
    PEK%XRESP_BIOMASS(:,5) = 0.0
  ENDWHERE
  !
ENDIF
!
IF (LHOOK) CALL DR_HOOK('VEGETATION_EVOL:INTERACTIVE_VEGETATION',1,ZHOOK_HANDLE)
!
END SUBROUTINE INTERACTIVE_VEGETATION
!
!
!-------------------------------------------------------------------------------
!
!
SUBROUTINE AGRI_OR_IRRIG
!
IMPLICIT NONE
!
LOGICAL, DIMENSION(SIZE(PEK%XLAI,1)) :: GINV       ! Check seed/reap inversion (T/F)
LOGICAL, DIMENSION(SIZE(PEK%XLAI,1)) :: GMASK_AGRI
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
IF (LHOOK) CALL DR_HOOK('VEGETATION_EVOL:AGRI_OR_IRRIG',0,ZHOOK_HANDLE)
!
! 1 Determining if there are inversions in seeding and reaping date (=winter crops) 
! -----------------------------------------------------------------------------------
!
GINV(:) = .FALSE.
WHERE ( PEK%TSEED(:)%TDATE%MONTH /= NUNDEF .AND. PEK%TREAP(:)%TDATE%MONTH /= NUNDEF .AND.  &
      ( PEK%TSEED(:)%TDATE%MONTH >  PEK%TREAP(:)%TDATE%MONTH .OR.                          &
        ( PEK%TSEED(:)%TDATE%MONTH==PEK%TREAP(:)%TDATE%MONTH .AND. PEK%TSEED(:)%TDATE%DAY>PEK%TREAP(:)%TDATE%DAY ) ) )
  GINV(:) = .TRUE.
ENDWHERE
!
!
! 2 Define the GMASK_AGRI for winter crops and summer crops 
! -----------------------------------------------------------
!
GMASK_AGRI(:) = .FALSE.
WHERE ( ( PEK%TSEED(:)%TDATE%MONTH /= NUNDEF .AND. PEK%TREAP(:)%TDATE%MONTH /= NUNDEF .AND.       & ! TSEED AND TREAD defines
        ( ( (.NOT.GINV(:)) .AND.                                                                  & ! If summer crop:
            ( ( (KMONTH == PEK%TSEED(:)%TDATE%MONTH .AND. KDAY < PEK%TSEED(:)%TDATE%DAY) .OR.     & !        - Before seeding
            KMONTH < PEK%TSEED(:)%TDATE%MONTH ) .OR. ( KMONTH > PEK%TREAP(:)%TDATE%MONTH .OR.     & !        - OR after reaping
            (KMONTH == PEK%TREAP(:)%TDATE%MONTH .AND. KDAY >= PEK%TREAP(:)%TDATE%DAY) ) ) )       &
          .OR. ( (GINV(:) ) .AND.                                                                 & ! If winter crop:
            ( ( (KMONTH == PEK%TSEED(:)%TDATE%MONTH .AND. KDAY < PEK%TSEED(:)%TDATE%DAY) .OR.     & !        - Before seeding
            KMONTH < PEK%TSEED(:)%TDATE%MONTH ) .AND. (KMONTH > PEK%TREAP(:)%TDATE%MONTH .OR.     & !        - AND after reaping
            (KMONTH == PEK%TREAP(:)%TDATE%MONTH .AND. KDAY >= PEK%TREAP(:)%TDATE%DAY) ) ) ) ) )   & 
       .AND. (PEK%TSEED(:)%TDATE%MONTH /= 1  .OR. PEK%TSEED(:)%TDATE%DAY /= 1 .OR.                & ! Ajout de l'exeption si toute l'année: pas de coupe !!
              PEK%TREAP(:)%TDATE%MONTH /= 12 .OR. PEK%TREAP(:)%TDATE%DAY /= 31) )
  GMASK_AGRI(:) = .TRUE.
ENDWHERE
!
! 3 Remove crops if not allowed (exept for trees)
! -----------------------------------------------------------
!
IF ( OAGRIP .AND. OBIOM_REAP ) THEN
   !
   WHERE (GMASK_AGRI(:))
      PEK%XLAI(:)         = PEK%XLAIMIN(:)
      ZBIOMASS_LEAF(:)    = PEK%XLAI(:) * ZBSLAI_NITRO(:)
   ENDWHERE
   !
   WHERE (GMASK_AGRI(:))
      PEK%XBIOMASS(:,1)       = 0.0
      PEK%XBIOMASS(:,2)       = 0.0
      PEK%XBIOMASS(:,3)       = 0.0
      PEK%XRESP_BIOMASS(:,2)  = 0.0
      PEK%XRESP_BIOMASS(:,3)  = 0.0
   ENDWHERE
   !
   IF (IO%CPHOTO == 'NCB') THEN
      WHERE (GMASK_AGRI(:)) 
        PEK%XBIOMASS(:,4)       = 0.0
        PEK%XBIOMASS(:,5)       = 0.0
        PEK%XBIOMASS(:,6)       = 0.0
        PEK%XRESP_BIOMASS(:,4)  = 0.0
      ENDWHERE
   ENDIF
   !
ENDIF
!
! 4 Check if there is more than one season for vegetation type and update the dates
! -----------------------------------------------------------
!
IF ( LMULTI_SEASON ) THEN
  !
  ! Update the date only when we are outside a irrigation season and only if there is differents season in this point
  !
  WHERE (GMASK_AGRI(:) .AND. PEK%MULTI_TSEED(:,2)%TDATE%MONTH /= NUNDEF )
      !
      ! Check if we are still in the prescendent season dates for oagrip/oirrigation
      !
      WHERE ( ( PEK%TREAP(:)%TDATE%MONTH <   KMONTH                                          .OR. &
               (PEK%TREAP(:)%TDATE%MONTH ==  KMONTH .AND. PEK%TREAP(:)%TDATE%DAY <= KDAY) )  .AND.& ! wrong season if after treap
              (.NOT.GINV(:).AND. .NOT.( PEK%MULTI_TREAP(:,1)%TDATE%MONTH==PEK%TREAP(:)%TDATE%MONTH .AND. & ! IF no inversion, execpt if the current date is after the reaping date of the last season: we have to be at the first season
                                        PEK%MULTI_TREAP(:,1)%TDATE%DAY  ==PEK%TREAP(:)%TDATE%DAY   .AND. & !
                                 ((PEK%MULTI_TSEED(:,3)%TDATE%MONTH==NUNDEF .AND. (PEK%MULTI_TSEED(:,2)%TDATE%MONTH<KMONTH .OR.    & ! check if we are not after the last season (without inversion in this last season)
                                   (PEK%MULTI_TSEED(:,2)%TDATE%MONTH==KMONTH .AND. PEK%MULTI_TSEED(:,2)%TDATE%DAY  <=KDAY))).OR.   &
                                  (PEK%MULTI_TSEED(:,3)%TDATE%MONTH/=NUNDEF .AND. (PEK%MULTI_TSEED(:,3)%TDATE%MONTH<KMONTH .OR.    &
                                   (PEK%MULTI_TSEED(:,3)%TDATE%MONTH==KMONTH .AND. PEK%MULTI_TSEED(:,3)%TDATE%DAY  <=KDAY)))))).OR.&
              ( GINV(:)      .AND. ((PEK%MULTI_TSEED(:,1)%TDATE%MONTH >   KMONTH)                                 .OR. & ! IF inversion, wrong season if it's moreover not during the last season
                               (PEK%MULTI_TSEED(:,1)%TDATE%MONTH ==  KMONTH .AND. PEK%MULTI_TSEED(:,1)%TDATE%DAY > KDAY)) ) )
        ! Check in witch season we are and them update dates in concequence
        WHERE ( GINV(:) ) ! In case of inversion and multiseason, the next season have to be the first one
          PEK%TSEED(:)%TDATE%MONTH = PEK%MULTI_TSEED(:,1)%TDATE%MONTH
          PEK%TSEED(:)%TDATE%DAY   = PEK%MULTI_TSEED(:,1)%TDATE%DAY
          PEK%TREAP(:)%TDATE%MONTH = PEK%MULTI_TREAP(:,1)%TDATE%MONTH
          PEK%TREAP(:)%TDATE%DAY   = PEK%MULTI_TREAP(:,1)%TDATE%DAY
        ELSEWHERE ( PEK%TSEED(:)%TDATE%MONTH == PEK%MULTI_TSEED(:,1)%TDATE%MONTH .AND. &
                    PEK%TSEED(:)%TDATE%DAY   == PEK%MULTI_TSEED(:,1)%TDATE%DAY )  ! In case that we was is the first season
          PEK%TSEED(:)%TDATE%MONTH = PEK%MULTI_TSEED(:,2)%TDATE%MONTH
          PEK%TSEED(:)%TDATE%DAY   = PEK%MULTI_TSEED(:,2)%TDATE%DAY
          PEK%TREAP(:)%TDATE%MONTH = PEK%MULTI_TREAP(:,2)%TDATE%MONTH
          PEK%TREAP(:)%TDATE%DAY   = PEK%MULTI_TREAP(:,2)%TDATE%DAY
        ELSEWHERE ( PEK%TSEED(:)%TDATE%MONTH == PEK%MULTI_TSEED(:,2)%TDATE%MONTH .AND. &
                    PEK%TSEED(:)%TDATE%DAY   == PEK%MULTI_TSEED(:,2)%TDATE%DAY   .AND. & 
                    PEK%MULTI_TSEED(:,3)%TDATE%MONTH /= NUNDEF )  ! In case that we was is the second season and there is a third one
          PEK%TSEED(:)%TDATE%MONTH = PEK%MULTI_TSEED(:,3)%TDATE%MONTH
          PEK%TSEED(:)%TDATE%DAY   = PEK%MULTI_TSEED(:,3)%TDATE%DAY
          PEK%TREAP(:)%TDATE%MONTH = PEK%MULTI_TREAP(:,3)%TDATE%MONTH
          PEK%TREAP(:)%TDATE%DAY   = PEK%MULTI_TREAP(:,3)%TDATE%DAY
        ELSEWHERE  ! In case that we was is the second season and there is not a third one
          PEK%TSEED(:)%TDATE%MONTH = PEK%MULTI_TSEED(:,1)%TDATE%MONTH
          PEK%TSEED(:)%TDATE%DAY   = PEK%MULTI_TSEED(:,1)%TDATE%DAY
          PEK%TREAP(:)%TDATE%MONTH = PEK%MULTI_TREAP(:,1)%TDATE%MONTH
          PEK%TREAP(:)%TDATE%DAY   = PEK%MULTI_TREAP(:,1)%TDATE%DAY
        ENDWHERE
        !
      ENDWHERE
      !
  ENDWHERE
  !
ENDIF
!
IF (LHOOK) CALL DR_HOOK('VEGETATION_EVOL:AGRI_OR_IRRIG',1,ZHOOK_HANDLE)
!
END SUBROUTINE AGRI_OR_IRRIG
!
!-----------------------------------------------------------------
!
END SUBROUTINE VEGETATION_EVOL
