!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!-----------------------------------------------------------------
!     #######################
!
      SUBROUTINE WRITE_DISCHARGE_FILE_SUB(HPROGRAM,HFILE,HFORM,KYEAR,KMONTH,KDAY,KH,KM)
!     #######################
!
!!****  *WRITE_DISCHARGE_FILE*  
!!
!!    PURPOSE
!!    -------
!     This routine aims at reading topographic files
!!**  METHOD
!!    ------
!
!!    EXTERNAL
!!    --------
!!
!!    none
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------ 
!!
!!    
!!    
!!
!!      
!!    REFERENCE
!!    ---------
!!
!!    
!!      
!!    AUTHOR
!!    ------
!!
!!      B. Vincendon    * Meteo-France *
!!
!!    MODIFICATIONS
!!    -------------
!!
!!      Original   11/2006
!-------------------------------------------------------------------------------
!
!*       0.     DECLARATIONS
!               ------------
!
USE MODD_TOPODYN, ONLY  : CCAT, NNCAT, NNB_TOPD_STEP, NNCAT_SUB, XQTOT_SUB
USE MODD_COUPLING_TOPD, ONLY : CSUBCAT
!
USE MODI_GET_LUOUT
USE MODI_OPEN_FILE
USE MODI_CLOSE_FILE
!
USE YOMHOOK   ,ONLY : LHOOK, DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
CHARACTER(LEN=*),      INTENT(IN)  :: HPROGRAM   !
CHARACTER(LEN=*),      INTENT(IN)  :: HFILE      ! File to be read
CHARACTER(LEN=*),      INTENT(IN)  :: HFORM      ! Format of the file to be read
INTEGER, DIMENSION(:), INTENT(IN)  :: KYEAR      ! Year of the beginning of the simulation.
INTEGER, DIMENSION(:), INTENT(IN)  :: KMONTH     ! Month of the beginning of the simulation.
INTEGER, DIMENSION(:), INTENT(IN)  :: KDAY       ! Day of the beginning of the simulation.
INTEGER, DIMENSION(:), INTENT(IN)  :: KH         ! Hour of the beginning of the simulation.
INTEGER, DIMENSION(:), INTENT(IN)  :: KM         ! Minutes of the beginning of the simulation.
!
!
!*      0.2    declarations of local variables
!
INTEGER                   :: JJ,JCAT ! loop control 
INTEGER                   :: IUNIT       ! Unit of the files
INTEGER                   :: ILUOUT      ! Unit of the files
!
CHARACTER(LEN=28) :: YFILE
CHARACTER(LEN=40) :: YFORM          ! Writing format
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('WRITE_DISCHARGE_FILE_SUB',0,ZHOOK_HANDLE)
!
!*       0.3    preparing file openning
!               ----------------------
!
CALL GET_LUOUT(HPROGRAM,ILUOUT)
YFORM='(I4,A1,I2,A1,I2,A1,I2,A1,I2,A1,F7.2)'
!
DO JCAT=1,NNCAT
  !
  YFILE = TRIM(CCAT(JCAT))//'_'//TRIM(HFILE)
  !
  CALL OPEN_FILE(HPROGRAM,IUNIT,YFILE,HFORM,HACTION='WRITE')
  !
  WRITE(IUNIT,*) 'YEAR;MO;DA;HO;MI;',CCAT(JCAT),CSUBCAT(JCAT,1:NNCAT_SUB(JCAT))
  DO JJ=1,NNB_TOPD_STEP
    WRITE(IUNIT,*) KYEAR(JJ),';',KMONTH(JJ),';',KDAY(JJ),';',&
                       KH(JJ)   ,';',KM(JJ)    ,';',&
                       XQTOT_SUB(JCAT,1:NNCAT_SUB(JCAT),JJ)
  ENDDO
  !
  CALL CLOSE_FILE(HPROGRAM,IUNIT)
  !
ENDDO
!
IF (LHOOK) CALL DR_HOOK('WRITE_DISCHARGE_FILE_SUB',1,ZHOOK_HANDLE)
!
END SUBROUTINE WRITE_DISCHARGE_FILE_SUB

