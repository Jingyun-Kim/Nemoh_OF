!--------------------------------------------------------------------------------------
!
!   NEMOH - March 2022
!   Contributors list:
!   - R. Kurnia
!--------------------------------------------------------------------------------------
MODULE CONSTANTS

      IMPLICIT NONE

      INTEGER, PARAMETER        :: ID_DP= 0     ! 0 for single precision 1 double
                                                ! set 1 if it is compile with option -r8
      REAL, PARAMETER           :: PI=4.*ATAN(1.), DPI=2.*PI, DPI2=2.*PI**2  
      COMPLEX, PARAMETER        :: II=CMPLX(0.,1.)
      REAL, PARAMETER           :: INFINITE_DEPTH=0.0
      INTEGER, PARAMETER        :: NO_Y_SYMMETRY=0, Y_SYMMETRY=1
      REAL, PARAMETER           :: EPS=0.0001, ZERO=0.
      COMPLEX, PARAMETER        :: CZERO=CMPLX(0.0,0.0)
      INTEGER, PARAMETER        :: DIFFRACTION_PROBLEM=1, RADIATION_PROBLEM=-1
END MODULE CONSTANTS
