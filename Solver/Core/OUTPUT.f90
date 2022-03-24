!--------------------------------------------------------------------------------------
!
!   Copyright 2014 Ecole Centrale de Nantes, 1 rue de la Noë, 44300 Nantes, France
!
!   Licensed under the Apache License, Version 2.0 (the "License");
!   you may not use this file except in compliance with the License.
!   You may obtain a copy of the License at
!
!       http://www.apache.org/licenses/LICENSE-2.0
!
!   Unless required by applicable law or agreed to in writing, software
!   distributed under the License is distributed on an "AS IS" BASIS,
!   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
!   See the License for the specific language governing permissions and
!   limitations under the License. 
!
!   Contributors list:
!   - A. Babarit 
!
!--------------------------------------------------------------------------------------
MODULE OUTPUT

  USE Constants
  USE MMesh

  IMPLICIT NONE

  !PUBLIC :: WRITE_DATA_ON_MESH

  PUBLIC 
  INTEGER, PARAMETER :: RAW_OUTPUT=0, TECPLOT_OUTPUT=1
  INTEGER :: output_format = TECPLOT_OUTPUT

CONTAINS
!
!     AC: add write sources for QTF
!    SUBROUTINE WRITE_SOURCES(ID,Pbnumber)
!    USE MIDENTIFICATION
!!     USE MMesh
!    USE COM_VAR
!    TYPE(TID) :: ID
!    INTEGER :: Pbnumber,i
!    CHARACTER*5 :: str
!    WRITE(str,'(I5)') Pbnumber
!    OPEN(11,FILE=ID%ID(1:ID%lID)//'/results/sources.'//str//'.dat')
!    
!    do i=1,NFA
!      WRITE(11,'(2(X,E14.7))') REAL(ZIGB(I)), IMAG(ZIGB(I))
!    end do
!    
!    do i=1,NFA
!      WRITE(11,'(2(X,E14.7))') REAL(ZIGS(I)), IMAG(ZIGS(I))
!    end do
!    
!    CLOSE(11)
!    END SUBROUTINE WRITE_SOURCES
    
!    SUBROUTINE WRITE_KOCHIN(ID,Pbnumber,HKochin,Ntheta,Theta)
!    USE MIDENTIFICATION
!    IMPLICIT NONE
!    TYPE(TID) :: ID
!    INTEGER :: Pbnumber
!    COMPLEX,DIMENSION(*) :: HKochin
!    REAL,DIMENSION(*) :: Theta
!    INTEGER :: Ntheta,i
!    CHARACTER*5 :: str
!    WRITE(str,'(I5)') Pbnumber
!    OPEN(11,FILE=ID%ID(1:ID%lID)//'/results/Kochin.'//str//'.dat')
!    DO i=1,Ntheta
!        WRITE(11,*) Theta(i),ABS(HKochin(i)),ATAN2(IMAG(HKochin(i)),REAL(HKochin(i)))
!    END DO
!    CLOSE(11)
!    END SUBROUTINE WRITE_KOCHIN
!
   ! SUBROUTINE WRITE_FS(ID,Pbnumber,PHI,MeshFS)
   ! USE MIDENTIFICATION
   ! USE MMesh
   ! IMPLICIT NONE
   ! INTEGER :: Pbnumber
   ! TYPE(TID) :: ID
   ! TYPE(TMesh) :: MeshFS
   ! INTEGER:: i
   ! COMPLEX,DIMENSION(*) :: PHI
   ! CHARACTER*5 :: str
   ! WRITE(str,'(I5)') Pbnumber
   ! OPEN(11,FILE=ID%ID(1:ID%lID)//'/results/freesurface.'//str//'.dat')
   ! WRITE(11,*) 'VARIABLES="X" "Y" "abs(eta) (m)" "angle(phi) (rad)" "PRE1" "PRE2"'
   ! WRITE(11,'(A,I7,A,I7,A)') 'ZONE N=',MeshFS%Npoints,' ,E = ',MeshFS%Npanels,' ,F=FEPOINT,ET=QUADRILATERAL'
   ! DO i=1,MeshFS%Npoints
!  !     	Adrien Combourieu: remove "-" before real and imaginary part to be consistent with amplitude and phase
   !     WRITE(11,'(6(X,E14.7))') MeshFS%X(1,i),MeshFS%X(2,i),ABS(PHI(i)),ATAN2(IMAG(PHI(i)),REAL(PHI(i))),REAL(PHI(i)),IMAG(PHI(i)) 
   ! END DO
   ! DO i=1,MeshFS%Npanels
   !     WRITE(11,*) MeshFS%P(1,i),MeshFS%P(2,i),MeshFS%P(3,i),MeshFS%P(4,i)
   ! END DO
   ! CLOSE(11)
   ! END SUBROUTINE WRITE_FS 
!
!    SUBROUTINE WRITE_POTENTIAL(ID,Pbnumber,PRESSURE)
!    USE MIDENTIFICATION
!    USE MMesh
!    USE COM_VAR
!    IMPLICIT NONE
!    INTEGER :: Pbnumber
!    TYPE(TID) :: ID
!    INTEGER:: i
!    COMPLEX,DIMENSION(*) :: PRESSURE
!    CHARACTER*5 :: str
!    WRITE(str,'(I5)') Pbnumber
!    OPEN(11,FILE=ID%ID(1:ID%lID)//'/results/pressure.'//str//'.dat')
!    WRITE(11,*) 'VARIABLES="X" "Y" "Z" "abs(p) (Pa)" "angle(p) (rad)"'
!    WRITE(11,'(A,I7,A,I7,A)') 'ZONE N=',4*NFA*2**NSYMY,',E = ',NFA*2**NSYMY,',F=FEPOINT,ET=QUADRILATERAL'
!    DO i=1,NFA
!        WRITE(11,*) X(M1(i)),Y(M1(i)),Z(M1(i)),ABS(PRESSURE(i)),ATAN2(IMAG(PRESSURE(i)),REAL(PRESSURE(i))) 
!        WRITE(11,*) X(M2(i)),Y(M2(i)),Z(M2(i)),ABS(PRESSURE(i)),ATAN2(IMAG(PRESSURE(i)),REAL(PRESSURE(i))) 
!        WRITE(11,*) X(M3(i)),Y(M3(i)),Z(M3(i)),ABS(PRESSURE(i)),ATAN2(IMAG(PRESSURE(i)),REAL(PRESSURE(i))) 
!        WRITE(11,*) X(M4(i)),Y(M4(i)),Z(M4(i)),ABS(PRESSURE(i)),ATAN2(IMAG(PRESSURE(i)),REAL(PRESSURE(i))) 
!    END DO
!
!  IF (NSYMY.EQ.1) THEN
!        DO i=1,NFA
!            WRITE(11,*) X(M1(i)),-Y(M1(i)),Z(M1(i)),ABS(PRESSURE(i+NFA)),ATAN2(IMAG(PRESSURE(i+NFA)),REAL(PRESSURE(i+NFA)))
!            WRITE(11,*) X(M2(i)),-Y(M2(i)),Z(M2(i)),ABS(PRESSURE(i+NFA)),ATAN2(IMAG(PRESSURE(i+NFA)),REAL(PRESSURE(i+NFA)))
!            WRITE(11,*) X(M3(i)),-Y(M3(i)),Z(M3(i)),ABS(PRESSURE(i+NFA)),ATAN2(IMAG(PRESSURE(i+NFA)),REAL(PRESSURE(i+NFA)))
!            WRITE(11,*) X(M4(i)),-Y(M4(i)),Z(M4(i)),ABS(PRESSURE(i+NFA)),ATAN2(IMAG(PRESSURE(i+NFA)),REAL(PRESSURE(i+NFA)))
!        END DO
!    END IF  
!    DO i=1,NFA
!        WRITE(11,*) 1+(i-1)*4,2+(i-1)*4,3+(i-1)*4,4+(i-1)*4
!    END DO
!    IF (NSYMY.EQ.1) THEN
!        DO i=1,NFA
!            WRITE(11,*) 1+(i-1)*4+4*NFA,2+(i-1)*4+4*NFA,3+(i-1)*4+4*NFA,4+(i-1)*4+4*NFA
!        END DO
!    END IF  
!!    DO i=1,NP
!!        WRITE(11,*) X(i),Y(i),Z(i),ABS(PRESSURE(i)),ATAN2(IMAG(PRESSURE(i)),REAL(PRESSURE(i))) 
!!    END DO
!!    IF (NSYMY.EQ.1) THEN
!!        DO i=1,NP
!!            WRITE(11,*) X(i),-Y(i),Z(i),ABS(PRESSURE(i+NFA)),ATAN2(IMAG(PRESSURE(i+NFA)),REAL(PRESSURE(i+NFA)))
!!        END DO
!!    END IF
!!    DO i=1,NFA
!!        WRITE(11,*) M1(i),M2(i),M3(i),M4(i)
!!    END DO
!!    IF (NSYMY.EQ.1) THEN
!!        DO i=1,NFA
!!            WRITE(11,*) M1(i)+NP,M2(i)+NP,M3(i)+NP,M4(i)+NP
!!        END DO
!!    END IF
!!    WRITE(11,'(A,I7)') 'ZONE t="Pressure",F=POINT,I=',NFA*2**NSYMY
!!    DO i=1,NFA
!!        WRITE(11,*) XG(i),YG(i),ZG(i),ABS(PRESSURE(i)),ATAN2(IMAG(PRESSURE(i)),REAL(PRESSURE(i))) 
!!    END DO
!!    IF (NSYMY.EQ.1) THEN
!!        DO i=1,NFA
!!            WRITE(11,*)  XG(i),-YG(i),ZG(i),ABS(PRESSURE(i+NFA)),ATAN2(IMAG(PRESSURE(i+NFA)),REAL(PRESSURE(i+NFA))) 
!!        END DO
!!    END IF
!    CLOSE(11) 
!    END SUBROUTINE WRITE_POTENTIAL

  SUBROUTINE WRITE_DATA_ON_MESH(Mesh, cdata, filename)
    ! Save in a file a field of complex values on a mesh.

    TYPE(TMesh),                                   INTENT(IN) :: Mesh
    COMPLEX, DIMENSION(Mesh%NPanels*2**Mesh%ISym), INTENT(IN) :: cdata
    CHARACTER(LEN=*),                              INTENT(IN) :: filename ! Output file
    
    ! Local variables
    INTEGER            :: u, i, j
    REAL, DIMENSION(3) :: x

    OPEN(NEWUNIT=u, FILE=filename, ACTION='WRITE')

    IF (output_format == TECPLOT_OUTPUT) THEN
      WRITE(u,*) 'VARIABLES="X" "Y" "Z" "abs(p) (Pa)" "angle(p) (rad)"'
      WRITE(u,'(A,I7,A,I7,A)') 'ZONE N=', Mesh%Npoints*2**Mesh%ISym,',E = ', Mesh%Npanels*2**Mesh%ISym,', F=FEPOINT, ET=QUADRILATERAL'
    END IF 

    DO i = 1, Mesh%NPanels
      DO j = 1, 4
        WRITE(u, *) Mesh%X(:, Mesh%P(j, i)), ABS(cdata(i)), ATAN2(AIMAG(cdata(i)),REAL(cdata(i)))
      END DO
    END DO

    IF (Mesh%ISym == Y_SYMMETRY) THEN
      DO i = Mesh%NPanels+1, 2*Mesh%NPanels
        DO j = 1, 4
          x(:) = Mesh%X(:, Mesh%P(j, i-Mesh%NPanels))
          x(2) = -x(2)
          WRITE(u, *) x(:), ABS(cdata(i)), ATAN2(AIMAG(cdata(i)),REAL(cdata(i)))
        END DO
      END DO
    END IF

    IF (output_format == TECPLOT_OUTPUT) THEN
      DO i=1,Mesh%Npanels
        WRITE(u, *) 1+(i-1)*4, 2+(i-1)*4, 3+(i-1)*4, 4+(i-1)*4
      END DO

      IF (Mesh%ISym == Y_SYMMETRY) THEN
        DO i=1,Mesh%Npanels
          WRITE(u, *) 1+(i-1)*4+4*Mesh%NPanels, 2+(i-1)*4+4*Mesh%NPanels, 3+(i-1)*4+4*Mesh%NPanels, 4+(i-1)*4+4*Mesh%NPanels
        END DO
      END IF
    END IF

    CLOSE(u)

  END SUBROUTINE WRITE_DATA_ON_MESH

END MODULE OUTPUT
