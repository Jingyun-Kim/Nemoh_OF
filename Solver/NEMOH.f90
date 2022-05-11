!--------------------------------------------------------------------------------------
!
!   NEMOH V1.0 - BVP solver - January 2014
!
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
!   - G. Delhommeau
!   - P. Guével
!   - J.C. Daubisse
!   - J. Singh
!   - A. Babarit
!
!--------------------------------------------------------------------------------------

PROGRAM Main

  USE Constants
  USE MMesh,                ONLY: TMesh,           ReadTMesh
  USE MEnvironment,         ONLY: TEnvironment,    ReadTEnvironment
  USE MBodyConditions,      ONLY: TBodyConditions, ReadTBodyConditions
  USE M_Solver,             ONLY: TSolver,         ReadTSolver, ID_GMRES
  USE MLogFile              ! 
  ! Preprocessing and initialization
  USE MFace,                ONLY: TVFace, Prepare_FaceMesh
  USE M_INITIALIZE_GREEN,   ONLY: TGREEN, INITIALIZE_GREEN
  USE Elementary_functions, ONLY: X0

  ! Resolution
  USE SOLVE_BEM_DIRECT,     ONLY: SOLVE_POTENTIAL_DIRECT
  ! Post processing and output
  USE OUTPUT,               ONLY: WRITE_DATA_ON_MESH,WRITE_SOURCES
  USE FORCES,               ONLY: COMPUTE_AND_WRITE_FORCES
  USE KOCHIN,               ONLY: COMPUTE_AND_WRITE_KOCHIN
  USE FREESURFACE,          ONLY: COMPUTE_AND_WRITE_FREE_SURFACE_ELEVATION

  IMPLICIT NONE

  CHARACTER(LEN=1000)   :: wd             ! Working directory path (max length: 1000 characters, increase if necessary)
  TYPE(TMesh)           :: Mesh           ! Mesh of the floating body
  TYPE(TBodyConditions) :: BodyConditions ! Physical conditions on the floating body
  TYPE(TEnvironment)    :: Env            ! Physical conditions of the environment
  TYPE(TSolver)         :: SolverOpt      ! Solver Option, specified by user in input_solver.txt 
  
  INTEGER                            :: i_problem          ! Index of the current problem
  REAL                               :: omega, wavenumber  ! Wave frequency and wavenumber
  COMPLEX, DIMENSION(:), ALLOCATABLE     :: ZIGB, ZIGS     ! Computed source distribution
  COMPLEX, DIMENSION(:,:,:), ALLOCATABLE :: V, Vinv,S      ! Influece coefficients
  COMPLEX, DIMENSION(:), ALLOCATABLE :: Potential          ! Computed potential

  TYPE(TVFACE)                       :: VFace              ! Face Mesh structure variable                   
  TYPE(TGREEN)                       :: IGreen             ! Initial Green variables
  
  REAL                               :: tcpu_start
  CHARACTER(LEN=1000)                :: LogTextToBeWritten

  ! Initialization ---------------------------------------------------------------------

  WRITE(*,*) ' '
  WRITE(*,'(A,$)') '  -> Initialisation '

  ! Get working directory from command line argument
  IF (COMMAND_ARGUMENT_COUNT() >= 1) THEN
    CALL GET_COMMAND_ARGUMENT(1, wd)
  ELSE
    wd = "."
  END IF

  CALL ReadTMesh(Mesh, TRIM(wd)//'/mesh/')
  ALLOCATE(ZIGB(Mesh%NPanels), ZIGS(Mesh%NPanels))
  ALLOCATE(Potential(Mesh%NPanels*2**Mesh%Isym))

  CALL ReadTBodyConditions            &
  ( BodyConditions,                   &
    Mesh%Npanels*2**Mesh%Isym,        &
    TRIM(wd)//'/Normalvelocities.dat' &
    )

  CALL ReadTEnvironment(Env, file=TRIM(wd)//'/Nemoh.cal')

  CALL ReadTSolver(SolverOpt,TRIM(wd))

  CALL Prepare_FaceMesh(Mesh,SolverOpt%NP_GQ,VFace)

  CALL INITIALIZE_GREEN(VFace,Mesh,Env%depth,IGreen)
  ALLOCATE(S(Mesh%NPanels,Mesh%NPanels,2**Mesh%Isym))
  ALLOCATE(V(Mesh%NPanels,Mesh%NPanels,2**Mesh%Isym))
  ALLOCATE(Vinv(Mesh%NPanels,Mesh%NPanels,2**Mesh%Isym))

  WRITE(*, *) ' '
  WRITE(LogTextToBeWritten,*) 'NP Gauss Quadrature Integ.: ', SolverOpt%NP_GQ
  CALL WRITE_LOGFILE(trim(wd)//'/logfile.txt',TRIM(LogTextToBeWritten),IdStartLog,IdprintTerm)
  WRITE(*, *) '. Done !'
  WRITE(*, *) ' '

  ! Solve BVPs and calculate forces ----------------------------------------------------
  WRITE(*, *) ' -> Solve BVPs and calculate forces '
  WRITE(LogTextToBeWritten,*) 'Linear Solver: ', SolverOpt%SNAME
  CALL WRITE_LOGFILE(trim(wd)//'/logfile.txt',TRIM(LogTextToBeWritten),IdAppend,IdprintTerm)
  CALL START_RECORD_TIME(tcpu_start,trim(wd)//'/logfile.txt',IdAppend)
  WRITE(*, *) ' '

  DO i_problem = 1, BodyConditions%Nproblems
    WRITE(*,'(A,I5,A,I5,A,A,$)') ' Problem ',i_problem,' / ',BodyConditions%Nproblems,' ',CHAR(13)

    omega = BodyConditions%omega(i_problem) ! Wave frequency
    ! Compute wave number
    IF ((Env%depth == INFINITE_DEPTH) .OR. (omega**2*Env%depth/Env%g >= 20)) THEN
      wavenumber = omega**2/Env%g
    ELSE
      wavenumber = X0(omega**2*Env%depth/Env%g)/Env%depth
      ! X0(y) returns the solution of y = x * tanh(x)
    END IF
    !===============
    ! BEM Resolution
    !===============
      CALL SOLVE_POTENTIAL_DIRECT                                              &
      !==========================
      ( VFace, Mesh, Env, omega, wavenumber,IGreen,                            &
        BodyConditions%NormalVelocity(1:Mesh%Npanels*2**Mesh%Isym, i_problem), &
        S,V,Vinv,ZIGB, ZIGS,                                                   &
        Potential(:),SolverOpt,trim(wd))

    !===========================
    ! Post processing and output
    !===========================

    CALL COMPUTE_AND_WRITE_FORCES            &
    !============================
    ( TRIM(wd)//'/mesh/Integration.dat',     &
      Mesh, Env%rho, omega, Potential,       &
      Bodyconditions%Switch_type(i_problem), &
      TRIM(wd)//'/results/Forces.dat'        &
      )
    IF (BodyConditions%Switch_Potential(i_problem) == 1) THEN
      ! Write pressure field on the floating body in file
      CALL WRITE_DATA_ON_MESH                                      &
      !=======================
      ( Mesh,                                                      &
        Env%rho*II*omega*Potential(:),                             &
        TRIM(wd)//'/results/pressure.'//string(i_problem)//'.dat'  &
        )
    END IF

    IF (BodyConditions%Switch_Kochin(i_problem) == 1) THEN
      CALL COMPUTE_AND_WRITE_KOCHIN                             &
      !============================
      ( TRIM(wd)//'/mesh/Kochin.dat',                           &
        Mesh, Env, wavenumber, ZIGB, ZIGS,                      &
        TRIM(wd)//'/results/Kochin.'//string(i_problem)//'.dat' &
        )
    END IF
    IF (BodyConditions%Switch_FreeSurface(i_problem) == 1) THEN
      CALL COMPUTE_AND_WRITE_FREE_SURFACE_ELEVATION                  &
      !============================================
      ( TRIM(wd)//'/mesh/Freesurface.dat', IGreen,VFace,             &
        Mesh, Env, omega, wavenumber, ZIGB, ZIGS,                    &
        TRIM(wd)//'/results/freesurface.'//string(i_problem)//'.dat' &
        )
    END IF
    
    IF (BodyConditions%Switch_SourceDistr(i_problem) == 1) THEN
      ! Write pressure field on the floating body in file
      CALL WRITE_SOURCES(ZIGB,ZIGS,Mesh%Npanels,                    &
       TRIM(wd)//'/results/sources/sources.'//string(i_problem)//'.dat')
    END IF

  END DO
  CALL END_RECORD_TIME(tcpu_start,trim(wd)//'/logfile.txt')
  WRITE(*,*) '. Done !'
  ! Finalize ---------------------------------------------------------------------------

  DEALLOCATE(ZIGB, ZIGS, Potential,S,V,Vinv)

CONTAINS

  FUNCTION string (i) result (s)
    ! For example 5 -> "00005"
    INTEGER :: i
    CHARACTER(LEN=5) :: s
    WRITE(s, '(I0.5)') i
  END FUNCTION
      
END PROGRAM Main
