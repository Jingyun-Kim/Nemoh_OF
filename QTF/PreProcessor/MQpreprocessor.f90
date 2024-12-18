!--------------------------------------------------------------------------------------
!
!    Copyright (C) 2022 - LHEEA Lab., Ecole Centrale de Nantes, UMR CNRS 6598
!
!    This program is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    This program is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with this program.  If not, see <http://www.gnu.org/licenses/>.
!--------------------------------------------------------------------------------------
!   Contributors list:
!   - Gerard Delhommeau (13/11/2014, d'apres LA VERSION 1.1 d'AVRIL 1991
!     PROGRAMME StOK LABORATOIRE D'HYDRODYNAMIQUE NAVALE DE L'E.N.S.M. DE NANTES & SIREHNA )
!   - Adrien Combourieu, INNOSEA (adrien.combourieu@innosea.fr) 2014
!   - Ruddy Kurnia (LHEEA,ECN) 2022
!--------------------------------------------------------------------------------------
!
!   NEMOH2 - second order (QTF) - May 2022
!
!--------------------------------------------------------------------------------------
Module MQpreprocessor
USE CONSTANTS
USE MReadInputFiles,            ONLY: TSource,TMeshFS
USE MMesh,                      ONLY: TMesh
USE MFace,                      ONLY: TVFace, TWLine
USE MEnvironment,               ONLY: TEnvironment,Fun_inverseDispersion, &
                                      COMPUTE_INC_POTENTIAL_VELOCITY
USE Elementary_functions,       ONLY: DOT_PRODUCT_COMPLEX,&
                                      X0 !invers of disp. relation

! Green functions
USE M_INITIALIZE_GREEN,         ONLY: TGREEN,LISC
USE MInfluenceMatrix,           ONLY: construct_influence_matrix, &
                                      construct_influence_matrix_FS
!
USE MFileDirectoryList
USE MLogFile
!
IMPLICIT NONE

INTEGER,parameter :: ID_DEBUG=0 ! for debugging potential will be saved

CONTAINS

  SUBROUTINE COMPUTE_POTENTIALS_AND_VELOCITIES_FS(wd,Iw,omega,Vbeta,Nbeta,Nradiation,&
                                              Env,Mesh,VFace,IGreen,         &
                                              SourceDistr,MotionIw,MeshFS)
        !Potential and velocities on free surface panels
        !direct computation for each flow point-Ipflow to avoid large memory used
        !INPUT/OUTPUT
        CHARACTER(LEN=*),               INTENT(IN)::wd
        INTEGER,                        INTENT(IN)::Iw,Nbeta,Nradiation
        REAL,                           INTENT(IN)::omega       !rad freq w(Iw)
        REAL,DIMENSION(Nbeta),          INTENT(IN)::Vbeta       !Angle vector
        TYPE(TEnvironment),             INTENT(IN)::Env
        TYPE(TMesh),                    INTENT(IN)::Mesh        !body mesh
        TYPE(TVFace),                   INTENT(IN)::VFace       !body face
        TYPE(TGREEN),                   INTENT(INOUT)::IGreen
        TYPE(TSource),                  INTENT(IN)::SourceDistr
        COMPLEX,DIMENSION(Nradiation,Nbeta),INTENT(IN):: MotionIw
        TYPE(TMeshFS),                  INTENT(IN)::MeshFS      !Freesurface mesh
        !LOCAL
        INTEGER                                ::NPFLOW,uFile,ILINE,Ipflow
        REAL                                   :: wavenumber
        COMPLEX, DIMENSION(:,:)  , ALLOCATABLE :: S       ! Inf. coef. Integ. of Green func.
        COMPLEX, DIMENSION(:,:,:), ALLOCATABLE :: GradS   ! Inf. coef. Integ. of Gradient Green func.
        INTEGER                                :: Ibeta,Ipanel,Irad
        COMPLEX,DIMENSION(Mesh%Npanels,Nbeta)  :: ZPGB,ZPGS ! Perturbation sources, B for Body
                                                            ! S for symmetric part
        COMPLEX,ALLOCATABLE,DIMENSION(:,:)       ::IncPotential  ! Incoming potential
        COMPLEX,ALLOCATABLE,DIMENSION(:,:,:)     ::IncVelocity   ! Incoming velocity
        COMPLEX,ALLOCATABLE,DIMENSION(:,:)       ::Potential    ! Total potential, Perturb+Incoming pot
        COMPLEX,ALLOCATABLE,DIMENSION(:,:,:)     ::Velocity     ! Total velocity
        COMPLEX,ALLOCATABLE,DIMENSION(:,:)       ::RadPotential ! Rad pot mod I
        COMPLEX,ALLOCATABLE,DIMENSION(:,:,:)     ::RadVelocity  ! Rad velocity mod I
        CHARACTER(LEN=1000)                      :: LogTextToBeWritten
        COMPLEX                                  :: PotIter,VxIter,VyIter,VzIter

        REAL                                     :: XM_I(3)

        NPFLOW=MeshFS%Mesh%NPanels+MeshFS%BdyLine%NWlineseg!computed points on FS
        ALLOCATE(S(Mesh%NPanels,2**Mesh%Isym))         ! same size as body panels
        ALLOCATE(GradS(Mesh%NPanels,3,2**Mesh%Isym))
        ALLOCATE(IncPotential(NPFLOW*2**Mesh%Isym,Nbeta)) ! same size as free surface panels
        ALLOCATE(IncVelocity(NPFLOW*2**Mesh%Isym,3,Nbeta))
        ALLOCATE(Potential(NPFLOW*2**Mesh%Isym,Nbeta)) ! same size as free surface panels
        ALLOCATE(Velocity(NPFLOW*2**Mesh%Isym,3,Nbeta))
        ALLOCATE(RadPotential(NPFLOW*2**Mesh%Isym,Nradiation))
        ALLOCATE(RadVelocity(NPFLOW*2**Mesh%Isym,3,Nradiation))

        wavenumber=Fun_inverseDispersion(omega,Env%depth,Env%g)



        DO Ibeta=1,Nbeta
           !CONSTRUCT PERTURBATION(Diff+Rad) SINGULAR DISTRIBUTION FOR EACH WAVE DIRECTION
           !assign the diffraction singular source distribution for all panels
           ZPGB(1:Mesh%Npanels,Ibeta)=SourceDistr%ZIGB(1:Mesh%Npanels,Nradiation+Ibeta)
           ZPGS(1:Mesh%Npanels,Ibeta)=SourceDistr%ZIGS(1:Mesh%Npanels,Nradiation+Ibeta)
              !sum the diffraction + radiation singular source distribution for all panels
              DO Irad=1,Nradiation
              ZPGB(1:Mesh%Npanels,Ibeta)=ZPGB(1:Mesh%Npanels,Ibeta)                             &
                                 -II*omega*SourceDistr%ZIGB(1:Mesh%Npanels,Irad)                &
                                  *MotionIw(Irad,Ibeta)
              ZPGS(1:Mesh%Npanels,Ibeta)=ZPGS(1:Mesh%Npanels,Ibeta)                             &
                                 -II*omega*SourceDistr%ZIGS(1:Mesh%Npanels,Irad)                &
                                  *MotionIw(Irad,Ibeta)
              ENDDO
              !--------------------------------------------------------------------------------
              ! Compute incoming Potential and velocity
              CALL COMPUTE_INC_POTENTIAL_VELOCITY(wavenumber,omega,Vbeta(Ibeta),                &
                                  MeshFS%VFace%XM,MeshFS%Mesh%Npanels,MeshFS%BdyLine%XM,        &
                                  MeshFS%BdyLine%NWlineseg,Env,MeshFS%Mesh%ISym,                &
                                  IncPotential(:,Ibeta),IncVelocity(:,:,Ibeta))
              !--------------------------------------------------------------------------------
        ENDDO

           IF ((.NOT.Env%depth == INFINITE_DEPTH)) THEN
           !for the green function
                IF (omega**2*Env%depth/Env%g.LT.20) THEN
                   CALL LISC(omega**2*Env%depth/Env%g, wavenumber*Env%depth,IGreen)
                ENDIF
           END IF

        DO Ipflow=1,NPFLOW
            ! compute the inf coef S and GradS for each row-Ipflow
            IF (Ipflow<=MeshFS%Mesh%NPanels) THEN
             XM_I=MeshFS%VFace%XM(Ipflow,:)
            ELSE
             XM_I=MeshFS%BdyLine%XM(Ipflow-MeshFS%Mesh%NPanels,:)
            ENDIF

            CALL CONSTRUCT_INFLUENCE_MATRIX_FS(omega,wavenumber,Env,IGreen,                                 &
                            Mesh,VFace,XM_I,S,GradS)
            !---------------------------------------------------------
            !Compute the physical variables for each flow point-Ipflow
            DO Ibeta=1,Nbeta
                  ! Compute total potential and velocity
              IF (Mesh%ISym==NO_Y_SYMMETRY) THEN
                Potential(Ipflow,Ibeta)=IncPotential(Ipflow,Ibeta)+                                         &
                     DOT_PRODUCT_COMPLEX(S(1:Mesh%Npanels,1),ZPGB(1:Mesh%Npanels,Ibeta),Mesh%Npanels)
                !Vx
                Velocity(Ipflow,1,Ibeta)=IncVelocity(Ipflow,1,Ibeta)+                                       &
                     DOT_PRODUCT_COMPLEX(gradS(1:Mesh%Npanels,1,1),ZPGB(1:Mesh%Npanels,Ibeta),Mesh%Npanels)
                !Vy
                Velocity(Ipflow,2,Ibeta)=IncVelocity(Ipflow,2,Ibeta)+                                       &
                     DOT_PRODUCT_COMPLEX(gradS(1:Mesh%Npanels,2,1),ZPGB(1:Mesh%Npanels,Ibeta),Mesh%Npanels)
                !Vz
                Velocity(Ipflow,3,Ibeta)=IncVelocity(Ipflow,3,Ibeta)+                                       &
                     DOT_PRODUCT_COMPLEX(gradS(1:Mesh%Npanels,3,1),ZPGB(1:Mesh%Npanels,Ibeta),Mesh%Npanels)

              ELSE IF (Mesh%ISym == Y_SYMMETRY) THEN
                Potential(Ipflow,Ibeta)=IncPotential(Ipflow,Ibeta)                                          &
                  +DOT_PRODUCT_COMPLEX(S(:, 1) + S(:, 2), ZPGB(:,Ibeta),Mesh%Npanels)/2                     &
                  +DOT_PRODUCT_COMPLEX(S(:, 1) - S(:, 2), ZPGS(:,Ibeta),Mesh%Npanels)/2
                Potential(NPFLOW+Ipflow,Ibeta) =IncPotential(NPFLOW+IPflow,Ibeta)+                          &
                  +DOT_PRODUCT_COMPLEX(S(:, 1) - S(:, 2), ZPGB(:,Ibeta),Mesh%Npanels)/2                     &
                  +DOT_PRODUCT_COMPLEX(S(:, 1) + S(:, 2), ZPGS(:,Ibeta),Mesh%Npanels)/2
                !Vx
                Velocity(Ipflow,1,Ibeta)=IncVelocity(Ipflow,1,Ibeta)+                                       &
                   +DOT_PRODUCT_COMPLEX(gradS(:, 1, 1) + gradS(:, 1, 2), ZPGB(:,Ibeta),Mesh%Npanels)/2      &
                   +DOT_PRODUCT_COMPLEX(gradS(:, 1, 1) - gradS(:, 1, 2), ZPGS(:,Ibeta),Mesh%Npanels)/2
                Velocity(NPFLOW+Ipflow,1,Ibeta) =IncVelocity(NPFLOW+Ipflow,1,Ibeta)                         &
                   +DOT_PRODUCT_COMPLEX(gradS(:, 1, 1) - gradS(:, 1, 2), ZPGB(:,Ibeta),Mesh%Npanels)/2      &
                   +DOT_PRODUCT_COMPLEX(gradS(:, 1, 1) + gradS(:, 1, 2), ZPGS(:,Ibeta),Mesh%Npanels)/2
                !Vy
                Velocity(Ipflow,2,Ibeta)=IncVelocity(Ipflow,2,Ibeta)+                                       &
                   +DOT_PRODUCT_COMPLEX(gradS(:, 2, 1) + gradS(:, 2, 2), ZPGB(:,Ibeta),Mesh%Npanels)/2      &
                   +DOT_PRODUCT_COMPLEX(gradS(:, 2, 1) - gradS(:, 2, 2), ZPGS(:,Ibeta),Mesh%Npanels)/2
                Velocity(NPFLOW+Ipflow,2,Ibeta) =IncVelocity(NPFLOW+Ipflow,2,Ibeta)                         &
                   -DOT_PRODUCT_COMPLEX(gradS(:, 2, 1) - gradS(:, 2, 2), ZPGB(:,Ibeta),Mesh%Npanels)/2      &
                   -DOT_PRODUCT_COMPLEX(gradS(:, 2, 1) + gradS(:, 2, 2), ZPGS(:,Ibeta),Mesh%Npanels)/2
                !Vz
                 Velocity(Ipflow,3,Ibeta)=IncVelocity(Ipflow,3,Ibeta)+                                      &
                   +DOT_PRODUCT_COMPLEX(gradS(:, 3, 1) + gradS(:, 3, 2), ZPGB(:,Ibeta),Mesh%Npanels)/2      &
                   +DOT_PRODUCT_COMPLEX(gradS(:, 3, 1) - gradS(:, 3, 2), ZPGS(:,Ibeta),Mesh%Npanels)/2
                 Velocity(NPFLOW+Ipflow,3,Ibeta) =IncVelocity(NPFLOW+Ipflow,3,Ibeta)                        &
                   +DOT_PRODUCT_COMPLEX(gradS(:, 3, 1) - gradS(:, 3, 2), ZPGB(:,Ibeta),Mesh%Npanels)/2      &
                   +DOT_PRODUCT_COMPLEX(gradS(:, 3, 1) + gradS(:, 3, 2), ZPGS(:,Ibeta),Mesh%Npanels)/2
              ENDIF


            ENDDO

            ! Compute radiation potential and velocity
            DO IRad=1,Nradiation
                   IF (Mesh%ISym==NO_Y_SYMMETRY) THEN
                    RadPotential(Ipflow,IRad)=                                                                       &
                      DOT_PRODUCT_COMPLEX(S(1:Mesh%Npanels,1),SourceDistr%ZIGB(:,Irad),Mesh%Npanels)
                    !Vx
                    RadVelocity(Ipflow,1,IRad)=                                                                      &
                      DOT_PRODUCT_COMPLEX(gradS(1:Mesh%Npanels,1,1),SourceDistr%ZIGB(:,Irad),Mesh%Npanels)
                    !Vy
                    RadVelocity(Ipflow,2,IRad)=                                                                      &
                      DOT_PRODUCT_COMPLEX(gradS(1:Mesh%Npanels,2,1),SourceDistr%ZIGB(:,Irad),Mesh%Npanels)
                    !Vz
                    RadVelocity(Ipflow,3,IRad)=                                                                      &
                      DOT_PRODUCT_COMPLEX(gradS(1:Mesh%Npanels,3,1),SourceDistr%ZIGB(:,Irad),Mesh%Npanels)
                  ELSE IF (Mesh%ISym == Y_SYMMETRY) THEN
                     RadPotential(Ipflow,IRad)=                                                                      &
                        DOT_PRODUCT_COMPLEX(S(:, 1) + S(:, 2),SourceDistr%ZIGB(:,Irad),Mesh%Npanels)/2               &
                       +DOT_PRODUCT_COMPLEX(S(:, 1) - S(:, 2),SourceDistr%ZIGS(:,Irad),Mesh%Npanels)/2
                     RadPotential(NPFLOW+Ipflow,IRad) =                                                              &
                        DOT_PRODUCT_COMPLEX(S(:, 1) - S(:, 2),SourceDistr%ZIGB(:,Irad),Mesh%Npanels)/2               &
                       +DOT_PRODUCT_COMPLEX(S(:, 1) + S(:, 2),SourceDistr%ZIGS(:,Irad),Mesh%Npanels)/2
                    !Vx
                    RadVelocity(Ipflow,1,IRad)=                                                                      &
                        DOT_PRODUCT_COMPLEX(gradS(:, 1, 1) + gradS(:, 1, 2),SourceDistr%ZIGB(:,Irad),Mesh%Npanels)/2 &
                       +DOT_PRODUCT_COMPLEX(gradS(:, 1, 1) - gradS(:, 1, 2),SourceDistr%ZIGS(:,Irad),Mesh%Npanels)/2
                    RadVelocity(NPFLOW+Ipflow,1,IRad) =                                                              &
                        DOT_PRODUCT_COMPLEX(gradS(:, 1, 1) - gradS(:, 1, 2),SourceDistr%ZIGB(:,Irad),Mesh%Npanels)/2 &
                       +DOT_PRODUCT_COMPLEX(gradS(:, 1, 1) + gradS(:, 1, 2),SourceDistr%ZIGS(:,Irad),Mesh%Npanels)/2
                    !Vy
                    RadVelocity(Ipflow,2,IRad)=                                                                      &
                        DOT_PRODUCT_COMPLEX(gradS(:, 2, 1) + gradS(:, 2, 2),SourceDistr%ZIGB(:,Irad),Mesh%Npanels)/2 &
                       +DOT_PRODUCT_COMPLEX(gradS(:, 2, 1) - gradS(:, 2, 2),SourceDistr%ZIGS(:,Irad),Mesh%Npanels)/2
                    RadVelocity(NPFLOW+Ipflow,2,IRad) =                                                              &
                       -DOT_PRODUCT_COMPLEX(gradS(:, 2, 1) - gradS(:, 2, 2),SourceDistr%ZIGB(:,Irad),Mesh%Npanels)/2 &
                       -DOT_PRODUCT_COMPLEX(gradS(:, 2, 1) + gradS(:, 2, 2),SourceDistr%ZIGS(:,Irad),Mesh%Npanels)/2
                    !Vz
                    RadVelocity(Ipflow,3,IRad)=                                                                      &
                        DOT_PRODUCT_COMPLEX(gradS(:, 3, 1) + gradS(:, 3, 2),SourceDistr%ZIGB(:,Irad),Mesh%Npanels)/2 &
                       +DOT_PRODUCT_COMPLEX(gradS(:, 3, 1) - gradS(:, 3, 2),SourceDistr%ZIGS(:,Irad),Mesh%Npanels)/2
                    RadVelocity(NPFLOW+Ipflow,3,IRad) =                                                              &
                        DOT_PRODUCT_COMPLEX(gradS(:, 3, 1) - gradS(:, 3, 2),SourceDistr%ZIGB(:,Irad),Mesh%Npanels)/2 &
                       +DOT_PRODUCT_COMPLEX(gradS(:, 3, 1) + gradS(:, 3, 2),SourceDistr%ZIGS(:,Irad),Mesh%Npanels)/2
                  ENDIF

            ENDDO

        ENDDO

        DO Ibeta=1,Nbeta

           ILINE=(Iw-1)*Nbeta+Ibeta
           CALL WRITE_POTENTIAL_FS_DATA(wd,IncPotFILE_FS,omega,wavenumber,          &
                   Vbeta(Ibeta), NPFLOW*2**Mesh%Isym,ILINE,IncPotential(:,Ibeta))
           CALL WRITE_POTENTIAL_FS_DATA(wd,TotPotFILE_FS,omega,wavenumber,          &
                   Vbeta(Ibeta), NPFLOW*2**Mesh%Isym,ILINE,Potential(:,Ibeta))
           CALL WRITE_VELOCITY_FS_DATA(wd,IncVelFILE_FS,omega,wavenumber,          &
                   Vbeta(Ibeta), NPFLOW*2**Mesh%Isym,ILINE,IncVelocity(:,:,Ibeta))
           CALL WRITE_VELOCITY_FS_DATA(wd,TotVelFILE_FS,omega,wavenumber,          &
                   Vbeta(Ibeta), NPFLOW*2**Mesh%Isym,ILINE,Velocity(:,:,Ibeta))

           IF (ID_DEBUG==1) THEN
           CALL WRITE_POTENTIAL_FS_DATA_DEBUG(wd,IncPotFILE_FSdbg,omega,wavenumber,          &
                   Vbeta(Ibeta), NPFLOW*2**Mesh%Isym,IncPotential(:,Ibeta))
           CALL WRITE_POTENTIAL_FS_DATA_DEBUG(wd,TotPotFILE_FSdbg,omega,wavenumber,          &
                   Vbeta(Ibeta), NPFLOW*2**Mesh%Isym,Potential(:,Ibeta))
           CALL WRITE_VELOCITY_FS_DATA_DEBUG(wd,IncVelFILE_FSdbg,omega,wavenumber,          &
                   Vbeta(Ibeta), NPFLOW*2**Mesh%Isym,IncVelocity(:,:,Ibeta))
           CALL WRITE_VELOCITY_FS_DATA_DEBUG(wd,TotVelFILE_FSdbg,omega,wavenumber,          &
                   Vbeta(Ibeta), NPFLOW*2**Mesh%Isym,Velocity(:,:,Ibeta))
           ENDIF


        ENDDO
        DEALLOCATE(Potential,Velocity)
        DEALLOCATE(IncPotential,IncVelocity)

        DO IRad=1,Nradiation

           ILINE=(Iw-1)*Nradiation+Irad
           CALL WRITE_RADPOTENTIAL_FS_DATA(wd,RadPotFILE_FS,omega,wavenumber,         &
                   IRad, NPFLOW*2**Mesh%Isym,ILINE,RadPotential(:,IRad))
           CALL WRITE_RADVELOCITY_FS_DATA(wd,RadVelFILE_FS,omega,wavenumber,          &
                   IRad, NPFLOW*2**Mesh%Isym,ILINE,RadVelocity(:,:,IRad))

           IF (ID_DEBUG==1) THEN
           CALL WRITE_RADPOTENTIAL_FS_DATA_DEBUG(wd,RadPotFILE_FSdbg,omega,wavenumber,         &
                   IRad, NPFLOW*2**Mesh%Isym,RadPotential(:,IRad))
           CALL WRITE_RADVELOCITY_FS_DATA_DEBUG(wd,RadVelFILE_FSdbg,omega,wavenumber,          &
                   IRad, NPFLOW*2**Mesh%Isym,RadVelocity(:,:,IRad))
           ENDIF

        ENDDO

        DEALLOCATE(RadPotential,RadVelocity)
        DEALLOCATE(S,GradS)
  END SUBROUTINE

   SUBROUTINE COMPUTE_POTENTIALS_AND_VELOCITIES(wd,Iw,omega,Vbeta,Nbeta,Nradiation,&
                                              Env,Mesh,VFace,WLine,IGreen,         &
                                              SourceDistr,MotionIw)
        !Potential and velocities on body panels
        !Matrix operation, as in NEMOH1
        !INPUT/OUTPUT
        CHARACTER(LEN=*),               INTENT(IN)::wd
        INTEGER,                        INTENT(IN)::Iw,Nbeta,Nradiation
        REAL,                           INTENT(IN)::omega       !rad freq w(Iw)
        REAL,DIMENSION(Nbeta),          INTENT(IN)::Vbeta       !Angle vector
        TYPE(TEnvironment),             INTENT(IN)::Env
        TYPE(TMesh),                    INTENT(IN)::Mesh
        TYPE(TVFace),                   INTENT(IN)::VFace
        TYPE(TWLine),                   INTENT(IN)::WLine
        TYPE(TGREEN),                   INTENT(INOUT)::IGreen
        TYPE(TSource),                  INTENT(IN)::SourceDistr
        COMPLEX,DIMENSION(Nradiation,Nbeta),INTENT(IN):: MotionIw
        !LOCAL
        INTEGER                                ::NPFLOW,uFile,ILINE
        REAL                                   :: wavenumber
        COMPLEX, DIMENSION(:,:,:)  , ALLOCATABLE :: S       ! Inf. coef. Integ. of Green func.
        COMPLEX, DIMENSION(:,:,:,:), ALLOCATABLE :: GradS   ! Inf. coef. Integ. of Gradient Green func.
        INTEGER                                :: Ibeta,Ipanel,Irad
        COMPLEX,DIMENSION(Mesh%Npanels,Nbeta)  :: ZPGB,ZPGS ! Perturbation sources, B for Body
                                                            ! S for symmetric part
        COMPLEX,ALLOCATABLE,DIMENSION(:)       ::Potential  ! Total Potential=Phi_InC+Phi_Perturb
        COMPLEX,ALLOCATABLE,DIMENSION(:,:)     ::Velocity   ! Total velocity
        COMPLEX,ALLOCATABLE,DIMENSION(:)       ::RadPotential ! Rad pot mod I
        COMPLEX,ALLOCATABLE,DIMENSION(:,:)     ::RadVelocity  ! Rad velocity mod I
        CHARACTER(LEN=1000)                    :: LogTextToBeWritten

        INTEGER                                :: RECLength, RECLength2

        NPFLOW=Mesh%NPanels+WLine%NWlineseg
        ALLOCATE(S(NPFLOW,Mesh%NPanels,2**Mesh%Isym))
        ALLOCATE(GradS(NPFLOW,Mesh%NPanels,3,2**Mesh%Isym))
        ALLOCATE(Potential(NPFLOW*2**Mesh%Isym))
        ALLOCATE(Velocity(NPFLOW*2**Mesh%Isym,3))
        wavenumber=Fun_inverseDispersion(omega,Env%depth,Env%g)

        CALL CONSTRUCT_INFLUENCE_MATRIX(omega,wavenumber,Env,IGreen,                  &
                        Mesh,VFace,WLine%XM,WLine%NWlineseg,S,GradS)

        DO Ibeta=1,Nbeta
           !CONSTRUCT PERTURBATION(Diff+Rad) SINGULAR DISTRIBUTION FOR EACH WAVE DIRECTION
           !assign the diffraction singular source distribution for all panels
           ZPGB(1:Mesh%Npanels,Ibeta)=SourceDistr%ZIGB(1:Mesh%Npanels,Nradiation+Ibeta)
           ZPGS(1:Mesh%Npanels,Ibeta)=SourceDistr%ZIGS(1:Mesh%Npanels,Nradiation+Ibeta)
              !sum the diffraction + radiation singular source distribution for all panels
              DO Irad=1,Nradiation
              ZPGB(1:Mesh%Npanels,Ibeta)=ZPGB(1:Mesh%Npanels,Ibeta)                             &
                                 -II*omega*SourceDistr%ZIGB(1:Mesh%Npanels,Irad)                &
                                  *MotionIw(Irad,Ibeta)
              ZPGS(1:Mesh%Npanels,Ibeta)=ZPGS(1:Mesh%Npanels,Ibeta)                             &
                                 -II*omega*SourceDistr%ZIGS(1:Mesh%Npanels,Irad)                &
                                  *MotionIw(Irad,Ibeta)
              ENDDO
              !--------------------------------------------------------------------------------
              ! Compute incoming Potential and velocity
              CALL COMPUTE_INC_POTENTIAL_VELOCITY(wavenumber,omega,Vbeta(Ibeta),                &
                                  VFace%XM,Mesh%Npanels,WLine%XM,WLine%NWlineseg,               &
                                  Env,Mesh%ISym,Potential,Velocity)
              !--------------------------------------------------------------------------------
              ! Compute total potential and velocity
              IF (Mesh%ISym==NO_Y_SYMMETRY) THEN
                Potential(1:NPFLOW)=Potential(1:NPFLOW)+                                        &
                        MATMUL(S(1:NPFLOW,1:Mesh%Npanels,1),ZPGB(1:Mesh%Npanels,Ibeta))
                !Vx
                Velocity(1:NPFLOW,1)=Velocity(1:NPFLOW,1)+                                      &
                        MATMUL(gradS(1:NPFLOW,1:Mesh%Npanels,1,1),ZPGB(1:Mesh%Npanels,Ibeta))
                !Vy
                Velocity(1:NPFLOW,2)=Velocity(1:NPFLOW,2)+                                      &
                        MATMUL(gradS(1:NPFLOW,1:Mesh%Npanels,2,1),ZPGB(1:Mesh%Npanels,Ibeta))
                !Vz
                Velocity(1:NPFLOW,3)=Velocity(1:NPFLOW,3)+                                      &
                        MATMUL(gradS(1:NPFLOW,1:Mesh%Npanels,3,1),ZPGB(1:Mesh%Npanels,Ibeta))

              ELSE IF (Mesh%ISym == Y_SYMMETRY) THEN
                 Potential(1:NPFLOW)=Potential(1:NPFLOW)+                                       &
                   MATMUL(S(:, :, 1) + S(:, :, 2), ZPGB(:,Ibeta))/2                             &
                   +MATMUL(S(:, :, 1) - S(:, :, 2), ZPGS(:,Ibeta))/2
                 Potential(NPFLOW+1:2*NPFLOW) =Potential(NPFLOW+1:2*NPFLOW)+                    &
                   MATMUL(S(:, :, 1) - S(:, :, 2), ZPGB(:,Ibeta))/2                             &
                   +MATMUL(S(:, :, 1) + S(:, :, 2), ZPGS(:,Ibeta))/2
                !Vx
                Velocity(1:NPFLOW,1)=Velocity(1:NPFLOW,1)+                                      &
                   MATMUL(gradS(:, :, 1, 1) + gradS(:, :, 1, 2), ZPGB(:,Ibeta))/2               &
                   +MATMUL(gradS(:, :, 1, 1) - gradS(:, :, 1, 2), ZPGS(:,Ibeta))/2
                 Velocity(NPFLOW+1:2*NPFLOW,1) =Velocity(NPFLOW+1:2*NPFLOW,1)+                  &
                   MATMUL(gradS(:, :, 1, 1) - gradS(:, :, 1, 2), ZPGB(:,Ibeta))/2               &
                   +MATMUL(gradS(:, :, 1, 1) + gradS(:, :, 1, 2), ZPGS(:,Ibeta))/2
                !Vy
                Velocity(1:NPFLOW,2)=Velocity(1:NPFLOW,2)+                                      &
                   MATMUL(gradS(:, :, 2, 1) + gradS(:, :, 2, 2), ZPGB(:,Ibeta))/2               &
                   +MATMUL(gradS(:, :, 2, 1) - gradS(:, :, 2, 2), ZPGS(:,Ibeta))/2
                 Velocity(NPFLOW+1:2*NPFLOW,2) =Velocity(NPFLOW+1:2*NPFLOW,2)-                  &
                   MATMUL(gradS(:, :, 2, 1) - gradS(:, :, 2, 2), ZPGB(:,Ibeta))/2               &
                   -MATMUL(gradS(:, :, 2, 1) + gradS(:, :, 2, 2), ZPGS(:,Ibeta))/2
                !Vz
                Velocity(1:NPFLOW,3)=Velocity(1:NPFLOW,3)+                                      &
                   MATMUL(gradS(:, :, 3, 1) + gradS(:, :, 3, 2), ZPGB(:,Ibeta))/2               &
                   +MATMUL(gradS(:, :, 3, 1) - gradS(:, :, 3, 2), ZPGS(:,Ibeta))/2
                 Velocity(NPFLOW+1:2*NPFLOW,3) =Velocity(NPFLOW+1:2*NPFLOW,3)+                  &
                   MATMUL(gradS(:, :, 3, 1) - gradS(:, :, 3, 2), ZPGB(:,Ibeta))/2               &
                   +MATMUL(gradS(:, :, 3, 1) + gradS(:, :, 3, 2), ZPGS(:,Ibeta))/2
               ! DO Ipanel=1,NPFLOW
               !         print*,Ipanel,Potential(Ipanel)
               ! ENDDO
               ! STOP
              ENDIF
              INQUIRE(iolength=RECLength) (omega)
              INQUIRE(iolength=RECLength2) (REAL(Potential(Ipanel)),Ipanel=1,NPFLOW*2**Mesh%Isym)
              OPEN(NEWUNIT=uFile, FILE=TRIM(wd)//'/'//PreprocDir//'/'//TotPotFILE,              &
                        STATUS='UNKNOWN',ACCESS='DIRECT',RECL=3*RECLength+2*RECLength2)
              ILINE=(Iw-1)*Nbeta+Ibeta
              WRITE(uFile,REC=ILINE) omega,wavenumber,Vbeta(Ibeta),                             &
                             (REAL(Potential(Ipanel)),Ipanel=1,NPFLOW*2**Mesh%Isym),            &
                             (AIMAG(Potential(Ipanel)),Ipanel=1,NPFLOW*2**Mesh%Isym)
              CLOSE(uFile)

              OPEN(NEWUNIT=uFile, FILE=TRIM(wd)//'/'//PreprocDir//'/'//TotVelFILE,              &
                      STATUS='UNKNOWN',ACCESS='DIRECT',RECL=3*RECLength+6*RECLength2)
              WRITE(uFile,REC=ILINE) omega,wavenumber,Vbeta(Ibeta),                             &
                               ( REAL(Velocity(Ipanel,1)),Ipanel=1,NPFLOW*2**Mesh%Isym),        &
                               (AIMAG(Velocity(Ipanel,1)),Ipanel=1,NPFLOW*2**Mesh%Isym),        &
                               ( REAL(Velocity(Ipanel,2)),Ipanel=1,NPFLOW*2**Mesh%Isym),        &
                               (AIMAG(Velocity(Ipanel,2)),Ipanel=1,NPFLOW*2**Mesh%Isym),        &
                               ( REAL(Velocity(Ipanel,3)),Ipanel=1,NPFLOW*2**Mesh%Isym),        &
                               (AIMAG(Velocity(Ipanel,3)),Ipanel=1,NPFLOW*2**Mesh%Isym)
              CLOSE(uFile)
        ENDDO
        DEALLOCATE(Potential,Velocity)

        ! Compute radiation potential and velocity
        ALLOCATE(RadPotential(NPFLOW*2**Mesh%Isym))
        ALLOCATE(RadVelocity(NPFLOW*2**Mesh%Isym,3))
        DO IRad=1,Nradiation
               IF (Mesh%ISym==NO_Y_SYMMETRY) THEN
                RadPotential(1:NPFLOW)=                                                         &
                  MATMUL(S(1:NPFLOW,1:Mesh%Npanels,1),SourceDistr%ZIGB(:,Irad))
                !Vx
                RadVelocity(1:NPFLOW,1)=                                                        &
                  MATMUL(gradS(1:NPFLOW,1:Mesh%Npanels,1,1),SourceDistr%ZIGB(:,Irad))
                !Vy
                RadVelocity(1:NPFLOW,2)=                                                        &
                  MATMUL(gradS(1:NPFLOW,1:Mesh%Npanels,2,1),SourceDistr%ZIGB(:,Irad))
                !Vz
                RadVelocity(1:NPFLOW,3)=                                                        &
                  MATMUL(gradS(1:NPFLOW,1:Mesh%Npanels,3,1),SourceDistr%ZIGB(:,Irad))
              ELSE IF (Mesh%ISym == Y_SYMMETRY) THEN
                 RadPotential(1:NPFLOW)=                                                        &
                   MATMUL(S(:, :, 1) + S(:, :, 2), SourceDistr%ZIGB(:,Irad))/2                  &
                   +MATMUL(S(:, :, 1) - S(:, :, 2),SourceDistr%ZIGS(:,Irad))/2
                 RadPotential(NPFLOW+1:2*NPFLOW) =                                              &
                   MATMUL(S(:, :, 1) - S(:, :, 2), SourceDistr%ZIGB(:,Irad))/2                  &
                   +MATMUL(S(:, :, 1) + S(:, :, 2),SourceDistr%ZIGS(:,Irad))/2
                !Vx
                RadVelocity(1:NPFLOW,1)=                                                        &
                   MATMUL(gradS(:, :, 1, 1) + gradS(:, :, 1, 2), SourceDistr%ZIGB(:,Irad))/2    &
                   +MATMUL(gradS(:, :, 1, 1) - gradS(:, :, 1, 2),SourceDistr%ZIGS(:,Irad))/2
                RadVelocity(NPFLOW+1:2*NPFLOW,1) =                                              &
                   MATMUL(gradS(:, :, 1, 1) - gradS(:, :, 1, 2), SourceDistr%ZIGB(:,Irad))/2    &
                   +MATMUL(gradS(:, :, 1, 1) + gradS(:, :, 1, 2),SourceDistr%ZIGS(:,Irad))/2
                !Vy
                RadVelocity(1:NPFLOW,2)=                                                        &
                   MATMUL(gradS(:, :, 2, 1) + gradS(:, :, 2, 2), SourceDistr%ZIGB(:,Irad))/2    &
                   +MATMUL(gradS(:, :, 2, 1) - gradS(:, :, 2, 2),SourceDistr%ZIGS(:,Irad))/2
                RadVelocity(NPFLOW+1:2*NPFLOW,2) =                                              &
                   -MATMUL(gradS(:, :, 2, 1) - gradS(:, :, 2, 2), SourceDistr%ZIGB(:,Irad))/2   &
                   -MATMUL(gradS(:, :, 2, 1) + gradS(:, :, 2, 2),SourceDistr%ZIGS(:,Irad))/2
                !Vz
                RadVelocity(1:NPFLOW,3)=                                                        &
                   MATMUL(gradS(:, :, 3, 1) + gradS(:, :, 3, 2), SourceDistr%ZIGB(:,Irad))/2    &
                   +MATMUL(gradS(:, :, 3, 1) - gradS(:, :, 3, 2),SourceDistr%ZIGS(:,Irad))/2
                RadVelocity(NPFLOW+1:2*NPFLOW,3) =                                              &
                   MATMUL(gradS(:, :, 3, 1) - gradS(:, :, 3, 2), SourceDistr%ZIGB(:,Irad))/2    &
                   +MATMUL(gradS(:, :, 3, 1) + gradS(:, :, 3, 2),SourceDistr%ZIGS(:,Irad))/2
              ENDIF

              OPEN(NEWUNIT=uFile, FILE=TRIM(wd)//'/'//PreprocDir//'/'//RadPotFILE,              &
                      STATUS='UNKNOWN',ACCESS='DIRECT',RECL=3*RECLength+2*RECLength2)
              ILINE=(Iw-1)*Nradiation+Irad
              WRITE(uFile,REC=ILINE) omega,wavenumber,Irad,                                     &
                             (REAL(RadPotential(Ipanel)),Ipanel=1,NPFLOW*2**Mesh%Isym),         &
                             (AIMAG(RadPotential(Ipanel)),Ipanel=1,NPFLOW*2**Mesh%Isym)
              CLOSE(uFile)


              OPEN(NEWUNIT=uFile, FILE=TRIM(wd)//'/'//PreprocDir//'/'//RadVelFILE,              &
                      STATUS='UNKNOWN',ACCESS='DIRECT',RECL=3*RECLength+6*RECLength2)
              WRITE(uFile,REC=ILINE) omega,wavenumber,Irad,                                     &
                               ( REAL(RadVelocity(Ipanel,1)),Ipanel=1,NPFLOW*2**Mesh%Isym),     &
                               (AIMAG(RadVelocity(Ipanel,1)),Ipanel=1,NPFLOW*2**Mesh%Isym),     &
                               ( REAL(RadVelocity(Ipanel,2)),Ipanel=1,NPFLOW*2**Mesh%Isym),     &
                               (AIMAG(RadVelocity(Ipanel,2)),Ipanel=1,NPFLOW*2**Mesh%Isym),     &
                               ( REAL(RadVelocity(Ipanel,3)),Ipanel=1,NPFLOW*2**Mesh%Isym),     &
                               (AIMAG(RadVelocity(Ipanel,3)),Ipanel=1,NPFLOW*2**Mesh%Isym)
              CLOSE(uFile)

        ENDDO

        DEALLOCATE(RadPotential,RadVelocity)
        DEALLOCATE(S,GradS)
  END SUBROUTINE

  SUBROUTINE WRITE_QTFLOGFILE(wd,beta,Nbeta,w,Nw,NP_GQ,eps_zmin,Nbodies,depth,NPanels,NFSpanels)
        CHARACTER(LEN=*),             INTENT(IN)::wd
        INTEGER,                      INTENT(IN)::Nbeta,Nw,NP_GQ,Nbodies
        INTEGER,                      INTENT(IN)::Npanels,NFSpanels
        REAL, DIMENSION(Nbeta),       INTENT(IN)::beta
        REAL, DIMENSION(Nw),          INTENT(IN)::w
        REAL,                         INTENT(IN)::depth,eps_zmin
        CHARACTER(LEN=1000)                     ::LogTextToBeWritten


        WRITE(*, *) ' '
        WRITE(LogTextToBeWritten,*) '----Pre-Processing (QTF Module)---'
        CALL WRITE_LOGFILE(TRIM(wd)//'/'//LogFILE,TRIM(LogTextToBeWritten),IdStartLog,IdprintTerm)
        WRITE(LogTextToBeWritten,'(A,I4,A,3(F7.3,A))') ' NFreq= ', Nw, ', omega = (', w(1),':',w(2)-w(1),':',w(Nw),') rad/s'
        CALL WRITE_LOGFILE(TRIM(wd)//'/'//LogFILE,TRIM(LogTextToBeWritten),IdAppend,IdprintTerm)
        IF (Nbeta>1) THEN
        WRITE(LogTextToBeWritten,'(A,I4,A,3(F7.2,A))') ' Nbeta= ', Nbeta, ', beta  = (',beta(1)*180/PI, &
                ':',(beta(2)-beta(1))*180/PI,':',beta(Nbeta)*180/PI,') deg'
        ELSE
         WRITE(LogTextToBeWritten,'(A,I4,A,F7.2,A)') ' Nbeta= ', Nbeta, ', beta  = ',beta(1)*180/PI,' deg'
        ENDIF
        CALL WRITE_LOGFILE(TRIM(wd)//'/'//LogFILE,TRIM(LogTextToBeWritten),IdAppend,IdprintTerm)
        WRITE(LogTextToBeWritten,'(A,I3,A)') ' NBodies=', Nbodies, ', NDOF=6'
        CALL WRITE_LOGFILE(TRIM(wd)//'/'//LogFILE,TRIM(LogTextToBeWritten),IdAppend,IdprintTerm)
        WRITE(LogTextToBeWritten,'(A,I7,A)') ' NPanels (Body)        =', NPanels
        CALL WRITE_LOGFILE(TRIM(wd)//'/'//LogFILE,TRIM(LogTextToBeWritten),IdAppend,IdprintTerm)
        WRITE(LogTextToBeWritten,'(A,I7,A)') ' NPanels (Free Surface)=', NFSPanels
        CALL WRITE_LOGFILE(TRIM(wd)//'/'//LogFILE,TRIM(LogTextToBeWritten),IdAppend,IdprintTerm)

        IF (depth==INFINITE_DEPTH) THEN
        WRITE(LogTextToBeWritten,*) 'Waterdepth= Infinite'
        ELSE
        WRITE(LogTextToBeWritten,'(A,F7.3,A)') ' Waterdepth= ', depth,' m'
        ENDIF
        CALL WRITE_LOGFILE(TRIM(wd)//'/'//LogFILE,TRIM(LogTextToBeWritten),IdAppend,IdprintTerm)
         WRITE(LogTextToBeWritten,'(A,I3)') ' NP Gauss Quadrature Integ.: ', NP_GQ
        CALL WRITE_LOGFILE(TRIM(wd)//'/'//LogFILE,TRIM(LogTextToBeWritten),IdAppend,IdprintTerm)
        WRITE(LogTextToBeWritten,*) 'EPS min z                 : ', eps_zmin
         CALL WRITE_LOGFILE(TRIM(wd)//'/'//LogFILE,TRIM(LogTextToBeWritten),IdAppend,IdprintTerm)

  END SUBROUTINE

  SUBROUTINE INITIALIZE_GREEN_FS(IGreen,IGreenFS)
    ! this soubroutine is only passing parameters from IGREEN to IGREENFS
    ! No computaion of the first term influence coefficients

    !INPUT/OUTPUT
                                                        !used in QPreprocessor
    TYPE(TGREEN),                 INTENT(IN)  :: IGreen
    TYPE(TGREEN),                 INTENT(OUT) :: IGreenFS

    ALLOCATE(IGreenFS%XR(IGreen%IR),IGreenFS%XZ(IGreen%JZ))
    ALLOCATE(IGreenFS%APD1X(IGreen%IR,IGreen%JZ),IGreenFS%APD2X(IGreen%IR,IGreen%JZ))
    ALLOCATE(IGreenFS%APD1Z(IGreen%IR,IGreen%JZ),IGreenFS%APD2Z(IGreen%IR,IGreen%JZ))

    !passing parameters
    IGreenFS%IR       =IGreen%IR
    IGreenFS%JZ       =IGreen%NPINTE
    IGreenFS%MAX_KR   =IGreen%MAX_KR
    IGreenFS%MIN_KZ   =IGreen%MIN_KZ
    IGreenFS%EPS_ZMIN =IGreen%EPS_ZMIN
    IGreenFS%NEXP     =IGreen%NEXP
    IGreenFS%AMBDA    =IGreen%AMBDA
    IGreenFS%AR       =IGreen%AR

    IGreenFS%XR       =IGreen%XR
    IGreenFS%XZ       =IGreen%XZ
    IGreenFS%APD1X    =IGreen%APD1X
    IGreenFS%APD1Z    =IGreen%APD1Z
    IGreenFS%APD2X    =IGreen%APD2X
    IGreenFS%APD2Z    =IGreen%APD2Z
   END SUBROUTINE

   SUBROUTINE INITIALIZE_POTVELFS_OUTPUT_FILES(wd,Npanels,XM,Isym,NBdyLine,XM_Line)
        CHARACTER(LEN=*),          INTENT(IN)::wd
        INTEGER,                   INTENT(IN)::Npanels,Isym,NbdyLine
        REAL,DIMENSION(3,Npanels), INTENT(IN)::XM
        REAL,DIMENSION(NBdyLine,3),INTENT(IN)::XM_Line
        INTEGER                              ::Ipanel,u,ITERFILE
        CHARACTER(LEN=100)        :: FileName

        DO IterFile=1,6
        IF (IterFile==1) FileName=IncPotFILE_FSdbg
        IF (IterFile==2) FileName=TotPotFILE_FSdbg
        IF (IterFile==3) FileName=RadPotFILE_FSdbg
        IF (IterFile==4) FileName=IncVelFILE_FSdbg
        IF (IterFile==5) FileName=TotVelFILE_FSdbg
        IF (IterFile==6) FileName=RadVelFILE_FSdbg

        OPEN(NEWUNIT=u,FILE=TRIM(wd)//'/'//PreprocDir//'/'//TRIM(FileName), &
                                                              ACTION='WRITE')

        WRITE(u,*) '#Npanels,NBdyLine'
        WRITE(u,*) '#[XM(1,1:Npanels),XM(2,1:Npanels),XM(3,1:Npanels)]'
        IF (IterFile<=3) THEN
        WRITE(u,'(4(X,A))') '#w_i','beta_i','Real(Pot_i)','Imag(Pot_i)'
        ELSE
        WRITE(u,'(8(X,A))')'#w_i','beta_i','Real(Velx_i)','Imag(Velx_i)',&
                'Real(Vely_i)','Imag(Vely_i)','Real(Velz_i)','Imag(Velz_i)'
        ENDIF
        WRITE(u,*) '#'
        IF (((3*Npanels).GT.10000).OR.((3*NBdyLine).GT.10000)) THEN
          PRINT*,'Increase the value used in the following write commands'
          STOP
        ENDIF
        IF (Isym==0) THEN
           WRITE(u,'(3(I6,X))') Npanels,NBdyLine,Isym
           IF (NBdyLine==0) THEN
           WRITE(u,'(10000(E14.7,X))')              &
                 (XM(1,Ipanel),Ipanel=1,Npanels),         &
                 (XM(2,Ipanel),Ipanel=1,Npanels),         &
                 (XM(3,Ipanel),Ipanel=1,Npanels)
           ELSE
          WRITE(u,'(10000(E14.7,X),10000(E14.7,X))') &
                (XM(1,Ipanel),Ipanel=1,Npanels),         &
                (XM(2,Ipanel),Ipanel=1,Npanels),         &
                (XM(3,Ipanel),Ipanel=1,Npanels),         &
                (XM_Line(Ipanel,1),Ipanel=1,NBdyLine),   &
                (XM_Line(Ipanel,2),Ipanel=1,NBdyLine),   &
                (XM_Line(Ipanel,3),Ipanel=1,NBdyLine)
           ENDIF
        ELSE
           WRITE(u,'(3(I6,X))') Npanels,NBdyLine,Isym
           IF (NBdyLine==0) THEN
           WRITE(u,'(20000(E14.7,X))')            &
               (XM(1,Ipanel),Ipanel=1,Npanels),           &
               (XM(2,Ipanel),Ipanel=1,Npanels),           &
               (XM(3,Ipanel),Ipanel=1,Npanels),           &
               (XM(1,Ipanel),Ipanel=1,Npanels),           &
               (-XM(2,Ipanel),Ipanel=1,Npanels),          &
               (XM(3,Ipanel),Ipanel=1,Npanels)
           ELSE
           WRITE(u,'(40000(E14.7,X))' )&
               (XM(1,Ipanel),Ipanel=1,Npanels),           &
               (XM(2,Ipanel),Ipanel=1,Npanels),           &
               (XM(3,Ipanel),Ipanel=1,Npanels),           &
               (XM_Line(Ipanel,1),Ipanel=1,NBdyLine),     &
               (XM_Line(Ipanel,2),Ipanel=1,NBdyLine),     &
               (XM_Line(Ipanel,3),Ipanel=1,NBdyLine),     &
               (XM(1,Ipanel),Ipanel=1,Npanels),           &
               (-XM(2,Ipanel),Ipanel=1,Npanels),          &
               (XM(3,Ipanel),Ipanel=1,Npanels),           &
               (XM_Line(Ipanel ,1),Ipanel=1,NBdyLine),     &
               (-XM_Line(Ipanel,2),Ipanel=1,NBdyLine),    &
               (XM_Line(Ipanel ,3),Ipanel=1,NBdyLine)
           ENDIF

        ENDIF
        CLOSE(u)
        ENDDO
   END SUBROUTINE

   SUBROUTINE WRITE_POTENTIAL_FS_DATA(wd,FileM,w,k,beta,Npanels,Iline,Potential)
         CHARACTER(LEN=*),    INTENT(IN)::wd,FileM
         INTEGER,             INTENT(IN)::Npanels,Iline
         REAL,                INTENT(IN)::w,k,beta
         COMPLEX,DIMENSION(Npanels),INTENT(IN):: Potential
         Integer :: Ipanel,uFile,RECLength

         INQUIRE(iolength=RECLength) w
         OPEN(NEWUNIT=uFile, FILE=TRIM(wd)//'/'//PreprocDir//'/'//FileM,  &
                 STATUS='UNKNOWN',ACCESS='DIRECT',RECL=(3+2*Npanels)*RECLength)
         WRITE(uFile,REC=ILINE) w,k,beta,                             &
                        (REAL(Potential(Ipanel)),Ipanel=1,Npanels), &
                        (AIMAG(Potential(Ipanel)),Ipanel=1,Npanels)
         CLOSE(uFile)
    END SUBROUTINE

    SUBROUTINE WRITE_VELOCITY_FS_DATA(wd,FileM,w,k,beta,Npanels,Iline,Velocity)
         CHARACTER(LEN=*),    INTENT(IN)::wd,FileM
         INTEGER,             INTENT(IN)::Npanels,Iline
         REAL,                INTENT(IN)::w,k,beta
         COMPLEX,DIMENSION(Npanels,3),INTENT(IN):: Velocity
         Integer :: Ipanel,uFile,RECLength

         INQUIRE(iolength=RECLength) w
         OPEN(NEWUNIT=uFile, FILE=TRIM(wd)//'/'//PreprocDir//'/'//FileM,   &
                   STATUS='UNKNOWN',ACCESS='DIRECT',RECL=(3+3*2*Npanels)*RECLength)

         WRITE(uFile,REC=ILINE) w,k,beta,                              &
                            ( REAL(Velocity(Ipanel,1)),Ipanel=1,Npanels),  &
                            (AIMAG(Velocity(Ipanel,1)),Ipanel=1,Npanels),  &
                            ( REAL(Velocity(Ipanel,2)),Ipanel=1,Npanels),  &
                            (AIMAG(Velocity(Ipanel,2)),Ipanel=1,Npanels),  &
                            ( REAL(Velocity(Ipanel,3)),Ipanel=1,Npanels),  &
                            (AIMAG(Velocity(Ipanel,3)),Ipanel=1,Npanels)
         CLOSE(uFile)
   END SUBROUTINE

   SUBROUTINE WRITE_RADPOTENTIAL_FS_DATA(wd,FileM,w,k,IRad,Npanels,Iline,Potential)
         CHARACTER(LEN=*),    INTENT(IN)::wd,FileM
         INTEGER,             INTENT(IN)::Npanels,Iline,IRad
         REAL,                INTENT(IN)::w,k
         COMPLEX,DIMENSION(Npanels),INTENT(IN):: Potential
         Integer :: Ipanel,uFile,RECLength

         INQUIRE(iolength=RECLength) w
         OPEN(NEWUNIT=uFile, FILE=TRIM(wd)//'/'//PreprocDir//'/'//FileM,  &
                 STATUS='UNKNOWN',ACCESS='DIRECT',RECL=(3+2*Npanels)*RECLength)
         WRITE(uFile,REC=ILINE) w,k,IRad,                             &
                        (REAL(Potential(Ipanel)),Ipanel=1,Npanels), &
                        (AIMAG(Potential(Ipanel)),Ipanel=1,Npanels)
         CLOSE(uFile)
    END SUBROUTINE

    SUBROUTINE WRITE_RADVELOCITY_FS_DATA(wd,FileM,w,k,IRad,Npanels,Iline,Velocity)
         CHARACTER(LEN=*),    INTENT(IN)::wd,FileM
         INTEGER,             INTENT(IN)::Npanels,Iline,IRad
         REAL,                INTENT(IN)::w,k
         COMPLEX,DIMENSION(Npanels,3),INTENT(IN):: Velocity
         Integer :: Ipanel,uFile,RECLength

         INQUIRE(iolength=RECLength) w
         OPEN(NEWUNIT=uFile, FILE=TRIM(wd)//'/'//PreprocDir//'/'//FileM,   &
                   STATUS='UNKNOWN',ACCESS='DIRECT',RECL=(3+3*2*Npanels)*RECLength)

         WRITE(uFile,REC=ILINE) w,k,IRad,                              &
                            ( REAL(Velocity(Ipanel,1)),Ipanel=1,Npanels),  &
                            (AIMAG(Velocity(Ipanel,1)),Ipanel=1,Npanels),  &
                            ( REAL(Velocity(Ipanel,2)),Ipanel=1,Npanels),  &
                            (AIMAG(Velocity(Ipanel,2)),Ipanel=1,Npanels),  &
                            ( REAL(Velocity(Ipanel,3)),Ipanel=1,Npanels),  &
                            (AIMAG(Velocity(Ipanel,3)),Ipanel=1,Npanels)
         CLOSE(uFile)
   END SUBROUTINE

   SUBROUTINE WRITE_POTENTIAL_FS_DATA_DEBUG(wd,FileM,w,k,beta,Npanels,PotDat)
         CHARACTER(LEN=*),    INTENT(IN)::wd,FileM
         INTEGER,             INTENT(IN)::Npanels
         REAL,                INTENT(IN)::w,k,beta
         COMPLEX,DIMENSION(Npanels),INTENT(IN):: PotDat
         Integer :: Ipanel,u1

         IF ((2*Npanels).GT.10000) THEN
           PRINT*,'Increase the value used in the following write commands'
           STOP
         ENDIF
         OPEN(NEWUNIT=u1, FILE=wd//'/'//PreprocDir//'/'//FileM,    &
                 ACTION='WRITE',POSITION='APPEND')
         WRITE(u1,'(3(F10.3,X),10000(E14.7,X))') w,k,beta,   &
         (REAL(PotDat(Ipanel)),Ipanel=1,Npanels),&
         (AIMAG(PotDat(Ipanel)),Ipanel=1,Npanels)
         CLOSE(u1)
   END SUBROUTINE

    SUBROUTINE WRITE_VELOCITY_FS_DATA_DEBUG(wd,FileM,w,k,beta,Npanels,VelDat)
         CHARACTER(LEN=*),    INTENT(IN)::wd,FileM
         INTEGER,             INTENT(IN)::Npanels
         REAL,                INTENT(IN)::w,k,beta
         COMPLEX,DIMENSION(Npanels,3),INTENT(IN):: VelDat
         Integer :: Ipanel,u1

         IF ((3*Npanels).GT.10000) THEN
           PRINT*,'Increase the value used in the following write commands'
           STOP
         ENDIF
         OPEN(NEWUNIT=u1, FILE=wd//'/'//PreprocDir//'/'//FileM,    &
                 ACTION='WRITE',POSITION='APPEND')
         WRITE(u1,'(3(F10.3,X),20000(E14.7,X))') w,k,beta,       &
         (REAL(VelDat(Ipanel,1)),Ipanel=1,Npanels),&
         (AIMAG(VelDat(Ipanel,1)),Ipanel=1,Npanels),&
         (REAL(VelDat(Ipanel,2)),Ipanel=1,Npanels),&
         (AIMAG(VelDat(Ipanel,2)),Ipanel=1,Npanels),&
         (REAL(VelDat(Ipanel,3)),Ipanel=1,Npanels),&
         (AIMAG(VelDat(Ipanel,3)),Ipanel=1,Npanels)
         CLOSE(u1)
   END SUBROUTINE

   SUBROUTINE WRITE_RADPOTENTIAL_FS_DATA_DEBUG(wd,FileM,w,k,Irad,Npanels,PotDat)
         CHARACTER(LEN=*),    INTENT(IN)::wd,FileM
         INTEGER,             INTENT(IN)::Npanels
         REAL,                INTENT(IN)::w,k
         INTEGER,             INTENT(IN)::Irad
         COMPLEX,DIMENSION(Npanels),INTENT(IN):: PotDat
         Integer :: Ipanel,u1

         IF ((2*Npanels).GT.10000) THEN
           PRINT*,'Increase the value used in the following write commands'
           STOP
         ENDIF
         OPEN(NEWUNIT=u1, FILE=wd//'/'//PreprocDir//'/'//FileM,    &
                 ACTION='WRITE',POSITION='APPEND')
         WRITE(u1,'(2(F10.3,X),(I3,X),10000(E14.7,X))') w,k,Irad,       &
         (REAL(PotDat(Ipanel)),Ipanel=1,Npanels),&
         (AIMAG(PotDat(Ipanel)),Ipanel=1,Npanels)
         CLOSE(u1)
   END SUBROUTINE

    SUBROUTINE WRITE_RADVELOCITY_FS_DATA_DEBUG(wd,FileM,w,k,Irad,Npanels,VelDat)
         CHARACTER(LEN=*),    INTENT(IN)::wd,FileM
         INTEGER,             INTENT(IN)::Npanels
         REAL,                INTENT(IN)::w,k
         INTEGER,             INTENT(IN):: Irad
         COMPLEX,DIMENSION(Npanels,3),INTENT(IN):: VelDat
         Integer :: Ipanel,u1
         IF ((3*Npanels).GT.10000) THEN
           PRINT*,'Increase the value used in the following write commands'
           STOP
         ENDIF
         OPEN(NEWUNIT=u1, FILE=wd//'/'//PreprocDir//'/'//FileM,    &
                 ACTION='WRITE',POSITION='APPEND')
         WRITE(u1,'(2(F10.3,X),(I3,X),20000(E14.7,X))') w,k,Irad,       &
         (REAL(VelDat(Ipanel,1)),Ipanel=1,Npanels),&
         (AIMAG(VelDat(Ipanel,1)),Ipanel=1,Npanels),&
         (REAL(VelDat(Ipanel,2)),Ipanel=1,Npanels),&
         (AIMAG(VelDat(Ipanel,2)),Ipanel=1,Npanels),&
         (REAL(VelDat(Ipanel,3)),Ipanel=1,Npanels),&
         (AIMAG(VelDat(Ipanel,3)),Ipanel=1,Npanels)
         CLOSE(u1)
   END SUBROUTINE



END MODULE
