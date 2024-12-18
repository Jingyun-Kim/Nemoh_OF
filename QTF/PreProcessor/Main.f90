!-----------------------------------------------------------------------------------
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
!-------------------------------------------------------------------------------------
!   Contributors list:
!   - Gerard Delhommeau (13/11/2014, d'apres LA VERSION 1.1 d'AVRIL 1991
!     PROGRAMME StOK LABORATOIRE D'HYDRODYNAMIQUE NAVALE DE L'E.N.S.M. DE NANTES & SIREHNA )
!   - Adrien Combourieu, INNOSEA (adrien.combourieu@innosea.fr)  Version 2014
!   - Ruddy Kurnia (LHEEA,ECN) 2022
!--------------------------------------------------------------------------------------
!
!   NEMOH2 - second order (QTF) - May 2022
!
!--------------------------------------------------------------------------------------

PROGRAM Main
!
USE MIdentification
USE MNemohCal,          ONLY:TNemCal,READ_TNEMOHCAL,Tqtfinput,               &
                             Discretized_Omega_and_Beta
USE MMesh
USE MFace,              ONLY:TVFace, Prepare_FaceMesh,TWLine,Prepare_Waterline
USE MReadInputFiles,    ONLY:Read_NP_GaussQuad,Read_Mechanical_Coefs,TMech,  &
                             Read_FirstOrderLoad,TLoad1,Read_Motion,TSource, &
                             Read_SourceDistribution,Read_Eps_Zmin,          &
                             TMeshFS,Read_Prepare_FreeSurface_Mesh
USE M_INITIALIZE_GREEN, ONLY: TGREEN, INITIALIZE_GREEN
USE MQpreprocessor
USE MLogFile


IMPLICIT NONE
!
!Declaration variables
!
        TYPE(TID)                       :: ID
        TYPE(TMesh)                     :: Mesh
        TYPE(TMeshFS)                   :: MeshFS
        TYPE(TNemCal)                   :: inpNEMOHCAL
        TYPE(Tqtfinput)                 :: QTFinputNem
        INTEGER                         :: Nw,Nbeta           ! Number of Freq,direction
        REAL, ALLOCATABLE,DIMENSION(:)  :: w,beta             ! vector of freq [rad/s],
                                                              ! direction angle [rad]
        TYPE(TWLine)  :: WLine          ! Waterline
        TYPE(TVFace)  :: VFace          ! Face of body panel
        TYPE(TMech)   :: MechCoef       ! MechCoef
        TYPE(TLoad1)  :: Forces1        ! First order forces
        INTEGER       :: NP_GQ          ! Number of point for Gauss Quad. Integration
        INTEGER       :: Nintegration   ! Number of excitation force integration
        INTEGER       :: Nradiation     ! Number of radiation problem
        COMPLEX,ALLOCATABLE,DIMENSION(:,:,:) :: Motion
        TYPE(TSource) :: SOURCEDISTR                            !First order NEMOH solution
        TYPE(TGREEN)                         :: IGreen          ! Initial Green variables for body
        TYPE(TGREEN)                         :: IGreenFS        ! Initial Green variables for free surface
        REAL                                 :: EPS_ZMIN
        INTEGER       :: I,J,uFile
        CHARACTER(LEN=1000)                :: LogTextToBeWritten
        REAL                               :: tcpu_start
!
!   --- Initialize and read input datas -----------------------------------------------
!
        CALL ReadTID(ID)
        CALL ReadTMesh(Mesh,TRIM(ID%ID)//'/mesh/')
        CALL READ_TNEMOHCAL(ID,InpNEMOHCAL)
        CALL Read_Mechanical_Coefs(TRIM(ID%ID),InpNEMOHCAL%Nbodies,MechCoef)
!
        Nw          =InpNEMOHCAL%waveinput%NFreq
        Nbeta       =InpNEMOHCAL%waveinput%NBeta
        Nintegration=InpNEMOHCAL%Nintegtot
        Nradiation  =InpNEMOHCAL%Nradtot
        QTFinputNem =InpNEMOHCAL%qtfinput
        CALL Read_FirstOrderLoad(TRIM(ID%ID),Nw,Nbeta,Nintegration,Nradiation,Forces1)
        ALLOCATE(Motion(Nw,Nradiation,Nbeta))
        CALL Read_Motion(TRIM(ID%ID),Nw,Nbeta,Nradiation,Motion)!RAO
!
        ALLOCATE(w(Nw),beta(Nbeta))
        CALL Discretized_Omega_and_Beta(1,InpNEMOHCAL%waveinput,Nw,Nbeta,w,beta)
!
        NP_GQ=Read_NP_GaussQuad(TRIM(ID%ID))
        EPS_ZMIN=Read_Eps_Zmin(TRIM(ID%ID))
!       Prepare Body Mesh
        CALL Prepare_FaceMesh(Mesh,NP_GQ,VFace)
        CALL Prepare_Waterline(VFace,EPS_ZMIN,Mesh%xy_diameter,Mesh%Npanels,WLine)
!
        CALL INITIALIZE_GREEN(VFace,Mesh,InpNEMOHCAL%Env%depth, &
                              WLine%XM,WLine%NWlineseg,EPS_ZMIN,IGreen)
!       -------------------------------------
!       Prepare Free-surface mesh
        IF  (InpNEMOHCAL%qtfinput%Ncontrib==3) THEN
          CALL Read_Prepare_FreeSurface_Mesh(TRIM(ID%ID),MeshFS,QTFinputNem)
          CALL INITIALIZE_GREEN_FS(IGreen,IGreenFS)
        ENDIF
!       --------------------------------------
!
        IF  (InpNEMOHCAL%qtfinput%Ncontrib==3) THEN
        CALL WRITE_QTFLOGFILE(TRIM(ID%ID),beta,Nbeta,w,Nw,NP_GQ,EPS_ZMIN,                      &
                InpNEMOHCAL%Nbodies,InpNEMOHCAL%Env%depth,Mesh%Npanels,MeshFS%Mesh%Npanels)
        ELSE
        CALL WRITE_QTFLOGFILE(TRIM(ID%ID),beta,Nbeta,w,Nw,NP_GQ,EPS_ZMIN,                      &
                InpNEMOHCAL%Nbodies,InpNEMOHCAL%Env%depth,Mesh%Npanels,0)
        ENDIF

        CALL START_RECORD_TIME(tcpu_start,TRIM(ID%ID)//'/'//LogFILE,IdAppend)
        WRITE(LogTextToBeWritten,*) '-------'
        CALL WRITE_LOGFILE(TRIM(ID%ID)//'/'//LogFILE,TRIM(LogTextToBeWritten),IdAppend,IdprintTerm)

        ALLOCATE(SOURCEDISTR%ZIGB(Mesh%Npanels,Nradiation+Nbeta))
        ALLOCATE(SOURCEDISTR%ZIGS(Mesh%Npanels,Nradiation+Nbeta))

        CALL make_directory(TRIM(ID%ID)//'/'//PreprocDir)

        IF  (InpNEMOHCAL%qtfinput%Ncontrib==3) THEN
             IF (ID_DEBUG==1) CALL INITIALIZE_POTVELFS_OUTPUT_FILES                         &
                                (ID%ID,MeshFS%Mesh%Npanels,MeshFS%Mesh%XM,MeshFS%Mesh%ISym, &
                                 MeshFS%BdyLine%NWlineseg,MeshFS%BdyLine%XM)
        ENDIF
! ------Computing potentials and velocities--------------------------------------------------
        DO I=1,Nw
            CALL Read_SourceDistribution(TRIM(ID%ID),I,Nw,Nradiation,Nbeta,                 &
                                         Mesh%Npanels,SourceDistr)
            ! Calc Body Potentials
            CALL COMPUTE_POTENTIALS_AND_VELOCITIES(TRIM(ID%ID),                             &
                                         I,w(I),beta,Nbeta,Nradiation,InpNEMOHCAL%Env,Mesh, &
                                         VFace,WLine,IGreen,SourceDistr,Motion(I,:,:))

           IF  (InpNEMOHCAL%qtfinput%Ncontrib==3) THEN
            ! Calc Free Surface Potentials
            CALL COMPUTE_POTENTIALS_AND_VELOCITIES_FS(TRIM(ID%ID),                          &
                                         I,w(I),beta,Nbeta,Nradiation,InpNEMOHCAL%Env,Mesh, &
                                         VFace,IGreenFS,SourceDistr,Motion(I,:,:),MeshFS)
           ENDIF

           IF ((InpNEMOHCAL%Env%depth == INFINITE_DEPTH).OR.                                &
                   (w(I)**2*InpNEMOHCAL%Env%depth/InpNEMOHCAL%Env%g.GE.20)) THEN
             WRITE(LogTextToBeWritten,'(A,F7.3,A)') 'Omega= ', w(I),' [rad/s],              &
                   Green Fun: Infinite-Depth. DONE!'
           ELSE
             WRITE(LogTextToBeWritten,'(A,F7.3,A)') 'Omega= ', w(I),' [rad/s],              &
                   Green Fun: Finite-Depth. DONE!'
           ENDIF
           CALL WRITE_LOGFILE(TRIM(ID%ID)//'/'//LogFILE,TRIM(LogTextToBeWritten),IdAppend,IdprintTerm)
        ENDDO
        WRITE(LogTextToBeWritten,*) '-------'
        CALL WRITE_LOGFILE(TRIM(ID%ID)//'/'//LogFILE,TRIM(LogTextToBeWritten),IdAppend,IdprintTerm)
        CALL END_RECORD_TIME(tcpu_start,TRIM(ID%ID)//'/'//LogFILE)
        WRITE(LogTextToBeWritten,*) '---- DONE ---'
        CALL WRITE_LOGFILE(TRIM(ID%ID)//'/'//LogFILE,TRIM(LogTextToBeWritten),IdAppend,IdprintTerm)
! ----- Finalize ---------------------------------------------------------------------------
!       DEALOCATING variables
        DEALLOCATE(w,beta)
        DO I=1,InpNEMOHCAL%Nbodies
                DEALLOCATE(inpNEMOHCAL%bodyinput(I)%RadCase)
                DEALLOCATE(inpNEMOHCAL%bodyinput(I)%IntCase)
        ENDDO
        DEALLOCATE(inpNEMOHCAL%bodyinput)

        DEALLOCATE(VFace%X,VFace%XM,VFace%N,VFace%A,VFace%tDis)
        DEALLOCATE(VFace%dXdXG_WGQ_per_A,VFace%XM_GQ)
        DEALLOCATE(WLine%XM,WLine%SegLength,WLine%IndexPanel)
        DEALLOCATE(MechCoef%MassMat,MechCoef%StiffMat)
        DEALLOCATE(MechCoef%StiffMat_EXT,MechCoef%DampCoefMat_EXT)
        DEALLOCATE(Forces1%addedmass,Forces1%dampcoef)
        DEALLOCATE(Forces1%excitation)
        DEALLOCATE(Motion)
        DEALLOCATE(SOURCEDISTR%ZIGB,SOURCEDISTR%ZIGS)
        DEALLOCATE(IGreen%FSP1,IGreen%FSM1,IGreen%VSP1,IGreen%VSM1)
        DEALLOCATE(IGreen%FSP1_INF,IGreen%FSM1_INF,IGreen%VSP1_INF,IGreen%VSM1_INF)
        DEALLOCATE(IGreen%XR,IGreen%XZ)
        DEALLOCATE(IGreen%APD1X,IGreen%APD2X,IGreen%APD1Z,IGreen%APD2Z)
        IF  (InpNEMOHCAL%qtfinput%Ncontrib==3) THEN
        DEALLOCATE(IGreenFS%XR,IGreenFS%XZ)
        DEALLOCATE(IGreenFS%APD1X,IGreenFS%APD2X,IGreenFS%APD1Z,IGreenFS%APD2Z)
        END IF
END PROGRAM Main
