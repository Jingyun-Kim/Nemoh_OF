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
!
!   Contributors list:
!   - R. Kurnia
!
!--------------------------------------------------------------------------------------

        MODULE MNemohCal

          USE Constants
          USE MIdentification
          USE MEnvironment, only: TEnvironment
          IMPLICIT NONE

          PUBLIC      :: READ_TNEMOHCAL

          INTEGER, PARAMETER :: IdRadFreq=1     !unit rad/s
          INTEGER, PARAMETER :: IdFreqHz =2     !unit Hz
          INTEGER, PARAMETER :: IdPeriod =3     !unit s

          TYPE Twaveinput
              INTEGER  :: FreqType
              INTEGER  :: NFreq,NBeta
              REAL     :: Freq1,Freq2,Beta1,Beta2
          END TYPE Twaveinput

          TYPE TIRF
              INTEGER  :: Switch
              REAL     :: time_step,duration
          END TYPE TIRF

          TYPE TKochin
              INTEGER  :: Switch,Ntheta
              REAL     :: min_theta,max_theta   !degree
          END TYPE TKochin

          TYPE TFreesurface
              INTEGER     :: Switch
              INTEGER     :: NX,NY      !Number of points
              REAL        :: Lx,LY      !Length of domain
          END TYPE TFreesurface

          TYPE TOptOutput
              INTEGER            :: Switch_POTENTIAL
              INTEGER            :: Switch_SourceDistr
              INTEGER            :: Switch_RAO
              INTEGER            :: FreqType
              TYPE(TIRF)         :: IRF
              TYPE(TKochin)      :: Kochin
              TYPE(TFreesurface) :: Freesurface
          END TYPE TOptOutput

          TYPE TBCase
              INTEGER :: ICase
              REAL,DIMENSION(3) :: Direction,Axis
          END TYPE TBCase

          TYPE Tbodyinput
               CHARACTER(LEN=100)       :: meshfile
               INTEGER                  :: Npoints,NPanels
               INTEGER                  :: NRadiation
               INTEGER                  :: NIntegration
               TYPE(TBCase),ALLOCATABLE :: RadCase(:)
               TYPE(TBCase),ALLOCATABLE :: IntCase(:)
          END TYPE Tbodyinput

          TYPE Tqtfinput
              INTEGER  :: switch_qtfp       !LQTFP
              REAL,DIMENSION(3) :: omega    !rad freq [Nw,wmin,wmax]
              REAL     :: body_forward_speed!body forward-speed
              INTEGER  :: bidirection       !bi-direction
              INTEGER  :: NContrib          !Contrib
              CHARACTER(LEN=100) :: FSmeshfile !Free surface meshfile
              REAL     :: FSRe              !Free surface exterior radius
              INTEGER  :: FSNRe             !Number of points in Free surface exterior radius
              INTEGER  :: FSNBessel         !Number of bessel function
              INTEGER  :: FreqTypeOutput    !Frq type output
              INTEGER  :: switch_quadHM     !quadratic hydrostatic and moment terms
              INTEGER  :: switch_qtfduok    !Loutduok
              INTEGER  :: switch_qtfhasbo   !Louthasbo
              INTEGER  :: switch_qtfhasfs   !Louthasfs
          END TYPE Tqtfinput

          TYPE TNemCal
              TYPE(Tenvironment)           :: Env
              INTEGER                      :: Nbodies,Nradtot,Nintegtot
              TYPE(Tbodyinput),ALLOCATABLE :: bodyinput(:)
              TYPE(Twaveinput)             :: waveinput
              TYPE(TOptOutput)             :: OptOUTPUT
              TYPE(Tqtfinput)              :: qtfinput
          END TYPE TNemCal

        CONTAINS

         SUBROUTINE READ_TNEMOHCAL(ID,InpNEMOHCAL)

          TYPE(TID),            INTENT(IN)      :: ID
          TYPE(TNemCal),        INTENT(OUT)     :: InpNEMOHCAL

          !Local var
          INTEGER ufile,I,J,K
          LOGICAL :: exist_dir

          OPEN(NEWUNIT=ufile, FILE=TRIM(ID%ID)//'/Nemoh.cal')
           READ(ufile,*) !--- Environment ----------------------------!
           READ(ufile,*) InpNEMOHCAL%Env%RHO
           READ(ufile,*) InpNEMOHCAL%Env%G
           READ(ufile,*) InpNEMOHCAL%Env%Depth
           READ(ufile,*) InpNEMOHCAL%Env%Xeff, InpNEMOHCAL%Env%Yeff
           READ(ufile,*) !--- Description of floating bodies-----------!
           READ(ufile,*) InpNEMOHCAL%Nbodies
           !! ALLOCATE body input variable
           ALLOCATE(InpNEMOHCAL%bodyinput(InpNEMOHCAL%Nbodies))
           !!
           InpNEMOHCAL%Nradtot=0
           InpNEMOHCAL%Nintegtot=0

           DO I=1,InpNEMOHCAL%Nbodies
             READ(ufile,*) !----Body(I)--------------------------------!
             READ(ufile,*) InpNEMOHCAL%bodyinput(I)%meshfile
             READ(ufile,*) InpNEMOHCAL%bodyinput(I)%Npoints,           &
                                  InpNEMOHCAL%bodyinput(I)%Npanels
             READ(ufile,*) InpNEMOHCAL%bodyinput(I)%NRadiation
             InpNEMOHCAL%Nradtot=InpNEMOHCAL%Nradtot+                  &
                                 InpNEMOHCAL%bodyinput(I)%NRadiation
             ALLOCATE(InpNEMOHCAL%bodyinput(I)%RadCase(                &
                              InpNEMOHCAL%bodyinput(I)%NRadiation))
             DO J=1,InpNEMOHCAL%bodyinput(I)%NRadiation
              READ(ufile,*)InpNEMOHCAL%bodyinput(I)%RadCase(J)%ICase,  &
              (InpNEMOHCAL%bodyinput(I)%RadCase(J)%Direction(K),K=1,3),&
              (InpNEMOHCAL%bodyinput(I)%RadCase(J)%Axis(K),K=1,3)
             END DO
             READ(ufile,*) InpNEMOHCAL%bodyinput(I)%NIntegration
             InpNEMOHCAL%Nintegtot=InpNEMOHCAL%Nintegtot+              &
                                 InpNEMOHCAL%bodyinput(I)%NIntegration
             ALLOCATE(InpNEMOHCAL%bodyinput(I)%IntCase(                &
                              InpNEMOHCAL%bodyinput(I)%NIntegration))
             DO J=1,InpNEMOHCAL%bodyinput(I)%NIntegration
              READ(ufile,*)InpNEMOHCAL%bodyinput(I)%IntCase(J)%ICase,  &
              (InpNEMOHCAL%bodyinput(I)%IntCase(J)%Direction(K),K=1,3),&
              (InpNEMOHCAL%bodyinput(I)%IntCase(J)%Axis(K),K=1,3)
             END DO
             READ(ufile,*) !For additional info, ie. for generalize mode, Not implemented yet
           END DO
             READ(ufile,*) !--- Load cases to be solved ---------------!
             READ(ufile,*) InpNEMOHCAL%waveinput%FreqType,             &
                           InpNEMOHCAL%waveinput%NFreq,                &
                           InpNEMOHCAL%waveinput%Freq1,                &
                           InpNEMOHCAL%waveinput%Freq2
             READ(ufile,*) InpNEMOHCAL%waveinput%NBeta,                &
                           InpNEMOHCAL%waveinput%Beta1,                &
                           InpNEMOHCAL%waveinput%Beta2
             READ(ufile,*) !--- Post Processing -----------------------!
             READ(ufile,*) InpNEMOHCAL%OptOUTPUT%IRF%Switch,           &
                           InpNEMOHCAL%OptOUTPUT%IRF%time_step,        &
                           InpNEMOHCAL%OptOUTPUT%IRF%duration
             READ(ufile,*) InpNEMOHCAL%OptOUTPUT%Switch_POTENTIAL
             READ(ufile,*) InpNEMOHCAL%OptOUTPUT%Kochin%Ntheta,        &
                           InpNEMOHCAL%OptOUTPUT%Kochin%min_theta,     &
                           InpNEMOHCAL%OptOUTPUT%Kochin%max_theta
             IF (InpNEMOHCAL%OptOUTPUT%Kochin%Ntheta.GT.0 ) THEN
                           InpNEMOHCAL%OptOUTPUT%Kochin%Switch=1
             ELSE
                           InpNEMOHCAL%OptOUTPUT%Kochin%Switch=0
             ENDIF

             READ(ufile,*) InpNEMOHCAL%OptOUTPUT%Freesurface%Nx,       &
                           InpNEMOHCAL%OptOUTPUT%Freesurface%Ny,       &
                           InpNEMOHCAL%OptOUTPUT%Freesurface%Lx,       &
                           InpNEMOHCAL%OptOUTPUT%Freesurface%Ly
             IF (InpNEMOHCAL%OptOUTPUT%Freesurface%Nx.GT.0 ) THEN
                           InpNEMOHCAL%OptOUTPUT%Freesurface%Switch=1
             ELSE
                           InpNEMOHCAL%OptOUTPUT%Freesurface%Switch=0
             ENDIF

             READ(ufile,*) InpNEMOHCAL%OptOUTPUT%Switch_RAO
             READ(ufile,*) InpNEMOHCAL%OptOUTPUT%FreqType

             READ(ufile,*)! ---QTF----
             READ(ufile,*)InpNEMOHCAL%OptOUTPUT%Switch_SourceDistr
             IF (InpNEMOHCAL%OptOUTPUT%Switch_SourceDistr==1) THEN
               READ(ufile,*)(InpNEMOHCAL%qtfinput%omega(K),K=1,3)
               READ(ufile,*) InpNEMOHCAL%qtfinput%bidirection
               !READ(ufile,*)InpNEMOHCAL%qtfinput%body_forward_speed
               InpNEMOHCAL%qtfinput%body_forward_speed=0 ! for now 0
               InpNEMOHCAL%qtfinput%switch_QTFP=1 ! always compute QTFP
               READ(ufile,*)InpNEMOHCAL%qtfinput%Ncontrib
               READ(ufile,*)InpNEMOHCAL%qtfinput%FSmeshfile
               READ(ufile,*)InpNEMOHCAL%qtfinput%FSRe,          &
                            InpNEMOHCAL%qtfinput%FSNRe,         &
                            InpNEMOHCAL%qtfinput%FSNBessel
               READ(ufile,*)InpNEMOHCAL%qtfinput%switch_quadHM
               READ(ufile,*)InpNEMOHCAL%qtfinput%FreqTypeOutput
               READ(ufile,*)InpNEMOHCAL%qtfinput%switch_qtfduok
               READ(ufile,*)InpNEMOHCAL%qtfinput%switch_qtfhasbo
               READ(ufile,*)InpNEMOHCAL%qtfinput%switch_qtfhasfs
             ENDIF
          CLOSE(ufile)

          IF(InpNEMOHCAL%OptOUTPUT%Switch_SourceDistr==1) THEN
            !INQUIRE (DIRECTORY=TRIM(ID%ID)//'/results/sources',           &
            !           EXIST=exist_dir) !this is Intel-specific
            INQUIRE (FILE=TRIM(ID%ID)//'/results/sources/.',           &
                       EXIST=exist_dir)
            IF (.NOT.exist_dir) CALL SYSTEM('mkdir '//TRIM(ID%ID)//       &
                                    FILESEP//'results'//FILESEP//'sources')
          END IF
         END SUBROUTINE

         FUNCTION IntegrationAXIS_FROM_MNEMOHCAL(InpNEMOHCAL)          &
                                                       result(IntegAxis)
          TYPE(TNemCal),       INTENT(IN)     :: InpNEMOHCAL
          REAL,DIMENSION(3,InpNEMOHCAL%Nintegtot)::IntegAxis
          INTEGER       :: I,J,Iinteg
          Iinteg=1
          DO I=1,InpNEMOHCAL%Nbodies
              DO J=1,InpNEMOHCAL%bodyinput(I)%NIntegration
              IntegAxis(:,Iinteg)=      &
                InpNEMOHCAL%bodyinput(I)%IntCase(J)%Axis(1:3)
                Iinteg=Iinteg+1
             END DO
          END DO

         END FUNCTION

         SUBROUTINE Discretized_Omega_and_Beta(IDQTF,waveinp,Nw,Nbeta, &
                                                w,beta)
           !input/output
           TYPE(Twaveinput),         INTENT(IN) :: waveinp
           INTEGER,                  INTENT(IN) :: IDQTF,Nw,Nbeta
           REAl,DIMENSION(Nw),       INTENT(OUT):: w
           REAl,DIMENSION(Nbeta),    INTENT(OUT):: beta
           !local variables
           INTEGER                              ::j
           REAL                                 ::wmin,wmax,dw,dwtemp
           REAL                                 ::betamin,betamax

            wmin        =waveinp%Freq1
            wmax        =waveinp%Freq2

            IF (Nw.GT.1) THEN
                dw=(wmax-wmin)/(Nw-1)
                IF (IDQTF==1) THEN
                   dwtemp = wmax/Nw
                   IF (abs(dwtemp-dw).LT.0.001) THEN
                   dw =wmin
                   END IF
           !        IF(dw.NE.wmin) THEN
           !           WRITE(*,*) &
           !  "ERROR: minimum frequency must be equal to frequency step!"
           !           STOP
           !        END IF

                END IF

                DO j=1,Nw
                    w(j)=wmin+dw*(j-1)
                END DO
            ELSE
                IF (IDQTF==1) THEN
                      WRITE(*,*) &
             "ERROR: QTF can not be run with 1 frequency only!"
                      STOP
                ENDIF
                w(1)=wmin
            END IF

            IF (waveinp%FreqType==IdFreqHz) w(:)=2*PI*w(:)
            IF (waveinp%FreqType==IdPeriod) w(:)=2*PI/w(:)

            betamin        =waveinp%Beta1
            betamax        =waveinp%Beta2
            IF (Nbeta.GT.1) THEN
                DO j=1,Nbeta
                    beta(j)=(betamin+(betamax-betamin)*(j-1)/(Nbeta-1))&
                                *PI/180.
                END DO
            ELSE
                beta(1)=betamin*PI/180.
            END IF

         END SUBROUTINE


        END MODULE
