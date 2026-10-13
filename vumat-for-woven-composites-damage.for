      SUBROUTINE VUMAT(
     1 NBLOCK, NDIR, NSHR, NSTATEV, NFIELDV, NPROPS, LANNEAL,
     2 STEPTIME, TOTALTIME, DT, CMNAME, COORDMP, CHARLENGTH,
     3 PROPS, DENSITY, STRAININC, RELSPININC,
     4 TEMPOLD, STRETCHOLD, DEFGRADOLD, FIELDOLD,
     5 STRESSOLD, STATEOLD, ENERINTERNOLD, ENERINELASOLD,
     6 TEMPNEW, STRETCHNEW, DEFGRADNEW, FIELDNEW,
     7 STRESSNEW, STATENEW, ENERINTERNNEW, ENERINELASNEW)
C
      INCLUDE 'VABA_PARAM.INC'
C
      CHARACTER*80 CMNAME
      DIMENSION COORDMP(NBLOCK,*), CHARLENGTH(NBLOCK),
     1 PROPS(NPROPS), DENSITY(NBLOCK),
     2 STRAININC(NBLOCK,NDIR+NSHR),
     3 RELSPININC(NBLOCK,NSHR), TEMPOLD(NBLOCK),
     4 STRETCHOLD(NBLOCK,NDIR+NSHR),
     5 DEFGRADOLD(NBLOCK,NDIR+NSHR+NSHR),
     6 FIELDOLD(NBLOCK,NFIELDV), STRESSOLD(NBLOCK,NDIR+NSHR),
     7 STATEOLD(NBLOCK,NSTATEV), ENERINTERNOLD(NBLOCK),
     8 ENERINELASOLD(NBLOCK), TEMPNEW(NBLOCK),
     9 STRETCHNEW(NBLOCK,NDIR+NSHR),
     1 DEFGRADNEW(NBLOCK,NDIR+NSHR+NSHR),
     2 FIELDNEW(NBLOCK,NFIELDV),
     3 STRESSNEW(NBLOCK,NDIR+NSHR), STATENEW(NBLOCK,NSTATEV),
     4 ENERINTERNNEW(NBLOCK), ENERINELASNEW(NBLOCK)
      INTEGER I, J, K
      DOUBLE PRECISION KF, WORK, DEPS(6), S0(6), S1(6)
      DOUBLE PRECISION SV0(12), SV1(12)
C
C     3D solid VUMAT. 
C     Required properties: 14.
C     PROPS: E1,E2,E3,G12,G13,G23,nu12,nu13,nu23,Xft,Xfc,Tm,
C            beta,gamma. Stress units must match strength units.
C     Abaqus VUMAT shear strain increments are tensor components.
C
      DO K=1,NBLOCK
         TEMPNEW(K)=TEMPOLD(K)
         ENERINTERNNEW(K)=ENERINTERNOLD(K)
         ENERINELASNEW(K)=ENERINELASOLD(K)
         DO J=1,NDIR+NSHR
            STRETCHNEW(K,J)=STRETCHOLD(K,J)
            STRESSNEW(K,J)=STRESSOLD(K,J)
         END DO
         DO J=1,NDIR+NSHR+NSHR
            DEFGRADNEW(K,J)=DEFGRADOLD(K,J)
         END DO
         DO J=1,NFIELDV
            FIELDNEW(K,J)=FIELDOLD(K,J)
         END DO
         DO J=1,NSTATEV
            STATENEW(K,J)=STATEOLD(K,J)
         END DO
      END DO
C
C     This implementation is restricted to 3D stress/strain.
      IF (NDIR.NE.3 .OR. NSHR.NE.3) RETURN

C
      KF=1.5D0
      IF (NPROPS.GE.15) KF=PROPS(15)
      IF (KF.LE.1.0D0) KF=1.5D0
C
      DO K=1,NBLOCK
         DO J=1,6
            DEPS(J)=STRAININC(K,J)
            S0(J)=STRESSOLD(K,J)
         END DO
         DO J=1,12
            SV0(J)=STATEOLD(K,J)
         END DO
         CALL MMF_UPDATE(PROPS,KF,DEPS,S0,SV0,S1,SV1)
         DO J=1,6
            STRESSNEW(K,J)=S1(J)
         END DO
         DO J=1,12
            STATENEW(K,J)=SV1(J)
         END DO
C        Work conjugacy for tensor shear strain components.
         WORK=0.0D0
         DO J=1,3
            WORK=WORK+0.5D0*(STRESSOLD(K,J)+STRESSNEW(K,J))
     1                 *STRAININC(K,J)
         END DO
         DO J=4,6
            WORK=WORK+(STRESSOLD(K,J)+STRESSNEW(K,J))
     1                 *STRAININC(K,J)
         END DO
         IF (DENSITY(K).GT.0.0D0) THEN
            ENERINTERNNEW(K)=ENERINTERNOLD(K)+WORK/DENSITY(K)
         END IF
      END DO
      RETURN
      END
C
C====================================================================
      SUBROUTINE MMF_UPDATE(PROPS,KF,DEPS,SOLD,SVOLD,SNEW,SVNEW)
      IMPLICIT NONE
      DOUBLE PRECISION PROPS(*), KF, DEPS(6), SOLD(6), SVOLD(*),
     1 SNEW(6), SVNEW(*)
      DOUBLE PRECISION E1,E2,E3,G12,G13,G23,V12,V13,V23
      DOUBLE PRECISION XFT,XFC,TM,BETA,GAMMA
      DOUBLE PRECISION M1Y(6,6),M1M(6,6),M2F(6,6),M2M(6,6)
      DOUBLE PRECISION EPS(6), SIGTR(6), SIGLOC(6), SY(6)
      DOUBLE PRECISION SMESO(6), SF(6), SMIC(6), SROT(6)
      DOUBLE PRECISION C(6,6), DF, DMM, DMI, DM
      DOUBLE PRECISION KFT,KFC,KMESO,KMIC,DMNEW
      DOUBLE PRECISION KM
      INTEGER I,J
C
C     SAF data are hard-coded in Fortran column-major order.

      DATA M1Y /
     1 2.160D0,0.250D0,0.150D0,0.175D0,0.0D0,0.0D0,
     2 0.270D0,1.860D0,0.120D0,0.160D0,0.0D0,0.0D0,
     3 0.150D0,0.120D0,1.450D0,0.080D0,0.0D0,0.0D0,
     4 0.102D0,0.109D0,0.850D0,0.150D0,0.0D0,0.0D0,
     5 0.0D0,0.0D0,0.0D0,0.0D0,0.850D0,0.0D0,
     6 0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.750D0 /
      DATA M1M /
     1 0.370D0,-0.27D0,-0.18D0,-0.01D0,0.0D0,0.0D0,
     2 -0.27D0,0.360D0,-0.22D0,0.0D0,0.0D0,0.0D0,
     3 -0.18D0,-0.22D0,0.420D0,0.0D0,0.0D0,0.0D0,
     4 0.0D0,0.0D0,0.0D0,0.130D0,0.0D0,0.0D0,
     5 0.0D0,0.0D0,0.0D0,0.0D0,0.150D0,0.0D0,
     6 0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.180D0 /
      DATA M2F /
     1 1.750D0,-0.200D0,-0.01D0,0.0D0,0.0D0,0.0D0,
     2 0.03D0,1.450D0,0.25D0,0.0D0,0.0D0,0.0D0,
     3 -0.01D0,0.25D0,1.350D0,0.0D0,0.0D0,0.0D0,
     4 0.0D0,0.0D0,0.0D0,-0.420D0,0.0D0,0.0D0,
     5 0.0D0,0.0D0,0.0D0,0.0D0,-0.380D0,0.0D0,
     6 0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,-0.350D0 /
      DATA M2M /
     1 0.800D0,0.05D0,0.03D0,0.0D0,0.0D0,0.0D0,
     2 0.0D0,0.670D0,0.15D0,0.0D0,0.0D0,0.0D0,
     3 0.03D0,0.15D0,0.450D0,0.0D0,0.0D0,0.0D0,
     4 0.0D0,0.0D0,0.0D0,0.200D0,0.0D0,0.0D0,
     5 0.0D0,0.0D0,0.0D0,0.0D0,0.250D0,0.0D0,
     6 0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.220D0 /
C
      E1=PROPS(1)
      E2=PROPS(2)
      E3=PROPS(3)
      G12=PROPS(4)
      G13=PROPS(5)
      G23=PROPS(6)
      V12=PROPS(7)
      V13=PROPS(8)
      V23=PROPS(9)
      XFT=PROPS(10)
      XFC=PROPS(11)
      TM=PROPS(12)
      BETA=PROPS(13)
      GAMMA=PROPS(14)
      IF (GAMMA.LE.0.0D0) GAMMA=1.5D0
C
C     State variables: 1:6 total tensor strain; 
C     9 fiber damage; 
C     12  matrix damage.
      DO I=1,12
         SVNEW(I)=SVOLD(I)
      END DO
      DO I=1,6
         EPS(I)=SVOLD(I)+DEPS(I)
      END DO
      DMM=MAX(0.0D0,MIN(0.99D0,SVOLD(7)))
      DMI=MAX(0.0D0,MIN(0.99D0,SVOLD(8)))
      DF=MAX(0.0D0,MIN(0.99D0,SVOLD(9)))
      DM=1.0D0-(1.0D0-DMM)*(1.0D0-DMI)
C
C     Trial stress from the current damaged orthotropic stiffness.
      CALL BUILD_C(E1,E2,E3,G12,G13,G23,V12,V13,V23,
     1 DF,DM,C)
      CALL MATVEC6(C,EPS,SIGTR)
C
C     Meso matrix pocket: matrix-pocket SAF in macro material axes.
      CALL MATVEC6(M1M,SIGTR,SMESO)
      CALL MATRIX_INDEX(SMESO,BETA,TM,KMESO)
C     Macro and meso axes are aligned. First obtain yarn stress in the
C     meso frame; then rotate it into warp and weft axes before M2 SAFs.
      KFT=0.0D0
      KMIC=0.0D0
      CALL MATVEC6(M1Y,SIGTR,SY)
      CALL ROTATE_Z(SY,0.0D0,SIGLOC)
      CALL MATVEC6(M2F,SIGLOC,SF)
      CALL MATVEC6(M2M,SIGLOC,SMIC)
      CALL FIBER_INDEX(SF,XFT,XFC,KFT)
      CALL MATRIX_INDEX(SMIC,BETA,TM,KM)
      KMIC=MAX(KMIC,KM)
      CALL ROTATE_Z(SY,90.0D0,SROT)
      CALL MATVEC6(M2F,SROT,SF)
      CALL MATVEC6(M2M,SROT,SMIC)
      CALL FIBER_INDEX(SF,XFT,XFC,KFC)
      KFT=MAX(KFT,KFC)
      CALL MATRIX_INDEX(SMIC,BETA,TM,KM)
      KMIC=MAX(KMIC,KM)
C
C     Brittle fiber failure and continuous normalized matrix evolution.
      IF (KFT.GE.1.0D0) DF=0.99D0
      CALL NORMAL_DAMAGE(KMESO,GAMMA,KF,DMNEW)
      DMM=MAX(DMM,DMNEW)
      CALL NORMAL_DAMAGE(KMIC,GAMMA,KF,DMNEW)
      DMI=MAX(DMI,DMNEW)
      DM=1.0D0-(1.0D0-DMM)*(1.0D0-DMI)
C
C     Damaged stiffness per manuscript: d1=d2=df, d3=dm,
C     d4=d5=d6=1-(1-df)(1-dm).
      CALL BUILD_C(E1,E2,E3,G12,G13,G23,V12,V13,V23,
     1 DF,DM,C)
      CALL MATVEC6(C,EPS,SNEW)
      DO I=1,6
         SVNEW(I)=EPS(I)
      END DO
      SVNEW(7)=DMM
      SVNEW(8)=DMI
      SVNEW(9)=DF
      SVNEW(12)=DM
      RETURN
      END
C
C====================================================================
      SUBROUTINE MATVEC6(A,X,Y)
      IMPLICIT NONE
      DOUBLE PRECISION A(6,6),X(6),Y(6)
      INTEGER I,J
      DO I=1,6
         Y(I)=0.0D0
         DO J=1,6
            Y(I)=Y(I)+A(I,J)*X(J)
         END DO
      END DO
      RETURN
      END
C
C====================================================================
      SUBROUTINE ROTATE_Z(SG,ANGLE,SL)
      IMPLICIT NONE
      DOUBLE PRECISION SG(6),SL(6),ANGLE,C,S,PI
      PI=3.141592653589793D0
      C=COS(ANGLE*PI/180.0D0)
      S=SIN(ANGLE*PI/180.0D0)
      SL(1)=C*C*SG(1)+S*S*SG(2)+2.0D0*C*S*SG(4)
      SL(2)=S*S*SG(1)+C*C*SG(2)-2.0D0*C*S*SG(4)
      SL(3)=SG(3)
      SL(4)=C*S*(SG(2)-SG(1))+(C*C-S*S)*SG(4)
      SL(5)=C*SG(5)+S*SG(6)
      SL(6)=-S*SG(5)+C*SG(6)
      RETURN
      END
C
C====================================================================
      SUBROUTINE FIBER_INDEX(SF,XFT,XFC,KFIB)
      IMPLICIT NONE
      DOUBLE PRECISION SF(6),XFT,XFC,KFIB
      KFIB=0.0D0
      IF (SF(1).GE.0.0D0 .AND. XFT.GT.0.0D0) THEN
         KFIB=SF(1)/XFT
      ELSE IF (SF(1).LT.0.0D0 .AND. XFC.GT.0.0D0) THEN
         KFIB=ABS(SF(1))/XFC
      END IF
      RETURN
      END
C
C====================================================================
      SUBROUTINE MATRIX_INDEX(S,BETA,TM,K)
      IMPLICIT NONE
      DOUBLE PRECISION S(6),BETA,TM,K,I1,VM,SEQ
      I1=S(1)+S(2)+S(3)
      VM=SQRT(0.5D0*((S(1)-S(2))**2+(S(2)-S(3))**2+
     1 (S(3)-S(1))**2+6.0D0*(S(4)**2+S(5)**2+S(6)**2)))
      SEQ=0.0D0
      IF (BETA.GT.0.0D0) THEN
         SEQ=((BETA-1.0D0)*I1+
     1   SQRT((BETA-1.0D0)**2*I1**2+4.0D0*BETA*VM**2))
     2   /(2.0D0*BETA)
      END IF
      K=0.0D0
      IF (TM.GT.0.0D0) K=MAX(0.0D0,SEQ/TM)
      RETURN
      END
C
C====================================================================
      SUBROUTINE NORMAL_DAMAGE(K,GAMMA,KF,D)
      IMPLICIT NONE
      DOUBLE PRECISION K,GAMMA,KF,D,DMAX,DEN,NUM
      DMAX=0.99D0
      D=0.0D0
      IF (K.LE.1.0D0) RETURN
      IF (K.GE.KF) THEN
         D=DMAX
         RETURN
      END IF
      DEN=1.0D0-EXP(-GAMMA*(KF-1.0D0))
      NUM=1.0D0-EXP(-GAMMA*(K-1.0D0))
      IF (ABS(DEN).GT.1.0D-14) THEN
         D=DMAX*NUM/DEN
      ELSE
         D=DMAX*(K-1.0D0)/(KF-1.0D0)
      END IF
      D=MAX(0.0D0,MIN(DMAX,D))
      RETURN
      END
C
C====================================================================
      SUBROUTINE BUILD_C(E1,E2,E3,G12,G13,G23,V12,V13,V23,
     1 DF,DM,C6)
      IMPLICIT NONE
      DOUBLE PRECISION E1,E2,E3,G12,G13,G23,V12,V13,V23
      DOUBLE PRECISION DF,DM,C6(6,6),S11,S22,S33,S12,S13,S23
      DOUBLE PRECISION NU21,NU31,NU32,DET,DS
      INTEGER I,J
      DO I=1,6
         DO J=1,6
            C6(I,J)=0.0D0
         END DO
      END DO
      IF (E1.LE.0.0D0 .OR. E2.LE.0.0D0 .OR. E3.LE.0.0D0)
     1 RETURN
      IF (G12.LE.0.0D0 .OR. G13.LE.0.0D0 .OR. G23.LE.0.0D0)
     1 RETURN
      NU21=V12*E2/E1
      NU31=V13*E3/E1
      NU32=V23*E3/E2
C     Damaged compliance follows the d1=d2=df, d3=dm rule.
      S11=1.0D0/((1.0D0-DF)*E1)
      S22=1.0D0/((1.0D0-DF)*E2)
      S33=1.0D0/((1.0D0-DM)*E3)
      S12=-V12/E1
      S13=-V13/E1
      S23=-V23/E2
      DET=S11*S22*S33+2.0D0*S12*S13*S23
     1    -S11*S23*S23-S22*S13*S13-S33*S12*S12
      IF (ABS(DET).LT.1.0D-30) RETURN
      C6(1,1)=(S22*S33-S23*S23)/DET
      C6(2,2)=(S11*S33-S13*S13)/DET
      C6(3,3)=(S11*S22-S12*S12)/DET
      C6(1,2)=(S13*S23-S12*S33)/DET
      C6(2,1)=C6(1,2)
      C6(1,3)=(S12*S23-S13*S22)/DET
      C6(3,1)=C6(1,3)
      C6(2,3)=(S12*S13-S11*S23)/DET
      C6(3,2)=C6(2,3)
C     VUMAT stores tensor shear strain; use 2G in stress update.
      DS=1.0D0-(1.0D0-DF)*(1.0D0-DM)
      C6(4,4)=2.0D0*G12*(1.0D0-DS)
      C6(5,5)=2.0D0*G13*(1.0D0-DS)
      C6(6,6)=2.0D0*G23*(1.0D0-DS)
      RETURN
      END
