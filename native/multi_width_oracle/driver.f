C Diagnostic O-only driver of exact pinned MULTI DAMP. No RT/EOS solve.
      PROGRAM WIDTH_PROBE
      INCLUDE 'PREC'
      INCLUDE 'PARAM'
      INCLUDE 'CATOM'
      INCLUDE 'CATMOS'
      INCLUDE 'CATMO2'
      INCLUDE 'CCONST'
      INCLUDE 'CINPUT'
      INCLUDE 'CLU'
      INCLUDE 'CTERM'
      REAL(8) GAMTOT(MDEP,MLINE),GAMRAD(MDEP,MLINE),
     *        GAMVDW(MDEP,MLINE),GAMSTK(MDEP,MLINE)
      COMMON/WDREC/ GAMTOT,GAMRAD,GAMVDW,GAMSTK
      READ(*,*) NK,NLINE,NDEP,AMASS
      IF(NK.LT.1.OR.NK.GT.MK.OR.NLINE.LT.1.OR.
     * NLINE.GT.MLINE.OR.NDEP.LT.1.OR.NDEP.GT.MDEP.OR.
     * AMASS.LE.0) STOP 2
      ATOMID='O'
      LJOBLO=6
      IWDAMP=1
      QNORM=10.0D0
      HCE=HH*CC/EE
      AWGT=AMASS*UU
      DO I=1,NK
        READ(*,*) E_CM,G(I),ION(I)
        EV(I)=E_CM*HCE
      ENDDO
      DO KR=1,NLINE
        READ(*,*) IRAD(KR),JRAD(KR),GA(KR),GW(KR),GQ(KR),
     *            ALAMB(KR)
        NTERM(KR)=0
      ENDDO
      DO K=1,NDEP
        READ(*,*) TEMP(K),NE(K),HTOTAL,VMIC
        IF(TEMP(K).LE.0.OR.NE(K).LE.0.OR.HTOTAL.LE.0.OR.
     *     VMIC.LT.0) STOP 3
C NH(1) is aggregate storage for DAMP's SUM NH(1:5) ONLY.
C This is not a hydrogen level-population model. ATOMID is hard O.
        DO J=1,6
          NH(J,K)=0
        ENDDO
        NH(1,K)=HTOTAL
        VTURB(K)=VMIC*1.0D-3/QNORM
      ENDDO
      CALL DAMP
      WRITE(*,'(A,6(ES25.16E3,1X))') 'CONSTANTS ',
     * BK,UU,HH,CC,EE,PI
      DO KR=1,NLINE
        DO K=1,NDEP
          WRITE(*,'(2(I5,1X),5(ES25.16E3,1X))')
     * KR,K,GAMTOT(K,KR),GAMRAD(K,KR),GAMVDW(K,KR),
     * GAMSTK(K,KR),ADAMP(K,KR)
        ENDDO
      ENDDO
      END

      SUBROUTINE WDAMP(K,KR,DDP,GT,GR,GV,GS)
      INCLUDE 'PREC'
      INCLUDE 'PARAM'
      REAL(8) GAMTOT(MDEP,MLINE),GAMRAD(MDEP,MLINE),
     *        GAMVDW(MDEP,MLINE),GAMSTK(MDEP,MLINE)
      COMMON/WDREC/ GAMTOT,GAMRAD,GAMVDW,GAMSTK
      IF(K.EQ.0.AND.KR.EQ.0) RETURN
      IF(K.LT.1.OR.K.GT.MDEP.OR.KR.LT.1.OR.KR.GT.MLINE)
     * STOP 4
      GAMTOT(K,KR)=GT
      GAMRAD(K,KR)=GR
      GAMVDW(K,KR)=GV
      GAMSTK(K,KR)=GS
      END

      SUBROUTINE STOP(TEXT)
      CHARACTER*(*) TEXT
      WRITE(*,*) TEXT
      STOP 5
      END
