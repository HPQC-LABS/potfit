c===========================================================================
      SUBROUTINE quadCORR(NCMM,MCMM,NCMMAX,MMLR,De,CmVAL,CmEFF)
c===========================================================================
c** subroutine to generate and print MLR CmEFF values incorporating
c  quadratic 'Dattani' corrections to Cm values for both standard 'linear'
c  and A-F diagonalized uLR(r) functions for MLR potentials
c** Return MCMM= NCMM+1  for C9{adj} term for m_1= 3 potentials
c===========================================================================
      INTEGER NCMM,MCMM,NCMMAX,MMLR(NCMMAX)
      REAL*8 De,CmVAL(NCMMAX),CmEFF(NCMMAX)
c----------------------------------------------------------------------
      IF(MMLR(1).GT.0) THEN
c** For 'normal' inverse-power sum MLR case, with or without damping,
c   set up Dattani's 'Quadratic-corrected' effective Cm values 
          IF((MMLR(1).EQ.6).AND.(NCMM.GE.4)) THEN
c... First, consider C6/C12adj(C14adj) for MMLR(m)={6,8,10,(11),12,14} case
              IF(MMLR(4).EQ.12) THEN             ! explicitly MMLR(4)=12
                  CmEFF(4)= CmVAL(4)+ 0.25D0*CmVAL(1)**2/De
                  WRITE(6,710) MMLR(4),MMLR(4),CmEFF(4)
                  ENDIF
              IF(NCMM.GE.5) THEN
                 IF(MMLR(4).EQ.11) THEN         ! implicitly MMLR(5)=12
                     CmEFF(5)= CmVAL(5) + 0.25D0*CmVAL(1)**2/De
                     WRITE(6,710) MMLR(5),MMLR(5),CmEFF(5)
                     IF(NCMM.GE.6) THEN         ! implicitly MMLR(6)=14
                         CmEFF(6)= CmVAL(6)+ 0.5D0*CmVAL(1)*CmVAL(2)/De
                         WRITE(6,710) MMLR(6),MMLR(6),CmEFF(6)
                         ENDIF
                     ENDIF
                 IF(MMLR(4).EQ.12) THEN           ! assuming MMLR(5)=14
                     CmEFF(5)= CmVAL(5) + 0.5D0*CmVAL(1)*CmVAL(2)/De
                     WRITE(6,710) MMLR(5),MMLR(5),CmEFF(5)
                     ENDIF
                 ENDIF
              ENDIF
          IF((MMLR(1).EQ.5).AND.(NCMM.GE.4)) THEN
c... Then, consider C5/C10adj + C12adj for MMLR(m)={5,6,8,10,12,14} cases
              CmEFF(4)= CmVAL(4) + 0.25D0*CmVAL(1)**2/De
              WRITE(6,710) MMLR(4),MMLR(4),CmEFF(4)
              IF(NCMM.GE.5) THEN                 ! introduce C12^{adj}
                  CmEFF(5)= CmVAL(5) + 0.25D0*CmVAL(2)**2/De
                  WRITE(6,710) MMLR(5),MMLR(5),CmEFF(5)
                  IF(NCMM.GE.6) THEN             ! introduce C14^{adj}
                      CmEFF(6)= CmVAL(6) + 0.5D0*CmVAL(2)*CmVAL(3)/De
                      WRITE(6,710) MMLR(6),MMLR(6),CmEFF(6)
                      ENDIF
                  ENDIF
              ENDIF
          IF((MMLR(1).EQ.4).AND.(NCMM.GE.3).and.(MMLR(3).EQ.8)) THEN
c... Then, consider C4/C8adj + C12adj for MMLR(m)={4,6,8,10,12,14} cases
              CmEFF(3)= CmVAL(3) + 0.25D0*CmVAL(1)**2/De
              WRITE(6,712) MMLR(3),MMLR(3),CmEFF(3)
              IF(NCMM.GE.4) THEN                 ! implicitly MMLR(4)=10
                  CmEFF(4)= CmVAL(4) + 0.5D0*CmVAL(1)*CmVAL(2)/De
                  WRITE(6,710) MMLR(4),MMLR(4),CmEFF(4)
                  IF(NCMM.GE.5) THEN             ! implicitly MMLR(5)=12
                      CmEFF(5)= CmVAL(5) + 0.5D0*CmVAL(1)*CmVAL(3)/De
     1                                       + 0.25D0*CmVAL(2)**2/De
                      WRITE(6,710) MMLR(5),MMLR(5),CmEFF(5)
                      IF(NCMM.GE.6) THEN         ! implicitly MMLR(6)=14
                          CmEFF(6)= CmVAL(6)+ 0.5D0*CmVAL(2)*CmVAL(3)/De
     1                                      + 0.5D0*CmVAL(1)*CmVAL(4)/De
                          WRITE(6,710) MMLR(6),MMLR(6),CmEFF(6)
                          ENDIF
                      ENDIF
                  ENDIF
                ENDIF                   !! consider no further adjustment
          IF((MMLR(1).EQ.4).AND.(NCMM.GE.4).and.(MMLR(4).EQ.8)) THEN
c... Then, consider C4/C8adj + C12adj for MMLR(m)={4,6,7,8} cases
                  CmEFF(4)= CmVAL(4) + 0.25D0*CmVAL(1)**2/De
                  WRITE(6,712) MMLR(4),MMLR(4),CmEFF(4)
                  ENDIF
          IF((MMLR(1).EQ.3).AND.(NCMM.GE.2)) THEN
c... Then, consider C3/C6adj + C9adj for MMLR(m)={3,6,8,(9),10,12,14} cases
              CmEFF(2)= CmVAL(2) + 0.25D0*CmVAL(1)**2/De 
              WRITE(6,712) MMLR(2),MMLR(2),CmEFF(2)
              IF(NCMM.GE.3) THEN              ! introduce C9adj & MMLR=9
                  MCMM= NCMM + 1
                  MMLR(MCMM)= 9 
                  CmEFF(MCMM)= 0.5d0*CmVAL(1)*CmEFF(2)/De
                  WRITE(6,714) MMLR(MCMM),CmEFF(MCMM)
                  IF(NCMM.GE.5) THEN             ! implicitly MMLR(5)=12
                      CmEFF(5)= CmVAL(5) + 0.5D0*CmVAL(1)*CmEFF(MCMM)/De
     1                                         + 0.25D0*CmEFF(2)**2/De
                      WRITE(6,710) MMLR(5),MMLR(5),CmEFF(5)
                      IF(NCMM.GE.6) THEN         ! implicitly MMLR(6)=14
                          CmEFF(6)= CmVAL(6)+ 0.5D0*CmEFF(2)*CmVAL(3)/De
                          WRITE(6,710) MMLR(6),MMLR(6),CmEFF(6)
                          ENDIF
                      ENDIF
                  ENDIF
              ENDIF
          ENDIF
c======================================================================= c
c** End of  CmEFF= Cm + CmADJ  setup for non-AF case ===================
  710 Format("  'Quadratic correction' for   C",I2,'(MLR)   yields',
     1  6x,'C',I2,'{adj}=',1PD15.8)
  712 Format("  'Quadratic correction' for   C",I1,'(MLR)    yields',
     1  7x,'C'I1,'{adj}=',1PD15.8)
  714 Format("  'Quadratic corrn' for  MLR(m_1=3)  introduces   C",
     1    I1,'(',A4,',adj) =',1PD15.8)
  716 Format("  'Quadratic correction' for  C",I1,'(Sigma)  yields   C',
     1    I1,'(Sigma,adj)=',1PD15.8)
  718 Format("  'Quadratic correction' for  C",I1,'(^3Pi)   yields   C',
     1    I1,'(^3Pi,adj) =',1PD15.8)
  720 Format("  'Quadratic correction' for  C",I1,'(^1Pi)   yields  C',
     1    I1,'(^1Pi,adj) =',1PD15.8)
c=========================================================================      
      IF(MMLR(1).LE.0) THEN
c** implement Quadratic 'Dattani' MLR corrections for AF cases         
          IF(MMLR(1).GE.-1) THEN         !! first for the 2x2 cases ...
              CmEFF(4)= CmVAL(4) + 0.25*CmVAL(2)**2/De
              CmEFF(5)= CmVAL(5) + 0.25*CmVAL(3)**2/De
              WRITE(6,716) MMLR(4),MMLR(4),CmEFF(4)
              WRITE(6,718) MMLR(5),MMLR(5),CmEFF(5)
c*  prepare C9{adj} coefficients for addition to chosen root
              MMLR(8)= 9               !! These terms added just
              MMLR(9)= 9               !! before exit from  AFdiag
              Cmeff(8)= 0.5*CmVAL(2)*CmEFF(4)/De   
              WRITE(6,714) MMLR(8),'Sigm',CmEFF(8)
              Cmeff(9)= 0.5*CmVAL(3)*CmEFF(5)/De
              WRITE(6,714) MMLR(9),'^3Pi',CmEFF(9)
              ENDIF
          IF(MMLR(1).LE.-2) THEN         !! now for the 3x3 cases ...
              CmEFF(5)= CmVAL(5) + 0.25*CmVAL(2)**2/De
              WRITE(6,716) MMLR(5),MMLR(5),CmEFF(5)
              CmEFF(6)= CmVAL(6) + 0.25*CmVAL(3)**2/De
              WRITE(6,720) MMLR(6),MMLR(6),CmEFF(6)
              CmEFF(7)= CmVAL(7) + 0.25*CmVAL(4)**2/De
              WRITE(6,718) MMLR(7),MMLR(7),CmEFF(7)
c*  prepare C9{adj} coefficients for addition to chosen root
              MMLR(11)= 9               !! These terms added just
              MMLR(12)= 9               !! before exit from  AFdiag
              MMLR(13)= 9
              Cmeff(11)= 0.5*CmVAL(2)*CmEFF(5)/De   
              IF(MMLR(1).EQ.-2) WRITE(6,714) MMLR(11),'Sigm',CmEFF(11)
              Cmeff(12)= 0.5*CmVAL(3)*CmEFF(6)/De
              IF(MMLR(1).EQ.-3) WRITE(6,714) MMLR(12),'^3Pi',CmEFF(12)
              Cmeff(13)= 0.5*CmVAL(4)*CmEFF(7)/De
              IF(MMLR(1).EQ.-4) WRITE(6,714) MMLR(13),'^1Pi',CmEFF(13)
              ENDIF
          ENDIF
      RETURN
      END
c23456789012345678901234567890123456789012345678901234567890123456789012
