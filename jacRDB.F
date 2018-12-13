        DATA Z100,Z1000/1.0D02,1.0D03/,FLAG/' ','*'/
        READ(4,12) NT,NVU,NVL,REDMAS,ER,TIT
   48   NLINE=N0
        NCALC=N0
        RESB=Z0
        I=N1
   40   READ(4,987,END=159) (HEAD(I,J),J=1,10)
  987   FORMAT(10A1)
        READ(4,*) J1,J2,NDELJ,NPRTYU
        READ(4,*) (NL(J),J=J1,J2)
        READ(4,*) (WORK1(J),J=J1,J2)
        IF(ER.EQ.Z0) READ(4,*) (ERR(J),J=J1,J2)
        IF(ER.GT.Z0) WT=Z1/ER**2
        THOU=Z0
        HUN=Z0
        DO 42 J=J1,J2
                WJ=WORK1(J)
                IF(WJ.EQ.Z0) GO TO 42
                IF(WJ.LT.Z100) GO TO 165
                IF(WJ.LT.Z1000) GO TO 162
                THOU=IDINT(WJ/Z1000)*Z1000
                HUN=IDINT((WJ-THOU)/Z100)*Z100
                GO TO 42
  162           IF(WJ.lt.1000.0.and.(WJ-HUN).GT.400.0) THOU=Z0
                HUN=IDINT(WJ/Z100)*Z100
                WORK1(J)=THOU+WJ
                GO TO 42
  165           WORK1(J)=THOU+HUN+WJ
   42           CONTINUE
