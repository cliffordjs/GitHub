MIXED
  spring_drinking_A BY fall_drinking_A fall_drinking_P spring_drinking_P with Y1F_dem_2_A
  /RANDOM = fall_drinking_A fall_drinking_P spring_drinking_P
  /PRINT = SOLUTION TESTCOV
  /REPEATED = partnum | SUBJECT(Fall_PairID) COVTYPE(CS) .


MIXED spring_drinking_A BY fall_hall_drinks spring_hall_drinks spring_drinking_P fall_drinking_P fall_drinking_A randomornot_A fall_HallID
    WITH Y1F_dem_2_A
  /CRITERIA=CIN(95) MXITER(100) MXSTEP(10) SCORING(1) SINGULAR(0.000000000001) HCONVERGE(0, 
    ABSOLUTE) LCONVERGE(0, ABSOLUTE) PCONVERGE(0.000001, ABSOLUTE)
  /METHOD=REML
  /PRINT=SOLUTION TESTCOV
  /RANDOM=fall_drinking_P(fall_HallID) spring_drinking_P(fall_HallID) fall_hall_drinks fall_drinking_A fall_drinking_A (fall_HallID)spring_hall_drinks spring_drinking_P 
    fall_drinking_P randomornot_A | COVTYPE(VC).


MIXED spring_drinking_A BY randomornot_A WITH fall_drinking_A spring_drinking_P fall_drinking_P
  /CRITERIA=CIN(95) MXITER(100) MXSTEP(10) SCORING(1) SINGULAR(0.000000000001) HCONVERGE(0, 
    ABSOLUTE) LCONVERGE(0, ABSOLUTE) PCONVERGE(0.000001, ABSOLUTE)
  /FIXED=| SSTYPE(3)
  /METHOD=REML
  /PRINT=SOLUTION
  /RANDOM=randomornot_A fall_drinking_A spring_drinking_P fall_drinking_P 
    randomornot_A*fall_drinking_A randomornot_A*spring_drinking_P randomornot_A*fall_drinking_P 
    fall_drinking_A*spring_drinking_P fall_drinking_A*fall_drinking_P spring_drinking_P*fall_drinking_P 
    randomornot_A*fall_drinking_A*spring_drinking_P randomornot_A*fall_drinking_A*fall_drinking_P 
    randomornot_A*spring_drinking_P*fall_drinking_P fall_drinking_A*spring_drinking_P*fall_drinking_P 
    randomornot_A*fall_drinking_A*spring_drinking_P*fall_drinking_P | SUBJECT(fall_PairID) COVTYPE(VC)
  /REPEATED=partnum | SUBJECT(fall_PairID) COVTYPE(CS).


MIXED spring_drinking_A BY randomornot_A WITH fall_drinking_A spring_drinking_P fall_drinking_P
  /CRITERIA=CIN(95) MXITER(100) MXSTEP(10) SCORING(1) SINGULAR(0.000000000001) HCONVERGE(0, 
    ABSOLUTE) LCONVERGE(0, ABSOLUTE) PCONVERGE(0.000001, ABSOLUTE)
  /FIXED=| SSTYPE(3)
  /METHOD=REML
  /PRINT=SOLUTION TESTCOV
  /RANDOM=randomornot_A fall_drinking_A spring_drinking_P fall_drinking_P 
    randomornot_A*fall_drinking_A randomornot_A*spring_drinking_P randomornot_A*fall_drinking_P 
    fall_drinking_A*spring_drinking_P fall_drinking_A*fall_drinking_P spring_drinking_P*fall_drinking_P 
    randomornot_A*fall_drinking_A*spring_drinking_P randomornot_A*fall_drinking_A*fall_drinking_P 
    randomornot_A*spring_drinking_P*fall_drinking_P fall_drinking_A*spring_drinking_P*fall_drinking_P 
    randomornot_A*fall_drinking_A*spring_drinking_P*fall_drinking_P INTERCEPT | SUBJECT(fall_HallID) 
    COVTYPE(VC)
  /RANDOM=randomornot_A fall_drinking_A spring_drinking_P fall_drinking_P 
    randomornot_A*fall_drinking_A randomornot_A*spring_drinking_P randomornot_A*fall_drinking_P 
    fall_drinking_A*spring_drinking_P fall_drinking_A*fall_drinking_P spring_drinking_P*fall_drinking_P 
    randomornot_A*fall_drinking_A*spring_drinking_P randomornot_A*fall_drinking_A*fall_drinking_P 
    randomornot_A*spring_drinking_P*fall_drinking_P fall_drinking_A*spring_drinking_P*fall_drinking_P 
    randomornot_A*fall_drinking_A*spring_drinking_P*fall_drinking_P INTERCEPT  | SUBJECT(fall_HallID*fall_PairID).



MIXED spring_drinking_A WITH fall_drinking_A fall_drinking_P Y1F_dem_2_A fall_HallID randomornot_A
   /FIXED = fall_drinking_A fall_drinking_P Y1F_dem_2_A fall_HallID randomornot_A
   /PRINT = SOLUTION TESTCOV
   /REPEATED = partnum | SUBJECT(fall_HallID) COVTYPE(CS)
   /REPEATED = partnum | SUBJECT(fall_PairID*fall_HallID) COVTYPE(CS).
