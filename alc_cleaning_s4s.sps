

RECODE Y1F_alc_7 (1=0) (2=1) (3=2) (4=3) (5=4) (6=5) (7=6) (8=7) (9=8) (10=9) (11=10) (12=11) 
    (13=12) (14=13) (15=14) (16=15) (17=16) (18=17) (19=18) (20=19) (21=20) (22=21) (23=22) (24=23) 
    (25=24) (26=25) (27=26) (28=27) (29=28) (30=29) (31=30) (-9=SYSMIS) (-99=SYSMIS) (-999=SYSMIS) INTO past30days.
EXECUTE.

FREQUENCIES past30days.

RECODE Y1F_alc_5 (1=0) (2=1) (3=4) (4=12) (5=16) (-9=SYSMIS) (-99=SYSMIS) (-999=SYSMIS) INTO falloften.
EXECUTE.

FREQUENCIES falloften.

RECODE Y1F_alc_6 (1=2) (2=4) (3=6) (4=8) (5=10) (-9=SYSMIS) (-99=SYSMIS) (-999=SYSMIS) INTO fallmany.
EXECUTE.

FREQUENCIES fallmany.

RECODE past30days (0 thru 1=1) (2 thru 4=2) (4 thru 12=3) (13 thru 16=4) (17 thru Highest=5) (-9=SYSMIS) (-99=SYSMIS) (-999=SYSMIS) INTO 
    often2011.
EXECUTE.

RECODE Y1F_alc_7 (0=0) (1 thru 2=2) (3 thru 4=4) (5 thru 6=6) (7 thru 9=8) (10 thru Highest=10) (-9=SYSMIS) (-99=SYSMIS) (-999=SYSMIS)
    INTO many2011.
EXECUTE.

RECODE Y1S_alc_5 (1=0) (2=1) (3=4) (4=12) (5=16) (-9=SYSMIS) (-99=SYSMIS) (-999=SYSMIS) INTO springoften.
EXECUTE.

RECODE Y1S_alc_6 (1=2) (2=4) (3=6) (4=8) (5=10) (-9=SYSMIS) (-99=SYSMIS) (-999=SYSMIS) INTO springmany.
EXECUTE.

DO IF (Y1F_alc_1=2).
RECODE falloften fallmany often2011 many2011 (SYSMIS=0).
END IF.
EXECUTE.

DO IF (Y1S_alc_1 = 2).
RECODE springoften springmany (SYSMIS=0).
END IF.
EXECUTE.

RECODE falloften fallmany springoften springmany often2011 many2011 (SYSMIS=-99).
EXECUTE.

COMPUTE spring_drinking=springoften * springmany.
EXECUTE.

COMPUTE fall_drinking=falloften * fallmany.
EXECUTE.

FREQUENCIES VARIABLES=spring_drinking fall_drinking
  /STATISTICS=STDDEV MEAN
  /HISTOGRAM NORMAL
  /ORDER=ANALYSIS.

COMPUTE delta_drinking=spring_drinking - fall_drinking.
EXECUTE.

FREQUENCIES VARIABLES=delta_drinking
  /HISTOGRAM NORMAL
  /ORDER=ANALYSIS.

COMMENT transitioning into roommate by roommate data

SORT CASES BY Spring_Pair_ID .
CASESTOVARS
  /ID=Spring_Pair_ID
  /GROUPBY=VARIABLE.



