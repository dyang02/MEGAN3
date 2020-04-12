#! /bin/csh -fx
########################################################################
## Common setups
#source setcase.csh

## Directory setups
setenv PRJ USA 
setenv PROMPTFLAG N

# Program directory
setenv PROG megan_test_v2
#setenv PROG megan_test
setenv EXE    ./$PROG

# Input directory
setenv METDIR ../inputs/
setenv INPDIR ../inputs/

# Output directory
setenv OUTDIR ../outputs/

# Log directory
setenv LOGDIR ./LOGS
mkdir -p $LOGDIR
########################################################################

# Grid definition
foreach dom ( 12US1 ) #36 12 )
setenv GRIDDESC $cwd/GRIDDESC
setenv GDNAM3D  ${dom} 
set scen = "J3"

setenv SDATE 2016002        #Episode start date
setenv STIME 0		    #start time
#setenv RLENG 240000	    # time step of meteorology files
set JDATE = $SDATE

########################################################################
# CANTYP
setenv CANTYP $INPDIR/CT3_$GDNAM3D.ncf

# LAIS46
setenv LAIS46 $INPDIR/LAI3_$GDNAM3D.ncf

# MGNMET
setenv  MGNMET   $METDIR/MET.MEGAN.$GDNAM3D.${JDATE}.ncf

# LDFILE
setenv LDFILE $INPDIR/LDF_$GDNAM3D.$scen.ncf

########################################################################
# Output, dailyMET
setenv DailyMET $OUTDIR/DAYMET.$GDNAM3D.$JDATE.test.ncf

# Output, CANMET
setenv CANMET $OUTDIR/CANMET.$GDNAM3D.${SDATE}.test.ncf

# Output, MGNSEA
setenv MGNSEA $OUTDIR/MGNSEA.$GDNAM3D.${SDATE}.test.ncf

# Output, MEGVEA
setenv MGNERS $OUTDIR/MGNERS.$GDNAM3D.$scen.${SDATE}.nostress.test.ncf

# LOG file
setenv LOGFILE  $LOGDIR/run.$PROG.$GDNAM3D.$SDATE.log

########################################################################

rm -f $LOGFILE
if ( ! -e $LOGDIR ) mkdir -p $LOGDIR
#$EXE | tee $LOGDIR/log.run.$PROG.$GDNAM3D.txt
$EXE

end # dom
