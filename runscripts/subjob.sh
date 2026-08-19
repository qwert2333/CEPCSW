#!/bin/bash

RUNDIR=/afs/ihep.ac.cn/users/h/huangliheng/work/CEPCSW/runscripts/jobs
RUN=read.sh
RUN1=Runrc.sh

echo "/afs/ihep.ac.cn/soft/common/sysgroup/hep_job/bin/hep_sub -g cms $RUN -argu \"%{ProcId}\" -n 10 -e /afs/ihep.ac.cn/users/h/huangliheng/pub/err -o /afs/ihep.ac.cn/users/h/huangliheng/pub/out" >$RUNDIR/$RUN1

echo 'procid=$1
scriptname=$(sed -n "$((procid + 1))p" /afs/ihep.ac.cn/users/h/huangliheng/work/CEPCSW/runscripts/jobs/scriptnames.txt)
bash /afs/ihep.ac.cn/users/h/huangliheng/work/CEPCSW/runscripts/jobs/$scriptname' >$RUNDIR/$RUN

chmod u+x $RUNDIR/$RUN
chmod u+x $RUNDIR/$RUN1

cd $RUNDIR
./$RUN1
