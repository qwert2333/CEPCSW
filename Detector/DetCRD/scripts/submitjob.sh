#!/bin/bash

Outdir=/afs/ihep.ac.cn/users/h/huangliheng/work/CEPCSW/Detector/DetCRD/scripts/jobs
if [ ! -d "$Outdir" ]; then
    mkdir -p "$Outdir"
fi

>"./jobs/scriptnames.txt"

while IFS= read -r filename; do
    scriptname="${filename%.py}.sh"
    scriptpath="${Outdir}/${scriptname}"
    echo "#!/bin/bash" >"$scriptpath"
    echo "cd /afs/ihep.ac.cn/users/h/huangliheng/work/CEPCSW/" >>"$scriptpath"
    echo "source setup.sh" >>"$scriptpath"
    echo "./run.sh /afs/ihep.ac.cn/users/h/huangliheng/work/CEPCSW/Detector/DetCRD/scripts/bashes/${filename}" >>"$scriptpath"
    echo "generate: ${scriptname}"
    chmod u+x $scriptpath
    echo "${scriptname}" >>"$Outdir/scriptnames.txt"
done <"/afs/ihep.ac.cn/users/h/huangliheng/work/CEPCSW/Detector/DetCRD/scripts/bashes/filenames.txt"
