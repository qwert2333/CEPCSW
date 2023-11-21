#! /bin/bash

CurrPATH="$PWD"
WorkPATH=$CurrPATH
StdHepPATH="/cefs/data/stdhep/CEPC240/higgs/Higgs_10M/data/E240.Pnnh_gg.e0.p0.whizard195"
#--------------------------------------------------------------
EventName="E240_nnHgg"

# -----------------------------------------------------------
# Simulating / Processing event number per file
Nevt=100

# -----------------------------------------------------------
# RecoDir=LayerCut${LayerCut}_Threshold${HCAL_Threshold_MIP}MIP_Calibration${HCAL_Calibration}
# RecoPATH=$EvtDir/Reco/$GeoDir/$RecoDir
# RecoPATH=$WorkPATH/${EventName}/InitialTest_50Evt/Reco_BaselineArbor   #!!!!!!!!!!!!!!!!!!!!
#SimuPATH=$WorkPATH/${EventName}/Sim
SimuPATH=/cefs/higgs/guofy/CEPCSW_v212/run/E240_nnHgg/Sim
RecoPATH=$WorkPATH/${EventName}_Reco/Reco   #!!!!!!!!!!!!!!!!!!!!
jobPATH=$WorkPATH/${EventName}_Reco/job
#if [ ! -d "$SimuPATH" ]
#then
#mkdir -p $SimuPATH
#fi
if [ ! -d "$RecoPATH" ]
then
mkdir -p $RecoPATH
fi
if [ ! -d "$jobPATH" ]
then
mkdir -p $jobPATH
fi
#------------------

#--------------------------------------------------------------
iFile=1
nFile=50
while [ "$iFile" -le "$nFile" ]
do
    FileID=`printf "%05d" ${iFile}`
    #--------------------------------------------------------------
    SampleName=${EventName}_${iFile}
    StdHepFile=$StdHepPATH/nnh_gg.e0.p0.${FileID}.stdhep

    SimuFile=$SimuPATH/Simu_${SampleName}.root
    DigiFile=$RecoPATH/Digi_${SampleName}.root
    RecoFile=$RecoPATH/Reco_${SampleName}.root
    # -----------------------------------------------------------

    simScriptFile=$jobPATH/Sim_${SampleName}.py
    recScriptFile=$jobPATH/Digi_${SampleName}.py
    shFile=$jobPATH/sub_${SampleName}.sh

    #---------------------------------
#    /bin/cp -fr $WorkPATH/template/detsim_nnHgg_temp.py ${simScriptFile}
#    sed -i "s#RNDMSEED#${iFile}#g" ${simScriptFile}
#    sed -i "s#READSTDHEP#${StdHepFile}#g" ${simScriptFile}
#    sed -i "s#SIMFILE#${SimuFile}#g" ${simScriptFile}
#    sed -i "s#NEVT#${Nevt}#g" ${simScriptFile}

    /bin/cp -fr $WorkPATH/template/digiRec_temp.py ${recScriptFile}
    sed -i "s#SIMFILE#${SimuFile}#g" ${recScriptFile}
    sed -i "s#DIGIFILE#${DigiFile}#g" ${recScriptFile}
    sed -i "s#RECFILE#${RecoFile}#g" ${recScriptFile}
    sed -i "s#NEVT#${Nevt}#g" ${recScriptFile}

    #---------------------------------
    echo \
    "
cd $jobPATH
#./run.sh ${simScriptFile}
./run.sh ${recScriptFile}

    " > ${shFile}

    chmod +x  ${shFile}
    #-------------------- run *.sh -----------------------------
    # . ${shFile}
    #-------------------- or hep_sub single *.sh -----------------------------
    # hep_sub ${shFile} -g higgs -o ${outFile} -e ${errFile}
    # -------------------------

    let "iFile+=1"
done
/bin/cp $WorkPATH/run.sh $jobPATH

# --- hep_sub in batch mode ------------------------------------------------------>
NSubJob=$iFile
echo "Total $NSubJob files will be hep_sub."

jobBatchName=$jobPATH/sub_${EventName}_"%{ProcId}"
echo $jobBatchName

# export PATH=/afs/ihep.ac.cn/soft/common/sysgroup/hep_job/bin:$PATH
hep_sub ${jobBatchName}.sh -g atlas -o ${jobBatchName}.out -e ${jobBatchName}.err -n $NSubJob

# hep_sub -wt test ${jobBatchName}.sh -g higgs -o ${jobBatchName}.out -e ${jobBatchName}.err -n $NSubJob
# --------------------------------------------------------------------------------
## Note: -wt default(10h)/test(5min)/short(30min)/mid(100h)/long(720h)
# --------------------------------------------------------------------------------
