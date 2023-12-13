#! /bin/bash

CurrPATH="$PWD"
WorkPATH=$CurrPATH

#--------------------------------------------------------------
EventName="SinglePi_30GeV"

# -----------------------------------------------------------
# Simulating / Processing event number per file
Nevt=200

# -----------------------------------------------------------
# RecoDir=LayerCut${LayerCut}_Threshold${HCAL_Threshold_MIP}MIP_Calibration${HCAL_Calibration}
# RecoPATH=$EvtDir/Reco/$GeoDir/$RecoDir
# RecoPATH=$WorkPATH/${EventName}/InitialTest_50Evt/Reco_BaselineArbor   #!!!!!!!!!!!!!!!!!!!!
SimuPATH=$WorkPATH/${EventName}/Sim
RecoPATH=$WorkPATH/${EventName}/Reco   #!!!!!!!!!!!!!!!!!!!!
jobPATH=$WorkPATH/${EventName}/job
if [ ! -d "$SimuPATH" ]
then
mkdir -p $SimuPATH
fi
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
iFile=0
nFile=10
while [ "$iFile" -lt "$nFile" ]
do
    FileID=`printf "%03d" ${iFile}`
    #--------------------------------------------------------------
    SampleName=${EventName}_${FileID}

    # GunParticle=\"gamma\",\"gamma\",\"gamma\",\"pi-\",\"pi-\"
    # GunParticleX=0.,0.,0.,0.,0.
    # GunParticleY=0.,0.,0.,0.,0.
    # GunParticleZ=0.,0.,0.,0.,0.
    # GunParticleEMin=5.,5.,5.,5.,5.
    # GunParticleEMax=5.,5.,5.,5.,5.
    # GunParticleThetaMin=80.,80.,80.,80.,80.
    # GunParticleThetaMax=100.,100.,100.,100.,100.
    # GunParticlePhiMin=-5.,-5.,-5.,-10.,-10.
    # GunParticlePhiMax=15.,15.,15.,10.,10.

    # GunParticle=\"pi-\"
    # GunParticleX=0.
    # GunParticleY=0.
    # GunParticleZ=0.
    # GunParticleEMin=30.
    # GunParticleEMax=30.
    # GunParticleThetaMin=80.
    # GunParticleThetaMax=100.
    # GunParticlePhiMin=0.
    # GunParticlePhiMax=360.

    SimuFile=$SimuPATH/Simu_${SampleName}.root
    ECALDigiFile=$RecoPATH/DigiEcal_${SampleName}.root
    HCALDigiFile=$RecoPATH/DigiHcal_${SampleName}.root
    RecoFile=$RecoPATH/Reco_${SampleName}.root
    # -----------------------------------------------------------

    simScriptFile=$jobPATH/Sim_${SampleName}.py
    recScriptFile=$jobPATH/Digi_${SampleName}.py
    shFile=$jobPATH/sub_${SampleName}.sh

    #---------------------------------
    # /bin/cp -fr $WorkPATH/template/detsim_FullDet_HCAL_temp.py ${simScriptFile}
    # sed -i "s#RNDMSEED#${iFile}#g" ${simScriptFile}
    # sed -i "s#SimPARTICLES#${GunParticle}#g" ${simScriptFile}
    # sed -i "s#SimPARTICLEX#${GunParticleX}#g" ${simScriptFile}
    # sed -i "s#SimPARTICLEY#${GunParticleY}#g" ${simScriptFile}
    # sed -i "s#SimPARTICLEZ#${GunParticleZ}#g" ${simScriptFile}
    # sed -i "s#SimPARTICLEEMIN#${GunParticleEMin}#g" ${simScriptFile}
    # sed -i "s#SimPARTICLEEMAX#${GunParticleEMax}#g" ${simScriptFile}
    # sed -i "s#SimPARTICLETHETAMIN#${GunParticleThetaMin}#g" ${simScriptFile}
    # sed -i "s#SimPARTICLETHETAMAX#${GunParticleThetaMax}#g" ${simScriptFile}
    # sed -i "s#SimPARTICLEPHIMIN#${GunParticlePhiMin}#g" ${simScriptFile}
    # sed -i "s#SimPARTICLEPHIMAX#${GunParticlePhiMax}#g" ${simScriptFile}
    # sed -i "s#SIMFILE#${SimuFile}#g" ${simScriptFile}
    # sed -i "s#NEVT#${Nevt}#g" ${simScriptFile}

    /bin/cp -fr $WorkPATH/template/digiRec_temp_truth_s1.py ${recScriptFile}
    sed -i "s#SIMFILE#${SimuFile}#g" ${recScriptFile}
    sed -i "s#ECALDIGIFILE#${ECALDigiFile}#g" ${recScriptFile}
    sed -i "s#HCALDIGIFILE#${HCALDigiFile}#g" ${recScriptFile}
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

jobBatchName=$jobPATH/sub_${EventName}_00"%{ProcId}"
echo $jobBatchName

# export PATH=/afs/ihep.ac.cn/soft/common/sysgroup/hep_job/bin:$PATH
hep_sub ${jobBatchName}.sh -g atlas -o ${jobBatchName}.out -e ${jobBatchName}.err -n $NSubJob

# hep_sub -wt test ${jobBatchName}.sh -g higgs -o ${jobBatchName}.out -e ${jobBatchName}.err -n $NSubJob
# --------------------------------------------------------------------------------
## Note: -wt default(10h)/test(5min)/short(30min)/mid(100h)/long(720h)
# --------------------------------------------------------------------------------
