* Digitization package for rotated crystal ECAL design (by Fangyi Guo [guofangyi@ihep.ac.cn])

Data model: 
	* CaloStep: G4 step contribution to a bar, used in digitization. 
	* CaloBar: Digitized crystal bars, with member: cellID, Q1, Q2, T1, T2. 
  * CaloCluster: A set of bars after reconstruction. With member: vector<CaloBar>, phi, Z, E.

Source code: FanEcalDigiAlg
  * Merge G4 steps into digitized bars (CaloBar). Readout charge (Q1, Q2) is the sum of all G4 steps energy with an attenuation, scale factor is 1. Readout time (T1, T2) is the transport time from a G4 step to the edge of bar. Details can refer to https://indico.ihep.ac.cn/event/13888/session/3/contribution/22/material/slides/0.pdf. 
  * Bar length and hardware energy threshold can be readin from script: CrystalBarLength, EnergyThreshold. 
  * A very simple neighbor clustering algorithm for reconstruction. Cluster phi is calculated with bar cellID [crystal] in src/CaloCluster.h. May need to change if the geometry config changes. 
  * NOTE: Should separate reconstruction part and digitization part in future. 
Script: script/read.py
  * Readin file: simulated events from podioevent.inputs["RCEcalSim_Mu_Phi180.root"]
  * Import the algorithm: EcalDigi

