#ifndef CalibECALPrototypeAlg_C
#define CalibECALPrototypeAlg_C

#include "CalibECALPrototypeAlg.h"

#include "edm4hep/Vector3f.h"

#include "DD4hep/Detector.h"
#include <DD4hep/Objects.h>
#include "DD4hep/DD4hepUnits.h"
#include "DDRec/Vector3D.h"

#include "GaudiKernel/INTupleSvc.h"
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/IRndmGen.h"
#include "GaudiKernel/IRndmGenSvc.h"
#include "GaudiKernel/RndmGenerators.h"

#include "edm4hep/SimCalorimeterHit.h"
#include "edm4hep/CalorimeterHit.h"
#include "edm4hep/Vector3f.h"

#include "TVector2.h"
#include "TFitResult.h"
#include <cmath>
#include <iostream>
#include <algorithm>
#include <map>
#include <random>
#include <vector>

using namespace std ;

DECLARE_COMPONENT( CalibECALPrototypeAlg )

CalibECALPrototypeAlg::CalibECALPrototypeAlg(const std::string& name, ISvcLocator* svcLoc)
: GaudiAlgorithm(name, svcLoc)
{
  // Input collections
  declareProperty("RawHGHitsCollection", m_inputRawHGHitHdl, "Handle of the input hit collection");
  declareProperty("RawLGHitsCollection", m_inputRawLGHitHdl, "Handle of the input hit collection");
  // Output collections
  declareProperty("CalibHitCollection", m_outputHitHdl, "Handle of the output hit collection");
}


StatusCode CalibECALPrototypeAlg::initialize()
{
  // set rand seed;
  m_randSvc = service<IRndmGenSvc>("RndmGenSvc");

  // Ntuple 
  if(_writeNtuple)
  {
    std::string s_outfile = _filename.value();
    //m_wfile = new TFile(s_outfile.c_str(), "recreate");
    //m_tree = new TTree("tree", "muon digi tree");


  }

	//Get calibration constant from input MetaData
	if(!m_metadata_path.value().empty()){
		podio::ROOTFrameReader m_reader;
		m_reader.openFile(m_metadata_path.value());
		auto frame = podio::Frame(m_reader.readEntry(podio::Category::Metadata, 0));
    m_metadata = std::make_unique<podio::GenericParameters>(frame.getParameters());
		m_cellID = m_metadata->getValue<std::vector<int>>("CellIDVector");
    info() << "Read calibration constants from metadata "<< endmsg;

		if(!m_dopedestal){ 
			m_pedestalHG = m_metadata->getValue<std::vector<float>>("PedestalHGVector");
			m_pedestalLG = m_metadata->getValue<std::vector<float>>("PedestalLGVector");
			m_pedestalHGsigma = m_metadata->getValue<std::vector<float>>("PedestalHGSigmaVector");
			m_pedestalLGsigma = m_metadata->getValue<std::vector<float>>("PedestalLGSigmaVector");

			if(m_cellID.size()!=m_pedestalHG.size() || m_cellID.size()!=m_pedestalLG.size() ){
				error() << "In MetaData: cellID size ("<<m_cellID.size()<<") != pedestal size ("<<m_pedestalHG.size()<<" / "<<m_pedestalLG.size()<<") "<<endmsg;	
			}
			else{
				for(int i=0; i<m_cellID.size(); i++){
					map_pedestalHG[m_cellID[i]] = m_pedestalHG[i];
					map_pedestalLG[m_cellID[i]] = m_pedestalLG[i];
					map_pedestalHGsigma[m_cellID[i]] = m_pedestalHGsigma[i];
					map_pedestalLGsigma[m_cellID[i]] = m_pedestalLGsigma[i];
				}
			}
		}

    if(!m_domip){
      m_MIPHG = m_metadata->getValue<std::vector<float>>("MIPHGVector");
      m_MIPLG = m_metadata->getValue<std::vector<float>>("MIPLGVector");

      if(m_cellID.size()!=m_MIPHG.size() || m_cellID.size()!=m_MIPLG.size() ){
        error() << "In MetaData: cellID size ("<<m_cellID.size()<<") != MIP size ("<<m_MIPHG.size()<<" / "<<m_MIPLG.size()<<") "<<endmsg;
      }
      else{
        for(int i=0; i<m_cellID.size(); i++){
          map_MIPHG[m_cellID[i]] = m_MIPHG[i];
          map_MIPLG[m_cellID[i]] = m_MIPLG[i];
        }
      }
    }

    if(!m_dogain){
      m_GainFactor = m_metadata->getValue<std::vector<float>>("GainIFVector");

      if(m_cellID.size()!=m_GainFactor.size()){
        error() << "In MetaData: cellID size ("<<m_cellID.size()<<") != Gain inter-factor size ("<<m_GainFactor.size()<<") "<<endmsg;
      }
      else{
        for(int i=0; i<m_cellID.size(); i++){
          map_GainFactor[m_cellID[i]] = m_GainFactor[i];
        }
      }
    }
		
	}
  
	if(m_dopedestal){
    info() << "Do Pedestal calibration "<< endmsg; 
    map_pedestalHG.clear();
    map_pedestalLG.clear();
    map_pedestalHGsigma.clear();
    map_pedestalLGsigma.clear();
  }
  if(m_domip){
    info() << "Do MIP calibration "<< endmsg; 
    map_MIPHG.clear();
    map_MIPLG.clear();
  }
  if(m_dogain){
    info() << "Do Gain calibration "<< endmsg; 
    map_GainFactor.clear();
  }
  if(m_doana){
    info() << "Do Data analysis "<< endmsg;
    if(map_pedestalHG.size()==0 || map_pedestalLG.size()==0 || map_pedestalHGsigma.size()==0 || map_pedestalLGsigma.size()==0){
      error() << "Lack of pedestal data! " << endmsg;
      return StatusCode::FAILURE;
    }
    if(map_MIPHG.size()==0 || map_MIPLG.size()==0){
      error() << "Lack of MIP data! " << endmsg;
      return StatusCode::FAILURE;
    }
    if(map_GainFactor.size()==0){
      error() << "Lack of Gain factor data! " << endmsg;
      return StatusCode::FAILURE;
    }
  }

  //Create cellID map
  std::ifstream fMapFile(m_channel_map.value());
  if (!fMapFile.is_open())
  {
      std::cerr << "Error: Failed to open file: " << m_channel_map.value() << std::endl;
      return StatusCode::FAILURE;
  }
  std::string line;
  for (Int_t i = 0; i < m_skipline_map; ++i)
  {
      if (!std::getline(fMapFile, line))
      {
          std::cerr << "Error: Unexpected end of file while skipping header lines in: " << m_channel_map.value() << std::endl;
          return StatusCode::FAILURE;
      }
  }
  while (std::getline(fMapFile, line))
  {
      std::istringstream iss(line);
      Int_t board = -1, channel = -1, layer = -1, index = -1, side = -1, crystalid = -1, cellid = -1;
      if (!(iss >> board >> channel >> layer >> index >> side >> crystalid >> cellid))
      {
          std::cerr << "Error: Failed to parse line: " << line << std::endl;
          continue;
      }
      std::vector<Int_t> cellidvec(4);
      ChannelMap[cellid] = {board, channel, layer, index, side};
  }
  fMapFile.close();
/*
std::cout<<"Print pedestal"<<std::endl;
for(auto iter : map_pedestalHG){
  int cellid = iter.first;
  std::cout<<ChannelMap[cellid].Board<<'\t'<<ChannelMap[cellid].Channel<<'\t'<<map_pedestalHG[cellid]<<'\t'<<map_pedestalLG[cellid]<<std::endl;
}

std::cout<<"Print Gain factor"<<std::endl;
for(auto iter : map_GainFactor){
  int cellid = iter.first;
  std::cout<<ChannelMap[cellid].Board<<'\t'<<ChannelMap[cellid].Channel<<'\t'<<map_GainFactor[cellid]<<std::endl;
}

std::cout<<"Print MIP calibration"<<std::endl;
for(auto iter : map_MIPHG){
  int cellid = iter.first;
  std::cout<<ChannelMap[cellid].Board<<'\t'<<ChannelMap[cellid].Channel<<'\t'<<map_MIPHG[cellid]<<'\t'<<map_MIPLG[cellid]<<std::endl;
}
*/
  //m_geosvc = service<IGeomSvc>("GeomSvc");
  //if(!m_geosvc)
  //{
  //  error() << "Failed to get the GeomSvc" << endmsg;
  //  return StatusCode::FAILURE;
  //}
  //auto m_dd4hep = m_geosvc->lcdd();
  //if(!m_dd4hep)
  //{
  //  error() << "Failed to get the Detector from GeomSvc" << endmsg;
  //  return StatusCode::FAILURE;
  //}
  //m_cellIDConverter = new dd4hep::rec::CellIDPositionConverter(*m_dd4hep);
  //if(!m_cellIDConverter)
  //{
  //  error() << "Failed to get the m_cellIDConverter." << endmsg;
  //  return StatusCode::FAILURE;
  //}

 
  info() << "CalibECALPrototypeAlg::initialized" << endmsg;
  return GaudiAlgorithm::initialize();
}


StatusCode CalibECALPrototypeAlg::execute()
{
  if(m_nEvt%1000==0) info() << "Processing " << m_nEvt << " event: " << endmsg;

	const edm4hep::RawCalorimeterHitCollection* m_inputRawHGHitCol = nullptr;
	const edm4hep::RawCalorimeterHitCollection* m_inputRawLGHitCol = nullptr;
  edm4hep::CalorimeterHitCollection* w_caloVec = m_outputHitHdl.createAndPut();


	try{
		m_inputRawHGHitCol = m_inputRawHGHitHdl.get();
		m_inputRawLGHitCol = m_inputRawLGHitHdl.get();
	}catch(GaudiException &e){
		error() << "Miss input collection! " << endmsg;
		return StatusCode::FAILURE;
	}

	//Read data into vector
  int ChannelAboveThreshold = 0;
	for(auto hit : *m_inputRawHGHitCol){
		int cellID = hit.getCellID();
    if(cellID<0) continue;

    map_HGVec[cellID].push_back(hit.getAmplitude());

    if(hit.getAmplitude()> map_pedestalHG[cellID] + MIPCut_OverPedestalSigma*map_pedestalHGsigma[cellID]){
      ChannelAboveThreshold++;
    }
	}

  for(auto hit : *m_inputRawLGHitCol){
    int cellID = hit.getCellID();
    if(cellID<0) continue;

    map_LGVec[cellID].push_back(hit.getAmplitude());
  }

//std::cout<<"Channel over threshold: "<<ChannelAboveThreshold<<std::endl;

  if(ChannelAboveThreshold >= HitNumberLimitMIP){

    for(auto hit : *m_inputRawHGHitCol){
      int cellID = hit.getCellID();
      if(cellID<0) continue;

      //Note: use channel over threshold.
      //TODO: use MIP position cut?
      if(hit.getAmplitude()> map_pedestalHG[cellID] + MIPCut_OverPedestalSigma*map_pedestalHGsigma[cellID]){
        float MIP_ADC_HG = hit.getAmplitude() -  map_pedestalHG[cellID];
        if(MIP_ADC_HG>0) map_HGMipVec[cellID].push_back(MIP_ADC_HG);
      }
    }

    for(auto hit : *m_inputRawLGHitCol){
      int cellID = hit.getCellID();
      if(cellID<0) continue;

      if(map_HGMipVec.find(cellID)!=map_HGMipVec.end()){
        float MIP_ADC_LG = hit.getAmplitude() -  map_pedestalLG[cellID];
        if(MIP_ADC_LG>0) map_LGMipVec[cellID].push_back(MIP_ADC_LG);
      }
    }
  }

	//Analysis data
	//Note: data analyse and calibration can not be performed simultaneously. 
  if(m_doana){
    for(auto hit_HG : *m_inputRawHGHitCol){
      int cellID = hit_HG.getCellID();
      if(cellID<0) continue;  
      edm4hep::RawCalorimeterHit hit_LG; 
      for(auto iter : *m_inputRawLGHitCol){
        if(iter.getCellID()==cellID){
          hit_LG = iter;
          break;
        }
      }
//std::cout<<"  cellID "<<cellID<<", raw HG "<<hit_HG.getAmplitude()<<", raw LG "<<hit_LG.getAmplitude()<<std::endl;

      int ADC_LG_This = 8*(int)(hit_LG.getAmplitude()/8);
      int ADC_HG_This = 8*(int)(hit_HG.getAmplitude()/8);
      bool kSwitch = (ADC_HG_This > GainSwitchPoint);

      float chEnergyHG = (ADC_HG_This - map_pedestalHG[cellID])/map_MIPHG[cellID]*EcalMIPEnergy/2.;
      float chEnergyLG = (ADC_LG_This - map_pedestalLG[cellID])/(map_MIPHG[cellID]*map_GainFactor[cellID])*EcalMIPEnergy/2.;
      float chEnergy = 0;
      if(!kSwitch){
        chEnergy = chEnergyHG;
      }
      else if(chEnergyLG>0.9*chEnergyHG){
        chEnergy = chEnergyLG;
      }
      else chEnergy = 0.;
//std::cout<<"  HG energy "<<chEnergyHG<<", LG energy "<<chEnergyLG<<", final energy "<<chEnergy<<std::endl;
      if(chEnergy<EcalMIPCut*EcalMIPEnergy/2.) continue;

      //get position from cellID
      double pos_x, pos_y, pos_z;
      pos_x = 1000. + (ChannelMap[cellID].Layer-1)* 20. + 10.;
      if(ChannelMap[cellID].Layer % 2 ==0){
        pos_z=0;
        pos_y = -60. + (ChannelMap[cellID].Index-1)*20. + 10;
      }
      else{
        pos_y=0;
        pos_z = -60. + (ChannelMap[cellID].Index-1)*20. + 10;
      }
      edm4hep::Vector3f m_pos(pos_x, pos_y, pos_z);

//std::cout<<"  position: ("<<pos_x<<", "<<pos_y<<", "<<pos_z<<" ) "<<std::endl;

      auto calohit = w_caloVec->create();
      calohit.setCellID(cellID);
      calohit.setEnergy(chEnergy/1000.);  //MeV to GeV
      calohit.setTime(hit_HG.getTimeStamp());
      calohit.setPosition(m_pos);      
    }
  }

  m_nEvt++;
  return StatusCode::SUCCESS;
}

StatusCode CalibECALPrototypeAlg::finalize()
{

	//Fit calibration 
	if(m_dopedestal){
		for(const auto& iter : map_HGVec){
      if(iter.second.size()<50){
        error() << "Error: Low histogram entry: "<<iter.second.size()<<", HG cellID " << iter.first<<endmsg;
        continue;
      }

      TH1D* h_hist = new TH1D(Form("HG_Cell%d",iter.first),"",140,5,145);
      for(int i=0; i<iter.second.size(); i++) h_hist->Fill(iter.second[i]);

			double PedMPVHG = h_hist->GetBinCenter(h_hist->GetMaximumBin());
			auto PedResultHG = h_hist->Fit("gaus", "ERSMQ", "", PedMPVHG - 20, PedMPVHG + 20);
      PedMPVHG = PedResultHG->Parameter(1);
      double PedSigmaHG = PedResultHG->Parameter(2); 
      if(PedMPVHG-2*PedSigmaHG<0) PedMPVHG = h_hist->GetBinCenter(h_hist->GetMaximumBin());

      PedResultHG = h_hist->Fit("gaus", "ERSMQ", "", PedMPVHG-2*PedSigmaHG, PedMPVHG+2*PedSigmaHG);
      PedMPVHG = PedResultHG->Parameter(1);
      PedSigmaHG = PedResultHG->Parameter(2); 
      if(PedMPVHG-2*PedSigmaHG<0) PedMPVHG = h_hist->GetBinCenter(h_hist->GetMaximumBin());

      PedResultHG = h_hist->Fit("gaus", "ERSMQ", "", PedMPVHG-2*PedSigmaHG, PedMPVHG+2*PedSigmaHG);
      PedMPVHG = PedResultHG->Parameter(1);
      PedSigmaHG = PedResultHG->Parameter(2); 
      if(PedMPVHG-2*PedSigmaHG<0) PedMPVHG = h_hist->GetBinCenter(h_hist->GetMaximumBin());

      PedResultHG = h_hist->Fit("gaus", "ERSMQ", "", PedMPVHG-2*PedSigmaHG, PedMPVHG+2*PedSigmaHG);			

			map_pedestalHG[iter.first] = PedResultHG->Parameter(1);
			map_pedestalHGsigma[iter.first] = PedResultHG->Parameter(2);
      delete h_hist;
		}

    for(const auto& iter : map_LGVec){
      if(iter.second.size()<50){
        error() << "Error: Low histogram entry: "<<iter.second.size()<<", LG cellID " << iter.first<<endmsg;
        continue;
      }

      TH1D* h_hist = new TH1D(Form("LG_Cell%d",iter.first),"",140,5,145);
      for(int i=0; i<iter.second.size(); i++) h_hist->Fill(iter.second[i]);

      double PedMPVLG = h_hist->GetBinCenter(h_hist->GetMaximumBin());
      auto PedResultLG = h_hist->Fit("gaus", "ERSMQ", "", PedMPVLG - 20, PedMPVLG + 20);
      PedMPVLG = PedResultLG->Parameter(1);
      double PedSigmaLG = PedResultLG->Parameter(2);
      if(PedMPVLG-2*PedSigmaLG<0) PedMPVLG = h_hist->GetBinCenter(h_hist->GetMaximumBin());

      PedResultLG = h_hist->Fit("gaus", "ERSMQ", "", PedMPVLG-2*PedSigmaLG, PedMPVLG+2*PedSigmaLG);
      PedMPVLG = PedResultLG->Parameter(1);
      PedSigmaLG = PedResultLG->Parameter(2);
      if(PedMPVLG-2*PedSigmaLG<0) PedMPVLG = h_hist->GetBinCenter(h_hist->GetMaximumBin());

      PedResultLG = h_hist->Fit("gaus", "ERSMQ", "", PedMPVLG-2*PedSigmaLG, PedMPVLG+2*PedSigmaLG);
      PedMPVLG = PedResultLG->Parameter(1);
      PedSigmaLG = PedResultLG->Parameter(2);
      if(PedMPVLG-2*PedSigmaLG<0) PedMPVLG = h_hist->GetBinCenter(h_hist->GetMaximumBin());

      PedResultLG = h_hist->Fit("gaus", "ERSMQ", "", PedMPVLG-2*PedSigmaLG, PedMPVLG+2*PedSigmaLG);

      map_pedestalLG[iter.first] = PedResultLG->Parameter(1);
      map_pedestalLGsigma[iter.first] = PedResultLG->Parameter(2);
      delete h_hist;
    }
    std::cout<<"End fitting pedestal "<<std::endl;


		if(m_cellID.size()==0){
			for(auto& iter : map_pedestalHG){
				m_cellID.push_back(iter.first);
			}
		}

		for(auto& iter : map_pedestalHG){
      m_pedestalHG.push_back(static_cast<float>(iter.second));
      m_pedestalHGsigma.push_back(static_cast<float>(map_pedestalHGsigma[iter.first]));
    }

    for(auto& iter : map_pedestalLG){
      m_pedestalLG.push_back(static_cast<float>(iter.second));
      m_pedestalLGsigma.push_back(static_cast<float>(map_pedestalLGsigma[iter.first]));
    }

    std::cout<<"End writing metadata"<<std::endl;
	}

  if(m_dogain){

    for(const auto& iter : map_HGVec){
      int cellID = iter.first;
      if(map_LGVec.find(cellID)==map_LGVec.end() || map_LGVec[cellID].size()!=iter.second.size() ){
        error() << "Error: HG and LG not match " << endmsg;
        continue;
      }
      if(iter.second.size()<50 || map_LGVec[cellID].size()<50){
        error() << "Error: Low histogram entry: "<<iter.second.size()<<" / "<<map_LGVec[cellID].size()<<", cellID " << cellID <<endmsg;
        continue;
      }

      TH2D* h_hist = new TH2D(Form("HGLG_Cell%d",cellID), "", 500, -1, 1.5*m_GainLimit.value(), 500, -1, 2500);
      for(int i=0; i<iter.second.size(); i++)
        h_hist->Fill(iter.second[i], map_LGVec[cellID][i]);

      //Fit gain ratio from 2D histograms
      auto FitLinear = new TF1("FitLinear", "[1]*x+[0]", 0, 0);
      h_hist->Fit(FitLinear, "ERSMQ", "", m_GainLimit.value()*0.3, m_GainLimit.value());

      map_GainFactor[cellID] = FitLinear->GetParameter(1);
      delete h_hist;
    }

    if(m_cellID.size()==0){
      for(auto& iter : map_GainFactor){
        m_cellID.push_back(iter.first);
      }
    }

    for(auto& iter : map_GainFactor){
      m_GainFactor.push_back(iter.second);
    }

  }

  if(m_domip){
    if(map_pedestalHG.size()==0 || map_pedestalLG.size()==0)
      error() << "Can not do MIP calibration w/o pedestal. " << endmsg;
    else{

      for(const auto& iter : map_HGMipVec){
        if(iter.second.size()<50){
          error() << "Error: Low histogram entry: "<<iter.second.size()<<", HG cellID " << iter.first<<endmsg;
          continue;
        }
        int cellID = iter.first;

        TH1D* h_hist = new TH1D(Form("MIPHG_%d", cellID), "", 200, 200, 3000);
        for(int i=0; i<iter.second.size(); i++)
          h_hist->Fill(iter.second[i]);

        double MPVHG = h_hist->GetBinCenter(h_hist->GetMaximumBin());
        //auto ResultHG = h_hist->Fit("landau", "ERSMQ", "", MPVHG*0.8 , h_hist->GetXaxis()->GetXmax());
        //MPVHG = ResultHG->Parameter(1);

        map_MIPHG[cellID] = MPVHG;
        delete h_hist;
      }

      for(const auto& iter : map_LGMipVec){
        if(iter.second.size()<50){
          error() << "Error: Low histogram entry: "<<iter.second.size()<<", LG cellID " << iter.first<<endmsg;
          continue;
        }
        int cellID = iter.first;

        TH1D* h_hist = new TH1D(Form("MIPLG_%d", cellID), "", 200, 200, 3000);
        for(int i=0; i<iter.second.size(); i++)
          h_hist->Fill(iter.second[i]);

        double MPVLG = h_hist->GetBinCenter(h_hist->GetMaximumBin());
        //auto ResultLG = h_hist->Fit("landau", "ERSMQ", "", MPVLG*0.8 , h_hist->GetXaxis()->GetXmax());
        //MPVLG = ResultLG->Parameter(1);

        map_MIPLG[cellID] = MPVLG;
        delete h_hist;
      }
    }

    if(m_cellID.size()==0){
      for(auto& iter : map_MIPHG){
        m_cellID.push_back(iter.first);
      }
    }

    for(auto& iter : map_MIPHG){
      m_MIPHG.push_back(static_cast<float>(iter.second));
    }

    for(auto& iter : map_MIPLG){
      m_MIPLG.push_back(static_cast<float>(iter.second));
    }
  }

	m_pedestalHGHdl.put(m_pedestalHG);
	m_pedestalHGSigmaHdl.put(m_pedestalHGsigma);
	m_pedestalLGHdl.put(m_pedestalLG);
	m_pedestalLGSigmaHdl.put(m_pedestalLGsigma);
  m_GainFactorHdl.put(m_GainFactor);
  m_cellIDMetaDataHdl.put(m_cellID);
  m_MIPHGHdl.put(m_MIPHG);
  m_MIPLGHdl.put(m_MIPLG);

  if(_writeNtuple)
  {
    //m_wfile->cd();
    //m_tree->Write();

    //if(m_dopedestal){
    //  for(const auto& iter : map_HGVec)
    //    iter.second->Write();
    //  for(const auto& iter : map_LGVec)
    //    iter.second->Write();
    //}
    //m_wfile->Close();
    //delete m_tree;
    //delete m_wfile;
  }

  info() << "Processed " << m_nEvt << " events " << endmsg;
  return GaudiAlgorithm::finalize();
}

#endif
