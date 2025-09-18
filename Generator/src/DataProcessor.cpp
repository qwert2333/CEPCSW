#include "DataProcessor.h"

#include "GaudiKernel/IEventProcessor.h"
#include "GaudiKernel/IAppMgrUI.h"
#include "GaudiKernel/GaudiException.h"

#include <iostream>
#include <vector>
#include <fstream>


DECLARE_COMPONENT(DataProcessor)

DataProcessor::DataProcessor(const std::string& name, ISvcLocator* pSvcLocator): GaudiAlgorithm(name, pSvcLocator) {
		declareProperty("RawHGCaloHit", rawHGhit_hdl, "High-gain Raw hit collection name");
		declareProperty("RawLGCaloHit", rawLGhit_hdl, "Low-gain Raw hit collection name");
    m_evtid = 0;

}

StatusCode
DataProcessor::initialize() {

		//Read raw data file
		std::ifstream DataFile(m_input_file);
		if(!DataFile.is_open()){
			std::cout<<"Error: Raw data file does not exist. "<<std::endl;
			return StatusCode::FAILURE;
		}

    std::string DataLine;
		//Skip first 9 lines for ASIC information
    for (Int_t i = 0; i < m_skipline_data; ++i){
        std::getline(DataFile, DataLine);
				//std::cout<<"Header line: "<<DataLine<<std::endl;
		}

		EventCol.clear();
    Double_t Tstamp_us = 0;
    Int_t fTriggerID = 0;
    Int_t NHits = 0;
    Int_t fBoardID = -1;
    Int_t fChannelID = -1;
    Int_t fADC_LG = -1;
    Int_t fADC_HG = -1;
		Int_t fTriggerID_tmp = -1;
		Double_t Tstamp_us_tmp = -1;
    while (std::getline(DataFile, DataLine))
    {
        std::istringstream iss(DataLine);
        TString TempStr;

        Int_t NumCount = 0;
        while (iss >> TempStr)
            NumCount++;

        iss.clear();
				iss.seekg(0);
        if ( NumCount == 7 ) //New board
        {
          iss >> fBoardID;
          iss >> fChannelID;
          iss >> fADC_LG;
          iss >> fADC_HG;
          iss >> Tstamp_us;
          iss >> fTriggerID;
          iss >> NHits;

					fTriggerID_tmp = fTriggerID;
					Tstamp_us_tmp = Tstamp_us; 
					auto iter = std::find_if(EventCol.begin(), EventCol.end(), [&fTriggerID](const Event& evt) { return evt.TriggerID == fTriggerID;});
					if(iter==EventCol.end()){
						Event m_newEvt;
						m_newEvt.TriggerID = fTriggerID;
						m_newEvt.TStamp = Tstamp_us;
						m_newEvt.fBorderID.push_back(fBoardID);
						m_newEvt.fChannelID.push_back(fChannelID);
						m_newEvt.fADC_HG.push_back(fADC_HG);
						m_newEvt.fADC_LG.push_back(fADC_LG);
						m_newEvt.fToA_LSB.push_back(-1);
						m_newEvt.fToT_LSB.push_back(-1);
						EventCol.push_back(m_newEvt);
					}
					else{
            iter->TriggerID = fTriggerID;
            iter->TStamp = Tstamp_us;
            iter->fBorderID.push_back(fBoardID);
            iter->fChannelID.push_back(fChannelID);
            iter->fADC_HG.push_back(fADC_HG);
            iter->fADC_LG.push_back(fADC_LG);
            iter->fToA_LSB.push_back(-1);
            iter->fToT_LSB.push_back(-1);
					}
				}
				else if ( NumCount == 4 ){
          iss >> fBoardID;
          iss >> fChannelID;
          iss >> fADC_LG;
          iss >> fADC_HG;


          auto iter = std::find_if(EventCol.begin(), EventCol.end(), [&fTriggerID_tmp](const Event& evt) { return evt.TriggerID == fTriggerID_tmp;});
					if(iter==EventCol.end()){
						std::cerr<<"Error: unknown trigger ID: "<<fTriggerID_tmp<<std::endl;
						continue;
					}
					else{
            //iter->TriggerID = fTriggerID;
            //iter->TStamp = Tstamp_us;
            iter->fBorderID.push_back(fBoardID);
            iter->fChannelID.push_back(fChannelID);
            iter->fADC_HG.push_back(fADC_HG);
            iter->fADC_LG.push_back(fADC_LG);
            iter->fToA_LSB.push_back(-1);
            iter->fToT_LSB.push_back(-1);
					}

				}
				else{
					std::cerr<<" Unknown data size in line: "<<DataLine<<std::endl;
					continue;
				}	
    }		
		DataFile.close();

		std::cout<<"Event size: "<<EventCol.size()<<std::endl;

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
        ChannelMap[std::make_pair(board, channel)] = {layer, index, side, cellid};
    }
    fMapFile.close();
    
    return StatusCode::SUCCESS;
}

StatusCode
DataProcessor::execute() {
		edm4hep::RawCalorimeterHitCollection *rawHGhitcol = rawHGhit_hdl.createAndPut();
		edm4hep::RawCalorimeterHitCollection *rawLGhitcol = rawLGhit_hdl.createAndPut();

	  if(m_evtid>=EventCol.size()){
	  	//std::cout<<"Finished all events in file: "<<m_evtid<<std::endl;
      warning() << "Have read all events, stop now." << endmsg;

      auto ep = serviceLocator()->as<IEventProcessor>();
      if ( !ep ) {
          error() << "Cannot get IEventProcessor" << endmsg;
          return StatusCode::FAILURE;
      }
      ep->stopRun();
      return StatusCode::SUCCESS;
	  }

		//Construct event
		if(m_evtid%100==0) std::cout<<"  Processing event #"<<m_evtid<<std::endl;

		auto event = EventCol[m_evtid];

		for(int ihit=0; ihit<event.fBorderID.size(); ihit++){
			//Derive cell ID
			auto key = std::make_pair(event.fBorderID[ihit], event.fChannelID[ihit]);
			auto it = ChannelMap.find(key);
			if(it==ChannelMap.end()){
          std::cerr << "Error: Key not found in map for BoardID: " << event.fBorderID[ihit] << ", ChannelID: " << event.fChannelID[ihit] << std::endl;
          continue;
			}

			//Create HG hit
			auto rawhit_HG = rawHGhitcol->create();
			rawhit_HG.setCellID(it->second.CellID);
			rawhit_HG.setAmplitude(event.fADC_HG[ihit]);
			rawhit_HG.setTimeStamp(event.TStamp);

			//Create LG hit
      auto rawhit_LG = rawLGhitcol->create();
      rawhit_LG.setCellID(it->second.CellID);
      rawhit_LG.setAmplitude(event.fADC_LG[ihit]);
      rawhit_LG.setTimeStamp(event.TStamp);
		}

    m_evtid++;

    return StatusCode::SUCCESS;
}

StatusCode
DataProcessor::finalize() {
		std::cout<<"Finalizing data processing"<<std::endl;

		
		EventCol.clear();
		ChannelMap.clear();
    return StatusCode::SUCCESS;
}
