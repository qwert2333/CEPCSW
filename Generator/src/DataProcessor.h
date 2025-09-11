#ifndef DataProcessor_h
#define DataProcessor_h


#include <GaudiKernel/Algorithm.h>
#include <Gaudi/Property.h>
#include <GaudiKernel/ToolHandle.h>

#include "GaudiAlg/GaudiAlgorithm.h"
#include "k4FWCore/DataHandle.h"

#include "edm4hep/RawCalorimeterHitCollection.h"
#include "edm4hep/MCParticleCollection.h"


class Event{
public:

	int TriggerID;
	double TStamp;
	
	std::vector<int> fBorderID;
	std::vector<int> fChannelID;
	std::vector<int> fADC_HG;
	std::vector<int> fADC_LG;
	std::vector<int> fToA_LSB;
	std::vector<int> fToT_LSB;
};

struct ChannelMapEntry
{
    Int_t Layer;
    Int_t Index;
    Int_t Side;
    Int_t CellID;
};

class DataProcessor: public GaudiAlgorithm {

public:
    DataProcessor(const std::string& name, ISvcLocator* pSvcLocator);

    virtual StatusCode initialize();
    virtual StatusCode execute();
    virtual StatusCode finalize();

private:
    Gaudi::Property<std::string> m_input_file{this, "InputRawDataFile", "NULL"};
		Gaudi::Property<std::string> m_channel_map{this, "ChannelMap", "NULL"};
    //Gaudi::Property<std::string> m_output_file{this, "OutputRootFile", "NULL"};
		Gaudi::Property<int> m_skipline_data{this, "RawDataHeaderLine", 9};
		Gaudi::Property<int> m_skipline_map{this, "MapFileHeaderLine", 8};

    int m_evtid;                               
    DataHandle<edm4hep::RawCalorimeterHitCollection> rawhit_hdl{"RawCaloHit", Gaudi::DataHandle::Writer, this};
		std::vector<Event> EventCol; 
		std::map<std::pair<Int_t, Int_t>, ChannelMapEntry> ChannelMap; 

};

#endif
