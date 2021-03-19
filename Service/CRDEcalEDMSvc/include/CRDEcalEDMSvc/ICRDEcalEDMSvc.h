#ifndef ICRDECALEDMSvc_h
#define ICRDECALEDMSvc_h

#include "GaudiKernel/IService.h"
#include "CRDCaloBar.h"

using namespace CRDEcalEDM;

class ICRDEcalEDMSvc: virtual public IInterface {
public:
    DeclareInterfaceID(ICRDEcalEDMSvc, 0, 1); // major/minor version
    virtual ~ICRDEcalEDMSvc() = default;

    //In Digi: save the DigiSystem(DigiBlock) in service
    virtual void setDigiSystem( std::vector<CRDEcalEDM::DigiBlock> blockVec ) = 0;
    
    //In Recon: get the DigiSystem from service
    virtual std::vector<CRDEcalEDM::DigiBlock> getDigiSystem() = 0;

private:
    std::vector<CRDEcalEDM::DigiBlock> m_DigiBlockVec;

};

#endif
