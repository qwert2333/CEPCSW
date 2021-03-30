#ifndef ICRDECALEDMSvc_h
#define ICRDECALEDMSvc_h

#include "GaudiKernel/IService.h"
#include "CRDCaloBlock.h"

using namespace CRDEcalEDM;

class ICRDEcalEDMSvc: virtual public IInterface {
public:
    DeclareInterfaceID(ICRDEcalEDMSvc, 0, 1); // major/minor version
    virtual ~ICRDEcalEDMSvc() = default;

    //In Digi: save the DigiSystem(CRDCaloBlock) in service
    virtual void setDigiSystem( std::vector<CRDEcalEDM::CRDCaloBlock>& blockVec ) = 0;
    
    //In Recon: get the DigiSystem from service
    virtual void getDigiSystem( std::vector<CRDEcalEDM::CRDCaloBlock>& blockVec) = 0;

    virtual void ClearSystem() = 0;

private:
    std::vector<CRDEcalEDM::CRDCaloBlock> m_CaloBlockVec;

};

#endif
