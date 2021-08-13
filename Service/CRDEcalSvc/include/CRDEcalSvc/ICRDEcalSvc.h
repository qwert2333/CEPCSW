#ifndef ICRDECALSvc_h
#define ICRDECALSvc_h

#include "GaudiKernel/IService.h"
#include "Objects/CRDCaloBlock.h"

using namespace CRDEcalEDM;

class ICRDEcalSvc: virtual public IInterface {
public:
    DeclareInterfaceID(ICRDEcalSvc, 0, 1); // major/minor version
    virtual ~ICRDEcalSvc() = default;

    //In Digi: save the DigiSystem(CRDCaloBlock) in service
    virtual void setDigiSystem( std::vector<CRDEcalEDM::CRDCaloBlock>& blockVec ) = 0;
    
    //In Recon: get the DigiSystem from service
    virtual void getDigiSystem( std::vector<CRDEcalEDM::CRDCaloBlock>& blockVec) = 0;

    virtual void ClearSystem() = 0;

private:
    std::vector<CRDEcalEDM::CRDCaloBlock> m_CaloBlockVec;

};

#endif
