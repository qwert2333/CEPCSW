#ifndef ICRDECALSvc_h
#define ICRDECALSvc_h

#include "GaudiKernel/IService.h"
#include "Objects/CRDCaloBar.h"

using namespace CRDEcalEDM;

class ICRDEcalSvc: virtual public IInterface {
public:
    DeclareInterfaceID(ICRDEcalSvc, 0, 1); // major/minor version
    virtual ~ICRDEcalSvc() = default;

    //In Digi: save the DigiHits(CRDCaloBar) in service
    virtual void setDigiHits( std::vector<CRDEcalEDM::CRDCaloBar>& digiHits ) = 0;
    
    //In Recon: get the DigiHits from service
    virtual void getDigiHits( std::vector<CRDEcalEDM::CRDCaloBar>& digiHits) = 0;

    virtual void ClearSystem() = 0;

private:
    std::vector<CRDEcalEDM::CRDCaloBar> m_digiHits;

};

#endif
