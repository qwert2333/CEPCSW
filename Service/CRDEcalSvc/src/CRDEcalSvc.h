#ifndef CRDECALSvc_h
#define CRDECALSvc_h

#include <GaudiKernel/Service.h>
#include "CRDEcalSvc/ICRDEcalSvc.h"

class CRDEcalSvc: public extends<Service, ICRDEcalSvc> {
public:
    using extends::extends;

    StatusCode initialize() override;
    StatusCode finalize() override;

    void setDigiHits( std::vector<CRDEcalEDM::CRDCaloBar>& _hits ) { m_digiHits = _hits; }
    void getDigiHits( std::vector<CRDEcalEDM::CRDCaloBar>& _hits ) { _hits = m_digiHits; }
    void ClearSystem() { m_digiHits.clear(); }

private:
    std::vector<CRDEcalEDM::CRDCaloBar> m_digiHits;

};

#endif
