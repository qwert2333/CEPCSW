#include <iostream>
#include <cmath>
#include "edm4hep/TrackState.h"
#include "edm4hep/TrackCollection.h"
#include "edm4hep/ReconstructedParticleCollection.h"
#include "edm4hep/SimTrackerHitCollection.h"
#include "edm4hep/TrackerHitCollection.h"

#include <TVector3.h>

#include "TrackSystemSvc/LCIOTrackPropagators.h"
#include "DataHelper/HelixClass.h"

using namespace std;
using namespace edm4hep;


class MuonExtrapolator {
public:
    MuonExtrapolator( double BfieldRBound_, double BfieldZBound_ ) : BfieldRBound(BfieldRBound_), BfieldZBound(BfieldZBound_) {}

    TrackState extrap_Simple( const Track& track);
    TrackState extrap_CalCorr( const Track& track, const ReconstructedParticle& pfo);

    std::vector<double> angles(TrackState st, const TrackerHitCollection& hits, int variable);

private:
    double BfieldRBound;
    double BfieldZBound;
    
    double ECALBarrelInnerRBound=1830;
    double ECALBarrelOuterRBound=2130;
    double HCALBarrelInnerRBound=2130;
    double HCALBarrelOuterRBound=3535;

    double ECALEndcapInnerZBound=2900;
    double ECALEndcapOuterZBound=3230;
    double HCALEndcapInnerZBound=3230;
    double HCALEndcapOuterZBound=3260+1315;
    int Location( double x, double y, double z);
    
};