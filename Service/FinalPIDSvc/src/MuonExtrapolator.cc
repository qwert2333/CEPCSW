#include "FinalPIDSvc/MuonExtrapolator.h"

TrackState MuonExtrapolator::extrap_Simple( const Track& track) {

    for ( const auto& st: track.getTrackStates() ){
            
        if ( st.location == TrackState::AtLastHit ){
            
            TrackState st2b = st;
            TrackState st2e = st;

            LCIOTrackPropagators::PropagateLCIOToCylinder( st2b, BfieldRBound, 0, 0, 1, 1e-3 );
            LCIOTrackPropagators::PropagateLCIOToZPlane( st2e, st.referencePoint.z>0?BfieldZBound:-BfieldRBound );

            double st2b_refz = st2b.referencePoint.z;
            double st2e_refz = st2e.referencePoint.z;

            TrackState newst = fabs(st2b_refz) < fabs(st2e_refz) ? st2b : st2e;

            return newst;
        }
    }

    return {};

}

// Extrapolate to B field with momentum correction of calorimeters
TrackState MuonExtrapolator::extrap_CalCorr( const Track& track, const ReconstructedParticle& pfo) {

    if ( isnan(pfo.getMomentum()[0]) || isnan(pfo.getMomentum()[1]) || isnan(pfo.getMomentum()[2]) || pfo.getEnergy() == 0) {
        return {};
    }

    double M2=0.105658*0.105658;
    double P_original=sqrt(pfo.getMomentum()[0]*pfo.getMomentum()[0]+pfo.getMomentum()[1]*pfo.getMomentum()[1]+pfo.getMomentum()[2]*pfo.getMomentum()[2]);
    double E_original=sqrt(P_original*P_original+M2);

    std::vector<std::array<double, 4>> clusters_E, clusters_H, clusters;
    clusters_E.clear();
    clusters_H.clear();
    clusters.clear();
    for (unsigned i=0;i<pfo.clusters_size();i++) {
        int loc=Location(pfo.getClusters(i).getPosition().x,pfo.getClusters(i).getPosition().y,pfo.getClusters(i).getPosition().z);
        if (loc==1 || loc==2) {
            clusters_E.push_back({pfo.getClusters(i).getPosition().x,pfo.getClusters(i).getPosition().y,pfo.getClusters(i).getPosition().z,pfo.getClusters(i).getEnergy()});
        }
        else if (loc==3 || loc==4) {
            clusters_H.push_back({pfo.getClusters(i).getPosition().x,pfo.getClusters(i).getPosition().y,pfo.getClusters(i).getPosition().z,pfo.getClusters(i).getEnergy()});
        }
    }

    std::sort(clusters_E.begin(), clusters_E.end(), [](const std::array<double, 4>& a, const std::array<double, 4>& b) {
        double normA = a[0] * a[0] + a[1] * a[1] + a[2] * a[2];
        double normB = b[0] * b[0] + b[1] * b[1] + b[2] * b[2];
        return normA < normB; 
    });
    std::sort(clusters_H.begin(), clusters_H.end(), [](const std::array<double, 4>& a, const std::array<double, 4>& b) {
        double normA = a[0] * a[0] + a[1] * a[1] + a[2] * a[2];
        double normB = b[0] * b[0] + b[1] * b[1] + b[2] * b[2];
        return normA < normB; 
    });
    
    std::copy(clusters_E.begin(), clusters_E.end(), std::back_inserter(clusters));
    std::copy(clusters_H.begin(), clusters_H.end(), std::back_inserter(clusters));

    for ( const auto& st: track.getTrackStates() ){
            
            if ( st.location == TrackState::AtLastHit ){
                
                TrackState st_tmp=st;
                double E_tmp=E_original;
                double P_tmp=P_original;

                for (size_t ic=0;ic<clusters.size();ic++) {
                    std::array<double, 4> arr=clusters.at(ic);
                    int loc=Location(arr[0],arr[1],arr[2]);
                    double e_loss=arr[3];
                    if (loc==1 || loc==2) e_loss/=1.26;
                    if (loc==1 || loc==3) {
                        if (arr[0]*arr[0]+arr[1]*arr[1]>st_tmp.referencePoint[0]*st_tmp.referencePoint[0]+st_tmp.referencePoint[1]*st_tmp.referencePoint[1])
                            LCIOTrackPropagators::PropagateLCIOToCylinder( st_tmp, sqrt(arr[0]*arr[0]+arr[1]*arr[1]), 0, 0, 1, 1e-3 );
                    }
                    else if (loc==2 || loc==4) {
                        if (fabs(arr[2])>fabs(st_tmp.referencePoint[1]))
                            LCIOTrackPropagators::PropagateLCIOToZPlane( st_tmp, fabs(arr[2]) );
                    }

                    E_tmp-=e_loss;
                    if (E_tmp*E_tmp<=M2) E_tmp+=e_loss;
                    double SF=sqrt(E_tmp*E_tmp-M2)/P_tmp;
                    // cout << st_tmp.D0 <<","<<st_tmp.Z0 <<","<<st_tmp.tanLambda <<","<<st_tmp.omega <<","<<st_tmp.phi<<","<<st_tmp.referencePoint[0]<<","<<st_tmp.referencePoint[1]<<","<<st_tmp.referencePoint[2]<<endl;
                    // cout << SF << "," << loc << endl;
                    st_tmp.omega=st_tmp.omega/SF;
                    P_tmp*=SF;
                }

                TrackState st2b = st_tmp;
                TrackState st2e = st_tmp;

                LCIOTrackPropagators::PropagateLCIOToCylinder( st2b, BfieldRBound, 0, 0, 1, 1e-3 );
                LCIOTrackPropagators::PropagateLCIOToZPlane( st2e, BfieldZBound );

                double st2b_refz = st2b.referencePoint.z;
                double st2e_refz = st2e.referencePoint.z;

                TrackState newst = fabs(st2b_refz) < fabs(st2e_refz) ? st2b : st2e;

                return newst;
    
            }
    }

    return {};

}

std::vector<double> MuonExtrapolator::angles(TrackState st, const TrackerHitCollection& hits, int variable) {
    double refx = st.referencePoint.x;
    double refy = st.referencePoint.y;
    double refz = st.referencePoint.z;

    HelixClass helix;
    helix.Initialize_Canonical( st.phi, st.D0, st.Z0, st.omega, st.tanLambda, 3.0 );
    double dirx = helix.getMomentum()[0];
    double diry = helix.getMomentum()[1];
    double dirz = helix.getMomentum()[2];

    vector<double> dRs;
    for ( const auto& hit: hits ){

        double hitx = hit.getPosition().x;
        double hity = hit.getPosition().y;
        double hitz = hit.getPosition().z;

        double dx = hitx - refx;
        double dy = hity - refy;
        double dz = hitz - refz;

        TVector3 p_track(dirx,diry,dirz);
        TVector3 x_hit(dx,dy,dz);

        double dR;
        if (variable==0) dR = x_hit.DeltaR(p_track);
        else if (variable==1) dR = x_hit.DeltaPhi(p_track);
        else dR = x_hit.Eta() - p_track.Eta();
        dRs.push_back(dR);
    }
    return dRs;
}

int MuonExtrapolator::Location(double x, double y, double z) {
    if (x*x+y*y<ECALBarrelInnerRBound*ECALBarrelInnerRBound && fabs(z)<ECALEndcapInnerZBound) return 0;//tracker
    else if (x*x+y*y<HCALBarrelInnerRBound*HCALBarrelInnerRBound && fabs(z)<ECALEndcapInnerZBound) return 1;//EB
    else if (x*x+y*y<HCALBarrelInnerRBound*HCALBarrelInnerRBound && fabs(z)<HCALEndcapInnerZBound) return 2;//EE
    else if (x*x+y*y<HCALBarrelOuterRBound*HCALBarrelOuterRBound && fabs(z)<HCALEndcapInnerZBound) return 3;//HB
    else if (x*x+y*y<HCALBarrelOuterRBound*HCALBarrelOuterRBound && fabs(z)<BfieldZBound) return 4;//HE
    else return -1;//Outside
}