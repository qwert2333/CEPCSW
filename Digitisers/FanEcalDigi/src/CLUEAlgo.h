#ifndef CLUEAlgo_h
#define CLUEAlgo_h

// C/C++ headers
#include <set>
#include <string>
#include <vector>
#include <iostream>
#include <fstream>

#include <chrono>

#include "LayerTiles.h"
#include "Points.h"
#include "TTree.h"
#include "TFile.h"
#include <unordered_map>

using namespace std;
class CLUEAlgo{

  public:
    // constructor
    CLUEAlgo(float dc, float deltao, float deltac, float rhoc, bool verbose=false ){ 
      dc_ = dc; 
      deltao_ = deltao; 
      deltac_ = deltac; 
      rhoc_ = rhoc;
      dm_ = std::min(deltao_, deltac_);//max
      verbose_ = verbose;
    
    }
    // distrcutor
    ~CLUEAlgo(){} 
    
    // public variables
    float dc_, dm_, deltao_, deltac_, rhoc_;
    bool verbose_;
    Points points_;

    // public methods
    void setPoints(int n, float* x, float* y, int* layer, float* weight) {
      points_.clear();
      points_.n = n;
      // input variables
      points_.x.assign(x, x + n);
      points_.y.assign(y, y + n);
      points_.layer.assign(layer, layer + n);
      points_.weight.assign(weight, weight + n);
      // result variables
      points_.rho.resize(n,0);
      points_.delta.resize(n,std::numeric_limits<float>::max());
      points_.nearestHigher.resize(n,-1);
      points_.isSeed.resize(n,0);
      points_.followers.resize(n);
      points_.clusterIndex.resize(n,-1);
    }

    void setGroup( Points& groups, int n, float* phi, int* clusterID, float* weight){
      groups.clear();
      groups.n = n;
      groups.phi.assign(phi, phi + n);
      groups.clusterIndex.assign(clusterID, clusterID + n);
      groups.weight.assign(weight, weight + n);
    }

    void clearPoints(){ points_.clear(); }

    void makeClusters();

    void verboseToCluster(Points& AlgoClus, vector<CaloCluster>& m_clus, int nVerbose = -1){

	m_clus.clear();
        AlgoClus = points_;

        if (nVerbose ==-1) nVerbose=points_.n;

	// verbose to cluster
        std::unordered_map<int, vector<int>> dataMap;

        for(size_t i=0;i<points_.clusterIndex.size();++i){
	   int currentID = points_.clusterIndex[i];
	   dataMap[currentID].push_back(i);
	}
	for(const auto& pair: dataMap){
	   int groupID = pair.first;
	   const vector<int>& groupPoints = pair.second;
//	   cout<<"Gruop "<<groupID<<"contains points: ";

	   CaloCluster nClus; nClus.Clear();
	   for(int index : groupPoints){
	      CaloBar newbar;
	      newbar.setcellID(1,points_.x[index],points_.y[index]);
              newbar.setQ(points_.weight[index],points_.weight[index]);
              newbar.setClusterID(points_.clusterIndex[index]);
              nClus.AddBar(newbar);
//	      cout<<points_.x[index]<< " "<<points_.y[index]<<" ";
	   }
	   m_clus.push_back(nClus);
//	   cout<<endl;
	}
    }

    void verboseResults( Points& AlgoClus, std::string outputFileType = "cout", std::string outputFileName = "cout", int nVerbose = -1){ 
      
      if (verbose_) {
        
	AlgoClus = points_;
      
        if (nVerbose ==-1) nVerbose=points_.n;

        // verbose to screen
        if (outputFileType.compare("cout") == 0 )  {
          std::cout << "index,x,y,layer,weight,rho,delta,nh,isSeed,clusterId"<< std::endl;
          for(int i = 0; i < nVerbose; i++) {
            std::cout << i << ","<<points_.x[i]<< ","<<points_.y[i]<< ","<<points_.layer[i] << ","<<points_.weight[i];
            std::cout << "," << points_.rho[i];
            if (points_.delta[i] <= 999) 
              std::cout << ","<<points_.delta[i];
            else
              std::cout << ",999"; // convert +inf to 999 in verbose
            std::cout << ","<<points_.nearestHigher[i];
            std::cout << "," << points_.isSeed[i];
            std::cout << ","<<points_.clusterIndex[i];
            std::cout << std::endl;
          }
        }

        // verbose to file
        if(outputFileType.compare("csvFile") == 0) {

          std::string outputFileCSV = outputFileName+".csv";
          std::ofstream oFile(outputFileCSV);
          oFile << "index    x    y    layer    weight    rho    delta    nh    isSeed    clusterId\n";
          for(int i = 0; i < nVerbose; i++) {
            oFile << i << "    "<<points_.x[i]<< "    "<<points_.y[i]<< "    "<<points_.layer[i] << "    "<<points_.weight[i];
            oFile << "    " << points_.rho[i];
            if (points_.delta[i] <= 999) 
              oFile << "    "<<points_.delta[i];
            else
              oFile << "    999"; // convert +inf to 999 in verbose
            oFile << "    "<<points_.nearestHigher[i];
            oFile << "    " << points_.isSeed[i];
            oFile << "    clusterIndex: "<<points_.clusterIndex[i];
            oFile << "\n";
          }
          oFile.close();
        }
        //verbose to tree
        if(outputFileType.compare("tree") == 0){
          vector<float> crystal;    crystal.clear();
          vector<float> module;     module.clear();
          vector<float> E;          E.clear();
          vector<int> clusterID;    clusterID.clear();
          std::string outputFileRoot = outputFileName+".root";
          char outputRoot[50];
          strcpy(outputRoot,outputFileRoot.c_str());

          TFile *w_file = new TFile(outputRoot, "recreate");
          TTree *t_RecClus = new TTree("RecClus", "RecClus");
          t_RecClus->Branch("crystal", &crystal);
          t_RecClus->Branch("module", &module);
          t_RecClus->Branch("E", &E);
          t_RecClus->Branch("clusterID", &clusterID);

          for(int i = 0; i < nVerbose; i++) {
            crystal.push_back(points_.x[i]);
            module.push_back(points_.y[i]);
            E.push_back(points_.weight[i]);
            clusterID.push_back(points_.clusterIndex[i]);
          }
          t_RecClus->Fill();
          t_RecClus->Write();
        }
	
      }// end of if verbose_
        
    }
   
    
    double getClusEnergy(Points& AlgoClus, int ClusIndex = 0){
      AlgoClus = points_;
      double ClusE = 0;
      for(int i=0;i<AlgoClus.n;i++){
        if(/*AlgoClus.clusterIndex[i] != -1 && */AlgoClus.clusterIndex[i] == ClusIndex){
           ClusE += AlgoClus.weight[i];
        }
      }
      return ClusE;
    }

    void getAveClusPhi(Points& AlgoClus, int& ClusIndex, int& ClusSize, double& ClusPhi, double& clusE){

      AlgoClus = points_;
      ClusSize=0, ClusPhi=0;
      clusE = 0;
      clusE = getClusEnergy(AlgoClus,ClusIndex);
  //    cout<<"Here list cluster "<<ClusIndex<<endl;
      for(int i=0;i<AlgoClus.n;i++){
        if(/*AlgoClus.clusterIndex[i] != -1 && */AlgoClus.clusterIndex[i] == ClusIndex){
           ClusPhi += (double) AlgoClus.x[i]*AlgoClus.weight[i]/clusE;
//	   cout<<AlgoClus.x[i]<<" "<<AlgoClus.y[i]<<" "<<AlgoClus.weight[i]<<" "<<ClusIndex<<"      "<<AlgoClus.clusterIndex[i]<<endl;
           ClusSize +=1;
        }
      }
    }
/*
    Points GetPoints(Points& evtClus){
      evtClus.clear();
      Points groups; groups.clear();
      vector<float> _phi;     _phi.clear();
      vector<int> _clusID;    _clusID.clear();
      vector<float> _weight;  _weight.clear();

      evtClus = points_;
      int nOddClus = evtClus.NumberofCluster();
      int size=0;
      double phi=0, weight=0;
      for(int j=-1;j<nOddClus-1;j++){
          getAveClusPhi(evtClus, j, size, phi, weight);
              _phi.push_back(phi);
              _clusID.push_back(j);
	      _weight.push_back(weight);
//	      cout<<"clusID: "<<j<<"  ,size: "<<size<<"  ,phi: "<<phi<<"E: "<<weight<<endl;
      }
      setGroup(groups,_phi.size(),&_phi[0],&_clusID[0],&_weight[0]);
      cout<<"size: "<<_phi.size()<<endl;
      return groups;
   }
*/
  private:

    // private variables

    // private member methods
    void prepareDataStructures(std::array<LayerTiles, NLAYERS> & );
    void calculateLocalDensity(std::array<LayerTiles, NLAYERS> & );
    void calculateDistanceToHigher(std::array<LayerTiles, NLAYERS> & );
    void findAndAssignClusters();
    inline float distance(int , int) const ;
};

#endif
