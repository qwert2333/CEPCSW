#include "PandoraPlusDataCol.h"

using namespace std;
void PandoraPlusDataCol::PrintLayer(){

  cout<<"------------------------------------------------"<<endl;
  cout<<"-----------Print out Layer Collection-----------"<<endl;
  cout<<"------------------------------------------------"<<endl;

  cout<<"Layer collection: Number = "<<LayerCol.size()<<endl;
  cout<<"Loop print all layers: "<<endl;
  for(int i=0;i<LayerCol.size();i++){
  double totE;

  cout<<'\t'<<"BarX collection in layer #"<<i<<endl;
  totE=0;
  for(int j=0;j<LayerCol[i].barXCol.size();j++){
    printf(" \t Position: (%.2f, %.2f, %.2f), and Energy: %.4f \n", 
            LayerCol[i].barXCol[j].getPosition().x(), 
            LayerCol[i].barXCol[j].getPosition().y(), 
            LayerCol[i].barXCol[j].getPosition().z(), 
            LayerCol[i].barXCol[j].getEnergy() );   
    totE += LayerCol[i].barXCol[j].getEnergy();
  }
  cout<<'\t'<<"Total energy in BarX: "<<totE<<endl;
  totE=0;
  cout<<'\t'<<"BarY collection in layer #"<<i<<endl;
  for(int j=0;j<LayerCol[i].barYCol.size();j++){
    printf(" \t Position: (%.2f, %.2f, %.2f), and Energy: %.4f \n",
            LayerCol[i].barYCol[j].getPosition().x(),
            LayerCol[i].barYCol[j].getPosition().y(),
            LayerCol[i].barYCol[j].getPosition().z(),
            LayerCol[i].barYCol[j].getEnergy() );
    totE += LayerCol[i].barYCol[j].getEnergy();
  }
  cout<<'\t'<<"Total energy in BarY: "<<totE<<endl;
  cout<<'\t'<<"------------------------------------------"<<endl;


  cout<<'\t'<<"ShowerX collection in layer #"<<i<<endl;
  totE=0;
  for(int j=0;j<LayerCol[i].barShowerXCol.size();j++){
    printf(" \t Position: (%.2f, %.2f, %.2f), and Energy: %.4f \n",
            LayerCol[i].barShowerXCol[j].getPos().x(),
            LayerCol[i].barShowerXCol[j].getPos().y(),
            LayerCol[i].barShowerXCol[j].getPos().z(),
            LayerCol[i].barShowerXCol[j].getE() );
    totE += LayerCol[i].barShowerXCol[j].getE();
  }
  cout<<'\t'<<"Total energy in ShowerX: "<<totE<<endl;
  totE=0;
  cout<<'\t'<<"ShowerY collection in layer #"<<i<<endl;
  for(int j=0;j<LayerCol[i].barShowerYCol.size();j++){
    printf(" \t Position: (%.2f, %.2f, %.2f), and Energy: %.4f \n",
            LayerCol[i].barShowerYCol[j].getPos().x(),
            LayerCol[i].barShowerYCol[j].getPos().y(),
            LayerCol[i].barShowerYCol[j].getPos().z(),
            LayerCol[i].barShowerYCol[j].getE() );
    totE += LayerCol[i].barShowerYCol[j].getE();
  }
  cout<<'\t'<<"Total energy in ShowerY: "<<totE<<endl;
  cout<<'\t'<<"------------------------------------------"<<endl;

  }
  cout<<"------------------------------------------------"<<endl;
  cout<<"--------------------End print-------------------"<<endl;
}

void PandoraPlusDataCol::PrintShower(){

  cout<<"------------------------------------------------"<<endl;
  cout<<"----------Print out Shower Collection-----------"<<endl;
  cout<<"------------------------------------------------"<<endl;

  cout<<"------------------------------------------------"<<endl;
  cout<<"2D shower collection:  Number = "<<Shower2DCol.size()<<endl;
  cout<<"Loop print in layers: "<<endl;
  cout<<'\t'<<"Layer"<<'\t'<<"position (x, y, z)"<<'\t'<<"Energy"<<endl;
  double totE=0;
  for(int i=0;i<Shower2DCol.size();i++){
    cout<<'\t'<<Shower2DCol[i].getDlayer()<<'\t';
    printf( "(%.2f, %.2f, %.2f) \t%.4f \n", Shower2DCol[i].getPos().x(), Shower2DCol[i].getPos().y(), Shower2DCol[i].getPos().z(), Shower2DCol[i].getShowerE() );
    totE+=Shower2DCol[i].getShowerE();
  }
  cout<<"Total Energy: "<<totE<<endl;
  cout<<"------------------------------------------------"<<endl;
  cout<<"--------------------End print-------------------"<<endl;
}

void PandoraPlusDataCol::Print3DClus(){

  cout<<"------------------------------------------------"<<endl;
  cout<<"---------Print out Cluster Collection-----------"<<endl;
  cout<<"------------------------------------------------"<<endl;

  cout<<"Good 3D cluster Number = "<<GoodClus3DCol.size()<<endl;
  cout<<"Bad 3D cluster Number = " <<BadClus3DCol.size()<<endl;
  cout<<"Merged 3D showers: Number = "<<Clus3DCol.size()<<endl;
  cout<<endl;

  cout<<"------------------------------------------------"<<endl;
  cout<<"Loop good 3Dshowers: "<<endl;
  for(int i=0;i<GoodClus3DCol.size();i++){
    cout<<'\t'<<"Shower #"<<i;
    printf("\t (%2f, %2f, %2f), Energy=%2f, axis (%.2f, %.2f, %.2f) \n", 
            GoodClus3DCol[i].getShowerCenter().x(), GoodClus3DCol[i].getShowerCenter().y(), GoodClus3DCol[i].getShowerCenter().z(),  
            GoodClus3DCol[i].getShowerE(), 
            GoodClus3DCol[i].getAxis().x(), GoodClus3DCol[i].getAxis().y(), GoodClus3DCol[i].getAxis().z());
    cout<<'\t'<<'\t'<<"2D shower depth: ";
    for(int j=0;j<GoodClus3DCol[i].get2DShowers().size();j++) cout<<GoodClus3DCol[i].get2DShowers()[j].getPos().x()<<'\t';
    cout<<endl;
    cout<<'\t'<<'\t'<<"2D shower z: ";
    for(int j=0;j<GoodClus3DCol[i].get2DShowers().size();j++) cout<<GoodClus3DCol[i].get2DShowers()[j].getPos().z()<<'\t';
    cout<<endl;
  }

  cout<<"------------------------------------------------"<<endl;
  cout<<"Loop bad 3Dshowers: "<<endl;
  for(int i=0;i<BadClus3DCol.size();i++){
    cout<<'\t'<<"Shower #"<<i;
    printf("\t (%2f, %2f, %2f), Energy=%2f, axis (%.2f, %.2f, %.2f) \n", 
            BadClus3DCol[i].getShowerCenter().x(), BadClus3DCol[i].getShowerCenter().y(), BadClus3DCol[i].getShowerCenter().z(),  
            BadClus3DCol[i].getShowerE(), 
            BadClus3DCol[i].getAxis().x(), BadClus3DCol[i].getAxis().y(), BadClus3DCol[i].getAxis().z()  );
    cout<<'\t'<<'\t'<<"2D shower depth: ";
    for(int j=0;j<BadClus3DCol[i].get2DShowers().size();j++) cout<<BadClus3DCol[i].get2DShowers()[j].getPos().x()<<'\t';
    cout<<endl;
    cout<<'\t'<<'\t'<<"2D shower z: ";
    for(int j=0;j<BadClus3DCol[i].get2DShowers().size();j++) cout<<BadClus3DCol[i].get2DShowers()[j].getPos().z()<<'\t';
    cout<<endl;
  }

  cout<<"------------------------------------------------"<<endl;
  cout<<"Loop merged 3Dshowers: "<<endl;
  for(int i=0;i<Clus3DCol.size();i++){
    cout<<'\t'<<"Shower #"<<i;
    printf("\t (%2f, %2f, %2f), Energy=%2f, axis (%.2f, %.2f, %.2f) \n", 
            Clus3DCol[i].getShowerCenter().x(), Clus3DCol[i].getShowerCenter().y(), Clus3DCol[i].getShowerCenter().z(),  
            Clus3DCol[i].getShowerE(),
            Clus3DCol[i].getAxis().x(), Clus3DCol[i].getAxis().y(), Clus3DCol[i].getAxis().z()  );
    cout<<'\t'<<'\t'<<"2D shower depth: ";
    for(int j=0;j<Clus3DCol[i].get2DShowers().size();j++) cout<<Clus3DCol[i].get2DShowers()[j].getPos().x()<<'\t';
    cout<<endl;
    cout<<'\t'<<'\t'<<"2D shower z: ";
    for(int j=0;j<Clus3DCol[i].get2DShowers().size();j++) cout<<Clus3DCol[i].get2DShowers()[j].getPos().z()<<'\t';
    cout<<endl;
  }
  
  cout<<"------------------------------------------------"<<endl;
  cout<<"--------------------End print-------------------"<<endl;

}


void PandoraPlusDataCol::PrintArborTree(){
  cout<<"------------------------------------------------"<<endl;
  cout<<"-------------Print out ArborTrees---------------"<<endl;
  cout<<"------------------------------------------------"<<endl;

  cout<<"Arbor TreeX number: "<<ArborTreeCol.size()<<endl;
  cout<<"------------------------------------------------"<<endl;
  cout<<"Loop TreeX:  "<<endl;
  cout<<"  (Barycenter, Nnode):"<<endl;
  for(int i=0; i<ArborTreeCol.size(); i++){
    printf("     (%.2f, %.2f, %.2f, %d) \n", ArborTreeCol[i].GetBarycenter().x(), ArborTreeCol[i].GetBarycenter().y(), ArborTreeCol[i].GetBarycenter().z(), ArborTreeCol[i].GetNodes().size());
  }
  cout<<"------------------------------------------------"<<endl;
  cout<<"--------------------End print-------------------"<<endl;

}

void PandoraPlusDataCol::ClearTempCol(){
  BlockVec_raw.clear();

  BlockVec_iter0.clear();
  LayerCol_iter0.clear();
  Shower2DCol_iter0.clear();
  GoodClus3DCol_iter0.clear();
  BadClus3DCol_iter0.clear();
  Clus3DCol_iter0.clear();

  BlockVec_iter1.clear();
  LayerCol_iter1.clear();
  Shower2DCol_iter1.clear();
  GoodClus3DCol_iter1.clear();
  BadClus3DCol_iter1.clear();
  Clus3DCol_iter1.clear();

  BlockVec_iter2.clear();
  LayerCol_iter2.clear();
  Shower2DCol_iter2.clear();
  GoodClus3DCol_iter2.clear();
  BadClus3DCol_iter2.clear();
  Clus3DCol_iter2.clear();
}

void PandoraPlusDataCol::Clear(){
  ClearBlock();  
  ClearLayer();  
  ClearShower(); 
  ClearCluster();
  ClearTrack();
  ClearPFO();    
  ClearTempCol();
  ClearArbor();
}
