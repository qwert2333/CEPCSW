#include "PandoraPlusDataCol.h"

using namespace std;
void PandoraPlusDataCol::PrintLayer(){

  cout<<"------------------------------------------------"<<endl;
  cout<<"-----------Print out Layer Collection-----------"<<endl;
  cout<<"------------------------------------------------"<<endl;

  cout<<"Layer collection: Number = "<<LayerCol.size()<<endl;
  cout<<"Loop print all layers: "<<endl;
  for(int i=0;i<LayerCol.size();i++){

  cout<<'\t'<<"BarX collection in layer #"<<i<<endl;
  double totE=0;
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
  cout<<"2D shower collection:  Number = "<<shower2DCol.size()<<endl;
  cout<<"Loop print in layers: "<<endl;
  cout<<'\t'<<"Layer"<<'\t'<<"position (x, y, z)"<<'\t'<<"Energy"<<endl;
  double totE=0;
  for(int i=0;i<shower2DCol.size();i++){
    cout<<'\t'<<shower2DCol[i].getDlayer()<<'\t';
    cout<<shower2DCol[i].getPos().x()<<'\t'<<shower2DCol[i].getPos().y()<<'\t'<<shower2DCol[i].getPos().z()<<'\t';
    cout<<shower2DCol[i].getShowerE()<<endl;
    totE+=shower2DCol[i].getShowerE();
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
  }
  
  cout<<"------------------------------------------------"<<endl;
  cout<<"--------------------End print-------------------"<<endl;

}

void PandoraPlusDataCol::Clear(){
  Flag_Iter = false; 
  ClearBlock();  
  ClearLayer();  
  ClearShower(); 
  ClearCluster();
  ClearPFO();    

}
