#ifndef Points_h
#define Points_h


struct Points {
  
  std::vector<float> x;
  std::vector<float> y;
  std::vector<int> layer;
  std::vector<float> weight;
 
  std::vector<float> phi;
 
  std::vector<float> rho;
  std::vector<float> delta;
  std::vector<int> nearestHigher;
  std::vector<int> clusterIndex;
  std::vector<std::vector<int>> followers;
  std::vector<int> isSeed;
  // why use int instead of bool?
  // https://en.cppreference.com/w/cpp/container/vector_bool
  // std::vector<bool> behaves similarly to std::vector, but in order to be space efficient, it:
  // Does not necessarily store its elements as a contiguous array (so &v[0] + n != &v[n])

  int n;

  void clear() {
    x.clear();
    y.clear();
    layer.clear();
    weight.clear();

    phi.clear();

    rho.clear();
    delta.clear();
    nearestHigher.clear();
    clusterIndex.clear();
    followers.clear();
    isSeed.clear();

    n = 0;
  }

  int NumberofCluster(){ 
    int nClus = 0;

    sort(clusterIndex.begin(),clusterIndex.end());
    clusterIndex.erase(unique(clusterIndex.begin(),clusterIndex.end()),clusterIndex.end());
    nClus = clusterIndex.size();  
    return nClus; 
  }

  void EraseElementInPoint(int i){
     phi.erase(phi.begin()+i);
     clusterIndex.erase(clusterIndex.begin()+i);
     weight.erase(weight.begin()+i);
     if(phi.size()<1) n = 0;
      n = phi.size();
     // std::cout<<"erase n: "<<n<<endl;
  }

/*  void ToClusters(vector<Cluster> vec_clus){
     vec_clus.clear();
     int N = NumberofCluster();
     for(int index=0;index<N;index++){
       Cluster nclus;
       for(int i=0;i<n;i++){
         if(clusterIndex[i] == index ){
	   Bar nbar(x[i],y[i],clusterIndex[i],weight[i]);
	   nclus.AddBar(nbar);	
         }
       }
       vec_clus.push_back(nclus);
     }
  }*/

};
#endif
