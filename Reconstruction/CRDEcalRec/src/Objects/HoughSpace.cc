#ifndef HOUGHSPACE_C
#define HOUGHSPACE_C

#include "Objects/HoughSpace.h"
#include <iostream>
namespace PandoraPlus{

  int HoughSpace::getAlphaBin(double alpha) const{ 
    if(alpha<alpha_low){
      // cout << "int HoughSpace::getAlphaBin(double alpha): wrong alpha input: " << alpha << ", return 1" << endl;
      return 1;
    }
    if(alpha>alpha_high){
      // cout << "int HoughSpace::getAlphaBin(double alpha): wrong alpha input: " << alpha << ", return Nbins_alpha" << endl;
      return Nbins_alpha;
    }

    int alpha_bin = ceil((alpha-alpha_low)/bin_width_alpha);
    return alpha_bin;
  }

  double HoughSpace::getAlphaBinCenter(int alpha_bin) const{ 
    if(alpha_bin<1){
      // cout << "double HoughSpace::getAlphaBinCenter(int alpha_bin): wrong alpha_bin input: " << alpha_bin << ", return bincenter when alpha_bin=1" << endl;
      return alpha_low + (bin_width_alpha * 0.5);
    }
    if(alpha_bin>Nbins_alpha){
      // cout << "double HoughSpace::getAlphaBinCenter(int alpha_bin): wrong alpha_bin input: " << alpha_bin << ", return bincenter when alpha_bin=Nbins_alpha" << endl;
      return alpha_high - (bin_width_alpha * 0.5);
    }

    double alpha_bin_center = alpha_low + (alpha_bin*1.*bin_width_alpha) - (bin_width_alpha*0.5) ;
    return alpha_bin_center;
  }

  double HoughSpace::getAlphaBinLowEdge(int alpha_bin) const{ 
    if(alpha_bin<1){
      // cout << "double HoughSpace::getAlphaBinLowEdge(int alpha_bin): wrong alpha_bin input: " << alpha_bin << ", return bin low edge when alpha_bin=1" << endl;
      return alpha_low;
    }
    if(alpha_bin>Nbins_alpha){
      // cout << "double HoughSpace::getAlphaBinLowEdge(int alpha_bin): wrong alpha_bin input: " << alpha_bin << ", return bin low edge when alpha_bin=Nbins_alpha" << endl;
      return alpha_high - bin_width_alpha;
    }

    double alpha_bin_low_edge = alpha_low + (alpha_bin*1.*bin_width_alpha) - bin_width_alpha;
    return alpha_bin_low_edge;
  }

  double HoughSpace::getAlphaBinUpEdge(int alpha_bin) const{
    if(alpha_bin<1){
      // cout << "double HoughSpace::getAlphaBinUpEdge(int alpha_bin): wrong alpha_bin input: " << alpha_bin << ", return bin up edge when alpha_bin=1" << endl;
      return alpha_low + bin_width_alpha;
    }
    if(alpha_bin>Nbins_alpha){
      // cout << "double HoughSpace::getAlphaBinUpEdge(int alpha_bin): wrong alpha_bin input: " << alpha_bin << ", return bin up edge when alpha_bin=Nbins_alpha" << endl;
      return alpha_high;
    }
    
    double alpha_bin_up_edge = alpha_low + (alpha_bin*1.*bin_width_alpha);
    return alpha_bin_up_edge;
  }

  int HoughSpace::getRhoBin(double rho) const{
    if(rho<rho_low){
      // cout << "int HoughSpace::getRhoBin(double rho): wrong rho input: " << rho << ", return 1" << endl;
      return 1;
    }
    if(rho>rho_high){
      // cout << "int HoughSpace::getRhoBin(double rho): wrong rho input: " << rho << ", return Nbins_rho" << endl;
      return Nbins_rho;
    }

    int rho_bin = ceil((rho-rho_low)/bin_width_rho);
    return rho_bin;
  }

  double HoughSpace::getRhoBinCenter(int rho_bin) const{
    if(rho_bin<1){
      // cout << "double HoughSpace::getRhoBinCenter(int rho_bin): wrong rho_bin input: " << rho_bin << ", return bincenter when rho_bin=1" << endl;
      return rho_low + (bin_width_rho * 0.5);
    }
    if(rho_bin>Nbins_rho){
      // cout << "double HoughSpace::getRhoBinCenter(int rho_bin): wrong rho_bin input: " << rho_bin << ", return bincenter when rho_bin=Nbins_rho" << endl;
      return rho_high - (bin_width_rho * 0.5);
    }

    double rho_bin_center = rho_low + (rho_bin*1.*bin_width_rho) - (bin_width_rho*0.5) ;
    return rho_bin_center;
  }

  double HoughSpace::getRhoBinLowEdge(int rho_bin) const{
    if(rho_bin<1){
      // cout << "double HoughSpace::getRhoBinLowEdge(int rho_bin): wrong rho_bin input: " << rho_bin << ", return bin low edge when rho_bin=1" << endl;
      return rho_low;
    }
    if(rho_bin>Nbins_rho){
      // cout << "double HoughSpace::getRhoBinLowEdge(int rho_bin): wrong rho_bin input: " << rho_bin << ", return bin low edge when rho_bin=Nbins_rho" << endl;
      return rho_high - bin_width_rho;
    }

    double rho_bin_low_edge = rho_low + (rho_bin*1.*bin_width_rho) - bin_width_rho;
    return rho_bin_low_edge;
  }

  double HoughSpace::getRhoBinUpEdge(int rho_bin) const{
    if(rho_bin<1){
      // cout << "double HoughSpace::getRhoBinUpEdge(int rho_bin): wrong rho_bin input: " << rho_bin << ", return bin up edge when rho_bin=1" << endl;
      return rho_low + bin_width_rho;
    }
    if(rho_bin>Nbins_rho){
      // cout << "double HoughSpace::getRhoBinUpEdge(int rho_bin): wrong rho_bin input: " << rho_bin << ", return bin up edge when rho_bin=Nbins_rho" << endl;
      return rho_high;
    }
    
    double rho_bin_up_edge = rho_low + (rho_bin*1.*bin_width_rho);
    return rho_bin_up_edge;
  }

  void HoughSpace::AddBinHobj(int bin_alpha, int bin_rho, int index_Hobj){
    pair<int, int> bin(bin_alpha, bin_rho);
    Hough_bins[bin].insert(index_Hobj);
  }
 


};
#endif
