#ifndef CALOHOUGHSPACE_H
#define CALOHOUGHSPACE_H

#include<map>
#include<set>
#include "TH2.h"

using namespace std;

namespace PandoraPlus {
  class HoughSpace{
  public:
    HoughSpace(double _alpha_low, double _alpha_high, double _bin_width_alpha, int _Nbins_alpha, double _rho_low, double _rho_high, double _bin_width_rho, int _Nbins_rho){
      alpha_low = _alpha_low;
      alpha_high = _alpha_high;
      bin_width_alpha = _bin_width_alpha;
      Nbins_alpha = _Nbins_alpha;
      rho_low = _rho_low;
      rho_high = _rho_high;
      bin_width_rho = _bin_width_rho;
      Nbins_rho = _Nbins_rho;
    };

     // Functions
    int getAlphaBin(double alpha) const;
    double getAlphaBinCenter(int alpha_bin) const;
    double getAlphaBinLowEdge(int alpha_bin) const;
    double getAlphaBinUpEdge(int alpha_bin) const;
    int getRhoBin(double rho) const;
    double getRhoBinCenter(int rho_bin) const;
    double getRhoBinLowEdge(int rho_bin) const;
    double getRhoBinUpEdge(int rho_bin) const;

    void AddBinHobj(int bin_alpha, int bin_rho, int index_Hobj);
    map< pair<int, int>, set<int> > getHoughBins() const { return Hough_bins; }

  private:
    double alpha_low;
    double alpha_high;
    double bin_width_alpha;
    int Nbins_alpha;

    double rho_low;
    double rho_high;
    double bin_width_rho;
    int Nbins_rho;

    map< pair<int, int>, set<int> > Hough_bins;
    // description of the above map:
    // key: the pair represents index of alpha bin and rho bin. All start from 1 to Nbins
    // value: a set of index of HoughObject
    // With this map, we know the relationship of bin in Hough space and Hough object

  };


};
#endif