/* 

  Copyright University of Michigan.
  This file is part of the Larmor software suite and is made available under license.
  University of Michigan (UM) TECHtransfer: phone: 734-763-0614 email: techtransfer@umich.edu.
  
  Author: Sean M. Law and Aaron T. Frank
     
*/

#ifndef BIN_H
#define BIN_H

#include <vector>

class Bin {
  private:
    unsigned int n;
    unsigned int inx;
    std::vector<double> label;
  
  public:
    Bin();
    void setN(const unsigned int &nin);
    void setInx(const unsigned int &inxin);
    void setLabel(const std::vector<double> &labelin);
    unsigned int getN();
    unsigned int getInx();
    std::vector<double> getLabel();
    std::vector<double>& getLabelVec();
    const std::vector<double>& getLabelVec() const;
};

#endif

