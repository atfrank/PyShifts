/* 

  Copyright University of Michigan.
  This file is part of the Larmor software suite and is made available under license.
  University of Michigan (UM) TECHtransfer: phone: 734-763-0614 email: techtransfer@umich.edu.
  
  Author: Sean M. Law and Aaron T. Frank
     
*/

#ifndef PRMTOP_H
#define PRMTOP_H

#include <string>
#include <map>

class Prmtop {
  private:
    //Maps are referenced by pair(resname, atmname)
    std::map<std::pair<std::string, std::string>, double> mass;
    std::map<std::pair<std::string, std::string>, double> charge;

  public:
    Prmtop(); //Constructor
    void readTopology(const std::string& topin);
    void readParameter(const std::string& prmin);
    double getMass(const std::string& resnamein, const std::string& atmnamein);
    double getCharge(const std::string& resnamein, const std::string& atmnamein);
};

#endif
