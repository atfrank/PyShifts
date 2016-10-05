/* 

  Copyright University of Michigan.
  This file is part of the Larmor software suite and is made available under license.
  University of Michigan (UM) TECHtransfer: phone: 734-763-0614 email: techtransfer@umich.edu.
  
  Author: Sean M. Law and Aaron T. Frank
     
*/

#include "Prmtop.hpp"

#include "Misc.hpp"

#include <fstream>

Prmtop::Prmtop (){
  mass.clear();
  charge.clear();
}

void Prmtop::readTopology(const std::string& topin){
  std::ifstream topFile;
  std::istream* topinp;
  std::string line;
  std::vector<std::string> s;
  double m; //mass
  double c; //charge
  std::map<std::string, double> massRef;
  std::string resname;
  std::string word;

  resname.clear();

  if (topin.length() > 0){
    topFile.open(topin.c_str(), std::ios::in);
    topinp=&topFile;

    while (topinp->good() && !(topinp->eof())){
      getline(*topinp, line);
      word=Misc::toupper(Misc::trim(line).substr(0,4));
      if (word.compare(0,4,"MASS") == 0){
        Misc::splitStr(line, " \t", s, false);
        if (s.size() >= 4){
          std::stringstream(s.at(3)) >> m;
          massRef.insert(std::make_pair(s.at(2), m));
        }
      }
      else if (word.compare(0,4,"RESI") == 0 || Misc::trim(line).compare(0,4,"PRES") == 0){
        Misc::splitStr(line, " \t", s, false);
        if (s.size() >= 3){
          resname=s.at(1);
        }
      }
      else if (word.compare(0,4,"ATOM") == 0){
        Misc::splitStr(line, " \t", s, false);
        if (s.size() >= 3){
          mass.insert(std::make_pair(std::make_pair(resname, s.at(1)), massRef[s.at(2)]));
          std::stringstream(s.at(3)) >> c;
          charge.insert(std::make_pair(std::make_pair(resname, s.at(1)), c));
        }
      }
      else{
        //Do nothing
      }
    }
  }

}


void Prmtop::readParameter(const std::string& prmin){

}


double Prmtop::getMass(const std::string& resnamein, const std::string& atmnamein){
  if (mass.find(std::make_pair(resnamein, atmnamein)) != mass.end()){
    return  mass[std::make_pair(resnamein, atmnamein)];
  }
  else{
    std::cerr << "Warning: Could not a find mass for atom " << resnamein;
    std::cerr << " " << atmnamein << " and was set to 1.0" << std::endl;
    return 1.0;
  }
}

double Prmtop::getCharge(const std::string& resnamein, const std::string& atmnamein){
  if (charge.find(std::make_pair(resnamein, atmnamein)) != charge.end()){
    return charge[std::make_pair(resnamein, atmnamein)];
  }
  else{
    std::cerr << "Warning: Could not find a charge for atom " << resnamein;
    std::cerr << " " << atmnamein << " and was set to 0.0." << std::endl;
    return 0.0;
  }
}

