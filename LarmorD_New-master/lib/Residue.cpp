/* 

  Copyright University of Michigan.
  This file is part of the Larmor software suite and is made available under license.
  University of Michigan (UM) TECHtransfer: phone: 734-763-0614 email: techtransfer@umich.edu.
  
  Author: Sean M. Law and Aaron T. Frank
     
*/

#include "Residue.hpp"

#include "Atom.hpp"

Residue::Residue (){
  atmVec.clear();
//  sel=true;
}

void Residue::reset(){
  atmVec.clear();
//  sel=true;
}

void Residue::addAtom(Atom* atmEntry){
  if (atmEntry->getAtmNum()){
    atmVec.push_back(atmEntry);
  }
}

int Residue::getResId(){
  return this->getAtom(0)->getResId();
}

std::string Residue::getResName(){
  return this->getAtom(0)->getResName();
}

std::string Residue::getChainId(){
  return this->getAtom(0)->getChainId();
}

Atom* Residue::getStart(){
  return this->getAtom(0);
}

Atom* Residue::getEnd(){
  return this->getAtom(atmVec.size()-1);
}

std::string Residue::getSegId(){
  return this->getAtom(0)->getSegId();
}

std::vector<Atom*>& Residue::getAtmVec(){
  return atmVec;
}

Atom* Residue::getAtom (const unsigned int &element){
  if (element >= atmVec.size()){
    return NULL;
  }
  else{
    return atmVec.at(element);
  }
}

unsigned int Residue::getAtmVecSize (){
  return atmVec.size();
}

//void Residue::setSel (bool selin){
//  sel=selin;
//}

//bool& Residue::getSel (){
//  return sel;
//}

void Residue::selAll(){
//  sel=true;
  for(unsigned int i=0; i< this->getAtmVecSize(); i++){
    this->getAtom(i)->setSel(true);
  }
}

void Residue::deselAll(){
//  sel=false;
  for(unsigned int i=0; i< this->getAtmVecSize(); i++){
    this->getAtom(i)->setSel(false);
  }
}

std::string Residue::aa321(const std::string &aa){

  if(aa.compare("ALA") == 0){return "A";}
  else if(aa.compare("CYS") == 0){return "C";}
  else if (aa.compare("ASP") == 0){return "D";}
  else if (aa.compare("GLU") == 0){return "E";}
  else if (aa.compare("PHE") == 0){return "F";}
  else if (aa.compare("GLY") == 0){return "G";}
  else if (aa.compare("HIS") == 0 || aa.compare("HSD") == 0 || aa.compare("HSE") == 0 || aa.compare("HSP") == 0){
    return "H";
  }
  else if (aa.compare("ILE") == 0){return "I";}
  else if (aa.compare("LYS") == 0){return "K";}
  else if (aa.compare("LEU") == 0){return "L";}
  else if (aa.compare("MET") == 0){return "M";}
  else if (aa.compare("ASN") == 0){return "N";}
  else if (aa.compare("PRO") == 0){return "P";}
  else if (aa.compare("GLN") == 0){return "Q";}
  else if (aa.compare("ARG") == 0){return "R";}
  else if (aa.compare("SER") == 0){return "S";}
  else if (aa.compare("THR") == 0){return "T";}
  else if (aa.compare("VAL") == 0){return "V";}
  else if (aa.compare("TRP") == 0){return "W";}
  else if (aa.compare("TYR") == 0){return "Y";}
  else{
    return "";
  }
}

std::string Residue::aa123(const std::string &aa){
  if(aa.compare("A") == 0){return "A";}
  else if(aa.compare("C") == 0){return "CYS";} 
  else if (aa.compare("D") == 0){return "ASP";}
  else if (aa.compare("E") == 0){return "GLY";}
  else if (aa.compare("F") == 0){return "PHE";}
  else if (aa.compare("G") == 0){return "GLY";}
  else if (aa.compare("H") == 0){return "HIS";}
  else if (aa.compare("I") == 0){return "ILE";}
  else if (aa.compare("K") == 0){return "LYS";}
  else if (aa.compare("L") == 0){return "LEU";}
  else if (aa.compare("M") == 0){return "MET";}
  else if (aa.compare("N") == 0){return "ASN";}
  else if (aa.compare("P") == 0){return "PRO";}
  else if (aa.compare("Q") == 0){return "GLN";}
  else if (aa.compare("R") == 0){return "ARG";}
  else if (aa.compare("S") == 0){return "SER";}
  else if (aa.compare("T") == 0){return "THR";}
  else if (aa.compare("V") == 0){return "VAL";}
  else if (aa.compare("W") == 0){return "TRP";}
  else if (aa.compare("Y") == 0){return "TYR";}
  else{
    return "";
  } 
}
