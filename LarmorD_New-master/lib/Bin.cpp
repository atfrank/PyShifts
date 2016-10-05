/* 

  Copyright University of Michigan.
  This file is part of the Larmor software suite and is made available under license.
  University of Michigan (UM) TECHtransfer: phone: 734-763-0614 email: techtransfer@umich.edu.
  
  Author: Sean M. Law and Aaron T. Frank
     
*/

#include "Bin.hpp"

Bin::Bin (){
  n=0;
  inx=0;
  label.clear();
}

void Bin::setN (const unsigned int &nin){
  n=nin;
}

void Bin::setInx (const unsigned int &inxin){
  inx=inxin;
}

void Bin::setLabel (const std::vector<double> &labelin){
  label=labelin;
}

unsigned int Bin::getN (){
  return n;
}

unsigned int Bin::getInx (){
  return inx;
}

std::vector<double> Bin::getLabel (){
  return label;
}

std::vector<double>& Bin::getLabelVec (){
  return label;
}

const std::vector<double>& Bin::getLabelVec () const{
  return label;
}
