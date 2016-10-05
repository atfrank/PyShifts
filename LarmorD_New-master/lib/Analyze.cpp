/* 

  Copyright University of Michigan.
  This file is part of the Larmor software suite and is made available under license.
  University of Michigan (UM) TECHtransfer: phone: 734-763-0614 email: techtransfer@umich.edu.
  
  Author: Sean M. Law and Aaron T. Frank
     
*/

#include "Analyze.hpp"

#include "Molecule.hpp"
#include "Chain.hpp"
#include "Residue.hpp"
#include "Atom.hpp"
#include "Coor.hpp"
#include "Constants.hpp"
#include "Misc.hpp"

#include <fstream>
#include <iomanip>

Analyze::Analyze (){
  sel.clear();
  mol.clear();
  resel=false;
  ifile.clear();
  ofile.clear();
  verbose=false;
}

Analyze::~Analyze (){
  //Do nothing
}


void Analyze::addSel(const std::string& selin){
  this->sel.push_back(selin);
}

std::string Analyze::getSel(const int& element){
  return this->sel.at(element);
}

unsigned int Analyze::getNSel(){
  return this->sel.size();
}

void Analyze::addMol(Molecule* molin){
  this->mol.push_back(molin);
}

void Analyze::setMol(const int& element, Molecule* molin){
  this->mol.at(element)=molin;
}

void Analyze::clearMol(){
  this->mol.clear();
}

void Analyze::resizeNMol(const int sizein){
  this->mol.resize(sizein);
}

void Analyze::setupMolSel(Molecule* molin){
  Molecule* tmpmol;

  molin->select(this->getSel(0));
  tmpmol=molin->copy();
  this->addMol(tmpmol);
}

void AnalyzeDistance::setupMolSel(Molecule* molin){
  unsigned int i;
  Molecule* tmpmol;

  for (i=0; i< this->getNSel(); i++){
    molin->select(this->getSel(i));
    tmpmol=molin->copy();
    this->addMol(tmpmol);
  }
}

Molecule* Analyze::getMol(const int& element){
  return this->mol.at(element);
}

unsigned int Analyze::getNMol(){
  return this->mol.size();
}


void Analyze::readTopology(Molecule* molin, std::string topin){
  if (topin.length() > 0){
    molin->readTopology(topin);
    molin->setMass();
    molin->setCharge();
  }
}


//All preAnalysis Functions
void Analyze::preAnalysis(){
  //Do nothing
}

void Analyze::preAnalysis(Molecule* molin, std::string topin){
  this->readTopology(molin, topin);
  this->setupMolSel(molin);
}

//All postAnalysis functions

void Analyze::postAnalysis(){
  //Do nothing
}


//Basic analysis functions
Coor Analyze::centerOfGeometry(Molecule *mol, bool selFlag){
  Coor cog=Coor(0.0, 0.0, 0.0);

  for (unsigned int i=0; i< mol->getAtmVecSize(); i++){
    if (selFlag == true && mol->getAtom(i)->getSel() == false){
      continue;
    }
    cog+=mol->getAtom(i)->getCoor();
  }

  if (selFlag == true){
    cog/=mol->getNAtomSelected();
  }
  else{
    cog/=mol->getNAtom();
  }

  return cog;
}

double Analyze::distance (const Coor& u, const Coor& v){
  Coor d=u-v;
  return d.norm();
}

double Analyze::rmsd (Molecule *cmpmol, Molecule *refmol){
  double RMSD;
  double E;
  unsigned int i, j;
  Atom *atm;
  std::vector<double> dx;
  std::vector<double> dy;
  std::vector<double> dz;

  //For optimal efficiency, atoms need to be pre-selected
  //for both *cmpmol and *refmol before rmsd() is called
  RMSD=-1.0;
  E=0.0;

  //Check selection sizes and resize matrices
  if (cmpmol->getNAtomSelected() != refmol->getNAtomSelected()){
    std::cerr << std::endl << "Error: Atom number mismatch in RMSD calculation" << std::endl;
    std::cerr << "CMP-NATOM: " << cmpmol->getNAtomSelected() << ", ";
    std::cerr << "REF-NATOM: " << refmol->getNAtomSelected() << std::endl;
    return -1.0;
  }

  for (i=0; i< cmpmol->getNAtom(); i++){
    atm=cmpmol->getAtom(i);
    if (atm->getSel() == false){
      continue;
    }
    dx.push_back(atm->getX());
    dy.push_back(atm->getY());
    dz.push_back(atm->getZ());
  }

  j=0;
  for (i=0; i< refmol->getNAtom(); i++){
    atm=refmol->getAtom(i);
    if (atm->getSel() == false){
      continue;
    }
    dx.at(j)=dx.at(j)-atm->getX();
    dy.at(j)=dy.at(j)-atm->getY();
    dz.at(j)=dz.at(j)-atm->getZ();
    E+=dx.at(j)*dx.at(j)+dy.at(j)*dy.at(j)+dz.at(j)*dz.at(j);
    j++;
  }

  RMSD=sqrt(E/cmpmol->getNAtomSelected());

  return RMSD;

}

void Analyze::pairwiseDistance(Molecule *mol, std::vector<std::vector<double> >& pdin){
  std::vector<Atom*>::iterator ai;
  std::vector<Atom*>::iterator aj;
  unsigned int natom;
  unsigned int aiInx;
  bool flag;

  natom=mol->getAtmVecSize();

  pdin.clear();
  pdin.resize(natom);

  for (ai=mol->getAtmVec().begin(); ai != mol->getAtmVec().end(); ++ai){
    aiInx=(*ai)->getAtmInx();
    pdin.at(aiInx).resize(natom);
    pdin.at(aiInx).at(aiInx)=0.0; //Zero diagonal
  }

  for (ai=mol->getAtmVec().begin(); ai != mol->getAtmVec().end(); ++ai){
    aiInx=(*ai)->getAtmInx();
    if ((*ai)->getX() < 9999.9){
      flag=true;
    }
    else{
      flag=false;
    }
    //Lower Triangle
    for (aj=mol->getAtmVec().begin(); aj != ai; ++aj){
      if (flag && (*aj)->getX() < 9999.9){
        pdin.at(aiInx).at((*aj)->getAtmInx())=Analyze::distance((*ai)->getCoor(), (*aj)->getCoor());
      }
      else{
        pdin.at(aiInx).at((*aj)->getAtmInx())=9999.9;
      }
    }
  
    //Upper Triangle
    for (aj=ai+1; aj != mol->getAtmVec().end(); ++aj){
      if (flag && (*aj)->getX() < 9999.9){
        pdin.at(aiInx).at((*aj)->getAtmInx())=Analyze::distance((*ai)->getCoor(), (*aj)->getCoor());
      }
      else{
        pdin.at(aiInx).at((*aj)->getAtmInx())=9999.9;
      }
    }
  }
}

//All runAnalysis Functions

void AnalyzeDistance::runAnalysis(){
  std::cout << std::fixed;
  std::cout << std::setw(9) << std::right << std::setprecision(3) << Analyze::distance(Analyze::centerOfGeometry(this->getMol(0)), Analyze::centerOfGeometry(this->getMol(1)));
}
