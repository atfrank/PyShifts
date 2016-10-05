/* 

  Copyright University of Michigan.
  This file is part of the Larmor software suite and is made available under license.
  University of Michigan (UM) TECHtransfer: phone: 734-763-0614 email: techtransfer@umich.edu.
  
  Author: Sean M. Law and Aaron T. Frank
     
*/

#ifndef CHAIN_H
#define CHAIN_H

#include <vector>
#include <string>

//Forward Declaration
class Residue;
class Atom;

class Chain {
  private:
    std::vector<Residue *> resVec; //Coor of residue pointers
    std::vector<Atom *> atmVec; //Coor of atom pointers
//    bool sel;

  public:
    Chain();

    void reset();

    void addResidue(Residue* resEntry);
    void addAtom(Atom* atmEntry);

    //Get Chain info
    Atom* getAtom(const unsigned int& element);
    Residue* getResidue(const unsigned int& element);
    std::string getChainId();
    unsigned int getAtmVecSize();
    unsigned int getResVecSize();
//    void setSel(bool selin);
//    bool& getSel();
    void selAll();
    void deselAll();
};

#endif
