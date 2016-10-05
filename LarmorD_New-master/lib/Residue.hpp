/* 

  Copyright University of Michigan.
  This file is part of the Larmor software suite and is made available under license.
  University of Michigan (UM) TECHtransfer: phone: 734-763-0614 email: techtransfer@umich.edu.
  
  Author: Sean M. Law and Aaron T. Frank
     
*/

#ifndef RESIDUE_H
#define RESIDUE_H

#include <vector>
#include <string>

//Forward Declaration
class Atom;

class Residue {
  private:
    std::vector<Atom*> atmVec;
//    bool sel;

  public:
    Residue();

    void reset();
    int getResId();
    std::string getResName();
    std::string getChainId();
    Atom* getStart();
    Atom* getEnd();
    std::string getSegId();
    void addAtom(Atom* atmEntry);
    std::vector<Atom*>& getAtmVec();
    Atom* getAtom (const unsigned int &element);
    unsigned int getAtmVecSize();
//    void setSel(bool selin);
//    bool& getSel();
    void selAll();
    void deselAll();

    static std::string aa321(const std::string &aa);
    static std::string aa123(const std::string &aa);
};

#endif
