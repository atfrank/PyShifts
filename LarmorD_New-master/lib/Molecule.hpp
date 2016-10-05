/* 

  Copyright University of Michigan.
  This file is part of the Larmor software suite and is made available under license.
  University of Michigan (UM) TECHtransfer: phone: 734-763-0614 email: techtransfer@umich.edu.
  
  Author: Sean M. Law and Aaron T. Frank
     
*/

#ifndef MOLECULE_H
#define MOLECULE_H

#include "Prmtop.hpp"
#include "Coor.hpp"
#include "Enum.hpp"

#include <vector>
#include <map>

//Forward Declaration
class Chain;
class Residue;
class Atom;

//Base class
class Molecule {
  private:
    std::vector<Chain*> chnVec;
    std::vector<Residue*> resVec;
    std::vector<Atom*> atmVec;
    bool copyFlag; //This molecule is a copy if true
    std::map< std::string, std::vector<bool> > storedSel;
    std::string remarks;
    bool iCodeFlag;
    Prmtop toppar;
    unsigned int year;
    std::string exp;
    std::string tag;

  public:
    Molecule(); //Constructor
    virtual ~Molecule();
    static Molecule* readPDB (const std::string ifile, const int model=0, const std::string format="", const bool hetFlag=true);
    static Molecule* readPDB (const std::string ifile, const std::string format, const bool hetFlag=true);
//    std::string writePDB (bool selFlag=true, bool print=true, bool chnFlag=false);
    std::string writePDB (bool selFlag, bool print, bool chnFlag);
    std::string writePDB (bool selFlag, bool print);
    std::string writePDB (bool chnFlag);
    std::string writePDB ();
    static Molecule* readMol2 (const std::string ifile, const std::string format="", const bool stopFlag=false);
    Molecule* clone(bool selFlag=true, bool keep=true);
    Molecule* copy(bool selFlag=true);
    void cat (Molecule* catmol, bool selFlag=true, bool keep=true);
    void addAtom(Atom* atmEntry);
    Atom* getAtom(const unsigned int& element);
    void renameAtom(const std::string &search, const std::string &replace);
    void renameAtom(const std::vector<std::string> &search, const std::string &replace); 
    void addChain(Chain* chnEntry);
    void addMissingChainIds();
    void addResidue(Residue* resEntry);
    Chain* getChain(const unsigned int& element);
    unsigned int getChnVecSize();
    unsigned int getResVecSize();
    unsigned int getAtmVecSize();
    std::vector<Atom*>& getAtmVec();
    std::vector<Atom*> getAtmVecClone();
    Residue* getResidue(const unsigned int& element);
    void renameRes(const std::string &search, const std::string &replace);
    void renameRes(const std::vector<std::string> &search, const std::string &replace);
    void renameRes(int resid, const std::string &search, const std::string &replace);
    void readTopology(const std::string& topin);
    void readParameter(const std::string& prmin);
    void setMass();
    void setCharge();
    void selAll();
    void deselAll();
    void select(std::string sel, bool dieFlag=true, bool verbose=true);
    unsigned int getNAtom();
    unsigned int getNAtomSelected(); //Determining this on the fly is a good safety measure
    void setCopyFlag(bool copyFlagIn=false);
    bool getCopyFlag();
    void storeSel(std::string key="tmp");
    void recallSel(std::string key="tmp");
    void eraseSel(std::string key="tmp");
    void zeroCoor();
    void addRemark(const std::string& remin);
    void clearRemark();
    std::string getRemark();
    bool checkICode();
    void setICodeFlag(bool iCodeFlagIn=false);
    bool getICodeFlag();
    void assignAtmInx();
    void resetAtmInx();
    void setYear(const unsigned int& yearin);
    unsigned int getYear();
    void setExp(const std::string& expin);
    std::string getExp();
    void setTag(const std::string& tagin);
    std::string getTag();

    double rmsd (Molecule *refmol);
    void recenter (Molecule *recmol);
    void translate (const double &dx, const double &dy, const double &dz);
    void translate (const Coor &u);
    void rotate (const double &r1c1, const double &r1c2, const double &r1c3, const double &r2c1, const double &r2c2, const double &r2c3, const double &r3c1, const double &r3c2, const double &r3c3);
    void center (bool selFlag=true);
    void modPseudoCenter();
    void renameHis();

    //Virtual functions
    virtual void format();
};

class MoleculeCHARMM: public Molecule{
  public:
    void format();
};

#endif
