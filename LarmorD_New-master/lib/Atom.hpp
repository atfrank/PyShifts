/* 

  Copyright University of Michigan.
  This file is part of the Larmor software suite and is made available under license.
  University of Michigan (UM) TECHtransfer: phone: 734-763-0614 email: techtransfer@umich.edu.
  
  Author: Sean M. Law and Aaron T. Frank
     
*/

#ifndef ATOM_H
#define ATOM_H

#include "Coor.hpp"

#include <vector>
#include <string>

//Forward declaration
class Chain;
class Residue;

class Atom {
  private:
    std::string tag;
    std::string recname; //Record name: "ATOM  ", "HETATM"
    int  atmnum; //Atom serial number
    std::string atmname; //Atom name
    std::string atmtype; //Mol2
    std::string alt; //Alternate location indicator
    Residue* res;
    std::string resname; //Residue name
    Chain* chn;
    std::string chainid; //Chain identifier, modified if realid is duplicated
    std::string realid; //Store original chainid
    int  resid; //Residue sequence number
    std::string icode; //Code for insertion of residues
    Coor coor; //X, Y, Z Coordinates
    double occu; //Occupancy
    double bfac; //B-factor or Temperature factor
    std::string segid; //Segment identifier
    bool sel; //Selection flag
    std::string summary; //A:GLY1.CA style summary
    std::string ss; //Secondary structure
    double mass;
    double charge;
    std::string elem;
    unsigned int atminx; //For easier 2-D lookup tables
    std::vector<double> data;
    std::vector<Atom*> bonds;
    //All additional fields must also be added to Atom::clone(), Atom::dummy() functions!!

  public:
    Atom(); //Constructor
    Atom(int atmnum, std::string atmname, std::string resname, int resnum, Coor vec, std::string seg=0); //Overload constructor

    void reset();
    void clone(Atom* atmin);
    void dummy();

    //Get atom info
    std::string& getTag();
    std::string& getRecName();
    int& getAtmNum();
    std::string& getAtmName();
    std::string& getAtmType();
    std::string& getAlt();
    Residue* getResidue();
    std::string& getResName();
    Chain* getChain();
    std::string& getChainId();
    std::string& getRealId();
    int& getResId();
    std::string& getICode();
    Coor& getCoor ();
    double& getX();
    double& getY();
    double& getZ();
    double& getOccu();
    double& getBFac();
    std::string& getElem();
    std::string& getSegId();
    bool& getSel();
    std::string& getSummary();
    std::string& getSS();
    double& getMass();
    double& getCharge();
    std::vector<double>& getData();
    double& getDataPoint(const unsigned int element);
    unsigned int& getAtmInx();
    unsigned int getDataSize();
    unsigned int getBondsSize();
    std::vector<Atom*>& getBonds();
    Atom* getBond(const unsigned int &element);
  
    //Set atom info
    void setTag(const std::string& tagin);
    void setRecName(const std::string& recnamein);
    void setAtmNum(const int& atmnumin);
    void setAtmName(const std::string& atmnamein);
    void setAtmName(); //Clear
    void setAtmType(const std::string& atmtypein);
    void setAtmType(); //Clear
    void setAlt(const std::string& altin);
    void setAlt(); //Clear
    void setResidue(); //Clear
    void setResidue(Residue* resin=NULL);
    void setResName(const std::string& resnamein);
    void setResName(); //Clear
    void setChain(Chain* chnin=NULL);
    void setChainId(const std::string& chainidin);
    void setRealId(const std::string& realidin);
    void setChainId(); //Clear
    void setRealId(); //Clear
    void setResId(const int& residin);
    void setICode(const std::string& icodein);
    void setICode(); //Clear
    void setX(const double &xin);
    void setY(const double &yin);
    void setZ(const double &zin);
    void setCoor(const Coor& coorin);
    void setOccu(const double& occuin);
    void setBFac(const double& bfacin);
    void setElem(const std::string& elemin);
    void setSegId(const std::string& segidin);
    void setSegId(); //Clear
    void setSel(const bool selin);
    void setSummary(const std::string& summaryin);
    void setSS(const std::string& ssin);
    void setMass(const double& massin);
    void setCharge(const double& chargein);
    void setAtmInx(const unsigned int& atminxin);
    void addData(const double& din);
    void clearData();
    void addBond(Atom* atmin);
    void clearBonds();
};

#endif
