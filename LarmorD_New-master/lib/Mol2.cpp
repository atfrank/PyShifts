/* 

  Copyright University of Michigan.
  This file is part of the Larmor software suite and is made available under license.
  University of Michigan (UM) TECHtransfer: phone: 734-763-0614 email: techtransfer@umich.edu.
  
  Author: Sean M. Law and Aaron T. Frank
     
*/


#include "Mol2.hpp"

#include "Molecule.hpp"
#include "Chain.hpp"
#include "Residue.hpp"
#include "Atom.hpp"
#include "Misc.hpp"
#include "Select.hpp"

#include <fstream>
#include <cstdlib>
Mol2::Mol2(){
  chnMap.clear();
}

Molecule* Mol2::readMol2(std::string ifile, std::string format){
  std::ifstream mol2File;
  std::istream* inp;
  std::string line;
  Molecule *mol;
  Chain *chnEntry=new Chain;
  Residue *resEntry=new Residue;
  Atom *atmEntry; //Created in heap by processAtomLine
  Atom *lastAtom;
  Mol2 mol2;
  std::map<unsigned int, Atom*> atmMap;
  std::vector<std::string> s;
  unsigned int anum1, anum2;
  bool readAtoms;
  bool readBonds;

  if (format.compare("CHARMM") == 0){
    mol=new MoleculeCHARMM;
  }
  else{
    mol=new Molecule;
  }

  atmEntry=NULL;
  lastAtom=NULL;
  mol->setCopyFlag(false);
  readAtoms=false;
  readBonds=false;

  if (ifile.compare("-") == 0){ //Input from pipe
    inp=&std::cin;
  }
  else{ //Input from file
    mol2File.open((ifile).c_str());
    inp=&mol2File;
    mol->setTag(ifile.substr(0,4));
  }

  while (inp->good()){

    getline(*inp,line);
    if (line.length() == 0){
      continue;
    }
    
    if (readAtoms == true && line.size() >= 44 ){
      //Atom
      atmEntry=mol2.processAtomLine(line, lastAtom);
      if (atmEntry == NULL){
        continue;
      }
      mol->addAtom(atmEntry);
      atmMap.insert(std::pair<unsigned int, Atom*>(atmEntry->getAtmNum(), atmEntry));

      //Residue/Chain
      if (lastAtom != NULL){
        if(atmEntry->getChainId().compare(lastAtom->getChainId()) != 0 || atmEntry->getResId() < lastAtom->getResId()){
          //Store last
          chnEntry->addResidue(resEntry);
          mol->addResidue(resEntry);
          mol->addChain(chnEntry);
          //Create new
          chnEntry=new Chain;
          resEntry=new Residue;
        }
        else if (lastAtom->getResId() != atmEntry->getResId()){
          chnEntry->addResidue(resEntry);
          mol->addResidue(resEntry);
          resEntry=new Residue;
          if (lastAtom->getResName().compare("NME") == 0 || lastAtom->getResName().compare("AHE") == 0 || atmEntry->getResName().compare("UNK") == 0){
            mol->addChain(chnEntry);
            chnEntry=new Chain;
          }
        }
        else{
          //Do Nothing
        }
      }
      else{
        //Do nothing
      }

      resEntry->addAtom(atmEntry);
      chnEntry->addAtom(atmEntry);

      //Update for next atom
      lastAtom=atmEntry;
    }
    else if (readBonds == true){
      Misc::splitStr(Misc::trim(line), " \t", s, false);
      if (s.size() >= 3){
        std::stringstream(s.at(1)) >> anum1;
        std::stringstream(s.at(2)) >> anum2;
        if (atmMap.find(anum1) != atmMap.end() && atmMap.find(anum2) != atmMap.end()){
          atmMap.at(anum1)->addBond(atmMap.at(anum2));
          atmMap.at(anum2)->addBond(atmMap.at(anum1));
        }
      }
    }
    else if (line.compare(0,9,"@<TRIPOS>") == 0){
      readAtoms=false;
      readBonds=false;
      if (line.compare(9,4,"ATOM") == 0){
        readAtoms=true;
      }
      else if (line.compare(9,4,"BOND") == 0){
        readBonds=true;
      }
      else{
        //Do Nothing
      }
    }
    else{
      continue;
    }
  }
  //Add remaining residues and chains
  if (resEntry->getAtmVecSize() > 0){
    chnEntry->addResidue(resEntry);
    mol->addResidue(resEntry);
    mol->addChain(chnEntry);
  }
  
  if (ifile.compare("-") != 0){
    mol2File.close();
  }

  if (mol->getAtmVecSize() == 0){
    std::cerr << std::endl << "Error: The Mol2 file \"" << ifile << "\" ";
    std::cerr << "did not contain any valid atoms" << std::endl << std::endl;
    exit(1);
  }

  mol->select("!:PROTEIN+NUCLEIC+ION+METAL+SOLVENT+TERMINI.", false, false);
  //Rename unknown residues to "UNK"
  for (unsigned int j=0; j< mol->getNAtom(); j++){
    if (mol->getAtom(j)->getSel() == true){
      mol->getAtom(j)->setResName("UNK");
    }
  }

  mol->selAll();

  return mol;
}

Atom* Mol2::processAtomLine (std::string line, Atom* lastAtom){
  double x,y,z;
  std::string recname; //Record name: "ATOM  ", "HETATM"
  int  atmnum; //Atom serial number
  std::string chainid;
  int  resid; //Residue sequence number
  std::string segid; //Segment identifier
  Atom *atmEntry=new Atom;
  std::vector<std::string> s;
  Select *sel=new Select;
  Molecule *mol=new Molecule;
  std::vector<Atom *> atmSel;
  double charge;

  sel->initKeys(mol); //mol is not actually used, needed to initialize selection keys

  Misc::splitStr(Misc::trim(line), " \t", s, false);
  std::stringstream(s.at(0)) >> atmnum;
  atmEntry->setAtmNum(atmnum);
  atmEntry->setAtmName(Misc::trim(s.at(1)));
  std::stringstream(s.at(2)) >> x;
  std::stringstream(s.at(3)) >> y;
  std::stringstream(s.at(4)) >> z;
  atmEntry->setCoor(Coor(x,y,z));

  atmEntry->setAtmType(Misc::trim(s.at(5)));

  if (s.size() >= 6){
    std::stringstream(s.at(6)) >> resid;
    atmEntry->setResId(resid);
  }
  if (s.size() >= 7){
    if (s.at(7).length() >= 3){
      atmEntry->setResName(Misc::trim(s.at(7).substr(0,3)));
    }
  }
  if (s.size() >= 8){
    std::stringstream(s.at(8)) >> charge;
    atmEntry->setCharge(charge);
  }

  atmEntry->setChainId(" ");

  std::stringstream() >> resid;
  atmEntry->setResId(resid);

  atmEntry->setSegId("");
  atmEntry->setSel(true);

  std::stringstream ss;
  ss << resid;

  atmEntry->setSummary(chainid+":"+atmEntry->getResName()+ss.str()+"."+Misc::trim(atmEntry->getAtmName()));

  return atmEntry;
}
