/* 

  Copyright University of Michigan.
  This file is part of the Larmor software suite and is made available under license.
  University of Michigan (UM) TECHtransfer: phone: 734-763-0614 email: techtransfer@umich.edu.
  
  Author: Sean M. Law and Aaron T. Frank
     
*/

#include "PDB.hpp"

#include "Molecule.hpp"
#include "Chain.hpp"
#include "Residue.hpp"
#include "Atom.hpp"
#include "Misc.hpp"

#include <fstream>
#include <iomanip>
#include <cstdlib>

PDB::PDB(){
  chnMap.clear();
  format.clear();
}

void PDB::writePDBFormat (Molecule* mol, std::ostringstream &out, bool selFlag, bool chnFlag){
  Chain *chn;
  Residue *res;
  Atom *atm;
  int natom=0;
  int catom=0;
  unsigned int i;
  std::map<char, int> mapIds;
  char currId;
  bool addIdFlag;

  out.clear();
  currId='A';
  addIdFlag=false;

  //Create map of chainIds
  for (i=0; i< mol->getChnVecSize(); i++){
    chn=mol->getChain(i);
    if (chn->getChainId().compare(" ") != 0){
      mapIds[chn->getChainId().at(0)]=1; 
    }
  }

  for (i=0; i< mol->getChnVecSize(); i++){
    chn=mol->getChain(i);
    catom=0;
    //if(!chn->getSel()){continue;}
    for (unsigned int j=0; j< chn->getResVecSize(); j++){
      res=chn->getResidue(j);
      //if(!res->getSel()){continue;}
      for (unsigned int k=0; k< res->getAtmVecSize(); k++){
        atm=res->getAtom(k);
        if (selFlag == true && atm->getSel() == false){
          continue;
        }
        if (atm->getAtmNum() <= 99999){
          out << std::setw(6) << std::left << atm->getRecName();
          out << std::setw(5) << std::right << atm->getAtmNum(); //PDB format
          out << " "; //PDB format
        }
        else{
          out << std::setw(5) << std::left << atm->getRecName();
          out << std::setw(6) << std::right << atm->getAtmNum();
          out << " ";
        }
        out << std::setw(4) << std::left << atm->getAtmName();
        out << std::setw(1) << std::left << atm->getAlt();
        if (atm->getResName().length() < 4){
          out << std::setw(3) << std::left << atm->getResName(); //PDB format
          out << " "; //PDB format
        }
        else{
          out << std::setw(4) << std::left << atm->getResName();
        }
        if (chnFlag == true && atm->getChainId().compare(" ") == 0){
          while (mapIds.find(currId) != mapIds.end()){
            currId++;
          }
          out << std::setw(1) << std::left << currId;
          addIdFlag=true;
        }
        else{
          out << std::setw(1) << std::left << atm->getChainId();
        }
        if (atm->getResId() <= 999){
          out << std::setw(4) << std::right << atm->getResId(); //PDB format
          out << std::setw(1) << std::left << atm->getICode(); //PDB format
          out << "   "; //PDB format
        }
//        else if (atm->getResId() <= 9999){
//          out << std::setw(5) << std::right << atm->getResId();
//          out << "   ";
//        }
        else if (atm->getResId() <= 99999){
          out << std::setw(5) << std::right << atm->getResId();
          out << "   ";
        }
        else{
          out << std::setw(6) << std::left << atm->getResId();
        }
        out << std::fixed; //For setting precision
        out << std::setw(8) << std::right << std::setprecision(3) << atm->getX();
        out << std::setw(8) << std::right << std::setprecision(3) << atm->getY();
        out << std::setw(8) << std::right << std::setprecision(3) << atm->getZ();
        out << std::setw(6) << std::right << std::setprecision(2) << atm->getOccu();
        out << std::setw(6) << std::right << std::setprecision(2) << atm->getBFac();
        out << "      ";
        out << std::setw(4) << std::left << atm->getSegId();
        out << " ";
        out << std::setw(1) << std::left << atm->getElem();
        out << std::endl;
        natom++;
        catom++;
      }
    }
    if (catom>0){
      out << "TER" << std::endl;
      if (addIdFlag == true){
        mapIds[currId]=1;
        addIdFlag=false;
      }
    }
  }

  if (natom > 0){
    out << "END" << std::endl;
  }

}

Molecule* PDB::readPDB(const std::string ifile, const int model, const std::string format, const bool hetFlag, const bool remFlag){
  std::ifstream pdbFile;
  std::istream* inp;
  std::string line;
  int currModel=0;
  bool modelFlag=true;
  bool firstModelFlag=false;
  Molecule *mol;
  Chain *chnEntry=new Chain;
  Residue *resEntry=new Residue;
  Atom *atmEntry; //Created in heap by processAtomLine
  Atom *lastAtom;
  PDB pdb;
  std::string PDBID;
  unsigned int year;

  if (format.compare("CHARMM") == 0){
    mol=new MoleculeCHARMM;
  }
  else{
    mol=new Molecule;
  }

  atmEntry=NULL;
  lastAtom=NULL;
  mol->setCopyFlag(false);

  if (ifile.compare("-") == 0){ //Input from pipe
    inp=&std::cin;
  }
  else{ //Input from file
    pdbFile.open((ifile).c_str(), std::ios::in);
    inp=&pdbFile;
    PDBID=Misc::toupper(ifile.substr(0,4));
    mol->setTag(PDBID);
  }

  while (inp->good() && !(inp->eof())){

    getline(*inp,line);
    if (line.compare(0,6,"HEADER") == 0){
      std::stringstream(line.substr(57,2)) >> year;
      if (year < 71){
        mol->setYear(year+2000);
      }
      else{
        mol->setYear(year+1900);
      }
      if (Misc::trim(line.substr(62,4)).length() == 4){
        PDBID=Misc::toupper(Misc::trim(line.substr(62,4)));
      }
    }
    else if (line.compare(0,6,"EXPDTA") == 0){
      mol->setExp(Misc::trim(line.substr(10,60)));
    }
    else if (line.size() > 6 && line.compare(0,6,"MODEL ") == 0){
      std::stringstream(line.substr(10,4)) >> currModel;
      if ((model==0 && firstModelFlag==false) || currModel==model){
        firstModelFlag=true; //Processed first model
        modelFlag=true;
      }
      else{
        modelFlag=false;
      }
    }
    else if (line.size() > 6 && line.compare(0,6,"REMARK") == 0 && remFlag == true){
      mol->addRemark(line);
    }
    else if (modelFlag && line.size() >= 54 && (line.compare(0,4,"ATOM") == 0 || (line.compare(0,6,"HETATM") == 0 && hetFlag == true) || (line.compare(0,5,"HETAT") == 0 && hetFlag == true))){
      //Atom
      atmEntry=pdb.processAtomLine(line, lastAtom);
      atmEntry->setTag(PDBID);
      if (lastAtom != NULL && lastAtom->getResId() == atmEntry->getResId()){
        if (lastAtom->getResName().compare(atmEntry->getResName()) == 0 && (lastAtom->getAlt().compare(0,1," ") == 0 || lastAtom->getAlt().compare(0,1,atmEntry->getAlt(),0,1) == 0)){
          //Eliminate other alternate locations besides first
          mol->addAtom(atmEntry);
        }
        else{
          continue;
        }
      }
      else{
        mol->addAtom(atmEntry);
      }

      //Residue/Chain
      if (lastAtom != NULL && atmEntry->getChainId().compare(lastAtom->getChainId()) != 0) {
        //Store last
        chnEntry->addResidue(resEntry);
        mol->addResidue(resEntry);
        mol->addChain(chnEntry);
        //Create new
        chnEntry=new Chain;
        resEntry=new Residue;
      }
      else if (lastAtom != NULL && lastAtom->getResId() != atmEntry->getResId()){
        chnEntry->addResidue(resEntry);
        mol->addResidue(resEntry);
        resEntry=new Residue;
      }
      else{
        //Do nothing
      }
      
      atmEntry->setResidue(resEntry);
      atmEntry->setChain(chnEntry);
      resEntry->addAtom(atmEntry);
      chnEntry->addAtom(atmEntry);

      //Update for next atom
      lastAtom=atmEntry;
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
  
  if (ifile.compare("-") != 0 && pdbFile.is_open()){
    pdbFile.close();
  }

  if (mol->getAtmVecSize() == 0){
    std::cerr << std::endl << "Error: The PDB file \"" << ifile << "\" ";
    std::cerr << "did not contain any valid atoms" << std::endl << std::endl;
    exit(1);
  }

  mol->setICodeFlag(mol->checkICode());

  return mol;
}

Atom* PDB::processAtomLine (std::string line, Atom* lastAtom){
  double x,y,z;
  std::string recname; //Record name: "ATOM  ", "HETATM"
  int  atmnum; //Atom serial number
  std::string chainid;
  int  resid; //Residue sequence number
  double occu; //Occupancy
  double bfac; //B-factor or Temperature factor
  std::string segid; //Segment identifier
  Atom *atmEntry=new Atom;

  //substr: first character is denoted by a value of 0 (not 1)
  if (Misc::isdigit(Misc::trim(line.substr(4,7)))){
    atmEntry->setRecName(line.substr(0,4));
    std::stringstream(line.substr(4,7)) >> atmnum;
    atmEntry->setAtmNum(atmnum);
  }
  else{
    atmEntry->setRecName(line.substr(0,6)); //PDB format
    std::stringstream(line.substr(6,5)) >> atmnum; //PDB format
    atmEntry->setAtmNum(atmnum); //PDB format
  }
  atmEntry->setAtmName(Misc::trim(line.substr(12,4)));
  atmEntry->setAlt(line.substr(16,1));
  //atmEntry->setResName(line.substr(17,3)); //PDB format
  atmEntry->setResName(Misc::trim(line.substr(17,4)));
  chainid=line.substr(21,1);
  if(lastAtom == NULL || lastAtom->getChainId().compare(chainid) != 0){
    //New chain, check if duplicate
    if (chnMap.find(chainid) == chnMap.end()){
      chnMap[chainid]=chnMap.size();
    }
    else{
      //Duplicate, likely HETATM ligand bound to same chain
      //Assign numeric chainid, this is valid PDB format
      //std::ostringstream ossChainId;
      //ossChainId<< chnMap.find(chainid)->second;
      //chainid=ossChainId.str();
      //Set chainid to blank
      chainid=" ";
    }
  }
  atmEntry->setChainId(chainid);
  atmEntry->setRealId(line.substr(21,1)); 
  //std::stringstream(line.substr(22,4)) >> resid; //PDB format
  std::stringstream(line.substr(22,6)) >> resid; //Automagically truncates anything that is not numeric, i.e. iCode
  atmEntry->setResId(resid);
  atmEntry->setICode(line.substr(26,1));
  std::stringstream(line.substr(30,8)) >> x;
  std::stringstream(line.substr(38,8)) >> y;
  std::stringstream(line.substr(46,8)) >> z;
  atmEntry->setCoor(Coor(x,y,z));
  if (line.size() >= 60){
    std::stringstream(line.substr(54,6)) >> occu;
    atmEntry->setOccu(occu);
  }
  if (line.size() >= 66){
    std::stringstream(line.substr(60,6)) >> bfac;
    atmEntry->setBFac(bfac);
  }
  if (line.size() >= 76){
    segid=line.substr(72,4);
    atmEntry->setSegId(segid);
  }
  if (line.size() >= 78){
    if (line.substr(77,1).compare(" ") == 0){
      atmEntry->setElem(atmEntry->getAtmName().substr(0,1));
    }
    else{
      atmEntry->setElem(line.substr(77,1));
    }
  }
  atmEntry->setSel(true);
  std::stringstream ss;
  ss << resid;
  atmEntry->setSummary(chainid+":"+atmEntry->getResName()+ss.str()+Misc::trim(atmEntry->getICode())+"."+Misc::trim(atmEntry->getAtmName()));

  return atmEntry;
}

std::string PDB::formatCHARMMResName (Atom* atmEntry){
  if (atmEntry->getResName().compare("HIE") == 0){
    return "HSE";
  }
  else if (atmEntry->getResName().compare("HID") == 0){
    return "HSD";
  }
  else if (atmEntry->getResName().compare("HIP") == 0){
    return "HSP";
  }
  else if (atmEntry->getResName().compare("AHE") == 0){
    return "CT2";
  }
  else if (atmEntry->getResName().compare("NME") == 0){
    return "CT3";
  }
  else if (atmEntry->getResName().compare("CYX") == 0){
    std::cerr << "Warning: " << atmEntry->getSummary() << " has a disulfide bond" << std::endl;
    return "CYS";
  }
  else if (atmEntry->getResName().compare("FOR") == 0 || atmEntry->getResName().compare("CSO") == 0 || atmEntry->getResName().compare("CME") == 0 ){
    std::cerr << "Warning: " << atmEntry->getSummary() << " has no matching residue name in CHARMM" << std::endl;
    return atmEntry->getResName();
  }
  else{
    return atmEntry->getResName();
  }
}

int PDB::formatCHARMMResId(Atom* atmEntry, Residue* lastRes, Residue* nextRes){
  if (atmEntry->getResName().compare("ACE") == 0 && nextRes != NULL){
    return nextRes->getResId();
  }
  else if ((atmEntry->getResName().compare("AHE") == 0 || atmEntry->getResName().compare("NME") == 0) && lastRes != NULL){
    return lastRes->getResId();
  }
  else{
    return atmEntry->getResId();
  }
}
