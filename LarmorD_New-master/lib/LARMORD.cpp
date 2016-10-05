/* 

  Copyright University of Michigan.
  This file is part of the Larmor software suite and is made available under license.
  University of Michigan (UM) TECHtransfer: phone: 734-763-0614 email: techtransfer@umich.edu.
  
  Author: Sean M. Law and Aaron T. Frank
     
*/

#include "LARMORD.hpp"
#include "Molecule.hpp"
#include "Misc.hpp"
#include "Coor.hpp"
#include "Residue.hpp"
#include "Atom.hpp"
#include "Coor.hpp"
#include "Select.hpp"
#include "Analyze.hpp"

#include <fstream>
#include <stdlib.h>
#include <math.h>


LARMORD::LARMORD (Molecule *mol, const std::string fchemshift, const std::string fparmfile, const std::string freffile, const std::string faccfile, const std::string fcorfile,  bool residueBased, bool residueBasedWeights, bool mismatchCheck, bool extractor, std::string MolType, bool ringCurrent, double cutoffRing)
{
    /* nuclei for which chemical shifts will be calculated */
    this->initializeShiftAtoms();    
    /* set random coil chemical shifts */
    if (freffile.length()==0){
        this->initializeRandomShifts();
    } else {
        this->loadRefFile(freffile);
        //std::cout << "done with reference shifts..." << std::endl;      
    }
    
    this->residueBasedLarmor = residueBased;
    this->residueBasedWeightsLarmor = residueBasedWeights;
    this->mismatchCheckLarmor = mismatchCheck;
    this->MolTypeLamor = MolType;
    
    /* the accuracy weights for calculating errors */
    if (faccfile.length()==0){
        this->initializeAccuracyWeight();
    } else {
        this->loadAccFile(faccfile);
        //std::cout << "done with accuracy weights..." << std::endl;
    }
    
    if (fparmfile.length()==0){
        if (!extractor){
					std::cerr << std::endl << "Error: Please provide an input parameter file" << std::endl << std::endl;
					this->initializeAlpha();
        } 
        else {
        	this->initializeAtomTypes();
        }
    } else {
        this->loadParmFile(fparmfile);
        //std::cout << "done with parameters..." << std::endl;
    }            
    
    /* rename residue */
    if (mol != NULL){
        this->renameRes(mol);
    }    

    /* load chemical shifts from file */
    if (fchemshift.length() > 0){
        this->loadCSFile(fchemshift,mol);
        //std::cout << "done with chemical shifts..." << std::endl;
    }
    
    /* load relative correlation weighting from file */
    if (fcorfile.length() > 0){
        this->loadCorrFile(fcorfile);
    }

		/* set up ring separation */
		ringSeparation = 1.28;

		/* set cutoff */
		cutoffLarmorRing=cutoffRing;
	
		/* setUp aromatics resname keys */
		this->setUpAromaticsResKeys();

		/* setUp ringIntensity map */
		ringCurrentLarmor = ringCurrent;
		if (ringCurrentLarmor == true){
			this->setUpRingIntensityMap();		
			//std::cout << "Setting up maps.." << std::endl;
		}								                
}

std::vector<double> LARMORD::getAlpha(const std::string &key)
{
		std::vector <double> alpha;
    if (this->alphas.find (key) == this->alphas.end()){
    		alpha.push_back(0.0);
        return alpha;
    } else {
        return (this->alphas.at(key));
    }
}

std::vector<int> LARMORD::getBeta(const std::string &key)
{
    std::vector <int> beta;
    if (this->betas.find (key) == this->betas.end()){
    		beta.push_back(-3);
        return beta;        
    } else {
        return (this->betas.at(key));
    }
}

double LARMORD::getRandomShift(const std::string &key)
{
    if (this->randomShifts.find (key) == this->randomShifts.end()){
        return 0.0;
    } else {
        return (this->randomShifts.at(key));
    }
}

void LARMORD::loadCSFile(const std::string fchemshift, Molecule *mol)
{
	std::ifstream csFile;
	std::istream* csinp;
	std::string line, resname, nucleus, resid, key;
	std::vector<std::string> s;
	
	if (fchemshift.length() > 0){
		csFile.open(fchemshift.c_str(), std::ios::in);
		csinp=&csFile;
		if (csFile.is_open()) {
			while (csinp->good() && !(csinp->eof())){
				getline(*csinp, line);
				Misc::splitStr(line, " ", s, true);											
				if (s.size() == 5){
					resname = Misc::trim(s.at(0));
					resid = Misc::trim(s.at(1));
					nucleus = Misc::trim(s.at(2));
					key = resid+":"+resname+":"+nucleus;
			
					if (this->mismatchCheckLarmor) {
						mol->select(":"+resid+".+:"+resname+".+:."+nucleus,false,false);
						if (mol->getNAtomSelected()==0){
							std::cerr << "Warning: Mismatch between sequence in Molecule and CS file: " << ":"+resid+".+:"+resname+".+:."+nucleus << std::endl;	
						}
						else {
							this->experimentalCS.insert(std::pair<std::string,double>(key,atof(Misc::trim(s.at(3)).c_str())));
							this->errorCS.insert(std::pair<std::string,double>(key,atof(Misc::trim(s.at(4)).c_str())));						
						}
					}	
					else {			
						this->experimentalCS.insert(std::pair<std::string,double>(key,atof(Misc::trim(s.at(3)).c_str())));
						this->errorCS.insert(std::pair<std::string,double>(key,atof(Misc::trim(s.at(4)).c_str())));
					}
				}
			}
		} 
		else {
			std::cerr << "Error: Unable to open file: " << fchemshift << std::endl;	
			exit(0);
		}
	}
}

void LARMORD::loadParmFile(const std::string fparmfile)
{
	std::ifstream parmFile;
	std::istream* parminp;
	std::string line, key;
	std::vector<std::string> s;
	double alpha;
	int beta;
	std::vector <double> alphas_;
	std::vector <int> betas_;
	if (fparmfile.length() > 0){
		parmFile.open(fparmfile.c_str(), std::ios::in);
		parminp=&parmFile;
		if (parmFile.is_open()){
			while (parminp->good() && !(parminp->eof()))
			{
				getline(*parminp, line);
				Misc::splitStr(line, " ", s, true);            
				if (s.size() ==  6){
					if (this->residueBasedLarmor){
						key = Misc::trim(s.at(0))+":"+Misc::trim(s.at(1))+":"+Misc::trim(s.at(2))+":"+Misc::trim(s.at(3));
					}
					else {
						key = Misc::trim(s.at(1))+":"+Misc::trim(s.at(2))+":"+Misc::trim(s.at(3));
					}
				
					alpha = atof(Misc::trim(s.at(4)).c_str());
					//std::cout << "LARMOR " << key << " " << alpha << std::endl;
					if (this->alphas.count(key) > 0 ){
						//std::cout << "I am here" << std::endl;
						alphas_.clear();
						this->alphas.at(key).push_back(alpha);
					}
					else {
						alphas_.clear();
						alphas_.push_back(alpha);
						this->alphas.insert(std::pair<std::string,std::vector<double> >(key,alphas_));
					}                

					if (this->residueBasedLarmor){
						key = Misc::trim(s.at(0))+":"+Misc::trim(s.at(1))+":"+Misc::trim(s.at(2))+":"+Misc::trim(s.at(3));
					}
					else {
						key = Misc::trim(s.at(1))+":"+Misc::trim(s.at(2))+":"+Misc::trim(s.at(3));
					}
					beta = atoi(Misc::trim(s.at(5)).c_str());
					if (this->betas.count(key) > 0 ){
						betas_.clear();
						this->betas.at(key).push_back(beta);
					}
					else {
						betas_.clear();
						betas_.push_back(beta);
						this->betas.insert(std::pair<std::string,std::vector<int> >(key,betas_));
					}                                
				}
				//std::cout << key << " alpha " << this->alphas.at(key).size() << std::endl;                      
				//std::cout << key << " beta " << this->betas.at(key).size() << std::endl;                      
			}
		}
		else {
			std::cerr << "Error: Unable to open file: " << fparmfile << std::endl;	
			exit(0);		
		}
	}
}

void LARMORD::loadRefFile(const std::string freffile)
{
	std::ifstream refFile;
	std::istream* refinp;
	std::string line;
	std::vector<std::string> s;
	if (freffile.length() > 0){
		refFile.open(freffile.c_str(), std::ios::in);
		refinp=&refFile;
		if (refFile.is_open()) {
			while (refinp->good() && !(refinp->eof()))
			{
				getline(*refinp, line);
				Misc::splitStr(line, " ", s, true);
				if (s.size() ==  3){
						this->randomShifts.insert(std::pair<std::string,double>(Misc::trim(s.at(0))+":"+Misc::trim(s.at(1)),atof(Misc::trim(s.at(2)).c_str())));
				}
			}
		}
		else {
			std::cerr << "Error: Unable to open file: " << freffile << std::endl;	
			exit(0);		                
		}
	}
}

void LARMORD::loadAccFile(const std::string faccfile)
{
	std::ifstream accFile;
	std::istream* accinp;
	std::string line;
	std::vector<std::string> s;
	if (faccfile.length() > 0){
		accFile.open(faccfile.c_str(), std::ios::in);
		accinp=&accFile;
		if (accFile.is_open()){
			while (accinp->good() && !(accinp->eof()))
			{
				getline(*accinp, line);
				Misc::splitStr(line, " ", s, true);
				if (s.size() ==  3){
					if (this->residueBasedLarmor || this->residueBasedWeightsLarmor){
						this->accuracy_weight.insert(std::pair<std::string,double>(Misc::trim(s.at(0))+":"+Misc::trim(s.at(1)),atof(Misc::trim(s.at(2)).c_str())));
						//std::cout << "Accuracy Weights " << Misc::trim(s.at(0))+":"+Misc::trim(s.at(1)) << " " << this->accuracy_weight.at(Misc::trim(s.at(0))) << std::endl;
					} 
					else {
						this->accuracy_weight.insert(std::pair<std::string,double>(Misc::trim(s.at(0)),atof(Misc::trim(s.at(2)).c_str())));
						//std::cout << "Accuracy Weights " << Misc::trim(s.at(0)) << " " << this->accuracy_weight.at(Misc::trim(s.at(0))) << std::endl;							
					}
				}
			}
		}
		else {
			std::cerr << "Error: Unable to open file: " << faccfile << std::endl;	
			exit(0);		                			
		}
	}
}

void LARMORD::loadCorrFile(const std::string fcorrfile)
{
	std::ifstream corrFile;
	std::istream* corrinp;
	std::string line;
	std::vector<std::string> s;
	if (fcorrfile.length() > 0){
		corrFile.open(fcorrfile.c_str(), std::ios::in);
		corrinp=&corrFile;
		if (corrFile.is_open()){
			while (corrinp->good() && !(corrinp->eof()))
			{
				getline(*corrinp, line);
				Misc::splitStr(line, " ", s, true);
				if (s.size() ==  2){
						this->correlation_weight.insert(std::pair<std::string,double>(Misc::trim(s.at(0)),atof(Misc::trim(s.at(1)).c_str())));
				}
			}
		}
		else {
			std::cerr << "Error: Unable to open file: " << fcorrfile << std::endl;	
			exit(0);		                			
		}
	}
}

int LARMORD::getNShiftAtoms()
{
    return (this->shiftAtoms.size());
}

bool LARMORD::getShiftAtom(const std::string &key)
{
    if (this->shiftAtoms.find (key) != this->shiftAtoms.end()){
        return true;
    } else {
        return false;
    }
}

double LARMORD::getExperimentalCS(const std::string &key)
{
    if (this->experimentalCS.find (key) != this->experimentalCS.end()){
        return this->experimentalCS.at(key);
    } else {
        return 0.0;
    }
}

double LARMORD::getErrorCS(const std::string &key)
{
    if (this->errorCS.find (key) != this->errorCS.end()){
        return this->errorCS.at(key);
    } else {
        return 0.0;
    }
}

double LARMORD::getAccuracyWeight(const std::string &key)
{
    if (this->accuracy_weight.find (key) != this->accuracy_weight.end()){
        return this->accuracy_weight.at(key);
    } else {
        return 1.0;
    }
}

double LARMORD::getCorrelationWeight(const std::string &key)
{
    if (this->correlation_weight.find (key) != this->correlation_weight.end()){
        return this->correlation_weight.at(key);
    } else {
        return 1.0;
    }
}

void LARMORD::initializeAccuracyWeight()
{
		this->accuracy_weight.insert(std::pair<std::string,double>("CA",	0.533608491));
		this->accuracy_weight.insert(std::pair<std::string,double>("CB",	0.645544554));
		this->accuracy_weight.insert(std::pair<std::string,double>("C",	0.685740236));
		this->accuracy_weight.insert(std::pair<std::string,double>("HA",	2.420338983));
		this->accuracy_weight.insert(std::pair<std::string,double>("H",	1.247379455));
		this->accuracy_weight.insert(std::pair<std::string,double>("N",	0.236764222));
}

void LARMORD::renameRes(Molecule *mol)
{
    /* rename resname in mol to be consistent with LarmorD parameter file */
    /* GUA */
    mol->renameRes("RG","GUA");
    mol->renameRes("RG3","GUA");
    mol->renameRes("RG5","GUA");
    mol->renameRes("G","GUA");  
    /* ADE */
    mol->renameRes("RA","ADE");
    mol->renameRes("RA3","ADE");
    mol->renameRes("RA5","ADE");
    mol->renameRes("A","ADE");  
    /* CYT */
    mol->renameRes("RC","CYT");
    mol->renameRes("RC3","CYT");
    mol->renameRes("RC5","CYT");
    mol->renameRes("C","CYT");  
    /* URA */
    mol->renameRes("RU","URA");
    mol->renameRes("RU3","URA");
    mol->renameRes("RU5","URA");
    mol->renameRes("U","URA");  
    /* HIS */
    mol->renameRes("HID","HIS");
    mol->renameRes("HIE","HIS");
    mol->renameRes("HIP","HIS");
    mol->renameRes("HSD","HIS");
    mol->renameRes("HSE","HIS");
    mol->renameRes("HSP","HIS");
    /* rename resname in mol to be consistent with LarmorD parameter file */
    /* H2'1 => H2' */
    mol->renameAtom("H2'1","H2'");
    /* H2'2 => H2'' */
    mol->renameAtom("H2'2","H2''");
    /* H5'1 => H5' */
    mol->renameAtom("H5'1","H5'");
    /* H5'2 => H5'' */
    mol->renameAtom("H5'2","H5''");
    /* HN => H */
    mol->renameAtom("HN","H");
    //mol->writePDB();
}

void LARMORD::initializeShiftAtoms()
{
    this->shiftAtoms.insert(std::pair<std::string,bool>("P",true));
    this->shiftAtoms.insert(std::pair<std::string,bool>("O3'",true));
    this->shiftAtoms.insert(std::pair<std::string,bool>("O5'",true));
    this->shiftAtoms.insert(std::pair<std::string,bool>("O4'",true));    
    this->shiftAtoms.insert(std::pair<std::string,bool>("H1'",true));
    this->shiftAtoms.insert(std::pair<std::string,bool>("H2'",true));
    this->shiftAtoms.insert(std::pair<std::string,bool>("H3'",true));
    this->shiftAtoms.insert(std::pair<std::string,bool>("H4'",true));
    this->shiftAtoms.insert(std::pair<std::string,bool>("H5'",true));
    this->shiftAtoms.insert(std::pair<std::string,bool>("H5''",true));
    this->shiftAtoms.insert(std::pair<std::string,bool>("H2",true));
    this->shiftAtoms.insert(std::pair<std::string,bool>("H5",true));
    this->shiftAtoms.insert(std::pair<std::string,bool>("H6",true));
    this->shiftAtoms.insert(std::pair<std::string,bool>("H8",true));
    this->shiftAtoms.insert(std::pair<std::string,bool>("H1",true));
    this->shiftAtoms.insert(std::pair<std::string,bool>("H3",true));
    this->shiftAtoms.insert(std::pair<std::string,bool>("C1'",true));
    this->shiftAtoms.insert(std::pair<std::string,bool>("C2'",true));
    this->shiftAtoms.insert(std::pair<std::string,bool>("C3'",true));
    this->shiftAtoms.insert(std::pair<std::string,bool>("C4'",true));
    this->shiftAtoms.insert(std::pair<std::string,bool>("C5'",true));
    this->shiftAtoms.insert(std::pair<std::string,bool>("C2",true));
    this->shiftAtoms.insert(std::pair<std::string,bool>("C5",true));
    this->shiftAtoms.insert(std::pair<std::string,bool>("C6",true));
    this->shiftAtoms.insert(std::pair<std::string,bool>("C8",true));
    this->shiftAtoms.insert(std::pair<std::string,bool>("N1",true));
    this->shiftAtoms.insert(std::pair<std::string,bool>("N3",true));
    this->shiftAtoms.insert(std::pair<std::string,bool>("CA",true));
    this->shiftAtoms.insert(std::pair<std::string,bool>("CB",true));
    this->shiftAtoms.insert(std::pair<std::string,bool>("C",true));
    this->shiftAtoms.insert(std::pair<std::string,bool>("N",true));
    this->shiftAtoms.insert(std::pair<std::string,bool>("HA",true));
    this->shiftAtoms.insert(std::pair<std::string,bool>("H",true));
}

void LARMORD::initializeRandomShifts()
{
    /* generated using: awk '{print "this->randomShifts.insert(std::pair<std::string,double>(\""$1":"$2"\","$3"));"}' randomcoil.dat */
    this->randomShifts.insert(std::pair<std::string,double>("GUA:H5''",3.58));
    this->randomShifts.insert(std::pair<std::string,double>("GUA:C3'",73.9));
    this->randomShifts.insert(std::pair<std::string,double>("GUA:C1'",93.1));
    this->randomShifts.insert(std::pair<std::string,double>("GUA:C5'",66.3));
    this->randomShifts.insert(std::pair<std::string,double>("GUA:C4'",83.1));
    this->randomShifts.insert(std::pair<std::string,double>("GUA:C2'",76.3));
    this->randomShifts.insert(std::pair<std::string,double>("GUA:N1",146.3));
    this->randomShifts.insert(std::pair<std::string,double>("GUA:H2'",4.02));
    this->randomShifts.insert(std::pair<std::string,double>("GUA:H1'",5.05));
    this->randomShifts.insert(std::pair<std::string,double>("GUA:H4'",4.0));
    this->randomShifts.insert(std::pair<std::string,double>("GUA:H8",7.69));
    this->randomShifts.insert(std::pair<std::string,double>("GUA:H5'",3.78));
    this->randomShifts.insert(std::pair<std::string,double>("GUA:H3'",4.14));
    this->randomShifts.insert(std::pair<std::string,double>("GUA:H1",12.36));
    this->randomShifts.insert(std::pair<std::string,double>("GUA:C8",137.9));
    this->randomShifts.insert(std::pair<std::string,double>("ADE:H4'",4.0));
    this->randomShifts.insert(std::pair<std::string,double>("ADE:H8",8.6));
    this->randomShifts.insert(std::pair<std::string,double>("ADE:H5'",3.78));
    this->randomShifts.insert(std::pair<std::string,double>("ADE:H3'",4.14));
    this->randomShifts.insert(std::pair<std::string,double>("ADE:H2",8.68));
    this->randomShifts.insert(std::pair<std::string,double>("ADE:C8",141.2));
    this->randomShifts.insert(std::pair<std::string,double>("ADE:C3'",73.9));
    this->randomShifts.insert(std::pair<std::string,double>("ADE:C1'",93.1));
    this->randomShifts.insert(std::pair<std::string,double>("ADE:C5'",66.3));
    this->randomShifts.insert(std::pair<std::string,double>("ADE:H5''",3.58));
    this->randomShifts.insert(std::pair<std::string,double>("ADE:C4'",83.1));
    this->randomShifts.insert(std::pair<std::string,double>("ADE:C2'",76.3));
    this->randomShifts.insert(std::pair<std::string,double>("ADE:C2",154.3));
    this->randomShifts.insert(std::pair<std::string,double>("ADE:H2'",4.02));
    this->randomShifts.insert(std::pair<std::string,double>("ADE:N6",80.5));
    this->randomShifts.insert(std::pair<std::string,double>("ADE:H1'",5.03));
    this->randomShifts.insert(std::pair<std::string,double>("URA:H4'",4.0));
    this->randomShifts.insert(std::pair<std::string,double>("URA:H5'",3.78));
    this->randomShifts.insert(std::pair<std::string,double>("URA:H3'",4.14));
    this->randomShifts.insert(std::pair<std::string,double>("URA:H3",13.15));
    this->randomShifts.insert(std::pair<std::string,double>("URA:C3'",73.9));
    this->randomShifts.insert(std::pair<std::string,double>("URA:C1'",93.1));
    this->randomShifts.insert(std::pair<std::string,double>("URA:H5",5.95));
    this->randomShifts.insert(std::pair<std::string,double>("URA:C6",142.0));
    this->randomShifts.insert(std::pair<std::string,double>("URA:C5'",66.3));
    this->randomShifts.insert(std::pair<std::string,double>("URA:H5''",3.58));
    this->randomShifts.insert(std::pair<std::string,double>("URA:C4'",83.1));
    this->randomShifts.insert(std::pair<std::string,double>("URA:C2'",76.3));
    this->randomShifts.insert(std::pair<std::string,double>("URA:H6",7.85));
    this->randomShifts.insert(std::pair<std::string,double>("URA:N2",160.5));
    this->randomShifts.insert(std::pair<std::string,double>("URA:N3",158.82));
    this->randomShifts.insert(std::pair<std::string,double>("URA:H2'",4.02));
    this->randomShifts.insert(std::pair<std::string,double>("URA:C5",103.5));
    this->randomShifts.insert(std::pair<std::string,double>("URA:H1'",5.6));
    this->randomShifts.insert(std::pair<std::string,double>("CYT:H4'",4.0));
    this->randomShifts.insert(std::pair<std::string,double>("CYT:H5'",3.78));
    this->randomShifts.insert(std::pair<std::string,double>("CYT:H3'",4.14));
    this->randomShifts.insert(std::pair<std::string,double>("CYT:H6",7.8));
    this->randomShifts.insert(std::pair<std::string,double>("CYT:C3'",73.9));
    this->randomShifts.insert(std::pair<std::string,double>("CYT:C1'",93.1));
    this->randomShifts.insert(std::pair<std::string,double>("CYT:H5",6.2));
    this->randomShifts.insert(std::pair<std::string,double>("CYT:C5'",66.3));
    this->randomShifts.insert(std::pair<std::string,double>("CYT:H5''",3.58));
    this->randomShifts.insert(std::pair<std::string,double>("CYT:C4'",83.1));
    this->randomShifts.insert(std::pair<std::string,double>("CYT:C2'",76.3));
    this->randomShifts.insert(std::pair<std::string,double>("CYT:H2'",4.02));
    this->randomShifts.insert(std::pair<std::string,double>("CYT:N4",97.8));
    this->randomShifts.insert(std::pair<std::string,double>("CYT:C6",142.2));
    this->randomShifts.insert(std::pair<std::string,double>("CYT:C5",98.2));
    this->randomShifts.insert(std::pair<std::string,double>("CYT:H1'",5.28));
}

void LARMORD::initializeAlpha()
{
    this->alphas.insert(std::pair<std::string,double>("H5'':GUA:C1'",3.67815480447));
    this->alphas.insert(std::pair<std::string,double>("H5'':GUA:C2'",1.63100734262));
    this->alphas.insert(std::pair<std::string,double>("H5'':GUA:C3'",2.89281770448));
    this->alphas.insert(std::pair<std::string,double>("H5'':GUA:C4'",2.3548429102));
    this->alphas.insert(std::pair<std::string,double>("H5'':GUA:C5'",0.0321335101989));
    this->alphas.insert(std::pair<std::string,double>("H5'':GUA:P",-0.171617457218));
    this->alphas.insert(std::pair<std::string,double>("H5'':GUA:O5'",1.57702037737));
    this->alphas.insert(std::pair<std::string,double>("H5'':GUA:O3'",0.326168827858));
    this->alphas.insert(std::pair<std::string,double>("H5'':GUA:C2",-6.64414148375));
    this->alphas.insert(std::pair<std::string,double>("H5'':GUA:C4",-14.2514277865));
    this->alphas.insert(std::pair<std::string,double>("H5'':GUA:C5",-18.4752206521));
    this->alphas.insert(std::pair<std::string,double>("H5'':GUA:C6",-23.3660033205));
    this->alphas.insert(std::pair<std::string,double>("H5'':GUA:C8",11.2256815971));
    this->alphas.insert(std::pair<std::string,double>("H5'':GUA:N1",-15.0357687755));
    this->alphas.insert(std::pair<std::string,double>("H5'':GUA:N2",18.8071249973));
    this->alphas.insert(std::pair<std::string,double>("H5'':GUA:N3",7.49438758708));
    this->alphas.insert(std::pair<std::string,double>("H5'':GUA:N7",6.47382380489));
    this->alphas.insert(std::pair<std::string,double>("H5'':GUA:N9",-1.33876233081));
    this->alphas.insert(std::pair<std::string,double>("H5'':GUA:O6",13.3592994749));
    this->alphas.insert(std::pair<std::string,double>("H5'':ADE:C1'",-6.87501637071));
    this->alphas.insert(std::pair<std::string,double>("H5'':ADE:C2'",-0.125715293229));
    this->alphas.insert(std::pair<std::string,double>("H5'':ADE:C3'",0.934354550638));
    this->alphas.insert(std::pair<std::string,double>("H5'':ADE:C4'",0.256912097145));
    this->alphas.insert(std::pair<std::string,double>("H5'':ADE:C5'",0.571155032105));
    this->alphas.insert(std::pair<std::string,double>("H5'':ADE:P",-0.433798841887));
    this->alphas.insert(std::pair<std::string,double>("H5'':ADE:O5'",1.5556382895));
    this->alphas.insert(std::pair<std::string,double>("H5'':ADE:O3'",0.852536186953));
    this->alphas.insert(std::pair<std::string,double>("H5'':ADE:C2",19.465835991));
    this->alphas.insert(std::pair<std::string,double>("H5'':ADE:C4",-2.63457000853));
    this->alphas.insert(std::pair<std::string,double>("H5'':ADE:C5",-6.20431462114));
    this->alphas.insert(std::pair<std::string,double>("H5'':ADE:C6",-5.52564589385));
    this->alphas.insert(std::pair<std::string,double>("H5'':ADE:C8",3.72779586352));
    this->alphas.insert(std::pair<std::string,double>("H5'':ADE:N1",8.97415774656));
    this->alphas.insert(std::pair<std::string,double>("H5'':ADE:N3",-3.65500648528));
    this->alphas.insert(std::pair<std::string,double>("H5'':ADE:N6",-21.6733704635));
    this->alphas.insert(std::pair<std::string,double>("H5'':ADE:N7",0.39668878159));
    this->alphas.insert(std::pair<std::string,double>("H5'':ADE:N9",1.89834133117));
    this->alphas.insert(std::pair<std::string,double>("H5'':URA:C1'",-7.8361401694));
    this->alphas.insert(std::pair<std::string,double>("H5'':URA:C2'",-1.05475191362));
    this->alphas.insert(std::pair<std::string,double>("H5'':URA:C3'",-1.23351633568));
    this->alphas.insert(std::pair<std::string,double>("H5'':URA:C4'",4.8561240664));
    this->alphas.insert(std::pair<std::string,double>("H5'':URA:C5'",0.0));
    this->alphas.insert(std::pair<std::string,double>("H5'':URA:P",1.86617611577));
    this->alphas.insert(std::pair<std::string,double>("H5'':URA:O5'",1.01063675972));
    this->alphas.insert(std::pair<std::string,double>("H5'':URA:O3'",0.949438715561));
    this->alphas.insert(std::pair<std::string,double>("H5'':URA:C2",3.44987727435));
    this->alphas.insert(std::pair<std::string,double>("H5'':URA:C4",-7.50237209988));
    this->alphas.insert(std::pair<std::string,double>("H5'':URA:C5",-14.857005414));
    this->alphas.insert(std::pair<std::string,double>("H5'':URA:C6",3.44656331225));
    this->alphas.insert(std::pair<std::string,double>("H5'':URA:N1",2.71399930928));
    this->alphas.insert(std::pair<std::string,double>("H5'':URA:N3",19.9073925365));
    this->alphas.insert(std::pair<std::string,double>("H5'':URA:O4",23.4150612805));
    this->alphas.insert(std::pair<std::string,double>("H5'':CYT:C1'",0.706640452452));
    this->alphas.insert(std::pair<std::string,double>("H5'':CYT:C2'",0.166171606825));
    this->alphas.insert(std::pair<std::string,double>("H5'':CYT:C3'",-3.1035073294));
    this->alphas.insert(std::pair<std::string,double>("H5'':CYT:C4'",-1.57037444807));
    this->alphas.insert(std::pair<std::string,double>("H5'':CYT:C5'",1.32626298398));
    this->alphas.insert(std::pair<std::string,double>("H5'':CYT:P",-0.440649244862));
    this->alphas.insert(std::pair<std::string,double>("H5'':CYT:O5'",-2.63138895134));
    this->alphas.insert(std::pair<std::string,double>("H5'':CYT:O3'",1.92397549643));
    this->alphas.insert(std::pair<std::string,double>("H5'':CYT:C2",8.20045630532));
    this->alphas.insert(std::pair<std::string,double>("H5'':CYT:C4",1.29679400155));
    this->alphas.insert(std::pair<std::string,double>("H5'':CYT:C5",4.77736668512));
    this->alphas.insert(std::pair<std::string,double>("H5'':CYT:C6",1.05351000601));
    this->alphas.insert(std::pair<std::string,double>("H5'':CYT:N1",2.78643038779));
    this->alphas.insert(std::pair<std::string,double>("H5'':CYT:N3",6.38091065052));
    this->alphas.insert(std::pair<std::string,double>("H5'':CYT:N4",-12.0144974986));
    this->alphas.insert(std::pair<std::string,double>("H5'':CYT:O2",-2.71778957123));
    this->alphas.insert(std::pair<std::string,double>("C3':GUA:C1'",-3.00751015521));
    this->alphas.insert(std::pair<std::string,double>("C3':GUA:C2'",-1.54074163816));
    this->alphas.insert(std::pair<std::string,double>("C3':GUA:C3'",-27.7836942588));
    this->alphas.insert(std::pair<std::string,double>("C3':GUA:C4'",-1.27318800607));
    this->alphas.insert(std::pair<std::string,double>("C3':GUA:C5'",9.59221362253));
    this->alphas.insert(std::pair<std::string,double>("C3':GUA:P",71.3126488556));
    this->alphas.insert(std::pair<std::string,double>("C3':GUA:O5'",-20.5766593968));
    this->alphas.insert(std::pair<std::string,double>("C3':GUA:O3'",-0.606352937576));
    this->alphas.insert(std::pair<std::string,double>("C3':GUA:C2",-22.8316282596));
    this->alphas.insert(std::pair<std::string,double>("C3':GUA:C4",48.9774249661));
    this->alphas.insert(std::pair<std::string,double>("C3':GUA:C5",113.012349032));
    this->alphas.insert(std::pair<std::string,double>("C3':GUA:C6",40.5761254389));
    this->alphas.insert(std::pair<std::string,double>("C3':GUA:C8",-118.759576817));
    this->alphas.insert(std::pair<std::string,double>("C3':GUA:N1",-87.8860464611));
    this->alphas.insert(std::pair<std::string,double>("C3':GUA:N2",-40.4439078282));
    this->alphas.insert(std::pair<std::string,double>("C3':GUA:N3",67.9037015402));
    this->alphas.insert(std::pair<std::string,double>("C3':GUA:N7",64.6678187857));
    this->alphas.insert(std::pair<std::string,double>("C3':GUA:N9",-64.2137999992));
    this->alphas.insert(std::pair<std::string,double>("C3':GUA:O6",26.2516229808));
    this->alphas.insert(std::pair<std::string,double>("C3':ADE:C1'",8.67950913628));
    this->alphas.insert(std::pair<std::string,double>("C3':ADE:C2'",-0.0906529556619));
    this->alphas.insert(std::pair<std::string,double>("C3':ADE:C3'",80.788456279));
    this->alphas.insert(std::pair<std::string,double>("C3':ADE:C4'",-0.57773676838));
    this->alphas.insert(std::pair<std::string,double>("C3':ADE:C5'",7.82457724041));
    this->alphas.insert(std::pair<std::string,double>("C3':ADE:P",66.854253122));
    this->alphas.insert(std::pair<std::string,double>("C3':ADE:O5'",-22.4203583793));
    this->alphas.insert(std::pair<std::string,double>("C3':ADE:O3'",-0.622749715407));
    this->alphas.insert(std::pair<std::string,double>("C3':ADE:C2",-60.6454854662));
    this->alphas.insert(std::pair<std::string,double>("C3':ADE:C4",-7.70758976477));
    this->alphas.insert(std::pair<std::string,double>("C3':ADE:C5",26.7850361983));
    this->alphas.insert(std::pair<std::string,double>("C3':ADE:C6",54.0703855716));
    this->alphas.insert(std::pair<std::string,double>("C3':ADE:C8",-42.6236884114));
    this->alphas.insert(std::pair<std::string,double>("C3':ADE:N1",-19.3843799162));
    this->alphas.insert(std::pair<std::string,double>("C3':ADE:N3",-38.2503385293));
    this->alphas.insert(std::pair<std::string,double>("C3':ADE:N6",92.6338030812));
    this->alphas.insert(std::pair<std::string,double>("C3':ADE:N7",-13.4529123783));
    this->alphas.insert(std::pair<std::string,double>("C3':ADE:N9",-63.6874750901));
    this->alphas.insert(std::pair<std::string,double>("C3':URA:C1'",16.5836400001));
    this->alphas.insert(std::pair<std::string,double>("C3':URA:C2'",0.0738201348987));
    this->alphas.insert(std::pair<std::string,double>("C3':URA:C3'",205.965675527));
    this->alphas.insert(std::pair<std::string,double>("C3':URA:C4'",0.981957814507));
    this->alphas.insert(std::pair<std::string,double>("C3':URA:C5'",18.3972140189));
    this->alphas.insert(std::pair<std::string,double>("C3':URA:P",55.2118583147));
    this->alphas.insert(std::pair<std::string,double>("C3':URA:O5'",-12.9658867755));
    this->alphas.insert(std::pair<std::string,double>("C3':URA:O3'",1.15410456364));
    this->alphas.insert(std::pair<std::string,double>("C3':URA:C2",-56.8055224661));
    this->alphas.insert(std::pair<std::string,double>("C3':URA:C4",-4.54827872323));
    this->alphas.insert(std::pair<std::string,double>("C3':URA:C5",71.5910685132));
    this->alphas.insert(std::pair<std::string,double>("C3':URA:C6",-55.673454124));
    this->alphas.insert(std::pair<std::string,double>("C3':URA:N1",-152.855079302));
    this->alphas.insert(std::pair<std::string,double>("C3':URA:N3",-82.8401158362));
    this->alphas.insert(std::pair<std::string,double>("C3':URA:O4",-10.5439074089));
    this->alphas.insert(std::pair<std::string,double>("C3':CYT:C1'",8.61734600999));
    this->alphas.insert(std::pair<std::string,double>("C3':CYT:C2'",-0.817553475382));
    this->alphas.insert(std::pair<std::string,double>("C3':CYT:C3'",74.8750720266));
    this->alphas.insert(std::pair<std::string,double>("C3':CYT:C4'",-3.04428988812));
    this->alphas.insert(std::pair<std::string,double>("C3':CYT:C5'",5.36871431506));
    this->alphas.insert(std::pair<std::string,double>("C3':CYT:P",53.5397757171));
    this->alphas.insert(std::pair<std::string,double>("C3':CYT:O5'",-17.7328713715));
    this->alphas.insert(std::pair<std::string,double>("C3':CYT:O3'",-1.43184594746));
    this->alphas.insert(std::pair<std::string,double>("C3':CYT:C2",-34.8814718485));
    this->alphas.insert(std::pair<std::string,double>("C3':CYT:C4",12.7609629701));
    this->alphas.insert(std::pair<std::string,double>("C3':CYT:C5",50.4665629912));
    this->alphas.insert(std::pair<std::string,double>("C3':CYT:C6",-35.4908266074));
    this->alphas.insert(std::pair<std::string,double>("C3':CYT:N1",-36.0724446145));
    this->alphas.insert(std::pair<std::string,double>("C3':CYT:N3",-74.1084337498));
    this->alphas.insert(std::pair<std::string,double>("C3':CYT:N4",-3.22867035327));
    this->alphas.insert(std::pair<std::string,double>("C3':CYT:O2",-69.7971671866));
    this->alphas.insert(std::pair<std::string,double>("C1':GUA:C1'",114.82714112));
    this->alphas.insert(std::pair<std::string,double>("C1':GUA:C2'",-3.0472483545));
    this->alphas.insert(std::pair<std::string,double>("C1':GUA:C3'",2.56562882759));
    this->alphas.insert(std::pair<std::string,double>("C1':GUA:C4'",2.55213521656));
    this->alphas.insert(std::pair<std::string,double>("C1':GUA:C5'",18.0599695162));
    this->alphas.insert(std::pair<std::string,double>("C1':GUA:P",-24.4303586772));
    this->alphas.insert(std::pair<std::string,double>("C1':GUA:O5'",-15.0621200083));
    this->alphas.insert(std::pair<std::string,double>("C1':GUA:O3'",-88.0697364179));
    this->alphas.insert(std::pair<std::string,double>("C1':GUA:C2",-17.444879619));
    this->alphas.insert(std::pair<std::string,double>("C1':GUA:C4",0.328054956697));
    this->alphas.insert(std::pair<std::string,double>("C1':GUA:C5",-3.81465601003));
    this->alphas.insert(std::pair<std::string,double>("C1':GUA:C6",-39.5873330747));
    this->alphas.insert(std::pair<std::string,double>("C1':GUA:C8",11.248000247));
    this->alphas.insert(std::pair<std::string,double>("C1':GUA:N1",7.35925865133));
    this->alphas.insert(std::pair<std::string,double>("C1':GUA:N2",18.2584850105));
    this->alphas.insert(std::pair<std::string,double>("C1':GUA:N3",2.23613444361));
    this->alphas.insert(std::pair<std::string,double>("C1':GUA:N7",20.1672098418));
    this->alphas.insert(std::pair<std::string,double>("C1':GUA:N9",-0.444676697704));
    this->alphas.insert(std::pair<std::string,double>("C1':GUA:O6",5.30464135861));
    this->alphas.insert(std::pair<std::string,double>("C1':ADE:C1'",26.7803620224));
    this->alphas.insert(std::pair<std::string,double>("C1':ADE:C2'",0.0643894906118));
    this->alphas.insert(std::pair<std::string,double>("C1':ADE:C3'",0.428229417609));
    this->alphas.insert(std::pair<std::string,double>("C1':ADE:C4'",-5.22500751805));
    this->alphas.insert(std::pair<std::string,double>("C1':ADE:C5'",-6.10048656029));
    this->alphas.insert(std::pair<std::string,double>("C1':ADE:P",-19.2046858444));
    this->alphas.insert(std::pair<std::string,double>("C1':ADE:O5'",59.139553853));
    this->alphas.insert(std::pair<std::string,double>("C1':ADE:O3'",-66.3828700687));
    this->alphas.insert(std::pair<std::string,double>("C1':ADE:C2",14.8223131853));
    this->alphas.insert(std::pair<std::string,double>("C1':ADE:C4",2.45322773921));
    this->alphas.insert(std::pair<std::string,double>("C1':ADE:C5",0.948741218553));
    this->alphas.insert(std::pair<std::string,double>("C1':ADE:C6",-9.45700712442));
    this->alphas.insert(std::pair<std::string,double>("C1':ADE:C8",-7.33934779384));
    this->alphas.insert(std::pair<std::string,double>("C1':ADE:N1",-28.1405418759));
    this->alphas.insert(std::pair<std::string,double>("C1':ADE:N3",20.8037714018));
    this->alphas.insert(std::pair<std::string,double>("C1':ADE:N6",-3.07801267047));
    this->alphas.insert(std::pair<std::string,double>("C1':ADE:N7",-4.27413038584));
    this->alphas.insert(std::pair<std::string,double>("C1':ADE:N9",-3.75134428442));
    this->alphas.insert(std::pair<std::string,double>("C1':URA:C1'",85.4492359609));
    this->alphas.insert(std::pair<std::string,double>("C1':URA:C2'",-0.500419726263));
    this->alphas.insert(std::pair<std::string,double>("C1':URA:C3'",6.04410355845));
    this->alphas.insert(std::pair<std::string,double>("C1':URA:C4'",2.20333963079));
    this->alphas.insert(std::pair<std::string,double>("C1':URA:C5'",14.492473609));
    this->alphas.insert(std::pair<std::string,double>("C1':URA:P",35.1018638346));
    this->alphas.insert(std::pair<std::string,double>("C1':URA:O5'",12.5044100403));
    this->alphas.insert(std::pair<std::string,double>("C1':URA:O3'",-132.934409968));
    this->alphas.insert(std::pair<std::string,double>("C1':URA:C2",5.87483874579));
    this->alphas.insert(std::pair<std::string,double>("C1':URA:C4",-9.97429015969));
    this->alphas.insert(std::pair<std::string,double>("C1':URA:C5",6.4669395485));
    this->alphas.insert(std::pair<std::string,double>("C1':URA:C6",4.74730058919));
    this->alphas.insert(std::pair<std::string,double>("C1':URA:N1",0.0771970231211));
    this->alphas.insert(std::pair<std::string,double>("C1':URA:N3",1.84090735735));
    this->alphas.insert(std::pair<std::string,double>("C1':URA:O4",-18.2846256053));
    this->alphas.insert(std::pair<std::string,double>("C1':CYT:C1'",72.5621426102));
    this->alphas.insert(std::pair<std::string,double>("C1':CYT:C2'",-0.496183941411));
    this->alphas.insert(std::pair<std::string,double>("C1':CYT:C3'",-5.30736061501));
    this->alphas.insert(std::pair<std::string,double>("C1':CYT:C4'",2.14590174078));
    this->alphas.insert(std::pair<std::string,double>("C1':CYT:C5'",21.416401533));
    this->alphas.insert(std::pair<std::string,double>("C1':CYT:P",45.0099682134));
    this->alphas.insert(std::pair<std::string,double>("C1':CYT:O5'",-11.2962402263));
    this->alphas.insert(std::pair<std::string,double>("C1':CYT:O3'",-91.128074134));
    this->alphas.insert(std::pair<std::string,double>("C1':CYT:C2",3.30000689038));
    this->alphas.insert(std::pair<std::string,double>("C1':CYT:C4",6.72735787053));
    this->alphas.insert(std::pair<std::string,double>("C1':CYT:C5",6.35239516973));
    this->alphas.insert(std::pair<std::string,double>("C1':CYT:C6",-6.80540362903));
    this->alphas.insert(std::pair<std::string,double>("C1':CYT:N1",1.0097976237));
    this->alphas.insert(std::pair<std::string,double>("C1':CYT:N3",22.3065471995));
    this->alphas.insert(std::pair<std::string,double>("C1':CYT:N4",0.453608361434));
    this->alphas.insert(std::pair<std::string,double>("C1':CYT:O2",8.73664288121));
    this->alphas.insert(std::pair<std::string,double>("C5':GUA:C1'",16.1493679335));
    this->alphas.insert(std::pair<std::string,double>("C5':GUA:C2'",-21.0820561809));
    this->alphas.insert(std::pair<std::string,double>("C5':GUA:C3'",0.185243662952));
    this->alphas.insert(std::pair<std::string,double>("C5':GUA:C4'",0.133595380896));
    this->alphas.insert(std::pair<std::string,double>("C5':GUA:C5'",99.1730582404));
    this->alphas.insert(std::pair<std::string,double>("C5':GUA:P",-9.23902438957));
    this->alphas.insert(std::pair<std::string,double>("C5':GUA:O5'",0.0));
    this->alphas.insert(std::pair<std::string,double>("C5':GUA:O3'",-13.456658537));
    this->alphas.insert(std::pair<std::string,double>("C5':GUA:C2",-27.4045018974));
    this->alphas.insert(std::pair<std::string,double>("C5':GUA:C4",73.5029165562));
    this->alphas.insert(std::pair<std::string,double>("C5':GUA:C5",45.4035136869));
    this->alphas.insert(std::pair<std::string,double>("C5':GUA:C6",-8.06359109244));
    this->alphas.insert(std::pair<std::string,double>("C5':GUA:C8",-37.6815071775));
    this->alphas.insert(std::pair<std::string,double>("C5':GUA:N1",-104.53349979));
    this->alphas.insert(std::pair<std::string,double>("C5':GUA:N2",-21.6711379422));
    this->alphas.insert(std::pair<std::string,double>("C5':GUA:N3",93.5781163699));
    this->alphas.insert(std::pair<std::string,double>("C5':GUA:N7",-7.79109679095));
    this->alphas.insert(std::pair<std::string,double>("C5':GUA:N9",18.1875985544));
    this->alphas.insert(std::pair<std::string,double>("C5':GUA:O6",-9.65187163721));
    this->alphas.insert(std::pair<std::string,double>("C5':ADE:C1'",85.0499175912));
    this->alphas.insert(std::pair<std::string,double>("C5':ADE:C2'",1.38512954512));
    this->alphas.insert(std::pair<std::string,double>("C5':ADE:C3'",10.8117562756));
    this->alphas.insert(std::pair<std::string,double>("C5':ADE:C4'",-0.222493348115));
    this->alphas.insert(std::pair<std::string,double>("C5':ADE:C5'",158.038652504));
    this->alphas.insert(std::pair<std::string,double>("C5':ADE:P",-35.0427061927));
    this->alphas.insert(std::pair<std::string,double>("C5':ADE:O5'",-0.40607070889));
    this->alphas.insert(std::pair<std::string,double>("C5':ADE:O3'",-32.006828376));
    this->alphas.insert(std::pair<std::string,double>("C5':ADE:C2",-75.6750732327));
    this->alphas.insert(std::pair<std::string,double>("C5':ADE:C4",55.6709169916));
    this->alphas.insert(std::pair<std::string,double>("C5':ADE:C5",34.9845590279));
    this->alphas.insert(std::pair<std::string,double>("C5':ADE:C6",5.1415296815));
    this->alphas.insert(std::pair<std::string,double>("C5':ADE:C8",-56.8649710038));
    this->alphas.insert(std::pair<std::string,double>("C5':ADE:N1",-26.5128323222));
    this->alphas.insert(std::pair<std::string,double>("C5':ADE:N3",0.545911414774));
    this->alphas.insert(std::pair<std::string,double>("C5':ADE:N6",-11.6463603549));
    this->alphas.insert(std::pair<std::string,double>("C5':ADE:N7",-15.8195414579));
    this->alphas.insert(std::pair<std::string,double>("C5':ADE:N9",30.9536830257));
    this->alphas.insert(std::pair<std::string,double>("C5':URA:C1'",59.024605134));
    this->alphas.insert(std::pair<std::string,double>("C5':URA:C2'",53.3203542733));
    this->alphas.insert(std::pair<std::string,double>("C5':URA:C3'",-0.13508670821));
    this->alphas.insert(std::pair<std::string,double>("C5':URA:C4'",-2.69011187991));
    this->alphas.insert(std::pair<std::string,double>("C5':URA:C5'",102.966073736));
    this->alphas.insert(std::pair<std::string,double>("C5':URA:P",-22.139267743));
    this->alphas.insert(std::pair<std::string,double>("C5':URA:O5'",-1.77225149476));
    this->alphas.insert(std::pair<std::string,double>("C5':URA:O3'",-39.1065543578));
    this->alphas.insert(std::pair<std::string,double>("C5':URA:C2",33.023481697));
    this->alphas.insert(std::pair<std::string,double>("C5':URA:C4",-25.5913753432));
    this->alphas.insert(std::pair<std::string,double>("C5':URA:C5",107.801823979));
    this->alphas.insert(std::pair<std::string,double>("C5':URA:C6",-41.014986907));
    this->alphas.insert(std::pair<std::string,double>("C5':URA:N1",-4.02590822327));
    this->alphas.insert(std::pair<std::string,double>("C5':URA:N3",-45.9609877183));
    this->alphas.insert(std::pair<std::string,double>("C5':URA:O4",-102.138531804));
    this->alphas.insert(std::pair<std::string,double>("C5':CYT:C1'",33.8730095875));
    this->alphas.insert(std::pair<std::string,double>("C5':CYT:C2'",6.7506764108));
    this->alphas.insert(std::pair<std::string,double>("C5':CYT:C3'",-6.48897785418));
    this->alphas.insert(std::pair<std::string,double>("C5':CYT:C4'",-0.427028899144));
    this->alphas.insert(std::pair<std::string,double>("C5':CYT:C5'",28.1193528894));
    this->alphas.insert(std::pair<std::string,double>("C5':CYT:P",-26.4436785779));
    this->alphas.insert(std::pair<std::string,double>("C5':CYT:O5'",-2.21510693602));
    this->alphas.insert(std::pair<std::string,double>("C5':CYT:O3'",-10.9233688901));
    this->alphas.insert(std::pair<std::string,double>("C5':CYT:C2",13.2326178755));
    this->alphas.insert(std::pair<std::string,double>("C5':CYT:C4",-5.62022587458));
    this->alphas.insert(std::pair<std::string,double>("C5':CYT:C5",144.158519149));
    this->alphas.insert(std::pair<std::string,double>("C5':CYT:C6",-11.8528825557));
    this->alphas.insert(std::pair<std::string,double>("C5':CYT:N1",12.7977908631));
    this->alphas.insert(std::pair<std::string,double>("C5':CYT:N3",-131.233604726));
    this->alphas.insert(std::pair<std::string,double>("C5':CYT:N4",-87.9010855548));
    this->alphas.insert(std::pair<std::string,double>("C5':CYT:O2",26.6207521659));
    this->alphas.insert(std::pair<std::string,double>("C4':GUA:C1'",8.71060263931));
    this->alphas.insert(std::pair<std::string,double>("C4':GUA:C2'",-25.4513319283));
    this->alphas.insert(std::pair<std::string,double>("C4':GUA:C3'",0.764317392138));
    this->alphas.insert(std::pair<std::string,double>("C4':GUA:C4'",-19.5284692834));
    this->alphas.insert(std::pair<std::string,double>("C4':GUA:C5'",0.789558244056));
    this->alphas.insert(std::pair<std::string,double>("C4':GUA:P",12.3535043115));
    this->alphas.insert(std::pair<std::string,double>("C4':GUA:O5'",17.4371935171));
    this->alphas.insert(std::pair<std::string,double>("C4':GUA:O3'",1.95279982285));
    this->alphas.insert(std::pair<std::string,double>("C4':GUA:C2",-38.932298418));
    this->alphas.insert(std::pair<std::string,double>("C4':GUA:C4",59.6548632631));
    this->alphas.insert(std::pair<std::string,double>("C4':GUA:C5",57.7096288411));
    this->alphas.insert(std::pair<std::string,double>("C4':GUA:C6",42.8760300646));
    this->alphas.insert(std::pair<std::string,double>("C4':GUA:C8",-79.5036965416));
    this->alphas.insert(std::pair<std::string,double>("C4':GUA:N1",-61.0498362907));
    this->alphas.insert(std::pair<std::string,double>("C4':GUA:N2",-41.6879098573));
    this->alphas.insert(std::pair<std::string,double>("C4':GUA:N3",20.5937728869));
    this->alphas.insert(std::pair<std::string,double>("C4':GUA:N7",9.90036774758));
    this->alphas.insert(std::pair<std::string,double>("C4':GUA:N9",-2.4024367298));
    this->alphas.insert(std::pair<std::string,double>("C4':GUA:O6",66.9450257745));
    this->alphas.insert(std::pair<std::string,double>("C4':ADE:C1'",8.54166699918));
    this->alphas.insert(std::pair<std::string,double>("C4':ADE:C2'",-5.79698822206));
    this->alphas.insert(std::pair<std::string,double>("C4':ADE:C3'",0.895149158984));
    this->alphas.insert(std::pair<std::string,double>("C4':ADE:C4'",7.16229324557));
    this->alphas.insert(std::pair<std::string,double>("C4':ADE:C5'",0.864054599871));
    this->alphas.insert(std::pair<std::string,double>("C4':ADE:P",29.95885311));
    this->alphas.insert(std::pair<std::string,double>("C4':ADE:O5'",-0.487494999503));
    this->alphas.insert(std::pair<std::string,double>("C4':ADE:O3'",3.42045251706));
    this->alphas.insert(std::pair<std::string,double>("C4':ADE:C2",-65.3180175192));
    this->alphas.insert(std::pair<std::string,double>("C4':ADE:C4",3.72334355816));
    this->alphas.insert(std::pair<std::string,double>("C4':ADE:C5",0.102200641083));
    this->alphas.insert(std::pair<std::string,double>("C4':ADE:C6",44.5850882716));
    this->alphas.insert(std::pair<std::string,double>("C4':ADE:C8",-51.4540856351));
    this->alphas.insert(std::pair<std::string,double>("C4':ADE:N1",22.0351469565));
    this->alphas.insert(std::pair<std::string,double>("C4':ADE:N3",-26.7430678812));
    this->alphas.insert(std::pair<std::string,double>("C4':ADE:N6",51.2443732262));
    this->alphas.insert(std::pair<std::string,double>("C4':ADE:N7",31.668157955));
    this->alphas.insert(std::pair<std::string,double>("C4':ADE:N9",-12.8965616378));
    this->alphas.insert(std::pair<std::string,double>("C4':URA:C1'",7.38104945511));
    this->alphas.insert(std::pair<std::string,double>("C4':URA:C2'",2.64752138057));
    this->alphas.insert(std::pair<std::string,double>("C4':URA:C3'",0.154941388675));
    this->alphas.insert(std::pair<std::string,double>("C4':URA:C4'",58.2666686251));
    this->alphas.insert(std::pair<std::string,double>("C4':URA:C5'",-0.162059428093));
    this->alphas.insert(std::pair<std::string,double>("C4':URA:P",33.6088515479));
    this->alphas.insert(std::pair<std::string,double>("C4':URA:O5'",-2.65782323069));
    this->alphas.insert(std::pair<std::string,double>("C4':URA:O3'",28.7847347721));
    this->alphas.insert(std::pair<std::string,double>("C4':URA:C2",-56.8682858838));
    this->alphas.insert(std::pair<std::string,double>("C4':URA:C4",12.9320342444));
    this->alphas.insert(std::pair<std::string,double>("C4':URA:C5",89.0470742678));
    this->alphas.insert(std::pair<std::string,double>("C4':URA:C6",-63.8804979266));
    this->alphas.insert(std::pair<std::string,double>("C4':URA:N1",-42.3510867693));
    this->alphas.insert(std::pair<std::string,double>("C4':URA:N3",-112.279335177));
    this->alphas.insert(std::pair<std::string,double>("C4':URA:O4",-51.0512635834));
    this->alphas.insert(std::pair<std::string,double>("C4':CYT:C1'",1.21707478528));
    this->alphas.insert(std::pair<std::string,double>("C4':CYT:C2'",-1.67697617788));
    this->alphas.insert(std::pair<std::string,double>("C4':CYT:C3'",-0.821268148493));
    this->alphas.insert(std::pair<std::string,double>("C4':CYT:C4'",51.4214069772));
    this->alphas.insert(std::pair<std::string,double>("C4':CYT:C5'",0.231907011312));
    this->alphas.insert(std::pair<std::string,double>("C4':CYT:P",-18.6868420828));
    this->alphas.insert(std::pair<std::string,double>("C4':CYT:O5'",7.55997157808));
    this->alphas.insert(std::pair<std::string,double>("C4':CYT:O3'",1.32967840758));
    this->alphas.insert(std::pair<std::string,double>("C4':CYT:C2",-21.2862655685));
    this->alphas.insert(std::pair<std::string,double>("C4':CYT:C4",-33.1256894794));
    this->alphas.insert(std::pair<std::string,double>("C4':CYT:C5",135.421098181));
    this->alphas.insert(std::pair<std::string,double>("C4':CYT:C6",11.6942300127));
    this->alphas.insert(std::pair<std::string,double>("C4':CYT:N1",3.38143972674));
    this->alphas.insert(std::pair<std::string,double>("C4':CYT:N3",-143.324597908));
    this->alphas.insert(std::pair<std::string,double>("C4':CYT:N4",-84.0763306727));
    this->alphas.insert(std::pair<std::string,double>("C4':CYT:O2",-36.9103183744));
    this->alphas.insert(std::pair<std::string,double>("C2':GUA:C1'",2.56921049828));
    this->alphas.insert(std::pair<std::string,double>("C2':GUA:C2'",-19.4615237082));
    this->alphas.insert(std::pair<std::string,double>("C2':GUA:C3'",-2.55886379337));
    this->alphas.insert(std::pair<std::string,double>("C2':GUA:C4'",2.35609081846));
    this->alphas.insert(std::pair<std::string,double>("C2':GUA:C5'",-14.1194328648));
    this->alphas.insert(std::pair<std::string,double>("C2':GUA:P",0.91469105538));
    this->alphas.insert(std::pair<std::string,double>("C2':GUA:O5'",-12.1478215755));
    this->alphas.insert(std::pair<std::string,double>("C2':GUA:O3'",13.825145053));
    this->alphas.insert(std::pair<std::string,double>("C2':GUA:C2",18.8081097792));
    this->alphas.insert(std::pair<std::string,double>("C2':GUA:C4",-0.600485584917));
    this->alphas.insert(std::pair<std::string,double>("C2':GUA:C5",-2.41395757402));
    this->alphas.insert(std::pair<std::string,double>("C2':GUA:C6",-7.45067598699));
    this->alphas.insert(std::pair<std::string,double>("C2':GUA:C8",-7.84252192089));
    this->alphas.insert(std::pair<std::string,double>("C2':GUA:N1",20.2041435183));
    this->alphas.insert(std::pair<std::string,double>("C2':GUA:N2",24.7018790467));
    this->alphas.insert(std::pair<std::string,double>("C2':GUA:N3",-17.1827741209));
    this->alphas.insert(std::pair<std::string,double>("C2':GUA:N7",-20.1751620649));
    this->alphas.insert(std::pair<std::string,double>("C2':GUA:N9",3.26056131916));
    this->alphas.insert(std::pair<std::string,double>("C2':GUA:O6",-28.2982420406));
    this->alphas.insert(std::pair<std::string,double>("C2':ADE:C1'",0.602941104848));
    this->alphas.insert(std::pair<std::string,double>("C2':ADE:C2'",33.410581859));
    this->alphas.insert(std::pair<std::string,double>("C2':ADE:C3'",1.14500663022));
    this->alphas.insert(std::pair<std::string,double>("C2':ADE:C4'",2.06977973328));
    this->alphas.insert(std::pair<std::string,double>("C2':ADE:C5'",0.279592919088));
    this->alphas.insert(std::pair<std::string,double>("C2':ADE:P",-34.0469456495));
    this->alphas.insert(std::pair<std::string,double>("C2':ADE:O5'",-11.3333487803));
    this->alphas.insert(std::pair<std::string,double>("C2':ADE:O3'",8.16840735035));
    this->alphas.insert(std::pair<std::string,double>("C2':ADE:C2",-7.70162427119));
    this->alphas.insert(std::pair<std::string,double>("C2':ADE:C4",1.31551830428));
    this->alphas.insert(std::pair<std::string,double>("C2':ADE:C5",-7.36437922948));
    this->alphas.insert(std::pair<std::string,double>("C2':ADE:C6",-4.77353015653));
    this->alphas.insert(std::pair<std::string,double>("C2':ADE:C8",-2.31940645258));
    this->alphas.insert(std::pair<std::string,double>("C2':ADE:N1",-7.55219557282));
    this->alphas.insert(std::pair<std::string,double>("C2':ADE:N3",-15.9979176294));
    this->alphas.insert(std::pair<std::string,double>("C2':ADE:N6",22.1344871732));
    this->alphas.insert(std::pair<std::string,double>("C2':ADE:N7",-2.05059128996));
    this->alphas.insert(std::pair<std::string,double>("C2':ADE:N9",3.04580838472));
    this->alphas.insert(std::pair<std::string,double>("C2':URA:C1'",-0.011087024962));
    this->alphas.insert(std::pair<std::string,double>("C2':URA:C2'",68.0096523683));
    this->alphas.insert(std::pair<std::string,double>("C2':URA:C3'",0.317643844016));
    this->alphas.insert(std::pair<std::string,double>("C2':URA:C4'",1.26246119176));
    this->alphas.insert(std::pair<std::string,double>("C2':URA:C5'",-12.9852561071));
    this->alphas.insert(std::pair<std::string,double>("C2':URA:P",-22.4884946359));
    this->alphas.insert(std::pair<std::string,double>("C2':URA:O5'",-26.8758043936));
    this->alphas.insert(std::pair<std::string,double>("C2':URA:O3'",29.9291837084));
    this->alphas.insert(std::pair<std::string,double>("C2':URA:C2",-0.734971911213));
    this->alphas.insert(std::pair<std::string,double>("C2':URA:C4",-19.7395960619));
    this->alphas.insert(std::pair<std::string,double>("C2':URA:C5",32.6151150471));
    this->alphas.insert(std::pair<std::string,double>("C2':URA:C6",-25.4846471712));
    this->alphas.insert(std::pair<std::string,double>("C2':URA:N1",-3.47800665755));
    this->alphas.insert(std::pair<std::string,double>("C2':URA:N3",-30.4374466821));
    this->alphas.insert(std::pair<std::string,double>("C2':URA:O4",-3.29918848118));
    this->alphas.insert(std::pair<std::string,double>("C2':CYT:C1'",0.659442759842));
    this->alphas.insert(std::pair<std::string,double>("C2':CYT:C2'",39.0671676325));
    this->alphas.insert(std::pair<std::string,double>("C2':CYT:C3'",0.430309211905));
    this->alphas.insert(std::pair<std::string,double>("C2':CYT:C4'",3.77395532579));
    this->alphas.insert(std::pair<std::string,double>("C2':CYT:C5'",4.15690325159));
    this->alphas.insert(std::pair<std::string,double>("C2':CYT:P",-30.3347493869));
    this->alphas.insert(std::pair<std::string,double>("C2':CYT:O5'",-9.07279503151));
    this->alphas.insert(std::pair<std::string,double>("C2':CYT:O3'",16.2007194464));
    this->alphas.insert(std::pair<std::string,double>("C2':CYT:C2",-1.09775035107));
    this->alphas.insert(std::pair<std::string,double>("C2':CYT:C4",-3.47225855282));
    this->alphas.insert(std::pair<std::string,double>("C2':CYT:C5",-46.344005682));
    this->alphas.insert(std::pair<std::string,double>("C2':CYT:C6",-7.36817278071));
    this->alphas.insert(std::pair<std::string,double>("C2':CYT:N1",1.43600565073));
    this->alphas.insert(std::pair<std::string,double>("C2':CYT:N3",-12.6599239566));
    this->alphas.insert(std::pair<std::string,double>("C2':CYT:N4",57.3736531401));
    this->alphas.insert(std::pair<std::string,double>("C2':CYT:O2",-15.603598269));
    this->alphas.insert(std::pair<std::string,double>("N1:GUA:N7",-60.9541652271));
    this->alphas.insert(std::pair<std::string,double>("N1:GUA:C2",-0.754297767631));
    this->alphas.insert(std::pair<std::string,double>("N1:GUA:C4'",-71.4258300201));
    this->alphas.insert(std::pair<std::string,double>("N1:GUA:N9",-9.09889937229));
    this->alphas.insert(std::pair<std::string,double>("N1:GUA:O5'",140.625554782));
    this->alphas.insert(std::pair<std::string,double>("N1:GUA:C3'",33.0880670226));
    this->alphas.insert(std::pair<std::string,double>("N1:GUA:C1'",78.8131712947));
    this->alphas.insert(std::pair<std::string,double>("N1:GUA:C5",4.35280475387));
    this->alphas.insert(std::pair<std::string,double>("N1:GUA:C5'",-162.554182881));
    this->alphas.insert(std::pair<std::string,double>("N1:GUA:P",113.986552671));
    this->alphas.insert(std::pair<std::string,double>("N1:GUA:C8",-87.3441084355));
    this->alphas.insert(std::pair<std::string,double>("N1:GUA:C2'",21.3319027006));
    this->alphas.insert(std::pair<std::string,double>("N1:GUA:O6",2.12352800307));
    this->alphas.insert(std::pair<std::string,double>("N1:GUA:N1",2.41975323989));
    this->alphas.insert(std::pair<std::string,double>("N1:GUA:N2",0.262076632023));
    this->alphas.insert(std::pair<std::string,double>("N1:GUA:N3",0.00293001430038));
    this->alphas.insert(std::pair<std::string,double>("N1:GUA:C6",0.540220402937));
    this->alphas.insert(std::pair<std::string,double>("N1:GUA:O3'",-58.9629654665));
    this->alphas.insert(std::pair<std::string,double>("N1:GUA:C4",0.622051087229));
    this->alphas.insert(std::pair<std::string,double>("N1:ADE:N7",7.29486589511));
    this->alphas.insert(std::pair<std::string,double>("N1:ADE:N6",35.5809756633));
    this->alphas.insert(std::pair<std::string,double>("N1:ADE:N9",-53.1021041115));
    this->alphas.insert(std::pair<std::string,double>("N1:ADE:C1'",-65.9718112612));
    this->alphas.insert(std::pair<std::string,double>("N1:ADE:C3'",70.3869855475));
    this->alphas.insert(std::pair<std::string,double>("N1:ADE:O5'",-8.52794185159));
    this->alphas.insert(std::pair<std::string,double>("N1:ADE:C5",6.03840316685));
    this->alphas.insert(std::pair<std::string,double>("N1:ADE:C5'",63.2995970145));
    this->alphas.insert(std::pair<std::string,double>("N1:ADE:P",4.14678459581));
    this->alphas.insert(std::pair<std::string,double>("N1:ADE:C4'",6.48882075418));
    this->alphas.insert(std::pair<std::string,double>("N1:ADE:C2'",37.5887061074));
    this->alphas.insert(std::pair<std::string,double>("N1:ADE:N1",-29.1317483515));
    this->alphas.insert(std::pair<std::string,double>("N1:ADE:C2",-24.1365784558));
    this->alphas.insert(std::pair<std::string,double>("N1:ADE:C8",-9.28688096445));
    this->alphas.insert(std::pair<std::string,double>("N1:ADE:N3",14.4036270493));
    this->alphas.insert(std::pair<std::string,double>("N1:ADE:C6",11.8280366599));
    this->alphas.insert(std::pair<std::string,double>("N1:ADE:O3'",-11.8485086423));
    this->alphas.insert(std::pair<std::string,double>("N1:ADE:C4",-16.7723751035));
    this->alphas.insert(std::pair<std::string,double>("N1:URA:O4",-29.5254573925));
    this->alphas.insert(std::pair<std::string,double>("N1:URA:O5'",107.971883488));
    this->alphas.insert(std::pair<std::string,double>("N1:URA:C3'",-70.6614307649));
    this->alphas.insert(std::pair<std::string,double>("N1:URA:O3'",-22.795837646));
    this->alphas.insert(std::pair<std::string,double>("N1:URA:C1'",63.8406338391));
    this->alphas.insert(std::pair<std::string,double>("N1:URA:C5'",20.8968068256));
    this->alphas.insert(std::pair<std::string,double>("N1:URA:P",59.6075179566));
    this->alphas.insert(std::pair<std::string,double>("N1:URA:C4'",113.791188942));
    this->alphas.insert(std::pair<std::string,double>("N1:URA:C2'",-107.902262722));
    this->alphas.insert(std::pair<std::string,double>("N1:URA:N1",-27.0485507383));
    this->alphas.insert(std::pair<std::string,double>("N1:URA:C2",19.1323928576));
    this->alphas.insert(std::pair<std::string,double>("N1:URA:N3",46.0535506017));
    this->alphas.insert(std::pair<std::string,double>("N1:URA:C6",-127.652849914));
    this->alphas.insert(std::pair<std::string,double>("N1:URA:C5",-48.9498074799));
    this->alphas.insert(std::pair<std::string,double>("N1:URA:C4",49.8764830447));
    this->alphas.insert(std::pair<std::string,double>("N1:CYT:O5'",89.343230742));
    this->alphas.insert(std::pair<std::string,double>("N1:CYT:C3'",-17.455325413));
    this->alphas.insert(std::pair<std::string,double>("N1:CYT:O3'",95.5097221533));
    this->alphas.insert(std::pair<std::string,double>("N1:CYT:C1'",-88.2322742213));
    this->alphas.insert(std::pair<std::string,double>("N1:CYT:O2",27.4454321209));
    this->alphas.insert(std::pair<std::string,double>("N1:CYT:C5'",-31.9647652047));
    this->alphas.insert(std::pair<std::string,double>("N1:CYT:P",45.9856156017));
    this->alphas.insert(std::pair<std::string,double>("N1:CYT:C4'",11.6047433939));
    this->alphas.insert(std::pair<std::string,double>("N1:CYT:C2'",-54.7501251966));
    this->alphas.insert(std::pair<std::string,double>("N1:CYT:N1",-15.7354433186));
    this->alphas.insert(std::pair<std::string,double>("N1:CYT:C2",33.6150942453));
    this->alphas.insert(std::pair<std::string,double>("N1:CYT:N3",31.9686989868));
    this->alphas.insert(std::pair<std::string,double>("N1:CYT:N4",-33.8448398728));
    this->alphas.insert(std::pair<std::string,double>("N1:CYT:C6",7.98261343783));
    this->alphas.insert(std::pair<std::string,double>("N1:CYT:C5",6.7847540724));
    this->alphas.insert(std::pair<std::string,double>("N1:CYT:C4",24.2090811976));
    this->alphas.insert(std::pair<std::string,double>("H2':GUA:C1'",0.443213913089));
    this->alphas.insert(std::pair<std::string,double>("H2':GUA:C2'",-0.00109234454735));
    this->alphas.insert(std::pair<std::string,double>("H2':GUA:C3'",1.22756367961));
    this->alphas.insert(std::pair<std::string,double>("H2':GUA:C4'",1.02783065306));
    this->alphas.insert(std::pair<std::string,double>("H2':GUA:C5'",5.03768899238));
    this->alphas.insert(std::pair<std::string,double>("H2':GUA:P",2.09032241937));
    this->alphas.insert(std::pair<std::string,double>("H2':GUA:O5'",-0.976924815948));
    this->alphas.insert(std::pair<std::string,double>("H2':GUA:O3'",-1.59206906818));
    this->alphas.insert(std::pair<std::string,double>("H2':GUA:C2",0.549545469052));
    this->alphas.insert(std::pair<std::string,double>("H2':GUA:C4",-0.342907679656));
    this->alphas.insert(std::pair<std::string,double>("H2':GUA:C5",-0.213933165529));
    this->alphas.insert(std::pair<std::string,double>("H2':GUA:C6",-0.535284229993));
    this->alphas.insert(std::pair<std::string,double>("H2':GUA:C8",1.9560296341));
    this->alphas.insert(std::pair<std::string,double>("H2':GUA:N1",-20.3271232254));
    this->alphas.insert(std::pair<std::string,double>("H2':GUA:N2",3.04349846521));
    this->alphas.insert(std::pair<std::string,double>("H2':GUA:N3",5.28584984899));
    this->alphas.insert(std::pair<std::string,double>("H2':GUA:N7",0.199316154869));
    this->alphas.insert(std::pair<std::string,double>("H2':GUA:N9",1.03234348523));
    this->alphas.insert(std::pair<std::string,double>("H2':GUA:O6",9.35543099661));
    this->alphas.insert(std::pair<std::string,double>("H2':ADE:C1'",-0.102701478619));
    this->alphas.insert(std::pair<std::string,double>("H2':ADE:C2'",0.289431267668));
    this->alphas.insert(std::pair<std::string,double>("H2':ADE:C3'",0.839207071517));
    this->alphas.insert(std::pair<std::string,double>("H2':ADE:C4'",2.25206233142));
    this->alphas.insert(std::pair<std::string,double>("H2':ADE:C5'",1.97864562812));
    this->alphas.insert(std::pair<std::string,double>("H2':ADE:P",1.20547544542));
    this->alphas.insert(std::pair<std::string,double>("H2':ADE:O5'",0.539669033025));
    this->alphas.insert(std::pair<std::string,double>("H2':ADE:O3'",-0.559613204168));
    this->alphas.insert(std::pair<std::string,double>("H2':ADE:C2",-2.85210822007));
    this->alphas.insert(std::pair<std::string,double>("H2':ADE:C4",-0.127242057563));
    this->alphas.insert(std::pair<std::string,double>("H2':ADE:C5",-3.55121297143));
    this->alphas.insert(std::pair<std::string,double>("H2':ADE:C6",-6.83341172749));
    this->alphas.insert(std::pair<std::string,double>("H2':ADE:C8",2.52383001923));
    this->alphas.insert(std::pair<std::string,double>("H2':ADE:N1",1.08311829787));
    this->alphas.insert(std::pair<std::string,double>("H2':ADE:N3",0.92975307245));
    this->alphas.insert(std::pair<std::string,double>("H2':ADE:N6",-0.287897152995));
    this->alphas.insert(std::pair<std::string,double>("H2':ADE:N7",2.10462660808));
    this->alphas.insert(std::pair<std::string,double>("H2':ADE:N9",-0.14839311818));
    this->alphas.insert(std::pair<std::string,double>("H2':URA:C1'",1.54974921579));
    this->alphas.insert(std::pair<std::string,double>("H2':URA:C2'",-0.0130021000405));
    this->alphas.insert(std::pair<std::string,double>("H2':URA:C3'",-0.052255551788));
    this->alphas.insert(std::pair<std::string,double>("H2':URA:C4'",1.6742177488));
    this->alphas.insert(std::pair<std::string,double>("H2':URA:C5'",0.917824138431));
    this->alphas.insert(std::pair<std::string,double>("H2':URA:P",5.1099478725));
    this->alphas.insert(std::pair<std::string,double>("H2':URA:O5'",-1.92484224319));
    this->alphas.insert(std::pair<std::string,double>("H2':URA:O3'",1.30455722532));
    this->alphas.insert(std::pair<std::string,double>("H2':URA:C2",0.454503591558));
    this->alphas.insert(std::pair<std::string,double>("H2':URA:C4",-3.19148236974));
    this->alphas.insert(std::pair<std::string,double>("H2':URA:C5",4.86517666667));
    this->alphas.insert(std::pair<std::string,double>("H2':URA:C6",-0.324998742596));
    this->alphas.insert(std::pair<std::string,double>("H2':URA:N1",-0.295025124752));
    this->alphas.insert(std::pair<std::string,double>("H2':URA:N3",-4.31179803052));
    this->alphas.insert(std::pair<std::string,double>("H2':URA:O4",2.46710999318));
    this->alphas.insert(std::pair<std::string,double>("H2':CYT:C1'",-0.000363426973778));
    this->alphas.insert(std::pair<std::string,double>("H2':CYT:C2'",-0.0));
    this->alphas.insert(std::pair<std::string,double>("H2':CYT:C3'",-0.103864924564));
    this->alphas.insert(std::pair<std::string,double>("H2':CYT:C4'",-0.425476300534));
    this->alphas.insert(std::pair<std::string,double>("H2':CYT:C5'",-1.66067411242));
    this->alphas.insert(std::pair<std::string,double>("H2':CYT:P",8.17750077138));
    this->alphas.insert(std::pair<std::string,double>("H2':CYT:O5'",0.814616687963));
    this->alphas.insert(std::pair<std::string,double>("H2':CYT:O3'",-0.262832618758));
    this->alphas.insert(std::pair<std::string,double>("H2':CYT:C2",3.56451512396));
    this->alphas.insert(std::pair<std::string,double>("H2':CYT:C4",-2.21568840444));
    this->alphas.insert(std::pair<std::string,double>("H2':CYT:C5",0.12297745991));
    this->alphas.insert(std::pair<std::string,double>("H2':CYT:C6",0.199977535576));
    this->alphas.insert(std::pair<std::string,double>("H2':CYT:N1",1.12189381095));
    this->alphas.insert(std::pair<std::string,double>("H2':CYT:N3",-0.240168309961));
    this->alphas.insert(std::pair<std::string,double>("H2':CYT:N4",-9.90396375057));
    this->alphas.insert(std::pair<std::string,double>("H2':CYT:O2",1.51301098691));
    this->alphas.insert(std::pair<std::string,double>("N3:GUA:N7",6.58458596584));
    this->alphas.insert(std::pair<std::string,double>("N3:GUA:C2",-1.020229594));
    this->alphas.insert(std::pair<std::string,double>("N3:GUA:C4'",13.3695838149));
    this->alphas.insert(std::pair<std::string,double>("N3:GUA:N9",13.6851718278));
    this->alphas.insert(std::pair<std::string,double>("N3:GUA:O5'",-21.9209550798));
    this->alphas.insert(std::pair<std::string,double>("N3:GUA:C3'",151.992664753));
    this->alphas.insert(std::pair<std::string,double>("N3:GUA:C1'",74.6879487798));
    this->alphas.insert(std::pair<std::string,double>("N3:GUA:C5",10.3167521438));
    this->alphas.insert(std::pair<std::string,double>("N3:GUA:C5'",-86.8145406829));
    this->alphas.insert(std::pair<std::string,double>("N3:GUA:P",-156.918066634));
    this->alphas.insert(std::pair<std::string,double>("N3:GUA:C8",-3.22806610962));
    this->alphas.insert(std::pair<std::string,double>("N3:GUA:C2'",244.481297931));
    this->alphas.insert(std::pair<std::string,double>("N3:GUA:O6",-45.4255244858));
    this->alphas.insert(std::pair<std::string,double>("N3:GUA:N1",-10.1119544361));
    this->alphas.insert(std::pair<std::string,double>("N3:GUA:N2",48.3369135141));
    this->alphas.insert(std::pair<std::string,double>("N3:GUA:N3",-27.6187012117));
    this->alphas.insert(std::pair<std::string,double>("N3:GUA:C6",-22.0260300556));
    this->alphas.insert(std::pair<std::string,double>("N3:GUA:O3'",89.5336193189));
    this->alphas.insert(std::pair<std::string,double>("N3:GUA:C4",-9.07730721415));
    this->alphas.insert(std::pair<std::string,double>("N3:ADE:N7",-48.562077098));
    this->alphas.insert(std::pair<std::string,double>("N3:ADE:N6",-81.5183649793));
    this->alphas.insert(std::pair<std::string,double>("N3:ADE:N9",66.5560541678));
    this->alphas.insert(std::pair<std::string,double>("N3:ADE:C1'",-12.4726805814));
    this->alphas.insert(std::pair<std::string,double>("N3:ADE:C3'",-110.854262208));
    this->alphas.insert(std::pair<std::string,double>("N3:ADE:O5'",150.316169826));
    this->alphas.insert(std::pair<std::string,double>("N3:ADE:C5",34.6729195227));
    this->alphas.insert(std::pair<std::string,double>("N3:ADE:C5'",-24.5349238782));
    this->alphas.insert(std::pair<std::string,double>("N3:ADE:P",38.0076953551));
    this->alphas.insert(std::pair<std::string,double>("N3:ADE:C4'",-22.7069365658));
    this->alphas.insert(std::pair<std::string,double>("N3:ADE:C2'",-128.562350308));
    this->alphas.insert(std::pair<std::string,double>("N3:ADE:N1",54.5039553683));
    this->alphas.insert(std::pair<std::string,double>("N3:ADE:C2",2.57298688532));
    this->alphas.insert(std::pair<std::string,double>("N3:ADE:C8",4.42282476843));
    this->alphas.insert(std::pair<std::string,double>("N3:ADE:N3",18.0879831343));
    this->alphas.insert(std::pair<std::string,double>("N3:ADE:C6",-2.12990387137));
    this->alphas.insert(std::pair<std::string,double>("N3:ADE:O3'",-90.3587295911));
    this->alphas.insert(std::pair<std::string,double>("N3:ADE:C4",100.756259226));
    this->alphas.insert(std::pair<std::string,double>("N3:URA:O4",8.74121667938));
    this->alphas.insert(std::pair<std::string,double>("N3:URA:O5'",187.450114456));
    this->alphas.insert(std::pair<std::string,double>("N3:URA:C3'",-59.959167756));
    this->alphas.insert(std::pair<std::string,double>("N3:URA:O3'",123.48063276));
    this->alphas.insert(std::pair<std::string,double>("N3:URA:C1'",-0.667365340646));
    this->alphas.insert(std::pair<std::string,double>("N3:URA:C5'",-119.685169009));
    this->alphas.insert(std::pair<std::string,double>("N3:URA:P",93.8217539026));
    this->alphas.insert(std::pair<std::string,double>("N3:URA:C4'",-10.3095050415));
    this->alphas.insert(std::pair<std::string,double>("N3:URA:C2'",-12.4752420875));
    this->alphas.insert(std::pair<std::string,double>("N3:URA:N1",0.24318638595));
    this->alphas.insert(std::pair<std::string,double>("N3:URA:C2",-2.87705813296));
    this->alphas.insert(std::pair<std::string,double>("N3:URA:N3",36.8388488999));
    this->alphas.insert(std::pair<std::string,double>("N3:URA:C6",-23.019933208));
    this->alphas.insert(std::pair<std::string,double>("N3:URA:C5",-13.7796780268));
    this->alphas.insert(std::pair<std::string,double>("N3:URA:C4",2.33458669607));
    this->alphas.insert(std::pair<std::string,double>("N3:CYT:O5'",-19.6336193154));
    this->alphas.insert(std::pair<std::string,double>("N3:CYT:C3'",-47.9601595526));
    this->alphas.insert(std::pair<std::string,double>("N3:CYT:O3'",-82.5436265534));
    this->alphas.insert(std::pair<std::string,double>("N3:CYT:C1'",47.1833039321));
    this->alphas.insert(std::pair<std::string,double>("N3:CYT:O2",154.02013224));
    this->alphas.insert(std::pair<std::string,double>("N3:CYT:C5'",57.7609729657));
    this->alphas.insert(std::pair<std::string,double>("N3:CYT:P",0.510190378477));
    this->alphas.insert(std::pair<std::string,double>("N3:CYT:C4'",61.5831729679));
    this->alphas.insert(std::pair<std::string,double>("N3:CYT:C2'",-2.09906095257));
    this->alphas.insert(std::pair<std::string,double>("N3:CYT:N1",-41.4674631272));
    this->alphas.insert(std::pair<std::string,double>("N3:CYT:C2",17.1776008859));
    this->alphas.insert(std::pair<std::string,double>("N3:CYT:N3",-30.2446420309));
    this->alphas.insert(std::pair<std::string,double>("N3:CYT:N4",-49.2058116702));
    this->alphas.insert(std::pair<std::string,double>("N3:CYT:C6",-47.9050413804));
    this->alphas.insert(std::pair<std::string,double>("N3:CYT:C5",32.9805864051));
    this->alphas.insert(std::pair<std::string,double>("N3:CYT:C4",33.8260380518));
    this->alphas.insert(std::pair<std::string,double>("H1':GUA:C1'",0.428401203243));
    this->alphas.insert(std::pair<std::string,double>("H1':GUA:C2'",0.982419947699));
    this->alphas.insert(std::pair<std::string,double>("H1':GUA:C3'",2.76813085461));
    this->alphas.insert(std::pair<std::string,double>("H1':GUA:C4'",1.73102062865));
    this->alphas.insert(std::pair<std::string,double>("H1':GUA:C5'",3.71722447508));
    this->alphas.insert(std::pair<std::string,double>("H1':GUA:P",-1.12256934598));
    this->alphas.insert(std::pair<std::string,double>("H1':GUA:O5'",-0.678871448769));
    this->alphas.insert(std::pair<std::string,double>("H1':GUA:O3'",1.98327820674));
    this->alphas.insert(std::pair<std::string,double>("H1':GUA:C2",2.6954631793));
    this->alphas.insert(std::pair<std::string,double>("H1':GUA:C4",0.171315095916));
    this->alphas.insert(std::pair<std::string,double>("H1':GUA:C5",-15.3453413609));
    this->alphas.insert(std::pair<std::string,double>("H1':GUA:C6",-8.21106569206));
    this->alphas.insert(std::pair<std::string,double>("H1':GUA:C8",4.55036045992));
    this->alphas.insert(std::pair<std::string,double>("H1':GUA:N1",2.3550934953));
    this->alphas.insert(std::pair<std::string,double>("H1':GUA:N2",0.744104022029));
    this->alphas.insert(std::pair<std::string,double>("H1':GUA:N3",0.83095222479));
    this->alphas.insert(std::pair<std::string,double>("H1':GUA:N7",-5.89263198531));
    this->alphas.insert(std::pair<std::string,double>("H1':GUA:N9",1.91857659357));
    this->alphas.insert(std::pair<std::string,double>("H1':GUA:O6",6.94542079849));
    this->alphas.insert(std::pair<std::string,double>("H1':ADE:C1'",1.51229359219));
    this->alphas.insert(std::pair<std::string,double>("H1':ADE:C2'",1.72348040483));
    this->alphas.insert(std::pair<std::string,double>("H1':ADE:C3'",0.288474060452));
    this->alphas.insert(std::pair<std::string,double>("H1':ADE:C4'",4.27459478411));
    this->alphas.insert(std::pair<std::string,double>("H1':ADE:C5'",-10.2068987335));
    this->alphas.insert(std::pair<std::string,double>("H1':ADE:P",4.16606990943));
    this->alphas.insert(std::pair<std::string,double>("H1':ADE:O5'",-0.850448133504));
    this->alphas.insert(std::pair<std::string,double>("H1':ADE:O3'",-2.88093697619));
    this->alphas.insert(std::pair<std::string,double>("H1':ADE:C2",1.68212803064));
    this->alphas.insert(std::pair<std::string,double>("H1':ADE:C4",0.390061499576));
    this->alphas.insert(std::pair<std::string,double>("H1':ADE:C5",-28.5950918785));
    this->alphas.insert(std::pair<std::string,double>("H1':ADE:C6",-13.9985976343));
    this->alphas.insert(std::pair<std::string,double>("H1':ADE:C8",2.63612395693));
    this->alphas.insert(std::pair<std::string,double>("H1':ADE:N1",6.41895773111));
    this->alphas.insert(std::pair<std::string,double>("H1':ADE:N3",-0.374343440468));
    this->alphas.insert(std::pair<std::string,double>("H1':ADE:N6",9.45033832712));
    this->alphas.insert(std::pair<std::string,double>("H1':ADE:N7",1.58825189805));
    this->alphas.insert(std::pair<std::string,double>("H1':ADE:N9",-0.0486718470351));
    this->alphas.insert(std::pair<std::string,double>("H1':URA:C1'",-0.0));
    this->alphas.insert(std::pair<std::string,double>("H1':URA:C2'",-0.189115925066));
    this->alphas.insert(std::pair<std::string,double>("H1':URA:C3'",1.41917076395));
    this->alphas.insert(std::pair<std::string,double>("H1':URA:C4'",-0.00304489198373));
    this->alphas.insert(std::pair<std::string,double>("H1':URA:C5'",-1.89108479616));
    this->alphas.insert(std::pair<std::string,double>("H1':URA:P",-3.06345113005));
    this->alphas.insert(std::pair<std::string,double>("H1':URA:O5'",-11.7629160567));
    this->alphas.insert(std::pair<std::string,double>("H1':URA:O3'",14.5957500162));
    this->alphas.insert(std::pair<std::string,double>("H1':URA:C2",-0.184076072202));
    this->alphas.insert(std::pair<std::string,double>("H1':URA:C4",-3.4364381447));
    this->alphas.insert(std::pair<std::string,double>("H1':URA:C5",9.13129336492));
    this->alphas.insert(std::pair<std::string,double>("H1':URA:C6",-3.26480276004));
    this->alphas.insert(std::pair<std::string,double>("H1':URA:N1",-0.467575736342));
    this->alphas.insert(std::pair<std::string,double>("H1':URA:N3",5.29094212663));
    this->alphas.insert(std::pair<std::string,double>("H1':URA:O4",-2.73269759854));
    this->alphas.insert(std::pair<std::string,double>("H1':CYT:C1'",0.0223761113264));
    this->alphas.insert(std::pair<std::string,double>("H1':CYT:C2'",0.468384222507));
    this->alphas.insert(std::pair<std::string,double>("H1':CYT:C3'",2.66350021597));
    this->alphas.insert(std::pair<std::string,double>("H1':CYT:C4'",2.52513704899));
    this->alphas.insert(std::pair<std::string,double>("H1':CYT:C5'",-1.89097112757));
    this->alphas.insert(std::pair<std::string,double>("H1':CYT:P",-7.1344290962));
    this->alphas.insert(std::pair<std::string,double>("H1':CYT:O5'",-6.3274276126));
    this->alphas.insert(std::pair<std::string,double>("H1':CYT:O3'",4.93832856133));
    this->alphas.insert(std::pair<std::string,double>("H1':CYT:C2",0.403091552837));
    this->alphas.insert(std::pair<std::string,double>("H1':CYT:C4",-2.34344412576));
    this->alphas.insert(std::pair<std::string,double>("H1':CYT:C5",2.16880387864));
    this->alphas.insert(std::pair<std::string,double>("H1':CYT:C6",2.26705259774));
    this->alphas.insert(std::pair<std::string,double>("H1':CYT:N1",1.35760636487));
    this->alphas.insert(std::pair<std::string,double>("H1':CYT:N3",-7.07094179203));
    this->alphas.insert(std::pair<std::string,double>("H1':CYT:N4",3.43108123584));
    this->alphas.insert(std::pair<std::string,double>("H1':CYT:O2",0.453256398308));
    this->alphas.insert(std::pair<std::string,double>("H4':GUA:C1'",0.709970102502));
    this->alphas.insert(std::pair<std::string,double>("H4':GUA:C2'",1.75726958963));
    this->alphas.insert(std::pair<std::string,double>("H4':GUA:C3'",0.820067830994));
    this->alphas.insert(std::pair<std::string,double>("H4':GUA:C4'",-0.0));
    this->alphas.insert(std::pair<std::string,double>("H4':GUA:C5'",1.19834177254));
    this->alphas.insert(std::pair<std::string,double>("H4':GUA:P",1.41704358137));
    this->alphas.insert(std::pair<std::string,double>("H4':GUA:O5'",0.214611193547));
    this->alphas.insert(std::pair<std::string,double>("H4':GUA:O3'",2.55098124965));
    this->alphas.insert(std::pair<std::string,double>("H4':GUA:C2",-4.34378277938));
    this->alphas.insert(std::pair<std::string,double>("H4':GUA:C4",-13.0612443102));
    this->alphas.insert(std::pair<std::string,double>("H4':GUA:C5",-5.82963520485));
    this->alphas.insert(std::pair<std::string,double>("H4':GUA:C6",-9.0166408166));
    this->alphas.insert(std::pair<std::string,double>("H4':GUA:C8",7.36187099));
    this->alphas.insert(std::pair<std::string,double>("H4':GUA:N1",-1.57108231666));
    this->alphas.insert(std::pair<std::string,double>("H4':GUA:N2",15.0809389007));
    this->alphas.insert(std::pair<std::string,double>("H4':GUA:N3",-6.01448166863));
    this->alphas.insert(std::pair<std::string,double>("H4':GUA:N7",3.39001108825));
    this->alphas.insert(std::pair<std::string,double>("H4':GUA:N9",0.122165194123));
    this->alphas.insert(std::pair<std::string,double>("H4':GUA:O6",-1.94562554953));
    this->alphas.insert(std::pair<std::string,double>("H4':ADE:C1'",0.78744848204));
    this->alphas.insert(std::pair<std::string,double>("H4':ADE:C2'",0.561501546578));
    this->alphas.insert(std::pair<std::string,double>("H4':ADE:C3'",1.64313775127));
    this->alphas.insert(std::pair<std::string,double>("H4':ADE:C4'",0.300548491127));
    this->alphas.insert(std::pair<std::string,double>("H4':ADE:C5'",0.294734414741));
    this->alphas.insert(std::pair<std::string,double>("H4':ADE:P",0.0159862029319));
    this->alphas.insert(std::pair<std::string,double>("H4':ADE:O5'",0.0995858064853));
    this->alphas.insert(std::pair<std::string,double>("H4':ADE:O3'",0.506875271718));
    this->alphas.insert(std::pair<std::string,double>("H4':ADE:C2",2.93010658659));
    this->alphas.insert(std::pair<std::string,double>("H4':ADE:C4",-5.17376747269));
    this->alphas.insert(std::pair<std::string,double>("H4':ADE:C5",-0.480965062922));
    this->alphas.insert(std::pair<std::string,double>("H4':ADE:C6",1.2166392753));
    this->alphas.insert(std::pair<std::string,double>("H4':ADE:C8",-2.10344066894));
    this->alphas.insert(std::pair<std::string,double>("H4':ADE:N1",2.82222793768));
    this->alphas.insert(std::pair<std::string,double>("H4':ADE:N3",-4.72870593571));
    this->alphas.insert(std::pair<std::string,double>("H4':ADE:N6",-1.83851947014));
    this->alphas.insert(std::pair<std::string,double>("H4':ADE:N7",-3.0654907133));
    this->alphas.insert(std::pair<std::string,double>("H4':ADE:N9",-0.347280799401));
    this->alphas.insert(std::pair<std::string,double>("H4':URA:C1'",-0.0734787771009));
    this->alphas.insert(std::pair<std::string,double>("H4':URA:C2'",-0.0418791891531));
    this->alphas.insert(std::pair<std::string,double>("H4':URA:C3'",4.4195603373));
    this->alphas.insert(std::pair<std::string,double>("H4':URA:C4'",-0.0));
    this->alphas.insert(std::pair<std::string,double>("H4':URA:C5'",-0.758779742566));
    this->alphas.insert(std::pair<std::string,double>("H4':URA:P",1.88158819111));
    this->alphas.insert(std::pair<std::string,double>("H4':URA:O5'",-0.159475685151));
    this->alphas.insert(std::pair<std::string,double>("H4':URA:O3'",-0.66012676606));
    this->alphas.insert(std::pair<std::string,double>("H4':URA:C2",-0.898713697159));
    this->alphas.insert(std::pair<std::string,double>("H4':URA:C4",1.90197590148));
    this->alphas.insert(std::pair<std::string,double>("H4':URA:C5",-5.09621978092));
    this->alphas.insert(std::pair<std::string,double>("H4':URA:C6",2.88734341205));
    this->alphas.insert(std::pair<std::string,double>("H4':URA:N1",-0.238694211403));
    this->alphas.insert(std::pair<std::string,double>("H4':URA:N3",6.04312923997));
    this->alphas.insert(std::pair<std::string,double>("H4':URA:O4",0.568652667817));
    this->alphas.insert(std::pair<std::string,double>("H4':CYT:C1'",1.60494866062));
    this->alphas.insert(std::pair<std::string,double>("H4':CYT:C2'",1.39553501382));
    this->alphas.insert(std::pair<std::string,double>("H4':CYT:C3'",-0.0825747696457));
    this->alphas.insert(std::pair<std::string,double>("H4':CYT:C4'",-0.0));
    this->alphas.insert(std::pair<std::string,double>("H4':CYT:C5'",2.59564025402));
    this->alphas.insert(std::pair<std::string,double>("H4':CYT:P",1.62848792608));
    this->alphas.insert(std::pair<std::string,double>("H4':CYT:O5'",-1.42189755891));
    this->alphas.insert(std::pair<std::string,double>("H4':CYT:O3'",-0.310238368481));
    this->alphas.insert(std::pair<std::string,double>("H4':CYT:C2",-4.88594277688));
    this->alphas.insert(std::pair<std::string,double>("H4':CYT:C4",0.612360383042));
    this->alphas.insert(std::pair<std::string,double>("H4':CYT:C5",-7.61626186985));
    this->alphas.insert(std::pair<std::string,double>("H4':CYT:C6",-1.44130880698));
    this->alphas.insert(std::pair<std::string,double>("H4':CYT:N1",2.10265014085));
    this->alphas.insert(std::pair<std::string,double>("H4':CYT:N3",3.01762359772));
    this->alphas.insert(std::pair<std::string,double>("H4':CYT:N4",20.7989807661));
    this->alphas.insert(std::pair<std::string,double>("H4':CYT:O2",-4.12741813733));
    this->alphas.insert(std::pair<std::string,double>("H8:GUA:C1'",4.00482565816));
    this->alphas.insert(std::pair<std::string,double>("H8:GUA:C2'",-8.12091216026));
    this->alphas.insert(std::pair<std::string,double>("H8:GUA:C3'",-14.3896391287));
    this->alphas.insert(std::pair<std::string,double>("H8:GUA:C4'",2.31794287133));
    this->alphas.insert(std::pair<std::string,double>("H8:GUA:C5'",8.16037339359));
    this->alphas.insert(std::pair<std::string,double>("H8:GUA:P",3.45533761065));
    this->alphas.insert(std::pair<std::string,double>("H8:GUA:O5'",0.205181849438));
    this->alphas.insert(std::pair<std::string,double>("H8:GUA:O3'",3.56158512944));
    this->alphas.insert(std::pair<std::string,double>("H8:GUA:C2",-2.69736866523));
    this->alphas.insert(std::pair<std::string,double>("H8:GUA:C4",-1.278034993));
    this->alphas.insert(std::pair<std::string,double>("H8:GUA:C5",0.0567206506642));
    this->alphas.insert(std::pair<std::string,double>("H8:GUA:C6",-2.47397542994));
    this->alphas.insert(std::pair<std::string,double>("H8:GUA:C8",0.0745025487613));
    this->alphas.insert(std::pair<std::string,double>("H8:GUA:N1",-10.1495445705));
    this->alphas.insert(std::pair<std::string,double>("H8:GUA:N2",10.4355893394));
    this->alphas.insert(std::pair<std::string,double>("H8:GUA:N3",-3.09616040374));
    this->alphas.insert(std::pair<std::string,double>("H8:GUA:N7",2.16610153393));
    this->alphas.insert(std::pair<std::string,double>("H8:GUA:N9",0.832700368469));
    this->alphas.insert(std::pair<std::string,double>("H8:GUA:O6",-0.466839945921));
    this->alphas.insert(std::pair<std::string,double>("H8:ADE:C1'",3.52704450421));
    this->alphas.insert(std::pair<std::string,double>("H8:ADE:C2'",-3.75374008142));
    this->alphas.insert(std::pair<std::string,double>("H8:ADE:C3'",-11.0907022456));
    this->alphas.insert(std::pair<std::string,double>("H8:ADE:C4'",-0.478774580883));
    this->alphas.insert(std::pair<std::string,double>("H8:ADE:C5'",3.77839459615));
    this->alphas.insert(std::pair<std::string,double>("H8:ADE:P",7.67424558464));
    this->alphas.insert(std::pair<std::string,double>("H8:ADE:O5'",0.297596357271));
    this->alphas.insert(std::pair<std::string,double>("H8:ADE:O3'",-4.57043202239));
    this->alphas.insert(std::pair<std::string,double>("H8:ADE:C2",2.2290120441));
    this->alphas.insert(std::pair<std::string,double>("H8:ADE:C4",-5.44125700448));
    this->alphas.insert(std::pair<std::string,double>("H8:ADE:C5",-0.899917093831));
    this->alphas.insert(std::pair<std::string,double>("H8:ADE:C6",-0.181914952374));
    this->alphas.insert(std::pair<std::string,double>("H8:ADE:C8",-0.270472242826));
    this->alphas.insert(std::pair<std::string,double>("H8:ADE:N1",4.42622379352));
    this->alphas.insert(std::pair<std::string,double>("H8:ADE:N3",-10.8721776816));
    this->alphas.insert(std::pair<std::string,double>("H8:ADE:N6",1.99181605224));
    this->alphas.insert(std::pair<std::string,double>("H8:ADE:N7",0.995630718092));
    this->alphas.insert(std::pair<std::string,double>("H8:ADE:N9",1.03744918127));
    this->alphas.insert(std::pair<std::string,double>("H8:URA:C1'",8.7558424199));
    this->alphas.insert(std::pair<std::string,double>("H8:URA:C2'",-7.14921251737));
    this->alphas.insert(std::pair<std::string,double>("H8:URA:C3'",-4.37090163167));
    this->alphas.insert(std::pair<std::string,double>("H8:URA:C4'",22.0578397163));
    this->alphas.insert(std::pair<std::string,double>("H8:URA:C5'",-33.6272663834));
    this->alphas.insert(std::pair<std::string,double>("H8:URA:P",5.16592469422));
    this->alphas.insert(std::pair<std::string,double>("H8:URA:O5'",17.1185798436));
    this->alphas.insert(std::pair<std::string,double>("H8:URA:O3'",-3.17826835622));
    this->alphas.insert(std::pair<std::string,double>("H8:URA:C2",2.99060807107));
    this->alphas.insert(std::pair<std::string,double>("H8:URA:C4",1.29933542607));
    this->alphas.insert(std::pair<std::string,double>("H8:URA:C5",8.99822223527));
    this->alphas.insert(std::pair<std::string,double>("H8:URA:C6",0.188915420359));
    this->alphas.insert(std::pair<std::string,double>("H8:URA:N1",2.15678635216));
    this->alphas.insert(std::pair<std::string,double>("H8:URA:N3",-15.9963722626));
    this->alphas.insert(std::pair<std::string,double>("H8:URA:O4",-9.54813957384));
    this->alphas.insert(std::pair<std::string,double>("H8:CYT:C1'",9.3648822583));
    this->alphas.insert(std::pair<std::string,double>("H8:CYT:C2'",-2.83366318774));
    this->alphas.insert(std::pair<std::string,double>("H8:CYT:C3'",-3.40834060478));
    this->alphas.insert(std::pair<std::string,double>("H8:CYT:C4'",-3.92155591952));
    this->alphas.insert(std::pair<std::string,double>("H8:CYT:C5'",25.6156743029));
    this->alphas.insert(std::pair<std::string,double>("H8:CYT:P",15.257012093));
    this->alphas.insert(std::pair<std::string,double>("H8:CYT:O5'",-11.0162561204));
    this->alphas.insert(std::pair<std::string,double>("H8:CYT:O3'",-8.54731883901));
    this->alphas.insert(std::pair<std::string,double>("H8:CYT:C2",-3.41419878573));
    this->alphas.insert(std::pair<std::string,double>("H8:CYT:C4",-4.4524182262));
    this->alphas.insert(std::pair<std::string,double>("H8:CYT:C5",4.01507634801));
    this->alphas.insert(std::pair<std::string,double>("H8:CYT:C6",2.39771595124));
    this->alphas.insert(std::pair<std::string,double>("H8:CYT:N1",0.00828300837738));
    this->alphas.insert(std::pair<std::string,double>("H8:CYT:N3",-12.8018546124));
    this->alphas.insert(std::pair<std::string,double>("H8:CYT:N4",0.108888464294));
    this->alphas.insert(std::pair<std::string,double>("H8:CYT:O2",-0.469277322514));
    this->alphas.insert(std::pair<std::string,double>("H5':GUA:C1'",-8.44071572782));
    this->alphas.insert(std::pair<std::string,double>("H5':GUA:C2'",1.8752436953));
    this->alphas.insert(std::pair<std::string,double>("H5':GUA:C3'",-4.48117610879));
    this->alphas.insert(std::pair<std::string,double>("H5':GUA:C4'",0.234220002026));
    this->alphas.insert(std::pair<std::string,double>("H5':GUA:C5'",0.857103839465));
    this->alphas.insert(std::pair<std::string,double>("H5':GUA:P",4.11313783673));
    this->alphas.insert(std::pair<std::string,double>("H5':GUA:O5'",-0.231174055532));
    this->alphas.insert(std::pair<std::string,double>("H5':GUA:O3'",1.18745571032));
    this->alphas.insert(std::pair<std::string,double>("H5':GUA:C2",-2.76054023465));
    this->alphas.insert(std::pair<std::string,double>("H5':GUA:C4",-14.9411124745));
    this->alphas.insert(std::pair<std::string,double>("H5':GUA:C5",-29.3784338021));
    this->alphas.insert(std::pair<std::string,double>("H5':GUA:C6",-11.8425565839));
    this->alphas.insert(std::pair<std::string,double>("H5':GUA:C8",8.81887564713));
    this->alphas.insert(std::pair<std::string,double>("H5':GUA:N1",-1.63161700471));
    this->alphas.insert(std::pair<std::string,double>("H5':GUA:N2",-6.06136003562));
    this->alphas.insert(std::pair<std::string,double>("H5':GUA:N3",30.853675298));
    this->alphas.insert(std::pair<std::string,double>("H5':GUA:N7",5.3530752714));
    this->alphas.insert(std::pair<std::string,double>("H5':GUA:N9",-3.97021274524));
    this->alphas.insert(std::pair<std::string,double>("H5':GUA:O6",36.7719994209));
    this->alphas.insert(std::pair<std::string,double>("H5':ADE:C1'",-13.9938375685));
    this->alphas.insert(std::pair<std::string,double>("H5':ADE:C2'",-8.96147392926));
    this->alphas.insert(std::pair<std::string,double>("H5':ADE:C3'",0.502416888463));
    this->alphas.insert(std::pair<std::string,double>("H5':ADE:C4'",4.86108013938));
    this->alphas.insert(std::pair<std::string,double>("H5':ADE:C5'",0.102754464532));
    this->alphas.insert(std::pair<std::string,double>("H5':ADE:P",0.317818936758));
    this->alphas.insert(std::pair<std::string,double>("H5':ADE:O5'",1.13950631352));
    this->alphas.insert(std::pair<std::string,double>("H5':ADE:O3'",3.21320608882));
    this->alphas.insert(std::pair<std::string,double>("H5':ADE:C2",8.95627600854));
    this->alphas.insert(std::pair<std::string,double>("H5':ADE:C4",-18.7134373485));
    this->alphas.insert(std::pair<std::string,double>("H5':ADE:C5",-21.2412508963));
    this->alphas.insert(std::pair<std::string,double>("H5':ADE:C6",-15.5030849437));
    this->alphas.insert(std::pair<std::string,double>("H5':ADE:C8",46.1911451803));
    this->alphas.insert(std::pair<std::string,double>("H5':ADE:N1",29.621936072));
    this->alphas.insert(std::pair<std::string,double>("H5':ADE:N3",-6.73405453253));
    this->alphas.insert(std::pair<std::string,double>("H5':ADE:N6",-2.25602938268));
    this->alphas.insert(std::pair<std::string,double>("H5':ADE:N7",-41.9721845811));
    this->alphas.insert(std::pair<std::string,double>("H5':ADE:N9",12.6886019538));
    this->alphas.insert(std::pair<std::string,double>("H5':URA:C1'",-2.16842104664));
    this->alphas.insert(std::pair<std::string,double>("H5':URA:C2'",13.6293292753));
    this->alphas.insert(std::pair<std::string,double>("H5':URA:C3'",-6.95405795781));
    this->alphas.insert(std::pair<std::string,double>("H5':URA:C4'",-0.996712744262));
    this->alphas.insert(std::pair<std::string,double>("H5':URA:C5'",0.742522506903));
    this->alphas.insert(std::pair<std::string,double>("H5':URA:P",5.70688301664));
    this->alphas.insert(std::pair<std::string,double>("H5':URA:O5'",0.422165968021));
    this->alphas.insert(std::pair<std::string,double>("H5':URA:O3'",-6.10740160504));
    this->alphas.insert(std::pair<std::string,double>("H5':URA:C2",-21.3944056028));
    this->alphas.insert(std::pair<std::string,double>("H5':URA:C4",18.4725724635));
    this->alphas.insert(std::pair<std::string,double>("H5':URA:C5",-15.8585287763));
    this->alphas.insert(std::pair<std::string,double>("H5':URA:C6",1.77552510105));
    this->alphas.insert(std::pair<std::string,double>("H5':URA:N1",-7.26652083368));
    this->alphas.insert(std::pair<std::string,double>("H5':URA:N3",22.528584822));
    this->alphas.insert(std::pair<std::string,double>("H5':URA:O4",35.647471832));
    this->alphas.insert(std::pair<std::string,double>("H5':CYT:C1'",-11.00003478));
    this->alphas.insert(std::pair<std::string,double>("H5':CYT:C2'",-2.50643664735));
    this->alphas.insert(std::pair<std::string,double>("H5':CYT:C3'",1.30151430032));
    this->alphas.insert(std::pair<std::string,double>("H5':CYT:C4'",1.8052910683));
    this->alphas.insert(std::pair<std::string,double>("H5':CYT:C5'",0.138541957319));
    this->alphas.insert(std::pair<std::string,double>("H5':CYT:P",7.20877204834));
    this->alphas.insert(std::pair<std::string,double>("H5':CYT:O5'",1.74523712771));
    this->alphas.insert(std::pair<std::string,double>("H5':CYT:O3'",2.35480878038));
    this->alphas.insert(std::pair<std::string,double>("H5':CYT:C2",-11.1129874154));
    this->alphas.insert(std::pair<std::string,double>("H5':CYT:C4",-31.9628760684));
    this->alphas.insert(std::pair<std::string,double>("H5':CYT:C5",37.4266675818));
    this->alphas.insert(std::pair<std::string,double>("H5':CYT:C6",-1.99245750623));
    this->alphas.insert(std::pair<std::string,double>("H5':CYT:N1",-12.4403658348));
    this->alphas.insert(std::pair<std::string,double>("H5':CYT:N3",-6.12722995569));
    this->alphas.insert(std::pair<std::string,double>("H5':CYT:N4",-5.33821382262));
    this->alphas.insert(std::pair<std::string,double>("H5':CYT:O2",15.278789451));
    this->alphas.insert(std::pair<std::string,double>("H3':GUA:C1'",-1.64377881421));
    this->alphas.insert(std::pair<std::string,double>("H3':GUA:C2'",1.44163564589));
    this->alphas.insert(std::pair<std::string,double>("H3':GUA:C3'",-0.0));
    this->alphas.insert(std::pair<std::string,double>("H3':GUA:C4'",0.0188538222202));
    this->alphas.insert(std::pair<std::string,double>("H3':GUA:C5'",1.95436198594));
    this->alphas.insert(std::pair<std::string,double>("H3':GUA:P",3.25098142074));
    this->alphas.insert(std::pair<std::string,double>("H3':GUA:O5'",-0.286634061208));
    this->alphas.insert(std::pair<std::string,double>("H3':GUA:O3'",1.27785305986));
    this->alphas.insert(std::pair<std::string,double>("H3':GUA:C2",-7.84153735907));
    this->alphas.insert(std::pair<std::string,double>("H3':GUA:C4",-0.251897570061));
    this->alphas.insert(std::pair<std::string,double>("H3':GUA:C5",-8.65116941654));
    this->alphas.insert(std::pair<std::string,double>("H3':GUA:C6",-17.0470697079));
    this->alphas.insert(std::pair<std::string,double>("H3':GUA:C8",0.172357243887));
    this->alphas.insert(std::pair<std::string,double>("H3':GUA:N1",-8.58599641342));
    this->alphas.insert(std::pair<std::string,double>("H3':GUA:N2",6.36963166123));
    this->alphas.insert(std::pair<std::string,double>("H3':GUA:N3",17.9855649403));
    this->alphas.insert(std::pair<std::string,double>("H3':GUA:N7",4.54704703156));
    this->alphas.insert(std::pair<std::string,double>("H3':GUA:N9",-2.67152176475));
    this->alphas.insert(std::pair<std::string,double>("H3':GUA:O6",7.73167920399));
    this->alphas.insert(std::pair<std::string,double>("H3':ADE:C1'",0.175979103481));
    this->alphas.insert(std::pair<std::string,double>("H3':ADE:C2'",0.0149092600303));
    this->alphas.insert(std::pair<std::string,double>("H3':ADE:C3'",0.727303502292));
    this->alphas.insert(std::pair<std::string,double>("H3':ADE:C4'",-0.0283881674575));
    this->alphas.insert(std::pair<std::string,double>("H3':ADE:C5'",1.916417227));
    this->alphas.insert(std::pair<std::string,double>("H3':ADE:P",3.89152201245));
    this->alphas.insert(std::pair<std::string,double>("H3':ADE:O5'",-1.18825286383));
    this->alphas.insert(std::pair<std::string,double>("H3':ADE:O3'",-1.57936473722));
    this->alphas.insert(std::pair<std::string,double>("H3':ADE:C2",2.6542847575));
    this->alphas.insert(std::pair<std::string,double>("H3':ADE:C4",-2.47786178306));
    this->alphas.insert(std::pair<std::string,double>("H3':ADE:C5",-7.56136599687));
    this->alphas.insert(std::pair<std::string,double>("H3':ADE:C6",-10.7998748174));
    this->alphas.insert(std::pair<std::string,double>("H3':ADE:C8",3.20450872743));
    this->alphas.insert(std::pair<std::string,double>("H3':ADE:N1",16.111075801));
    this->alphas.insert(std::pair<std::string,double>("H3':ADE:N3",-14.6180472161));
    this->alphas.insert(std::pair<std::string,double>("H3':ADE:N6",2.03932479979));
    this->alphas.insert(std::pair<std::string,double>("H3':ADE:N7",-0.474990961129));
    this->alphas.insert(std::pair<std::string,double>("H3':ADE:N9",-0.2072750359));
    this->alphas.insert(std::pair<std::string,double>("H3':URA:C1'",-0.00081685001849));
    this->alphas.insert(std::pair<std::string,double>("H3':URA:C2'",1.23848728448));
    this->alphas.insert(std::pair<std::string,double>("H3':URA:C3'",0.340917141751));
    this->alphas.insert(std::pair<std::string,double>("H3':URA:C4'",-1.41405156342));
    this->alphas.insert(std::pair<std::string,double>("H3':URA:C5'",-1.58027394554));
    this->alphas.insert(std::pair<std::string,double>("H3':URA:P",1.31124533661));
    this->alphas.insert(std::pair<std::string,double>("H3':URA:O5'",0.625543527478));
    this->alphas.insert(std::pair<std::string,double>("H3':URA:O3'",-0.221857256205));
    this->alphas.insert(std::pair<std::string,double>("H3':URA:C2",4.15515479807));
    this->alphas.insert(std::pair<std::string,double>("H3':URA:C4",0.473065849415));
    this->alphas.insert(std::pair<std::string,double>("H3':URA:C5",0.904693261332));
    this->alphas.insert(std::pair<std::string,double>("H3':URA:C6",-0.291646081039));
    this->alphas.insert(std::pair<std::string,double>("H3':URA:N1",0.380766469758));
    this->alphas.insert(std::pair<std::string,double>("H3':URA:N3",-4.25137541707));
    this->alphas.insert(std::pair<std::string,double>("H3':URA:O4",2.16879657134));
    this->alphas.insert(std::pair<std::string,double>("H3':CYT:C1'",1.26028054892));
    this->alphas.insert(std::pair<std::string,double>("H3':CYT:C2'",0.279194094177));
    this->alphas.insert(std::pair<std::string,double>("H3':CYT:C3'",-0.00574197745572));
    this->alphas.insert(std::pair<std::string,double>("H3':CYT:C4'",0.14377589193));
    this->alphas.insert(std::pair<std::string,double>("H3':CYT:C5'",1.04994706915));
    this->alphas.insert(std::pair<std::string,double>("H3':CYT:P",2.01267368136));
    this->alphas.insert(std::pair<std::string,double>("H3':CYT:O5'",0.232390738304));
    this->alphas.insert(std::pair<std::string,double>("H3':CYT:O3'",-0.0040041348529));
    this->alphas.insert(std::pair<std::string,double>("H3':CYT:C2",1.36946286595));
    this->alphas.insert(std::pair<std::string,double>("H3':CYT:C4",1.7800018514));
    this->alphas.insert(std::pair<std::string,double>("H3':CYT:C5",-1.73139214366));
    this->alphas.insert(std::pair<std::string,double>("H3':CYT:C6",0.0996681263081));
    this->alphas.insert(std::pair<std::string,double>("H3':CYT:N1",-0.232624902118));
    this->alphas.insert(std::pair<std::string,double>("H3':CYT:N3",-4.02995634));
    this->alphas.insert(std::pair<std::string,double>("H3':CYT:N4",7.17115037889));
    this->alphas.insert(std::pair<std::string,double>("H3':CYT:O2",-0.686987553812));
    this->alphas.insert(std::pair<std::string,double>("H2:GUA:C1'",-5.46113414466));
    this->alphas.insert(std::pair<std::string,double>("H2:GUA:C2'",1.17993849123));
    this->alphas.insert(std::pair<std::string,double>("H2:GUA:C3'",-2.00802642559));
    this->alphas.insert(std::pair<std::string,double>("H2:GUA:C4'",3.04491844088));
    this->alphas.insert(std::pair<std::string,double>("H2:GUA:C5'",-9.19886010936));
    this->alphas.insert(std::pair<std::string,double>("H2:GUA:P",11.5341223322));
    this->alphas.insert(std::pair<std::string,double>("H2:GUA:O5'",-5.36025008063));
    this->alphas.insert(std::pair<std::string,double>("H2:GUA:O3'",70.8349884692));
    this->alphas.insert(std::pair<std::string,double>("H2:GUA:C2",-7.36349565976));
    this->alphas.insert(std::pair<std::string,double>("H2:GUA:C4",-8.0327985787));
    this->alphas.insert(std::pair<std::string,double>("H2:GUA:C5",1.11125948335));
    this->alphas.insert(std::pair<std::string,double>("H2:GUA:C6",3.29992572077));
    this->alphas.insert(std::pair<std::string,double>("H2:GUA:C8",-3.87839267324));
    this->alphas.insert(std::pair<std::string,double>("H2:GUA:N1",-0.547023548503));
    this->alphas.insert(std::pair<std::string,double>("H2:GUA:N2",8.35498770847));
    this->alphas.insert(std::pair<std::string,double>("H2:GUA:N3",-4.09246766531));
    this->alphas.insert(std::pair<std::string,double>("H2:GUA:N7",-10.9096339759));
    this->alphas.insert(std::pair<std::string,double>("H2:GUA:N9",4.02422373348));
    this->alphas.insert(std::pair<std::string,double>("H2:GUA:O6",2.82493428675));
    this->alphas.insert(std::pair<std::string,double>("H2:ADE:C1'",-5.43934389942));
    this->alphas.insert(std::pair<std::string,double>("H2:ADE:C2'",7.03643407835));
    this->alphas.insert(std::pair<std::string,double>("H2:ADE:C3'",13.6097531761));
    this->alphas.insert(std::pair<std::string,double>("H2:ADE:C4'",-2.24993943313));
    this->alphas.insert(std::pair<std::string,double>("H2:ADE:C5'",13.8789137533));
    this->alphas.insert(std::pair<std::string,double>("H2:ADE:P",7.71139170047));
    this->alphas.insert(std::pair<std::string,double>("H2:ADE:O5'",-31.7626125005));
    this->alphas.insert(std::pair<std::string,double>("H2:ADE:O3'",39.5023918732));
    this->alphas.insert(std::pair<std::string,double>("H2:ADE:C2",-1.96191523492));
    this->alphas.insert(std::pair<std::string,double>("H2:ADE:C4",-4.86346526461));
    this->alphas.insert(std::pair<std::string,double>("H2:ADE:C5",-3.42166269582));
    this->alphas.insert(std::pair<std::string,double>("H2:ADE:C6",10.286586097));
    this->alphas.insert(std::pair<std::string,double>("H2:ADE:C8",-14.6259086489));
    this->alphas.insert(std::pair<std::string,double>("H2:ADE:N1",5.87569946522));
    this->alphas.insert(std::pair<std::string,double>("H2:ADE:N3",1.84448378874));
    this->alphas.insert(std::pair<std::string,double>("H2:ADE:N6",2.11839966655));
    this->alphas.insert(std::pair<std::string,double>("H2:ADE:N7",-3.81308269133));
    this->alphas.insert(std::pair<std::string,double>("H2:ADE:N9",-19.2435336057));
    this->alphas.insert(std::pair<std::string,double>("H2:URA:C1'",4.21339339235));
    this->alphas.insert(std::pair<std::string,double>("H2:URA:C2'",-2.66299578913));
    this->alphas.insert(std::pair<std::string,double>("H2:URA:C3'",7.05506416845));
    this->alphas.insert(std::pair<std::string,double>("H2:URA:C4'",18.6329978406));
    this->alphas.insert(std::pair<std::string,double>("H2:URA:C5'",-4.34814159235));
    this->alphas.insert(std::pair<std::string,double>("H2:URA:P",10.8486757079));
    this->alphas.insert(std::pair<std::string,double>("H2:URA:O5'",29.814938725));
    this->alphas.insert(std::pair<std::string,double>("H2:URA:O3'",-33.2445387743));
    this->alphas.insert(std::pair<std::string,double>("H2:URA:C2",-1.48697436872));
    this->alphas.insert(std::pair<std::string,double>("H2:URA:C4",-0.11124084618));
    this->alphas.insert(std::pair<std::string,double>("H2:URA:C5",-9.63916566524));
    this->alphas.insert(std::pair<std::string,double>("H2:URA:C6",-18.4271913979));
    this->alphas.insert(std::pair<std::string,double>("H2:URA:N1",6.08515775436));
    this->alphas.insert(std::pair<std::string,double>("H2:URA:N3",-4.40163709982));
    this->alphas.insert(std::pair<std::string,double>("H2:URA:O4",-0.42523909654));
    this->alphas.insert(std::pair<std::string,double>("H2:CYT:C1'",8.53170872151));
    this->alphas.insert(std::pair<std::string,double>("H2:CYT:C2'",-26.9155883409));
    this->alphas.insert(std::pair<std::string,double>("H2:CYT:C3'",-9.15973565482));
    this->alphas.insert(std::pair<std::string,double>("H2:CYT:C4'",52.7155007838));
    this->alphas.insert(std::pair<std::string,double>("H2:CYT:C5'",29.896488049));
    this->alphas.insert(std::pair<std::string,double>("H2:CYT:P",-34.2560535023));
    this->alphas.insert(std::pair<std::string,double>("H2:CYT:O5'",-20.5309613086));
    this->alphas.insert(std::pair<std::string,double>("H2:CYT:O3'",-11.5449767971));
    this->alphas.insert(std::pair<std::string,double>("H2:CYT:C2",0.490358986334));
    this->alphas.insert(std::pair<std::string,double>("H2:CYT:C4",5.92594650161));
    this->alphas.insert(std::pair<std::string,double>("H2:CYT:C5",-4.8600581481));
    this->alphas.insert(std::pair<std::string,double>("H2:CYT:C6",-13.9469618462));
    this->alphas.insert(std::pair<std::string,double>("H2:CYT:N1",0.651819666906));
    this->alphas.insert(std::pair<std::string,double>("H2:CYT:N3",2.43450454067));
    this->alphas.insert(std::pair<std::string,double>("H2:CYT:N4",-0.330564985395));
    this->alphas.insert(std::pair<std::string,double>("H2:CYT:O2",-0.986902148517));
    this->alphas.insert(std::pair<std::string,double>("H3:GUA:N7",-9.90557456176));
    this->alphas.insert(std::pair<std::string,double>("H3:GUA:C2",13.2184695744));
    this->alphas.insert(std::pair<std::string,double>("H3:GUA:C4'",26.8577693898));
    this->alphas.insert(std::pair<std::string,double>("H3:GUA:N9",-14.9610360713));
    this->alphas.insert(std::pair<std::string,double>("H3:GUA:O5'",-52.5675912336));
    this->alphas.insert(std::pair<std::string,double>("H3:GUA:C3'",64.9374646361));
    this->alphas.insert(std::pair<std::string,double>("H3:GUA:C1'",-34.8620045721));
    this->alphas.insert(std::pair<std::string,double>("H3:GUA:C5",1.28356132507));
    this->alphas.insert(std::pair<std::string,double>("H3:GUA:C5'",52.9434434976));
    this->alphas.insert(std::pair<std::string,double>("H3:GUA:P",-3.51554075412));
    this->alphas.insert(std::pair<std::string,double>("H3:GUA:C8",-1.9581656465));
    this->alphas.insert(std::pair<std::string,double>("H3:GUA:C2'",40.156485581));
    this->alphas.insert(std::pair<std::string,double>("H3:GUA:O6",-6.11961833588));
    this->alphas.insert(std::pair<std::string,double>("H3:GUA:N1",3.7015664332));
    this->alphas.insert(std::pair<std::string,double>("H3:GUA:N2",5.71348570537));
    this->alphas.insert(std::pair<std::string,double>("H3:GUA:N3",5.10371087346));
    this->alphas.insert(std::pair<std::string,double>("H3:GUA:C6",5.49416381763));
    this->alphas.insert(std::pair<std::string,double>("H3:GUA:O3'",83.9347577362));
    this->alphas.insert(std::pair<std::string,double>("H3:GUA:C4",-9.28211830372));
    this->alphas.insert(std::pair<std::string,double>("H3:ADE:N7",12.5607242793));
    this->alphas.insert(std::pair<std::string,double>("H3:ADE:N6",9.55860596126));
    this->alphas.insert(std::pair<std::string,double>("H3:ADE:N9",-29.3796756729));
    this->alphas.insert(std::pair<std::string,double>("H3:ADE:C1'",3.17083792606));
    this->alphas.insert(std::pair<std::string,double>("H3:ADE:C3'",18.0090182513));
    this->alphas.insert(std::pair<std::string,double>("H3:ADE:O5'",-61.7384941189));
    this->alphas.insert(std::pair<std::string,double>("H3:ADE:C5",-1.51869271528));
    this->alphas.insert(std::pair<std::string,double>("H3:ADE:C5'",-31.3665002257));
    this->alphas.insert(std::pair<std::string,double>("H3:ADE:P",-61.9664635691));
    this->alphas.insert(std::pair<std::string,double>("H3:ADE:C4'",-8.46929838487));
    this->alphas.insert(std::pair<std::string,double>("H3:ADE:C2'",91.2258366886));
    this->alphas.insert(std::pair<std::string,double>("H3:ADE:N1",8.91518393344));
    this->alphas.insert(std::pair<std::string,double>("H3:ADE:C2",11.169832639));
    this->alphas.insert(std::pair<std::string,double>("H3:ADE:C8",-15.2230780589));
    this->alphas.insert(std::pair<std::string,double>("H3:ADE:N3",-5.67780526186));
    this->alphas.insert(std::pair<std::string,double>("H3:ADE:C6",-2.64625568141));
    this->alphas.insert(std::pair<std::string,double>("H3:ADE:O3'",-4.35718452693));
    this->alphas.insert(std::pair<std::string,double>("H3:ADE:C4",-25.5239290312));
    this->alphas.insert(std::pair<std::string,double>("H3:URA:O4",-1.6693842823));
    this->alphas.insert(std::pair<std::string,double>("H3:URA:O5'",-6.61503969182));
    this->alphas.insert(std::pair<std::string,double>("H3:URA:C3'",-11.077494508));
    this->alphas.insert(std::pair<std::string,double>("H3:URA:O3'",-71.3546165555));
    this->alphas.insert(std::pair<std::string,double>("H3:URA:C1'",27.3576794288));
    this->alphas.insert(std::pair<std::string,double>("H3:URA:C5'",-22.8684351677));
    this->alphas.insert(std::pair<std::string,double>("H3:URA:P",-25.6165547676));
    this->alphas.insert(std::pair<std::string,double>("H3:URA:C4'",9.96800153072));
    this->alphas.insert(std::pair<std::string,double>("H3:URA:C2'",5.24120407275));
    this->alphas.insert(std::pair<std::string,double>("H3:URA:N1",-0.312154189752));
    this->alphas.insert(std::pair<std::string,double>("H3:URA:C2",-2.62167860044));
    this->alphas.insert(std::pair<std::string,double>("H3:URA:N3",-0.1652688998));
    this->alphas.insert(std::pair<std::string,double>("H3:URA:C6",4.96159495346));
    this->alphas.insert(std::pair<std::string,double>("H3:URA:C5",1.3644318216));
    this->alphas.insert(std::pair<std::string,double>("H3:URA:C4",-5.14269223728));
    this->alphas.insert(std::pair<std::string,double>("H3:CYT:O5'",-6.79123882521));
    this->alphas.insert(std::pair<std::string,double>("H3:CYT:C3'",-47.4025922922));
    this->alphas.insert(std::pair<std::string,double>("H3:CYT:O3'",-114.895562462));
    this->alphas.insert(std::pair<std::string,double>("H3:CYT:C1'",50.3553017232));
    this->alphas.insert(std::pair<std::string,double>("H3:CYT:O2",0.0676760120205));
    this->alphas.insert(std::pair<std::string,double>("H3:CYT:C5'",-52.4860640466));
    this->alphas.insert(std::pair<std::string,double>("H3:CYT:P",-4.73356132565));
    this->alphas.insert(std::pair<std::string,double>("H3:CYT:C4'",-18.6572747597));
    this->alphas.insert(std::pair<std::string,double>("H3:CYT:C2'",-10.5829975514));
    this->alphas.insert(std::pair<std::string,double>("H3:CYT:N1",12.1223539563));
    this->alphas.insert(std::pair<std::string,double>("H3:CYT:C2",-2.13126511077));
    this->alphas.insert(std::pair<std::string,double>("H3:CYT:N3",0.17403809902));
    this->alphas.insert(std::pair<std::string,double>("H3:CYT:N4",-4.97461156671));
    this->alphas.insert(std::pair<std::string,double>("H3:CYT:C6",5.36643494926));
    this->alphas.insert(std::pair<std::string,double>("H3:CYT:C5",-0.671220957929));
    this->alphas.insert(std::pair<std::string,double>("H3:CYT:C4",2.70947652494));
    this->alphas.insert(std::pair<std::string,double>("H1:GUA:N7",-5.60107259381));
    this->alphas.insert(std::pair<std::string,double>("H1:GUA:C2",-1.87842557685));
    this->alphas.insert(std::pair<std::string,double>("H1:GUA:C4'",-6.65853692202));
    this->alphas.insert(std::pair<std::string,double>("H1:GUA:N9",4.0420177122));
    this->alphas.insert(std::pair<std::string,double>("H1:GUA:O5'",-24.9804818878));
    this->alphas.insert(std::pair<std::string,double>("H1:GUA:C3'",15.7349194051));
    this->alphas.insert(std::pair<std::string,double>("H1:GUA:C1'",19.6585044563));
    this->alphas.insert(std::pair<std::string,double>("H1:GUA:C5",-6.23791509446));
    this->alphas.insert(std::pair<std::string,double>("H1:GUA:C5'",7.28801756652));
    this->alphas.insert(std::pair<std::string,double>("H1:GUA:P",11.3643133938));
    this->alphas.insert(std::pair<std::string,double>("H1:GUA:C8",1.42853498487));
    this->alphas.insert(std::pair<std::string,double>("H1:GUA:C2'",20.600074305));
    this->alphas.insert(std::pair<std::string,double>("H1:GUA:O6",-4.14044333419));
    this->alphas.insert(std::pair<std::string,double>("H1:GUA:N1",-0.00330494503579));
    this->alphas.insert(std::pair<std::string,double>("H1:GUA:N2",-1.17265247106));
    this->alphas.insert(std::pair<std::string,double>("H1:GUA:N3",-16.2529366656));
    this->alphas.insert(std::pair<std::string,double>("H1:GUA:C6",-0.419380221405));
    this->alphas.insert(std::pair<std::string,double>("H1:GUA:O3'",59.5008667949));
    this->alphas.insert(std::pair<std::string,double>("H1:GUA:C4",-8.62230013167));
    this->alphas.insert(std::pair<std::string,double>("H1:ADE:N7",7.57502382674));
    this->alphas.insert(std::pair<std::string,double>("H1:ADE:N6",6.54784290165));
    this->alphas.insert(std::pair<std::string,double>("H1:ADE:N9",-4.0533908781));
    this->alphas.insert(std::pair<std::string,double>("H1:ADE:C1'",-22.8998376328));
    this->alphas.insert(std::pair<std::string,double>("H1:ADE:C3'",4.44027812873));
    this->alphas.insert(std::pair<std::string,double>("H1:ADE:O5'",22.1099680175));
    this->alphas.insert(std::pair<std::string,double>("H1:ADE:C5",-12.409874331));
    this->alphas.insert(std::pair<std::string,double>("H1:ADE:C5'",-1.36370540826));
    this->alphas.insert(std::pair<std::string,double>("H1:ADE:P",-79.8642046724));
    this->alphas.insert(std::pair<std::string,double>("H1:ADE:C4'",13.3808815755));
    this->alphas.insert(std::pair<std::string,double>("H1:ADE:C2'",4.56193495301));
    this->alphas.insert(std::pair<std::string,double>("H1:ADE:N1",-5.00480035898));
    this->alphas.insert(std::pair<std::string,double>("H1:ADE:C2",6.4831698135));
    this->alphas.insert(std::pair<std::string,double>("H1:ADE:C8",12.5268593142));
    this->alphas.insert(std::pair<std::string,double>("H1:ADE:N3",-2.70930018182));
    this->alphas.insert(std::pair<std::string,double>("H1:ADE:C6",-16.5674250417));
    this->alphas.insert(std::pair<std::string,double>("H1:ADE:O3'",36.3377913618));
    this->alphas.insert(std::pair<std::string,double>("H1:ADE:C4",-4.85371230369));
    this->alphas.insert(std::pair<std::string,double>("H1:URA:O4",5.68829352631));
    this->alphas.insert(std::pair<std::string,double>("H1:URA:O5'",5.84749620117));
    this->alphas.insert(std::pair<std::string,double>("H1:URA:C3'",-64.8064305308));
    this->alphas.insert(std::pair<std::string,double>("H1:URA:O3'",-80.4701932907));
    this->alphas.insert(std::pair<std::string,double>("H1:URA:C1'",26.7790568615));
    this->alphas.insert(std::pair<std::string,double>("H1:URA:C5'",42.0227783969));
    this->alphas.insert(std::pair<std::string,double>("H1:URA:P",14.2425936528));
    this->alphas.insert(std::pair<std::string,double>("H1:URA:C4'",62.3810706766));
    this->alphas.insert(std::pair<std::string,double>("H1:URA:C2'",-50.0850561192));
    this->alphas.insert(std::pair<std::string,double>("H1:URA:N1",22.6608263052));
    this->alphas.insert(std::pair<std::string,double>("H1:URA:C2",-4.83958991585));
    this->alphas.insert(std::pair<std::string,double>("H1:URA:N3",6.76234239923));
    this->alphas.insert(std::pair<std::string,double>("H1:URA:C6",1.53670457669));
    this->alphas.insert(std::pair<std::string,double>("H1:URA:C5",-22.0630596079));
    this->alphas.insert(std::pair<std::string,double>("H1:URA:C4",4.27526502637));
    this->alphas.insert(std::pair<std::string,double>("H1:CYT:O5'",-27.2662440308));
    this->alphas.insert(std::pair<std::string,double>("H1:CYT:C3'",-11.6178011132));
    this->alphas.insert(std::pair<std::string,double>("H1:CYT:O3'",-22.8744140273));
    this->alphas.insert(std::pair<std::string,double>("H1:CYT:C1'",-16.0469250337));
    this->alphas.insert(std::pair<std::string,double>("H1:CYT:O2",7.36914971981));
    this->alphas.insert(std::pair<std::string,double>("H1:CYT:C5'",7.9089870615));
    this->alphas.insert(std::pair<std::string,double>("H1:CYT:P",31.0524779388));
    this->alphas.insert(std::pair<std::string,double>("H1:CYT:C4'",0.817384702947));
    this->alphas.insert(std::pair<std::string,double>("H1:CYT:C2'",2.21206512961));
    this->alphas.insert(std::pair<std::string,double>("H1:CYT:N1",-1.19338518376));
    this->alphas.insert(std::pair<std::string,double>("H1:CYT:C2",3.47781892176));
    this->alphas.insert(std::pair<std::string,double>("H1:CYT:N3",8.26972546243));
    this->alphas.insert(std::pair<std::string,double>("H1:CYT:N4",5.52406558242));
    this->alphas.insert(std::pair<std::string,double>("H1:CYT:C6",-14.4542303866));
    this->alphas.insert(std::pair<std::string,double>("H1:CYT:C5",2.51226125817));
    this->alphas.insert(std::pair<std::string,double>("H1:CYT:C4",4.10136603248));
    this->alphas.insert(std::pair<std::string,double>("H6:GUA:C1'",6.01988641997));
    this->alphas.insert(std::pair<std::string,double>("H6:GUA:C2'",-5.08285644647));
    this->alphas.insert(std::pair<std::string,double>("H6:GUA:C3'",-4.52091866903));
    this->alphas.insert(std::pair<std::string,double>("H6:GUA:C4'",-0.433142999718));
    this->alphas.insert(std::pair<std::string,double>("H6:GUA:C5'",-4.56930402953));
    this->alphas.insert(std::pair<std::string,double>("H6:GUA:P",12.0796873257));
    this->alphas.insert(std::pair<std::string,double>("H6:GUA:O5'",-11.1597511502));
    this->alphas.insert(std::pair<std::string,double>("H6:GUA:O3'",2.16616730987));
    this->alphas.insert(std::pair<std::string,double>("H6:GUA:C2",-0.854393723056));
    this->alphas.insert(std::pair<std::string,double>("H6:GUA:C4",2.64061323232));
    this->alphas.insert(std::pair<std::string,double>("H6:GUA:C5",-1.57346543305));
    this->alphas.insert(std::pair<std::string,double>("H6:GUA:C6",-3.79979604374));
    this->alphas.insert(std::pair<std::string,double>("H6:GUA:C8",-1.47086756181));
    this->alphas.insert(std::pair<std::string,double>("H6:GUA:N1",-2.69602059671));
    this->alphas.insert(std::pair<std::string,double>("H6:GUA:N2",7.70719883091));
    this->alphas.insert(std::pair<std::string,double>("H6:GUA:N3",-1.27913733627));
    this->alphas.insert(std::pair<std::string,double>("H6:GUA:N7",-6.11746616629));
    this->alphas.insert(std::pair<std::string,double>("H6:GUA:N9",4.93050939031));
    this->alphas.insert(std::pair<std::string,double>("H6:GUA:O6",5.59238905239));
    this->alphas.insert(std::pair<std::string,double>("H6:ADE:C1'",-5.91643267378));
    this->alphas.insert(std::pair<std::string,double>("H6:ADE:C2'",-6.10255984373));
    this->alphas.insert(std::pair<std::string,double>("H6:ADE:C3'",0.799240407977));
    this->alphas.insert(std::pair<std::string,double>("H6:ADE:C4'",-1.77001130172));
    this->alphas.insert(std::pair<std::string,double>("H6:ADE:C5'",-5.94890184729));
    this->alphas.insert(std::pair<std::string,double>("H6:ADE:P",-17.0430204847));
    this->alphas.insert(std::pair<std::string,double>("H6:ADE:O5'",24.1732869202));
    this->alphas.insert(std::pair<std::string,double>("H6:ADE:O3'",-0.406367266508));
    this->alphas.insert(std::pair<std::string,double>("H6:ADE:C2",1.50825225401));
    this->alphas.insert(std::pair<std::string,double>("H6:ADE:C4",-6.13911568603));
    this->alphas.insert(std::pair<std::string,double>("H6:ADE:C5",-8.42599651672));
    this->alphas.insert(std::pair<std::string,double>("H6:ADE:C6",-2.83989718182));
    this->alphas.insert(std::pair<std::string,double>("H6:ADE:C8",8.10557608112));
    this->alphas.insert(std::pair<std::string,double>("H6:ADE:N1",20.2564641975));
    this->alphas.insert(std::pair<std::string,double>("H6:ADE:N3",-4.96756009057));
    this->alphas.insert(std::pair<std::string,double>("H6:ADE:N6",-3.1175293936));
    this->alphas.insert(std::pair<std::string,double>("H6:ADE:N7",-5.83934904549));
    this->alphas.insert(std::pair<std::string,double>("H6:ADE:N9",1.62032655081));
    this->alphas.insert(std::pair<std::string,double>("H6:URA:C1'",-0.0443361842129));
    this->alphas.insert(std::pair<std::string,double>("H6:URA:C2'",0.739379964858));
    this->alphas.insert(std::pair<std::string,double>("H6:URA:C3'",0.680086079535));
    this->alphas.insert(std::pair<std::string,double>("H6:URA:C4'",0.369043427392));
    this->alphas.insert(std::pair<std::string,double>("H6:URA:C5'",0.310549733109));
    this->alphas.insert(std::pair<std::string,double>("H6:URA:P",0.537440949349));
    this->alphas.insert(std::pair<std::string,double>("H6:URA:O5'",0.436050179417));
    this->alphas.insert(std::pair<std::string,double>("H6:URA:O3'",1.92352892091));
    this->alphas.insert(std::pair<std::string,double>("H6:URA:C2",-0.421147677745));
    this->alphas.insert(std::pair<std::string,double>("H6:URA:C4",-1.6198251992));
    this->alphas.insert(std::pair<std::string,double>("H6:URA:C5",-0.0249710559498));
    this->alphas.insert(std::pair<std::string,double>("H6:URA:C6",-0.0));
    this->alphas.insert(std::pair<std::string,double>("H6:URA:N1",-0.1007391836));
    this->alphas.insert(std::pair<std::string,double>("H6:URA:N3",-0.24909401143));
    this->alphas.insert(std::pair<std::string,double>("H6:URA:O4",-3.65887960987));
    this->alphas.insert(std::pair<std::string,double>("H6:CYT:C1'",-1.26768796735));
    this->alphas.insert(std::pair<std::string,double>("H6:CYT:C2'",-2.54564164504));
    this->alphas.insert(std::pair<std::string,double>("H6:CYT:C3'",-0.218257039925));
    this->alphas.insert(std::pair<std::string,double>("H6:CYT:C4'",0.29015387966));
    this->alphas.insert(std::pair<std::string,double>("H6:CYT:C5'",1.82577223504));
    this->alphas.insert(std::pair<std::string,double>("H6:CYT:P",1.58272498989));
    this->alphas.insert(std::pair<std::string,double>("H6:CYT:O5'",-0.288497257476));
    this->alphas.insert(std::pair<std::string,double>("H6:CYT:O3'",1.02862778963));
    this->alphas.insert(std::pair<std::string,double>("H6:CYT:C2",-0.555191466415));
    this->alphas.insert(std::pair<std::string,double>("H6:CYT:C4",-0.497986542593));
    this->alphas.insert(std::pair<std::string,double>("H6:CYT:C5",0.0669004703608));
    this->alphas.insert(std::pair<std::string,double>("H6:CYT:C6",-0.0));
    this->alphas.insert(std::pair<std::string,double>("H6:CYT:N1",-0.134744785147));
    this->alphas.insert(std::pair<std::string,double>("H6:CYT:N3",-0.38350657439));
    this->alphas.insert(std::pair<std::string,double>("H6:CYT:N4",4.29792738726));
    this->alphas.insert(std::pair<std::string,double>("H6:CYT:O2",0.387219300775));
    this->alphas.insert(std::pair<std::string,double>("H5:GUA:C1'",18.801493806));
    this->alphas.insert(std::pair<std::string,double>("H5:GUA:C2'",-10.5783213371));
    this->alphas.insert(std::pair<std::string,double>("H5:GUA:C3'",-33.3880729555));
    this->alphas.insert(std::pair<std::string,double>("H5:GUA:C4'",8.92435628157));
    this->alphas.insert(std::pair<std::string,double>("H5:GUA:C5'",10.8711364539));
    this->alphas.insert(std::pair<std::string,double>("H5:GUA:P",32.5455713216));
    this->alphas.insert(std::pair<std::string,double>("H5:GUA:O5'",-14.062766619));
    this->alphas.insert(std::pair<std::string,double>("H5:GUA:O3'",9.23250907634));
    this->alphas.insert(std::pair<std::string,double>("H5:GUA:C2",-3.33387600872));
    this->alphas.insert(std::pair<std::string,double>("H5:GUA:C4",0.804869301065));
    this->alphas.insert(std::pair<std::string,double>("H5:GUA:C5",-7.10892726493));
    this->alphas.insert(std::pair<std::string,double>("H5:GUA:C6",-3.65921586603));
    this->alphas.insert(std::pair<std::string,double>("H5:GUA:C8",-1.78019970189));
    this->alphas.insert(std::pair<std::string,double>("H5:GUA:N1",-2.80023020613));
    this->alphas.insert(std::pair<std::string,double>("H5:GUA:N2",6.79295378938));
    this->alphas.insert(std::pair<std::string,double>("H5:GUA:N3",3.23066286604));
    this->alphas.insert(std::pair<std::string,double>("H5:GUA:N7",-3.18580917165));
    this->alphas.insert(std::pair<std::string,double>("H5:GUA:N9",-4.43800178048));
    this->alphas.insert(std::pair<std::string,double>("H5:GUA:O6",5.66332648425));
    this->alphas.insert(std::pair<std::string,double>("H5:ADE:C1'",30.9955458732));
    this->alphas.insert(std::pair<std::string,double>("H5:ADE:C2'",-26.1209464829));
    this->alphas.insert(std::pair<std::string,double>("H5:ADE:C3'",-11.4258001662));
    this->alphas.insert(std::pair<std::string,double>("H5:ADE:C4'",-18.2729331255));
    this->alphas.insert(std::pair<std::string,double>("H5:ADE:C5'",9.60967236879));
    this->alphas.insert(std::pair<std::string,double>("H5:ADE:P",9.61799771197));
    this->alphas.insert(std::pair<std::string,double>("H5:ADE:O5'",18.6818451634));
    this->alphas.insert(std::pair<std::string,double>("H5:ADE:O3'",18.3685150052));
    this->alphas.insert(std::pair<std::string,double>("H5:ADE:C2",-2.71809366523));
    this->alphas.insert(std::pair<std::string,double>("H5:ADE:C4",-8.46864292171));
    this->alphas.insert(std::pair<std::string,double>("H5:ADE:C5",-3.32150035471));
    this->alphas.insert(std::pair<std::string,double>("H5:ADE:C6",2.31142583423));
    this->alphas.insert(std::pair<std::string,double>("H5:ADE:C8",-1.54418719293));
    this->alphas.insert(std::pair<std::string,double>("H5:ADE:N1",0.424499878363));
    this->alphas.insert(std::pair<std::string,double>("H5:ADE:N3",3.00401069893));
    this->alphas.insert(std::pair<std::string,double>("H5:ADE:N6",-3.35640604925));
    this->alphas.insert(std::pair<std::string,double>("H5:ADE:N7",-3.4790728058));
    this->alphas.insert(std::pair<std::string,double>("H5:ADE:N9",-10.3338738397));
    this->alphas.insert(std::pair<std::string,double>("H5:URA:C1'",1.41039102634));
    this->alphas.insert(std::pair<std::string,double>("H5:URA:C2'",6.78957512396));
    this->alphas.insert(std::pair<std::string,double>("H5:URA:C3'",-1.18134245297));
    this->alphas.insert(std::pair<std::string,double>("H5:URA:C4'",-8.33795737378));
    this->alphas.insert(std::pair<std::string,double>("H5:URA:C5'",7.20836242285));
    this->alphas.insert(std::pair<std::string,double>("H5:URA:P",-1.06492991733));
    this->alphas.insert(std::pair<std::string,double>("H5:URA:O5'",0.219470504227));
    this->alphas.insert(std::pair<std::string,double>("H5:URA:O3'",6.28036303578));
    this->alphas.insert(std::pair<std::string,double>("H5:URA:C2",-2.48085097083));
    this->alphas.insert(std::pair<std::string,double>("H5:URA:C4",-1.51388939805));
    this->alphas.insert(std::pair<std::string,double>("H5:URA:C5",0.0));
    this->alphas.insert(std::pair<std::string,double>("H5:URA:C6",-1.08246145427));
    this->alphas.insert(std::pair<std::string,double>("H5:URA:N1",1.29744038724));
    this->alphas.insert(std::pair<std::string,double>("H5:URA:N3",0.0797087574651));
    this->alphas.insert(std::pair<std::string,double>("H5:URA:O4",-0.781669192046));
    this->alphas.insert(std::pair<std::string,double>("H5:CYT:C1'",1.35733798405));
    this->alphas.insert(std::pair<std::string,double>("H5:CYT:C2'",1.07862490209));
    this->alphas.insert(std::pair<std::string,double>("H5:CYT:C3'",-3.89790387925));
    this->alphas.insert(std::pair<std::string,double>("H5:CYT:C4'",6.162170197));
    this->alphas.insert(std::pair<std::string,double>("H5:CYT:C5'",16.315493606));
    this->alphas.insert(std::pair<std::string,double>("H5:CYT:P",-10.7684961704));
    this->alphas.insert(std::pair<std::string,double>("H5:CYT:O5'",-4.46958202252));
    this->alphas.insert(std::pair<std::string,double>("H5:CYT:O3'",0.615786764817));
    this->alphas.insert(std::pair<std::string,double>("H5:CYT:C2",-3.41771410735));
    this->alphas.insert(std::pair<std::string,double>("H5:CYT:C4",-0.583805143394));
    this->alphas.insert(std::pair<std::string,double>("H5:CYT:C5",-0.0300448003023));
    this->alphas.insert(std::pair<std::string,double>("H5:CYT:C6",-1.87953523909));
    this->alphas.insert(std::pair<std::string,double>("H5:CYT:N1",-1.28006908324));
    this->alphas.insert(std::pair<std::string,double>("H5:CYT:N3",-0.560651207854));
    this->alphas.insert(std::pair<std::string,double>("H5:CYT:N4",-1.32780847621));
    this->alphas.insert(std::pair<std::string,double>("H5:CYT:O2",-0.264565186748));
    this->alphas.insert(std::pair<std::string,double>("C8:GUA:C1'",44.5836909036));
    this->alphas.insert(std::pair<std::string,double>("C8:GUA:C2'",-20.4995752332));
    this->alphas.insert(std::pair<std::string,double>("C8:GUA:C3'",-117.507595367));
    this->alphas.insert(std::pair<std::string,double>("C8:GUA:C4'",-0.697878358411));
    this->alphas.insert(std::pair<std::string,double>("C8:GUA:C5'",30.6329616365));
    this->alphas.insert(std::pair<std::string,double>("C8:GUA:P",31.9196493542));
    this->alphas.insert(std::pair<std::string,double>("C8:GUA:O5'",-33.8697598383));
    this->alphas.insert(std::pair<std::string,double>("C8:GUA:O3'",68.4786809726));
    this->alphas.insert(std::pair<std::string,double>("C8:GUA:C2",-37.3301009036));
    this->alphas.insert(std::pair<std::string,double>("C8:GUA:C4",8.12333897377));
    this->alphas.insert(std::pair<std::string,double>("C8:GUA:C5",7.76948561507));
    this->alphas.insert(std::pair<std::string,double>("C8:GUA:C6",-11.3280238481));
    this->alphas.insert(std::pair<std::string,double>("C8:GUA:C8",-33.0385995825));
    this->alphas.insert(std::pair<std::string,double>("C8:GUA:N1",-24.8369930374));
    this->alphas.insert(std::pair<std::string,double>("C8:GUA:N2",49.4841047282));
    this->alphas.insert(std::pair<std::string,double>("C8:GUA:N3",-10.3826215332));
    this->alphas.insert(std::pair<std::string,double>("C8:GUA:N7",1.78177372878));
    this->alphas.insert(std::pair<std::string,double>("C8:GUA:N9",1.12918894151));
    this->alphas.insert(std::pair<std::string,double>("C8:GUA:O6",-2.73346245549));
    this->alphas.insert(std::pair<std::string,double>("C8:ADE:C1'",34.3921682982));
    this->alphas.insert(std::pair<std::string,double>("C8:ADE:C2'",-3.29602701258));
    this->alphas.insert(std::pair<std::string,double>("C8:ADE:C3'",-46.0546472291));
    this->alphas.insert(std::pair<std::string,double>("C8:ADE:C4'",1.24652297363));
    this->alphas.insert(std::pair<std::string,double>("C8:ADE:C5'",-14.2723263023));
    this->alphas.insert(std::pair<std::string,double>("C8:ADE:P",-8.73876261467));
    this->alphas.insert(std::pair<std::string,double>("C8:ADE:O5'",4.83743393635));
    this->alphas.insert(std::pair<std::string,double>("C8:ADE:O3'",23.8483011795));
    this->alphas.insert(std::pair<std::string,double>("C8:ADE:C2",-33.7596652463));
    this->alphas.insert(std::pair<std::string,double>("C8:ADE:C4",5.4049704572));
    this->alphas.insert(std::pair<std::string,double>("C8:ADE:C5",5.57909571565));
    this->alphas.insert(std::pair<std::string,double>("C8:ADE:C6",4.43101804192));
    this->alphas.insert(std::pair<std::string,double>("C8:ADE:C8",-50.644828533));
    this->alphas.insert(std::pair<std::string,double>("C8:ADE:N1",7.29492027133));
    this->alphas.insert(std::pair<std::string,double>("C8:ADE:N3",-51.1547012056));
    this->alphas.insert(std::pair<std::string,double>("C8:ADE:N6",25.0100674445));
    this->alphas.insert(std::pair<std::string,double>("C8:ADE:N7",-0.642837099892));
    this->alphas.insert(std::pair<std::string,double>("C8:ADE:N9",2.42758464536));
    this->alphas.insert(std::pair<std::string,double>("C8:URA:C1'",81.7578486382));
    this->alphas.insert(std::pair<std::string,double>("C8:URA:C2'",17.4294108492));
    this->alphas.insert(std::pair<std::string,double>("C8:URA:C3'",-88.9191472766));
    this->alphas.insert(std::pair<std::string,double>("C8:URA:C4'",28.7364916186));
    this->alphas.insert(std::pair<std::string,double>("C8:URA:C5'",-40.4885933504));
    this->alphas.insert(std::pair<std::string,double>("C8:URA:P",-10.8021390336));
    this->alphas.insert(std::pair<std::string,double>("C8:URA:O5'",-9.45976341132));
    this->alphas.insert(std::pair<std::string,double>("C8:URA:O3'",49.5772211538));
    this->alphas.insert(std::pair<std::string,double>("C8:URA:C2",-16.7728371808));
    this->alphas.insert(std::pair<std::string,double>("C8:URA:C4",30.1408035608));
    this->alphas.insert(std::pair<std::string,double>("C8:URA:C5",31.508888927));
    this->alphas.insert(std::pair<std::string,double>("C8:URA:C6",-43.0261189302));
    this->alphas.insert(std::pair<std::string,double>("C8:URA:N1",1.33916456425));
    this->alphas.insert(std::pair<std::string,double>("C8:URA:N3",-67.8633706605));
    this->alphas.insert(std::pair<std::string,double>("C8:URA:O4",-14.5255924169));
    this->alphas.insert(std::pair<std::string,double>("C8:CYT:C1'",-14.2506021238));
    this->alphas.insert(std::pair<std::string,double>("C8:CYT:C2'",66.3818935843));
    this->alphas.insert(std::pair<std::string,double>("C8:CYT:C3'",-106.848152305));
    this->alphas.insert(std::pair<std::string,double>("C8:CYT:C4'",49.8185467829));
    this->alphas.insert(std::pair<std::string,double>("C8:CYT:C5'",140.080222459));
    this->alphas.insert(std::pair<std::string,double>("C8:CYT:P",-15.6311466019));
    this->alphas.insert(std::pair<std::string,double>("C8:CYT:O5'",-92.9275463373));
    this->alphas.insert(std::pair<std::string,double>("C8:CYT:O3'",-79.5908585758));
    this->alphas.insert(std::pair<std::string,double>("C8:CYT:C2",13.0050608114));
    this->alphas.insert(std::pair<std::string,double>("C8:CYT:C4",-27.1143823659));
    this->alphas.insert(std::pair<std::string,double>("C8:CYT:C5",4.98837073311));
    this->alphas.insert(std::pair<std::string,double>("C8:CYT:C6",37.9931162836));
    this->alphas.insert(std::pair<std::string,double>("C8:CYT:N1",15.7476512819));
    this->alphas.insert(std::pair<std::string,double>("C8:CYT:N3",-37.0543891699));
    this->alphas.insert(std::pair<std::string,double>("C8:CYT:N4",-64.1410009143));
    this->alphas.insert(std::pair<std::string,double>("C8:CYT:O2",-39.0763148731));
    this->alphas.insert(std::pair<std::string,double>("C2:GUA:C1'",-32.9807788245));
    this->alphas.insert(std::pair<std::string,double>("C2:GUA:C2'",100.280554763));
    this->alphas.insert(std::pair<std::string,double>("C2:GUA:C3'",7.42641155626));
    this->alphas.insert(std::pair<std::string,double>("C2:GUA:C4'",19.6982440956));
    this->alphas.insert(std::pair<std::string,double>("C2:GUA:C5'",25.9130237914));
    this->alphas.insert(std::pair<std::string,double>("C2:GUA:P",61.2244337979));
    this->alphas.insert(std::pair<std::string,double>("C2:GUA:O5'",6.64288430188));
    this->alphas.insert(std::pair<std::string,double>("C2:GUA:O3'",-81.8932477836));
    this->alphas.insert(std::pair<std::string,double>("C2:GUA:C2",-0.79870959906));
    this->alphas.insert(std::pair<std::string,double>("C2:GUA:C4",-4.6479002722));
    this->alphas.insert(std::pair<std::string,double>("C2:GUA:C5",-37.1263534285));
    this->alphas.insert(std::pair<std::string,double>("C2:GUA:C6",-5.68411803005));
    this->alphas.insert(std::pair<std::string,double>("C2:GUA:C8",-67.0271760218));
    this->alphas.insert(std::pair<std::string,double>("C2:GUA:N1",32.3964154536));
    this->alphas.insert(std::pair<std::string,double>("C2:GUA:N2",4.166294997));
    this->alphas.insert(std::pair<std::string,double>("C2:GUA:N3",3.55856574403));
    this->alphas.insert(std::pair<std::string,double>("C2:GUA:N7",-4.28685433872));
    this->alphas.insert(std::pair<std::string,double>("C2:GUA:N9",8.70037493757));
    this->alphas.insert(std::pair<std::string,double>("C2:GUA:O6",30.9923242491));
    this->alphas.insert(std::pair<std::string,double>("C2:ADE:C1'",22.005709057));
    this->alphas.insert(std::pair<std::string,double>("C2:ADE:C2'",22.9473951297));
    this->alphas.insert(std::pair<std::string,double>("C2:ADE:C3'",24.4977028163));
    this->alphas.insert(std::pair<std::string,double>("C2:ADE:C4'",34.6235509957));
    this->alphas.insert(std::pair<std::string,double>("C2:ADE:C5'",1.57070574818));
    this->alphas.insert(std::pair<std::string,double>("C2:ADE:P",-65.1384503611));
    this->alphas.insert(std::pair<std::string,double>("C2:ADE:O5'",-93.1025905854));
    this->alphas.insert(std::pair<std::string,double>("C2:ADE:O3'",136.936702582));
    this->alphas.insert(std::pair<std::string,double>("C2:ADE:C2",-1.61872208456));
    this->alphas.insert(std::pair<std::string,double>("C2:ADE:C4",0.92579195053));
    this->alphas.insert(std::pair<std::string,double>("C2:ADE:C5",-0.0764013184191));
    this->alphas.insert(std::pair<std::string,double>("C2:ADE:C6",3.97491046118));
    this->alphas.insert(std::pair<std::string,double>("C2:ADE:C8",-52.251130646));
    this->alphas.insert(std::pair<std::string,double>("C2:ADE:N1",1.14263231343));
    this->alphas.insert(std::pair<std::string,double>("C2:ADE:N3",0.20003701758));
    this->alphas.insert(std::pair<std::string,double>("C2:ADE:N6",18.0354506772));
    this->alphas.insert(std::pair<std::string,double>("C2:ADE:N7",-23.4826642866));
    this->alphas.insert(std::pair<std::string,double>("C2:ADE:N9",-21.2404133268));
    this->alphas.insert(std::pair<std::string,double>("C2:URA:C1'",-31.5671632289));
    this->alphas.insert(std::pair<std::string,double>("C2:URA:C2'",-39.5833190126));
    this->alphas.insert(std::pair<std::string,double>("C2:URA:C3'",17.9285088806));
    this->alphas.insert(std::pair<std::string,double>("C2:URA:C4'",-10.2746994221));
    this->alphas.insert(std::pair<std::string,double>("C2:URA:C5'",73.7430209149));
    this->alphas.insert(std::pair<std::string,double>("C2:URA:P",-18.2481272786));
    this->alphas.insert(std::pair<std::string,double>("C2:URA:O5'",-22.4652214107));
    this->alphas.insert(std::pair<std::string,double>("C2:URA:O3'",116.446276977));
    this->alphas.insert(std::pair<std::string,double>("C2:URA:C2",-12.2479894627));
    this->alphas.insert(std::pair<std::string,double>("C2:URA:C4",-6.64858952374));
    this->alphas.insert(std::pair<std::string,double>("C2:URA:C5",-2.69210436717));
    this->alphas.insert(std::pair<std::string,double>("C2:URA:C6",4.8235118848));
    this->alphas.insert(std::pair<std::string,double>("C2:URA:N1",24.1762733166));
    this->alphas.insert(std::pair<std::string,double>("C2:URA:N3",-16.2503020503));
    this->alphas.insert(std::pair<std::string,double>("C2:URA:O4",-6.99027921664));
    this->alphas.insert(std::pair<std::string,double>("C2:CYT:C1'",63.588160415));
    this->alphas.insert(std::pair<std::string,double>("C2:CYT:C2'",17.6662312072));
    this->alphas.insert(std::pair<std::string,double>("C2:CYT:C3'",-113.625102818));
    this->alphas.insert(std::pair<std::string,double>("C2:CYT:C4'",63.1870368469));
    this->alphas.insert(std::pair<std::string,double>("C2:CYT:C5'",-161.507199625));
    this->alphas.insert(std::pair<std::string,double>("C2:CYT:P",39.796423685));
    this->alphas.insert(std::pair<std::string,double>("C2:CYT:O5'",-6.53706291975));
    this->alphas.insert(std::pair<std::string,double>("C2:CYT:O3'",-96.3538521485));
    this->alphas.insert(std::pair<std::string,double>("C2:CYT:C2",14.4354812796));
    this->alphas.insert(std::pair<std::string,double>("C2:CYT:C4",-14.2471075668));
    this->alphas.insert(std::pair<std::string,double>("C2:CYT:C5",-4.75181771335));
    this->alphas.insert(std::pair<std::string,double>("C2:CYT:C6",-2.20116305901));
    this->alphas.insert(std::pair<std::string,double>("C2:CYT:N1",-10.3975596808));
    this->alphas.insert(std::pair<std::string,double>("C2:CYT:N3",24.7102359667));
    this->alphas.insert(std::pair<std::string,double>("C2:CYT:N4",-15.2565119241));
    this->alphas.insert(std::pair<std::string,double>("C2:CYT:O2",-34.9105915644));
    this->alphas.insert(std::pair<std::string,double>("C6:GUA:C1'",49.0262369405));
    this->alphas.insert(std::pair<std::string,double>("C6:GUA:C2'",-70.3926470372));
    this->alphas.insert(std::pair<std::string,double>("C6:GUA:C3'",-37.7569424084));
    this->alphas.insert(std::pair<std::string,double>("C6:GUA:C4'",-2.2719201484));
    this->alphas.insert(std::pair<std::string,double>("C6:GUA:C5'",25.9133964552));
    this->alphas.insert(std::pair<std::string,double>("C6:GUA:P",75.9731996192));
    this->alphas.insert(std::pair<std::string,double>("C6:GUA:O5'",-82.8518141642));
    this->alphas.insert(std::pair<std::string,double>("C6:GUA:O3'",50.883605954));
    this->alphas.insert(std::pair<std::string,double>("C6:GUA:C2",-1.42561004533));
    this->alphas.insert(std::pair<std::string,double>("C6:GUA:C4",3.54431594504));
    this->alphas.insert(std::pair<std::string,double>("C6:GUA:C5",-16.9692255309));
    this->alphas.insert(std::pair<std::string,double>("C6:GUA:C6",-15.4850897301));
    this->alphas.insert(std::pair<std::string,double>("C6:GUA:C8",-16.3509591815));
    this->alphas.insert(std::pair<std::string,double>("C6:GUA:N1",-19.5220251761));
    this->alphas.insert(std::pair<std::string,double>("C6:GUA:N2",3.17886913036));
    this->alphas.insert(std::pair<std::string,double>("C6:GUA:N3",5.34404554519));
    this->alphas.insert(std::pair<std::string,double>("C6:GUA:N7",-28.4108818465));
    this->alphas.insert(std::pair<std::string,double>("C6:GUA:N9",28.8191698418));
    this->alphas.insert(std::pair<std::string,double>("C6:GUA:O6",21.0042473357));
    this->alphas.insert(std::pair<std::string,double>("C6:ADE:C1'",49.0346171667));
    this->alphas.insert(std::pair<std::string,double>("C6:ADE:C2'",-75.4548652827));
    this->alphas.insert(std::pair<std::string,double>("C6:ADE:C3'",-30.4393215796));
    this->alphas.insert(std::pair<std::string,double>("C6:ADE:C4'",115.939339352));
    this->alphas.insert(std::pair<std::string,double>("C6:ADE:C5'",30.3226209489));
    this->alphas.insert(std::pair<std::string,double>("C6:ADE:P",-129.373039712));
    this->alphas.insert(std::pair<std::string,double>("C6:ADE:O5'",-11.420581455));
    this->alphas.insert(std::pair<std::string,double>("C6:ADE:O3'",10.1960749571));
    this->alphas.insert(std::pair<std::string,double>("C6:ADE:C2",-8.35072903732));
    this->alphas.insert(std::pair<std::string,double>("C6:ADE:C4",-22.6175823206));
    this->alphas.insert(std::pair<std::string,double>("C6:ADE:C5",-1.18103414415));
    this->alphas.insert(std::pair<std::string,double>("C6:ADE:C6",-5.79009572491));
    this->alphas.insert(std::pair<std::string,double>("C6:ADE:C8",2.14389961275));
    this->alphas.insert(std::pair<std::string,double>("C6:ADE:N1",21.551045154));
    this->alphas.insert(std::pair<std::string,double>("C6:ADE:N3",-43.4466396718));
    this->alphas.insert(std::pair<std::string,double>("C6:ADE:N6",1.93812342348));
    this->alphas.insert(std::pair<std::string,double>("C6:ADE:N7",-5.2926036866));
    this->alphas.insert(std::pair<std::string,double>("C6:ADE:N9",30.3167461698));
    this->alphas.insert(std::pair<std::string,double>("C6:URA:C1'",5.85689730713));
    this->alphas.insert(std::pair<std::string,double>("C6:URA:C2'",24.0889300202));
    this->alphas.insert(std::pair<std::string,double>("C6:URA:C3'",-66.2695328947));
    this->alphas.insert(std::pair<std::string,double>("C6:URA:C4'",18.935228391));
    this->alphas.insert(std::pair<std::string,double>("C6:URA:C5'",55.9917463903));
    this->alphas.insert(std::pair<std::string,double>("C6:URA:P",42.5633071009));
    this->alphas.insert(std::pair<std::string,double>("C6:URA:O5'",-25.2838030211));
    this->alphas.insert(std::pair<std::string,double>("C6:URA:O3'",23.7457273192));
    this->alphas.insert(std::pair<std::string,double>("C6:URA:C2",1.02023901636));
    this->alphas.insert(std::pair<std::string,double>("C6:URA:C4",2.28285925841));
    this->alphas.insert(std::pair<std::string,double>("C6:URA:C5",0.349109558141));
    this->alphas.insert(std::pair<std::string,double>("C6:URA:C6",-44.8033816552));
    this->alphas.insert(std::pair<std::string,double>("C6:URA:N1",-0.0));
    this->alphas.insert(std::pair<std::string,double>("C6:URA:N3",-1.18309475013));
    this->alphas.insert(std::pair<std::string,double>("C6:URA:O4",4.65496651495));
    this->alphas.insert(std::pair<std::string,double>("C6:CYT:C1'",11.7427652254));
    this->alphas.insert(std::pair<std::string,double>("C6:CYT:C2'",-3.80144284162));
    this->alphas.insert(std::pair<std::string,double>("C6:CYT:C3'",-51.7352686786));
    this->alphas.insert(std::pair<std::string,double>("C6:CYT:C4'",6.33367441138));
    this->alphas.insert(std::pair<std::string,double>("C6:CYT:C5'",80.768478577));
    this->alphas.insert(std::pair<std::string,double>("C6:CYT:P",5.85126619552));
    this->alphas.insert(std::pair<std::string,double>("C6:CYT:O5'",-10.883411394));
    this->alphas.insert(std::pair<std::string,double>("C6:CYT:O3'",-10.4079749642));
    this->alphas.insert(std::pair<std::string,double>("C6:CYT:C2",0.954981214108));
    this->alphas.insert(std::pair<std::string,double>("C6:CYT:C4",2.41800906721));
    this->alphas.insert(std::pair<std::string,double>("C6:CYT:C5",0.127509708064));
    this->alphas.insert(std::pair<std::string,double>("C6:CYT:C6",-21.1752142054));
    this->alphas.insert(std::pair<std::string,double>("C6:CYT:N1",0.1152975529));
    this->alphas.insert(std::pair<std::string,double>("C6:CYT:N3",-1.04437655782));
    this->alphas.insert(std::pair<std::string,double>("C6:CYT:N4",-14.8995817002));
    this->alphas.insert(std::pair<std::string,double>("C6:CYT:O2",2.47657965074));
    this->alphas.insert(std::pair<std::string,double>("C5:GUA:C1'",24.593911474));
    this->alphas.insert(std::pair<std::string,double>("C5:GUA:C2'",-75.282270576));
    this->alphas.insert(std::pair<std::string,double>("C5:GUA:C3'",-78.1636680519));
    this->alphas.insert(std::pair<std::string,double>("C5:GUA:C4'",30.9682728125));
    this->alphas.insert(std::pair<std::string,double>("C5:GUA:C5'",96.1140588603));
    this->alphas.insert(std::pair<std::string,double>("C5:GUA:P",52.4276691491));
    this->alphas.insert(std::pair<std::string,double>("C5:GUA:O5'",-27.0274676105));
    this->alphas.insert(std::pair<std::string,double>("C5:GUA:O3'",-41.3197959468));
    this->alphas.insert(std::pair<std::string,double>("C5:GUA:C2",-8.1618947817));
    this->alphas.insert(std::pair<std::string,double>("C5:GUA:C4",19.7370997353));
    this->alphas.insert(std::pair<std::string,double>("C5:GUA:C5",15.2929550907));
    this->alphas.insert(std::pair<std::string,double>("C5:GUA:C6",0.753212563409));
    this->alphas.insert(std::pair<std::string,double>("C5:GUA:C8",-38.958014046));
    this->alphas.insert(std::pair<std::string,double>("C5:GUA:N1",-85.7215474874));
    this->alphas.insert(std::pair<std::string,double>("C5:GUA:N2",-8.80753963559));
    this->alphas.insert(std::pair<std::string,double>("C5:GUA:N3",47.0096018819));
    this->alphas.insert(std::pair<std::string,double>("C5:GUA:N7",10.9322023381));
    this->alphas.insert(std::pair<std::string,double>("C5:GUA:N9",-22.8355322127));
    this->alphas.insert(std::pair<std::string,double>("C5:GUA:O6",1.10242747223));
    this->alphas.insert(std::pair<std::string,double>("C5:ADE:C1'",-27.6093794852));
    this->alphas.insert(std::pair<std::string,double>("C5:ADE:C2'",-37.8285572835));
    this->alphas.insert(std::pair<std::string,double>("C5:ADE:C3'",-55.3795978436));
    this->alphas.insert(std::pair<std::string,double>("C5:ADE:C4'",-0.152016580486));
    this->alphas.insert(std::pair<std::string,double>("C5:ADE:C5'",26.3251158306));
    this->alphas.insert(std::pair<std::string,double>("C5:ADE:P",77.5961992304));
    this->alphas.insert(std::pair<std::string,double>("C5:ADE:O5'",69.9227778171));
    this->alphas.insert(std::pair<std::string,double>("C5:ADE:O3'",-31.9531855583));
    this->alphas.insert(std::pair<std::string,double>("C5:ADE:C2",23.8801308624));
    this->alphas.insert(std::pair<std::string,double>("C5:ADE:C4",7.8868974429));
    this->alphas.insert(std::pair<std::string,double>("C5:ADE:C5",0.227915764595));
    this->alphas.insert(std::pair<std::string,double>("C5:ADE:C6",-20.1984026786));
    this->alphas.insert(std::pair<std::string,double>("C5:ADE:C8",-3.70834974414));
    this->alphas.insert(std::pair<std::string,double>("C5:ADE:N1",-9.87726352451));
    this->alphas.insert(std::pair<std::string,double>("C5:ADE:N3",11.9331759598));
    this->alphas.insert(std::pair<std::string,double>("C5:ADE:N6",-58.6139621643));
    this->alphas.insert(std::pair<std::string,double>("C5:ADE:N7",-3.50889163414));
    this->alphas.insert(std::pair<std::string,double>("C5:ADE:N9",-2.73016060981));
    this->alphas.insert(std::pair<std::string,double>("C5:URA:C1'",1.78969552925));
    this->alphas.insert(std::pair<std::string,double>("C5:URA:C2'",-62.6199245926));
    this->alphas.insert(std::pair<std::string,double>("C5:URA:C3'",-64.0533854471));
    this->alphas.insert(std::pair<std::string,double>("C5:URA:C4'",10.9457122978));
    this->alphas.insert(std::pair<std::string,double>("C5:URA:C5'",-16.0525734919));
    this->alphas.insert(std::pair<std::string,double>("C5:URA:P",94.8153756517));
    this->alphas.insert(std::pair<std::string,double>("C5:URA:O5'",-8.59217308486));
    this->alphas.insert(std::pair<std::string,double>("C5:URA:O3'",21.4928097155));
    this->alphas.insert(std::pair<std::string,double>("C5:URA:C2",3.16361626679));
    this->alphas.insert(std::pair<std::string,double>("C5:URA:C4",1.7614091899));
    this->alphas.insert(std::pair<std::string,double>("C5:URA:C5",-31.0903874324));
    this->alphas.insert(std::pair<std::string,double>("C5:URA:C6",1.53742590466));
    this->alphas.insert(std::pair<std::string,double>("C5:URA:N1",8.65730423799));
    this->alphas.insert(std::pair<std::string,double>("C5:URA:N3",7.29672365987));
    this->alphas.insert(std::pair<std::string,double>("C5:URA:O4",11.2491331868));
    this->alphas.insert(std::pair<std::string,double>("C5:CYT:C1'",80.5514386717));
    this->alphas.insert(std::pair<std::string,double>("C5:CYT:C2'",184.314401196));
    this->alphas.insert(std::pair<std::string,double>("C5:CYT:C3'",-179.702276224));
    this->alphas.insert(std::pair<std::string,double>("C5:CYT:C4'",-9.07448098252));
    this->alphas.insert(std::pair<std::string,double>("C5:CYT:C5'",200.060772639));
    this->alphas.insert(std::pair<std::string,double>("C5:CYT:P",-33.0382485855));
    this->alphas.insert(std::pair<std::string,double>("C5:CYT:O5'",-191.121444667));
    this->alphas.insert(std::pair<std::string,double>("C5:CYT:O3'",-18.386965361));
    this->alphas.insert(std::pair<std::string,double>("C5:CYT:C2",-2.78213004972));
    this->alphas.insert(std::pair<std::string,double>("C5:CYT:C4",1.08479962586));
    this->alphas.insert(std::pair<std::string,double>("C5:CYT:C5",-46.6666595343));
    this->alphas.insert(std::pair<std::string,double>("C5:CYT:C6",0.505706111039));
    this->alphas.insert(std::pair<std::string,double>("C5:CYT:N1",12.2757220355));
    this->alphas.insert(std::pair<std::string,double>("C5:CYT:N3",0.389718934631));
    this->alphas.insert(std::pair<std::string,double>("C5:CYT:N4",15.9492324187));
    this->alphas.insert(std::pair<std::string,double>("C5:CYT:O2",-48.2342535574));
}

void LARMORD::initializeAtomTypes()
{
	if(this->MolTypeLamor==Misc::tolower("protein")){
		this->atomTypes.clear();
		this->atomTypes.push_back("CYS:C");
		this->atomTypes.push_back("CYS:CB");
		this->atomTypes.push_back("CYS:CA");
		this->atomTypes.push_back("CYS:O");
		this->atomTypes.push_back("CYS:N");
		this->atomTypes.push_back("CYS:SG");
		this->atomTypes.push_back("ASP:C");
		this->atomTypes.push_back("ASP:CB");
		this->atomTypes.push_back("ASP:CA");
		this->atomTypes.push_back("ASP:CG");
		this->atomTypes.push_back("ASP:O");
		this->atomTypes.push_back("ASP:N");
		this->atomTypes.push_back("ASP:OD1");
		this->atomTypes.push_back("ASP:OD2");
		this->atomTypes.push_back("SER:C");
		this->atomTypes.push_back("SER:OG");
		this->atomTypes.push_back("SER:CB");
		this->atomTypes.push_back("SER:CA");
		this->atomTypes.push_back("SER:O");
		this->atomTypes.push_back("SER:N");
		this->atomTypes.push_back("GLN:C");
		this->atomTypes.push_back("GLN:CB");
		this->atomTypes.push_back("GLN:CA");
		this->atomTypes.push_back("GLN:CG");
		this->atomTypes.push_back("GLN:O");
		this->atomTypes.push_back("GLN:N");
		this->atomTypes.push_back("GLN:CD");
		this->atomTypes.push_back("GLN:NE2");
		this->atomTypes.push_back("GLN:OE1");
		this->atomTypes.push_back("LYS:C");
		this->atomTypes.push_back("LYS:CB");
		this->atomTypes.push_back("LYS:CA");
		this->atomTypes.push_back("LYS:CG");
		this->atomTypes.push_back("LYS:CE");
		this->atomTypes.push_back("LYS:CD");
		this->atomTypes.push_back("LYS:NZ");
		this->atomTypes.push_back("LYS:O");
		this->atomTypes.push_back("LYS:N");
		this->atomTypes.push_back("ILE:C");
		this->atomTypes.push_back("ILE:CG2");
		this->atomTypes.push_back("ILE:CB");
		this->atomTypes.push_back("ILE:CA");
		this->atomTypes.push_back("ILE:O");
		this->atomTypes.push_back("ILE:N");
		this->atomTypes.push_back("ILE:CD1");
		this->atomTypes.push_back("ILE:CG1");
		this->atomTypes.push_back("ILE:CD");
		this->atomTypes.push_back("PRO:C");
		this->atomTypes.push_back("PRO:CB");
		this->atomTypes.push_back("PRO:CA");
		this->atomTypes.push_back("PRO:CG");
		this->atomTypes.push_back("PRO:O");
		this->atomTypes.push_back("PRO:N");
		this->atomTypes.push_back("PRO:CD");
		this->atomTypes.push_back("THR:C");
		this->atomTypes.push_back("THR:CB");
		this->atomTypes.push_back("THR:CA");
		this->atomTypes.push_back("THR:OG1");
		this->atomTypes.push_back("THR:O");
		this->atomTypes.push_back("THR:N");
		this->atomTypes.push_back("THR:CG2");
		this->atomTypes.push_back("PHE:C");
		this->atomTypes.push_back("PHE:CE1");
		this->atomTypes.push_back("PHE:CB");
		this->atomTypes.push_back("PHE:CA");
		this->atomTypes.push_back("PHE:CG");
		this->atomTypes.push_back("PHE:O");
		this->atomTypes.push_back("PHE:N");
		this->atomTypes.push_back("PHE:CZ");
		this->atomTypes.push_back("PHE:CD1");
		this->atomTypes.push_back("PHE:CD2");
		this->atomTypes.push_back("PHE:CE2");
		this->atomTypes.push_back("ALA:CB");
		this->atomTypes.push_back("ALA:C");
		this->atomTypes.push_back("ALA:CA");
		this->atomTypes.push_back("ALA:O");
		this->atomTypes.push_back("ALA:N");
		this->atomTypes.push_back("GLY:C");
		this->atomTypes.push_back("GLY:CA");
		this->atomTypes.push_back("GLY:O");
		this->atomTypes.push_back("GLY:N");
		this->atomTypes.push_back("HIS:C");
		this->atomTypes.push_back("HIS:CE1");
		this->atomTypes.push_back("HIS:CB");
		this->atomTypes.push_back("HIS:CA");
		this->atomTypes.push_back("HIS:CG");
		this->atomTypes.push_back("HIS:O");
		this->atomTypes.push_back("HIS:N");
		this->atomTypes.push_back("HIS:CD2");
		this->atomTypes.push_back("HIS:ND1");
		this->atomTypes.push_back("HIS:NE2");
		this->atomTypes.push_back("GLU:C");
		this->atomTypes.push_back("GLU:CB");
		this->atomTypes.push_back("GLU:CA");
		this->atomTypes.push_back("GLU:CG");
		this->atomTypes.push_back("GLU:O");
		this->atomTypes.push_back("GLU:N");
		this->atomTypes.push_back("GLU:OE2");
		this->atomTypes.push_back("GLU:CD");
		this->atomTypes.push_back("GLU:OE1");
		this->atomTypes.push_back("LEU:C");
		this->atomTypes.push_back("LEU:CB");
		this->atomTypes.push_back("LEU:CA");
		this->atomTypes.push_back("LEU:CG");
		this->atomTypes.push_back("LEU:O");
		this->atomTypes.push_back("LEU:N");
		this->atomTypes.push_back("LEU:CD1");
		this->atomTypes.push_back("LEU:CD2");
		this->atomTypes.push_back("ARG:C");
		this->atomTypes.push_back("ARG:CB");
		this->atomTypes.push_back("ARG:CA");
		this->atomTypes.push_back("ARG:CG");
		this->atomTypes.push_back("ARG:NE");
		this->atomTypes.push_back("ARG:O");
		this->atomTypes.push_back("ARG:CD");
		this->atomTypes.push_back("ARG:CZ");
		this->atomTypes.push_back("ARG:NH1");
		this->atomTypes.push_back("ARG:NH2");
		this->atomTypes.push_back("ARG:N");
		this->atomTypes.push_back("TRP:C");
		this->atomTypes.push_back("TRP:CD1");
		this->atomTypes.push_back("TRP:CZ2");
		this->atomTypes.push_back("TRP:CB");
		this->atomTypes.push_back("TRP:CA");
		this->atomTypes.push_back("TRP:CG");
		this->atomTypes.push_back("TRP:O");
		this->atomTypes.push_back("TRP:N");
		this->atomTypes.push_back("TRP:CH2");
		this->atomTypes.push_back("TRP:CE3");
		this->atomTypes.push_back("TRP:CE2");
		this->atomTypes.push_back("TRP:CD2");
		this->atomTypes.push_back("TRP:CZ3");
		this->atomTypes.push_back("TRP:NE1");
		this->atomTypes.push_back("VAL:C");
		this->atomTypes.push_back("VAL:CB");
		this->atomTypes.push_back("VAL:CA");
		this->atomTypes.push_back("VAL:O");
		this->atomTypes.push_back("VAL:N");
		this->atomTypes.push_back("VAL:CG1");
		this->atomTypes.push_back("VAL:CG2");
		this->atomTypes.push_back("ASN:C");
		this->atomTypes.push_back("ASN:CB");
		this->atomTypes.push_back("ASN:CA");
		this->atomTypes.push_back("ASN:CG");
		this->atomTypes.push_back("ASN:O");
		this->atomTypes.push_back("ASN:N");
		this->atomTypes.push_back("ASN:OD1");
		this->atomTypes.push_back("ASN:OT1");
		this->atomTypes.push_back("ASN:ND2");
		this->atomTypes.push_back("TYR:C");
		this->atomTypes.push_back("TYR:CE1");
		this->atomTypes.push_back("TYR:OH");
		this->atomTypes.push_back("TYR:CB");
		this->atomTypes.push_back("TYR:CA");
		this->atomTypes.push_back("TYR:CG");
		this->atomTypes.push_back("TYR:O");
		this->atomTypes.push_back("TYR:N");
		this->atomTypes.push_back("TYR:CZ");
		this->atomTypes.push_back("TYR:CD1");
		this->atomTypes.push_back("TYR:CD2");
		this->atomTypes.push_back("TYR:CE2");
		this->atomTypes.push_back("MET:C");
		this->atomTypes.push_back("MET:CB");
		this->atomTypes.push_back("MET:CA");
		this->atomTypes.push_back("MET:CG");
		this->atomTypes.push_back("MET:CE");
		this->atomTypes.push_back("MET:N");
		this->atomTypes.push_back("MET:O");
		this->atomTypes.push_back("MET:SD");
	}
	else if(this->MolTypeLamor==Misc::tolower("rna")){
		this->atomTypes.push_back("GUA:C1'");
		this->atomTypes.push_back("GUA:C2'");
		this->atomTypes.push_back("GUA:C3'");
		this->atomTypes.push_back("GUA:C4'");
		this->atomTypes.push_back("GUA:C5'");
		this->atomTypes.push_back("GUA:P");
		this->atomTypes.push_back("GUA:O5'");
		this->atomTypes.push_back("GUA:O3'");
		this->atomTypes.push_back("GUA:C2");
		this->atomTypes.push_back("GUA:C4");
		this->atomTypes.push_back("GUA:C5");
		this->atomTypes.push_back("GUA:C6");
		this->atomTypes.push_back("GUA:C8");
		this->atomTypes.push_back("GUA:N1");
		this->atomTypes.push_back("GUA:N2");
		this->atomTypes.push_back("GUA:N3");
		this->atomTypes.push_back("GUA:N7");
		this->atomTypes.push_back("GUA:N9");
		this->atomTypes.push_back("GUA:O6");
		this->atomTypes.push_back("ADE:C1'");
		this->atomTypes.push_back("ADE:C2'");
		this->atomTypes.push_back("ADE:C3'");
		this->atomTypes.push_back("ADE:C4'");
		this->atomTypes.push_back("ADE:C5'");
		this->atomTypes.push_back("ADE:P");
		this->atomTypes.push_back("ADE:O5'");
		this->atomTypes.push_back("ADE:O3'");
		this->atomTypes.push_back("ADE:C2");
		this->atomTypes.push_back("ADE:C4");
		this->atomTypes.push_back("ADE:C5");
		this->atomTypes.push_back("ADE:C6");
		this->atomTypes.push_back("ADE:C8");
		this->atomTypes.push_back("ADE:N1");
		this->atomTypes.push_back("ADE:N3");
		this->atomTypes.push_back("ADE:N6");
		this->atomTypes.push_back("ADE:N7");
		this->atomTypes.push_back("ADE:N9");
		this->atomTypes.push_back("URA:C1'");
		this->atomTypes.push_back("URA:C2'");
		this->atomTypes.push_back("URA:C3'");
		this->atomTypes.push_back("URA:C4'");
		this->atomTypes.push_back("URA:C5'");
		this->atomTypes.push_back("URA:P");
		this->atomTypes.push_back("URA:O5'");
		this->atomTypes.push_back("URA:O3'");
		this->atomTypes.push_back("URA:C2");
		this->atomTypes.push_back("URA:C4");
		this->atomTypes.push_back("URA:C5");
		this->atomTypes.push_back("URA:C6");
		this->atomTypes.push_back("URA:N1");
		this->atomTypes.push_back("URA:N3");
		this->atomTypes.push_back("URA:O4");
		this->atomTypes.push_back("CYT:C1'");
		this->atomTypes.push_back("CYT:C2'");
		this->atomTypes.push_back("CYT:C3'");
		this->atomTypes.push_back("CYT:C4'");
		this->atomTypes.push_back("CYT:C5'");
		this->atomTypes.push_back("CYT:P");
		this->atomTypes.push_back("CYT:O5'");
		this->atomTypes.push_back("CYT:O3'");
		this->atomTypes.push_back("CYT:C2");
		this->atomTypes.push_back("CYT:C4");
		this->atomTypes.push_back("CYT:C5");
		this->atomTypes.push_back("CYT:C6");
		this->atomTypes.push_back("CYT:N1");
		this->atomTypes.push_back("CYT:N3");
		this->atomTypes.push_back("CYT:N4");
		this->atomTypes.push_back("CYT:O2");	
	}
}

double LARMORD::ellintk (double m)
{
    double mr1 = 1.0 - m;
    double mr2 = mr1*mr1;
    double mr3 = mr2*mr1;
    double mr4 = mr3*mr1;
    double ad = 1.38629436112 + 0.09666344259*mr1 + 0.03590092383*mr2 +
               0.03742563713*mr3 + 0.01451196212*mr4;
    double bd = 0.5 + 0.12498593597*mr1 + 0.06880248576*mr2 +
               0.03328355346*mr3 + 0.00441787012*mr4;
    return (ad + bd * log (1.0 / mr1));
}

double LARMORD::ellinte (double m)
{
    /* ellinte, calculates the E factor of the elliptical integral
       expects 0<=m<1 as input. This function was difined from formula 17.3.36
       Abramowitz and Segun, Handbook of Mathematical Functions, 
       Dover publications 1965
    */
    double mr1 = 1.0 - m;
    double mr2 = mr1*mr1;
    double mr3 = mr2*mr1;
    double mr4 = mr3*mr1;
    double ad = 1.0 + 0.44325141463*mr1 + 0.06260601220*mr2 +
    0.04757383546*mr3 + 0.01736506451*mr4;
    double bd = 0.24998368310*mr1 + 0.09200180037*mr2 +
    0.04069697526*mr3 + 0.00526449639*mr4;
    return (ad + bd * log (1.0 / mr1));
}

Coor LARMORD::base_plane_normal(Molecule *base)
{
    /* check if number of residues is get than one */
    Coor a;
    Coor b;
    a = base->getAtom(2)->getCoor() - base->getAtom(0)->getCoor();
    b = base->getAtom(4)->getCoor() - base->getAtom(0)->getCoor();
    return (b.cross(a));
}

double LARMORD::ringCurrentCompute(const Coor point)
{
  double rc = 0.0;
  double intensity;
  Coor com1;
  Coor norm;
  Coor f;
  /* loop over all ring */
  for (unsigned int i=0;i<this->ringCenterOfMass.size();i++){   
    /* get vector between point and com of ring */
    com1 = this->ringCenterOfMass.at(i);
    norm = this->ringNormal.at(i);
    intensity = this->ringIntensity.at(i);    
    f = point-com1;
    
    //std::cout << sqrt(f.dot(f)) << std::endl;
    if(sqrt(f.dot(f)) > this->cutoffLarmorRing) continue;
    
    rc += intensity*jb_geo(f,norm); 
  }
  return (rc);  
}

double LARMORD::jb_geo(const Coor point, const Coor normal)
{
  /* z, distance above the plane */
  /* rho, projected distance on plane */
  double jb;
  double k;
  double z,z2;
  double rho,rho2;

  z = fabs(point.dot(normal))+(ringSeparation/2);
  z2 = z*z;
  
  rho = (point-(normal*(point.dot(normal)))).norm();
  rho2 = rho*rho;
  k = sqrt((4*rho/(pow((1+rho),2)+z2)));
  jb = -1*(sqrt((1/(pow((1+rho),2)+z2)))*(this->ellintk(k)+this->ellinte(k)*((1-rho2-z2)/(pow((1-rho),2)+z2))));
  
  //std::cout << "z = " << z << std::endl;
  //std::cout << "rho = " << rho << std::endl;
  //std::cout << "k = " << k << std::endl;
  //std::cout << "ellintk = " << this->ellintk(k) << std::endl;
  //std::cout << "ellinte = " << this->ellinte(k) << std::endl;
  //create a unit test
  //coor1 =  Coor(1.000,2.000,3.000);
    //coor2 =  Coor(4.000,5.000,6.000);
    //this->jb_geo(coor1,coor2)      
    //answer should be: -0.000310637
  return (jb);
}

void LARMORD::setUpRings(Molecule *mol)
{
  std::stringstream resid;
  std::string sele1,sele2;
  std::string resname;  
  Molecule *tmpmol, *cmol;
  Residue *res;
  Coor coor, cog; 
  double rc;
  for (unsigned int h=0; h < this->aromaticsResKeys.size();h++){
    sele1 = ":"+this->aromaticsResKeys.at(h)+".";
    mol->select(sele1,false,false);    
    if (mol->getNAtomSelected()>0){
      tmpmol= mol->copy();
      for (unsigned int i=0; i < tmpmol->getResVecSize();i++){
        resid.str(""); //Clear stringstream
        res=tmpmol->getResidue(i);
        resid << res->getResId();
        resname = res->getResName();
        rc = this->getRingIntensity(this->aromaticsResKeys.at(h)+":RING6");
        if(rc!=0.0){
          sele2 = ":"+resid.str()+".RING6";
          tmpmol->select(sele2);
          cmol= tmpmol->copy();
          this->ringNormal.push_back(this->base_plane_normal(cmol));
          this->ringCenterOfMass.push_back(Analyze::centerOfGeometry(cmol));
          this->ringIntensity.push_back(rc);  
        }
        rc = this->getRingIntensity(this->aromaticsResKeys.at(h)+":RING5");     
        if(rc!=0.0){
          sele2 = ":"+resid.str()+".RING5";
          tmpmol->select(sele2);
          cmol= tmpmol->copy();
          this->ringNormal.push_back(this->base_plane_normal(cmol));
          this->ringCenterOfMass.push_back(Analyze::centerOfGeometry(cmol));
          this->ringIntensity.push_back(rc);  
        }
      }
    } 
  }
}

void LARMORD::setUpAromaticsResKeys()
{
  this->aromaticsResKeys.push_back("ADENINES");
  this->aromaticsResKeys.push_back("URIDINES");
  this->aromaticsResKeys.push_back("CYTOSINES");
  this->aromaticsResKeys.push_back("GUANINES");
  this->aromaticsResKeys.push_back("THYMINES");
  this->aromaticsResKeys.push_back("PHENYLALANINES");
  this->aromaticsResKeys.push_back("TYROSINES");
  this->aromaticsResKeys.push_back("TRYPTOPHANS");
  this->aromaticsResKeys.push_back("HISTIDINES");
}

void LARMORD::setUpRingIntensityMap()
{
  //this->ringIntensityMap.insert(std::pair<std::string,double>("ADENINES:RING6",0.90));
  //this->ringIntensityMap.insert(std::pair<std::string,double>("URIDINES:RING6",0.11));
  //this->ringIntensityMap.insert(std::pair<std::string,double>("CYTOSINES:RING6",0.28));
  //this->ringIntensityMap.insert(std::pair<std::string,double>("GUANINES:RING6",0.30));
  //this->ringIntensityMap.insert(std::pair<std::string,double>("THYMINES:RING6",0.11));
  //this->ringIntensityMap.insert(std::pair<std::string,double>("PHENYLALANINES:RING6",1.00));
  //this->ringIntensityMap.insert(std::pair<std::string,double>("TYROSINES:RING6",0.94));
  //this->ringIntensityMap.insert(std::pair<std::string,double>("TRYPTOPHANS:RING6",1.04));
  //this->ringIntensityMap.insert(std::pair<std::string,double>("HISTIDINES:RING6",0.53));    
  //this->ringIntensityMap.insert(std::pair<std::string,double>("ADENINES:RING5",0.66));
  //this->ringIntensityMap.insert(std::pair<std::string,double>("GUANINES:RING5",0.64));
  //this->ringIntensityMap.insert(std::pair<std::string,double>("TRYPTOPHANS:RING5",0.56));   

  /* Case, 1995 */  
  this->ringIntensityMap.insert(std::pair<std::string,double>("ADENINES:RING6",0.83));
  this->ringIntensityMap.insert(std::pair<std::string,double>("URIDINES:RING6",0.24));
  this->ringIntensityMap.insert(std::pair<std::string,double>("CYTOSINES:RING6",0.31));
  this->ringIntensityMap.insert(std::pair<std::string,double>("GUANINES:RING6",0.49));
  this->ringIntensityMap.insert(std::pair<std::string,double>("THYMINES:RING6",0.28));
  this->ringIntensityMap.insert(std::pair<std::string,double>("PHENYLALANINES:RING6",1.27));
  this->ringIntensityMap.insert(std::pair<std::string,double>("TYROSINES:RING6",1.10));
  this->ringIntensityMap.insert(std::pair<std::string,double>("TRYPTOPHANS:RING6",1.27));
  this->ringIntensityMap.insert(std::pair<std::string,double>("HISTIDINES:RING6",1.40));    
  this->ringIntensityMap.insert(std::pair<std::string,double>("ADENINES:RING5",0.95));
  this->ringIntensityMap.insert(std::pair<std::string,double>("GUANINES:RING5",0.81));
  this->ringIntensityMap.insert(std::pair<std::string,double>("TRYPTOPHANS:RING5",1.02));       
}

double LARMORD::getRingIntensity(const std::string &key)
{
    if (this->ringIntensityMap.find (key) == this->ringIntensityMap.end()){
        return 0.0;
    } else {
        return (this->ringIntensityMap.at(key));
    }
}

void LARMORD::clearRings()
{
	 this->ringIntensity.clear();
	 this->ringNormal.clear();
	 this->ringCenterOfMass.clear();
}

