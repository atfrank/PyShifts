/* 

  Copyright University of Michigan.
  This file is part of the Larmor software suite and is made available under license.
  University of Michigan (UM) TECHtransfer: phone: 734-763-0614 email: techtransfer@umich.edu.
  
  Author: Aaron T. Frank
     
*/

#include "Molecule.hpp"
#include "Residue.hpp"
#include "Atom.hpp"
#include "Misc.hpp"
#include "Analyze.hpp"
#include "LARMORD.hpp"
#include "Trajectory.hpp"

#include <iostream>
#include <cstdlib>
#include <fstream>
#include <limits>
#include <sstream>
#include <math.h>
#include <algorithm>


void usage(){
  std::cerr << "====================================================" << std::endl;
  std::cerr << "LARMORD Chemical Shift Predictor v1.00" << std::endl;
  std::cerr << "(c) 2014 Aaron T. Frank and University of Michigan." << std::endl;
  std::cerr << "====================================================" << std::endl;
  std::cerr << "Usage:   larmord [-options] <PDBfile>" << std::endl;
  std::cerr << "Options: [-csfile CSfile]" << std::endl;
  std::cerr << "         [-cutoff CUToff]" << std::endl;
  std::cerr << "         [-beta BETA]" << std::endl;
  std::cerr << "         [-rna]" << std::endl;
  std::cerr << "         [-residueSelection]" << std::endl;
  std::cerr << "         [-nucleusSelection]" << std::endl;
  std::cerr << "         [-ring] [-cutoffRing (9999.0 Å)]" << std::endl;     
  std::cerr << "         [-trj TRAJfile]" << std::endl;
  std::cerr << "         [-skip frames] [-start frame] [-stop frame]" << std::endl;  
  std::cerr << "         [-identification ID]" << std::endl;
  std::cerr << std::endl;
  exit(0);
}

int main (int argc, char **argv){
  int i;
  unsigned int f;
  unsigned int ainx;
  std::stringstream resid;
  std::vector<std::string> pdbs;
  std::string currArg;
  std::string fchemshift;
  std::string nucleus;
  std::string resname;
  std::string resnameCode;
  std::string atomname;
  std::string key, residID;
  std::string identification;
  std::string residue_select;
  std::string nucleus_select;
  std::vector<int> selected_residues;
  std::vector<std::string> selected_nuclei;
  bool isResidue;
  bool isNucleus;
  bool mismatchCheck;
  std::map<std::string, double> histo;
  std::string moltype;
  
  std::vector<std::string> trajs;
  int start;
  int stop=std::numeric_limits<int>::max();
  int skip;
  bool startFlag=false;
  unsigned int itrj;
  std::ifstream trjin;
  Trajectory *ftrjin;
  unsigned int nframe;
  unsigned int process;
  unsigned int k;

  start=0;
  skip=0;
  nframe=0;
  process=0;
  
  int beta;
  double dist;
  double randcs;
  double expcs;
  double cutoff;
  
  std::vector<std::vector<double> > neighborDistances;
 	bool ringCurrent=false;
  double cutoffRing=9999.9;
  double ringc;
 
  Molecule *neighbormol;
  neighbormol=NULL;
  
  Atom *ai, *aj;  
  ai=NULL;
  aj=NULL;
  fchemshift="";
  identification="None";
  cutoff=99999.9;
  beta=-3;
  residue_select = "";
  nucleus_select = "";
  moltype = "protein";
  
 
  LARMORD *larm;
  larm=NULL;

  pdbs.clear();
  selected_residues.clear();
  selected_nuclei.clear();
  isResidue = true;
  mismatchCheck = false;
  histo.clear();

  for (i=1; i<argc; i++){
    currArg=argv[i];
    if (currArg.compare("-h") == 0 || currArg.compare("-help") == 0)
    {
      usage();
    }
    else if (currArg.compare("-csfile") == 0 || currArg.compare("-c") == 0 )
    {
      currArg=argv[++i];
      fchemshift=currArg;
    }  
    else if (currArg.compare("-mismatchCheck") == 0 )
    {
				mismatchCheck=true;
    }    
    else if (currArg.compare("-residueSelection") == 0 )
    {
      currArg=argv[++i];    
			residue_select=currArg;			
    }      
    else if (currArg.compare("-nucleusSelection") == 0 )
    {
      currArg=argv[++i];    
			nucleus_select=currArg;
    }                
    else if (currArg.compare("-identification") == 0 || currArg.compare("-i") == 0 )
    {
      currArg=argv[++i];
      identification=currArg;
    }      
    else if (currArg.compare("-cutoff") == 0)
    {
      currArg=argv[++i];
      std::stringstream(currArg) >> cutoff;
    }
    else if (currArg.compare("-cutoffRing") == 0){
      currArg=argv[++i];
      cutoffRing = atof(currArg.c_str());
    }   
    else if (currArg.compare("-ring") == 0){
      ringCurrent = true;
    }    
    else if (currArg.compare("-beta") == 0)
    {
      currArg=argv[++i];
      std::stringstream(currArg) >> beta;
    }
    else if (currArg.compare("-rna") == 0)
    {
      moltype="rna";
    }
    else if (currArg.compare("-trj") == 0 || currArg.compare("-traj") == 0)
    {
      currArg=argv[++i];
      trajs.push_back(currArg);
    }
    else if (currArg.compare("-skip") == 0)
    {
      currArg=argv[++i];
      std::stringstream(currArg) >> skip;
    }
    else if (currArg.compare("-start") == 0)
    {
      currArg=argv[++i];
      std::stringstream(currArg) >> start;
      start--;
      startFlag=true;
    }
    else if (currArg.compare("-stop") == 0)
    {
      currArg=argv[++i];
      std::stringstream(currArg) >> stop;
    }    
    else if (currArg.compare(0,1,"-") == 0)
    {
      std::cerr << "Warning: Skipping unknown option \"" << currArg << "\"" << std::endl;
    }
    else{
      pdbs.push_back(currArg);
    }
  }
  if (pdbs.size() == 0)
  {
    std::cerr << std::endl << "Error: Please provide an input file" << std::endl << std::endl;
    usage();
  }
  
  Molecule *mol=NULL;

  if (residue_select.length() != 0 || nucleus_select.length() != 0){
		Molecule *cmol;
		cmol = Molecule::readPDB(pdbs.at(0));
		
		if (residue_select.length() != 0){
			cmol->select(residue_select);
			neighbormol=cmol->copy(true);	
			Residue *res;
			selected_residues.clear();
	
			for (unsigned int j=0; j< neighbormol->getResVecSize(); j++)
			{
				res=neighbormol->getResidue(j);
				selected_residues.push_back(res->getResId());
			}
		}

		if (nucleus_select.length() != 0){
			cmol->select(nucleus_select);
			neighbormol=cmol->copy(true);	
			Atom *atm;
			selected_nuclei.clear();			
			for (unsigned int j=0; j< neighbormol->getAtmVecSize(); j++)
			{
				atm=neighbormol->getAtom(j);
				if (std::find(selected_nuclei.begin(),selected_nuclei.end(),atm->getAtmName()) == selected_nuclei.end()){
					selected_nuclei.push_back(atm->getAtmName());
				}					
			}
		}
	}

  
  if (trajs.size() > 0)
  {
    if (pdbs.size() > 1)
    {
      std::cerr << std::endl << "Warning: Only the first PDB structure is used for trajectory analysis" << std::endl << std::endl;
    }
    /* Trajectory analysis */
    /* instantiate LARMORD */
    mol=Molecule::readPDB(pdbs.at(0));
    mol->selAll();
    larm = new LARMORD(mol,fchemshift,"","","","",false,false,false,true,moltype,ringCurrent,cutoffRing);
        
    /* Process trajectories */
    for (itrj=0; itrj< trajs.size(); itrj++)
    {
      trjin.open(trajs.at(itrj).c_str(), std::ios::binary);
      if (trjin.is_open())
      {
        ftrjin=new Trajectory;
        ftrjin->setMolecule(mol);
        if (ftrjin->findFormat(trjin) == true)
        {
          ftrjin->readHeader(trjin);
          if (skip > 0 && startFlag == false)
          {
            start=skip;
          }
          /* Loop through desired frames */
          for (i=start; i< ftrjin->getNFrame() && i< stop; i=i+1+skip)
          {
            if( ftrjin->readFrame(trjin, i) == false)
            {
              std::cerr << "Warning: EOF found before the next frame could be read" << std::endl;
              break;
            }
            nframe++;
            /* get distance matrix */
            mol->assignAtmInx();
            mol->selAll();
            Analyze::pairwiseDistance(mol, neighborDistances);
            
            mol->select(":.HEAVY");
            neighbormol= mol->copy();

						if(ringCurrent == true){
							larm->clearRings();	
							larm->setUpRings(mol);	
						}
         
						for (unsigned int j=0; j< mol->getAtmVecSize(); j++)
						{
							ai = mol->getAtom(j);
							nucleus = ai->getAtmName();
							if (larm->getShiftAtom(nucleus)==true)
							{
								resname = ai->getResName();
								if(resname == "GUA")
									resnameCode = "1 0 0 0";
								if(resname == "ADE")
									resnameCode = "0 1 0 0";
								if(resname == "CYT")
									resnameCode = "0 0 1 0";
								if(resname == "URA")
									resnameCode = "0 0 0 1";

								resid << ai->getResId();
								residID = resid.str();
								key = resid.str()+":"+resname+":"+nucleus;
								resid.str("");            
								expcs = larm->getExperimentalCS(key);

								// ring current effect
								ringc = 0.0;
								if(ringCurrent == true){
									ringc = larm->ringCurrentCompute(ai->getCoor());
								}          								
								//std::cout << key << expcs << std::endl;
								isResidue = std::find(selected_residues.begin(),selected_residues.end(),ai->getResId())!= selected_residues.end();
								isNucleus = std::find(selected_nuclei.begin(),selected_nuclei.end(),ai->getAtmName()) != selected_nuclei.end();
								if( (fchemshift.length() == 0 || expcs != 0.0) && (selected_residues.size() == 0 || isResidue) && (selected_nuclei.size() == 0 || isNucleus) )
								{
									key = resname+":"+nucleus;
									randcs = larm->getRandomShift(key);
									histo.clear();        																	
									if(true)
									{
										ainx = ai->getAtmInx();
										for (unsigned int l=0; l < neighbormol->getAtmVecSize(); l++)
										{
											aj = neighbormol->getAtom(l);
											dist = neighborDistances.at(ainx).at(aj->getAtmInx());
								
											if(ai!=aj && dist < cutoff)
											{
												key = aj->getResName()+":"+aj->getAtmName();												
												if (histo.find(key) != histo.end()){
													histo.at(key)+=pow(dist,beta);
												}
												else{
													histo.insert(std::pair<std::string, double>(key, pow(dist,beta)));
												}
											}
										}
									}
									process++;
									// print header
									if(process==1)
									{
										std:: cout << "frame ID resid nucleus expCS randCS resnameG resnameA resnameC resnameU ringCurrent";
										for (k=0; k< larm->atomTypes.size(); k++){
											key=larm->atomTypes.at(k);
											std::cout << key;
											if (k != larm->atomTypes.size()-1 ){
												std::cout << "   ";
											}
										}
										std::cout << std::endl;  
									} 
									// print histogram
									std:: cout << nframe << "   " << identification << "   " << residID << "   " <<  nucleus << "   " << expcs << "   "  << randcs << "   " << resnameCode << "   " << ringc << "   "; 
									for (k=0; k< larm->atomTypes.size(); k++)
									{
										key=larm->atomTypes.at(k);
										if (histo.find(key) != histo.end()){
											std::cout << histo.at(key);
										}
										else{
											std::cout << "0";
										}
										if (j != larm->atomTypes.size()-1){
											std::cout << "   ";
										}
									}
									std::cout << std::endl;															
								}
							}        
						}
            
          }
        }
        else
        {
          std::cerr << "Warning: Skipping unknown trajectory format \"";
          std::cerr << trajs.at(itrj) << "\"" << std::endl;
        }
        if (ftrjin != NULL){
          delete ftrjin;
        }
      }
      trjin.close();
    }
  }
  else 
  { 
    /* instantiate LARMORD */
    for (f=0; f< pdbs.size(); f++)
    {  
      mol=Molecule::readPDB(pdbs.at(f));
      larm = new LARMORD(mol,fchemshift,"","","","",false,false,false,true,moltype,ringCurrent,cutoffRing);
      
      //std::cerr << "Processing file \"" << pdbs.at(f) << "..." << std::endl;
      /* get distance matrix */
      mol->assignAtmInx();
      mol->selAll();
      Analyze::pairwiseDistance(mol, neighborDistances);

      mol->select(":.HEAVY");
      neighbormol= mol->copy();
 
			if(ringCurrent == true){
				larm->clearRings();	
				larm->setUpRings(mol);	
			}
			            
      for (unsigned int j=0; j< mol->getAtmVecSize(); j++)
      {
        ai = mol->getAtom(j);
        nucleus = ai->getAtmName();
        if (larm->getShiftAtom(nucleus)==true)
        {
          resname = ai->getResName();
					if(resname == "GUA")
						resnameCode = "1 0 0 0";
					if(resname == "ADE")
						resnameCode = "0 1 0 0";
					if(resname == "CYT")
						resnameCode = "0 0 1 0";
					if(resname == "URA")
						resnameCode = "0 0 0 1";          
          resid << ai->getResId();
          residID = resid.str();
          key = resid.str()+":"+resname+":"+nucleus;
          resid.str("");            
          expcs = larm->getExperimentalCS(key);

					// ring current effect
					ringc = 0.0;
					if(ringCurrent == true){
						ringc = larm->ringCurrentCompute(ai->getCoor());
					}          
          //std::cout << key << expcs << std::endl;
          isResidue = std::find(selected_residues.begin(),selected_residues.end(),ai->getResId())!= selected_residues.end();
          isNucleus = std::find(selected_nuclei.begin(),selected_nuclei.end(),ai->getAtmName()) != selected_nuclei.end();
          if( (fchemshift.length() == 0 || expcs != 0.0) && (selected_residues.size() == 0 || isResidue) && (selected_nuclei.size() == 0 || isNucleus) )
          {
						key = resname+":"+nucleus;
						randcs = larm->getRandomShift(key);
						histo.clear();        																	
						if(true)
						{
							ainx = ai->getAtmInx();							
							for (unsigned int l=0; l < neighbormol->getAtmVecSize(); l++)
							{
								aj = neighbormol->getAtom(l);
								dist = neighborDistances.at(ainx).at(aj->getAtmInx());
								
								if(ai!=aj && dist < cutoff)
								{
									key = aj->getResName()+":"+aj->getAtmName();												
									if (histo.find(key) != histo.end()){
										histo.at(key)+=pow(dist,beta);
									}
									else{
										histo.insert(std::pair<std::string, double>(key, pow(dist,beta)));
									}
								}
							}
						}
						process++;
						// print header
						if(process==1)
						{
							std:: cout << "frame ID resid nucleus expCS randCS resnameG resnameA resnameC resnameU ringCurrent ";
							for (k=0; k< larm->atomTypes.size(); k++){
								key=larm->atomTypes.at(k);
								std::cout << key;
								if (k != larm->atomTypes.size()-1 ){
									std::cout << "   ";
								}
							}
							std::cout << std::endl;  
						} 
						// print histogram
						std:: cout << f << "   " << identification << "   " << resname << "   " << residID << "   " <<  nucleus << "   " << expcs << "   " << randcs << "   " << resnameCode <<  "   " << ringc << "   "; 
						for (k=0; k< larm->atomTypes.size(); k++)
						{
							key=larm->atomTypes.at(k);
							if (histo.find(key) != histo.end()){
								std::cout << histo.at(key);
							}
							else{
								std::cout << "0";
							}
							if (j != larm->atomTypes.size()-1){
								std::cout << "   ";
							}
						}
						std::cout << std::endl;															
          }
        }        
      }
      delete mol;
    }
  }  
  if(larm!=NULL)
    delete larm;
  return 0;
}
