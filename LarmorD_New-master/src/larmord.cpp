/* 

  Copyright University of Michigan.
  This file is part of the Larmor software suite and is made available under license.
  University of Michigan (UM) TECHtransfer: phone: 734-763-0614 email: techtransfer@umich.edu.
  
  Author: Sean M. Law and Aaron T. Frank
     
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
  std::cerr << "         [-parmfile PARAMfile]" << std::endl;
  std::cerr << "         [-reffile REFfile]" << std::endl;
  std::cerr << "         [-accfile ACCfile]" << std::endl;
  std::cerr << "         [-corrfile CORRfile]" << std::endl;  
  std::cerr << "         [-cutoff CUToff]" << std::endl;
  std::cerr << "         [-printError]" << std::endl;
  std::cerr << "         [-residueSelection]" << std::endl;
  std::cerr << "         [-nucleusSelection]" << std::endl;
  std::cerr << "         [-residueBasedWeights]" << std::endl;  
  std::cerr << "         [-mismatchCheck]" << std::endl; 
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
  unsigned int counter = 0;
  std::stringstream resid;
  std::vector<std::string> pdbs;
  std::string currArg;
  std::string fchemshift;
  std::string fparmfile;
  std::string freffile;
  std::string faccfile;
  std::string fcorfile;
  std::string nucleus;
  std::string resname;
  std::string atomname;
  std::string key;
  std::string identification;
  bool print_error;
  bool residue_based;
  bool residue_based_weights;
  std::string residue_select;
  std::string nucleus_select;
  std::vector<int> selected_residues;
  std::vector<std::string> selected_nuclei;
  bool isResidue;
  bool isNucleus;
  bool mismatchCheck;

  std::vector<std::string> trajs;
  int start;
  int stop=std::numeric_limits<int>::max();
  int skip;
  bool startFlag=false;
  unsigned int itrj;
  std::ifstream trjin;
  Trajectory *ftrjin;
  unsigned int nframe;

  start=0;
  skip=0;
  nframe=0;
  
  std::vector<double> alpha;
  std::vector<int> beta;
  double dist;
  double cspred;
  double randcs;
  double expcs;
  double cutoff;
  double error= 0.0;
  double error_mae=0.0;
  double error_rmse=0.0;
  double error_wmae=0.0;
  double error_wrmse=0.0;
  double error_flat_chi2=0.0;
  double chi2_c=2.98;
  double weight=0.0;
  double errorCS=0.0;
  double tau;
  std::vector<double> cs_C1,cs_C2,cs_CA1,cs_CA2,cs_CB1,cs_CB2,cs_N1,cs_N2,cs_HA1,cs_HA2,cs_H1,cs_H2;
  std::vector<std::vector<double> > neighborDistances;
  bool ringCurrent=false;
  double cutoffRing=9999.9;
    

  Molecule *neighbormol;
  neighbormol=NULL;
  
  Atom *ai, *aj;  
  ai=NULL;
  aj=NULL;
  fchemshift="";
  fparmfile="";
  freffile="";  
  faccfile="";  
  fcorfile="";
  identification="None";
  cutoff=99999.9;
  print_error = false;
  residue_based = false;
  residue_based_weights =false;
  residue_select = "";
  nucleus_select = "";
  cs_C1.clear();
  cs_C2.clear();
  cs_CA1.clear();
  cs_CA2.clear();
  cs_CB1.clear();
  cs_CB2.clear();
  cs_N1.clear();
  cs_N2.clear();
  cs_HA1.clear();
  cs_HA2.clear();
  cs_H1.clear();
  cs_H2.clear();

  LARMORD *larm;
  larm=NULL;

  pdbs.clear();
  selected_residues.clear();
  selected_nuclei.clear();
  isResidue = true;
  mismatchCheck = false;

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
    else if (currArg.compare("-parmfile") == 0 || currArg.compare("-p") == 0 )
    {
      currArg=argv[++i];
      fparmfile=currArg;
    }
    else if (currArg.compare("-reffile") == 0 || currArg.compare("-r") == 0 )
    {
      currArg=argv[++i];
      freffile=currArg;
    }
    else if (currArg.compare("-corrfile") == 0 || currArg.compare("-r") == 0 )
    {
      currArg=argv[++i];
      fcorfile=currArg;
    }
    else if (currArg.compare("-residueBased") == 0 )
    {
				residue_based=true;
    }    
    else if (currArg.compare("-residueBasedWeights") == 0 )
    {
				residue_based_weights=true;
    }    
    else if (currArg.compare("-mismatchCheck") == 0 )
    {
				mismatchCheck=true;
    }    
    else if (currArg.compare("-printError") == 0 )
    {
			if(fchemshift.length() > 0)
			{
					print_error=true;
			} else 
			{
					std::cerr << "Warning: -printError ignored" << std::endl;
					std::cerr << "Warning: -printError is valid only when a chemical shifts data file is suppled via --csfile " << std::endl;
					print_error=false;
			}      
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
    else if (currArg.compare("-accfile") == 0 || currArg.compare("-a") == 0 )
    {
      currArg=argv[++i];
      faccfile=currArg;
    }
    else if (currArg.compare("-identification") == 0 || currArg.compare("-i") == 0 )
    {
      currArg=argv[++i];
      identification=currArg;
    }      
    else if (currArg.compare("-trj") == 0 || currArg.compare("-traj") == 0)
    {
      currArg=argv[++i];
      trajs.push_back(currArg);
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
    larm = new LARMORD(mol,fchemshift,fparmfile,freffile,faccfile,fcorfile,residue_based,residue_based_weights,mismatchCheck);
        
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

						if (print_error)
						{
							errorCS = 0.0;
							error_mae = 0.0;
							error_rmse = 0.0;
							error_wmae = 0.0;
							error_wrmse = 0.0;
							error_flat_chi2 = 0.0;
							cs_C1.clear();
							cs_C2.clear();
							cs_CA1.clear();
							cs_CA2.clear();
							cs_CB1.clear();
							cs_CB2.clear();
							cs_N1.clear();
							cs_N2.clear();
							cs_HA1.clear();
							cs_HA2.clear();
							cs_H1.clear();
							cs_H2.clear();							
							counter = 0;
						}
            
            for (unsigned int j=0; j< mol->getAtmVecSize(); j++)
            {
              ai = mol->getAtom(j);
              nucleus = ai->getAtmName();
              if (larm->getShiftAtom(nucleus)==true)
              {
                resname = ai->getResName();
                resid << ai->getResId();
                key = resid.str()+":"+resname+":"+nucleus;
                resid.str("");
                cspred = 0.0;
                expcs = larm->getExperimentalCS(key);
                errorCS = larm->getErrorCS(key);                
                
								isResidue = std::find(selected_residues.begin(),selected_residues.end(),ai->getResId())!= selected_residues.end();
								isNucleus = std::find(selected_nuclei.begin(),selected_nuclei.end(),ai->getAtmName()) != selected_nuclei.end();
								if( (fchemshift.length() == 0 || expcs != 0.0) && (selected_residues.size() == 0 || isResidue) && (selected_nuclei.size() == 0 || isNucleus) )
								{
									key = resname+":"+nucleus;
									randcs = larm->getRandomShift(key);
									if(randcs != 0.0)
									{
										ainx = ai->getAtmInx();
										for (unsigned int l=0; l < neighbormol->getAtmVecSize(); l++)
										{
											aj = neighbormol->getAtom(l);
											if(ai!=aj){
												if (residue_based)
												{
													key = resname+":"+nucleus+":"+aj->getResName()+":"+aj->getAtmName();
												}
												else {
													key = nucleus+":"+aj->getResName()+":"+aj->getAtmName();
												}
												alpha = larm->getAlpha(key);
												beta = larm->getBeta(key);                    
												dist = neighborDistances.at(ainx).at(aj->getAtmInx());
												if (dist < cutoff){
													for (unsigned int m = 0; m < alpha.size(); m++)
													{
														cspred = cspred + alpha.at(m)*pow(dist,beta.at(m));
													}
												}                          
											}
										}
										cspred += randcs;
										if(print_error)
										{
											if (residue_based || residue_based_weights)
											{
												key = nucleus+":"+resname;
											}
											else 
											{
												key = nucleus;
											}
											weight = larm->getAccuracyWeight(key);
											error = (cspred - expcs) + errorCS;
											error_mae += fabs(error);
											error_rmse += (error*error);
											error_wmae += weight*fabs(error);
											error_wrmse += weight*weight*(error*error);
											error_flat_chi2 += (1 - exp(-(pow((weight*error/chi2_c),2))));											
											if(nucleus=="C")
											{
												cs_C1.push_back(cspred-randcs);
												cs_C2.push_back(expcs-randcs);	
											}
											else if (nucleus=="CA")
											{
												cs_CA1.push_back(cspred-randcs);
												cs_CA2.push_back(expcs-randcs);									
											}	
											else if (nucleus=="CB")
											{
												cs_CB1.push_back(cspred-randcs);
												cs_CB2.push_back(expcs-randcs);									
											}	
											else if (nucleus=="N")
											{
												cs_N1.push_back(cspred-randcs);
												cs_N2.push_back(expcs-randcs);									
											}	
											else if (nucleus=="HA")
											{
												cs_HA1.push_back(cspred-randcs);
												cs_HA2.push_back(expcs-randcs);									
											}	
											else if (nucleus=="H")
											{
												cs_H1.push_back(cspred-randcs);
												cs_H2.push_back(expcs-randcs);									
											}															
											
											counter++; 											
										} 
										else 
										{                               
											std::cout << 0 << " " << i << " " << ai->getResId() << " " << ai->getResName() << " " << nucleus << " " << cspred << " " << expcs << " " <<  randcs << " " << identification << std::endl;
										}                      
									}
                }
              }
            }
						if (print_error)
						{
							tau = 0.0; 
							if(cs_C1.size() > 0)
								tau += larm->getCorrelationWeight("C")*(1.0-Misc::kendall(cs_C1,cs_C2));
							if(cs_CA1.size() > 0)
								tau += larm->getCorrelationWeight("CA")*(1.0-Misc::kendall(cs_CA1,cs_CA2));
							if(cs_CB1.size() > 0)
								tau += larm->getCorrelationWeight("CB")*(1.0-Misc::kendall(cs_CB1,cs_CB2));
							if(cs_N1.size() > 0)
								tau += larm->getCorrelationWeight("N")*(1.0-Misc::kendall(cs_N1,cs_N2));
							if(cs_HA1.size() > 0)
								tau += larm->getCorrelationWeight("HA")*(1.0-Misc::kendall(cs_HA1,cs_HA2));
							if(cs_H1.size() > 0)
								tau += larm->getCorrelationWeight("H")*(1.0-Misc::kendall(cs_H1,cs_H2));
							std::cout << 0 << " " << i << " " << error_mae/counter << " " << sqrt(error_rmse/counter) << " " << error_wmae/counter << " " <<  sqrt(error_wrmse/counter)<< " " <<  chi2_c*chi2_c*(error_flat_chi2/counter) << " " << tau  << " " << identification << std::endl;
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
      larm = new LARMORD(mol,fchemshift,fparmfile,freffile,faccfile,fcorfile,residue_based,residue_based_weights,mismatchCheck);
      
      //std::cerr << "Processing file \"" << pdbs.at(f) << "..." << std::endl;
      /* get distance matrix */
      mol->assignAtmInx();
      mol->selAll();
      Analyze::pairwiseDistance(mol, neighborDistances);

      mol->select(":.HEAVY");
      neighbormol= mol->copy();
      
      if (print_error)
      {
				errorCS = 0.0;
				error_mae = 0.0;
				error_rmse = 0.0;
				error_wmae = 0.0;
				error_wrmse = 0.0;
				error_flat_chi2 = 0.0;
				cs_C1.clear();
				cs_C2.clear();
				cs_CA1.clear();
				cs_CA2.clear();
				cs_CB1.clear();
				cs_CB2.clear();
				cs_N1.clear();
				cs_N2.clear();
				cs_HA1.clear();
				cs_HA2.clear();
				cs_H1.clear();
				cs_H2.clear();
				
				counter = 0;
      }
      
      for (unsigned int j=0; j< mol->getAtmVecSize(); j++)
      {
        ai = mol->getAtom(j);
        nucleus = ai->getAtmName();
        if (larm->getShiftAtom(nucleus)==true)
        {
          resname = ai->getResName();
          resid << ai->getResId();
          key = resid.str()+":"+resname+":"+nucleus;  
          resid.str("");        
          cspred = 0.0;
          expcs = larm->getExperimentalCS(key);
          errorCS = larm->getErrorCS(key);
          isResidue = std::find(selected_residues.begin(),selected_residues.end(),ai->getResId())!= selected_residues.end();
          isNucleus = std::find(selected_nuclei.begin(),selected_nuclei.end(),ai->getAtmName()) != selected_nuclei.end();
          if( (fchemshift.length() == 0 || expcs != 0.0) && (selected_residues.size() == 0 || isResidue) && (selected_nuclei.size() == 0 || isNucleus) )
          {
						key = resname+":"+nucleus;
						randcs = larm->getRandomShift(key);
						//std::cout << key << " " << randcs << std::endl;
						if(randcs != 0.0)
						{
							ainx = ai->getAtmInx();
							for (unsigned int l=0; l < neighbormol->getAtmVecSize(); l++){
								aj = neighbormol->getAtom(l);
								if(ai!=aj){
									if (residue_based)
									{
										key = resname+":"+nucleus+":"+aj->getResName()+":"+aj->getAtmName();
									}
									else {
										key = nucleus+":"+aj->getResName()+":"+aj->getAtmName();
									}
									alpha = larm->getAlpha(key);
									beta = larm->getBeta(key);									
									dist = neighborDistances.at(ainx).at(aj->getAtmInx());
									if (dist < cutoff)
									{
										for (unsigned int m = 0; m < alpha.size(); m++)
										{
											//std::cout << key << " " << alpha.at(m) << " " << beta.at(m) << std::endl;                    
											//std::cout << nucleus+":"+aj->getResName()+":"+aj->getAtmName() << " " << alpha.at(m) << " " << beta.at(m) << " " << randcs << std::endl;
											cspred = cspred + alpha.at(m)*pow(dist,beta.at(m));
										}
									}
								}								
							}
							cspred += randcs;							
							if(print_error)
							{
								if (residue_based || residue_based_weights)
								{
									key = nucleus+":"+resname;
								}
								else 
								{
									key = nucleus;
								}
								weight = larm->getAccuracyWeight(key);							
								error = (cspred - expcs ) + errorCS;
								error_mae += fabs(error);
								error_rmse += (error*error);
								error_wmae += weight*fabs(error);
								error_wrmse += weight*weight*(error*error);
								error_flat_chi2 += (1 - exp(-(pow((weight*error/chi2_c),2))));
								if(nucleus=="C")
								{
									cs_C1.push_back(cspred-randcs);
									cs_C2.push_back(expcs-randcs);	
								}
								else if (nucleus=="CA")
								{
									cs_CA1.push_back(cspred-randcs);
									cs_CA2.push_back(expcs-randcs);									
								}	
								else if (nucleus=="CB")
								{
									cs_CB1.push_back(cspred-randcs);
									cs_CB2.push_back(expcs-randcs);									
								}	
								else if (nucleus=="N")
								{
									cs_N1.push_back(cspred-randcs);
									cs_N2.push_back(expcs-randcs);									
								}	
								else if (nucleus=="HA")
								{
									cs_HA1.push_back(cspred-randcs);
									cs_HA2.push_back(expcs-randcs);									
								}	
								else if (nucleus=="H")
								{
									cs_H1.push_back(cspred-randcs);
									cs_H2.push_back(expcs-randcs);									
								}															
								counter++; 
							} 
							else 
							{                               
								std::cout << 0 << " " << f+1 << " " << ai->getResId() << " " << ai->getResName() << " " << nucleus << " " << cspred << " " << expcs << " " <<  randcs << " " << identification << std::endl;
							}
						}
          }
        }        
      }
      if (print_error)
      {
				tau = 0.0;
				if(cs_C1.size() > 0)
					tau += larm->getCorrelationWeight("C")*(1.0-Misc::kendall(cs_C1,cs_C2));
				if(cs_CA1.size() > 0)
					tau += larm->getCorrelationWeight("CA")*(1.0-Misc::kendall(cs_CA1,cs_CA2));
				if(cs_CB1.size() > 0)
					tau += larm->getCorrelationWeight("CB")*(1.0-Misc::kendall(cs_CB1,cs_CB2));
				if(cs_N1.size() > 0)
					tau += larm->getCorrelationWeight("N")*(1.0-Misc::kendall(cs_N1,cs_N2));
				if(cs_HA1.size() > 0)
					tau += larm->getCorrelationWeight("HA")*(1.0-Misc::kendall(cs_HA1,cs_HA2));
				if(cs_H1.size() > 0)
					tau += larm->getCorrelationWeight("H")*(1.0-Misc::kendall(cs_H1,cs_H2));
      	std::cout << 0 << " " << f+1 << " " << error_mae/counter << " " << sqrt(error_rmse/counter) << " " << error_wmae/counter << " " <<  sqrt(error_wrmse/counter) << " " <<  chi2_c*chi2_c*(error_flat_chi2/counter) << " " << tau << " "  << identification << std::endl;
      }
      delete mol;
    }
  }  
  if(larm!=NULL)
    delete larm;
  return 0;
}
