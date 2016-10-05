/* 

  Copyright University of Michigan.
  This file is part of the Larmor software suite and is made available under license.
  University of Michigan (UM) TECHtransfer: phone: 734-763-0614 email: techtransfer@umich.edu.
  
  Author: Sean M. Law and Aaron T. Frank
     
*/

#ifndef LARMORD_H
#define LARMORD_H

#include <string>
#include <vector>
#include <map>


//#include "Molecule.hpp"

/* Forward Declaration  (only valid for pointers and references) */
 class Molecule;
 class Coor;

class LARMORD {
    private:
         std::map<std::string,bool> shiftAtoms;
         std::map<std::string,double> randomShifts;
         std::map<std::string, std::vector < double>  > alphas;
         std::map<std::string,std::vector < int> > betas;
         std::map<std::string,double> experimentalCS; 
         std::map<std::string,double> errorCS; 
         std::map<std::string,double> accuracy_weight;
         std::map<std::string,double> correlation_weight;         
         bool residueBasedLarmor;
         bool residueBasedWeightsLarmor;
         bool mismatchCheckLarmor; 
         std::string MolTypeLamor;
         
         std::map<std::string,double> ringIntensityMap;
         std::vector<double> ringIntensity;
         std::vector<Coor> ringNormal;
         std::vector<Coor> ringCenterOfMass;
         std::vector<std::string> aromaticsResKeys;
         double ringSeparation;
         bool ringCurrentLarmor;
         double cutoffLarmorRing;
         
    public:
        LARMORD (Molecule *mol=NULL, const std::string fchemshift="",const std::string fparmfile="",const std::string freffile="",const std::string faccfile="", const std::string fcorfile="", bool residueBased=false, bool residueBasedWeights=false, bool mismatchCheck=false, bool extractor=false, const std::string MolType="protein", bool ringCurrent=false, double cutoffRing=9999.9);
        void initializeShiftAtoms();
        void initializeRandomShifts();
        void initializeAlpha();
        void initializeAccuracyWeight();
        void initializeAtomTypes();
        bool getShiftAtom(const std::string &key);
        double getRandomShift(const std::string &key);
        std::vector<double> getAlpha (const std::string &key);
        std::vector<int> getBeta (const std::string &key);
        double getExperimentalCS(const std::string &key);
        double getErrorCS(const std::string &key);
        double getAccuracyWeight(const std::string &key);
        double getCorrelationWeight(const std::string &key);
        int getNShiftAtoms();
        void renameRes(Molecule *mol);
        void loadCSFile(const std::string fchemshift, Molecule *mol);
        void loadParmFile(const std::string fparmfile);
        void loadRefFile(const std::string freffile);
        void loadAccFile(const std::string faccfile);
        void loadCorrFile(const std::string fcorrfile);
        std::vector<std::string> atomTypes;

        /* for calculating ring current effects */
        void setUpAromaticsResKeys();
        void setUpRingIntensityMap();
        void setUpRings(Molecule *mol);
        double ellintk(const double m);
        double ellinte(const double m);
        Coor base_plane_normal(Molecule *base);
        double ringCurrentCompute(const Coor point);
        double jb_geo(const Coor point, const Coor normal);  
        double getRingIntensity(const std::string &key);
        void clearRings();                            
};
#endif
