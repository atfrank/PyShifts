/* 

  Copyright University of Michigan.
  This file is part of the Larmor software suite and is made available under license.
  University of Michigan (UM) TECHtransfer: phone: 734-763-0614 email: techtransfer@umich.edu.
  
  Author: Sean M. Law and Aaron T. Frank
     
*/

#ifndef ANALYZE_H
#define ANALYZE_H

#include "Enum.hpp"
#include <vector>
#include <string>


//Forward Declaration
class Molecule;
class Coor;
class DTree;

//Abstract base class (cannot create instance of it!)
class Analyze {
  private:
    //Since this is an abstract base class
    //all members need to be accessed via a function or passed
    //directly to the analysis function or set the members as protected
    std::vector<std::string> sel;
    std::vector<Molecule*> mol;
    std::vector<double> tdata; //Time dependent data, maybe for averaging
    std::vector<std::vector<double> > fdata; //Frame data, cleared after each frame
    std::vector<unsigned int> modes;
    int ndata; //Total number of datapoints
    bool resel; //Re-do selection for each analysis, not implemented yet
    std::string ifile;
    std::string ofile;
    bool verbose;

  public:
    Analyze ();
    virtual ~Analyze ();
    void addSel(const std::string& selin);
    std::string getSel(const int& element);
    unsigned int getNSel();
    void addMol(Molecule* molin);
    void setMol(const int& element, Molecule* molin);
    void clearMol();
    void resizeNMol(const int sizein);
    Molecule* getMol(const int& element);
    unsigned int getNMol();
        
    //Virtual functions
    virtual void readTopology(Molecule* molin, std::string topin="");
    virtual void setupMolSel(Molecule* molin);
    virtual void preAnalysis(Molecule* molin, std::string topin=""); 
    virtual void preAnalysis();
    virtual void runAnalysis() =0; //Pure virtual function
    virtual void postAnalysis();

    //Analysis functions
    static Coor centerOfGeometry(Molecule* mol, bool selFlag=true);
    static double rmsd (Molecule* cmpmol, Molecule* refmol);
    static void rmsf (Molecule* cmpmol, Molecule* refmol, std::vector<double> &tdataIO, int &ndataIO);
    static double distance (const Coor& u, const Coor& v);
    static double distance (Molecule* sel1, Molecule* sel2, bool selFlag=true);
    static void pairwiseDistance(Molecule *mol, std::vector<std::vector<double> >& pdin);
};

//Derived classes

class AnalyzeCOG: public Analyze {
  public:
    void runAnalysis();
};

class AnalyzeRMSD: public Analyze {
  public:
    void preAnalysis(Molecule* molin, std::string topin="");
    void runAnalysis();
};


class AnalyzeDistance: public Analyze {
  public:
    void setupMolSel(Molecule* molin);
    void runAnalysis();
};

class AnalyzePairwiseDistance: public Analyze {
  public:
    void runAnalysis();
};

#endif
