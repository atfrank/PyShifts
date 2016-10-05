/* 

  Copyright University of Michigan.
  This file is part of the Larmor software suite and is made available under license.
  University of Michigan (UM) TECHtransfer: phone: 734-763-0614 email: techtransfer@umich.edu.
  
  Author: Sean M. Law and Aaron T. Frank
     
*/

#ifndef TRAJECTORY_H
#define TRAJECTORY_H

#include <string>
#include <vector>

//Forward Declaration
class Molecule;

class Trajectory {
  private:
    std::string format;
    bool swab; //Swap bytes
    Molecule *mol;
    bool show;
    bool scan;
    unsigned int iframe;
    unsigned int hdrSize;
    unsigned int frameSize;
    unsigned int lastFrame;
    std::streampos lastPos;

    std::string hdr;
    int nframe; //ICNTRL[1], Number of frames
    int npriv; //ICNTRL[2], Number of previous integration steps
    int nsavc; //ICNTRL[3], Frequency for saving of frames
    int nstep; //ICNTRL[4],Number of steps in the run that created this file
    int qvelocity; //ICNTRL[5], Velocity flag
    
    int nblock; //ICNTRL[7], Total number of blocks in lambda dynamics
    int dof; //ICNTRL[8], Degrees of freedom
    int nfixed; //ICNTRL[9], Number of fixed atoms
    float tstep; //ICNTRL[10], AKMA units
    int qcrystal; //ICNTRL[11], Crystal lattice/Periodic boundaries flag
    int q4d; //ICNTRL[12], 4D trajectory flag
    int qcharge; //ICNTRL[13], Fluctuating charges flag
    int qcheck; //ICNTRL[14], Consistency check flag, See "NOCHeck" keyword

    int version; //ICNTRL[20], CHARMM version
  
    std::vector<std::string> title;
    int natom;
    std::vector<int> fixinx;

    std::vector<double> pb;
    double pbx;
    double pby;
    double pbz;
    double pbalpha;
    double pbbeta;
    double pbgamma;
    
    std::vector<float> x;
    std::vector<float> y;
    std::vector<float> z;

    std::string endian;

    typedef union {
      unsigned int ui;
      int i;
      char c[4];
      float f;
    } binbuf;

  public:
    Trajectory ();
    bool findFormat(std::ifstream &trjin);
    template <class BinBuf> 
      BinBuf* readFortran(std::ifstream &trjin, BinBuf *buffer, int &length);
    template <class BinBuf>
      void writeFortran(std::ofstream &trjout, BinBuf *buffer, int &length);
    void clearHeader();
    void readHeader(std::ifstream &trjin);
    void writeHeader(std::ofstream &trjout);
    void showHeader();
    void cloneHeader(Trajectory *ftrjin);
    bool readFrame(std::ifstream &trjin, unsigned int frame);
    void scanFrame(std::ifstream &trjin, unsigned int frame);
    void writeFrame(std::ofstream &trjout, Trajectory *ftrjin=NULL);
    std::string getHeader(){return hdr;};
    void setMolecule(Molecule *molin);

    //Get
    std::string getFormat();
    unsigned int getHdrSize();
    unsigned int getFrameSize();
    unsigned int getLastFrame();
    std::streampos getLastPos();
    bool getShow();
    bool getScan();
    Molecule* getMolecule();
    std::string getHdr();
    int getNFrame();
    int getNPriv(); 
    int getNSavc();
    int getNStep();
    int getQVelocity();

    int getDOF(); 
    int getNFixed(); 
    float getTStepAKMA();
    double getTStepPS();
    int getQCrystal(); 
    int getQ4D(); 
    int getQCharge(); 
    int getQCheck(); 

    int getVersion(); 

    std::vector<std::string> getTitle();
    std::string getTitle(const unsigned int &element);
    int getNAtom();
    unsigned int getFixInxVecSize();
    int getFixInx(int element); 

    //Set
    void setFormat(const std::string &formatin);
    void setDefaultHeader();
    void setHdrSize(const unsigned int &size);
    void setFrameSize(const unsigned int &size);
    void setLastFrame(const unsigned int &frame);
    void setLastPos(const std::streampos &pos);
    void setShow(const bool &val);
    void setScan(const bool &val);
    void setHdr(const std::string &hdrin);
    void setNFrame(const int &nframein);
    void setNPriv(const int &nprivin);
    void setNSavc(const int &nsavcin);
    void setNStep(const int &nstepin);
    void setQVelocity(const int &qvelocityin);

    void setDOF(const int &dofin);
    void setNFixed(const int &nfixedin);
    void setTStep(const double &tstepin); //Picosecond input units
    void setTStep(const float &tstepin); //AKMA input units
    void setQCrystal(const int &qcrystalin);
    void setQ4D(const int &q4din);
    void setQCharge(const int &qchargein);
    void setQCheck(const int &qcheckin);

    void setVersion(const int &versionin);

    void setTitle(const std::vector<std::string> &titlein);
    void setTitle(const std::string &titlein, const unsigned int &element=0);
    void clearTitle();
    void addTitle(const std::string &titlein);
    void setNAtom(const int &natomin);
    void addFixInx(const int &elementin);
    void clearFixInx();
};


#endif
