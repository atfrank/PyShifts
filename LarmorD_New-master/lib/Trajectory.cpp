/* 

  Copyright University of Michigan.
  This file is part of the Larmor software suite and is made available under license.
  University of Michigan (UM) TECHtransfer: phone: 734-763-0614 email: techtransfer@umich.edu.
  
  Author: Sean M. Law and Aaron T. Frank
     
*/

#include "Trajectory.hpp"

#include "Constants.hpp"
#include "Molecule.hpp"
#include "Atom.hpp"

#include <fstream>
#include <iomanip>
#include <iostream>
#include <cstdlib>
#include <limits>

//Deal with Endianess

Trajectory::Trajectory (){
  swab=false;
  mol=NULL;
  show=false;
  scan=false;
  iframe=0;
  hdrSize=0;
  frameSize=0;
  lastFrame=0;
  lastPos=0;
  hdr.clear();
  nframe=0;
  npriv=0;
  nsavc=0;
  nstep=0;
  qvelocity=false;
  dof=0;
  nfixed=0;
  tstep=0.0;
  qcrystal=false;
  q4d=false;
  qcharge=false;
  qcheck=false;
  version=0;
  title.clear();
  natom=0;
  fixinx.clear();
}

void Trajectory::setDefaultHeader(){
  this->setFormat("CHARMM");
  this->setHdr("CORD");
  this->setNFrame(0);
  this->setNPriv(1);
  this->setNSavc(1);
  this->setNStep(1);
  this->setQVelocity(false);
  this->setDOF(0);
  this->setNFixed(0);
  this->setTStep(1.0);
  this->setQCrystal(false);
  this->setQ4D(false);
  this->setQCharge(false);
  this->setQCheck(false);
  this->setVersion(0);
  title.clear();
  this->setNAtom(0);
  fixinx.clear();
}

bool Trajectory::findFormat(std::ifstream &trjin){
  binbuf *buffer;
  int length;

  buffer=NULL;

  trjin.seekg(0, std::ios::beg);

  buffer=readFortran(trjin, buffer, length);

  hdr.assign(buffer[0].c,4);

  if (length == 84 && hdr.compare("CORD") == 0){
    format="CHARMM";
    swab=false;
  }

  if (buffer != NULL){
    delete buffer;
  }

  trjin.seekg(0, std::ios::beg);

  clearHeader();

  if (format.length() > 0){
    return true;
  }
  else{
    return false;
  }
}

template <class BinBuf> 
BinBuf* Trajectory::readFortran(std::ifstream &trjin, BinBuf *buffer, int &length){
  int recStart;
  int recEnd;
  BinBuf *binOut;

  binOut=NULL;

  //Read Fortran record lengths and buffer
  trjin.read(reinterpret_cast<char*>(&recStart), sizeof(int));
  binOut = new BinBuf [recStart];
  trjin.read(reinterpret_cast<char*>(binOut), recStart);
  trjin.read(reinterpret_cast<char*>(&recEnd), sizeof(int));

  //Check Fortran record length mismatch
  if (recStart == recEnd){
    length=recStart;
  }
  else{
    std::cerr << "Error: Fortran record length marker mismatch" << std::endl;
    std::exit(0);
  }

  return binOut;
}
template Trajectory::binbuf* Trajectory::readFortran<Trajectory::binbuf> (std::ifstream&, Trajectory::binbuf*, int&);
template double* Trajectory::readFortran<double> (std::ifstream&, double*, int&);
template float* Trajectory::readFortran<float> (std::ifstream&, float*, int&);
template char* Trajectory::readFortran<char> (std::ifstream&, char*, int&);

template <class BinBuf>
void Trajectory::writeFortran(std::ofstream &trjout, BinBuf *buffer, int &length){
  
  //Write Fortran record lengths and buffer
  trjout.write(reinterpret_cast<char*>(&length), sizeof(int));
  trjout.write(reinterpret_cast<char*>(buffer), length);
  trjout.write(reinterpret_cast<char*>(&length), sizeof(int));
}
template void Trajectory::writeFortran<Trajectory::binbuf> (std::ofstream&, Trajectory::binbuf*, int&);
template void Trajectory::writeFortran<double> (std::ofstream&, double*, int&);
template void Trajectory::writeFortran<float> (std::ofstream&, float*, int&);
template void Trajectory::writeFortran<char> (std::ofstream&, char*, int&);

void Trajectory::clearHeader(){
  hdr.clear();
  nframe=0;
  npriv=0;
  nsavc=0;
  nstep=0;
  qvelocity=0;
  dof=0;
  nfixed=0;
  tstep=0.0;
  qcrystal=0;
  q4d=0;
  qcharge=0;
  qcheck=0;
  version=0;
  title.clear();
  natom=0;
  fixinx.clear();
}

void Trajectory::readHeader(std::ifstream &trjin){
  binbuf *buffer;
  char *cbuffer;
  int length;
  int i;

  buffer=NULL;
  cbuffer=NULL;

  trjin.seekg(0, std::ios::end);
  this->setLastPos(trjin.tellg());

  trjin.seekg(0, std::ios::beg);

  if (format.compare("CHARMM") == 0){
    buffer=readFortran(trjin, buffer, length);
    this->setHdrSize(this->getHdrSize()+length+sizeof(int)*2);

    //HDR
    hdr.assign(buffer[0].c,4);

    //ICNTRL 1-20
    nframe=buffer[1].i;
    npriv=buffer[2].i;
    nsavc=buffer[3].i;
    nstep=buffer[4].i;
    qvelocity=buffer[5].i;
    if (qvelocity > 0){
      std::cerr << "Error: Velocity reading has yet to be implemented" << std::endl;
    }
  
    dof=buffer[8].i;
    nfixed=buffer[9].i;
  
    tstep=buffer[10].f;
    
    qcrystal=buffer[11].i;
    q4d=buffer[12].i;
    qcharge=buffer[13].i;
    qcheck=buffer[14].i;
    
    version=buffer[20].i;

    if (buffer != NULL){
      delete buffer;
    }

    //Title
    cbuffer=readFortran(trjin, cbuffer, length);
    this->setHdrSize(this->getHdrSize()+length+sizeof(int)*2);
    int ntitle=cbuffer[0];
    title.resize(ntitle);
    for (i=0; i< ntitle; i++){
      title[i].assign(cbuffer+sizeof(int)+i*80,80);
    }

    if (cbuffer != NULL){
      delete cbuffer;
    }
    

    //NATOM
    buffer=readFortran(trjin, buffer, length);
    this->setHdrSize(this->getHdrSize()+length+sizeof(int)*2);
    natom=buffer[0].i;
    if (mol != NULL && static_cast<int>(mol->getNAtom()) != natom){
      std::cerr << "Error: Atom number mismatch!" << std::endl;
    }

    if (buffer != NULL){
      delete buffer;
    }

    //FIXED
    if (nfixed > 0){
      buffer=readFortran(trjin, buffer, length);
      this->setHdrSize(this->getHdrSize()+length+sizeof(int)*2);
      std::cerr << "Warning: Fixed atoms has yet to be implemented" << std::endl;
      //for (i=0; i< length; i++){
      //  std::cerr << buffer[i].i << ":";
      //  fixinx.push_back(buffer[i].i);
      //}

      if (buffer != NULL){
        delete buffer;
      }
    }

    if (this->getShow() == true){
      this->showHeader();
    } 
  }
  else if (format.compare("AMBER") == 0){
    
  }
  else{
    //Do nothing
  }

}

void Trajectory::writeHeader(std::ofstream &trjout){
  if (format.compare("CHARMM") == 0){
    trjout.seekp(0, std::ios::beg);

    binbuf icntrl[21];
    binbuf *buffer=reinterpret_cast<binbuf *>(&icntrl[0]);
    int *ibuffer;
    int length;
    unsigned int i;

    for (i=0; i< 21; i++){
      buffer[i].i=0;
    }

    buffer[0].c[0]='C';
    buffer[0].c[1]='O';
    buffer[0].c[2]='R';
    buffer[0].c[3]='D';

    buffer[1].i=getNFrame();
    buffer[2].i=getNPriv();
    buffer[3].i=getNSavc();
    buffer[4].i=getNStep();
    icntrl[5].i=getQVelocity();

    buffer[8].i=getDOF();
    buffer[9].i=getNFixed();
    buffer[10].f=getTStepAKMA();
    buffer[11].i=getQCrystal();
    buffer[12].i=getQ4D();
    buffer[13].i=getQCharge();
    buffer[14].i=getQCheck();

    buffer[20].i=getVersion();

    length=sizeof(int)*21;

    writeFortran(trjout, buffer, length);

   
    length=sizeof(int)+title.size()*80;
    char *otitle=new char [length];
    unsigned int *ntitle=reinterpret_cast<unsigned int *>(otitle);
    *ntitle=title.size();
    for (i=0; i< title.size(); i++){
      //We need to use a temporary/intermediate string because
      //getTitle(i) returns a copy of the title. Multiple calls
      //to the same element returns two separate copies with different
      //begin() and end()
      std::string tmpTitle=getTitle(i); 
      std::copy(tmpTitle.begin(), tmpTitle.end() , otitle+sizeof(int)+i*80);
    }

    writeFortran(trjout, otitle, length);
    delete otitle;

    //binbuf natm=natom;
    if (mol != NULL){
      natom=mol->getNAtom();
    }
    ibuffer=&natom;
    length=sizeof(int);
    writeFortran(trjout, ibuffer, length);

    if (nfixed > 0){
      std::cerr << "Warning: Fixed atoms has yet to be implemented" << std::endl;
      //for (i=0; i< length; i++){
      //  std::cerr << buffer[i].i << ":";
      //  fixinx.push_back(buffer[i].i);
      //}
    }

  }
}

void Trajectory::showHeader(){
  for (unsigned int i=0; i< title.size(); i++){
    std::cout << title[i] << std::endl;
  }
  std::cout << std::fixed;
  std::cout << std::setw(25) << std::left << "Atoms" << ": " << natom << std::endl;
  if (nframe == 0){
    std::cout << std::setw(25) << std::left << "Frames" << ": " << nframe << " ***Warning: Zero frames detected!***" << std::endl;
  }
  else{
    std::cout << std::setw(25) << std::left << "Frames" << ": " << nframe << std::endl;
  }
  std::cout << std::setw(25) << std::left << "Start Frame" << ": " << npriv << std::endl;
  std::cout << std::setw(25) << std::left << "Save Frequency" << ": " << nsavc << std::endl;
  std::cout << std::setw(25) << std::left << "Dynamics Steps" << ": " << nstep << std::endl;
  std::cout << std::setw(25) << std::left << "Degrees of Freedom" << ": " << dof << std::endl;
  std::cout << std::setw(25) << std::left << "Number of Fixed" << ": " << nfixed << std::endl;
  std::cout << std::setw(25) << std::left << "Time Step (ps)" << ": " << this->getTStepPS() << std::endl;
  std::cout << std::setw(25) << std::left << "Start Time (ps)" << ": " << npriv*this->getTStepPS()/nsavc << std::endl;
  std::cout << std::setw(25) << std::left << "Periodic Boundaries" << ": " << qcrystal << std::endl;
  std::cout << std::setw(25) << std::left << "4D Trajectory" << ": " << q4d << std::endl;
  std::cout << std::setw(25) << std::left << "Fluctuating Charges" << ": " << qcharge << std::endl;
  std::cout << std::setw(25) << std::left << "Consistency Check" << ": " << qcheck << std::endl;
  std:: cout << std::setw(25) << std::left << "Version" << ": " << version << std::endl;
}

void Trajectory::cloneHeader(Trajectory *ftrjin){
  format=ftrjin->getFormat();
  hdr=ftrjin->getHdr();
  nframe=ftrjin->getNFrame();
  npriv=ftrjin->getNPriv();
  nsavc=ftrjin->getNSavc();
  nstep=ftrjin->getNStep();
  qvelocity=ftrjin->getQVelocity();

  dof=ftrjin->getDOF();
  nfixed=ftrjin->getNFixed();
  tstep=ftrjin->getTStepAKMA();
  qcrystal=ftrjin->getQCrystal();
  q4d=ftrjin->getQ4D();
  qcharge=ftrjin->getQCharge();
  qcheck=ftrjin->getQCheck();

  version=ftrjin->getVersion();

  title=ftrjin->getTitle();
  natom=ftrjin->getNAtom();
  for (unsigned int i=0; i< ftrjin->getFixInxVecSize(); i++){
    fixinx.push_back(ftrjin->getFixInx(i));
  }
}

bool Trajectory::readFrame(std::ifstream &trjin, unsigned int frame){
  double *dbuffer; //Needed for crystal!
  float *xbuffer;
  float *ybuffer;
  float *zbuffer;
  int length;
  int i;
  unsigned int j;
  unsigned int maxSeekFrame;

  if (this->getFrameSize() == 0){
    this->setFrameSize((this->getNAtom()*sizeof(float)+sizeof(int)*2)*3);
    if (qcrystal > 0){
      this->setFrameSize(this->getFrameSize()+sizeof(double)*6+2*sizeof(int));
    }
  }

  maxSeekFrame=(std::numeric_limits<int>::max()/this->getFrameSize())-10; //10 buffer

  //Assume header has been read
  //Seek frames
  if (frame-this->getLastFrame() > maxSeekFrame){
    //Seek frame by frame to be safe
    for (j=0; j< frame-this->getLastFrame(); j++){
      trjin.seekg(this->getFrameSize(), std::ios::cur);
    }
  }
  else{
    trjin.seekg((frame-this->getLastFrame())*this->getFrameSize(), std::ios::cur);
  }
  this->setLastFrame(frame+1); 
  //The "+1" positions the pointer at the end of this frame after it has been read


  if (trjin.tellg() == this->getLastPos()){
    //End-of-file, note that trjin.good()/trjin.eof() wouldn't work!
    return false;
  }

  this->iframe++;
  if (this->getShow() == true || this->getScan() == true){
    std::cout << std::fixed;
    std::cout << "Frame : " << this->iframe << " ( " << frame+1 << " / " << getNFrame() << " ) " << std::endl; 
  }

  dbuffer=NULL;
  xbuffer=NULL;
  ybuffer=NULL;
  zbuffer=NULL;

  if (qcrystal > 0){
    dbuffer=readFortran(trjin, dbuffer, length);

    if (this->getScan() == false){ //
      pb.resize(6);
      for (j=0; j< pb.size(); j++){
        pb.at(j)=dbuffer[j];
      }

      //Box dimensions
      pbx=sqrt(pb.at(0)*pb.at(0)+pb.at(1)*pb.at(1)+pb.at(3)*pb.at(3));
      pby=sqrt(pb.at(1)*pb.at(1)+pb.at(2)*pb.at(2)+pb.at(4)*pb.at(4));
      pbz=sqrt(pb.at(3)*pb.at(3)+pb.at(4)*pb.at(4)+pb.at(5)*pb.at(5));

      //Box angles
      pbalpha=acos(pb.at(4)*(pb.at(2)+pb.at(5))+pb.at(1)*pb.at(3)/(pby*pbz))*180.0/PI;
      pbbeta=acos(pb.at(3)*(pb.at(0)+pb.at(5))+pb.at(1)*pb.at(4)/(pbz*pbx))*180.0/PI;
      pbgamma=acos(pb.at(1)*(pb.at(0)+pb.at(2))+pb.at(3)*pb.at(4)/(pbx*pby))*180.0/PI;

      if (this->getShow() == true){ 
        std::cout << "Periodic Box :  ";
        std::cout << pbx << " x " << pby << " x " << pbz << " ";
        std::cout << "( " << pbalpha << " x " << pbbeta << " x " << pbgamma << " )" << std::endl;
      }
    }

    if (dbuffer != NULL){
      delete dbuffer;
    }
  }

  //Coordinates
  xbuffer=readFortran(trjin, xbuffer, length);
  ybuffer=readFortran(trjin, ybuffer, length);
  zbuffer=readFortran(trjin, zbuffer, length);

  if (this->getScan() == false){
    if (x.size() == 0 || y.size() == 0 || z.size() == 0){
      x.resize(natom);
      y.resize(natom);
      z.resize(natom);
    }
    if (natom != static_cast<int>(x.size()) || 
      natom != static_cast<int>(y.size()) || 
      natom != static_cast<int>(z.size())){
        std::cerr << "Error: The number of atoms \"natom\" does not match the number of atoms in each frame" << std::endl;
    }

    if (nfixed > 0){
      for (i=0; i< nfixed; i++){
        //Do something
      }
      std::cerr << "Warning: Fixed atoms has yet to be implemented"<< std::endl;
    }
    else{
      if (mol != NULL && length/sizeof(float) != mol->getNAtom()){
        std::cerr << "Error: The natom mismatch!" << std::endl;
      }
      for (i=0; i< natom; i++){
        if (mol == NULL){
          x.at(i)=xbuffer[i];
          y.at(i)=ybuffer[i];
          z.at(i)=zbuffer[i];
        }
        else{
          mol->getAtom(i)->setCoor(Coor(static_cast<double>(xbuffer[i]), static_cast<double>(ybuffer[i]), static_cast<double>(zbuffer[i])));
        }
  
        if (this->getShow() == true){
          std::cout << std::fixed;
          std::cout << "coor" << std::setw(10) << std::right << i+1;
          std::cout << std::setw(14) << xbuffer[i] << " ";
          std::cout << std::setw(14) << ybuffer[i] << " ";
          std::cout << std::setw(14) << zbuffer[i] << std::endl;
        } 
      }
    }
  }

  if (xbuffer != NULL){
    delete xbuffer;
  }
  if (ybuffer != NULL){
    delete ybuffer;
  }
  if (zbuffer != NULL){
    delete zbuffer;
  }

  return true; 
}

void Trajectory::writeFrame(std::ofstream &trjout, Trajectory *ftrjin){
  double *dbuffer;
  float *xbuffer;
  float *ybuffer;
  float *zbuffer;
  unsigned int i;
  int length;
  Molecule *mol;

  //crystal
  if (this->getQCrystal() > 0){
    double pbc[ftrjin->pb.size()];
    for(i=0; i< ftrjin->pb.size(); i++){
      pbc[i]=ftrjin->pb.at(i);
    }
    dbuffer=&pbc[0];
    length=sizeof(double)*ftrjin->pb.size();
    writeFortran(trjout, dbuffer, length);
  }

  if (this->getMolecule() == NULL){
    float x[ftrjin->x.size()];
    float y[ftrjin->y.size()];
    float z[ftrjin->z.size()];

    for(i=0; i< ftrjin->x.size(); i++){
      x[i]=ftrjin->x.at(i);
      y[i]=ftrjin->y.at(i);
      z[i]=ftrjin->z.at(i);
    } 

    length=sizeof(float)*ftrjin->x.size();

    xbuffer=&x[0];
    ybuffer=&y[0];
    zbuffer=&z[0];

    writeFortran(trjout, xbuffer, length);
    writeFortran(trjout, ybuffer, length);
    writeFortran(trjout, zbuffer, length);
  
  }
  else{
    mol=this->getMolecule();
    float x[mol->getNAtom()];
    float y[mol->getNAtom()];
    float z[mol->getNAtom()];

    for(i=0; i< mol->getNAtom(); i++){
      x[i]=mol->getAtom(i)->getX();
      y[i]=mol->getAtom(i)->getY();
      z[i]=mol->getAtom(i)->getZ();
    }

    length=sizeof(float)*mol->getNAtom();

    xbuffer=&x[0];
    ybuffer=&y[0];
    zbuffer=&z[0];

    writeFortran(trjout, xbuffer, length);
    writeFortran(trjout, ybuffer, length);
    writeFortran(trjout, zbuffer, length);

  }
}

void Trajectory::setMolecule(Molecule *molin){
  mol=molin;
}

unsigned int Trajectory::getHdrSize(){
  return hdrSize;
}

unsigned int Trajectory::getFrameSize(){
  return frameSize;
}

unsigned int Trajectory::getLastFrame(){
  return lastFrame;
}

std::streampos Trajectory::getLastPos(){
  return lastPos;
}

std::string Trajectory::getFormat(){
  return format;
}

bool Trajectory::getShow(){
  return show;
}

bool Trajectory::getScan(){
  return scan;
}

Molecule* Trajectory::getMolecule(){
  return mol;
}


std::string Trajectory::getHdr(){
  return hdr;
}

int Trajectory::getNFrame(){
  return nframe;
}

int Trajectory::getNPriv(){
  return npriv;
}

int Trajectory::getNSavc(){
  return nsavc;
}

int Trajectory::getNStep(){
  return nstep;
}

int Trajectory::getQVelocity(){
  return qvelocity;
}

int Trajectory::getDOF(){
  return dof;
}
int Trajectory::getNFixed(){
  return nfixed;
}

float Trajectory::getTStepAKMA(){
  return tstep;
}

double Trajectory::getTStepPS(){
  return  static_cast<double>(static_cast<int>(static_cast<double>(tstep)*AKMATPS*nsavc*100000.0+0.5))/100000.0;
}

int Trajectory::getQCrystal(){
  return qcrystal;
}

int Trajectory::getQ4D(){
  return q4d;
}

int Trajectory::getQCharge(){
  return qcharge;
}

int Trajectory::getQCheck(){
  return qcheck;
}

int Trajectory::getVersion(){
  return version;
}

std::vector<std::string> Trajectory::getTitle(){
  for (unsigned int i=0; i< title.size(); i++){
    if (title[i].find("*") != 0){
      title[i].insert(0,"*");
    }
    title[i].resize(80, ' ');
  }
  return title;
}

std::string Trajectory::getTitle(const unsigned int &element){
  return title.at(element);
}

int Trajectory::getNAtom(){
  return natom;
}

unsigned int Trajectory::getFixInxVecSize(){
  return fixinx.size();
}

int Trajectory::getFixInx(int element){
  return fixinx.at(element);
}

//Set

void Trajectory::setFormat(const std::string &formatin){
  format=formatin;
}

void Trajectory::setHdrSize(const unsigned int &size){
  hdrSize=size;
}

void Trajectory::setFrameSize(const unsigned int &size){
  frameSize=size;
}

void Trajectory::setLastFrame(const unsigned int &frame){
  lastFrame=frame;
}

void Trajectory::setLastPos(const std::streampos &pos){
  lastPos=pos;
}

void Trajectory::setShow(const bool &val){
  show=val;
}

void Trajectory::setScan(const bool &val){
  scan=val;
}

void Trajectory::setHdr(const std::string &hdrin){
  hdr=hdrin;
}
void Trajectory::setNFrame(const int &nframein){
  nframe=nframein;
}

void Trajectory::setNPriv(const int &nprivin){
  npriv=nprivin;
}

void Trajectory::setNSavc(const int &nsavcin){
  nsavc=nsavcin;
}

void Trajectory::setNStep(const int &nstepin){
  nstep=nstepin;
}

void Trajectory::setQVelocity(const int &qvelocityin){
  qvelocity=qvelocityin;
}

void Trajectory::setDOF(const int &dofin){
  dof=dofin;
}

void Trajectory::setNFixed(const int &nfixedin){
  nfixed=nfixedin;
}

void Trajectory::setTStep(const double &tstepin){
  //Convert from picoseconds double to AKMA float
  tstep=static_cast<float>((static_cast<double>(static_cast<int>(tstepin*100000.0))-0.5)/100000.0/nsavc/AKMATPS);
}

void Trajectory::setTStep(const float &tstepin){
  //Float input is already in AKMA units
  tstep=tstepin;
}

void Trajectory::setQCrystal(const int &qcrystalin){
  qcrystal=qcrystalin;
}

void Trajectory::setQ4D(const int &q4din){
  q4d=q4din;
}

void Trajectory::setQCharge(const int &qchargein){
  qcharge=qchargein;
}

void Trajectory::setQCheck(const int &qcheckin){
  qcheck=qcheckin;
}

void Trajectory::setVersion(const int &versionin){
  version=versionin;
}

void Trajectory::setTitle(const std::vector<std::string> &titlein){
  title=titlein;
  if (title.size() == 1){
    title.resize(2);
  }
  for (unsigned int i=0; i< title.size(); i++){
    if (title[i].find("*") != 0){
      title[i].insert(0,"*");
    }
    title[i].resize(80, ' ');
  }
}

void Trajectory::setTitle(const std::string &titlein, const unsigned int &element){
  if (title.size() < element){
    title.resize(element);
  }
  title[element]=titlein;
  if (title.size() == 1){
    title.resize(2);
  }
  for (unsigned int i=0; i< title.size(); i++){
    if (title[i].find("*") != 0){
      title[i].insert(0,"*");
    }
    title[i].resize(80, ' ');
  }
}

void Trajectory::clearTitle(){
  title.clear();
}

void Trajectory::addTitle(const std::string &titlein){
  title.push_back(titlein);
  for (unsigned int i=0; i< title.size(); i++){
    if (title[i].find("*") != 0){
      title[i].insert(0,"*");
    }
    title[i].resize(80, ' ');
  }
}

void Trajectory::setNAtom(const int &natomin){
  natom=natomin;
}
void Trajectory::addFixInx(const int &elementin){
  fixinx.push_back(elementin);
}

void Trajectory::clearFixInx(){
  fixinx.clear();
}

