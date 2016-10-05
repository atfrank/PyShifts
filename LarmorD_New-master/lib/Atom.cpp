//Sean M. Law
//Aaron T. Frank
  
/*
*/

#include "Atom.hpp"

#include "Constants.hpp"

#include <iostream>
#include <limits>

Atom::Atom(){
  tag="NoID";
  recname="ATOM";
  atmnum=0;
  atmname="    ";
  atmtype="  ";
  alt=" ";
  res=NULL;
  resname="   ";
  chn=NULL;
  chainid=" ";
  realid=" ";
  resid=0;
  coor=Coor(0.0, 0.0, 0.0);
  occu=0.0;
  bfac=0.0;
  elem=" ";
  segid=" ";
  sel=true;
  summary="";
  ss="";
  mass=1.0;
  charge=0.0;
  atminx=std::numeric_limits<unsigned int>::max();
  data.clear();
  bonds.clear();
}

Atom::Atom(int atmnumin, std::string atmnamein, std::string resnamein, int residin, Coor coorin, std::string segidin){
  tag="NoID";
  recname="ATOM";
  atmnum=atmnumin;
  atmname=atmnamein;
  atmname="  ";
  alt=" ";
  res=NULL;
  resname=resnamein;
  chn=NULL;
  chainid=" ";
  realid=" ";
  resid=residin;
  coor=coorin;
  occu=0.0;
  bfac=0.0;
  elem=" ";
  segid=segidin;
  sel=true;
  summary="";
  ss="";
  mass=1.0;
  charge=0.0;
  atminx=std::numeric_limits<unsigned int>::max();
  data.clear();
  bonds.clear();
}

void Atom::reset(){
  tag="NoID";
  recname="ATOM";
  atmnum=0;
  atmname="    ";
  atmtype="  ";
  alt=" ";
  res=NULL;
  resname="   ";
  chn=NULL;
  chainid=" ";
  realid=" ";
  resid=0;
  coor=Coor(0.0, 0.0, 0.0);
  occu=0.0;
  bfac=0.0;
  elem=" ";
  segid="  ";
  sel=true;
  summary="";
  ss="";
  mass=1.0;
  charge=0.0;
  data.clear();
  bonds.clear();
}

void Atom::clone(Atom* atmin){
  tag=atmin->getTag();
  recname=atmin->getRecName();
  atmnum=atmin->getAtmNum();
  atmname=atmin->getAtmName();
  atmtype=atmin->getAtmType();
  alt=atmin->getAlt();
  res=atmin->getResidue();
  resname=atmin->getResName();
  chn=atmin->getChain();
  chainid=atmin->getChainId();
  realid=atmin->getRealId();
  resid=atmin->getResId();
  icode=atmin->getICode();
  coor=atmin->getCoor();
  occu=atmin->getOccu();
  bfac=atmin->getBFac();
  elem=atmin->getElem();
  segid=atmin->getSegId();
  sel=atmin->getSel();
  summary=atmin->getSummary();
  ss=atmin->getSS();
  mass=atmin->getMass();
  charge=atmin->getCharge();
  atminx=atmin->getAtmInx();
  data=atmin->getData();
  bonds=atmin->getBonds();
}

void Atom::dummy(){
  tag="NoID";
  recname="ATOM";
  atmnum=1;
  atmname="CA";
  atmtype="C.3";
  alt=" ";
  res=NULL;
  resname="ALA";
  chn=NULL;
  chainid="+";
  realid="+";
  resid=999;
  coor=Coor(0.0, 0.0, 0.0);
  occu=0.0;
  bfac=0.0;
  elem=" ";
  segid="DUMY";
  sel=true;
  summary="";
  ss="";
  mass=1.0;
  charge=0.0;
  atminx=std::numeric_limits<unsigned int>::max();
  data.clear();
  bonds.clear();
}

//Get atom info
std::string& Atom::getTag(){
  return tag;
}

std::string& Atom::getRecName(){
  return recname;
}

int& Atom::getAtmNum(){
  return atmnum;
}

std::string& Atom::getAtmName(){
  return atmname;
}

std::string& Atom::getAtmType(){
  return atmtype;
}

std::string& Atom::getAlt(){
  return alt;
}

Residue* Atom::getResidue(){
  return res;
}

std::string& Atom::getResName(){
  return resname;
}

Chain* Atom::getChain(){
  return chn;
}

std::string& Atom::getChainId(){
  return chainid;
}

std::string& Atom::getRealId(){
  return realid;
}

int& Atom::getResId(){
  return resid;
}

std::string& Atom::getICode(){
  return icode;
}

Coor& Atom::getCoor () {
  return coor;
}

double& Atom::getX () {
  return coor.x();
}

double& Atom::getY () {
  return coor.y();
}

double& Atom::getZ () {
  return coor.z();
}

double& Atom::getOccu(){
  return occu;
}

double& Atom::getBFac(){
  return bfac;
}

std::string& Atom::getElem(){
  return elem;
}

std::string& Atom::getSegId(){
  return segid;
}

bool& Atom::getSel(){
  return sel;
}

std::string& Atom::getSummary(){
  return summary;
}

std::string& Atom::getSS(){
  return ss;
}

double& Atom::getMass(){
  return mass;
}

double& Atom::getCharge(){
  return charge;
}

unsigned int& Atom::getAtmInx(){
  return atminx;
}


std::vector<double>& Atom::getData(){
  return data;
}

double& Atom::getDataPoint(const unsigned int element){
  if (element >= data.size()){
    std::cerr << "Error: Atom::getDataPoint Out of Range" << std::endl;
  }
  return data.at(element);
}

unsigned int Atom::getDataSize(){
  return data.size();
}

unsigned int Atom::getBondsSize(){
  return bonds.size();
}

std::vector<Atom*>& Atom::getBonds(){
  return bonds;
}

Atom* Atom::getBond(const unsigned int &element){
  if (element >= bonds.size()){
    std::cerr << "Error: Atom::getBond Out of Range" << std::endl;
  }
  return bonds.at(element);
}

//Set atom info
void Atom::setTag(const std::string& tagin){
  tag=tagin;
}

void Atom::setRecName(const std::string& recnamein){
  recname=recnamein;
}

void Atom::setAtmNum(const int& atmnumin){
  atmnum=atmnumin;
}

void Atom::setAtmName(const std::string& atmnamein){
  atmname=atmnamein;
}

void Atom::setAtmName(){
  atmname="    ";
}

void Atom::setAtmType(const std::string& atmtypein){
  atmtype=atmtypein;
}

void Atom::setAtmType(){
  atmtype="  ";
}

void Atom::setAlt(const std::string& altin){
  alt=altin;
}

void Atom::setAlt(){
  alt=" ";
}

void Atom::setResidue(Residue* resin){
  res=resin;
}

void Atom::setResName(const std::string& resnamein){
  resname=resnamein;
}

void Atom::setResName(){
  resname=" ";
}

void Atom::setChain(Chain* chnin){
  chn=chnin;  
}

void Atom::setChainId(const std::string& chainidin){
  //Modified chain id if duplicated
  chainid=chainidin;
}

void Atom::setRealId(const std::string& realidin){
  //Original chain id from PDB
  realid=realidin;
}

void Atom::setChainId(){
  chainid=" ";
}

void Atom::setRealId(){
  realid=" ";
}

void Atom::setResId(const int& residin){
  resid=residin;
}

void Atom::setICode(const std::string& icodein){
  icode=icodein;
}

void Atom::setICode(){
  icode=" ";
}

void Atom::setX (const double& xin){
  coor.x()=xin;
}

void Atom::setY (const double& yin){
  coor.y()=yin;
}

void Atom::setZ (const double& zin){
  coor.z()=zin;
}

void Atom::setCoor (const Coor& coorin){
  coor=coorin;
}

void Atom::setOccu(const double& occuin){
  occu=occuin;
}

void Atom::setBFac(const double& bfacin){
  bfac=bfacin;
}

void Atom::setElem(const std::string& elemin){
  elem=elemin;
}

void Atom::setSegId(const std::string& segidin){
  segid=segidin;
}

void Atom::setSegId(){
  segid="  ";
}

void Atom::setSel(const bool selin){
  sel=selin;
}

void Atom::setSummary(const std::string& summaryin){
  summary=summaryin;
}

void Atom::setSS(const std::string& ssin){
  ss=ssin;
}

void Atom::setMass(const double& massin){
  mass=massin;
}

void Atom::setCharge(const double& chargein){
  charge=chargein;
}

void Atom::setAtmInx(const unsigned int& atminxin){
  atminx=atminxin;
}

void Atom::addData(const double& din){
  data.push_back(din);
}

void Atom::clearData(){
  data.clear();
}

void Atom::addBond(Atom* atmin){
  bonds.push_back(atmin);
}

void Atom::clearBonds(){
  bonds.clear();
}

