/* 

  Copyright University of Michigan.
  This file is part of the Larmor software suite and is made available under license.
  University of Michigan (UM) TECHtransfer: phone: 734-763-0614 email: techtransfer@umich.edu.
  
  Author: Sean M. Law and Aaron T. Frank
     
*/


#include "Misc.hpp"

#include <cmath>
#include <algorithm>

void Misc::splitStr (const std::string &str, const std::string &delim, std::vector<std::string> &out, const bool repeat){
  size_t p0=0;
  size_t p1=std::string::npos;
  size_t plast=std::string::npos;
  out.clear();

  //"repeat" = true means that a blank string is added when there are
  //back-to-back delimiters. Otherwise, repeat=false ignores back-to-back delimiters.

  plast=str.find_last_not_of(delim);
  p1=str.find_first_of(delim,p0);

  while (p1 != std::string::npos){
    if (p1-p0 > 0){
      out.push_back(str.substr(p0,p1-p0));
    }
    else{
      if (repeat){
        out.push_back(str.substr(p0,p1-p0));
      }
    }
    p0=p1+1;
    p1=str.find_first_of(delim, p0);
  }
  //After last delimiter
  if (plast != std::string::npos && plast >= p0 ){
    out.push_back(str.substr(p0,p1-p0));
  }
  else{
    if (repeat){
      out.push_back(str.substr(p0,p1-p0));
    }
  }
  //std::cerr << out.size() << std::endl;
}

template <class SplitVec>
void Misc::splitNum (const std::string &str, const std::string &delim, std::vector<SplitVec> &out, const bool repeat){
  out.clear();
  out.reserve(500); //Still need resizing but reduces moving of data in memory
  size_t p0=0;
  size_t p1=std::string::npos;
  size_t plast=std::string::npos;
  SplitVec tmp;

  //Only use this if all columns are needed as a number!

  //"repeat" = true means that a blank string is added when there are
  //back-to-back delimiters. Otherwise, repeat=false ignores back-to-back delimiters.

  plast=str.find_last_not_of(delim);
  p1=str.find_first_of(delim,p0);

  while (p1 != std::string::npos){
    if (p1-p0 > 0){
      std::stringstream(str.substr(p0,p1-p0)) >> tmp;
      out.push_back(tmp);
    }
    else{
      if (repeat){
        std::stringstream(str.substr(p0,p1-p0)) >> tmp;
        out.push_back(tmp);
      }
    }
    p0=p1+1;
    p1=str.find_first_of(delim, p0);
  }
  //After last delimiter
  if (plast != std::string::npos && plast >= p0){
    std::stringstream(str.substr(p0,p1-p0)) >> tmp;
    out.push_back(tmp);
  }
  else{
    if (repeat){
      std::stringstream(str.substr(p0,p1-p0)) >> tmp;
      out.push_back(tmp);
    }
  }
  //std::cerr << out.size() << std::endl;
}

template void Misc::splitNum<int> (const std::string&, const std::string&, std::vector<int>&, const bool);

template void Misc::splitNum<unsigned int> (const std::string&, const std::string&, std::vector<unsigned int>&, const bool);

template void Misc::splitNum<double> (const std::string&, const std::string&, std::vector<double>&, const bool);

template void Misc::splitNum<float> (const std::string&, const std::string&, std::vector<float>&, const bool);

std::string Misc::replace (const std::string &str, const std::string searchStr, const std::string replaceStr, const bool globalFlag){
  std::string out;

  out=str;
  for (std::string::size_type pos=0; (pos=out.find(searchStr, pos)) != std::string::npos;){
    out.replace(pos, searchStr.length(), replaceStr);
    pos += replaceStr.length()-searchStr.length()+1;
    if (globalFlag == false){
      break;
    }
  }

  return out;
}

bool Misc::isdigit (const std::string &str){
  return str.find_first_not_of("0123456789") == std::string::npos;
}

bool Misc::isdouble (const std::string &str){
  unsigned int nDecimals=0; 
  unsigned int negPos=0;

  for (unsigned int i=0; i< str.size(); i++){
    switch (str[i]){
      case '-': negPos=i;
        continue;
      case '.': nDecimals++;
        continue;
      case '0':
      case '1':
      case '2':
      case '3':
      case '4':
      case '5':
      case '6':
      case '7':
      case '8':
      case '9':
        continue;
      default:
        return false;
    }
  }
  if (nDecimals > 1 || negPos != 0){
    return false;
  }
  return true;
}

bool Misc::isfloat (const std::string &str){
  return isdouble(str);
}

bool Misc::isalpha (const std::string &str){
  return str.find_first_not_of("ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz") == std::string::npos;
}

bool Misc::isrange (const std::string &str){

  if (str.find_first_not_of("0123456789-") != std::string::npos){
    return false;
  }
  if (str.find("-") == std::string::npos){
    return false;
  }
  if (str.find_first_not_of("-") != 0 || str.find_last_not_of("-") < str.find("-")){
    return false;
  }

  return true;
}

std::string Misc::trim (const std::string &str, const std::string t){
  std::string out=str;

  size_t pos = out.find_last_not_of(t);
  if (pos != std::string::npos){
    if (out.length() != pos+1){
      out.erase(pos+1);
    }
    pos=out.find_first_not_of(t);
    if (pos != 0){
      out.erase(0, pos);
    }
  }
  else{
    out="";
  }

  return out;
}

std::string Misc::processRange (const std::string &start, const std::string &end){
  std::stringstream ss;
  int i, istart, iend;
  
  std::stringstream(start) >> istart;
  std::stringstream(end) >> iend;

  ss << istart;
  istart++;

  for (i=istart; i<=iend; i++){
    ss << "+" << i;
  }

  return ss.str();
}

std::string Misc::toupper (const std::string &str){
  std::string out=str;
  std::transform(out.begin(), out.end(), out.begin(), ::toupper);
  return out;
}

std::string Misc::tolower (const std::string &str){
  std::string out=str;
  std::transform(out.begin(), out.end(), out.begin(), ::tolower);
  return out;
}

int Misc::atoi (std::string &str, const unsigned int offset){
  int val;
  const char *a;
  val=0;

  if (Misc::isalpha(str.substr(offset, 1))){
    a=str.substr(offset,1).c_str();
    val=static_cast<int>(a[0]-'A')+1; //A = 1, B = 2, etc
  }
  return val;
}

double Misc::hypot (const double &a, const double &b){
  if (a !=0){
    return a*sqrt(1+(b/a)*(b/a));
  }
  else if (b !=0){
    return b*sqrt(1+(a/b)*(a/b));
  }
  else{
    return 0.0;
  }
}

template <class First, class Second>
bool Misc::sortPairFirst(const std::pair<First, Second> &a, const std::pair<First, Second> &b){
  return (a.first < b.first);
}

template bool Misc::sortPairFirst<double, std::string>(const std::pair<double, std::string> &a, const std::pair<double, std::string> &b);

template <class First, class Second>
bool Misc::sortPairSecond(const std::pair<First, Second> &a, const std::pair<First, Second> &b){
  return (a.second < b.second);
}

template bool Misc::sortPairSecond<int, int>(const std::pair<int, int> &a, const std::pair<int, int> &b);

template <class First, class Second>
bool Misc::findUniqueFirst(const std::pair<First, Second> &a, const std::pair<First, Second> &b){
  return (a.first == b.first);
}

template bool Misc::findUniqueFirst<int, int>(const std::pair<int, int> &a, const std::pair<int, int> &b);


template <class First, class Second>
bool Misc::findUniqueSecond(const std::pair<First, Second> &a, const std::pair<First, Second> &b){
  return (a.second == b.second);
}

template bool Misc::findUniqueSecond(const std::pair<int, int> &a, const std::pair<int, int> &b);


double Misc::rmse (std::vector<double> errorVec){
    double error=0.0;
    if(errorVec.size()>0){    
      for (unsigned int i=0; i<errorVec.size();i++){
          error += errorVec.at(i)*errorVec.at(i);
      }
      return sqrt(error/errorVec.size());
    } else {
      return 0.0;
    }
}

double Misc::mae (std::vector<double> errorVec){
    double error=0.0;
    if(errorVec.size()>0){
      for (unsigned int i=0; i<errorVec.size();i++){
          error += fabs(errorVec.at(i));
      }
      return error/errorVec.size();
    } else {
      return 0.0;
    }
}

double Misc::kendall (std::vector<double> data1, std::vector<double> data2)
{
	int is=0,j,k,n2=0,n1=0;
	double aa,a2,a1,tau;
	int n=data1.size();
	for (j=0;j<n-1;j++) {
		for (k=j+1;k<n;k++) {
			a1=data1[j]-data1[k];
			a2=data2[j]-data2[k];
			aa=a1*a2;
			if (aa != 0.0) {
				++n1;
				++n2;
				aa > 0.0 ? ++is : --is;
			} else {
				if (a1 != 0.0) ++n1;
				if (a2 != 0.0) ++n2;
			}
		}
	}
	tau=is/(sqrt((double)n1)*sqrt((double)n2));
	return tau;
}
