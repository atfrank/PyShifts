/* 

  Copyright University of Michigan.
  This file is part of the Larmor software suite and is made available under license.
  University of Michigan (UM) TECHtransfer: phone: 734-763-0614 email: techtransfer@umich.edu.
  
  Author: Sean M. Law and Aaron T. Frank
     
*/

#ifndef MISC_H
#define MISC_H

#include <vector>
#include <string>
#include <sstream>
#include <iostream>
//#include <cctype>
//#include <cmath>
//#include <algorithm>
//#include <typeinfo>

class Misc {
  private:

  public:
    static void splitStr (const std::string &str, const std::string &delim, std::vector<std::string> &out, const bool repeat=true);
    template <class SplitVec>
      static void splitNum (const std::string &str, const std::string &delim, std::vector<SplitVec> &out, const bool repeat=true);
    static std::string replace (const std::string &str, const std::string search=" ", const std::string replace="_", const bool globalFlag=false);
    static bool isdigit (const std::string &str);
    static bool isdouble (const std::string &str);
    static bool isfloat (const std::string &str);
    static bool isalpha (const std::string &str);
    static bool isrange (const std::string &str);
    static std::string trim (const std::string &str, const std::string t=" ");
    static std::string processRange (const std::string &start, const std::string &end);
    static std::string toupper (const std::string &str);
    static std::string tolower (const std::string &str);
    static int atoi (std::string &str, const unsigned int offset=0);
    static double hypot (const double &a, const double &b);
    template <class First, class Second>
      static bool sortPairFirst(const std::pair<First, Second> &a, const std::pair<First, Second> &b);
    template <class First, class Second>
      static bool sortPairSecond(const std::pair<First, Second> &a, const std::pair<First, Second> &b);
    template <class First, class Second>
      static bool findUniqueFirst(const std::pair<First, Second> &a, const std::pair<First, Second> &b);
    template <class First, class Second>
      static bool findUniqueSecond(const std::pair<First, Second> &a, const std::pair<First, Second> &b);
    static double rmse(std::vector<double> errorVec);
    static double mae(std::vector<double> errorVec);
    static double kendall(std::vector<double> data1, std::vector<double> data2);
};

#endif
