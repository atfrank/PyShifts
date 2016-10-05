/* 

  Copyright University of Michigan.
  This file is part of the Larmor software suite and is made available under license.
  University of Michigan (UM) TECHtransfer: phone: 734-763-0614 email: techtransfer@umich.edu.
  
  Author: Sean M. Law and Aaron T. Frank
     
*/

#ifndef VECTOR_H
#define VECTOR_H

//#include <cmath>
//#include <iostream>

class Coor {
  private:
    double xcoor;
    double ycoor;
    double zcoor;

  public:
    Coor();
    Coor(double xcoorin, double ycoorin, double zcoorin); //Constructor
    Coor(const Coor& vec); //Overload Constructor 

    Coor& operator= (const Coor& vec);
    Coor& operator= (const double val);
    //Addition
    Coor operator+ (const Coor& vec) const;
    Coor& operator+= (const Coor& vec);
    Coor operator+ (const double val) const;
    Coor& operator+= (const double val);
    //Subtraction
    Coor operator- (const Coor& vec) const;
    Coor& operator-= (const Coor& vec);
    Coor operator- (const double val) const;
    Coor& operator-= (const double val);
    //Multiplication
    Coor operator* (const Coor& vec) const;
    Coor& operator*= (const Coor& vec);
    Coor operator* (const double val) const;
    Coor& operator*= (const double val);
    //Division
    Coor operator/ (const Coor& vec) const;
    Coor& operator/= (const Coor& vec);
    Coor operator/ (const double val) const;
    Coor& operator/= (const double val);

    double& x(){return xcoor;};
    double& y(){return ycoor;};
    double& z(){return zcoor;};

    Coor operator- () const;
    double dot (const Coor& vec) const; //Dot Product
    Coor cross (const Coor& vec) const; //Cross Product
    double norm () const;
};

#endif
