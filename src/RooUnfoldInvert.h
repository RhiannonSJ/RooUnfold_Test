#ifndef ROOUNFOLDINVERT_H_
#define ROOUNFOLDINVERT_H_

#include "RooUnfold.h"

class RooUnfoldResponse;
class TH1;
class TH1D;
class TH2D;

class RooUnfoldInvert : public RooUnfold {

public:


 RooUnfoldInvert(); // default constructor
  RooUnfoldInvert (const char*    name, const char*    title); // named constructor
  RooUnfoldInvert (const TString& name, const TString& title); // named constructor
  RooUnfoldInvert (const RooUnfoldInvert& rhs); // copy constructor
  virtual ~RooUnfoldInvert(); // destructor
  RooUnfoldInvert& operator= (const RooUnfoldInvert& rhs); // assignment operator
  virtual RooUnfoldInvert* Clone (const char* newname= 0) const;
RooUnfoldInvert (const RooUnfoldResponse* res, const TH1* meas, const char* name=0, const char* title=0);


protected:
virtual void Unfold();
virtual void GetCov();
virtual void GetSettings();

private:
  TMatrixD Hres_i;
  TVectorD Hmeas_m;
  
public:

  ClassDef (RooUnfoldInvert, 0) 
};

inline RooUnfoldInvert::RooUnfoldInvert()                                           : RooUnfold()           {Init();}
inline RooUnfoldInvert::RooUnfoldInvert (const char* name, const char* title)       : RooUnfold(name,title) {Init();}
inline RooUnfoldInvert::RooUnfoldInvert (const TString& name, const TString& title) : RooUnfold(name,title) {Init();}
inline RooUnfoldInvert& RooUnfoldInvert::operator= (const RooUnfoldInvert& rhs) {Assign(rhs); return *this;}

#endif /*ROOUNFOLDINVERT_H_*/
