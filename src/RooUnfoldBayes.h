//=====================================================================-*-C++-*-
// File and Version Information:
//      $Id: RooUnfoldBayes.h,v 1.5 2010-01-19 15:33:47 adye Exp $
//
// Description:
//      Bayesian unfolding. Just an interface to RooUnfoldBayesImpl.
//
// Author: Tim Adye <T.J.Adye@rl.ac.uk>
//
//==============================================================================

#ifndef ROOUNFOLDBAYES_HH
#define ROOUNFOLDBAYES_HH

#include "RooUnfold.h"
#include <vector>
using std::vector;

class RooUnfoldResponse;
class TH1;
class TH2D;
class RooUnfoldBayesImpl;
class Array2D;

class RooUnfoldBayes : public RooUnfold {

public:

  // Standard methods

  RooUnfoldBayes(); // default constructor
  RooUnfoldBayes (const char*    name, const char*    title); // named constructor
  RooUnfoldBayes (const TString& name, const TString& title); // named constructor
  RooUnfoldBayes (const RooUnfoldBayes& rhs); // copy constructor
  RooUnfoldBayes& operator= (const RooUnfoldBayes& rhs); // assignment operator

  // Special constructors
  RooUnfoldBayes (const RooUnfoldResponse* res, const TH1* meas, Int_t niter= 4, Bool_t smoothit= false,
                  const char* name= 0, const char* title= 0);

  // Set up an existing object
  virtual RooUnfoldBayes& Clear ();
  virtual RooUnfoldBayes& Setup (const RooUnfoldBayes& rhs);
  virtual RooUnfoldBayes& Setup (const RooUnfoldResponse* res, const TH1* meas, Int_t niter= 4, Bool_t smoothit= false);

  virtual void Print (Option_t* option= "") const;

  static vector<Double_t>& H2VD (const TH1*  h, vector<Double_t>& v);
  static Array2D&          H2AD (const TH2D* h, Array2D& m, const TH1* norm= 0);
  static TVectorD&         VD2V (const vector<Double_t>& vd, TVectorD& v);
  static TMatrixD&         AD2M (const Array2D& ad, TMatrixD& m);

protected:

  virtual RooUnfoldBayes& Setup();
  virtual RooUnfoldBayes& Setup (Int_t niter, Bool_t smoothit);
  virtual Int_t unfold (vector<Double_t>& causes);
  virtual Int_t getCovariance() const;
  virtual void GetCov() const;  // actually updates mutable _cov

  // instance variables
  RooUnfoldBayesImpl* _bayes;
  Int_t _niter;
  Int_t _smoothit;

public:

  ClassDef (RooUnfoldBayes, 0) // Bayesian Unfolding
};

// Inline method definitions

inline RooUnfoldBayes::RooUnfoldBayes()                                           : RooUnfold()           {Setup();}
inline RooUnfoldBayes::RooUnfoldBayes (const char* name, const char* title)       : RooUnfold(name,title) {Setup();}
inline RooUnfoldBayes::RooUnfoldBayes (const TString& name, const TString& title) : RooUnfold(name,title) {Setup();}
inline RooUnfoldBayes& RooUnfoldBayes::operator= (const RooUnfoldBayes& rhs) {Assign(rhs); return *this;}

#endif
