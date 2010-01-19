//=====================================================================-*-C++-*-
// File and Version Information:
//      $Id: RooUnfoldBayes.cxx,v 1.10 2010-01-19 23:30:59 adye Exp $
//
// Description:
//      Bayesian unfolding. Just an interface to RooUnfoldBayesImpl.
//      Note that we use 1D distributions in RooUnfoldBayesImpl, even if we are
//      unfolding multi-dimensional distributions: RooUnfold already converted
//      to 1D. Except for smoothing (which RooUnfoldBayesImpl doesn't implement)
//      this is just a matter of bookkeeping.
//
// Author: Tim Adye <T.J.Adye@rl.ac.uk>
//
//==============================================================================

#include "RooUnfoldBayes.h"

#include <iostream>
#include <vector>
#include <math.h>

#include "TNamed.h"
#include "TH1.h"
#include "TH2.h"

#include "RooUnfoldResponse.h"
#include "Array2D.h"
#include "RooUnfoldBayesImpl.h"

using std::vector;
using std::cerr;
using std::endl;

ClassImp (RooUnfoldBayes);

RooUnfoldBayes::RooUnfoldBayes (const RooUnfoldBayes& rhs)
  : RooUnfold (rhs.GetName(), rhs.GetTitle())
{
  Setup ();
  Setup (rhs);
}

RooUnfoldBayes::RooUnfoldBayes (const RooUnfoldResponse* res, const TH1* meas, Int_t niter, Bool_t smoothit,
                                      const char* name, const char* title)
  : RooUnfold (name, title)
{
  Setup();
  Setup (res, meas, niter, smoothit);
}

RooUnfoldBayes&
RooUnfoldBayes::Clear()
{
  delete _bayes;
  RooUnfold::Clear();
  return *this;
}

RooUnfoldBayes&
RooUnfoldBayes::Setup()
{
  _bayes= 0;
  _niter= 0;
  _smoothit= false;
  return *this;
}

RooUnfoldBayes&
RooUnfoldBayes::Setup (const RooUnfoldBayes& rhs)
{
  return Setup (rhs._niter, rhs._smoothit);
}

RooUnfoldBayes&
RooUnfoldBayes::Setup (const RooUnfoldResponse* res, const TH1* meas, Int_t niter, Bool_t smoothit)
{
  RooUnfold::Setup (res, meas);
  return Setup (niter, smoothit);
}

RooUnfoldBayes&
RooUnfoldBayes::Setup (Int_t niter, Bool_t smoothit)
{
  _niter= niter;
  _smoothit= smoothit;

  _bayes= new RooUnfoldBayesImpl (GetName(), GetTitle());

  _bayes->build (1, vector<Int_t>(1,_nt), vector<Double_t>(1,0.0), vector<Double_t>(1,1.0),
                 1, vector<Int_t>(1,_nm), vector<Double_t>(1,0.0), vector<Double_t>(1,1.0));

  if (verbose() >= 2) Print();

  vector<Double_t> vtruth(_nt), vtrain(_nm), vmeasured(_nm);
  Array2D aresp(_nt,_nm);
  _bayes->setupTrain  (H2VD (_res->Htruth(),    vtruth),
                       H2VD (_res->Hmeasured(), vtrain),
                       H2AD (_res->Hresponse(), aresp));
  _bayes->setupUnfold (H2VD (_meas,             vmeasured));

  if (verbose() >= 2) Print();

  vector<Double_t> causes;
  unfold (causes);

  if (verbose() >= 2) Print();

  _rec.ResizeTo (_nm);
  VD2V (causes, _rec);
  _haveCov= false;

  return *this;
}

void
RooUnfoldBayes::GetCov() const
{
  _cov.ResizeTo (_nm, _nm);
  getCovariance();
  if (_bayes->error() != 0.0) {
    AD2M (_bayes->covariance(), _cov);
  } else {
    cerr << "Covariance matrix not calculated - fill errors with sqrt(N)" << endl;
    for (Int_t i= 0; i < _nm; i++)
      _cov(i,i)= sqrt (fabs (_rec(i)));
  }
  _haveCov= true;
}

// routines that can be overridden by RooUnfoldBinByBin.
Int_t RooUnfoldBayes::unfold (vector<Double_t>& causes) { return _bayes->unfold (causes, _niter, _smoothit); }
Int_t RooUnfoldBayes::getCovariance() const             { return _bayes->getCovariance();                    }

void
RooUnfoldBayes::Print(Option_t* option) const
{
  RooUnfold::Print (option);
  _bayes->info();
}

vector<Double_t>&
RooUnfoldBayes::H2VD (const TH1* h, vector<Double_t>& v)
{
  if (!h) return v;
  Int_t nb= v.size();
  for (size_t i= 0; i < nb; i++)
    v[i]= h->GetBinContent (RooUnfoldResponse::GetBin(h,i));
  return v;
}

Array2D&
RooUnfoldBayes::H2AD (const TH2D* h, Array2D& m, const TH1* norm)
{
  if (!h) return m;
  Int_t nx= m.GetNrows(), ny= m.GetNcols();
  for (size_t j= 0; j < ny; j++) {
    Double_t nTrue= norm ? norm->GetBinContent (j+1) : 1.0;
    if (nTrue == 0.0) {
      for (size_t i= 0; i < nx; i++)
        m.Set (j, i, 0.0);
    } else {
      for (size_t i= 0; i < nx; i++)
        m.Set (j, i, h->GetBinContent(i+1,j+1) / nTrue);
    }
  }
  return m;
}

TVectorD&
RooUnfoldBayes::VD2V (const vector<Double_t>& vd, TVectorD& v)
{
  Int_t nb= v.GetNrows();
  for (size_t i= 0; i < nb; i++)
    v(i)= vd[i];
  return v;
}

TMatrixD&
RooUnfoldBayes::AD2M (const Array2D& ad, TMatrixD& m)
{
  Int_t nx= m.GetNrows(), ny= m.GetNcols();
  for (size_t i= 0; i < nx; i++)
    for (size_t j= 0; j < ny; j++)
      m(i,j)= ad.Get(j,i);
  return m;
}
