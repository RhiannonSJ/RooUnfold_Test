//=====================================================================-*-C++-*-
// File and Version Information:
//      $Id: RooUnfoldBayes.cxx,v 1.16 2010-07-19 21:45:10 adye Exp $
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
  : RooUnfold (rhs)
{
  Init();
  CopyData (rhs);
}

RooUnfoldBayes::RooUnfoldBayes (const RooUnfoldResponse* res, const TH1* meas, Int_t niter, Bool_t smoothit,
                                const char* name, const char* title)
  : RooUnfold (res, meas, name, title), _niter(niter), _smoothit(smoothit)
{
  Init();
}

RooUnfoldBayes*
RooUnfoldBayes::Clone (const char* newname) const
{
  RooUnfoldBayes* unfold= new RooUnfoldBayes(*this);
  if (newname && strlen(newname)) unfold->SetName(newname);
  return unfold;
}

void
RooUnfoldBayes::Init()
{
  _bayes= 0;
}

void
RooUnfoldBayes::Destroy()
{
  delete _bayes;
}

void
RooUnfoldBayes::Reset()
{
  Destroy();
  Init();
  RooUnfold::Reset();
}

void
RooUnfoldBayes::Assign (const RooUnfoldBayes& rhs)
{
  RooUnfold::Assign (rhs);
  CopyData (rhs);
}

void
RooUnfoldBayes::CopyData (const RooUnfoldBayes& rhs)
{
  _niter=    rhs._niter;
  _smoothit= rhs._smoothit;
}

TObject*
RooUnfoldBayes::Impl()
{
  return _bayes;
}

void
RooUnfoldBayes::Unfold()
{
  _bayes= new RooUnfoldBayesImpl (GetName(), GetTitle());

  Int_t overflow= (_overflow ? 2 : 0), nt= _nt+overflow, nm= _nm+overflow;
  _bayes->build (1, vector<Int_t>(1,nt), vector<Double_t>(1,0.0), vector<Double_t>(1,1.0),
                 1, vector<Int_t>(1,nm), vector<Double_t>(1,0.0), vector<Double_t>(1,1.0));

  _bayes->setDebug (verbose());
  if (verbose() >= 2) Print();

  vector<Double_t> vtruth(nt), vtrain(nm), vmeasured(nm);
  Array2D aresp(nt,nm);
  _bayes->setupTrain  (H2VD (_res->Htruth(),    vtruth, _overflow),
                       H2VD (_res->Hmeasured(), vtrain, _overflow),
                       H2AD (_res->Hresponse(), aresp,  0, _overflow));
  _bayes->setupUnfold (H2VD (_meas,             vmeasured, _overflow));

  if (verbose() >= 2) Print();

  vector<Double_t> causes;
  unfold (causes);

  if (verbose() >= 2) Print();

  _rec.ResizeTo (nt);
  VD2V (causes, _rec);
  _unfolded= true;
  _haveCov=  false;
}

void
RooUnfoldBayes::GetCov()
{
  if (!_unfolded) Unfold();
  Int_t nt= _nt + (_overflow ? 2 : 0);
  _cov.ResizeTo (nt, nt);
  getCovariance();
  if (_bayes->error() != 0.0) {
    AD2M (_bayes->covariance(), _cov);
  } else {
    cerr << "Covariance matrix not calculated - fill errors with sqrt(N)" << endl;
    for (Int_t i= 0; i < nt; i++)
      _cov(i,i)= fabs (_rec(i));
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
  if (_bayes) _bayes->info();
}

vector<Double_t>&
RooUnfoldBayes::H2VD (const TH1* h, vector<Double_t>& v, Bool_t overflow)
{
  if (!h) return v;
  size_t nb= v.size();
  for (size_t i= 0; i < nb; i++)
    v[i]= h->GetBinContent (RooUnfoldResponse::GetBin(h,i,overflow));
  return v;
}

Array2D&
RooUnfoldBayes::H2AD (const TH2D* h, Array2D& m, const TH1* norm, Bool_t overflow)
{
  if (!h) return m;
  Int_t first= overflow ? 0 : 1;
  Int_t nt= m.GetNrows(), nm= m.GetNcols();
  for (Int_t j= 0; j < nt; j++) {
    Double_t nTrue= norm ? norm->GetBinContent (j+first) : 1.0;
    if (nTrue == 0.0) {
      for (Int_t i= 0; i < nm; i++)
        m.Set (j, i, 0.0);
    } else {
      for (Int_t i= 0; i < nm; i++)
        m.Set (j, i, h->GetBinContent(i+first,j+first) / nTrue);
    }
  }
  return m;
}

TVectorD&
RooUnfoldBayes::VD2V (const vector<Double_t>& vd, TVectorD& v)
{
  Int_t nb= v.GetNrows();
  for (Int_t i= 0; i < nb; i++)
    v(i)= vd[i];
  return v;
}

TMatrixD&
RooUnfoldBayes::AD2M (const Array2D& ad, TMatrixD& m)
{
  Int_t nm= m.GetNrows(), nt= m.GetNcols();
  for (Int_t i= 0; i < nm; i++)
    for (Int_t j= 0; j < nt; j++)
      m(i,j)= ad.Get(j,i);
  return m;
}
