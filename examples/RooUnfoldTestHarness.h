//=====================================================================-*-C++-*-
// File and Version Information:
//      $Id: RooUnfoldTestHarness.h,v 1.14 2010-01-26 00:53:15 adye Exp $
//
// Description:
//      Test Harness class for the RooUnfold package using toy MC generated
//      according to PDFs defined in RooUnfoldTestPdf.icc or RooUnfoldTestPdfRooFit.icc.
//
// Authors: Tim Adye <T.J.Adye@rl.ac.uk> and Fergus Wilson <fwilson@slac.stanford.edu>
//
//==============================================================================

#ifndef ROOUNFOLDTESTHARNESS_HH
#define ROOUNFOLDTESTHARNESS_HH

#if !defined(__CINT__) || defined(__MAKECINT__)
#include "TNamed.h"
#if ROOT_VERSION_CODE >= ROOT_VERSION(5,0,0)
#include "TVectorDfwd.h"
#else
class TVectorD;
#endif
#endif

#ifdef __CINT__
#include "TPad.h"    // Why is this necessary (at least in ROOT 5.18+)?
#include "ArgVars.h"
#else
class ArgVars;
#endif

class TCanvas;
class TPostScript;
class TH1;
class TH1D;
class TH2D;
class RooUnfoldResponse;
class RooUnfold;

class RooUnfoldTestHarness : public TNamed {
public:
  // Parameters
  Int_t    method, stage, ftrainx, ftestx, ntx, ntest, ntrain;
  Int_t    regparm, ntoys, nmx, onepage, doerror, dim, dosmear, nbPDF;
  Double_t xlo, xhi, mtrainx, wtrainx, btrainx, mtestx, wtestx, btestx;
  Double_t effxlo, effxhi, xbias, xsmear;

  // Data
  Int_t              error, ipad, ntbins, nmbins;
  TCanvas*           canvas;
  TPostScript*       ps;
  TH1                *hTrain, *hTrainTrue, *hTrue, *hMeas, *hReco, *hRes, *hPulls;
  TH1D               *hPDFx, *hTestPDFx;
  TH2D               *hResmat;
  RooUnfoldResponse* response;
  RooUnfold*         unfold;

  // Constructors
  RooUnfoldTestHarness (const char* name= "RooUnfoldTest");
  RooUnfoldTestHarness (const char* name, const char* args);
  RooUnfoldTestHarness (const char* name, int argc, const char* const* argv);
  virtual ~RooUnfoldTestHarness();

  // Methods and functions
  virtual void     Parms (ArgVars& args);
  virtual Int_t    Run();
  virtual void     SetupCanvas();
  virtual Int_t    RunTests();
  virtual Int_t    Train();
  virtual Int_t    Test();
  virtual Int_t    Unfold();
  virtual void     ShowTest();
  virtual void     Results();
  virtual TH1D*    Generate (TVectorD& x, const char* name, const char* title, Int_t nt, Int_t fpdf,
                             Int_t nx, Double_t xlo, Double_t xhi, Double_t bkg, Double_t mean, Double_t width);
  virtual bool     Eff      (Double_t x,           Double_t xlo, Double_t xhi, Double_t efflo, Double_t effhi) const;
  virtual Double_t Smear    (Double_t x, Int_t nt, Double_t xlo, Double_t xhi, Double_t bias,  Double_t smear) const;
  virtual void     Reset();
  virtual void     SetDefaults();
  virtual int      SetArgs  (int argc, const char* const* argv, bool split= false);
  virtual void     Init();
  virtual Int_t    CheckParms();
  virtual void     PrintParms (std::ostream& o) const;
  static  void     setmax   (TH1* h, const TH1* h1= 0, const TH1* h2= 0, const TH1* h3= 0,
                             const TH1* h4= 0, const TH1* h5= 0, const TH1* h6= 0);
};

#ifndef NOINLINE
#include "RooUnfoldTestHarness.icc"
#endif

#endif
