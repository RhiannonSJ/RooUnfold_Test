

#include "RooUnfoldTUnfold.h"

#include <iostream>
#include <iomanip>
#include <cmath>

#include "TNamed.h"
#include "TH1.h"
#include "TH2.h"
#include "TVectorD.h"
#include "TMatrixD.h"
#include "TUnfold.h"
#include "TGraph.h"

#include "RooUnfoldResponse.h"

using std::cout;
using std::cerr;
using std::endl;
using std::sqrt;

ClassImp (RooUnfoldTUnfold);


RooUnfoldTUnfold::RooUnfoldTUnfold (const RooUnfoldResponse* res, const TH1* meas, Int_t kterm,
                            const char* name, const char* title)
  : RooUnfold (res, meas, name, title), _kterm(kterm)
{
  Init();
}

RooUnfoldTUnfold::~RooUnfoldTUnfold()
{
	//Destructor
	Reset();
	delete ematrix;
	delete _unf;
}

RooUnfoldTUnfold*
RooUnfoldTUnfold::Clone (const char* newname) const
{
	//Clones object
  RooUnfoldTUnfold* unfold= new RooUnfoldTUnfold(*this);
  if (newname && strlen(newname)) unfold->SetName(newname);
  return unfold;
}

void
RooUnfoldTUnfold::Reset()
{
	//Resets all values
  Init();
  RooUnfold::Reset();
}


void
RooUnfoldTUnfold::Init()
{
	//Sets error matrix
	tau_set=false;
	_tau=0;
	ematrix=0;
	_unf=0;
}

TObject*
RooUnfoldTUnfold::Impl()
{
	return _unf;
}


void
RooUnfoldTUnfold::Unfold()
{
	/* Does the unfolding
	   The _kterm value decides the regularisation method.
	 k=0: No regularisation
	 k=1: Minimise x-x0
	 k=2: Minimise 1st derivative of (x-x0)
	 k=3: Minimise 2nd derivative of (x-x0)
	 Uses TUnfold.ScanLCurve to do unfolding. 
	 Creates covariance matrix. 
	 */
  const TH2D* Hres=_res->Hresponse();
  if (_fail) return;
  TUnfold::ERegMode regmode=TUnfold::kRegModeNone;
  switch (_kterm){
  	case 0:
  	regmode=TUnfold::kRegModeNone;
  	break;
  	case 1:
  	regmode=TUnfold::kRegModeSize;
  	break;
  	case 2:
  	regmode=TUnfold::kRegModeDerivative;
  	break;
  	case 3:
  	regmode=TUnfold::kRegModeCurvature;
  	break;
  	default:
  	regmode=TUnfold::kRegModeSize;
  }
  TH2D* Hresc=CopyOverflow(Hres);
  TH2D* Hres_flipped=Flip2D(Hresc);
  _unf= new TUnfold(Hres_flipped,TUnfold::kHistMapOutputHoriz,regmode);
  Int_t nScan=30;
  // use automatic L-curve scan: start with taumin=taumax=0.0
  Double_t tauMin=0.0;
  Double_t tauMax=0.0;
  Int_t iBest;
  TSpline *logTauX,*logTauY;
  TGraph *lCurve;
  // this method scans the parameter tau and finds the kink in the L curve
  // finally, the unfolding is done for the best choice of tau
#if ROOT_VERSION_CODE >= ROOT_VERSION(5,23,0)  /* TUnfold v6 (included in ROOT 5.22) didn't have setInput return value */
  if(_unf->SetInput(_meas,0.0)>=10000) {
    cerr<<"Unfolding result may be wrong\n";
  }
#else
  _unf->SetInput(_meas,0.0);
#endif
  _unf->SetInput(_meas);
  if (!tau_set){
  iBest=_unf->ScanLcurve(nScan,tauMin,tauMax,&lCurve,&logTauX,&logTauY);
  }
  else{
  	_unf->DoUnfold(_tau);
  }
  TH1D* reco=_unf->GetOutput("_rec","reconstructed dist",0,0);
  if (_overflow){
	  _rec.ResizeTo (reco->GetNbinsX());
	  for (int i=1;i<reco->GetNbinsX();i++){
	  	_rec(i)=(reco->GetBinContent(i));
	  }
  }
  else{
  	_rec.ResizeTo (reco->GetNbinsX());
	  for (int i=0;i<reco->GetNbinsX();i++){
	  	_rec(i)=reco->GetBinContent(i+1);
	  }
  }
  ematrix=_unf->GetEmatrix("ematrix","error matrix",0,0);
  delete Hresc;
  delete reco;
  delete Hres_flipped;
  _unfolded= true;
  _haveCov=  false;
}

void
RooUnfoldTUnfold::GetCov()
{
	//Gets Covariance matrix
	Int_t nt=_meas->GetNbinsX();
	if (_overflow){nt+=2;}
	if (!_unfolded) Unfold();
	if (_fail) return;
	_cov.ResizeTo (nt,nt);
	for (Int_t i= 0; i<nt; i++) {
		for (Int_t j= 0; j<nt; j++) {
			_cov (i,j)= ematrix->GetBinContent(i+1,j+1);
		}
	}
	delete ematrix;
	_haveCov= true;
}


TH2D*
RooUnfoldTUnfold::Flip2D(const TH2D* h)
{
	//Returns the transpose of a matrix expressed as a TH2D
	//Inefficiencies are returned in the underflow bin
  Int_t nx= h->GetNbinsX(), ny= h->GetNbinsY();
  Double_t xlo= h->GetXaxis()->GetXmin(), xhi= h->GetXaxis()->GetXmax();
  Double_t ylo= h->GetYaxis()->GetXmin(), yhi= h->GetYaxis()->GetXmax();
  TH2D* hx= new TH2D ("h_flipped", h->GetTitle(), ny, ylo, yhi, nx, xlo, xhi);
  if (nx<ny){
  	cerr<<"Warning: fewer x bins than y bins. Unfolding may not work correctly"<<endl;
  } 
  for (Int_t i= 1; i <=ny; i++) {
  	Double_t sumineff=0;
    for (Int_t j= 1; j <=nx; j++) {
      hx->SetBinContent (j, i, h->GetBinContent (i, j));
      hx->SetBinError   (j, i, h->GetBinError   (i, j));
      sumineff+=h->GetBinContent(i,j); 
    }
    hx->SetBinContent(i,0,_res->Htruth()->GetBinContent(i)-sumineff);
  }
  return hx;
}

TH2D*
RooUnfoldTUnfold::CopyOverflow (const TH2D* h) const
{
  if (!_overflow) return (TH2D*)h->Clone();
  Int_t nx= h->GetNbinsX(), ny= h->GetNbinsX();
  Double_t xlo= h->GetXaxis()->GetXmin(), xhi= h->GetXaxis()->GetXmax(), xb= (xhi-xlo)/nx;
  Double_t ylo= h->GetYaxis()->GetXmin(), yhi= h->GetYaxis()->GetXmax(), yb= (yhi-ylo)/ny;
  nx += 2; ny += 2;
  TH2D* hx= new TH2D (h->GetName(), h->GetTitle(), nx, xlo-xb, xhi+xb, ny, ylo-yb, yhi+yb);
  for (Int_t i= 0; i < nx; i++) {
    for (Int_t j= 0; j < ny; j++) {
      hx->SetBinContent (i+1, j+1, h->GetBinContent (i, j));
      hx->SetBinError   (i+1, j+1, h->GetBinError   (i, j));
    }
  }
  return hx;
}

void 
RooUnfoldTUnfold::FixTau(Double_t tau)
{
	_tau=tau;
	tau_set=true;
}


void
RooUnfoldTUnfold::OptimiseTau()
{
	tau_set=false;
}

void
RooUnfoldTUnfold::Get_settings()
{
	_minparm=0;
	_maxparm=2;
	_stepsizeparm=5e-2;
	_defaultparm=1;
}
