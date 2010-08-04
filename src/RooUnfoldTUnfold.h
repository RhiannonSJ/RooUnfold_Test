#ifndef ROOUNFOLDTUNFOLD_H_
#define ROOUNFOLDTUNFOLD_H_

#include "RooUnfold.h"

class RooUnfoldResponse;
class TH1;
class TH1D;
class TH2D;
class TUnfold;

class RooUnfoldTUnfold : public RooUnfold {

public:

  virtual ~RooUnfoldTUnfold(); // destructor
  virtual RooUnfoldTUnfold* Clone (const char* newname= 0) const;
  RooUnfoldTUnfold (const RooUnfoldResponse* res, const TH1* meas, Int_t kterm= 1,
                const char* name= 0, const char* title= 0);
	void Reset();
	
	TObject* Impl();
	void FixTau(Double_t tau);
	void OptimiseTau();
	virtual void SetRegParm(Double_t parm){FixTau(parm);}
	Double_t GetTau() const { return _tau;    }
	virtual Double_t GetRegParm() const { return GetTau(); }
	virtual void Get_settings();
	
protected:
  TH2D* ematrix; //matrix of errors
  void Init();
  virtual void Unfold();
  void GetCov();
  Int_t _kterm; //Regularisation parameter
  TUnfold* _unf; //TUnfold object
  TH2D* Flip2D(const TH2D* Hres);
  TH2D* CopyOverflow(const TH2D* h)const;
  Bool_t tau_set;
  Double_t _tau;
  
public:

  ClassDef (RooUnfoldTUnfold, 0) 
};

#endif /*ROOUNFOLDTUNFOLD_H_*/
