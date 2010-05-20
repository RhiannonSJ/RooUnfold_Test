//=====================================================================-*-C++-*-
// File and Version Information:
//      $Id: ArgVars.h,v 1.6 2010-05-20 22:50:57 adye Exp $
//
// Description:
//      Parse argument list for parameter settings
//
// Author: Tim Adye <T.J.Adye@rl.ac.uk>
//
//==============================================================================

#ifndef ARGVARS_H
#define ARGVARS_H

#if !defined(__CINT__) || defined(__MAKECINT__)
#include "TROOT.h"
#include <iostream>
#include "TObject.h"
#include "TList.h"
#endif

class ArgVar : public TNamed {
public:
  const char *defhelp;
  Int_t      *ivar, idef;
  Double_t   *fvar, fdef;
  bool       setdef;
  ArgVar (const char *n,    Int_t *v) : TNamed (n,""), ivar(v), idef(0), fvar(0), fdef(0), setdef(false) {}
  ArgVar (const char *n, Double_t *v) : TNamed (n,""), ivar(0), idef(0), fvar(v), fdef(0), setdef(false) {}
  ArgVar (const char *n,    Int_t *v, Int_t    d, const char* h=0, const char* dh=0)
    : TNamed(n,h), defhelp(dh), ivar(v), idef(d), fvar(0), fdef(0), setdef(true) {}
  ArgVar (const char *n, Double_t *v, Double_t d, const char* h=0, const char* dh=0)
    : TNamed(n,h), defhelp(dh), ivar(0), idef(0), fvar(v), fdef(d), setdef(true) {}
};

class ArgVars : public TObject {
private:
  TList lst;
  static bool CmpOpt (const char* p, const char*  opt, const char* s);
  ArgVar* Find (const char* name) const;
  ArgVars& Add (ArgVar* arg);
public:
  ArgVars() {}
  ~ArgVars();
  ArgVars& Add (const char* name,    Int_t* var) { return Add (new ArgVar (name, var)); }
  ArgVars& Add (const char* name, Double_t* var) { return Add (new ArgVar (name, var)); }
  ArgVars& Add (const char* name,    Int_t* var, Int_t    def, const char* help=0, const char* defhelp=0)
    { return Add (new ArgVar (name, var, def, help, defhelp)); }
  ArgVars& Add (const char* name, Double_t* var, Double_t def, const char* help=0, const char* defhelp=0)
    { return Add (new ArgVar (name, var, def, help, defhelp)); }
  ArgVars& Add (const ArgVars& args);
  ArgVars& SetDefault (const char* name, Int_t    def);
  ArgVars& SetDefault (const char* name, Double_t def);
  Int_t SetArgs (int argc, const char* const* argv, bool split= false) const;
  void  SetDefaults() const;
  virtual void  Print (std::ostream& o, const char* sep= " ") const;
  virtual void  Print (const char* sep= " ") const { Print (std::cout, sep); }
  void  Usage (const char* prog) const;
  void  ArgHelp (std::ostream& o) const;
};

#ifndef NOINLINE
#include "ArgVars.icc"
#endif

#endif
