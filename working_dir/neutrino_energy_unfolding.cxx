//=====================================================================-*-C++-*-
// File and Version Information:
//      $Id$
//
// Description:
//      Simple example usage of the RooUnfold package using toy MC.
//
// Authors: Tim Adye <T.J.Adye@rl.ac.uk> and Fergus Wilson <fwilson@slac.stanford.edu>
//
//==============================================================================

#include <iostream>
using std::cout;
using std::endl;

#include "roo_unfold.h"

#include "../src/RooUnfoldResponse.h"
#include "../src/RooUnfoldBayes.h"
//#include "RooUnfoldSvd.h"
//#include "RooUnfoldTUnfold.h"

//==============================================================================
// Calculate the reconstructed energy
//==============================================================================

void EReco( vector< double > E_mu, vector< double > p_mu, vector< double > cos_mu, vector< double > &E_nu_reco ){

    // Calculate the reconstructed energy and fill a vector of doubles with the result
    // Define variables
    double m_N  = 0.93828; // Nucleon madd, GeV
    double m_mu = 0.10566; // Muon mass, GeV

    // Loop over the vectors and do the calculation
    
    if ( E_mu.size() != p_mu.size() 
      || E_mu.size() != cos_mu.size() 
      || p_mu.size() != cos_mu.size() ){
      
        cerr << " Vectors need to be the same size " << endl;
        exit(1);
    
    }

    int vect_size = E_mu.size();

    for ( int i = 0; i < vect_size; ++i ){
        double e_reco;

        e_reco = ( 1 / ( 1 - ( ( 1 / m_N ) * ( E_mu[i] - p_mu[i] * cos_mu[i] ) ) ) ) * ( E_mu[i] - ( ( 1 / ( 2 * m_N ) ) * m_mu * m_mu ) );

        E_nu_reco.push_back( e_reco );
    }

}

//==============================================================================
// The main function
//==============================================================================
void neutrino_energy_unfolding() { 
    //==============================================================================
    // Reading in the root file 
    //==============================================================================
    TFile f("/hepstore/rjones/Exercises/Flavours/Default+MEC/sbnd/1M/gntp.10000.gst.root");
    if(f.IsZombie()){
        cerr << " Error opening file " << endl;
        exit(1);
    }
    else{
        cout << "============================== Default + MEC open =============================" << endl;
    }

    //==============================================================================
    // Get everything from the tree
    //==============================================================================

    TTree *gst = (TTree*) f.Get("gst");

    TBranch *b_nu_e  = gst->GetBranch("Ev");
    TBranch *b_mu_e  = gst->GetBranch("El");
    TBranch *b_mu_p  = gst->GetBranch("pl");
    TBranch *b_theta = gst->GetBranch("cthl");
    TBranch *b_nfpi0 = gst->GetBranch("nfpi0");
    TBranch *b_nfpip = gst->GetBranch("nfpip");
    TBranch *b_nfpim = gst->GetBranch("nfpim");
    TBranch *b_cc    = gst->GetBranch("cc");
    TBranch *b_nc    = gst->GetBranch("nc");
    TBranch *b_pdgf  = gst->GetBranch("pdgf");
    TBranch *b_fspl  = gst->GetBranch("fspl");

    // Define and fill vectors for the muon energy, momentum and opening angle
    // Define an empty vector to hold the reconstructed energy
    vector< double > mu_e;
    vector< double > mu_p;
    vector< double > mu_cth;
    vector< double > nu_e;
    vector< double > nu_e_reco;

    int n_entries = gst->GetEntries();

    // Fill the vectors
    for ( int i = 0; i < n_entries; ++i ){
    
        gst->GetEntry(i);
    
        // Fill vectors if primary lepton is a muon and it is a cc0pi event
        if ( b_fspl->GetLeaf( "fspl" )->GetValue() == 13
          && b_nfpi0->GetLeaf( "nfpi0" )->GetValue() == 0 
          && b_nfpip->GetLeaf( "nfpip" )->GetValue() == 0 
          && b_nfpim->GetLeaf( "nfpim" )->GetValue() == 0
          && b_cc->GetLeaf( "cc" )->GetValue() == 1 ){
    
            nu_e.push_back(   b_nu_e->GetLeaf( "Ev" )->GetValue() );
            mu_e.push_back(   b_mu_e->GetLeaf( "El" )->GetValue() );
            mu_p.push_back(   b_mu_p->GetLeaf( "pl" )->GetValue() );
            mu_cth.push_back( b_theta->GetLeaf( "cthl" )->GetValue() );
    
        }
    }

    // Get the reconstructed energies
    EReco( mu_e, mu_p, mu_cth, nu_e_reco );

    //==============================================================================
    // Define the histograms
    //==============================================================================
    
    TCanvas *c    = new TCanvas( "c", "Unfolding", 800, 600 );

    TLegend *l    = new TLegend( 0.68, 0.68, 0.88, 0.88 );

    TH1D *hTrue   = new TH1D( "hTrue",   "Incoming neutrino energy, CC0#pi", 200./4., 0., 10./4. );
    TH1D *hTrue1  = new TH1D( "hTrue1",  "Incoming neutrino energy, training, CC0#pi", 200./4., 0., 10./4. );
    TH1D *hTrue2  = new TH1D( "hTrue2",  "Incoming neutrino energy, unfolding, CC0#pi", 200./4., 0., 10./4. );

    TH1D *hReco   = new TH1D( "hReco",   "Reconstructed neutrino energy, CC0#pi", 200./4., 0., 10./4. );
    TH1D *hReco1  = new TH1D( "hReco1",  "Reconstructed neutrino energy, training, CC0#pi", 200./4., 0., 10./4. );
    TH1D *hReco2  = new TH1D( "hReco2",  "Reconstructed neutrino energy, unfolding, CC0#pi", 200./4., 0., 10./4. );
    
    //==============================================================================
    // Train the unfolding algorithm and get the response matrix
    //==============================================================================

    int n_half = TMath::Floor( double( nu_e_reco.size() ) / 2. );
    int n_all  = nu_e_reco.size();
    
    cout << "==================================== TRAIN ====================================" << endl;
    RooUnfoldResponse response ( 200./4., 0., 10./4. );

    for ( int i = 0; i < n_half; ++i ){

        // Fill the histograms
        hTrue1->Fill( nu_e[i] );
        hReco1->Fill( nu_e_reco[i] );
        
        // Fill the response histograms
        response.Fill( nu_e_reco[i], nu_e[i] );
    }

    for ( int i = n_half; i < n_all; ++i ){

        // Fill the histograms
        hTrue2->Fill( nu_e[i] );
        hReco2->Fill( nu_e_reco[i] );
    }

    cout << "==================================== UNFOLD ===================================" << endl;
    RooUnfoldBayes    unfold   ( &response, hReco2, 1 ); // Try different numbers of iterations

    // Unfold
    // Histogram output
    TH1D *hUnfold =  (TH1D*) unfold.Hreco();
    
    // Vector and covariance matrix output
    TVectorD unfolded_dist = unfold.Vreco();
    TMatrixD unfolded_errs = unfold.Ereco;

    l->AddEntry( hUnfold, "Unfolded", "p" );
    l->AddEntry( hReco2, "Reconstructed", "l" );
    l->AddEntry( hTrue2, "True", "l" );

    
    unfold.PrintTable (cout, hTrue2);
   
    // double norm = hUnfold->Integral();

    hUnfold->GetYaxis()->SetTitleOffset(1.5);
    hUnfold->SetStats(kFALSE);
    hUnfold->GetXaxis()->SetTitle("E_{#nu}");
    hUnfold->GetYaxis()->SetTitle("Number of events");
    hUnfold->SetTitle("Unfolded reconstructed E_{#nu}");
    hUnfold->SetMarkerStyle(2);
    //hUnfold->Scale(1/norm);
    hUnfold->Draw();
    
    hReco2->SetLineColor( kRed + 2 );
    //hReco2->Scale(1/norm);
    hReco2->Draw("SAME");
    
    hTrue2->SetLineColor( kGreen + 2 );
    //hTrue2->Scale(1/norm);
    hTrue2->Draw("SAME");
   
    l->Draw();

    c->SaveAs( "my_work/unfolded_distributions/1D_unfolding_ex.png" );

    delete hUnfold;

    delete hReco;
    delete hReco1;
    delete hReco2;

    delete hTrue;
    delete hTrue1;
    delete hTrue2;
    
    delete c;
    delete l;

} 
