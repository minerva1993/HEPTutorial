#include <TFitResult.h>
#include <TMath.h>

void zMassFit(){

  //read root file
  TFile *_file0 = TFile::Open("hist_data.root");
  TH1D *h_Mmumu=(TH1D*)_file0->Get("Mmumu");

  TFile *_file1 = TFile::Open("hist_stack.root");
  TH1F *h_Mmumumc=(TH1F*)_file1->Get("Mmumu");

  //create clone of hgg dist
  TH1F *h_mmumu = (TH1F*)h_Mmumu->Clone();
  TH1F *h_mmumumc = (TH1F*)h_Mmumumc->Clone();


  //define gaus fit fct
  TF1 *mygaus = new TF1("mygaus", "[0]/(sqrt(2.0*3.14151927)*[2])*exp(-(x-[1])*(x-[1])/(2.0*[2]*[2]))", 30.0, 150.0);
    mygaus->SetParName(0, "Normalize");
    mygaus->SetParName(1, "Mean");
    mygaus->SetParName(2, "Width");
    double initval1[3] = { 13000, 90.0, 30};
    mygaus->SetParameters(initval1);

  //create canvas
  TCanvas *c=new TCanvas("c","DiMuon(Z) Mass",600,600);

  //store fit result
  TFitResultPtr r= h_mmumu->Fit(mygaus, "S");
  Double_t par0   = r->Parameter(0);
  Double_t par1   = r->Parameter(1);
  Double_t par2   = r->Parameter(2);
  double initval3[3] = {par0,par1,par2};

  //from TFitResult, we create bkg2 which will appear on our canvas
  TF1*Mmumu = new TF1("Mmumu", "[8]*exp(-[9]*(x/100)+[10]*(x/100)*(x/100))", 30.0, 150.0);
  Mmumu->SetParameters(initval3);

  //draw
  c->cd();
  h_mmumu->Draw("P");
  h_mmumu->SetMarkerStyle(8);
  h_mmumu->SetMarkerSize(0.6);
  h_mmumu->GetYaxis()->SetRangeUser(0,3000);
  h_mmumu->GetXaxis()->SetRangeUser(0,400);
  h_mmumu->SetTitle("DiMuon Mass;Mass (GeV);Events");
  h_mmumu->GetYaxis()->SetTitleOffset(1.5);
  h_mmumumc->Draw("hist same");
  h_mmumumc->SetFillStyle(0);
  h_mmumumc->SetLineColor(3);

  gStyle->SetOptFit(1);
  h_mmumu->SetStats(0);  // kill legend of h0

  c->Print("fitZ.pdf");
}
     
