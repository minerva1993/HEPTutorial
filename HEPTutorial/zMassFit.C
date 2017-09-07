#include <TFitResult.h>
#include <TMath.h>
Double_t mybw(Double_t* x, Double_t* par);

void zMassFit(){

  //read root file
  TFile *_file0 = TFile::Open("hist_data.root");
  TH1D *h_Mmumu=(TH1D*)_file0->Get("Mmumu");

  TFile *_file1 = TFile::Open("hist_stack.root");
  TH1F *h_Mmumumc=(TH1F*)_file1->Get("Mmumu");

  //create clone of hgg dist
  TH1F *h_mmumu = (TH1F*)h_Mmumu->Clone();
  TH1F *h_mmumu2 = (TH1F*)h_Mmumu->Clone();
  TH1F *h_mmumumc = (TH1F*)h_Mmumumc->Clone();

  //define gaus fit fct
  TF1 *mygaus = new TF1("mygaus", "[0]/(sqrt(2.0*3.14151927)*[2])*exp(-(x-[1])*(x-[1])/(2.0*[2]*[2]))", 0.0, 200.0);
    mygaus->SetParName(0, "Normalize");
    mygaus->SetParName(1, "Mean");
    mygaus->SetParName(2, "Width");
    double initval1[3] = { 13000, 90.0, 1};
    mygaus->SetParameters(initval1);

  TF1 *bw = new TF1("mybw", mybw, 0, 150, 3);
    bw->SetParName(0, "Constant");
    bw->SetParName(1, "Sigma");
    bw->SetParName(2, "Mean");
    double initval2[3] = { 1.0, 5.0, 90.0};
    bw->SetParameters(initval2);

  //create canvas
  TCanvas *c=new TCanvas("c","DiMuon(Z) Mass",600,600);

  //store fit result
  TFitResultPtr r= h_mmumu->Fit(mygaus, "S");
  Double_t par0   = r->Parameter(0);
  Double_t par1   = r->Parameter(1);
  Double_t par2   = r->Parameter(2);
  double initval3[3] = {par0,par1,par2};

  TF1*Mmumu = new TF1("Mmumu", "[0]/(sqrt(2.0*3.14151927)*[2])*exp(-(x-[1])*(x-[1])/(2.0*[2]*[2]))", 30.0, 150.0);
  Mmumu->SetParameters(initval3);

  TFitResultPtr r2= h_mmumu2->Fit(bw, "S");
  Double_t par3   = r2->Parameter(0);
  Double_t par4   = r2->Parameter(1);
  Double_t par5   = r2->Parameter(2);
  double initval4[3] = {par3,par4,par5};

  TF1 *Mmumu2 = new TF1("Mmumu2", mybw, 0, 150, 3);
  Mmumu2->SetParameters(initval4);

  //draw
  c->cd();
  h_mmumu->Draw("P");
  h_mmumu->SetMarkerStyle(8);
  h_mmumu->SetMarkerSize(0.6);
  h_mmumu->GetYaxis()->SetRangeUser(0,3000);
  h_mmumu->GetXaxis()->SetRangeUser(0,400);
  h_mmumu->SetTitle("DiMuon Mass;Mass (GeV);Events");
  h_mmumu->GetYaxis()->SetTitleOffset(1.5);
  h_mmumu2->Draw("p same");
  h_mmumu2->SetMarkerSize(0.);
  h_mmumumc->Draw("hist same");
  h_mmumumc->SetFillStyle(0);
  h_mmumumc->SetLineColor(3);

  gStyle->SetOptFit(1);
  h_mmumu->SetStats(0);  // kill legend of h0

  c->Print("fitZ.pdf");
}

Double_t mybw(Double_t* x, Double_t* par){
  Double_t arg1 = 14.0/22.0; // 2 over pi
  Double_t arg2 = par[1]*par[1]*par[2]*par[2]; //Gamma=par[1]  M=par[2]
  Double_t arg3 = ((x[0]*x[0]) - (par[2]*par[2]))*((x[0]*x[0]) - (par[2]*par[2]));
  Double_t arg4 = x[0]*x[0]*x[0]*x[0]*((par[1]*par[1])/(par[2]*par[2]));
  return par[0]*arg1*arg2/(arg3 + arg4);
}  
