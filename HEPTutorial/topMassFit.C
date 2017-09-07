#include <TFitResult.h>
#include <TMath.h>

void topMassFit(){

  //read root file
  TFile *_file0 = TFile::Open("hist_data.root");
  TH1F *h_LepT=(TH1F*)_file0->Get("LepTopMass_S5");
  TH1F *h_HadT=(TH1F*)_file0->Get("HadTopMass_S5");

  TFile *_file1 = TFile::Open("hist_stack.root");
  TH1F *h_LepTmc=(TH1F*)_file1->Get("LepTopMass_S5");
  TH1F *h_HadTmc=(TH1F*)_file1->Get("HadTopMass_S5");

  //create clone of hgg dist
  TH1F *h_lepT = (TH1F*)h_LepT->Clone("h_lepT");
  TH1F *h_hadT = (TH1F*)h_HadT->Clone();
  TH1F *h_lepTmc = (TH1F*)h_LepTmc->Clone("h_lepTmc");
  TH1F *h_hadTmc = (TH1F*)h_HadTmc->Clone();


  //define gaus fit fct
  TF1 *mygaus1 = new TF1("mygaus1", "[0]/(sqrt(2.0*3.14151927)*[2])*exp(-(x-[1])*(x-[1])/(2.0*[2]*[2]))", 100.0, 240.0);
    mygaus1->SetParName(0, "Normalize");
    mygaus1->SetParName(1, "Mean");
    mygaus1->SetParName(2, "Width");
    double initval1[3] = { 100, 172.0, 10};
    mygaus1->SetParameters(initval1);

  TF1 *mygaus2 = new TF1("mygaus2", "[3]/(sqrt(2.0*3.14151927)*[5])*exp(-(x-[4])*(x-[4])/(2.0*[5]*[5]))", 100.0, 240.0);
    mygaus2->SetParName(3, "Normalize");
    mygaus2->SetParName(4, "Mean");
    mygaus2->SetParName(5, "Width");
    double initval2[3] = { 100.0, 172, 20};
    mygaus2->SetParameters(initval2);

  //create canvas
  TCanvas *c1=new TCanvas("c1","Leptonic Top Mass",600,600);
  TCanvas *c2=new TCanvas("c2","Hadronic Top Mass",600,600);

  //store fit result
  TFitResultPtr r1= h_lepT->Fit(mygaus1, "S");
  Double_t par0   = r1->Parameter(0);
  Double_t par1   = r1->Parameter(1);
  Double_t par2   = r1->Parameter(2);
	double initval3[3] = {par0,par1,par2};

  TFitResultPtr r2= h_hadT->Fit(mygaus2, "S");
  Double_t par3   = r2->Parameter(0);
  Double_t par4   = r2->Parameter(1);
  Double_t par5   = r2->Parameter(2);
  double initval4[3] = {par3,par4,par5};

  //from TFitResult, we create bkg2 which will appear on our canvas
  TF1*lepT2 = new TF1("lepT2", "[8]/(sqrt(2.0*3.14151927)*[10])*exp(-(x-[9])*(x-[9])/(2.0*[10]*[10]))", 100.0, 250.0);
	lepT2->SetParameters(initval3);

  TF1*hadT2 = new TF1("hadT2", "[11]/(sqrt(2.0*3.14151927)*[13])*exp(-(x-[12])*(x-[12])/(2.0*[13]*[13]))", 200.0, 250.0);
  hadT2->SetParameters(initval4);

  //draw
  c1->cd();
  h_lepT->Draw("P");
  h_lepT->SetMarkerStyle(8);
  h_lepT->SetMarkerSize(0.6);
  h_lepT->GetYaxis()->SetRangeUser(0,30);
  h_lepT->GetXaxis()->SetRangeUser(0,400);
  h_lepT->SetTitle("Lep Top Mass;Mass (GeV);Events");
  h_lepT->GetYaxis()->SetTitleOffset(1.5); 

  h_lepTmc->Draw("hist same");
  h_lepTmc->SetFillStyle(0);
  h_lepTmc->SetLineColor(3);

  gStyle->SetOptFit(1);
  h_lepT->SetStats(0);  // kill legend of h0

  c2->cd();
  h_hadT->Draw("P");
  h_hadT->SetMarkerStyle(8);
  h_hadT->SetMarkerSize(0.6);
  h_hadT->GetYaxis()->SetRangeUser(0,120);
  h_hadT->GetXaxis()->SetRangeUser(0,400);
  h_hadT->SetTitle("Had Top Mass;Mass (GeV);Events");
  h_hadT->GetYaxis()->SetTitleOffset(1.5);

  h_hadTmc->Draw("hist same");
  h_hadTmc->SetFillStyle(0);
  h_hadTmc->SetLineColor(3);

  gStyle->SetOptFit(1);
  h_hadT->SetStats(0);  // kill legend of h0

  //Print
  c1->Print("fitLep.pdf");
  c2->Print("fitHad.pdf");
}
