#include "MyAnalysis.h"
#include "Plotter.h"
#include <iostream>
#include <TChain.h>
#include <TGraphAsymmErrors.h>
#include <string>

int main() {
   float lumi = 50.;
   
   MyAnalysis *A = new MyAnalysis();
   TChain* ch = new TChain("events");
   ch->Add("files/data.root");
   ch->Process(A);
   
   MyAnalysis *B = new MyAnalysis();
   TChain* ch2 = new TChain("events");
   ch2->Add("files/ttbar.root");
   ch2->Process(B);
   
   MyAnalysis *C = new MyAnalysis();
   TChain* ch3 = new TChain("events");
   ch3->Add("files/wjets.root");
   ch3->Process(C);
   
   MyAnalysis *D = new MyAnalysis();
   TChain* ch4 = new TChain("events");
   ch4->Add("files/dy.root");
   ch4->Process(D);
   
   MyAnalysis *E = new MyAnalysis();
   TChain* ch5 = new TChain("events");
   ch5->Add("files/ww.root");
   ch5->Process(E);

   MyAnalysis *F = new MyAnalysis();
   TChain* ch6 = new TChain("events");
   ch6->Add("files/wz.root");
   ch6->Process(F);

   MyAnalysis *G = new MyAnalysis();
   TChain* ch7 = new TChain("events");
   ch7->Add("files/zz.root");
   ch7->Process(G);

   MyAnalysis *H = new MyAnalysis();
   TChain* ch8 = new TChain("events");
   ch8->Add("files/qcd.root");
   ch8->Process(H);
   
   MyAnalysis *I = new MyAnalysis();
   TChain* ch9 = new TChain("events");
   ch9->Add("files/single_top.root");
   ch9->Process(I);

	Plotter P;
	P.SetData(A->histograms, std::string("Data"));
	P.AddBg(B->histograms, std::string("TTbar"));
	P.AddBg(C->histograms, std::string("Wjets"));
	P.AddBg(D->histograms, std::string("DY"));
	P.AddBg(E->histograms, std::string("WW"));
	P.AddBg(F->histograms, std::string("WZ"));
	P.AddBg(G->histograms, std::string("ZZ"));
	P.AddBg(H->histograms, std::string("QCD"));
	P.AddBg(I->histograms, std::string("single Top"));
   
	P.Plot(string("results.pdf"));
   
	Plotter P_MC;
	P_MC.AddBg(B->histograms_MC, std::string("TTbar"));
	P_MC.AddBg(C->histograms_MC, std::string("Wjets"));
	P_MC.AddBg(D->histograms_MC, std::string("DY"));
	P_MC.AddBg(E->histograms_MC, std::string("WW"));
	P_MC.AddBg(F->histograms_MC, std::string("WZ"));
	P_MC.AddBg(G->histograms_MC, std::string("ZZ"));
	P_MC.AddBg(H->histograms_MC, std::string("QCD"));
	P_MC.AddBg(I->histograms_MC, std::string("single Top"));
  P_MC.Plot(string("results_MC.pdf"));

   for ( int i = 0 ; i < A->histograms.size() ; i++)
   {
   double N_data = A->histograms[i]->Integral();
   double ttbar_mc = B->histograms[i]->Integral();
   double Wjets_mc = C->histograms[i]->Integral();
   double DY_mc = D->histograms[i]->Integral();
   double WW_mc = E->histograms[i]->Integral();
   double WZ_mc = F->histograms[i]->Integral();
   double ZZ_mc = G->histograms[i]->Integral();
   double QCD_mc = H->histograms[i]->Integral();
   double singleTop_mc = I->histograms[i]->Integral();
   double total_mc = ttbar_mc + Wjets_mc + DY_mc + WW_mc + WZ_mc + ZZ_mc + QCD_mc + singleTop_mc;
	 double total_bkg = Wjets_mc + DY_mc + WW_mc + WZ_mc + ZZ_mc + QCD_mc + singleTop_mc;
   double sig = ttbar_mc / sqrt(total_mc);
	 double cx = (N_data - total_bkg) / (50*(ttbar_mc / B->histograms[4]->Integral() ) );
   cout << "h_" << i+1 << endl;
   cout << "data : " << N_data << endl;
   cout << "ttbar : " << ttbar_mc << endl;
   cout << "Wjets : " << Wjets_mc << endl;
   cout << "DY : " << DY_mc << endl;
   cout << "WW : " << WW_mc << endl;
   cout << "WZ : " << WZ_mc << endl;
   cout << "ZZ : " << ZZ_mc << endl;
   cout << "QCD : " << QCD_mc << endl;
   cout << "singleTop : " << singleTop_mc << endl;
	 cout << "Total mc : " << total_mc << endl;
   cout << "Total Background with out ttbar_mc : " << total_bkg << endl;
   cout << "Significance : " << sig << endl;
	 cout << "Efficiency : " << ttbar_mc / B->histograms[4]->Integral() << endl;
	 cout << "Cross section : " << cx << endl;
	 cout << "error : " << cx / sqrt(N_data) << endl;
	 cout << endl; 
  }

}
