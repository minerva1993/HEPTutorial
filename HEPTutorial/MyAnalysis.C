#define MyAnalysis_cxx

#include "MyAnalysis.h"
#include <iostream>
#include <TH1F.h>
#include <TLatex.h>

using namespace std;

void MyAnalysis::BuildEvent() {

	Muons.clear();
	for (int i = 0; i < NMuon; ++i) {
	MyMuon muon(Muon_Px[i], Muon_Py[i], Muon_Pz[i], Muon_E[i]);
	muon.SetIsolation(Muon_Iso[i]);
	muon.SetCharge(Muon_Charge[i]);
	Muons.push_back(muon);
	}

	Electrons.clear();
	for (int i = 0; i < NElectron; ++i) {
	MyElectron electron(Electron_Px[i], Electron_Py[i], Electron_Pz[i], Electron_E[i]);
	electron.SetIsolation(Electron_Iso[i]);
	electron.SetCharge(Electron_Charge[i]);
	Electrons.push_back(electron);
	}

	Photons.clear();
	for (int i = 0; i < NPhoton; ++i) {
	MyPhoton photon(Photon_Px[i], Photon_Py[i], Photon_Pz[i], Photon_E[i]);
	photon.SetIsolation(Photon_Iso[i]);
	Photons.push_back(photon);
	}

	Jets.clear();
	for (int i = 0; i < NJet; ++i) {
	MyJet jet(Jet_Px[i], Jet_Py[i], Jet_Pz[i], Jet_E[i]);
	jet.SetBTagDiscriminator(Jet_btag[i]);
	jet.SetJetID(Jet_ID[i]);
	Jets.push_back(jet);
	}

	hadB.SetXYZM(MChadronicBottom_px, MChadronicBottom_py, MChadronicBottom_pz, 4.8);
	lepB.SetXYZM(MCleptonicBottom_px, MCleptonicBottom_py, MCleptonicBottom_pz, 4.8);
	hadWq.SetXYZM(MChadronicWDecayQuark_px, MChadronicWDecayQuark_py, MChadronicWDecayQuark_pz, 0.0);
	hadWqb.SetXYZM(MChadronicWDecayQuarkBar_px, MChadronicWDecayQuarkBar_py, MChadronicWDecayQuarkBar_pz, 0.0);
	lepWl.SetXYZM(MClepton_px, MClepton_py, MClepton_pz, 0.0);
	lepWn.SetXYZM(MCneutrino_px, MCneutrino_py, MCneutrino_pz, 0.0);
	met.SetXYZM(MET_px, MET_py, 0., 0.);

	EventWeight *= weight_factor;
   
}

void MyAnalysis::Begin(TTree * /*tree*/) {
   
  TString option = GetOption();
   
}

void MyAnalysis::SlaveBegin(TTree * /*tree*/) {
   
	TString option = GetOption();

	h_Mmumu = new TH1F("Mmumu", "Invariant di-muon mass", 60, 60, 120);
	h_Mmumu->SetXTitle("m_{#mu#mu}");
	h_Mmumu->Sumw2();
	histograms.push_back(h_Mmumu);
	histograms_MC.push_back(h_Mmumu);

	h_NMuon = new TH1F("NMuon", "Number of muons", 7, 0, 7);
	h_NMuon->SetXTitle("No. Muons");
	h_NMuon->Sumw2();
	histograms.push_back(h_NMuon);
	histograms_MC.push_back(h_NMuon);

	h_Mee = new TH1F("Mee", "eletcron mass", 60, 60, 120);
	h_Mee->SetXTitle("m_{ee}");
	h_Mee->Sumw2();
	histograms.push_back(h_Mee);
	histograms_MC.push_back(h_Mee);

	h_NElectron = new TH1F("NElectron", "Number of Electron", 7, 0, 7);
	h_NElectron->SetXTitle("No.Electron");
	h_NElectron->Sumw2();
	histograms.push_back(h_NElectron);
	histograms_MC.push_back(h_NElectron);

	for(int i=0 ; i < 9; i++){

  h_NJet_S[i] = new TH1F(Form("NJet_S%i_%s",i,option.Data()), "Number of jets", 12, 0, 12);
  h_NJet_S[i]->SetXTitle("Jet Multplicity");
  h_NJet_S[i]->Sumw2();
  histograms.push_back(h_NJet_S[i]);
  histograms_MC.push_back(h_NJet_S[i]);

  h_JetMass_S[i] = new TH1F(Form("JetMass_S%i_%s",i,option.Data()), "Mass of jets", 60, 60, 120);
  h_JetMass_S[i]->SetXTitle("m_{jj}");
  h_JetMass_S[i]->Sumw2();
  histograms.push_back(h_JetMass_S[i]);
  histograms_MC.push_back(h_JetMass_S[i]);

  h_NBJet_S[i] = new TH1F(Form("NBJet_S%i_%s",i,option.Data()), "Number of b jets", 5, 0, 5);
  h_NBJet_S[i]->SetXTitle("b-tagged Jet Multiplicity");
  h_NBJet_S[i]->Sumw2();
  histograms.push_back(h_NBJet_S[i]);
  histograms_MC.push_back(h_NBJet_S[i]);
	}
}

Bool_t MyAnalysis::Process(Long64_t entry) {
  
   ++TotalEvents;
   
   GetEntry(entry);
   
   if (TotalEvents % 10000 == 0)
      cout << "Next event -----> " << TotalEvents << endl;
   
   BuildEvent();
   
   double MuonPtCut = 30.;
   double MuonRelIsoCut = 0.10;
   double ElectronPtCut = 20;
   double ElectronRelIsoCut = 0.3;
   double JetPtCut = 35.;
   double JetRelIsoCut = 0;

/////Muon   
	 int N_IsoMuon = 0;
	 MyMuon *muon1, *muon2;

	 for (vector<MyMuon>::iterator muon = Muons.begin(); muon != Muons.end(); ++muon) {
		 if (muon->IsIsolated(MuonRelIsoCut) && abs(muon->Eta()) < 2.4) {
			 ++N_IsoMuon;
			 if (N_IsoMuon == 1) muon1 = &(*muon);
			 if (N_IsoMuon == 2) muon2 = &(*muon);
		 }
	 }

	 h_NMuon->Fill(N_IsoMuon, EventWeight);

	 if (N_IsoMuon > 1 && triggerIsoMu24) {
		 if (muon1->Pt()>MuonPtCut) {
			 h_Mmumu->Fill((*muon1 + *muon2).M(), EventWeight);
		 }
	 } 

////Electron
	 int N_IsoElectron = 0 ;
	 MyElectron *electron1, *electron2;

	 for (vector<MyElectron>::iterator elec = Electrons.begin(); elec != Electrons.end(); ++elec) {
		 if (elec->IsIsolated(ElectronRelIsoCut) && abs(elec->Eta()) < 2.4) {
			 ++N_IsoElectron;
			 if (N_IsoElectron == 1) electron1 = &(*elec);
			 if (N_IsoElectron == 2) electron2 = &(*elec);
		 }
	 }

	 h_NElectron->Fill(N_IsoElectron, EventWeight);

	 if (N_IsoElectron > 1 && triggerIsoMu24) {
		 if (electron1->Pt()>ElectronPtCut) {
			 h_Mee->Fill((*electron1 + *electron2).M(), EventWeight);
		 }
	 }

////Jets
	 MyJet *jet_1, *jet_2, *jet_3, *jet_4;
	 MyJet *bjet_1, *bjet_2;

		int N_Jets = 0;
		int N_BJets = 0;
		for (vector<MyJet>::iterator jt = Jets.begin(); jt != Jets.end(); ++jt) {
			if (abs(jt->Eta()) < 2.4 && (jt->Pt()) > JetPtCut) ++N_Jets;
			if (N_Jets == 1) jet_1 = &(*jt);
      if (N_Jets == 2) jet_2 = &(*jt);
      if (N_Jets == 3) jet_3 = &(*jt);
      if (N_Jets == 4) jet_4 = &(*jt);

			if(jt->IsBTagged()){ ++N_BJets;
				if (N_BJets == 1) bjet_1 = &(*jt);
				if (N_BJets == 2) bjet_2 = &(*jt);
			}
		}

			h_NJet_S[0]->Fill(N_Jets, EventWeight);
			h_NBJet_S[0]->Fill(N_BJets, EventWeight);

			if (N_Jets > 1){
      	if (jet_1->Pt()>JetPtCut){
        	h_JetMass_S[0]->Fill((*jet_1 + *jet_2).M(),EventWeight);
      	}
			}

////S1
if (triggerIsoMu24){

      h_NJet_S[1]->Fill(N_Jets, EventWeight);
      h_NBJet_S[1]->Fill(N_BJets, EventWeight);

      if (N_Jets > 1){
        if (jet_1->Pt()>JetPtCut){
          h_JetMass_S[1]->Fill((*jet_1 + *jet_2).M(),EventWeight);
        }
      }
////S2
 
	 	if(N_IsoMuon == 1){
			h_NJet_S[2]->Fill(N_Jets, EventWeight);
    	h_NBJet_S[2]->Fill(N_BJets, EventWeight);

			if (N_Jets > 1){
				if (jet_1->Pt()>JetPtCut){
					h_JetMass_S[2]->Fill((*jet_1 + *jet_2).M(),EventWeight);
				}
			}
		}

////S3
    if (N_Jets >= 4 && N_IsoMuon == 1) {
      h_NJet_S[3]->Fill(N_Jets, EventWeight);
      h_NBJet_S[3]->Fill(N_BJets, EventWeight);

      if (N_Jets >= 4){
        if (jet_1->Pt()>JetPtCut){
          h_JetMass_S[3]->Fill((*jet_1 + *jet_2).M(),EventWeight);
        }
      }
    }

////S4
    if (N_IsoMuon == 1 &&  N_BJets >= 1) {
      h_NJet_S[4]->Fill(N_Jets, EventWeight);
      h_NBJet_S[4]->Fill(N_BJets, EventWeight);

      if (N_Jets > 1){
        if (jet_1->Pt()>JetPtCut){
          h_JetMass_S[4]->Fill((*jet_1 + *jet_2).M(),EventWeight);
        }
      }
    }

////S5
		if (N_Jets >= 4 && N_IsoMuon == 1 &&  N_BJets >= 1) {
			h_NJet_S[5]->Fill(N_Jets, EventWeight);
			h_NBJet_S[5]->Fill(N_BJets, EventWeight);

			if (N_Jets >= 4){
      	if (jet_1->Pt()>JetPtCut){
        	h_JetMass_S[5]->Fill((*jet_1 + *jet_2).M(),EventWeight);
      	}
			}
		}

////S6
		if (N_IsoMuon == 1 && N_BJets >= 2){
			h_NJet_S[6]->Fill(N_Jets, EventWeight);
   		h_NBJet_S[6]->Fill(N_BJets, EventWeight);

			if (N_Jets > 1){
       	if (jet_1->Pt()>JetPtCut){
         	h_JetMass_S[6]->Fill((*jet_1 + *jet_2).M(),EventWeight);
       	}
     	}
   	}

////S7
    if (N_Jets >=4 && N_IsoMuon == 1 && N_BJets >= 2){
      h_NJet_S[7]->Fill(N_Jets, EventWeight);
      h_NBJet_S[7]->Fill(N_BJets, EventWeight);

      if (N_Jets >= 4){
        if (jet_1->Pt()>JetPtCut){
          h_JetMass_S[7]->Fill((*jet_1 + *jet_2).M(),EventWeight);
        }
      }
    }

///S8
  	if(N_IsoMuon == 2){
			h_NJet_S[8]->Fill(N_Jets, EventWeight);
    	h_NBJet_S[8]->Fill(N_BJets, EventWeight);

    	if (N_Jets > 1){
      	if (jet_1->Pt()>JetPtCut){
        	h_JetMass_S[8]->Fill((*jet_1 + *jet_2).M(),EventWeight);
      	}
    	}
		}
/////
	}
//////////////////////////////
   
   return kTRUE;
}

void MyAnalysis::SlaveTerminate() {
   // The SlaveTerminate() function is called after all entries or objects
   // have been processed. When running with PROOF SlaveTerminate() is called
   // on each slave server.

  
}

void MyAnalysis::Terminate() {
   // The Terminate() function is the last function to be called during
   // a query. It always runs on the client, it can be used to present
   // the results graphically or save the results to file.
   
}
