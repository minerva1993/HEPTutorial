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

  h_NMuon = new TH1F("NMuon", "Number of muons", 7, 0, 7);
  h_NMuon->SetXTitle("No. Muons");
  h_NMuon->Sumw2();
  histograms.push_back(h_NMuon);
  histograms_MC.push_back(h_NMuon);

  h_NElectron = new TH1F("NElectron", "Number of Electron", 7, 0, 7);
  h_NElectron->SetXTitle("No.Electron");
  h_NElectron->Sumw2();
  histograms.push_back(h_NElectron);
  histograms_MC.push_back(h_NElectron);

  for(int i=0 ; i < 9; i++){

    h_NJet_S[i] = new TH1F(Form("NJet_S%i_%s",i,option.Data()), "Number of jets", 12, 0, 12);
    h_NJet_S[i]->SetXTitle(Form("Jet Multplicity (S%i)",i));
    h_NJet_S[i]->Sumw2();
    histograms.push_back(h_NJet_S[i]);
    histograms_MC.push_back(h_NJet_S[i]);

    h_NBJet_S[i] = new TH1F(Form("NBJet_S%i_%s",i,option.Data()), "Number of b jets", 5, 0, 5);
    h_NBJet_S[i]->SetXTitle("b-tagged Jet Multiplicity");
    h_NBJet_S[i]->Sumw2();
    histograms.push_back(h_NBJet_S[i]);
    histograms_MC.push_back(h_NBJet_S[i]);

    h_metPt[i] = new TH1F(Form("metPt_S%i_%s",i,option.Data()), "MET Pt", 40, 0, 200);
    h_metPt[i]->SetXTitle("MET Pt (GeV)");
    h_metPt[i]->Sumw2();
    histograms.push_back(h_metPt[i]);
    histograms_MC.push_back(h_metPt[i]);

    h_LepWMass[i] = new TH1F(Form("LepWMass_S%i_%s",i,option.Data()), "W Mass (Lep Top)", 29 ,10 ,300);
    h_LepWMass[i]->SetXTitle("W Mass (Lep Top) (GeV)");
    h_LepWMass[i]->Sumw2();
    histograms.push_back(h_LepWMass[i]);
    histograms_MC.push_back(h_LepWMass[i]);

    h_LepTopMass[i] = new TH1F(Form("LepTopMass_S%i_%s",i,option.Data()), "Top Mass (Leptonic)", 23 ,40 ,500);
    h_LepTopMass[i]->SetXTitle("Top Mass (Leptonic) (GeV)");
    h_LepTopMass[i]->Sumw2();
    histograms.push_back(h_LepTopMass[i]);
    histograms_MC.push_back(h_LepTopMass[i]);

    h_HadWMass[i] = new TH1F(Form("HAdWMass_S%i_%s",i,option.Data()), "W Mass (Had Top)", 29 ,10 ,300);
    h_HadWMass[i]->SetXTitle("W Mass (Had Top) (GeV)");
    h_HadWMass[i]->Sumw2();
    histograms.push_back(h_HadWMass[i]);
    histograms_MC.push_back(h_HadWMass[i]);

    h_HadTopMass[i] = new TH1F(Form("HadTopMass_S%i_%s",i,option.Data()), "Top Mass (Hadronic)", 23 ,40 ,500);
    h_HadTopMass[i]->SetXTitle("Top Mass (Hadronic) (GeV)");
    h_HadTopMass[i]->Sumw2();
    histograms.push_back(h_HadTopMass[i]);
    histograms_MC.push_back(h_HadTopMass[i]);

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

  TLorentzVector MET;
  MET.SetXYZM(MET_px, MET_py, 0, 0);

/////Muon   
  int N_IsoMuon = 0;
  MyMuon *singleMu;

  for (vector<MyMuon>::iterator muon = Muons.begin(); muon != Muons.end(); ++muon){
    if (muon->IsIsolated(MuonRelIsoCut) && abs(muon->Eta()) < 2.4){
      ++N_IsoMuon;
      if (N_IsoMuon == 1) singleMu = &(*muon);
    }
  }

  h_NMuon->Fill(N_IsoMuon, EventWeight);

  ////Electron
  int N_IsoElectron = 0 ;
  MyElectron *singleEl;

  for (vector<MyElectron>::iterator elec = Electrons.begin(); elec != Electrons.end(); ++elec) {
    if (elec->IsIsolated(ElectronRelIsoCut) && abs(elec->Eta()) < 2.4) {
      ++N_IsoElectron;
      if (N_IsoElectron == 1) singleEl = &(*elec);
    }
  }

  h_NElectron->Fill(N_IsoElectron, EventWeight);

  ////Jets
  int N_Jets = 0;
  int N_BJets = 0;

  vector<TLorentzVector> v_jets;
  vector<TLorentzVector> v_bjets;

  for (vector<MyJet>::iterator jt = Jets.begin(); jt != Jets.end(); ++jt) {
    if (abs(jt->Eta()) < 2.4 && (jt->Pt()) > JetPtCut){
       ++N_Jets;
      TLorentzVector tmpjet;
      tmpjet.SetPtEtaPhiE(jt->Pt(), jt->Eta(), jt->Phi(), jt->E());
      v_jets.push_back(tmpjet);
    
      if (jt->IsBTagged()){
        ++N_BJets;
        v_bjets.push_back(tmpjet);
      }
    }
  }

  TLorentzVector Jets[N_Jets];
  TLorentzVector BJets[N_BJets];

  for(int i = 0; i < N_Jets; ++i){
    Jets[i].SetPxPyPzE(Jet_Px[i], Jet_Py[i], Jet_Pz[i], Jet_E[i]);
  }

  TLorentzVector lepton;
  double transverseM;

  h_NJet_S[0]->Fill(N_Jets, EventWeight);
  h_NBJet_S[0]->Fill(N_BJets, EventWeight);
  h_metPt[0]->Fill(MET.Pt(),EventWeight);


////S1
  if (triggerIsoMu24){

    h_NJet_S[1]->Fill(N_Jets, EventWeight);
    h_NBJet_S[1]->Fill(N_BJets, EventWeight);

      ////S2
    if (N_IsoElectron == 1 || N_IsoMuon == 1){
      if (N_IsoElectron == 1) lepton = *singleMu;
      if (N_IsoMuon == 1) lepton = *singleEl;
      transverseM = transverseMass(lepton, MET);

      double tmplepbDR = 99;
      double lepbDR = 99;
      double leptop;

      TLorentzVector lepb;

      if (N_BJets >= 1){
        for(int i = 0; i < N_BJets; i++){
          lepb = v_bjets[i];
          tmplepbDR = lepton.DeltaR(lepb);
          
          if(tmplepbDR < lepbDR){
            lepbDR = tmplepbDR;
            leptop = (lepton + MET + lepb).M();
          }
        }
      }

      h_NJet_S[2]->Fill(N_Jets, EventWeight);
      h_NBJet_S[2]->Fill(N_BJets, EventWeight);
      h_LepWMass[2]->Fill(transverseM, EventWeight);
      h_LepTopMass[2]->Fill(leptop, EventWeight);

      ////S3
      if (N_Jets >= 4) {
        h_NJet_S[3]->Fill(N_Jets, EventWeight);
        h_NBJet_S[3]->Fill(N_BJets, EventWeight);
        h_LepWMass[3]->Fill(transverseM, EventWeight);
        h_LepTopMass[3]->Fill(leptop, EventWeight);
      }

      ////S4
      if (N_BJets >= 1) {
        h_NJet_S[4]->Fill(N_Jets, EventWeight);
        h_NBJet_S[4]->Fill(N_BJets, EventWeight);
        h_LepWMass[4]->Fill(transverseM, EventWeight);
        h_LepTopMass[4]->Fill(leptop, EventWeight);
      }

      ////S5
      if (N_Jets >= 4 &&  N_BJets >= 1) {
        h_NJet_S[5]->Fill(N_Jets, EventWeight);
        h_NBJet_S[5]->Fill(N_BJets, EventWeight);
        h_LepWMass[5]->Fill(transverseM, EventWeight);
        h_LepTopMass[5]->Fill(leptop, EventWeight);
      }

      ////S6
      if (N_BJets >= 2){
        h_NJet_S[6]->Fill(N_Jets, EventWeight);
        h_NBJet_S[6]->Fill(N_BJets, EventWeight);
        h_LepWMass[6]->Fill(transverseM, EventWeight);
        h_LepTopMass[6]->Fill(leptop, EventWeight);
      }


      ////S7
      if (N_Jets >=4 && N_BJets >= 2){
        h_NJet_S[7]->Fill(N_Jets, EventWeight);
        h_NBJet_S[7]->Fill(N_BJets, EventWeight);
        h_LepWMass[7]->Fill(transverseM, EventWeight);
        h_LepTopMass[7]->Fill(leptop, EventWeight);
      }
/////
    }//pass lepton
  }//pass triger
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

double MyAnalysis::transverseMass( const TLorentzVector & lepton, const TLorentzVector & met){

  TLorentzVector leptonT(lepton.Px(),lepton.Py(),0.,lepton.E()*TMath::Sin(lepton.Theta()));
  TLorentzVector metT(met.Px(), met.Py(), 0, met.E());

  TLorentzVector sumT=leptonT+metT;
  double out = TMath::Sqrt( sumT.M2() );

  return out;

}

