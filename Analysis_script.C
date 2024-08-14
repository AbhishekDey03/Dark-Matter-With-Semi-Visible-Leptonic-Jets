/*
 * Analysis Code for isolating SVJl Signal.
 * 
 * Author: Abhishek Ranjan Dey
 * Created: 14/08/2023
 * 
 * Dependencies:
 * - ROOT (TLorentzVector, TFile, TTree, TH1D, TH1F, TStyle, TVector2, TVector3)
 * - FastJet (ClusterSequence, PseudoJet)
 * 
 * Example Compile Line:
 * g++ `root-config --cflags` `root-config --libs` -o Analysisscript Analysis_script.C `fastjet-install/bin/fastjet-config --cxxflags --libs --plugins`
 * 
 * The data files are not provided and therefore the analysis code requires you change the input and output file paths
 */


#include "TLorentzVector.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TVector2.h"
#include "TVector3.h"
#include "TH1F.h"
#include "TStyle.h"
#include <vector>
#include <iostream>
#include <set>
#include "fastjet/ClusterSequence.hh"

using namespace fastjet;
using namespace std;

#ifdef MAKECINT
#pragma link C++ class vector+;
#endif

float_t calculate_relative_isolation(float_t topoetcone20, float_t mu_pt) {
    return fabs(topoetcone20 / mu_pt);//relative isolation is calculated using etcone20
}

void setup_histogram(TH1F* hist, const char* xTitle, const char* yTitle) {
     //parameters that are set for each hsitorram
    hist->GetXaxis()->SetTitle(xTitle);
    hist->GetYaxis()->SetTitle(yTitle);
    hist->SetLineWidth(2);
}

int main(int argc, char *argv[]) {
    // Set the gStyle
    gStyle->SetLineWidth(0.5); 
    gStyle->SetTitleFont(42, "XYZ"); 
    gStyle->SetLabelFont(42, "XYZ"); 
    gStyle->SetTitleSize(0.05, "XYZ"); 
    gStyle->SetLabelSize(0.04, "XYZ"); 
    gStyle->SetHistLineWidth(0.5);

    // Retrieve the TTree for the input data
    TFile *inputFile = new TFile("/path/to/input/data.root", "READ");
    TTree *inputTree = (TTree*)inputFile->Get("nominal");

    //set the output file.
    TFile *outputFile = new TFile("/path/to/output/data.root", "RECREATE");
    outputFile->cd();

    // Set up the addresses for reading out the branches
    UInt_t runNumber;
    ULong64_t eventNumber;
    Float_t weight_jvt;
    Float_t weight_pileup;
    Float_t weight_mc;
    Float_t weight_leptonSF;
    Float_t totalEventsWeighted = 4000; //4000 is used for the signal data
    Float_t EventWT = 1;//initialise the event weight.
    Float_t met_met;
    Float_t met_phi;

    //initialise observable containers.
    vector<float> *jet_pt = 0;
    vector<float> *jet_phi = 0;
    vector<float> *jet_eta = 0;
    vector<float> *jet_e = 0;
    vector<char> *jet_isbtagged_DL1r_77 = 0;
    vector<float> *mu_pt = 0;
    vector<float> *mu_eta = 0;
    vector<float> *mu_phi = 0;
    vector<float> *mu_topoetcone20 = 0;
    vector<float> *mu_charge = 0;
    vector<float> *mu_e = 0;

    // Associate the addresses with the TTree
    inputTree->SetBranchAddress("runNumber", &runNumber);
    inputTree->SetBranchAddress("eventNumber", &eventNumber);
    inputTree->SetBranchAddress("weight_mc", &weight_mc);
    inputTree->SetBranchAddress("weight_jvt", &weight_jvt);
    inputTree->SetBranchAddress("weight_leptonSF", &weight_leptonSF);
    inputTree->SetBranchAddress("met_met", &met_met);
    inputTree->SetBranchAddress("met_phi", &met_phi);
    inputTree->SetBranchAddress("weight_pileup", &weight_pileup);

    inputTree->SetBranchAddress("jet_pt", &jet_pt);
    inputTree->SetBranchAddress("jet_phi", &jet_phi);
    inputTree->SetBranchAddress("jet_eta", &jet_eta);
    inputTree->SetBranchAddress("jet_e", &jet_e);

    inputTree->SetBranchAddress("mu_pt", &mu_pt);
    inputTree->SetBranchAddress("mu_eta", &mu_eta);
    inputTree->SetBranchAddress("mu_phi", &mu_phi);
    inputTree->SetBranchAddress("mu_charge", &mu_charge);
    inputTree->SetBranchAddress("mu_topoetcone20", &mu_topoetcone20);
    inputTree->SetBranchAddress("mu_e", &mu_e);

    // Declare CutFlow Histogram
    TH1F *CutFlow_sel = new TH1F("CutFlow_sel", "cutflow selections", 13, -0.5, 12.5);
    setup_histogram(CutFlow_sel, "Cutflow", "Events");

    // Declare histograms for cutflow steps 3, 7, and 11
    vector<TH1F*> histograms[3];

    //iterate over these and change the names using a simple iterative variable
    for (int step = 0; step < 3; ++step) {
        histograms[step].push_back(new TH1F(Form("HT_cutflow_%d", step), "HT", 40, 0, 4000));
        histograms[step].push_back(new TH1F(Form("MET_cutflow_%d", step), "MET", 40, 0, 3500));
        histograms[step].push_back(new TH1F(Form("RT_cutflow_%d", step), "RT", 30, 0, 3.0));
        histograms[step].push_back(new TH1F(Form("MT_cutflow_%d", step), "MT (dijet + met)", 35, 0, 3500));
        histograms[step].push_back(new TH1F(Form("mjj_cutflow_%d", step), "mjj", 35, 0, 3500));
        histograms[step].push_back(new TH1F(Form("mll_cutflow_%d", step), "mll", 50, 0, 25));
        histograms[step].push_back(new TH1F(Form("leading_jet_pt_cutflow_%d", step), "Leading Jet pT", 15, 0, 1500));
        histograms[step].push_back(new TH1F(Form("subleading_jet_pt_cutflow_%d", step), "Subleading Jet pT", 15, 0, 1500));
        histograms[step].push_back(new TH1F(Form("MINphi_cutflow_%d", step), "MINphi (closest jet, MET)", 35, 0, 3.5));
        histograms[step].push_back(new TH1F(Form("MAXphi_cutflow_%d", step), "MAXphi (farthest jet, MET)", 35, 0, 3.5));
        histograms[step].push_back(new TH1F(Form("MAX_MINphi_cutflow_%d", step), "MAX-MINphi", 35, 0, 3.5));
        histograms[step].push_back(new TH1F(Form("inter_iso_cutflow_%d", step), "Inter-iso", 10, 0, 1));
        histograms[step].push_back(new TH1F(Form("rel_iso_cutflow_%d", step), "Rel-iso", 10, 0, 1));
        histograms[step].push_back(new TH1F(Form("n_muons_non_iso_cutflow_%d", step), "Number of non-iso muons", 10, .5, 10.5));

        for (auto hist : histograms[step]) {
            hist->SetDirectory(outputFile);
        }
    }

    // Setup histogram properties
    for (auto& step_histograms : histograms) {
        for (auto hist : step_histograms) {
            setup_histogram(hist, hist->GetName(), "Events");
        }
    }

    // Start the event loop
    for (int eventIndex = 0; eventIndex < inputTree->GetEntries(); eventIndex++) {
        // Read the event info
        inputTree->GetEvent(eventIndex);

        EventWT = (totalEventsWeighted * weight_mc * weight_pileup * weight_jvt * weight_leptonSF);

        // Fill CutFlow for all events before cuts
        CutFlow_sel->Fill(1, EventWT); // CutFlow bin 1

        // Create booleans for each cut
        //largely used to avoid mismanaging counting events in loops
        bool cutflow_2 = false;
        bool cutflow_3 = false;
        bool cutflow_4 = false;
        bool cutflow_5 = false;
        bool cutflow_6 = false;
        bool cutflow_7 = false;
        bool cutflow_8 = false;
        bool cutflow_9 = false;
        bool cutflow_10 = false;
        bool cutflow_11 = false;

        // Initialise muon containers and non-isolated muon indices
        vector<TLorentzVector> allMuons;
        vector<size_t> nonIsoMuonsIndices;

        //recluster jets into a fastjet::PseudoJet container
        vector<PseudoJet> myjets;
        int smalljet_size = (*jet_pt).size();
        for (unsigned int j = 0; j < smalljet_size; j++) {
            if (fabs((*jet_eta)[j]) < 2.8 && ((*jet_pt)[j] / 1000) > 30) {
                Double_t pt = (*jet_pt)[j];
                Double_t eta = (*jet_eta)[j];
                Double_t phi = (*jet_phi)[j];
                Double_t energy = (*jet_e)[j];
                Double_t px = pt * cos(phi);
                Double_t py = pt * sin(phi);
                Double_t pz = pt * sinh(eta);
                myjets.emplace_back(PseudoJet(px, py, pz, energy));
            }
        }

        //define jets with R=1.0 and use the anti-k_t clustering algorithm.
        Double_t R = 1.0;
        JetDefinition jet_def(antikt_algorithm, R);
        ClusterSequence cs(myjets, jet_def);
        vector<PseudoJet> newradiusjets = sorted_by_pt(cs.inclusive_jets());
        
        int inclusiveJets_size = newradiusjets.size();
        
        // Convert FastJet jets to TLorentzVector for use throughout with other TLorentzVector objexts
        vector<TLorentzVector> selectedJets;
        for (auto &jet : newradiusjets) {
            TLorentzVector lv;
            lv.SetPtEtaPhiM(jet.pt(), jet.eta(), jet.phi(), jet.m());
            selectedJets.push_back(lv);
        }
        
        // Process muons and store them in TLorentzVector format
        for (size_t i = 0; i < mu_pt->size(); ++i) {
            TLorentzVector muon;
            muon.SetPtEtaPhiE((*mu_pt)[i], (*mu_eta)[i], (*mu_phi)[i], (*mu_e)[i]);
            allMuons.push_back(muon);
        }

        // Check jet multiplicity
        if (inclusiveJets_size > 1) {
            cutflow_2 = true;
            // Fill CutFlow for jets >= 2
            CutFlow_sel->Fill(2, EventWT);

            // Check leading and sub-leading jet pT and |η|
            if (selectedJets[0].Pt() > 200e3 && fabs(selectedJets[0].Eta()) < 2.4 &&
                selectedJets[1].Pt() > 200e3 && fabs(selectedJets[1].Eta()) < 2.4) {
                cutflow_3 = true;

                double leadPt = selectedJets[0].Pt();
                double leadEta = selectedJets[0].Eta();
                double subLeadPt = selectedJets[1].Pt();
                double subLeadEta = selectedJets[1].Eta();

                float scalarSumJetPt = 0;   
                for (const auto& jet : selectedJets) {
                    scalarSumJetPt += jet.Pt();
                }

                // Fill CutFlow for pT and η cut
                CutFlow_sel->Fill(3, EventWT); // CutFlow bin 3

                TLorentzVector metVector;
                metVector.SetPtEtaPhiM(met_met, 0.0, met_phi, 0.0);

                float dphiLeadMet = fabs(selectedJets[0].Phi() - metVector.Phi());
                float dphiSubLeadMet = fabs(selectedJets[1].Phi() - metVector.Phi());
                float dphiMin = 999.0;
                float dphiMax = -999.0;

                //find PHI_min/max.(Met_jet)
                for (const auto& jet : selectedJets) {
                    float dphiCheck = fabs(TVector2::Phi_mpi_pi(jet.Phi() - metVector.Phi()));
                    if (dphiCheck < dphiMin) {
                        dphiMin = dphiCheck;
                    }
                    if (dphiCheck > dphiMax) {
                        dphiMax = dphiCheck;
                    }
                }

                TLorentzVector jet1 = selectedJets[0];
                TLorentzVector jet2 = selectedJets[1];
                
                double dijetEt = sqrt(pow((jet1 + jet2).M(), 2) + pow((jet1 + jet2).Pt(), 2));
                double MT2 = sqrt(pow(dijetEt + met_met, 2) - pow((jet1 + jet2 + metVector).Pt(), 2));

                histograms[0][0]->Fill(scalarSumJetPt / 1e3, EventWT); // HT cutflow step 3
                histograms[0][1]->Fill(met_met / 1e3, EventWT); // MET cutflow step 3
                histograms[0][2]->Fill(met_met / MT2, EventWT); // RT cutflow step 3
                histograms[0][3]->Fill(MT2 / 1e3, EventWT); // MT (dijet + met) cutflow step 3
                histograms[0][4]->Fill((jet1 + jet2).M() / 1e3, EventWT); // mjj cutflow step 3
                // mll cutflow step 3 not filled here
                histograms[0][6]->Fill(leadPt / 1e3, EventWT); // Leading Jet pT cutflow step 3
                histograms[0][7]->Fill(subLeadPt / 1e3, EventWT); // Subleading Jet pT cutflow step 3
                histograms[0][8]->Fill(dphiMin, EventWT); // MINphi cutflow step 3
                histograms[0][9]->Fill(dphiMax, EventWT); // MAXphi cutflow step 3
                histograms[0][10]->Fill(fabs(dphiMin - dphiMax), EventWT); // MAX-MINphi cutflow step 3
                //TRIGGER
                if (met_met > 200e3) {
                    cutflow_4 = true;
                    // Fill CutFlow for Δφ(j_any, MET) < 2.0
                    CutFlow_sel->Fill(4, EventWT); // CutFlow bin 4

                    if (dphiLeadMet < 2.0 || dphiSubLeadMet < 2.0) {
                        cutflow_5 = true;

                        CutFlow_sel->Fill(5, EventWT); // CutFlow bin 5

                        if ((met_met / MT2) > 0.1) {
                            cutflow_6 = true;
                            // Fill CutFlow for dijet R_T > 0.1
                            CutFlow_sel->Fill(6, EventWT); // CutFlow bin 6

                            // Process muons for isolation checks
                            int numMuonsNonIso = 0;
                            vector<size_t> isolatedMuonsIndices;

                            for (size_t i = 0; i < allMuons.size(); ++i) {
                                TLorentzVector muon = allMuons[i];

                                // Check if muon is within jet radius
                                bool inJet = false;
                                double_t R=1.;
                                for (const auto& jet : selectedJets) {
                                    if (muon.DeltaR(jet) < R) {
                                        inJet = true;
                                        break;
                                    }
                                }

                                if (inJet && i < mu_topoetcone20->size()) {
                                    float relIso = calculate_relative_isolation((*mu_topoetcone20)[i], muon.Pt());
                                    histograms[1][12]->Fill(relIso, EventWT); // Rel-iso cutflow step 7

                                    if (relIso >= 0.1) {
                                        nonIsoMuonsIndices.push_back(i);
                                        numMuonsNonIso++;
                                    } else {
                                        isolatedMuonsIndices.push_back(i);
                                    }
                                }
                            }

                            // Fill CutFlow for at least one non-isolated muon
                            
                            bool  cutflow_10 = false;
                            if (allMuons.size() >= 2) {
                                cutflow_7 = true;
                                // Fill CutFlow for number of muons >= 2
                                CutFlow_sel->Fill(7, EventWT); // CutFlow bin 7
                                if (numMuonsNonIso >= 1) { //check if there are any muons that pass the rel iso check
                                    cutflow_8 = true;
                                    CutFlow_sel->Fill(8, EventWT); // CutFlow bin 8
                            }

                                histograms[1][0]->Fill(scalarSumJetPt / 1e3, EventWT); // HT cutflow step 7
                                histograms[1][1]->Fill(met_met / 1e3, EventWT); // MET cutflow step 7
                                histograms[1][2]->Fill(met_met / MT2, EventWT); // RT cutflow step 7
                                histograms[1][3]->Fill(MT2 / 1e3, EventWT); // MT (dijet + met) cutflow step 7
                                histograms[1][4]->Fill((jet1 + jet2).M() / 1e3, EventWT); // mjj cutflow step 7
                                // mll cutflow step 7 not filled here
                                histograms[1][6]->Fill(selectedJets[0].Pt() / 1e3, EventWT); // Leading Jet pT cutflow step 7
                                histograms[1][7]->Fill(selectedJets[1].Pt() / 1e3, EventWT); // Subleading Jet pT cutflow step 7
                                histograms[1][8]->Fill(dphiMin, EventWT); // MINphi cutflow step 7
                                histograms[1][9]->Fill(dphiMax, EventWT); // MAXphi cutflow step 7
                                histograms[1][10]->Fill(fabs(dphiMin - dphiMax), EventWT); // MAX-MINphi cutflow step 7
                                histograms[1][13]->Fill(numMuonsNonIso, EventWT); // Number of non-iso muons cutflow step 7

                                if (numMuonsNonIso >= 2) {
                                    cutflow_9 = true;
                                    // Fill CutFlow for rel-iso <= 0.3 and num-iso >= 2
                                    CutFlow_sel->Fill(9, EventWT); // CutFlow bin 9

                                    vector<pair<size_t, size_t>> nonIsolatedMuonPairs;
                                    set<size_t> muonsAfterChargeCutSet;

                                    for (size_t i = 0; i < nonIsoMuonsIndices.size(); ++i) {
                                        for (size_t j = i + 1; j < nonIsoMuonsIndices.size(); ++j) {
                                            size_t index_i = nonIsoMuonsIndices[i];
                                            size_t index_j = nonIsoMuonsIndices[j];

                                            if ((*mu_charge)[index_i] != (*mu_charge)[index_j]) {
                                                
                                                 cutflow_10 = true;

                                                TLorentzVector muon1 = allMuons[index_i];
                                                TLorentzVector muon2 = allMuons[index_j];

                                                float relIso1 = calculate_relative_isolation((*mu_topoetcone20)[index_i], muon1.Pt());
                                                float relIso2 = calculate_relative_isolation((*mu_topoetcone20)[index_j], muon2.Pt());

                                                histograms[1][12]->Fill(relIso1, EventWT);
                                                histograms[1][12]->Fill(relIso2, EventWT); 
                                                float dR1_l = muon1.DeltaR(selectedJets[0]);
                                                float dR1_s = muon1.DeltaR(selectedJets[1]);
                                                float dR2_l = muon2.DeltaR(selectedJets[0]);
                                                float dR2_s = muon2.DeltaR(selectedJets[1]);

                                                if (relIso1 >= 0.1 && relIso2 >=0.1 &&  cutflow_10 && (dR1_l<1.0||dR1_s<1.0) && (dR2_l<1.0||dR2_s<1.0)) {
                                                     cutflow_11 = true;

                                                    float iso12_1 = -99.0;
                                                    float iso12_2 = -99.0;

                                                    //remove the inter-iso selection
                                                    muonsAfterChargeCutSet.insert(index_i);
                                                    muonsAfterChargeCutSet.insert(index_j);
                                                    nonIsolatedMuonPairs.push_back(make_pair(index_i, index_j));
                                                    
                                                    for (size_t k = 0; k < allMuons.size(); ++k) {
                                                        if (index_i == k) continue;
                                                        TLorentzVector lepk = allMuons[k];
                                                        float deltaR1 = muon1.DeltaR(lepk);
                                                        if (deltaR1 < 0.5) {
                                                            iso12_1 += lepk.Pt();
                                                        }
                                                    }
                                                    float interIso1 = iso12_1 / muon1.Pt();

                                                    for (size_t k = 0; k < allMuons.size(); ++k) {
                                                        if (index_j == k) continue;
                                                        TLorentzVector lepk = allMuons[k];
                                                        float deltaR2 = muon2.DeltaR(lepk);
                                                        if (deltaR2 < 0.5) {
                                                            iso12_2 += lepk.Pt();
                                                        }
                                                    }
                                                    float interIso2 = iso12_2 / muon2.Pt();

                                                    histograms[1][11]->Fill(interIso1, EventWT); 
                                                    histograms[1][11]->Fill(interIso2, EventWT); 
                                                   
                                                }
                                            }
                                        }
                                    }

                                    histograms[1][13]->Fill(muonsAfterChargeCutSet.size(), EventWT); 

                                    if ( cutflow_10) { CutFlow_sel->Fill(10, EventWT); } // CutFlow bin 10

                                   
                                    // Fill histograms related to opposite charge pairs
                                    
                                    for (const auto& pair : nonIsolatedMuonPairs) {
                                        size_t index_i = pair.first;
                                        size_t index_j = pair.second;
                                        
                                        TLorentzVector muon1 = allMuons[index_i];
                                        TLorentzVector muon2 = allMuons[index_j];

                                        

                                        TLorentzVector ll = muon1 + muon2;


                                        //Previous studies used the below dilepton mass criteria
                                        /*
                                        if (ll.M() > 19e3 || ll.M() < 12e3) {
                                            continue;
                                        }*/

                                        
                                        TLorentzVector jj;

                                        if (selectedJets.size() > 1) {
                                            jj = selectedJets[0] + selectedJets[1];
                                        }

                                        // Calculate transverse mass
                                        TVector3 metVector(0, 0, 0);
                                        metVector.SetPtEtaPhi(met_met, 0, met_phi);

                                        float M_ll = ll.M();
                                        float M_jj = jj.M();

                                        histograms[2][0]->Fill(scalarSumJetPt / 1e3, EventWT); // HT cutflow step 11
                                        histograms[2][1]->Fill(met_met / 1e3, EventWT); // MET cutflow step 11
                                        histograms[2][2]->Fill(met_met / MT2, EventWT); // RT cutflow step 11
                                        histograms[2][3]->Fill(MT2 / 1e3, EventWT); // MT (dijet + met) cutflow step 11
                                        histograms[2][4]->Fill(M_jj / 1e3, EventWT); // mjj cutflow step 11
                                        histograms[2][5]->Fill(M_ll / 1e3, EventWT); // mll cutflow step 11
                                        histograms[2][6]->Fill(selectedJets[0].Pt() / 1e3, EventWT); // Leading Jet pT cutflow step 11
                                        histograms[2][7]->Fill(selectedJets[1].Pt() / 1e3, EventWT); // Subleading Jet pT cutflow step 11
                                        histograms[2][8]->Fill(dphiMin, EventWT); // MINphi cutflow step 11
                                        histograms[2][9]->Fill(dphiMax, EventWT); // MAXphi cutflow step 11
                                        histograms[2][10]->Fill(fabs(dphiMin - dphiMax), EventWT); // MAX-MINphi cutflow step 11
                                    }
                                    if ( cutflow_11) {
                                        CutFlow_sel->Fill(11, EventWT); // CutFlow bin 11
                                        }
                                    
                                    
                                }
                            }
                        }
                    }
                }
            }
        }
    allMuons.clear();
    nonIsoMuonsIndices.clear();
    myjets.clear();
    selectedJets.clear();   
    } // end event loop
    
    // Finalize
    inputFile->Close();
    outputFile->cd();
    outputFile->Write();
    outputFile->Close();

    return 0;
}
