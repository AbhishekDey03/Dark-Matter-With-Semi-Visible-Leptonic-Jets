/*
 * A plotting macro to create ATLAS style plots. The plots are all log-y scale. The macro takes in a list of legend labels and root files.
 * The macro is dependent on all the files having the same structure, therefore you should ensure that you have hadded all the background DSIDs into one .root file
 * 
 * Author: Abhishek Ranjan Dey
 * Created: 14/08/2023
 * 
 * Dependencies:
 * - ROOT
 * 
 * 
 * Example Compile Line:
 * root -l -b -q plotting_macro.C
 * 
 * The data files are not provided and therefore the analysis code requires you change the input and output file paths
 */


#include "TFile.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TPaveText.h"
#include "TText.h"
#include "TPad.h"
#include "TLatex.h"
#include <map>
#include <string>
#include <vector>
#include <iostream>

// Function to add the ATLAS label to the plot
void add_atlas_label(TPad* pad, const char* text="") {
    pad->cd();
    TLatex t;
    t.SetNDC(); // Set coordinates relative to the pad (0 to 1)
    t.SetTextFont(42);
    t.SetTextSize(0.05);

    // Add standard ATLAS label and additional text if provided
    t.DrawLatex(0.15, 0.83, "#font[72]{#scale[1]{ATLAS Work in Progress}}");
    t.DrawLatex(0.14, 0.76, "#scale[1]{#scale[0.7]{#sqrt{s} = 13 TeV}}");
    t.DrawLatex(0.15, 0.70, "#scale[0.7]{jet R=1.0}");

    if (text && text[0] != '\0') {
        t.DrawLatex(0.15, 0.64, text);
    }
}

// Function to plot and compare histograms from multiple files
void macroplot(const std::vector<std::string>& file_paths, const std::vector<std::string>& labels) {
    if (file_paths.size() != labels.size()) {
        std::cerr << "Error: Number of file paths and labels must be the same!" << std::endl;
        return;
    }

    // Set ATLAS style preferences
    gStyle->SetOptStat(0);
    gStyle->SetTitleFont(42, "XYZ");
    gStyle->SetLabelFont(42, "XYZ");
    gStyle->SetTitleSize(0.05, "XYZ");
    gStyle->SetLabelSize(0.04, "XYZ"); 
    gStyle->SetHistLineWidth(2); 

    // Open the ROOT files
    std::vector<TFile*> inputFiles;
    for (const auto& file_path : file_paths) {
        TFile* inputFile = TFile::Open(file_path.c_str(), "READ");
        if (!inputFile) {
            std::cerr << "Error opening file: " << file_path << std::endl;
            return;
        }
        inputFiles.push_back(inputFile);
    }

    // Define the names of all histograms to be plotted
    const std::vector<std::string> hist_names = {
        "HT_cutflow_0", "MET_cutflow_0", "RT_cutflow_0", "MT_cutflow_0", "mjj_cutflow_0", "mll_cutflow_0", "leading_jet_pt_cutflow_0",
        "subleading_jet_pt_cutflow_0", "MINphi_cutflow_0", "MAXphi_cutflow_0", "MAX_MINphi_cutflow_0", "inter_iso_cutflow_0",
        "rel_iso_cutflow_0", "n_muons_non_iso_cutflow_0", "HT_cutflow_1", "MET_cutflow_1", "RT_cutflow_1", "MT_cutflow_1",
        "mjj_cutflow_1", "mll_cutflow_1", "leading_jet_pt_cutflow_1", "subleading_jet_pt_cutflow_1", "MINphi_cutflow_1",
        "MAXphi_cutflow_1", "MAX_MINphi_cutflow_1", "inter_iso_cutflow_1", "rel_iso_cutflow_1", "n_muons_non_iso_cutflow_1",
        "HT_cutflow_2", "MET_cutflow_2", "RT_cutflow_2", "MT_cutflow_2", "mjj_cutflow_2", "mll_cutflow_2", "leading_jet_pt_cutflow_2",
        "subleading_jet_pt_cutflow_2", "MINphi_cutflow_2", "MAXphi_cutflow_2", "MAX_MINphi_cutflow_2", "inter_iso_cutflow_2",
        "rel_iso_cutflow_2", "n_muons_non_iso_cutflow_2", "CutFlow_sel"
    };

    // Define the corresponding labels for the histograms
    const std::vector<std::string> hist_labels = {
        "HT (CutFlow step 3)", "MET (CutFlow step 3)", "RT (CutFlow step 3)", "MT (dijet + met) (CutFlow step 3)", "mjj (CutFlow step 3)",
        "mll (CutFlow step 3)", "Leading Jet pT (CutFlow step 3)", "Subleading Jet pT (CutFlow step 3)", "MINphi (closest jet, MET) (CutFlow step 3)",
        "MAXphi (farthest jet, MET) (CutFlow step 3)", "MAX-MINphi (CutFlow step 3)", "Inter-iso (CutFlow step 3)", "Rel-iso (CutFlow step 3)",
        "Number of non-iso muons (CutFlow step 3)", "HT (CutFlow step 7)", "MET (CutFlow step 7)", "RT (CutFlow step 7)", "MT (dijet + met) (CutFlow step 7)",
        "mjj (CutFlow step 7)", "mll (CutFlow step 7)", "Leading Jet pT (CutFlow step 7)", "Subleading Jet pT (CutFlow step 7)", "MINphi (closest jet, MET) (CutFlow step 7)",
        "MAXphi (farthest jet, MET) (CutFlow step 7)", "MAX-MINphi (CutFlow step 7)", "Inter-iso (CutFlow step 7)", "Rel-iso (CutFlow step 7)",
        "Number of non-iso muons (CutFlow step 7)", "HT (CutFlow step 11)", "MET (CutFlow step 11)", "RT (CutFlow step 11)", "MT (dijet + met) (CutFlow step 11)",
        "mjj (CutFlow step 11)", "mll (CutFlow step 11)", "Leading Jet pT (CutFlow step 11)", "Subleading Jet pT (CutFlow step 11)",
        "MINphi (closest jet, MET) (CutFlow step 11)", "MAXphi (farthest jet, MET) (CutFlow step 11)", "MAX-MINphi (CutFlow step 11)",
        "Inter-iso (CutFlow step 11)", "Rel-iso (CutFlow step 11)", "Number of non-iso muons (CutFlow step 11)", "CutFlow"
    };

    // Create a canvas for plotting and open a PDF file to save the plots
    TCanvas *canvas = new TCanvas("canvas", "canvas", 1600, 1200); // Increase canvas size
    canvas->Print("/eos/home-a/abdey/ntupleanalyser_output/Plots/background/signal_vs_background_MET_Trigger_no_mll.pdf[");

    // Set margins for the canvas
    canvas->SetLeftMargin(0.12);   // Increase left margin
    canvas->SetRightMargin(0.12);  // Increase right margin
    canvas->SetTopMargin(0.12);    // Increase top margin
    canvas->SetBottomMargin(0.12); // Increase bottom margin

    // Define a list of colors for the histograms
    std::vector<int> colors = {kRed, kBlue, kGreen, kMagenta, kCyan}; 

    // Loop through each histogram, retrieve it from each file, and plot
    for (size_t idx = 0; idx < hist_names.size(); ++idx) {
        bool histograms_exist = true;
        std::vector<TH1*> hists;

        // Retrieve histograms from each input file
        for (size_t i = 0; i < inputFiles.size(); ++i) {
            TH1* hist = dynamic_cast<TH1*>(inputFiles[i]->Get(hist_names[idx].c_str()));
            if (!hist || hist->GetEntries() == 0) { // Skip histograms with 0 entries
                std::cerr << "Histogram " << hist_names[idx] << " not found or has 0 entries in file " << file_paths[i] << "!" << std::endl;
                histograms_exist = false;
                break;
            }
            hists.push_back(hist);
        }

        if (!histograms_exist) continue;

        // Adjust axis title offsets and alignment for better visibility
        for (auto hist : hists) {
            hist->GetXaxis()->SetTitleOffset(1.15); // Move x-axis title further away
            hist->GetYaxis()->SetTitleOffset(1.15); // Move y-axis title further away
            hist->GetXaxis()->CenterTitle(true);
            hist->GetYaxis()->CenterTitle(true);

            // Remove plot titles and set axis titles using the histogram labels
            hist->SetTitle("");
            hist->GetXaxis()->SetTitle(hist_labels[idx].c_str());
            hist->GetYaxis()->SetTitle("Events");

            double maxY = 0;
            double minY = 1e4;
            for (auto hist : hists) {
                double histMax = hist->GetMaximum();
                double histMin = hist->GetMinimum();
                if (histMax > maxY) maxY = histMax;
                if (histMin < minY && histMin > 0) minY = histMin; 
            }

            for (auto hist : hists) hist->SetMaximum(maxY * 1e2); 
            for (auto hist : hists) hist->SetMinimum(minY / 1e2); 
        }

        // Determine maximum and minimum Y values for setting Y-axis limits
        double maxY = 0;
        double minY = 1e4;
        for (auto hist : hists) {
            double histMax = hist->GetMaximum();
            double histMin = hist->GetMinimum();
            if (histMax > maxY) maxY = histMax;
            if (histMin < minY && histMin > 0) minY = histMin; 
        }

        for (auto hist : hists) hist->SetMaximum(maxY * 1.2); 
        for (auto hist : hists) hist->SetMinimum(minY / 1.2); 

        // Draw histograms on the canvas with log scale on the y-axis
        canvas->SetLogy(true);
        TLegend *legend = new TLegend(0.58, 0.75, 0.87, 0.85);
        legend->SetTextSize(0.03);
        legend->SetBorderSize(0);

        for (size_t i = 0; i < hists.size(); ++i) {
            hists[i]->SetLineColor(colors[i % colors.size()]);
            hists[i]->SetFillColor(0); // Remove fill color for clarity
            hists[i]->SetLineStyle(2);
            if (i == 0) {
                hists[i]->Draw("hist");
            } else {
                hists[i]->Draw("hist same");
            }
            legend->AddEntry(hists[i], labels[i].c_str(), "l");
        }

        legend->Draw();
        add_atlas_label(canvas, " ");
        canvas->Print("/eos/home-a/abdey/ntupleanalyser_output/Plots/background/signal_vs_background_MET_Trigger_no_mll.pdf");

        canvas->Clear(); // Clear the canvas for the next plot
    }

    // Close the PDF file and clean up
    canvas->Print("/eos/home-a/abdey/ntupleanalyser_output/Plots/background/signal_vs_background_MET_Trigger_no_mll.pdf]");
    for (auto file : inputFiles) file->Close(); // Close all input files
    delete canvas;
}

// Main function that sets up the files and labels and calls the macroplot function
void plotting_macro() {
    // Define the paths to the ROOT files to be analyzed
    std::vector<std::string> file_paths = {
        "/eos/user/a/abdey/ntupleanalyser_output/SVJl_output_NoOLR_ROOT_23_07.root",
        "/eos/user/a/abdey/ntupleanalyser_output/background/Zmumu/trigger_selects_1/all_Z_mumu_25_07.root",
        "/eos/user/a/abdey/ntupleanalyser_output/background/ttbar/trigger_selects_1/all_ttbar.root"
    };

    // Define the labels corresponding to the ROOT files
    std::vector<std::string> labels = {
        "signal",
        "Z mu mu background",
        "ttbar background"
    };

    // Call the macroplot function to create and save the plots
    macroplot(file_paths, labels);
}
