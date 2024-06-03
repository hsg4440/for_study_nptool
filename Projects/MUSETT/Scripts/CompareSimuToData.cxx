#include <TFile.h>
#include <TH2.h>
#include <fstream>
#include <iostream>
#include <sstream>
#include <TH2.h>
#include <TMath.h>
#include <TCanvas.h>
#include <TH2.h>
#include <TH1.h>
#include <TStyle.h>
#include <TString.h>
#include <TText.h>
#include <TLatex.h>

double ComputeChi2_TH2(TH2* histData, TH2* histSimu) {
    // Ensure the histograms have the same binning
    if (histData->GetNbinsX() != histSimu->GetNbinsX() || histData->GetNbinsY() != histSimu->GetNbinsY()) {
        std::cerr << "Histograms have different binning!" << std::endl;
        return -1;
    }

    // Normalize the simulation histogram to the data histogram
    double dataEntries = histData->GetEntries();
    double simuEntries = histSimu->GetEntries();
    if (simuEntries == 0) {
        std::cerr << "Simulation histogram has no entries!" << std::endl;
        return -1;
    }
    histSimu->Scale(dataEntries / simuEntries);

    // Compute the Chi2
    double chi2 = 0;
    int nBinsX = histData->GetNbinsX();
    int nBinsY = histData->GetNbinsY();
    int nBins = 0;

    for (int i = 1; i <= nBinsX; ++i) {
        for (int j = 1; j <= nBinsY; ++j) {
            double observed = histData->GetBinContent(i, j);
            double expected = histSimu->GetBinContent(i, j);
            double error_data = histData->GetBinError(i, j); // Assuming Poisson statistics
            double error_simu = histSimu->GetBinError(i, j); // Assuming Poisson statistics
            double error_tot = TMath::Sqrt(error_data*error_data + error_simu*error_simu);
            //std::cout << "observed = " << observed << std::endl;
            //std::cout << "expected = " << expected << std::endl;
            //std::cout << "error = " << error_tot << std::endl;
            //std::cout << "sqrt(observed) = " << (double) TMath::Sqrt(observed) << std::endl;
            if (error_tot == 0) error_tot = 1; // Avoid division by zero, assume minimum error
            chi2 += TMath::Power((observed - expected) / error_tot, 2);
            nBins++;
        }
    }

    return chi2 / nBins; // Normalized Chi2
}

double ComputeChi2_TH1(TH1* histData, TH1* histSimu) {
    // Ensure the histograms have the same binning
    if (histData->GetNbinsX() != histSimu->GetNbinsX()) {
        std::cerr << "Histograms have different binning!" << std::endl;
        return -1;
    }

    // Normalize the simulation histogram to the data histogram
    double dataEntries = histData->GetEntries();
    double simuEntries = histSimu->GetEntries();

    // Compute the Chi2
    double chi2 = 0;
    int nBins = histData->GetNbinsX();

    for (int i = 1; i <= nBins; ++i) {
        double observed = histData->GetBinContent(i);
        double expected = histSimu->GetBinContent(i);
        double error_data = histData->GetBinError(i); // Assuming Poisson statistics
        double error_simu = histSimu->GetBinError(i); // Assuming Poisson statistics
        double error_tot = TMath::Sqrt(error_data * error_data + error_simu * error_simu);
        if (error_tot == 0) error_tot = 1; // Avoid division by zero, assume minimum error
        chi2 += TMath::Power((observed - expected) / error_tot, 2);
    }

    return chi2 / nBins; // Normalized Chi2
}



void PlotComparison(TH2* histData, TH2* histSimu, const std::string& outputFileName, double chi2, int detectorIndex) {
    // Create a canvas for the full hit pattern comparison
    TCanvas* canvasFull = new TCanvas(Form("canvasFull_%d", detectorIndex), "Comparison Canvas", 1400, 700);
    canvasFull->Divide(2, 1); // Divide the canvas into two pads

    // Plot the simulation histogram
    canvasFull->cd(1);
    histSimu->SetTitle(Form("DET%d, Simulation", detectorIndex));
    histSimu->Draw();
    histSimu->Rebin2D(2,2);
    gPad->SetTopMargin(0.1);
    gPad->SetRightMargin(0.05);
    gPad->SetLeftMargin(0.12);
    gPad->Modified();
    histSimu->SetStats(0); // Enable statistics box
    gStyle->SetStatX(0.4); // Set statistics box position
    gStyle->SetStatY(0.9);

    // Plot the data histogram
    canvasFull->cd(2);
    histData->SetTitle(Form("DET%d, Data", detectorIndex));
    histData->Draw();
    gPad->SetTopMargin(0.1);
    gPad->SetRightMargin(0.05);
    gPad->SetLeftMargin(0.12);
    gPad->Modified();

    histData->SetStats(0); // Enable statistics box
    histData->Rebin2D(2,2);
    gStyle->SetStatX(0.4); // Set statistics box position
    gStyle->SetStatY(0.9);

    canvasFull->cd();
    TLatex latex;
    latex.SetTextSize(0.03);
    latex.SetNDC();
    //latex.DrawLatex(0.5, 0.96, Form("#chi^{2} / N_{pixels}: %.4f", chi2));

    // Save the full hit pattern comparison canvas as a PNG
    std::string fullOutputFileName = outputFileName + "_full.png";
    canvasFull->SaveAs(fullOutputFileName.c_str());
    delete canvasFull;

    // Create a canvas for the TH1 projections
    TCanvas* canvasProjX = new TCanvas(Form("canvasProjX_%d", detectorIndex), "Projections Canvas X", 2400, 1200);
    canvasProjX->Divide(4, 2); // Divide the canvas into 8 pads

    int nBinsX = histData->GetNbinsX();

    // Plot X-axis projections
    for (int i = 0; i < 8; ++i) {
        canvasProjX->cd(i + 1);
        int binLow = i * (nBinsX / 8) + 1;
        int binHigh = (i + 1) * (nBinsX / 8);

        TH1D* projXData = histData->ProjectionX(Form("projXData_%d_%d", detectorIndex, i), binLow, binHigh);
        TH1D* projXSimu = histSimu->ProjectionX(Form("projXSimu_%d_%d", detectorIndex, i), binLow, binHigh);

        projXData->SetStats(0);
        projXSimu->SetStats(0);


        projXData->SetLineColor(kRed);
        projXSimu->SetLineColor(kBlue);

        projXData->SetFillColorAlpha(kRed, 0.5);
        projXSimu->SetFillColorAlpha(kBlue, 0.5);

        projXData->SetFillStyle(3007);
        projXSimu->SetFillStyle(3006);

        projXSimu->SetLineWidth(3);
        projXData->SetLineWidth(3);

        projXData->Rebin(4);
        projXSimu->Rebin(4);

        projXData->SetTitle(Form("X Proj Bins %d-%d", binLow, binHigh));


        double chi2ProjX = ComputeChi2_TH1(projXData, projXSimu);




        projXData->Draw("E HIST");
        projXSimu->Draw("E HIST SAME");


        // Add Chi2 value to the subplot
        TPaveText* paveTextX =  new TPaveText(0.4, 0.2, 0.6, 0.4, "NDC");
        paveTextX->AddText(Form("#chi^{2}: %.2f", chi2ProjX));
        paveTextX->SetFillColor(0);
        paveTextX->SetTextColor(kBlack);
        paveTextX->Draw("same");
        gStyle->SetLineScalePS(0.5);
        gPad->Modified();
    }

    std::string projXOutputFileName = outputFileName + "_projX.png";
    canvasProjX->SaveAs(projXOutputFileName.c_str());
    delete canvasProjX;

    TCanvas* canvasProjY = new TCanvas(Form("canvasProjY_%d", detectorIndex), "Projections Canvas Y", 2400, 1200);
    canvasProjY->Divide(4, 2); // Divide the canvas into 8 pads

    int nBinsY = histData->GetNbinsY();

    // Plot Y-axis projections
    for (int i = 0; i < 8; ++i) {
        canvasProjY->cd(i + 1);
        int binLow = i * (nBinsY / 8) + 1;
        int binHigh = (i + 1) * (nBinsY / 8);

        TH1D* projYData = histData->ProjectionY(Form("projYData_%d_%d", detectorIndex, i), binLow, binHigh);
        TH1D* projYSimu = histSimu->ProjectionY(Form("projYSimu_%d_%d", detectorIndex, i), binLow, binHigh);


        projYData->SetStats(0);
        projYSimu->SetStats(0);


        projYData->SetLineColor(kRed);
        projYSimu->SetLineColor(kBlue);
        projYData->SetFillColorAlpha(kRed, 0.9);
        projYSimu->SetFillColorAlpha(kBlue, 0.5);
        projYData->SetFillStyle(3003);
        projYSimu->SetFillStyle(3006);
        projYSimu->SetLineWidth(3);
        projYData->SetLineWidth(3);

        projYData->Rebin(4);
        projYSimu->Rebin(4);


        projYData->SetTitle(Form("Y Proj Bins %d-%d", binLow, binHigh));





        projYData->Draw("E HIST");
        projYSimu->Draw("E HIST SAME");

        // Compute Chi2 for the projection
        double chi2ProjY = ComputeChi2_TH1(projYData, projYSimu);
        // Add Chi2 value to the subplot
        TPaveText* paveTextY = new TPaveText(0.4, 0.2, 0.6, 0.4, "NDC");
        paveTextY->AddText(Form("#chi^{2} : %.2f", chi2ProjY));
        paveTextY->SetFillColor(0);
        paveTextY->SetTextColor(kBlack);
        paveTextY->Draw("same");
        gStyle->SetLineScalePS(0.5);
        gPad->Modified();
    }

    std::string projYOutputFileName = outputFileName + "_projY.png";
    canvasProjY->SaveAs(projYOutputFileName.c_str());
    delete canvasProjY;
}


void PlotEnergyCorrelation(TH2* histData, TH2* histSimu)
{
    // Create a canvas
    TCanvas *c1 = new TCanvas("c1", "Energy Correlation", 1200, 600);

    // Create two pads, one for each plot
    c1->Divide(2, 1);

    // On the left plot the TH2 of simulation
    c1->cd(1);
    gPad->SetLogz(); // Set log scale for the z-axis
    histSimu->GetXaxis()->SetRangeUser(0, 8); // Set X-axis range
    histSimu->GetYaxis()->SetRangeUser(0, 8); // Set Y-axis range
    gPad->SetRightMargin(0.15); // To avoid cutting off the color scale
    histSimu->SetTitle("SIMULATION"); // Set the title
    histSimu->GetXaxis()->SetTitle("E_{1} [MeV]"); // Set X-axis label
    histSimu->GetYaxis()->SetTitle("E_{2} [MeV]"); // Set Y-axis label
    histSimu->Draw("COLZ");

    // Print the number of entries for debugging
    std::cout << "Nentries (Simulation) = " << histSimu->GetEntries() << std::endl;
    std::cout << "Nentries (Data) = " << histData->GetEntries() << std::endl;

    // Rescale number of entries of Data to match those of Simulation
    double scaleFactor = histSimu->GetEntries() / histData->GetEntries();
    //histData->Scale(5*scaleFactor);

    // Plot data on the right with the same color scale
    c1->cd(2);
    gPad->SetRightMargin(0.15); // To avoid cutting off the color scale
    gPad->SetLogz(); // Set log scale for the z-axis
    histData->GetXaxis()->SetRangeUser(0, 8); // Set X-axis range
    histData->GetYaxis()->SetRangeUser(0, 8); // Set Y-axis range
    histData->SetStats(0);
    histData->SetTitle("DATA"); // Set the title
    histData->GetXaxis()->SetTitle("E_{1} [MeV]"); // Set X-axis label
    histData->GetYaxis()->SetTitle("E_{2} [MeV]"); // Set Y-axis label
    histData->Draw("COLZ");

    // Set the color scale the same for both histograms
    double min = std::min(histSimu->GetMinimum(), histData->GetMinimum());
    double max = std::max(histSimu->GetMaximum(), histData->GetMaximum());
    histSimu->SetStats(0);
    //histSimu->GetZaxis()->SetRangeUser(min, max);
    //histData->GetZaxis()->SetRangeUser(min, max);

    // Update the canvas to draw the histograms
    c1->Update();
    c1->Modified(); // Ensure the canvas is updated and modified
}


void CompareSimuToData() {
    // Open the simulation file
    std::string dataFileName = "/Users/lh270370/Software/nptool/Outputs/Analysis/Data/222Ra/histoSanity_R275.root";
    std::string simuFileName = "/Users/lh270370/Software/nptool/Outputs/Analysis/Simu/Isolde/full_test_beam/histoSanity_222Ra_beam_test_AlGrid__5e6offsetTargetXm3__offsetTargetY6_0__offsetBeamXm6__offsetBeamYm150__sigmaBeamX0_2__sigmaBeamY0_2__holderShift0_5_physics.root";
    TFile* simuFile = TFile::Open(simuFileName.c_str());
    if (!simuFile) {
        std::cerr << "Failed to open simulation file: " << simuFileName << std::endl;
        return;
    }

    TDirectory *simuGeometryDir = (TDirectory*)simuFile->Get("Geometry/HitPattern");


    // Open the data file
    TFile* dataFile = TFile::Open(dataFileName.c_str());
    if (!dataFile) {
        std::cerr << "Failed to open data file: " << dataFileName << std::endl;
        simuFile->Close();
        return;
    }

    TDirectory *dataGeometryDir = (TDirectory*)dataFile->Get("Geometry/HitPattern");

    TDirectory *dataEnergyCorrelationDir = (TDirectory*)dataFile->Get("Energy correlation");
    TDirectory *simuEnergyCorrelationDir = (TDirectory*)simuFile->Get("Energy correlation");


    TH2F *histSimuEnergy = (TH2F*)simuEnergyCorrelationDir->Get("E1E2_all_diffDet");
    TH2F *histDataEnergy = (TH2F*)dataEnergyCorrelationDir->Get("E1E2_all_diffDet");


    PlotEnergyCorrelation(histDataEnergy,histSimuEnergy);
    std::cout << "After plot " << std::endl;
    //std::ofstream outFile("Chi2_results_2.txt", std::ios_base::app);


    /*
    // Loop over each detector and compare histograms
    for (int i = 0; i < 4; ++i) {
        // Construct histogram names
        std::string histName = "HitPattern_Det" + std::to_string(i);
        std::string dataHistName = histName;

        // Retrieve histograms
        TH2 *histSimu = (TH2*)simuGeometryDir->Get(histName.c_str());
        TH2 *histData = (TH2*)dataGeometryDir->Get(dataHistName.c_str());
        if (!histSimu || !histData) {
            std::cerr << "Histograms not found for: " << histName << std::endl;
            continue;
        }

        // Calculate Chi2
        double chi2 = ComputeChi2_TH2(histData, histSimu);
        std::cout << "Chi2 for " << histName << ": " << chi2 << std::endl;

        // Log the Chi2 value
        //outFile << "Detector " << i << " Chi2: " << chi2 << std::endl;

        // Plot comparison
        std::string outputFileName = "Beam_simu_data_DET" + std::to_string(i) + "_comparison";
        PlotComparison(histData, histSimu, outputFileName, chi2, i);
    }
    */
    //outFile.close();
    //simuFile->Close();
    //dataFile->Close();
}
