//This code draws the histogram of pulseH for MinBias events(when triggerBits==4) and does the Gaussian Fit and extract the values of mean and sigma. Also prints the values of mean+3sigma.
// Which is later used for the data cut. We will use data cut conditions PMT Hit > 3 if the pulseH >mu+3sigma during data selection for Michel electron calibration analysis.
#include <iostream>
#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TF1.h>
#include <TCanvas.h>
#include <vector>
#include <TStyle.h>
#include <TLatex.h>
#include <TPaveStats.h> // Include TPaveStats header

using namespace std;

// Mapping of PMT channels
const int pmtChannelMap[12] = {0, 10, 7, 2, 6, 3, 8, 9, 11, 4, 5, 1};

// Function to define baseline and noise distribution
void DefineBaselineAndNoise(const char *fileName) {
    TFile *file = TFile::Open(fileName);
    if (!file || file->IsZombie()) {
        cerr << "Error opening file: " << fileName << endl;
        return;
    }

    TTree *tree = (TTree*)file->Get("tree");
    if (!tree) {
        cerr << "Error accessing TTree!" << endl;
        file->Close();
        return;
    }

    // Variables to read from the tree
    Double_t pulseH[23];
    Int_t triggerBits;

    tree->SetBranchAddress("pulseH", pulseH);
    tree->SetBranchAddress("triggerBits", &triggerBits);

    // Histograms for pulse heights of each PMT
    TH1F *histPulseH[12];
    for (int i = 0; i < 12; i++) {
        histPulseH[i] = new TH1F(Form("PMT%d", i + 1),
                                Form(";ADC Counts;Events"), // Empty title
                                40, -10, 30); // 40 bins, range from -10 to 30
    }

    // Loop over all events to fill histograms for min bias events
    Long64_t nEntries = tree->GetEntries();
    for (Long64_t entry = 0; entry < nEntries; entry++) {
        tree->GetEntry(entry);

        // Select min bias events (triggerBits == 4)
        if (triggerBits != 4) continue;

        // Fill pulse height histograms for each PMT
        for (int pmt = 0; pmt < 12; pmt++) {
            histPulseH[pmt]->Fill(pulseH[pmtChannelMap[pmt]]);
        }
    }

    // Create a master canvas for the combined plot
    TCanvas *masterCanvas = new TCanvas("MasterCanvas", "Combined PMT Energy Distributions", 4000, 3900);
    masterCanvas->Divide(3, 4, 0.01, 0.01); // 3 columns, 4 rows, with small spacing between subplots

    // Add a title to the combined canvas
    masterCanvas->cd(); // Move to the main canvas area
    TLatex title;
    title.SetTextSize(0.02);
    title.DrawLatexNDC(0.4, 0.98, "Min Bias PulseH");

    // Define the layout of PMT channels on the canvas
    int layout[4][3] = {
        {9, 3, 7},  // Row 1: PMT 10, PMT 4, PMT 8
        {5, 4, 8},  // Row 2: PMT 6, PMT 5, PMT 9
        {0, 6, 1},  // Row 3: PMT 1, PMT 7, PMT 2
        {10, 11, 2} // Row 4: PMT 11, PMT 12, PMT 3
    };

    // Fit histograms with Gaussian and extract mean (μ) and standard deviation (σ)
    Double_t baselineMean[12] = {0};
    Double_t baselineSigma[12] = {0};

    // Increase the size of the statistical box
    gStyle->SetOptStat(1111); // Display all statistics (entries, mean, RMS, etc.)
    gStyle->SetStatFontSize(0.12); // Increase font size of the statistical box
    gStyle->SetStatW(0.25); // Increase width of the statistical box
    gStyle->SetStatH(0.10); // Increase height of the statistical box
    gStyle->SetStatX(0.95); // Position the statistical box at 90% of the pad width
    gStyle->SetStatY(0.90); // Position the statistical box at 90% of the pad height

    for (int row = 0; row < 4; row++) {
        for (int col = 0; col < 3; col++) {
            int pmtIndex = layout[row][col]; // Get PMT index from layout
            int padIndex = row * 3 + col + 1; // Calculate pad index (1-12)

            masterCanvas->cd(padIndex); // Select the appropriate pad

            if (histPulseH[pmtIndex]->GetEntries() == 0) {
                cerr << "Empty histogram for PMT " << pmtIndex + 1 << endl;
                continue;
            }

            // Dynamically determine the fit range
            Double_t xmin = histPulseH[pmtIndex]->GetXaxis()->GetXmin(); // Minimum x-value with data
            Double_t xmax = histPulseH[pmtIndex]->GetXaxis()->GetXmax(); // Maximum x-value with data

            // Fit with a Gaussian function
            TF1 *gaussFit = new TF1("gaussFit", "gaus", xmin, xmax); // Define the Gaussian function
            histPulseH[pmtIndex]->Fit(gaussFit, "RQ"); // Fit the histogram

            // Extract mean (μ) and standard deviation (σ)
            baselineMean[pmtIndex] = gaussFit->GetParameter(1); // Mean (μ)
            baselineSigma[pmtIndex] = gaussFit->GetParameter(2); // Standard deviation (σ)

            // Calculate μ + 3σ
            Double_t muPlus3Sigma = baselineMean[pmtIndex] + 3 * baselineSigma[pmtIndex];

            // Print only PMT number and μ + 3σ
            cout << "PMT " << pmtIndex + 1 << ": μ + 3σ = " << muPlus3Sigma << endl;

            // Adjust histogram text size
            histPulseH[pmtIndex]->GetXaxis()->SetTitleSize(0.07); // Increase x-axis title size
            histPulseH[pmtIndex]->GetYaxis()->SetTitleSize(0.06); // Increase y-axis title size
            histPulseH[pmtIndex]->GetXaxis()->SetLabelSize(0.04); // Increase x-axis label size
            histPulseH[pmtIndex]->GetYaxis()->SetLabelSize(0.04); // Increase y-axis label size

            // Ensure Y-axis label is visible
            histPulseH[pmtIndex]->GetYaxis()->SetTitle("Events per 1 ADC");
            histPulseH[pmtIndex]->GetYaxis()->SetTitleOffset(1.2); // Move Y-axis title closer to the axis line

            // Ensure X-axis label is visible
            histPulseH[pmtIndex]->GetXaxis()->SetTitle("ADC Counts");

            // Adjust margins for each subplot
            gPad->SetLeftMargin(0.15);   // Increase left margin
            gPad->SetRightMargin(0.10);  // Adjust right margin
            gPad->SetBottomMargin(0.15); // Increase bottom margin
            gPad->SetTopMargin(0.15);    // Adjust top margin

            // Draw the histogram and fit
            histPulseH[pmtIndex]->Draw();
            gaussFit->Draw("same");

            // Access the statistical box and increase its font size
            TPaveStats *stats = (TPaveStats*)histPulseH[pmtIndex]->FindObject("stats");
            if (stats) {
                stats->SetTextSize(0.06); // Increase font size of the statistical box
                stats->SetTextFont(42);   // Set font type (42 is the default font)
                stats->Draw(); // Redraw the statistical box
            }

            // Force the pad to update
            gPad->Modified();
            gPad->Update();

            delete gaussFit;
        }
    }

    // Save the master canvas
    masterCanvas->SaveAs("CombinedPMTEnergyDistributions.png");

    // Cleanup
    delete masterCanvas;
    for (int i = 0; i < 12; i++) delete histPulseH[i];
    file->Close();
}

// Main function
int main(int argc, char* argv[]) {
    if (argc != 2) {
        cerr << "Usage: " << argv[0] << " <input_file.root>" << endl;
        return 1;
    }

    // Call the function to define baseline and noise
    DefineBaselineAndNoise(argv[1]);

    return 0;
}
