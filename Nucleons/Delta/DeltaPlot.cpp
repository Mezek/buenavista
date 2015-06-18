/**
 * $Date: $
 * $Revision: $
 * $Author: $
 * $HeadURL: $
 * $Id: $
 *
 * @file
 * @brief	Preview of functions and data
 */

#include "DeltaPlot.h"

DeltaPlot::DeltaPlot ( std::size_t p ): k(0), c(p) {

	// To set white color for graph
    TStyle *plain = new TStyle("Plain","Plain Style(no colors/fill areas)");
	plain->SetCanvasBorderMode(0);
    plain->SetPadBorderMode(0);
    plain->SetPadColor(0);
    plain->SetCanvasColor(0);
    plain->SetTitleColor(1);
    plain->SetStatColor(0);
	plain->SetTitleBorderSize(0);
	//gROOT->SetStyle(plain);
	plain->cd();
	title = "How does function look like?";
	x0 = 10;
	y0 = 10;
	s = 25;
	w = 900;
	h = 700;
}

/// One graph with optional title

void DeltaPlot::view (Int_t num, Double_t axisX[], Double_t axisY[], Char_t const* title) {

	c[k] = new TCanvas (uName("c",k), uName("Graph_",k), x0+k*s, y0+k*s, w, h);
	c[k]->SetLogy();
	
	TGraph *gr1 = new TGraph (num, axisX, axisY);
	
    gr1->Draw("AL");
    gr1->GetXaxis()->CenterTitle();
    gr1->GetYaxis()->CenterTitle();

	gr1->GetXaxis()->SetTitle("t [GeV^{2}]");
	gr1->GetYaxis()->SetTitle(title);
	//gr1->SetTitle(title);
	gr1->SetTitle();
	++k;
}

/// View data

void DeltaPlot::viewData (Int_t num, Double_t axisX[], Double_t axisY[], Char_t const* title) {

	c[k] = new TCanvas (uName("c",k), uName("Graph_",k), x0+k*s, y0+k*s, w, h);
    //c[k]->SetLogx(); // logarithmic scale
	
	TGraph *gr1 = new TGraph (num, axisX, axisY);
	gr1->Draw("AP");

	gr1->SetMarkerStyle(22);
	gr1->SetMarkerSize(1.1);
	gr1->SetMarkerColor(9);

	gr1->GetXaxis()->CenterTitle();
    gr1->GetYaxis()->CenterTitle();

	gr1->GetXaxis()->SetTitle("Q^2 [GeV^{2}]");
	gr1->GetYaxis()->SetTitle(title);
	gr1->SetTitle(title);
	//gr1->SetTitle();
	++k;
}

/// One graph with data points and title

void DeltaPlot::viewPlusData (Int_t num, Double_t axisX[], Double_t axisY[], Int_t numD, Double_t axisXD[], Double_t axisYD[], Char_t const* title) {

	c[k] = new TCanvas (uName("c",k), uName("Graph_",k), x0+k*s, y0+k*s, w, h);
	c[k]->SetLogy(); // logarithmic scale

	TGraph *gr1 = new TGraph (num, axisX, axisY);
	gr1->Draw("AL");
    //gr1->SetTitle("Function with experimental points");
    gr1->SetTitle(title);
    gr1->GetXaxis()->CenterTitle();
    gr1->GetYaxis()->CenterTitle();
	//gr1->GetXaxis()->SetLimits(0.5,0.);
	//tl = new TLatex();
	//tl->DrawLatex(0.22,0.15,"#sqrt{s} = 500(GeV)");

    TGraph *gr2 = new TGraph (numD, axisXD, axisYD);
	gr2->Draw("P");
	gr2->SetMarkerStyle(21);
	gr2->SetMarkerSize(.5);
	++k;	
}

/*
 * 	Double_t X1[] = {2.8, 4.0};
	Double_t Y1[] = {-2.0, -3.1};
	Double_t U1[] = {1.3, 1.3};
	Double_t D1[] = {1.3, 1.3};
	g[1] = new TGraph (2, X1, Y1);
	//gr1->SetTitle(title);
	g[1]->SetMarkerColor(kBlue);
	g[1]->SetMarkerStyle(21);
	mg->Add(g[1]);

	Double_t X2[] = {2.8, 4.0};
	Double_t Y2[] = {-2.0, -3.1};
	Double_t U2[] = {1.3, 1.3};
	Double_t D2[] = {1.3, 1.3};
	*/
