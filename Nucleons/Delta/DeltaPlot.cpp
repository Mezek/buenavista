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

#include "../modules/PlotGraph.h"

PlotGraph::PlotGraph ( std::size_t p ): k(0), c(p) {

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

/// One graph with defined title

void PlotGraph::view (Int_t num, Double_t axisX[], Double_t axisY[]) {

	c[k] = new TCanvas (uName("c",k), uName("Graph_",k), x0+k*s, y0+k*s, w, h);
    c[k]->SetLogy(); // logarithmic scale

    TGraph *gr1 = new TGraph (num, axisX, axisY);

	gr1->Draw("AL");
    gr1->SetTitle(title);
    gr1->GetXaxis()->CenterTitle();
    gr1->GetYaxis()->CenterTitle();
	++k;
}

/// One graph with optional title

void PlotGraph::view (Int_t num, Double_t axisX[], Double_t axisY[], Char_t const* title) {

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

void PlotGraph::viewData (Int_t num, Double_t axisX[], Double_t axisY[]) {

	c[k] = new TCanvas (uName("c",k), uName("Graph_",k), x0+k*s, y0+k*s, w, h);
    c[k]->SetLogx(); // logarithmic scale
	
	TGraph *gr1 = new TGraph (num, axisX, axisY);
	gr1->Draw("AP");
	gr1->SetMarkerStyle(22);
	gr1->SetMarkerSize(.9);
	gr1->SetMarkerColor(9);
	++k;	
}

/// One graph with data points

void PlotGraph::viewPlusData (Int_t num, Double_t axisX[], Double_t axisY[], Int_t numD, Double_t axisXD[], Double_t axisYD[]) {

	c[k] = new TCanvas (uName("c",k), uName("Graph_",k), x0+k*s, y0+k*s, w, h);
    //c[k]->SetLogy(); // logarithmic scale

    TGraph *gr1 = new TGraph (num, axisX, axisY);
	gr1->Draw("AL");
    //gr1->SetTitle("Function with experimental points");
    gr1->SetTitle("");
    gr1->GetXaxis()->CenterTitle();
    gr1->GetYaxis()->CenterTitle();

    TGraph *gr2 = new TGraph (numD, axisXD, axisYD);
	gr2->Draw("P");
	gr2->SetMarkerStyle(21);
	gr2->SetMarkerSize(.5);
	++k;	
}

/// One graph with data points and title

void PlotGraph::viewPlusData (Int_t num, Double_t axisX[], Double_t axisY[], Int_t numD, Double_t axisXD[], Double_t axisYD[], Char_t const* title) {

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

/// One graph with data points and errors

void PlotGraph::viewPlusDataE (Int_t num, Double_t axisX[], Double_t axisY[], Int_t numD, Double_t axisXD[], Double_t axisYD[], Double_t axisXED[], Double_t axisYED[], const Char_t* title) {

	c[k] = new TCanvas (uName("c",k), uName("Graph_",k), x0+k*s, y0+k*s, w, h);
    c[k]->SetLogy(); // logarithmic scale

    TGraph *gr1 = new TGraph (num, axisX, axisY);
	gr1->Draw("AL");
	gr1->GetXaxis()->SetTitle("t [GeV^{2}]");
	gr1->GetYaxis()->SetTitle(title);
    //gr1->SetTitle(title);
    gr1->SetTitle();
    gr1->GetXaxis()->CenterTitle();
    gr1->GetYaxis()->CenterTitle();
	gr1->SetMaximum(2.);

    TGraphErrors *gr2 = new TGraphErrors (numD, axisXD, axisYD, axisXED, axisYED);
	gr2->Draw("P");
	gr2->SetMarkerStyle(21);
	gr2->SetMarkerSize(.7);
	gr2->SetMarkerColor(2);
	++k;	
}

