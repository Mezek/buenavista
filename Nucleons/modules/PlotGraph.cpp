/**
 * $Date$
 * $Revision$
 * $Author$
 * $HeadURL: http://triglav/repos/BuenaVista/Nucleons/modules/PlotGraph.cpp $
 * $Id$
 *
 * @file
 * @brief	Preview of functions and data.
 */

#include "PlotGraph.h"

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

/// Two separate graphs

void PlotGraph::view2 (Int_t num, Double_t axisX[], Double_t axisY1[], Double_t axisY2[]) {

	c[k] = new TCanvas (uName("c",k), uName("Graph_",k), x0+k*s, y0+k*s, 1.5*w, h);
	c[k]->SetLogy(); // logarithmic scale
	c[k]->Divide(2,1);

	c[k]->cd(1);
	gPad->SetLogy();
	TGraph *gr1 = new TGraph (num, axisX, axisY1);
	gr1->SetTitle("Graph 1");
	gr1->Draw("AL");
	
	c[k]->cd(2);
	gPad->SetLogy();
	TGraph *gr2 = new TGraph (num, axisX, axisY2);
	gr2->SetTitle("Graph 2");
	gr2->Draw("AL");
	++k;	
}

/// View data

void PlotGraph::viewData (Int_t num, Double_t axisX[], Double_t axisY[]) {

	c[k] = new TCanvas (uName("c",k), uName("Graph_",k), x0+k*s, y0+k*s, w, h);
	//c[k]->SetLogy(); // logarithmic scale
	
	TGraph *gr1 = new TGraph (num, axisX, axisY);
	gr1->Draw("AP");
	gr1->SetMarkerStyle(21);
	gr1->SetMarkerSize(.5);
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


/// One graph with data points and asymmetric errors

void PlotGraph::viewPlusDataAE (Int_t num, Double_t axisX[], Double_t axisY[], Int_t numD, Double_t axisXD[], Double_t axisYD[], Double_t axisXExl[], Double_t axisXExh[], Double_t axisYEyl[], Double_t axisYEyh[], Char_t const* title) {

	c[k] = new TCanvas (uName("c",k), uName("Graph_",k), x0+k*s, y0+k*s, w, h);
	c[k]->SetLogy(); // logarithmic scale

	TMultiGraph *mgr1 = new TMultiGraph();

	TGraph *gr1 = new TGraph (num, axisX, axisY);
	TGraphAsymmErrors *gr2 = new TGraphAsymmErrors (numD, axisXD, axisYD, axisXExl, axisXExh, axisYEyl, axisYEyh);
	gr2->SetMarkerStyle(21);
	gr2->SetMarkerSize(.7);
	gr2->SetMarkerColor(2);

	mgr1->Add(gr1,"L");
	mgr1->Add(gr2,"P");
	mgr1->Draw("A");

	mgr1->GetXaxis()->SetTitle("Q^2 [GeV^{2}]");
	mgr1->GetYaxis()->SetTitle(title);
	mgr1->GetXaxis()->CenterTitle();
	mgr1->GetYaxis()->CenterTitle();

	// Change the axis limits
	gPad->Modified();
	mgr1->GetXaxis()->SetLimits(-7.5,15.);
	mgr1->SetMinimum(0.01);
	mgr1->SetMaximum(10000.);

	TLegend *lg1 = new TLegend (.5,.7,.87,.87);
	lg1->SetFillColor(kWhite);
 	lg1->AddEntry(gr1,title,"l");
	lg1->AddEntry(gr2,"World data with errors","ep");
	lg1->SetTextSize(0.03);
	lg1->Draw();
	
	//c[k]->Print("viewPlusDataAE.eps");
	++k;
}

/// Two graphs in one canvas

void PlotGraph::viewPrint (Int_t num, Double_t axisX[], Double_t axisY1[], Double_t axisY2[]) {

	c[k] = new TCanvas (uName("c",k), uName("Graph_",k), x0+k*s, y0+k*s, w, h);
	//c[k]->SetLogy(); // logarithmic scale
	TGraph *gr1 = new TGraph (num, axisX, axisY1);
	gr1->SetLineStyle(kDashed);
	gr1->Draw("AL");
	
	TGraph *gr2 = new TGraph (num, axisX, axisY2);
	gr2->Draw("L");
  
	gr2->GetXaxis()->SetTitle("x");
	gr2->GetYaxis()->SetTitle("y");
	gr2->SetTitle("T");
	gr2->GetXaxis()->CenterTitle();
	gr2->GetYaxis()->CenterTitle();
	//gr2->SetMaximum(1);
	gr2->GetXaxis()->SetLimits(0.,1.00001);

	TLegend *lg1 = new TLegend (.7,.77,.87,.87);
	lg1->SetFillColor(kWhite);
 	lg1->AddEntry(gr1,"curve I.","l");
	lg1->AddEntry(gr2,"curve II.","l");
	lg1->Draw();
	++k;	
}

/// Four separate graphs

void PlotGraph::view4 (Int_t num, Double_t axisX[], Double_t axisY1[], Double_t axisY2[], Double_t axisY3[], Double_t axisY4[]) {

	c[k] = new TCanvas (uName("c",k), uName("Graph_",k), x0+k*s, y0+k*s, 1.5*w, h);
	c[k]->Divide(2,2);
	
	c[k]->cd(1);
	gPad->SetLogy();
	TGraph *gr1 = new TGraph (num, axisX, axisY1);
	gr1->SetTitle("Graph 1");
	gr1->SetLineColor(kBlue+2);
	gr1->Draw("AL");

	c[k]->cd(2);
	gPad->SetLogy();	
	TGraph *gr2 = new TGraph (num, axisX, axisY2);
	gr2->SetTitle("Graph 2");
	gr2->Draw("AL");
	
	c[k]->cd(3);
	gPad->SetLogy();
	TGraph *gr3 = new TGraph (num, axisX, axisY3);
	gr3->SetTitle("Graph 3");
	gr3->Draw("AL");

	c[k]->cd(4);
	gPad->SetLogy();
	TGraph *gr4 = new TGraph (num, axisX, axisY4);
	gr4->SetTitle("Graph 4");
	gr4->Draw("AL");
	
	++k;
}

/// Four separate graphs with experimental data

void PlotGraph::view4Exp (Int_t num, Double_t axisX[], Double_t axisY1[], Double_t axisY2[], Double_t axisY3[], Double_t axisY4[], Int_t num1, Double_t axisEX1[], Double_t axisEY1[], Int_t num2, Double_t axisEX2[], Double_t axisEY2[], Int_t num3, Double_t axisEX3[], Double_t axisEY3[], Int_t num4, Double_t axisEX4[], Double_t axisEY4[]) {

	c[k] = new TCanvas (uName("c",k), uName("Graph_",k), x0+k*s, y0+k*s, 2*w, h);
	c[k]->Divide(2,2);

	c[k]->cd(1);
	gPad->SetLogy();
	TGraph *gr1 = new TGraph (num, axisX, axisY1);
	gr1->SetTitle("|G_{E}^{p}|");
	gr1->Draw("AL");
	gr1->SetLineColor(kAzure);
	gr1->GetXaxis()->SetTitle("t [GeV^{2}]");	
	TGraph *gr2 = new TGraph (num1, axisEX1, axisEY1);
	//TGraphErrors *gr2 = new TGraphErrors (num1, axisEX1, axisEY1, axisXED, axisYED);
	gr2->Draw("P");
	gr2->SetMarkerStyle(21);
	gr2->SetMarkerSize(.5);
	gr2->SetMarkerColor(kBlue+2);

	c[k]->cd(2);
	gPad->SetLogy();	
	TGraph *gr3 = new TGraph (num, axisX, axisY2);
	gr3->SetTitle("|G_{M}^{p}|");
	gr3->Draw("AL");
	gr3->SetLineColor(kAzure);
	gr3->GetXaxis()->SetTitle("t [GeV^{2}]");	
	TGraph *gr4 = new TGraph (num2, axisEX2, axisEY2);
	gr4->Draw("P");
	gr4->SetMarkerStyle(21);
	gr4->SetMarkerSize(.5);
	gr4->SetMarkerColor(kBlue+2);
	
	c[k]->cd(3);
	gPad->SetLogy();
	TGraph *gr5 = new TGraph (num, axisX, axisY3);
	gr5->SetTitle("|G_{E}^{n}|");
	gr5->Draw("AL");
	gr5->SetLineColor(kAzure);
	gr5->GetXaxis()->SetTitle("t [GeV^{2}]");	
	TGraph *gr6 = new TGraph (num3, axisEX3, axisEY3);
	gr6->Draw("P");
	gr6->SetMarkerStyle(21);
	gr6->SetMarkerSize(.5);
	gr6->SetMarkerColor(kBlue+2);
	
	c[k]->cd(4);
	gPad->SetLogy();
	TGraph *gr7 = new TGraph (num, axisX, axisY4);
	gr7->SetTitle("|G_{M}^{n}|");
	gr7->Draw("AL");
	gr7->SetLineColor(kAzure);
	gr7->GetXaxis()->SetTitle("t [GeV^{2}]");	
	TGraph *gr8 = new TGraph (num4, axisEX4, axisEY4);
	gr8->Draw("P");
	gr8->SetMarkerStyle(21);
	gr8->SetMarkerSize(.5);
	gr8->SetMarkerColor(kBlue+2);
	
	++k;	
}
