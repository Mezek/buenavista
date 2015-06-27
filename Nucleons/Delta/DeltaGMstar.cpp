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

using namespace ROOT::Minuit2;

void DeltaPlot::viewGMstar (Char_t const* title)
{
	FFactor GMS(12);
	GMS.LoadParameters(parametersFile);
	GMS.CheckParameters();
	//GMS.PrintParameters();

	const int nPoints = 2500;
	double qMin;
	double qMax;
	double qStep;
	double qA;

	double DX[nPoints], DY[nPoints], RY[nPoints];
	qMin = 0.004;
	qMax = 3.0;
	qStep = (qMax-qMin)/nPoints;
	qA = qMin;

	for (int i = 0; i < nPoints; i++) {
		qA = qMin + i*qStep;
		double qA2 = qA*qA;
		double gen = GMS.AbsGEN(-qA2);
		double gmn = GMS.AbsGMN(-qA2);
		double msq = massDi*massDi - massNucl*massNucl - qA2;
		double abq = sqrt((qA2 + msq*msq)/(4.*massDi*massDi));
		double jsc = (massDi+massNucl)/2./massNucl*(1.-qA2/4./massDi/massDi)*(1.-qA2/4./massNucl/massNucl);
		DX[i] = qA2;
		double mDip = 1. + qA2/0.71;
		double gD = 1./mDip/mDip;
		double gDmn = gD*(-muN);
		DY[i] = sqrt(2.)*2./3.*gDmn*jsc;
		double masst = qA2/(4.*massNucl*massNucl);
		double gmo = (GMS.ScalarOne(-qA2)-GMS.VectorOne(-qA2))/qA2+(GMS.ScalarTwo(-qA2)-GMS.VectorTwo(-qA2))/massN/massN/4.;
		RY[i] = sqrt(2.)*2./3.*gmn*jsc;
	}	
	
	c[k] = new TCanvas (uName("c",k), uName("Graph_",k), x0+k*s, y0+k*s, w, h);
	//c[k]->SetLogy(); // logarithmic scale

	TGraphErrors *g[5];
	TMultiGraph *mg = new TMultiGraph();
	Double_t EX0[100] = {0};

	// 1999-PRL-82-45_Frolov
	Double_t X1[] = {2.8, 4.0};
	Double_t Y1[] = {0.0859, 0.0402};
	Double_t U1[] = {0.0035, 0.0019};
	Double_t D1[] = {0.0035, 0.0019};
	g[1] = new TGraphErrors (2, X1, Y1, EX0, D1);
	g[1]->SetTitle("JLab/Hall C");
	g[1]->SetMarkerColor(2);
	g[1]->SetMarkerStyle(21);
	mg->Add(g[1]);

	// 2006-PRL-97-112003_Ungaro
	Double_t X2[] = {3.0, 3.5, 4.2, 5.0, 6.0};
	Double_t Y2[] = {0.0697, 0.0524, 0.0346, 0.0242, 0.0134};
	Double_t U2[] = {0.0010, 0.0011, 0.0012, 0.0014, 0.0014};
	Double_t D2[] = {0.0010, 0.0011, 0.0012, 0.0014, 0.0014};
	g[2] = new TGraphErrors (5, X2, Y2, EX0, D2);
	g[2]->SetTitle("JLaB/CLAS");
	g[2]->SetMarkerColor(4);
	g[2]->SetMarkerStyle(21);
	mg->Add(g[2]);

	// 1968-PL-28-148B_Bartel
	Double_t X3[] = {0.20, 0.30, 0.40, 0.47, 0.48, 0.50, 0.60, 0.63, 0.63, 0.77, 0.78, 0.79, 0.97, 0.98, 1.15, 1.34, 1.57, 2.34};
	Double_t Y3[] = {1.7700, 1.3800, 1.1700, 0.9780, 0.9610, 0.9640, 0.7660, 0.7350, 0.7190, 0.5700, 0.5720, 0.5530, 0.4460, 0.4460, 0.3260, 0.2690, 0.2090, 0.1020};
	Double_t U3[] = {0.0620, 0.0483, 0.0351, 0.0293, 0.0336, 0.0289, 0.0268, 0.0221, 0.0252, 0.0200, 0.0172, 0.0194, 0.0156, 0.0156, 0.0147, 0.0121, 0.0115, 0.0082};
	Double_t D3[] = {0.0620, 0.0483, 0.0351, 0.0293, 0.0336, 0.0289, 0.0268, 0.0221, 0.0252, 0.0200, 0.0172, 0.0194, 0.0156, 0.0156, 0.0147, 0.0121, 0.0115, 0.0082};
	g[3] = new TGraphErrors (18, X3, Y3, EX0, D3);
	g[3]->SetTitle("DESY");
	g[3]->SetMarkerColor(6);
	g[3]->SetMarkerStyle(22);
	mg->Add(g[3]);

	// 1975-PR-D12-1884_Stein
	Double_t X4[] = {0.09, 0.22, 0.46, 0.78, 1.17, 1.48, 1.82};
	Double_t Y4[] = {2.2448, 1.5824, 0.9147, 0.5007, 0.2708, 0.1728, 0.1122};
	Double_t U4[] = {0.0709, 0.0332, 0.0243, 0.0136, 0.0116, 0.0095, 0.0064};
	Double_t D4[] = {0.0709, 0.0332, 0.0243, 0.0136, 0.0116, 0.0095, 0.0064};
	g[4] = new TGraphErrors (7, X4, Y4, EX0, D4);
	g[4]->SetTitle("SLAC");
	g[4]->SetMarkerColor(4);
	g[4]->SetMarkerStyle(23);
	mg->Add(g[4]);

	for (int i=0; i<4; ++i)
	{
		g[i+1]->SetFillColor(0);
		//g[i+1]->SetLineColor(4);
		g[i+1]->SetMarkerSize(1.2);
	}

	// Formula D
	TGraph *gD = new TGraph (nPoints, DX, DY);
	gD->SetTitle("Dipole formulae");
	gD->SetFillColor(0);
	gD->SetLineWidth(3);
	gD->SetMarkerSize(0.3);
	gD->SetMarkerStyle(21);
	gD->SetMarkerColor(4);
	gD->SetLineColor(4);
	mg->Add(gD);
	
	// Formula A
	TGraph *gA = new TGraph (nPoints, DX, RY);
	gA->SetTitle("Our result");
	gA->SetFillColor(0);
	gA->SetLineWidth(3);
	gA->SetMarkerSize(0.3);
	gA->SetMarkerStyle(21);
	//gA->SetLineColor(3);
	mg->Add(gA);

	mg->Draw("AP");

	//TF1 *fg = new TF1 ("fg", "[1]*x + [0]");
	//mg->Fit("poly5","Fit"); // fg

	TAxis *aX = mg->GetXaxis();
	aX->SetTitle("Q^{2} [GeV^{2}]");
	//aX->SetLimits(0.,10.);
	aX->SetRangeUser(-0.01,7.0);
	aX->SetTitleOffset(1.2);
	aX->CenterTitle();

	TAxis *aY = mg->GetYaxis();
	aY->SetTitle("G_{M}^{*}");
	//aY->SetRangeUser(0.6,1.8);
	aY->SetTitleOffset(1.2);
	aY->CenterTitle();

	gPad->SetFillColor(kWhite);
	
	TLegend *leg = c[k]->BuildLegend();
	leg->SetFillStyle(0);

	c[k]->Modified();

	c[k]->SaveAs("imgGMstar.pdf");
	
	++k;	
}
