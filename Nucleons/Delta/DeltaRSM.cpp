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

void DeltaPlot::viewRSM (Char_t const* title)
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
	qMin = 0.03;
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
		DX[i] = qA2;
		double mDip = 1. + qA2/0.71;
		double gD = 1./mDip/mDip;
		double gDmn = gD*(-muN);
		//DY[i] = sqrt(2.)*2./3.*gDmn;
		double masst = qA2/(4.*massNucl*massNucl);
		double pa = 0.9;
		double pd = 1.75;
		DY[i] = 100.*pa*masst/(1.+pd*masst)*(-1.)*massNucl/qA/2.;
		//std::cout << i << " "<< qA2 << " " << gen << " " << gmn << std::endl;
		RY[i] = -100.*gen/gmn;
		//DY[i] = 100.*massNucl*massDi/qA2/2.*gen/gmn;
	}	
	
	c[k] = new TCanvas (uName("c",k), uName("Graph_",k), x0+k*s, y0+k*s, w, h);
	//c[k]->SetLogy(); // logarithmic scale

	static Int_t nG = 6;
	TGraphErrors *g[nG];
	TMultiGraph *mg = new TMultiGraph();
	Double_t EX0[100] = {0};

	// 1999-PRL-82-45_Frolov
	Double_t X1[] = {2.8, 4.0};
	Double_t Y1[] = {-11.2, -14.8};
	Double_t U1[] = {1.6401, 1.6401};
	Double_t D1[] = {1.6401, 1.6401};
	g[1] = new TGraphErrors (2, X1, Y1, EX0, D1);
	g[1]->SetTitle("JLab/Hall C");
	g[1]->SetMarkerColor(2);
	g[1]->SetMarkerStyle(21);
	mg->Add(g[1]);

	// 2006-PRL-97-112003_Ungaro
	Double_t X2[] = {3.0, 3.5, 4.2, 5.0, 6.0};
	Double_t Y2[] = {-11.5, -13.0, -16.4, -24.8, -24.8};
	Double_t U2[] = {2.0713, 1.3292, 1.8288, 3.8897, 6.0902};
	Double_t D2[] = {2.0713, 1.3292, 1.8288, 3.8897, 6.0902};
	g[2] = new TGraphErrors (5, X2, Y2, EX0, D2);
	g[2]->SetTitle("JLaB/CLAS");
	g[2]->SetMarkerColor(4);
	g[2]->SetMarkerStyle(22);
	mg->Add(g[2]);

	// 2002-PRL-88-122001_Joo
	Double_t X3[] = {0.40, 0.52, 0.65, 0.75, 0.90, 0.65, 0.75, 0.90, 1.15, 1.45, 1.80};
	Double_t Y3[] = {-5.6, -6.4, -6.9, -7.4, -8.4, -6.6, -6.0, -7.2, -7.9, -7.7, -11.6};
	Double_t U3[] = {0.7211, 0.6403, 0.7810, 0.9434, 0.9849, 0.4472, 0.4472, 0.4123, 0.6403, 1.1402, 2.1932};
	Double_t D3[] = {0.7211, 0.6403, 0.7810, 0.9434, 0.9849, 0.4472, 0.4472, 0.4123, 0.6403, 1.1402, 2.1932};
	g[3] = new TGraphErrors (11, X3, Y3, EX0, D3);
	g[3]->SetTitle("JLab/CLAS");
	g[3]->SetMarkerColor(2);
	g[3]->SetMarkerStyle(22);
	mg->Add(g[3]);

	// 2001-PRL-86-2963_Mertz
	// 2003-PL-B564-21_Kunz
	// 2005-PRL-94-022003_Sparveris
	Double_t X4[] = {0.126, 0.127, 0.127};
	Double_t Y4[] = {-6.5, -6.5, -6.1};
	Double_t U4[] = {2.5080, 0.3000, 0.5385};
	Double_t D4[] = {2.5080, 0.3000, 0.5385};
	g[4] = new TGraphErrors (3, X4, Y4, EX0, D4);
	g[4]->SetTitle("MIT-Bates");
	g[4]->SetMarkerColor(4);
	g[4]->SetMarkerStyle(23);
	mg->Add(g[4]);

	// 2001-PRL-86-2959_Pospischil
	// 2006-EPJ-A30-471_Stave
	// 2006-EPJ-A27-91_Elsner
	// 2007-PL-B651-102_Sparveris
	Double_t X5[] = {0.121, 0.060, 0.200, 0.200};
	Double_t Y5[] = {-6.40, -4.81, -5.45, -5.09};
	Double_t U5[] = {1.0630, 0.3748, 0.4200, 0.4104};
	Double_t D5[] = {1.0630, 0.3748, 0.4200, 0.4104};
	g[5] = new TGraphErrors (4, X5, Y5, EX0, D5);
	g[5]->SetTitle("MAMI");
	g[5]->SetMarkerColor(4);
	g[5]->SetMarkerStyle(24);
	mg->Add(g[5]);
	
	for (int i=1; i<nG; ++i)
	{
		g[i]->SetFillColor(0);
		//g[i]->SetLineColor(4);
		g[i]->SetMarkerSize(1.2);
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
	aX->SetTitleOffset(1.2);
	aX->CenterTitle();

	TAxis *aY = mg->GetYaxis();
	aY->SetTitle("R_{SM} [%]");
	//aY->SetRangeUser(0.6,1.8);
	aY->SetTitleOffset(1.2);
	aY->CenterTitle();

	gPad->SetFillColor(kWhite);
	
	TLegend *leg = c[k]->BuildLegend();
	leg->SetFillStyle(0);

	c[k]->Modified();

	c[k]->SaveAs("imgRSM.pdf");
	
	++k;	
}
