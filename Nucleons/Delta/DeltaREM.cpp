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

void DeltaPlot::viewREM (Char_t const* title)
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

	double DX[nPoints], DY[nPoints];
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
		DX[i] = qA2;
		double mDip = 1. + qA2/0.71;
		double gD = 1./mDip/mDip;
		double gDmn = gD*(-muN);
		DY[i] = sqrt(2.)*2./3.*gDmn;

	}	
	
	c[k] = new TCanvas (uName("c",k), uName("Graph_",k), x0+k*s, y0+k*s, w, h);
	//c[k]->SetLogy(); // logarithmic scale

	static Int_t nG = 7;
	TGraphErrors *g[nG];
	TMultiGraph *mg = new TMultiGraph();
	Double_t EX0[100] = {0};

	// 1999-PRL-82-45_Frolov
	Double_t X1[] = {2.8, 4.0};
	Double_t Y1[] = {-2.0, -3.1};
	Double_t U1[] = {1.30, 1.30};
	Double_t D1[] = {1.30, 1.30};
	g[1] = new TGraphErrors (2, X1, Y1, EX0, D1);
	g[1]->SetTitle("JLab/Hall C");
	g[1]->SetMarkerColor(2);
	g[1]->SetMarkerStyle(21);
	mg->Add(g[1]);

	// 2006-PRL-97-112003_Ungaro
	Double_t X2[] = {3.0, 3.5, 4.2, 5.0, 6.0};
	Double_t Y2[] = {-1.61, -1.07, -3.15, -3.23, -3.84};
	Double_t U2[] = {0.4478, 0.4805, 0.7280, 1.5456, 3.0325};
	Double_t D2[] = {0.4478, 0.4805, 0.7280, 1.5456, 3.0325};
	g[2] = new TGraphErrors (5, X2, Y2, EX0, D2);
	g[2]->SetTitle("JLaB/CLAS");
	g[2]->SetMarkerColor(4);
	g[2]->SetMarkerStyle(22);
	mg->Add(g[2]);

	// 2002-PRL-88-122001_Joo
	Double_t X3[] = {0.40, 0.52, 0.65, 0.75, 0.90, 0.65, 0.75, 0.90, 1.15, 1.45, 1.80};
	Double_t Y3[] = {-3.4, -1.6, -1.9, -2.1, -1.8, -2.0, -1.6, -1.8, -1.6, -2.4, -0.9};
	Double_t U3[] = {0.5657, 0.5657, 0.7071, 0.9220, 0.7211, 0.5657, 0.7071, 0.5000, 0.5831, 0.8062, 1.3038};
	Double_t D3[] = {0.5657, 0.5657, 0.7071, 0.9220, 0.7211, 0.5657, 0.7071, 0.5000, 0.5831, 0.8062, 1.3038};
	g[3] = new TGraphErrors (11, X3, Y3, EX0, D3);
	g[3]->SetTitle("JLab/CLAS");
	g[3]->SetMarkerColor(2);
	g[3]->SetMarkerStyle(22);
	mg->Add(g[3]);

	// 2001-PRL-86-2963_Mertz
	// 2003-PL-B564-21_Kunz
	// 2005-PRL-94-022003_Sparveris
	Double_t X4[] = {0.126, 0.127, 0.127};
	Double_t Y4[] = {-2.1, -2.2, -2.3};
	Double_t U4[] = {2.01, 0.90, 0.67};
	Double_t D4[] = {2.01, 0.90, 0.67};
	g[4] = new TGraphErrors (3, X4, Y4, EX0, D4);
	g[4]->SetTitle("MIT-Bates");
	g[4]->SetMarkerColor(4);
	g[4]->SetMarkerStyle(23);
	mg->Add(g[4]);

	// 2000-PR-C61-035204_Beck
	// 2006-EPJ-A30-471_Stave
	// 2007-PL-B651-102_Sparveris
	Double_t X5[] = {0.000, 0.060, 0.200};
	Double_t Y5[] = {-2.50, -2.28, -1.96};
	Double_t U5[] = {0.2236, 0.3523, 0.7940};
	Double_t D5[] = {0.2236, 0.3523, 0.7940};
	g[5] = new TGraphErrors (3, X5, Y5, EX0, D5);
	g[5]->SetTitle("MAMI");
	g[5]->SetMarkerColor(4);
	g[5]->SetMarkerStyle(24);
	mg->Add(g[5]);

	// 2001-PRC-64-025203_Blanpied
	Double_t X6[] = {0.000};
	Double_t Y6[] = {-3.07};
	Double_t U6[] = {0.3538};
	Double_t D6[] = {0.3538};
	g[6] = new TGraphErrors (1, X6, Y6, EX0, D6);
	g[6]->SetTitle("LEGS");
	g[6]->SetMarkerColor(2);
	g[6]->SetMarkerStyle(20);
	mg->Add(g[6]);
	
	// Dipole formula
	TGraph *gD = new TGraph (nPoints, DX, DY);
	gD->SetTitle("Dipole formulae");
	gD->SetFillColor(0);
	mg->Add(gD);

	for (int i=1; i<nG; ++i)
	{
		g[i]->SetFillColor(0);
		//g[i]->SetLineColor(4);
		g[i]->SetMarkerSize(1.2);
	}
	mg->Draw("AP");

	//TF1 *fg = new TF1 ("fg", "[1]*x + [0]");
	//mg->Fit("poly5","Fit"); // fg

	TAxis *aX = mg->GetXaxis();
	aX->SetTitle("Q^{2} [GeV^{2}]");
	//aX->SetLimits(0.,10.);
	aX->SetTitleOffset(1.2);
	aX->CenterTitle();

	TAxis *aY = mg->GetYaxis();
	aY->SetTitle("R_{EM} [%]");
	//aY->SetRangeUser(0.6,1.8);
	aY->SetTitleOffset(1.2);
	aY->CenterTitle();

	gPad->SetFillColor(kWhite);
	
	TLegend *leg = c[k]->BuildLegend();
	leg->SetFillStyle(0);

	c[k]->Modified();

	c[k]->SaveAs("imgREM.pdf");
	
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
