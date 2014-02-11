/**
 * $Date$
 * $Revision$
 * $Author$
 * $Header$
 * $Id$
 *
 * @file
 * @brief	Header for plotting.
 */

#ifndef PlotGraph_H_
#define PlotGraph_H_

// Root graphics
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TMultiGraph.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TStyle.h"
#include "TColor.h"
#include "TAxis.h"
#include "TLegend.h"
#include "TLatex.h"

// unique name
const char* uName ( std::string name,  int i ) {
	std::ostringstream o;
	o << i;
	return name.append(o.str()).c_str();
}

class PlotGraph {
  private:
	Char_t const* title;
	std::vector<TCanvas *> c;
	int k;
	Double_t x0, y0, s, w, h;
  public:
    PlotGraph ( std::size_t );
    void view (Int_t num, Double_t axisX[], Double_t axisY[]);
    void view (Int_t num, Double_t axisX[], Double_t axisY[], Char_t const*);
    void view2 (Int_t num, Double_t axisX[], Double_t axisY1[], Double_t axisY2[]);
    void viewData (Int_t num, Double_t axisX[], Double_t axisY[]);
	void viewPlusData (Int_t num, Double_t axisX[], Double_t axisY[], Int_t numD, Double_t axisXD[], Double_t axisYD[]);
    void viewPlusData (Int_t num, Double_t axisX[], Double_t axisY[], Int_t numD, Double_t axisXD[], Double_t axisYD[], Char_t const*);
    void viewPlusDataE (Int_t num, Double_t axisX[], Double_t axisY[], Int_t numD, Double_t axisXD[], Double_t axisYD[], Double_t axisXED[], Double_t axisYED[], Char_t const*);
    void viewPlusDataAE (Int_t num, Double_t axisX[], Double_t axisY[], Int_t numD, Double_t axisXD[], Double_t axisYD[], Double_t axisXExl[],  Double_t axisXExh[], Double_t axisYEyl[], Double_t axisYEyh[], Char_t const*);
    void viewPrint (Int_t num, Double_t axisX[], Double_t axisY1[], Double_t axisY2[]);
    void view4 (Int_t num, Double_t axisX[], Double_t axisY1[], Double_t axisY2[], Double_t axisY3[], Double_t axisY4[]);
    void view4Exp (Int_t num, Double_t axisX[], Double_t axisY1[], Double_t axisY2[], Double_t axisY3[], Double_t axisY4[], Int_t num1, Double_t axisEX1[], Double_t axisEY1[], Int_t num2, Double_t axisEX2[], Double_t axisEY2[], Int_t num3, Double_t axisEX3[], Double_t axisEY3[], Int_t num4, Double_t axisEX4[], Double_t axisEY4[]);
};

#endif // PlotGraph_H_
