/**
 * ============================================================================
 * $Date$
 * $Revision$
 * $Author$
 * $HeadURL: http://triglav/repos/BuenaVista/Nucleons/modules/PlotGraph.h $
 * $Id$
 *
 * @file
 * @brief	Header for plotting.
 */

#ifndef _DeltaPlot_H_
#define _DeltaPlot_H_

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

class DeltaPlot {
  private:
	Char_t const* title;
	std::vector<TCanvas *> c;
	int k;
	Double_t x0, y0, s, w, h;
  public:
	DeltaPlot ( std::size_t );
	void view (Int_t num, Double_t axisX[], Double_t axisY[], Char_t const*);
	void viewData (Int_t num, Double_t axisX[], Double_t axisY[], Char_t const*);
	void viewPlusData (Int_t num, Double_t axisX[], Double_t axisY[], Int_t numD, Double_t axisXD[], Double_t axisYD[], Char_t const*);
	void viewGMstar (Char_t const*);
	void viewREM (Char_t const*);
	void viewRSM (Char_t const*);
};

#endif // _DeltaPlot_H_
