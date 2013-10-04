#ifndef _lxyanalysis_
#define _lxyanalysis_

#include "UserCode/llvv_fwk/interface/SmartSelectionMonitor.h"
#include "UserCode/llvv_fwk/interface/DataEventSummaryHandler.h"

/* #include "FWCore/ServiceRegistry/interface/Service.h" */
/* #include "CommonTools/UtilAlgos/interface/TFileService.h" */

#include "TTree.h"
#include "TFile.h"

class LxyAnalysis
{

public:
  LxyAnalysis(SmartSelectionMonitor &mon,bool runSystematics);

  void analyze(data::PhysicsObjectCollection_t & leptons, 
	       data::PhysicsObjectCollection_t &jets,
	       /* data::PhysicsObject_t &met,  */
	       ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > &met,
	       int nvtx,
	       data::PhysicsObjectCollection_t &mctruth,
	       float weight);


  void initializeTreeVariables();
  void writeTree(TFile &f);


 private:
  
  TTree* tree_;

  // tree variables
  // floats
  float flxy;
  float flxyerr;
  float fsecvtxmass;
  float fsecvtxntrk;
  float fweight;
  int fnvtx;
  int fchannel;
  
  /* int fcharge; */

  void clearTreeVariables();



  SmartSelectionMonitor *mon_;
  


};

#endif
