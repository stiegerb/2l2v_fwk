#ifndef _triganalysis_
#define _triganalysis_

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "UserCode/llvv_fwk/interface/DataEventSummaryHandler.h"
#include "UserCode/llvv_fwk/interface/MuScleFitCorrector.h"

#include "TTree.h"
#include "TFile.h"

class TrigAnalysis
{

public:
  TrigAnalysis(edm::ParameterSet, DataEventSummaryHandler*);

  void analyze(float);


  void initializeTreeVariables();
  void writeTree(TFile &f);

  bool isGoodElectron(data::PhysicsObject_t, float=10., float=2.5);
  bool isGoodMuon(data::PhysicsObject_t, float=10., float=2.1);
  bool isGoodJet(data::PhysicsObject_t,
                 data::PhysicsObjectCollection_t,
                 float=20., float=2.5);
  void selectLeptons(data::PhysicsObjectCollection_t,
                     data::PhysicsObjectCollection_t&);
  int selectJets(data::PhysicsObjectCollection_t,
                 data::PhysicsObjectCollection_t,
                 data::PhysicsObjectCollection_t&);

  TString getRecoChannel(data::PhysicsObjectCollection_t);
  TString getGenChannel();
  int channelToBin(TString);

 private:

  bool fIsMC;
  TString fJECDir;
  TString fURL;
  MuScleFitCorrector *fMuCor;

  DataEventSummaryHandler *fEvSummary_;
  DataEventSummary fCurrentEvent_;
  TTree* tree_;

  void clearTreeVariables();

  // Tree variables
  float fweight;
  int fnvtx;
  int fchannel;

  int nleps;
  int lepid[4];
  float leppt[4];
  float lepeta[4];
  float lepphi[4];

  int njets;
  float jetpt[10];
  float jeteta[10];
  float jetphi[10];
  float jetcsv[10];

};

#endif
