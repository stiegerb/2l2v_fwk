#include <iostream>
#include <boost/shared_ptr.hpp>
#include <fstream>

#include "UserCode/llvv_fwk/interface/MacroUtils.h"
#include "UserCode/llvv_fwk/interface/SmartSelectionMonitor.h"
#include "UserCode/llvv_fwk/interface/DataEventSummaryHandler.h"
#include "UserCode/llvv_fwk/interface/LxyAnalysis.h"
#include "UserCode/llvv_fwk/interface/UEAnalysis.h"
#include "UserCode/llvv_fwk/interface/BTVAnalysis.h"
#include "UserCode/llvv_fwk/interface/RAnalysis.h"
#include "UserCode/llvv_fwk/interface/TopPtWeighter.h"
#include "UserCode/llvv_fwk/interface/LeptonEfficiencySF.h"
#include "UserCode/llvv_fwk/interface/MuScleFitCorrector.h"

#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"

#include "FWCore/FWLite/interface/AutoLibraryLoader.h"
#include "FWCore/PythonParameterSet/interface/MakeParameterSets.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"

#include "PhysicsTools/CondLiteIO/interface/RecordWriter.h"
#include "DataFormats/FWLite/interface/Record.h"
#include "DataFormats/FWLite/interface/EventSetup.h"
#include "DataFormats/FWLite/interface/ESHandle.h"
#include "CondFormats/PhysicsToolsObjects/interface/BinningPointByMap.h"
#include "RecoBTag/PerformanceDB/interface/BtagPerformance.h"
#include "PhysicsTools/FWLite/interface/CommandLineParser.h"

#include "TSystem.h"
#include "TFile.h"
#include "TChain.h"
#include "TChainElement.h"
#include "TCanvas.h"
#include "TString.h"
#include "TDirectory.h"
#include "TEventList.h"
#include "TRandom.h"

#include <iostream>

using namespace std;

//
int main(int argc, char* argv[])
{
	// load framework libraries
	gSystem->Load( "libFWCoreFWLite" );
	AutoLibraryLoader::enable();

	//check arguments
	if ( argc < 2 ) {
		std::cout << "Usage : " << argv[0] << " parameters_cfg.py" << std::endl;
		return 0;
	}

	//
	// configure
	//
	const edm::ParameterSet &runProcess = edm::readPSetsFrom(argv[1])->getParameter<edm::ParameterSet>("runProcess");
	std::vector<std::string> urls=runProcess.getParameter<std::vector<std::string> >("input");
	TString url = TString(urls[0]);
	TString baseDir     = runProcess.getParameter<std::string>("dirName");
	bool runSystematics = runProcess.getParameter<bool>("runSystematics");
	TString jecDir      = runProcess.getParameter<std::string>("jecDir");
	bool isMC           = runProcess.getParameter<bool>("isMC");
	int mcTruthMode     = runProcess.getParameter<int>("mctruthmode");
	double xsec         = runProcess.getParameter<double>("xsec");
	bool isV0JetsMC(isMC && (url.Contains("DYJetsToLL_50toInf") || url.Contains("WJets")));
	TString out          = runProcess.getParameter<std::string>("outdir");
	bool saveSummaryTree = runProcess.getParameter<bool>("saveSummaryTree");
	std::vector<string>  weightsFile = runProcess.getParameter<std::vector<string> >("weightsFile");


	//jet energy scale uncertainties
	gSystem->ExpandPathName(jecDir);
	FactorizedJetCorrector *jesCor        = utils::cmssw::getJetCorrector(jecDir,isMC);
	JetCorrectionUncertainty *totalJESUnc = new JetCorrectionUncertainty((jecDir+"/MC_Uncertainty_AK5PFchs.txt").Data());

	//muon energy scale and uncertainties
	MuScleFitCorrector *muCor=getMuonCorrector(jecDir,url);

	//
	// read the background scale factors
	//

	//re-scale w+jets
	std::map<TString,float> wjetsSFmap;
	if(weightsFile.size() && (url.Contains("WJets") || url.Contains("W1Jets") || url.Contains("W2Jets") || url.Contains("W3Jets") || url.Contains("W4Jets")) && !(url.Contains("TTWJets"))  && isMC){
		cout << "Reading W+jets scale factors..." << endl;
		TFile *wjetsF=TFile::Open(weightsFile[0].c_str());
		TH1* wjetssfH=(TH1 *)wjetsF->Get("wjetssf");
		for(int ibin=1; ibin<=wjetssfH->GetXaxis()->GetNbins(); ibin++) {
			wjetsSFmap[wjetssfH->GetXaxis()->GetBinLabel(ibin)]=wjetssfH->GetBinContent(ibin);
			cout << wjetssfH->GetXaxis()->GetBinLabel(ibin) << " " << wjetssfH->GetBinContent(ibin) << endl;
		}
		wjetsF->Close();
	}

	//re-scale qcd
	std::map<TString,float> qcdSFmap;
	if(weightsFile.size() && url.Contains("QCD") && isMC){
		TFile *qcdF=TFile::Open(weightsFile[0].c_str());
		TH1* qcdsfH=(TH1 *)qcdF->Get("qcdsf");
		for(int ibin=1; ibin<=qcdsfH->GetXaxis()->GetNbins(); ibin++) {
			qcdSFmap[qcdsfH->GetXaxis()->GetBinLabel(ibin)]=qcdsfH->GetBinContent(ibin);
			cout << qcdsfH->GetXaxis()->GetBinLabel(ibin) << " " << qcdsfH->GetBinContent(ibin) << endl;
		}
		qcdF->Close();
	}


	//
	// check input file
	//
	TFile *inF = TFile::Open(url);
	if(inF==0) return -1;
	if(inF->IsZombie()) return -1;
	TString proctag=gSystem->BaseName(url);
	Ssiz_t pos=proctag.Index(".root");
	proctag.Remove(pos,proctag.Length());

	// only filter electrons or muons
	bool filterOnlyE(false), filterOnlyMU(false);
	if(!isMC){
		if(url.Contains("SingleEle")) filterOnlyE=true;
		if(url.Contains("SingleMu"))  filterOnlyMU=true;
	}

	//
	// pileup reweighter
	//
	std::vector<double> dataPileupDistributionDouble = runProcess.getParameter< std::vector<double> >("datapileup");
	std::vector<float> dataPileupDistribution; for(unsigned int i=0;i<dataPileupDistributionDouble.size();i++){dataPileupDistribution.push_back(dataPileupDistributionDouble[i]);}
	std::vector<float> mcPileupDistribution;
	if(isMC){
		TString puDist(baseDir+"/pileup");
		TH1F* histo = (TH1F *) inF->Get(puDist);
		if(!histo)std::cout<<"pileup histogram is null!!!\n";
		for(int i=1;i<=histo->GetNbinsX();i++) mcPileupDistribution.push_back(histo->GetBinContent(i));
		delete histo;
	}
	while(mcPileupDistribution.size()<dataPileupDistribution.size())   mcPileupDistribution.push_back(0.0);
	while(mcPileupDistribution.size()>dataPileupDistribution.size()) dataPileupDistribution.push_back(0.0);
	gROOT->cd();  //THIS LINE IS NEEDED TO MAKE SURE THAT HISTOGRAM INTERNALLY PRODUCED IN LumiReWeighting ARE NOT DESTROYED WHEN CLOSING THE FILE
	edm::LumiReWeighting *LumiWeights= isMC ? new edm::LumiReWeighting(mcPileupDistribution,dataPileupDistribution): 0;
	utils::cmssw::PuShifter_t PuShifters;
	if(isMC) PuShifters=utils::cmssw::getPUshifters(dataPileupDistribution,0.05);

	//systematic variations for final selection
	std::vector<TString> systVars(1,"");
	if(isMC && runSystematics){
		systVars.push_back("jerdown");
		systVars.push_back("jerup");
		systVars.push_back("jesdown");
		systVars.push_back("jesup");
		systVars.push_back("umetdown");
		systVars.push_back("umetup");
		systVars.push_back("pudown");
		systVars.push_back("puup");
	}

	//
	// control histograms
	//
	SmartSelectionMonitor controlHistos;
	TH1F* Hhepup        = (TH1F*) controlHistos.addHistogram(new TH1F ("heupnup"    , "hepupnup"    ,20,0,20) ) ;
	TH1F* Hcutflow      = (TH1F*) controlHistos.addHistogram(new TH1F ("cutflow"    , "cutflow"    ,4,0,4) ) ;
	TH1F* Hoptim_systs  = (TH1F*) controlHistos.addHistogram(new TH1F ("optim_systs"    , ";syst;", systVars.size(),0,systVars.size()) );

	//vertex multiplicity
	controlHistos.addHistogram( new TH1F ("nvertices", "; Vertex multiplicity; Events", 50, 0.,50.) );

	//event selection histogram
	TString labels[]={"1 lepton", "loose lepton veto", "#geq 4 jets"};
	int nsteps=sizeof(labels)/sizeof(TString);
	TH1F *cutflowH = (TH1F*) controlHistos.addHistogram( new TH1F("evtflow",";Cutflow;Events",nsteps,0,nsteps) );
	for(int ibin=0; ibin<nsteps; ibin++) cutflowH->GetXaxis()->SetBinLabel(ibin+1,labels[ibin]);

	TString lxylabels[]={"1 lepton", "loose lepton veto", "#geq 3 jets", "#geq 1-btag"};
	int lxynsteps=sizeof(lxylabels)/sizeof(TString);
	TH1F *lxycutflowH = (TH1F*) controlHistos.addHistogram( new TH1F("lxyevtflow",";Cutflow;Events",lxynsteps,0,lxynsteps) );
	for(int ibin=0; ibin<lxynsteps; ibin++) lxycutflowH->GetXaxis()->SetBinLabel(ibin+1,lxylabels[ibin]);

	TString uelabels[]={"1 lepton", "loose lepton veto", "#geq 4 jets", "#geq 2-btags"};
	int nuesteps=sizeof(uelabels)/sizeof(TString);
	TH1F *uecutflowH = (TH1F*) controlHistos.addHistogram( new TH1F("ueevtflow",";Cutflow;Events",nuesteps,0,nuesteps) );
	for(int ibin=0; ibin<nuesteps; ibin++) uecutflowH->GetXaxis()->SetBinLabel(ibin+1,uelabels[ibin]);

	for(size_t ivar=0;ivar<systVars.size(); ivar++) {
		Hoptim_systs->GetXaxis()->SetBinLabel(ivar+1,systVars[ivar]);
		TH1F *finalCutflowH=new TH1F("finalevtflow"+systVars[ivar],";Category;Events",4,0,4);
		finalCutflowH->GetXaxis()->SetBinLabel(1,"=1 jets");
		finalCutflowH->GetXaxis()->SetBinLabel(2,"=2 jets");
		finalCutflowH->GetXaxis()->SetBinLabel(3,"=3 jets");
		finalCutflowH->GetXaxis()->SetBinLabel(4,"=4 jets");
		controlHistos.addHistogram( finalCutflowH );
	}

	// FIXME: control regions for the l+jet analysis????
	TString ctrlCats[]={"presel_","eq1jet_","geq3jet1btag_","geq4jet_","btag_","bveto_"};
	for(size_t k=0;k<sizeof(ctrlCats)/sizeof(TString); k++) {
		controlHistos.addHistogram( new TH1F(ctrlCats[k]+"charge",";lepton charge;Events",3,-1.5,1.5) );
		controlHistos.addHistogram( new TH1F(ctrlCats[k]+"emva", "; e-id MVA; Electrons", 50, 0.95,1.0) );
		controlHistos.addHistogram( new TH1F(ctrlCats[k]+"met",";Missing transverse energy [GeV];Events",50,0,500) );
		controlHistos.addHistogram( new TH1F(ctrlCats[k]+"iso",";Isolation variable;Events",50,0,1.0) );
		TH1F *h=(TH1F *)controlHistos.addHistogram( new TH1F(ctrlCats[k]+"njets",";Jet multiplicity;Events",6,0,6) );
		for(int ibin=1; ibin<=h->GetXaxis()->GetNbins(); ibin++) {
			TString label( ibin==h->GetXaxis()->GetNbins() ? "#geq" : "=");
			label += (ibin-1);
			label += " jets";
			h->GetXaxis()->SetBinLabel(ibin,label);

			if(ibin==1) continue;
			label="jet"; label+=(ibin-1);
			Float_t jetPtaxis[]={30,35,40,45,50,55,60,65,70,80,90,100,125,150,200,250,500};
			const size_t nJetPtBins=sizeof(jetPtaxis)/sizeof(Float_t)-1;
			controlHistos.addHistogram( new TH1F(ctrlCats[k]+"btag_"+label+"_pt",";Transverse momentum [GeV];Events",nJetPtBins,jetPtaxis) );
			controlHistos.addHistogram( new TH1F(ctrlCats[k]+"btag_"+label+"_eta",";Pseudo-rapidity;Events",25,0,2.5) );
			TH1F *flavH=(TH1F *)controlHistos.addHistogram( new TH1F(ctrlCats[k]+"btag_"+label+"_flav",";Flavor;Events",5,0,5) );
			for(int ibin=1; ibin<=5; ibin++) {
				TString label("unmatched");
				if(ibin==2) label="g";
				if(ibin==3) label="uds";
				if(ibin==4) label="c";
				if(ibin==5) label="b";
				flavH->GetXaxis()->SetBinLabel(ibin,label);
			}
		}
	}


	//lepton efficiencies
	LeptonEfficiencySF lepEff;

	// run one of the analyses...
	// UEAnalysis ueAn(controlHistos);
	// TNtuple *ueNtuple = ueAn.getSummaryTuple();

	//BTVAnalysis btvAn(controlHistos,runSystematics);
	LxyAnalysis lxyAn(controlHistos,runSystematics);
	lxyAn.initializeTreeVariables();


	///
	// process events file
	//
	DataEventSummaryHandler evSummary;
	if( !evSummary.attach( (TTree *) inF->Get(baseDir+"/data") ) )  { inF->Close();  return -1; }
	// const Int_t totalEntries=evSummary.getEntries();
	const Int_t totalEntries=10000;

	float cnorm=1.0;
	if(isMC){
		TH1F* cutflowH = (TH1F *) inF->Get(baseDir+"/cutflow");
		if(cutflowH) cnorm=cutflowH->GetBinContent(1);
	}
	Hcutflow->SetBinContent(1,cnorm);

	cout << "Processing: " << proctag << " @ " << url << endl
	<< "Initial number of events: " << cnorm << endl
	<< "Events in tree:           " << totalEntries << endl
	<< " xSec x BR:               " << xsec << endl;

	//check if a summary should be saved
	Float_t evSummaryWeight(1.0);
	Float_t xsecWeight(isMC ? xsec/cnorm : 1.0);
	TFile *spyFile=0;
	TDirectory *spyDir=0;
	DataEventSummaryHandler *spyEvents=0;
	if(saveSummaryTree){
		gSystem->Exec("mkdir -p " + out);
		gDirectory->SaveSelf();
		TString summaryName(out + "/" + proctag);
		if(mcTruthMode!=0) { summaryName += "_filt"; summaryName += mcTruthMode; }
		summaryName += "_summary.root";
		gSystem->ExpandPathName(summaryName);
		cout << "Creating event summary file @ " << summaryName << endl;

		//open file
		spyEvents = new DataEventSummaryHandler;
		spyFile = TFile::Open(summaryName,"RECREATE");
		spyFile->rmdir(proctag);
		spyDir = spyFile->mkdir("dataAnalyzer");
		TTree *outT = evSummary.getTree()->CloneTree(0);
		outT->SetTitle("Event summary");
		outT->SetDirectory(spyDir);
		outT->SetAutoSave(1000000);
		outT->Branch("weight",&evSummaryWeight,"weight/F");
		spyEvents->init(outT,false);

		//add also other summary tuples
		// ueNtuple->SetDirectory(spyDir);
	}



	//
	// Event loop
	//
	DuplicatesChecker duplicatesChecker;
	int nDuplicates(0);
	for (int inum=0; inum < totalEntries; ++inum){
		if(inum%500==0) { printf("\r [ %d/100 ]",int(100*float(inum)/float(totalEntries))); cout << flush; }
		evSummary.getEntry(inum);
		DataEventSummary &ev = evSummary.getEvent();
		if(!isMC && duplicatesChecker.isDuplicate( ev.run, ev.lumi, ev.event) ) { nDuplicates++; continue; }

		//pileup weight
		float weightNom(1.0),weightUp(1.0), weightDown(1.0);
		if(LumiWeights) {
			weightNom  = LumiWeights->weight(ev.ngenITpu);
			weightUp   = weightNom*PuShifters[utils::cmssw::PUUP]->Eval(ev.ngenITpu);
			weightDown = weightNom*PuShifters[utils::cmssw::PUDOWN]->Eval(ev.ngenITpu);
		}

		if(isV0JetsMC && ev.nup>5) continue;
		Hhepup->Fill(ev.nup,1);

		//MC truth
		data::PhysicsObjectCollection_t gen=evSummary.getPhysicsObject(DataEventSummaryHandler::GENPARTICLES);
		bool hasTop(false);
		int ngenLeptonsStatus3(0);
		if(isMC){
			for(size_t igen=0; igen<gen.size(); igen++){
				if(gen[igen].get("status")!=3) continue;
				int absid=abs(gen[igen].get("id"));
				if(absid==6) hasTop=true;
				if(absid!=11 && absid!=13 && absid!=15) continue;
				ngenLeptonsStatus3++;
			}

			if(mcTruthMode==1 && (ngenLeptonsStatus3!=1 || !hasTop)) continue; // mcTruthMode 1 will select events with exactly 1 status 3 lepton and a top
			if(mcTruthMode==2 && (ngenLeptonsStatus3==0 || !hasTop)) continue; // mcTruthMode 2 will select events with exactly 0 status 3 lepton and a top
		}

		Hcutflow->Fill(1,1);
		Hcutflow->Fill(2,weightNom);
		Hcutflow->Fill(3,weightUp);
		Hcutflow->Fill(4,weightDown);


		//trigger bits
		bool eTrigger   = ev.t_bits[13] || ev.t_bits[14];
		bool muTrigger  = ev.t_bits[6];
		if(filterOnlyE)   { muTrigger=false; }
		if(filterOnlyMU)  { eTrigger=false;  }


		//
		// Lepton selection
		//
		data::PhysicsObjectCollection_t leptons=evSummary.getPhysicsObject(DataEventSummaryHandler::LEPTONS);
		data::PhysicsObjectCollection_t selLeptons;
		data::PhysicsObjectCollection_t vetoLeptons;
		data::PhysicsObjectCollection_t antiIsoLeptons;
		for(size_t ilep=0; ilep<leptons.size(); ilep++){
			Int_t id=leptons[ilep].get("id");
			bool passKin(true),passId(true),passIso(true);             // variables for lepton selection
			bool passVetoKin(true),passVetoId(true),passVetoIso(true); // variables for additional lepton veto
			bool passAntiIso(true);                                    // variables for QCD control region selection
			if(abs(id)==11) {
				float sceta=leptons[ilep].getVal("sceta");
				Float_t gIso    = leptons[ilep].getVal("gIso03");
				Float_t chIso   = leptons[ilep].getVal("chIso03");
				//Float_t puchIso = leptons[ilep].getVal("puchIso03");
				Float_t nhIso   = leptons[ilep].getVal("nhIso03");
				float relIso=(TMath::Max(nhIso+gIso-ev.rho*utils::cmssw::getEffectiveArea(11,sceta),Float_t(0.))+chIso)/leptons[ilep].pt();
				// selection of veto electrons
				if(leptons[ilep].pt()<20)                      passVetoKin=false;
				if(fabs(leptons[ilep].eta())>2.5)              passVetoKin=false;
				if(leptons[ilep].getVal("tk_d0")>0.4)          passVetoId =false;  // FIXME: is this the correct value?
				if(leptons[ilep].getVal("mvatrig")<0.0)        passVetoId =false;
				if(relIso>0.15)                                passVetoIso=false;
				// selection of tight electrons
				if(leptons[ilep].pt()<30)                      passKin=false;
				if(fabs(leptons[ilep].eta())>2.5)              passKin=false;
				if(fabs(sceta)>1.4442 && fabs(sceta)<1.5660)   passKin=false;
				if(leptons[ilep].getFlag("isconv"))            passId =false;
				if(leptons[ilep].getVal("tk_d0")>0.2)          passId =false;  // FIXME: is this the correct value?
				if(leptons[ilep].getVal("tk_lostInnerHits")>0) passId =false;
				if(leptons[ilep].getVal("mvatrig")<0.5)        passId =false;
				if(relIso>0.1)                                 passIso=false;
				// selection of the QCD control region
				if(relIso<0.20)                                passAntiIso=false;

			}
			else if (abs(id)==13) {
				if(muCor) {
					TLorentzVector p4(leptons[ilep].px(),leptons[ilep].py(),leptons[ilep].pz(),leptons[ilep].energy());
					muCor->applyPtCorrection(p4 , id<0 ? -1 : 1 );
					if(isMC) muCor->applyPtSmearing(p4, id<0 ? -1 : 1, false);
					leptons[ilep].SetPxPyPzE(p4.Px(),p4.Py(),p4.Pz(),p4.E());
				}
				Int_t idbits    = leptons[ilep].get("idbits");
				bool isTight    = ((idbits >> 10) & 0x1);
				Float_t gIso    = leptons[ilep].getVal("gIso04");
				Float_t chIso   = leptons[ilep].getVal("chIso04");
				Float_t puchIso = leptons[ilep].getVal("puchIso04");
				Float_t nhIso   = leptons[ilep].getVal("nhIso04");
				Float_t relIso=(TMath::Max(nhIso+gIso-0.5*puchIso,0.)+chIso)/leptons[ilep].pt();
				// selection of veto muons (same selection as dilepton muon, except the lower pt)
				if(leptons[ilep].pt()<10)                      passVetoKin=false;
				if(fabs(leptons[ilep].eta())>2.5)              passVetoKin=false;
				// if(!isTight)                                   passVetoId=false;
				if(relIso>0.2)                                 passVetoIso=false;
				// selection of tight muons
				if(leptons[ilep].pt()<26)                                passKin=false;
				if(fabs(leptons[ilep].eta())>2.1)                        passKin=false;
				// FIXME: add the tight selection variables
				// if(leptons[ilep].getVal("innerTrackChi2")>10.)           passId=false;  // FIXME: is this the correct variable?
				// if(leptons[ilep].getVal("trkLayersWithMeasurement")<5.)  passId=false;
				// if(leptons[ilep].getVal("validMuonHits")<=0)             passId=false;
				// if(leptons[ilep].getVal("tk_d0")>0.2)                    passId=false;  // FIXME: is this the correct value?
				// if(leptons[ilep].getVal("tk_dz")>0.5)                    passId=false;  // FIXME: is this the correct value?
				// if(leptons[ilep].getVal("tk_validPixelHits")<=0)         passId=false;
				// if(leptons[ilep].getVal("mn_nMatchedStations")<1)        passId=false;
				if(!isTight)                                             passId=false;
				if(relIso>0.12)                                          passIso=false;
				// selection of the QCD control region
				if(relIso<0.20)                                          passAntiIso=false;
			}

			// if(!passVetoKin || !passVetoId || !(passVetoIso || passAntiIso)) continue;
			if(!passVetoKin || !passVetoId || !passVetoIso) continue;

			// veto leptons
			if((passVetoKin && passVetoId && passVetoIso) && !(passKin && passId && passIso)) {
				vetoLeptons.push_back(leptons[ilep]);
			}

			// if(!passKin || !passId || !(passIso || passAntiIso)) continue;
			if(!passKin || !passId || !passIso) continue;
			selLeptons.push_back(leptons[ilep]);
		}
		sort(selLeptons.begin(),selLeptons.end(),data::PhysicsObject_t::sortByPt);
		sort(vetoLeptons.begin(),vetoLeptons.end(),data::PhysicsObject_t::sortByPt);


		// Trigger selection
		if(!eTrigger && !muTrigger) continue;


		// At least one selected lepton
		if(selLeptons.size()<1) continue;


		// Determine the lepton+jets channel
		ev.cat=1;
		float lScaleFactor(1.0);
		ev.cat *= selLeptons[0].get("id");
		int id(abs(selLeptons[0].get("id")));
		lScaleFactor = isMC ? lepEff.getLeptonEfficiency( selLeptons[0].pt(), selLeptons[0].eta(), id,  id ==11 ? "loose" : "tight" ).first : 1.0; // FIXME: check!

		TString chName;
		if     (abs(ev.cat)==11 && eTrigger)   { chName="e";  }
		else if(abs(ev.cat)==13 && muTrigger)  { chName="mu"; }
		else                                   continue; // don't want e events triggered by mu, and vice versa
		std::vector<TString> ch(1,chName);

		// FIXME: what about the hardcoded numbers for the lepton scale factors below?

		/*
		//determine the dilepton channel
		ev.cat=1;
		float llScaleFactor(1.0);
		for(size_t ilep=0; ilep<2; ilep++)
		{
		ev.cat *= selLeptons[ilep].get("id");
		int id(abs(selLeptons[ilep].get("id")));
		llScaleFactor *= isMC ? lepEff.getLeptonEfficiency( selLeptons[ilep].pt(), selLeptons[ilep].eta(), id,  id ==11 ? "loose" : "tight" ).first : 1.0;
		}

		TString chName;
		bool isOS(ev.cat<0);
		bool isSameFlavor(false);
		if     (abs(ev.cat)==11*11 && eeTrigger)   { chName="ee";  isSameFlavor=true;  if(ngenLeptonsStatus3>=2) llScaleFactor*=0.972; }
		else if(abs(ev.cat)==11*13 && emuTrigger)  { chName="emu";                     if(ngenLeptonsStatus3>=2) llScaleFactor*=0.968; }
		else if(abs(ev.cat)==13*13 && mumuTrigger) { chName="mumu"; isSameFlavor=true; if(ngenLeptonsStatus3>=2) llScaleFactor*=0.955; }
		else                                       continue;
		std::vector<TString> ch(1,chName);
		if(isSameFlavor) ch.push_back("ll");
		*/

		//apply data/mc correction factors and update the event weight
		float weight = weightNom*lScaleFactor;

		//jet/met
		data::PhysicsObjectCollection_t recoMet = evSummary.getPhysicsObject(DataEventSummaryHandler::MET);
		data::PhysicsObjectCollection_t jets    = evSummary.getPhysicsObject(DataEventSummaryHandler::JETS);
		utils::cmssw::updateJEC(jets,jesCor,totalJESUnc,ev.rho,ev.nvtx,isMC);
		std::vector<LorentzVector> met = utils::cmssw::getMETvariations(recoMet[0],jets,selLeptons,isMC);

		//
		// Jet selection
		//
		data::PhysicsObjectCollection_t looseJets, selJets;
		int nbtags(0);
		for(size_t ijet=0; ijet<jets.size(); ijet++){
			//cross-clean with selected leptons
			double minDRlj(9999.);
			for(size_t ilep=0; ilep<selLeptons.size(); ilep++)
				minDRlj = TMath::Min( minDRlj, deltaR(jets[ijet],selLeptons[ilep]) );
			if(minDRlj<0.4) continue;

			//cross-clean with vetoed leptons //FIXME: should not matter as we anyway veto events with additional leptons...
			minDRlj = 9999.;
			for(size_t ilep=0; ilep<vetoLeptons.size(); ilep++) minDRlj = TMath::Min( minDRlj, deltaR(jets[ijet],vetoLeptons[ilep]) );
			if(minDRlj<0.4) continue;

			//require to pass the loose id
			Int_t idbits=jets[ijet].get("idbits");
			bool passPFloose( ((idbits>>0) & 0x1));
			if(!passPFloose) continue;

			//top candidate jets
			looseJets.push_back(jets[ijet]);
			if(jets[ijet].pt()<30 || fabs(jets[ijet].eta())>2.5 ) continue;
			selJets.push_back(jets[ijet]);
			nbtags += (jets[ijet].getVal("csv")>0.405);
		}
		sort(looseJets.begin(),looseJets.end(),data::PhysicsObject_t::sortByCSV);
		sort(selJets.begin(),  selJets.end(),  data::PhysicsObject_t::sortByCSV);


		//
		// EVENT SELECTION
		//
		if(selLeptons.size()<1) continue; // Redundant??
		controlHistos.fillHisto("evtflow",    ch, 0, weight);
		controlHistos.fillHisto("lxyevtflow", ch, 0, weight);
		controlHistos.fillHisto("ueevtflow",  ch, 0, weight);
		controlHistos.fillHisto("nvertices",  ch, ev.nvtx, weight);

		bool passSoftLeptonVeto(vetoLeptons.size()==0 && selLeptons.size()==1 ); // exactly one selected and zero vetoed leptons
		bool pass3JetSelection (selJets.size()>=3 );                             // three or more selected jets
		bool passJetSelection  (selJets.size()>=4 );                             // four or more selected jets
		bool pass1BtagSelection(nbtags>=1 );                                     // one or more b-tagged jets
		bool passBtagSelection (nbtags>=2 );                                     // two or more b-tagged jets
		bool passBtagVeto      (nbtags==0 );                                     // no b-tagged jets

		//
		// NOMINAL SELECTION CONTROL
		//
		std::vector<TString> ctrlCategs;
		float wjetsWeight(1.0),qcdWeight(1.0),ibtagdyWeight(1.0);
		if(passSoftLeptonVeto)                                             { ctrlCategs.push_back("presel_");         if(wjetsSFmap.find(chName+"_presel")!=wjetsSFmap.end())       wjetsWeight=wjetsSFmap[chName+"_presel"];}
		if(passSoftLeptonVeto && selJets.size()==1 )                       { ctrlCategs.push_back("eq1jet_");         if(wjetsSFmap.find(chName+"_eq1jet")!=wjetsSFmap.end())       wjetsWeight=wjetsSFmap[chName+"_eq1jet"];}
		if(passSoftLeptonVeto && pass3JetSelection && pass1BtagSelection)  { ctrlCategs.push_back("geq3jet1btag_");   if(wjetsSFmap.find(chName+"_geq3jet1btag")!=wjetsSFmap.end()) wjetsWeight=wjetsSFmap[chName+"_geq3jet1btag"];}
		if(passSoftLeptonVeto && passJetSelection )                        { ctrlCategs.push_back("geq4jet_");        if(wjetsSFmap.find(chName+"_geq4jet")!=wjetsSFmap.end())      wjetsWeight=wjetsSFmap[chName+"_geq4jet"];}
		if(passSoftLeptonVeto && passJetSelection  && passBtagSelection)   { ctrlCategs.push_back("btag_");           if(wjetsSFmap.find(chName+"_btag")!=wjetsSFmap.end())         wjetsWeight=wjetsSFmap[chName+"_btag"];}
		if(passSoftLeptonVeto && passJetSelection  && passBtagVeto)        { ctrlCategs.push_back("bveto_");          if(wjetsSFmap.find(chName+"_bveto")!=wjetsSFmap.end())        wjetsWeight=wjetsSFmap[chName+"_bveto"];}

		////////////////////////////////////////////////////////////////////////////
		// SELECTION REGIONS
		//
		// presel_        : One lepton, zero veto leptons
		// eq1jet_        : presel + exactly one jet
		// geq3jet1btag_  : presel + three or more jets, one or more b-tags
		// geq4jet_       : presel + four or more jets
		// btag_          : presel + four or more jets, two or more b-tags
		// bveto_         : presel + four or more jets, no b-tags
		//
		////////////////////////////////////////////////////////////////////////////

		// fill different control distributions
		float charge = 0;
		if(passSoftLeptonVeto) {
			controlHistos.fillHisto("njets", ch, selJets.size(), weight);
		}
		// charge: pdg id conventions, e.g. +11 = electron, -11 = positron etc.
		// FIXME: there seems to be a problem concerning the charge... have more W- than W+ which gives negative scale factors...
		charge = (-1.)*selLeptons[0].get("id")/abs(selLeptons[0].get("id"));

		// FIXME: WJets sample seem to have different convention than rest (e.g. single top), temporary fix:
		if((url.Contains("WJets") || url.Contains("W1Jets") || url.Contains("W2Jets") || url.Contains("W3Jets") || url.Contains("W4Jets")) && !(url.Contains("TTWJets"))){
			charge *= -1.;
		}
		controlHistos.addHistogram( new TH1F("chargeCheck",";recoCharge * mcCharge;Events",3,-1.5,1.5) );
		if(isMC){
			const data::PhysicsObject_t &genParton=selLeptons[0].getObject("gen");
			int genPartonId=genParton.info.find("id")->second;
			bool lepton_matched(abs(genPartonId) == abs(selLeptons[0].get("id")));
			if (lepton_matched){
				int mc_charge = (-1.)*genPartonId/abs(genPartonId);
				controlHistos.fillHisto("chargeCheck", ch, charge*mc_charge, weight);
			}
		}

		for(size_t icat=0; icat<ctrlCategs.size(); icat++){
			float iweight(weight);
			iweight *=wjetsWeight;

			for(size_t ilep=0; ilep<1; ilep++){
				if(!passSoftLeptonVeto) continue;
				// fill muon specific variables
				if(abs(selLeptons[ilep].get("id"))==13) {
					Float_t gIso    = selLeptons[ilep].getVal("gIso04");
					Float_t chIso   = selLeptons[ilep].getVal("chIso04");
					Float_t puchIso = selLeptons[ilep].getVal("puchIso04");
					Float_t nhIso   = selLeptons[ilep].getVal("nhIso04");
					Float_t relIso=(TMath::Max(nhIso+gIso-0.5*puchIso,0.)+chIso)/selLeptons[ilep].pt();
					controlHistos.fillHisto(ctrlCategs[icat]+"iso" , ch, relIso, iweight);
				}

				// fill electron specific variables
				if(abs(selLeptons[ilep].get("id"))==11) {
					controlHistos.fillHisto(ctrlCategs[icat]+"emva", ch, selLeptons[ilep].getVal("mvatrig"), iweight);

					float sceta=selLeptons[ilep].getVal("sceta");
					Float_t gIso    = selLeptons[ilep].getVal("gIso03");
					Float_t chIso   = selLeptons[ilep].getVal("chIso03");
					// Float_t puchIso = selLeptons[ilep].getVal("puchIso03");
					Float_t nhIso   = selLeptons[ilep].getVal("nhIso03");
					float relIso=(TMath::Max(nhIso+gIso-ev.rho*utils::cmssw::getEffectiveArea(11,sceta),Float_t(0.))+chIso)/selLeptons[ilep].pt();
					controlHistos.fillHisto(ctrlCategs[icat]+"iso" , ch, relIso, iweight);
				}
			}

			// fill control histograms
			// add corresponding control distributions
			controlHistos.fillHisto(ctrlCategs[icat]+"met",             ch, met[0].pt(),    iweight);
			controlHistos.fillHisto(ctrlCategs[icat]+"charge",          ch, charge,         iweight);
			if(ctrlCategs[icat]!="") controlHistos.fillHisto(ctrlCategs[icat]+"njets",  ch, selJets.size(), iweight);



			for(size_t ijet=0; ijet<selJets.size(); ijet++){
				TString label("jet");
				label+=(ijet+1);
				const data::PhysicsObject_t &genJet=selJets[ijet].getObject("genJet");
				int flavId=genJet.info.find("id")->second;
				if(abs(flavId)==5 || abs(flavId)==4 ) flavId=abs(flavId)-1;
				else if(abs(flavId)>6)                flavId=1;
				else if(abs(flavId)==0)               flavId=0;
				else                                  flavId=2;
				controlHistos.fillHisto(ctrlCategs[icat]+"btag_"+label+"_pt",        ch, selJets[ijet].pt(), iweight, true);
				controlHistos.fillHisto(ctrlCategs[icat]+"btag_"+label+"_eta",       ch, fabs(selJets[ijet].eta()), iweight);
				controlHistos.fillHisto(ctrlCategs[icat]+"btag_"+label+"_flav",      ch, abs(flavId), iweight);
				controlHistos.fillHisto(ctrlCategs[icat]+"btag_"+label+"_nobsmearpt",ch, abs(flavId)==5 ? selJets[ijet].pt() : selJets[ijet].getVal("jer"), iweight);
				controlHistos.fillHisto(ctrlCategs[icat]+"btag_"+label+"_smearpt",   ch,                                       selJets[ijet].getVal("jer"), iweight);
			}
		}

		if(passSoftLeptonVeto && pass3JetSelection  && pass1BtagSelection) lxyAn.analyze(selLeptons,selJets,met[0],ev.nvtx,gen,weight*wjetsWeight);


		//select the event
		if(!passSoftLeptonVeto) continue;
		controlHistos.fillHisto("evtflow", ch, 1, weight*wjetsWeight);
		controlHistos.fillHisto("lxyevtflow", ch, 1, weight*wjetsWeight);
		controlHistos.fillHisto("ueevtflow", ch, 1, weight*wjetsWeight);

		// lxy cutflow histogram
		if(pass3JetSelection) {
			controlHistos.fillHisto("lxyevtflow", ch, 2, weight*wjetsWeight);
			if(pass1BtagSelection)
				controlHistos.fillHisto("lxyevtflow", ch, 3, weight*wjetsWeight);
		}

		if(passJetSelection) {
			controlHistos.fillHisto("evtflow", ch, 2, weight*wjetsWeight);

			//UE event analysis with PF candidates
			controlHistos.fillHisto("ueevtflow", ch, 2, weight*wjetsWeight);
			if(looseJets[0].pt()>30 && looseJets[1].pt()>30 && fabs(looseJets[0].eta())<2.5 && fabs(looseJets[1].eta())<2.5) {
				if(looseJets[0].getVal("csv")>0.405 && looseJets[1].getVal("csv")>0.405) {
					// float weight( weight*ibtagdyWeight );
					float weight( weight*ibtagdyWeight ); // FIXME: check weight
					controlHistos.fillHisto("ueevtflow", ch, 3, weight*wjetsWeight);
					data::PhysicsObjectCollection_t pf = evSummary.getPhysicsObject(DataEventSummaryHandler::PFCANDIDATES);
					// ueAn.analyze(selLeptons,looseJets,met[0],pf,gen,ev.nvtx,weight*wjetsWeight); // FIXME: UE analysis need 2 leptons?
					// if(saveSummaryTree) ueAn.fillSummaryTuple(xsecWeight);
				}
			}
		}

		//
		// STATISTICAL ANALYSIS (with systs variations)
		//
		for(size_t ivar=0;ivar<systVars.size(); ivar++){
			TString var=systVars[ivar];

			//re-select the jets
			int njets(0),nlocalbtags(0);
			for(size_t ijet=0; ijet<jets.size(); ijet++){
				float pt(jets[ijet].pt());
				if(var=="jerup")     pt=jets[ijet].getVal("jerup");
				if(var=="jerdown")   pt=jets[ijet].getVal("jerdown");
				if(var=="jesup")     pt=jets[ijet].getVal("jesup");
				if(var=="jesdown")   pt=jets[ijet].getVal("jesdown");
				if(fabs(jets[ijet].eta())<2.5 && pt>30) {
					njets++;
					nlocalbtags += (jets[ijet].getVal("csv")>0.405);
				}
			}
			if(njets<1) continue;
			bool passLocalJetSelection(njets>1);

			// FIXME: check and apply to QCD estimate
			//re-select the MET
			int metIdx(0);
			if(var=="jerup")    metIdx=utils::cmssw::JERUP;   if(var=="jerdown")  metIdx=utils::cmssw::JERDOWN;
			if(var=="jesup")    metIdx=utils::cmssw::JESUP;   if(var=="jesdown")  metIdx=utils::cmssw::JESDOWN;
			if(var=="umetup")   metIdx=utils::cmssw::UMETUP;  if(var=="umetdown") metIdx=utils::cmssw::UMETDOWN;
			LorentzVector iMet=met[metIdx];
			bool passLocalMetSelection(iMet.pt()>40 );


			// FIXME: define corresponding categories... QCD estimate
			//event category and dy scale factor
			TString localCtrlCateg(""), btagCtrlCateg("");
			float idyWeight(1.0),ibtagdyWeight(1.0);

			// FIXME: !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

			// if(passSoftLeptonVeto && passJetSelection  )                      { ctrlCategs.push_back("");        }
			// if(passSoftLeptonVeto && selJets.size()==1 )                      { ctrlCategs.push_back("eq1jets"); }
			// if(passSoftLeptonVeto && passJetSelection  && passBtagSelection)  { ctrlCategs.push_back("btag");    }
			// if(passSoftLeptonVeto && passJetSelection  && passBtagVeto)       { ctrlCategs.push_back("bveto");   }

			// //event category and dy scale factor
			// TString localCtrlCateg(""), btagCtrlCateg("");
			// float idyWeight(1.0),ibtagdyWeight(1.0);
			// if(!passLocalMetSelection && !passLocalJetSelection) localCtrlCateg="eq1jetslowmet";
			// if(!passLocalMetSelection && passLocalJetSelection)  localCtrlCateg="lowmet";
			// if(passLocalMetSelection  && !passLocalJetSelection) { localCtrlCateg="eq1jets"; if(dySFmap.find(chName+"eq1jets")!=dySFmap.end()) idyWeight=dySFmap[chName+"eq1jets"];    }
			// if(passLocalMetSelection  && passLocalJetSelection)  { localCtrlCateg="";        if(dySFmap.find(chName)!=dySFmap.end())           idyWeight=dySFmap[chName];              }
			// if(passLocalJetSelection  && nlocalbtags>=2)         { btagCtrlCateg="osbtag";   if(dySFmap.find(chName+"osbtag")!=dySFmap.end())  ibtagdyWeight=dySFmap[chName+"osbtag"]; }
			// if(passLocalJetSelection  && nlocalbtags==0)         { btagCtrlCateg="osbveto"; }

			//re-assign event weight
			float iweight(weightNom);
			if(var=="puup")   iweight=weightUp;
			if(var=="pudown") iweight=weightDown;
			iweight *= lScaleFactor;
			if(passLocalMetSelection) iweight*=idyWeight;

			// float thetall=utils::cmssw::getArcCos<LorentzVector>(selLeptons[0],selLeptons[1]);
			// float mtsum=utils::cmssw::getMT<LorentzVector>(selLeptons[0],iMet)+utils::cmssw::getMT<LorentzVector>(selLeptons[1],iMet);
			// controlHistos.fillHisto(localCtrlCateg+"mtsum"+systVars[ivar],        ch, mtsum,          iweight);
			// controlHistos.fillHisto(localCtrlCateg+"dilarccosine"+systVars[ivar], ch, thetall,        iweight);
			// if(btagCtrlCateg!="")
			//   {
			//     float ibtagweight(iweight*ibtagdyWeight/idyWeight);
			//     controlHistos.fillHisto(btagCtrlCateg+"mtsum"+systVars[ivar],        ch, mtsum,          ibtagweight);
			//     controlHistos.fillHisto(btagCtrlCateg+"dilarccosine"+systVars[ivar], ch, thetall,        ibtagweight);
			//   }

			//final selection
			if(!passLocalMetSelection) continue;
			if(njets<1 || njets>4) continue;
			controlHistos.fillHisto("finalevtflow"+systVars[ivar], ch, njets-1, iweight);
		}

	}
	if(nDuplicates) cout << "[Warning] found " << nDuplicates << " duplicate events in this ntuple" << endl;




	//
	// close opened files
	//
	inF->Close();
	if(spyFile){
		spyDir->cd();
		spyEvents->getTree()->Write();
		// ueNtuple->Write();
		spyFile->Close();
	}

	//
	// finally, save histos to local file
	//
	TString outUrl(out);
	gSystem->ExpandPathName(outUrl);
	gSystem->Exec("mkdir -p " + outUrl);
	outUrl += "/";
	outUrl += proctag;
	if(mcTruthMode!=0) { outUrl += "_filt"; outUrl += mcTruthMode; }
	outUrl += ".root";
	TFile *file=TFile::Open(outUrl, "recreate");
	controlHistos.Write();
	lxyAn.writeTree((*file));
	file->Close();

	//that's all folks!
}
