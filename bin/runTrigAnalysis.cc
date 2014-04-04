#include <iostream>
#include <boost/shared_ptr.hpp>
#include <fstream>

#include "UserCode/llvv_fwk/interface/MacroUtils.h"
#include "UserCode/llvv_fwk/interface/SmartSelectionMonitor.h"
#include "UserCode/llvv_fwk/interface/DataEventSummaryHandler.h"
#include "UserCode/llvv_fwk/interface/LeptonEfficiencySF.h"
#include "UserCode/llvv_fwk/interface/MuScleFitCorrector.h"

#include "UserCode/llvv_fwk/interface/TrigAnalysis.h"

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
	// Load framework libraries
	gSystem->Load( "libFWCoreFWLite" );
	AutoLibraryLoader::enable();

	// Check arguments
	if ( argc < 2 ) {
		std::cout << "Usage : " << argv[0] << " parameters_cfg.py" << std::endl;
		return 0;
	}

	//
	// Configure
	//
	const edm::ParameterSet &runProcess = edm::readPSetsFrom(argv[1])->getParameter<edm::ParameterSet>("runProcess");
	std::vector<std::string> urls=runProcess.getParameter<std::vector<std::string> >("input");
	TString url = TString(urls[0]);
	TString baseDir     = runProcess.getParameter<std::string>("dirName");
	TString jecDir      = runProcess.getParameter<std::string>("jecDir");
	bool isMC           = runProcess.getParameter<bool>("isMC");
	int mcTruthMode     = runProcess.getParameter<int>("mctruthmode");
	double xsec         = runProcess.getParameter<double>("xsec");
	TString out         = runProcess.getParameter<std::string>("outdir");
	bool isV0JetsMC(isMC && (url.Contains("DYJetsToLL_50toInf") || url.Contains("WJets")));
	bool saveSummaryTree = runProcess.getParameter<bool>("saveSummaryTree");
	std::vector<string>  weightsFile = runProcess.getParameter<std::vector<string> >("weightsFile");
	int maxEvents     = runProcess.getParameter<int>("maxEvents");

	// Jet energy scale uncertainties
	gSystem->ExpandPathName(jecDir);
	FactorizedJetCorrector *jesCor        = utils::cmssw::getJetCorrector(jecDir,isMC);
	JetCorrectionUncertainty *totalJESUnc = new JetCorrectionUncertainty((jecDir+"/MC_Uncertainty_AK5PFchs.txt").Data());

	//
	// check input file
	//
	TChain *chain = new TChain(baseDir+"/data");
	for( size_t i = 0; i < urls.size(); ++i ){
		cout << "Adding file " << TString(urls[i]) << endl;
		chain->Add(TString(urls[i]));
	}

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

	//
	// control histograms
	//
	SmartSelectionMonitor controlHistos;
	TH1F* Hcutflow = (TH1F*) controlHistos.addHistogram(new TH1F ("cutflow", "cutflow", 4, 0, 4)) ;

	//vertex multiplicity
	// controlHistos.addHistogram( new TH1F ("nvertices", "; Vertex multiplicity; Events", 50, 0.,50.) );

	//lepton efficiencies
	LeptonEfficiencySF lepEff;

	//
	// process events file
	//
	DataEventSummaryHandler evSummary;
	if( !evSummary.attach( chain ) )  {
		cout << "Mooep" << endl;
		inF->Close();
		return -1;
	}
	// if( !evSummary.attach( (TTree *) inF->Get(baseDir+"/data") ) )  { inF->Close();  return -1; }
	Int_t entries_to_process = -1;
	if(maxEvents > 0) entries_to_process = maxEvents;
	else              entries_to_process = evSummary.getEntries();
	const Int_t totalEntries = entries_to_process;


	TrigAnalysis trigAn(runProcess, &evSummary);
	trigAn.initializeTreeVariables();

	float cnorm=1.0;

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

		Hcutflow->Fill(1,1);
		Hcutflow->Fill(2,weightNom);
		Hcutflow->Fill(3,weightUp);
		Hcutflow->Fill(4,weightDown);


		//trigger bits
		bool eTrigger   = ev.t_bits[13] || ev.t_bits[14];
		bool muTrigger  = ev.t_bits[6];
		if(filterOnlyE)   { muTrigger=false; }
		if(filterOnlyMU)  { eTrigger=false;  }


		// Trigger selection
		if(!eTrigger && !muTrigger) continue;

		// // Determine the lepton+jets channel
		// data::PhysicsObject_t firstLepton;
		// if (selLeptons.size() > 0) firstLepton = selLeptons[0];

		// ev.cat = 1;
		// float lScaleFactor(1.0);
		// ev.cat *= firstLepton.get("id");

		// int id(abs(firstLepton.get("id")));
		// lScaleFactor = isMC ? lepEff.getLeptonEfficiency( firstLepton.pt(), firstLepton.eta(), id,  id ==11 ? "loose" : "tight" ).first : 1.0; // FIXME: check!

		// TString chName;
		// if     (abs(ev.cat)==11 && eTrigger)   { chName="e";  }
		// else if(abs(ev.cat)==13 && muTrigger)  { chName="mu"; }
		// else                                   continue; // don't want e events triggered by mu, and vice versa
		// std::vector<TString> ch(1,chName);

		// Apply data/mc correction factors and update the event weight
		// float weight = weightNom*lScaleFactor;

		float weight = weightNom;

		// // Jet/MET
		// data::PhysicsObjectCollection_t recoMet = evSummary.getPhysicsObject(DataEventSummaryHandler::MET);
		// data::PhysicsObjectCollection_t jets    = evSummary.getPhysicsObject(DataEventSummaryHandler::JETS);
		// utils::cmssw::updateJEC(jets,jesCor,totalJESUnc,ev.rho,ev.nvtx,isMC);
		// std::vector<LorentzVector> missingEt = utils::cmssw::getMETvariations(recoMet[0],jets,selLeptons,isMC);

		// data::PhysicsObjectCollection_t looseJets, selJets;
		// int nbtags = applyJetSelection(jets, selLeptons, vetoLeptons, looseJets, selJets);

		trigAn.analyze(weight);

	}
	cout << endl;
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
	trigAn.writeTree((*file));
	file->Close();

	//that's all folks!
}
