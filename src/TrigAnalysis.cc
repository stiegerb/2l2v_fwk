#include "UserCode/llvv_fwk/interface/TrigAnalysis.h"
#include <Math/VectorUtil.h>
#include "TSystem.h"

#include "UserCode/llvv_fwk/interface/MuScleFitCorrector.h"

using namespace std;
typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >::BetaVector BetaVector;

//
TrigAnalysis::TrigAnalysis(edm::ParameterSet runProcess, DataEventSummaryHandler* evSummary){
	fEvSummary_ = evSummary;

	fJECDir = runProcess.getParameter<std::string>("jecDir");
	gSystem->ExpandPathName(fJECDir);
	std::vector<std::string> urls = runProcess.getParameter<std::vector<std::string> >("input");
	fURL = TString(urls[0]);
	fMuCor = getMuonCorrector(fJECDir,fURL);
	fIsMC = runProcess.getParameter<bool>("isMC");
}


bool TrigAnalysis::isGoodElectron(data::PhysicsObject_t ele, float minpt, float maxeta){
	// Kinematic cuts
	if( ele.pt() < minpt ) return false;
	if( fabs(ele.eta()) > maxeta ) return false;
	float sceta = ele.getVal("sceta");
	if( fabs(sceta) > 1.4442 && fabs(sceta) < 1.5660 ) return false;

	// Isolation
	Float_t gIso  = ele.getVal("gIso03");
	Float_t chIso = ele.getVal("chIso03");
	Float_t nhIso = ele.getVal("nhIso03");
	float relIso = (TMath::Max(nhIso+gIso-fCurrentEvent_.rho*utils::cmssw::getEffectiveArea(11,sceta),Float_t(0.))+chIso)/ele.pt();
	if( relIso > 0.1 ) return false;

	// ID
	bool passId = true;
	if( ele.getFlag("isconv") )              passId = false;
	if( ele.getVal("tk_d0")>0.2 )            passId = false; // FIXME: is this the correct value?
	if( ele.getVal("tk_lostInnerHits") > 0 ) passId = false;
	if( ele.getVal("mvatrig")<0.5 )          passId = false;
	if( !passId ) return false;

	return true;
}

bool TrigAnalysis::isGoodMuon(data::PhysicsObject_t mu, float minpt, float maxeta){
	// Muon energy scale and uncertainties
	Int_t id = mu.get("id");
	if( fMuCor ){
		TLorentzVector p4(mu.px(),mu.py(),mu.pz(),mu.energy());
		fMuCor->applyPtCorrection(p4 , id<0 ? -1 : 1 );
		if( fIsMC ) fMuCor->applyPtSmearing(p4, id<0 ? -1 : 1, false);
		mu.SetPxPyPzE(p4.Px(),p4.Py(),p4.Pz(),p4.E());
	}

	// Kinematic cuts
	if( mu.pt() < minpt ) return false;
	if( fabs(mu.eta()) > maxeta ) return false;

	// Isolation
	Float_t gIso    = mu.getVal("gIso04");
	Float_t chIso   = mu.getVal("chIso04");
	Float_t puchIso = mu.getVal("puchIso04");
	Float_t nhIso   = mu.getVal("nhIso04");
	Float_t relIso = ( TMath::Max(nhIso+gIso-0.5*puchIso,0.)+chIso ) / mu.pt();
	if( relIso > 0.12 ) return false;

	// ID
	Int_t idbits = mu.get("idbits");
	bool isTight = ((idbits >> 10) & 0x1);
	if(!isTight) return false;

	return true;
}

bool TrigAnalysis::isGoodJet(data::PhysicsObject_t jet,
	                         data::PhysicsObjectCollection_t leptons,
	                         float minpt, float maxeta){
	// Cross-clean with selected leptons
	double minDRlj(9999.);
	for( size_t ilep=0; ilep<leptons.size(); ilep++ )
		minDRlj = TMath::Min( minDRlj, deltaR(jet, leptons[ilep]) );
	if( minDRlj<0.4 ) return false;

	// Require to pass the loose id
	Int_t idbits = jet.get("idbits");
	bool passPFloose( ((idbits>>0) & 0x1));
	if(!passPFloose) return false;

	// Top candidate jets
	if( jet.pt() < minpt ) return false;
 	if( fabs(jet.eta()) > maxeta ) return false;

	return true;
}

int TrigAnalysis::selectJets(data::PhysicsObjectCollection_t jets,
	                         data::PhysicsObjectCollection_t leptons,
	                         data::PhysicsObjectCollection_t &selJets){
	for(size_t ijet=0; ijet<jets.size(); ijet++){
		if( isGoodJet(jets[ijet], leptons, 10., 5.0) )
			selJets.push_back(jets[ijet]);
	}
	sort(selJets.begin(), selJets.end(), data::PhysicsObject_t::sortByPt);

	int nbtags(0);
	for(size_t ijet=0; ijet<selJets.size(); ijet++){
		if( jets[ijet].getVal("csv") > 0.679 ) nbtags++; // Loose = 0.244 Med = 0.679 Tight = 0.898
	}

	return nbtags;
}

void TrigAnalysis::selectLeptons(data::PhysicsObjectCollection_t leptons,
	                             data::PhysicsObjectCollection_t &selLeptons){
	for(size_t ilep=0; ilep<leptons.size(); ilep++){
		Int_t id=leptons[ilep].get("id");
		if(  (abs(id)==11 && isGoodElectron(leptons[ilep], 10., 2.5)) ||
		     (abs(id)==13 && isGoodMuon(leptons[ilep], 10., 2.5)) ) {
			selLeptons.push_back(leptons[ilep]);
		}
	}
	sort(selLeptons.begin(), selLeptons.end(), data::PhysicsObject_t::sortByPt);
}

TString TrigAnalysis::getRecoChannel(data::PhysicsObjectCollection_t leptons){
	TString channel = "had";
	if( leptons.size() == 0 ) return channel;

	int lid1 = leptons[0].get("id");
	if     ( abs(lid1) == 11 ) channel = "e";
	else if( abs(lid1) == 13 ) channel = "m";

	if( leptons.size() == 1 ) return channel;

	int lid2 = leptons[1].get("id");
	if     ( abs(lid1)*abs(lid2) == 11*11 ) channel = "ee";
	else if( abs(lid1)*abs(lid2) == 11*13 ) channel = "em";
	else if( abs(lid1)*abs(lid2) == 13*13 ) channel = "mm";

	return channel;
}

TString TrigAnalysis::getGenChannel(){
	if( !fIsMC ) return "-";

	data::PhysicsObjectCollection_t genParticles =
	              fEvSummary_->getPhysicsObject(DataEventSummaryHandler::GENPARTICLES);
	data::PhysicsObjectCollection_t genLeptonsS3;

	for(size_t igen=0; igen<genParticles.size(); igen++){
		if(genParticles[igen].get("status")!=3) continue;

		int absid = abs(genParticles[igen].get("id"));

		if(absid!=11 && absid!=13 && absid!=15) continue;
		genLeptonsS3.push_back(genParticles[igen]);
	}

	TString channel = "had";
	if( genLeptonsS3.size() == 0 ) return channel;

	int lid1 = genLeptonsS3[0].get("id");
	if     ( abs(lid1) == 11 ) channel = "e";
	else if( abs(lid1) == 13 ) channel = "m";

	if( genLeptonsS3.size() == 1 ) return channel;

	int lid2 = genLeptonsS3[1].get("id");
	if     ( abs(lid1)*abs(lid2) == 11*11 ) channel = "ee";
	else if( abs(lid1)*abs(lid2) == 11*13 ) channel = "em";
	else if( abs(lid1)*abs(lid2) == 13*13 ) channel = "mm";

	if( genLeptonsS3.size() == 2 ) return channel;
	std::cout << " Moooooep " << std::endl;
	channel = "?";
	return channel;
}

int TrigAnalysis::channelToBin(TString channel){
	if( channel == "had") return 1;
	if( channel == "e")   return 2;
	if( channel == "m")   return 3;
	if( channel == "ee")  return 4;
	if( channel == "mm")  return 5;
	if( channel == "em")  return 6;
	return -1;
}

//
void TrigAnalysis::analyze(float weight){
	fCurrentEvent_ = fEvSummary_->getEvent();
	clearTreeVariables();

	data::PhysicsObjectCollection_t leptons, jets;
	selectLeptons(fEvSummary_->getPhysicsObject(DataEventSummaryHandler::LEPTONS), leptons);
	int nbtags(0);
	nbtags = selectJets(fEvSummary_->getPhysicsObject(DataEventSummaryHandler::JETS), leptons, jets);

	// Check channel
	TString ch = getRecoChannel(leptons);
	// TString gench = getGenChannel();

	// Set some tree variables
	fchannel    = channelToBin(ch);
	fnvtx       = fCurrentEvent_.nvtx;
	fweight     = weight;

	nleps = leptons.size();
	for (size_t i = 0; i < leptons.size(); ++i){
		if(i > 3) break; // no more than 4 leptons
		leppt[i] = leptons[i].pt();
		lepeta[i] = leptons[i].eta();
		lepphi[i] = leptons[i].phi();
		lepid[i] = leptons[i].get("id");
	}

	njets = jets.size();
	for (size_t i = 0; i < jets.size(); ++i){
		if(i > 9) break; // no more than 10 jets
		jetpt[i] = jets[i].pt();
		jeteta[i] = jets[i].eta();
		jetphi[i] = jets[i].phi();
		jetcsv[i] = jets[i].getVal("csv");
	}

	tree_->Fill();
}

void TrigAnalysis::initializeTreeVariables() {
	tree_ = new TTree("tree","tree");
	tree_->Branch("nVtx", &fnvtx, "nVtx/I");
	tree_->Branch("PUWeight", &fweight, "PUWeight/F");
	tree_->Branch("RChannel", &fchannel, "RChannel/I");

	tree_->Branch("RNLeps", &nleps, "RNLeps/I");
	tree_->Branch("RLepId", lepid, "RLepId[RNLeps]/I");
	tree_->Branch("RLepPt", leppt, "RLepPt[RNLeps]/F");
	tree_->Branch("RLepEta", lepeta, "RLepEta[RNLeps]/F");
	tree_->Branch("RLepPhi", lepphi, "RLepPhi[RNLeps]/F");

	tree_->Branch("RNJets", &njets, "RNJets/I");
	tree_->Branch("RJetPt", jetpt, "RJetPt[RNJets]/F");
	tree_->Branch("RJetEta", jeteta, "RJetEta[RNJets]/F");
	tree_->Branch("RJetPhi", jetphi, "RJetPhi[RNJets]/F");
	tree_->Branch("RJetCSV", jetcsv, "RJetCSV[RNJets]/F");
}


void TrigAnalysis::clearTreeVariables() {
	fnvtx       = -1;
	fweight     = -1.;
	fchannel    = -1;

	nleps = 0;
	for( int i = 0; i < 4; ++i ){
		leppt[i] = -999.;
		lepid[i] = -999;
		lepeta[i] = -999.;
		lepphi[i] = -999.;
	}

	njets = 0;
	for( int i = 0; i < 10; ++i ){
		jetpt[i] = -999.;
		jeteta[i] = -999.;
		jetphi[i] = -999.;
		jetcsv[i] = -999.;
	}
}


void TrigAnalysis::writeTree(TFile &f) {
	f.cd();
	// tree_->Print();
	tree_->Write();
}
