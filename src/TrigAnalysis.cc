#include "UserCode/llvv_fwk/interface/TrigAnalysis.h"
#include <Math/VectorUtil.h>
#include "TSystem.h"

#include "UserCode/llvv_fwk/interface/MuScleFitCorrector.h"

using namespace std;
typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >::BetaVector BetaVector;

//
TrigAnalysis::TrigAnalysis(edm::ParameterSet runProcess, DataEventSummaryHandler* evSummary, SmartSelectionMonitor &mon) : mon_(&mon){
	fEvSummary_ = evSummary;

	fJECDir = runProcess.getParameter<std::string>("jecDir");
	gSystem->ExpandPathName(fJECDir);
	std::vector<std::string> urls = runProcess.getParameter<std::vector<std::string> >("input");
	fURL = TString(urls[0]);
	fMuCor = getMuonCorrector(fJECDir,fURL);
	fIsMC = runProcess.getParameter<bool>("isMC");

	// Start monitoring histograms of this analysis
	TH1F* H_recoChannel = (TH1F*) mon_->addHistogram(
		new TH1F ("recoChannel", "recoChannel;Channel (reconstructed);Events", 6, 0, 6)) ;
	TH1F* H_genChannel = (TH1F*) mon_->addHistogram(
		new TH1F ("genChannel", "genChannel;Channel (generated);Events", 6, 0, 6)) ;
	TH2F* H_channel = (TH2F*) mon_->addHistogram(
		new TH2F ("channel", "channel;Channel (generated);Channel (reconstructed);Events",
			      6, 0, 6, 6, 0, 6)) ;

}


bool TrigAnalysis::isGoodElectron(data::PhysicsObject_t ele, float minpt, float maxeta){
	// Kinematic cuts
	if(ele.pt()<minpt) return false;
	if(fabs(ele.eta())>maxeta) return false;
	float sceta=ele.getVal("sceta");
	if(fabs(sceta)>1.4442 && fabs(sceta)<1.5660) return false;

	// Isolation
	Float_t gIso    = ele.getVal("gIso03");
	Float_t chIso   = ele.getVal("chIso03");
	Float_t nhIso   = ele.getVal("nhIso03");
	float relIso=(TMath::Max(nhIso+gIso-fCurrentEvent_.rho*utils::cmssw::getEffectiveArea(11,sceta),Float_t(0.))+chIso)/ele.pt();
	if(relIso>0.1) return false;

	// ID
	bool passId = true;
	if(ele.getFlag("isconv"))            passId =false;
	if(ele.getVal("tk_d0")>0.2)          passId =false; // FIXME: is this the correct value?
	if(ele.getVal("tk_lostInnerHits")>0) passId =false;
	if(ele.getVal("mvatrig")<0.5)        passId =false;
	if(!passId) return false;

	return true;
}

bool TrigAnalysis::isGoodMuon(data::PhysicsObject_t mu, float minpt, float maxeta){
	// Muon energy scale and uncertainties
	Int_t id=mu.get("id");
	if(fMuCor){
		TLorentzVector p4(mu.px(),mu.py(),mu.pz(),mu.energy());
		fMuCor->applyPtCorrection(p4 , id<0 ? -1 : 1 );
		if(fIsMC) fMuCor->applyPtSmearing(p4, id<0 ? -1 : 1, false);
		mu.SetPxPyPzE(p4.Px(),p4.Py(),p4.Pz(),p4.E());
	}

	// Kinematic cuts
	if(mu.pt()<minpt) return false;
	if(fabs(mu.eta())>maxeta) return false;

	// Isolation
	Float_t gIso    = mu.getVal("gIso04");
	Float_t chIso   = mu.getVal("chIso04");
	Float_t puchIso = mu.getVal("puchIso04");
	Float_t nhIso   = mu.getVal("nhIso04");
	Float_t relIso=(TMath::Max(nhIso+gIso-0.5*puchIso,0.)+chIso)/mu.pt();
	if(relIso>0.12) return false;

	// ID
	Int_t idbits    = mu.get("idbits");
	bool isTight    = ((idbits >> 10) & 0x1);
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
		if( isGoodJet(jets[ijet], leptons, 30., 2.5) )
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
		if(  (abs(id)==11 && isGoodElectron(leptons[ilep], 30., 2.5)) ||
		     (abs(id)==13 && isGoodMuon(leptons[ilep], 26., 2.1)) ) {
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
	if( channel == "had") return 0;
	if( channel == "e")   return 1;
	if( channel == "m")   return 2;
	if( channel == "ee")  return 3;
	if( channel == "mm")  return 4;
	if( channel == "em")  return 5;
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

	// std::cout << leptons.size() << " " << jets.size() << " " << nbtags;

	// Check channel
	TString ch = getRecoChannel(leptons);
	TString gench = getGenChannel();

	// std::cout << " " << ch << " " << gench << std::endl;

	if( gench == "?" ){
		std::cout << " Moooooep " << std::endl;
	}

	mon_->fillHisto("recoChannel", "all", channelToBin(ch),    weight);
	mon_->fillHisto("genChannel",  "all", channelToBin(gench), weight);
	mon_->fillHisto("channel",     "all", channelToBin(gench), channelToBin(ch), weight);

	// //select
	// bool passEventSelection(true);
	// if(leptons.size()>1) {
	// 	float mll=(leptons[0]+leptons[1]).mass();
	// 	bool isZcand((ch=="ee" || ch=="mumu") && fabs(mll-91)<15);
	// 	bool passMet( ch=="emu" || ((ch=="ee" || ch=="mumu") && met.pt()>40));
	// 	bool isOS( leptons[0].get("id")*leptons[1].get("id") < 0 );
	// 	bool passJet(jets.size()>2);
	// 	if(isZcand || !passMet || !isOS || !passJet) passEventSelection=false;
	// }
	// if(!passEventSelection) return;

	// //get leading Lxy jet
	// data::PhysicsObject_t *lxyJet=0;
	// float leadingLxy(-1), leadingLxyErr(-1), leadingSecVtxNtrk(-1), leadingSecVtxMass(-1);
	// for(size_t ijet=0; ijet<jets.size(); ijet++){
	// 	// Get reconstructed secondary vertex
	// 	const data::PhysicsObject_t &svx = jets[ijet].getObject("svx");
	// 	float lxy = svx.vals.find("lxy")->second;
	// 	if(lxy <= 0) continue;
	// 	if(lxy < leadingLxy) continue;

	// 	//use lxy if it's larger, prefer also more central jets
	// 	if(lxyJet==0 || (fabs(lxyJet->eta())>1.1 && fabs(jets[ijet].eta())<fabs(lxyJet->eta()))){
	// 		lxyJet=&(jets[ijet]);
	// 		leadingLxy        = lxy;
	// 		leadingLxyErr     = svx.vals.find("lxyErr")->second;
	// 		leadingSecVtxNtrk = svx.info.find("ntrk")->second;
	// 		leadingSecVtxMass = svx.M();
	// 	}
	// }
	// if(leadingLxy < 0) return; // Select events with at least one Lxy found
	// bool passFiducialCut(fabs(lxyJet->eta())<1.1);

	// // MC truth
	// // get flavor
	// const data::PhysicsObject_t &genJet    = lxyJet->getObject("genJet");
	// const data::PhysicsObject_t &genParton = lxyJet->getObject("gen");
	// int flavId = genJet.info.find("id")->second;
	// int genId  = genParton.info.find("id")->second;

	// //match to a b-hadron
	// data::PhysicsObject_t *bHad=0,*top=0,*antitop=0;
	// for(size_t imc=0; imc<mctruth.size(); imc++){
	// 	int id = mctruth[imc].get("id");
	// 	if     (id==6)  top     = &(mctruth[imc]);
	// 	else if(id==-6) antitop = &(mctruth[imc]);
	// 	else if(abs(id)<500) continue;

	// 	if(deltaR(mctruth[imc],*lxyJet)>0.5) continue;
	// 	bHad=&mctruth[imc];
	// 	break;
	// }
	// float genLxy(-1);
	// if(bHad) genLxy=bHad->getVal("lxy");

	// //resolution plots
	// if(genLxy>0 && abs(flavId)==5){
	// 	mon_->fillHisto("bjetlxyreseta",ch,leadingLxy-genLxy,fabs(lxyJet->eta()),weight);

	// 	if(passFiducialCut){
	// 		mon_->fillHisto("bjetlxyrespt", ch,leadingLxy-genLxy,lxyJet->pt(),weight);

	// 		int absid=abs(bHad->get("id"));
	// 		int bHadronBin(3),fullBHadronBin(4); // other not so relevant
	// 		if(     absid==511 || absid==10511 || absid==513 || absid==10513 || absid==20513 || absid==515) {bHadronBin=0; fullBHadronBin=0;} // B0 family
	// 		else if(absid==521 || absid==10521 || absid==523 || absid==10523 || absid==20523 || absid==525) {bHadronBin=1; fullBHadronBin=1;} // B+ family
	// 		else if(absid==531 || absid==10531 || absid==533 || absid==10533 || absid==20533 || absid==535) {bHadronBin=2; fullBHadronBin=2;} // B0s family
	// 		else if(absid==541 || absid==10541 || absid==543 || absid==10543 || absid==20543 || absid==545) {bHadronBin=3; fullBHadronBin=3;} // B+c family
	// 		mon_->fillHisto("bjetlxybhad", ch, leadingLxy,  bHadronBin, weight);
	// 		mon_->fillHisto("lxybhad",     ch, fullBHadronBin,          weight);

	// 		if(genParton.pt()>0){
	// 			float xb(bHad->pt()/genParton.pt());
	// 			float bBoost = genParton.BoostToCM().R();
	// 			float genLxyCM=sqrt(1-bBoost*bBoost)*genLxy;

	// 			TString pf(""); pf+=bHadronBin;
	// 			mon_->fillHisto("bjetlxygenb"+pf+"had", ch, leadingLxy,     genLxyCM,   weight);
	// 			mon_->fillHisto("xbvslxy",              ch, xb,             leadingLxy, weight);
	// 			mon_->fillHisto("ptbvslxy",             ch, genParton.pt(), leadingLxy, weight);
	// 		}
	// 	}
	// }

	// // Lxy distribution
	// if(!passFiducialCut) return;
	// mon_->fillHisto("jetlxy",     ch, leadingLxy,               weight); // this will be used for the mass measurement later
	// mon_->fillHisto("jetlxysig",  ch, leadingLxy/leadingLxyErr, weight);
	// mon_->fillHisto("secvtxmass", ch, leadingSecVtxMass,        weight);
	// mon_->fillHisto("secvtxntrk", ch, leadingSecVtxNtrk,        weight);

	// if(top && antitop && abs(genId) == 5){
	// 	float genpt( genId == -5 ? antitop->pt() : top->pt() );
	// 	mon_->fillHisto("topptvslxy", ch, genpt, leadingLxy, weight );
	// }

	// // set some tree variables
	// flxy        = leadingLxy;
	// flxyerr     = leadingLxyErr;
	// fchannel    = abs(lid1*lid2);
	// fsecvtxmass = leadingSecVtxMass;
	// fsecvtxntrk = leadingSecVtxNtrk;
	// fnvtx       = nvtx;
	// fweight     = weight;

	tree_->Fill();
}

void TrigAnalysis::initializeTreeVariables() {
	tree_ = new TTree("tree","tree");
	tree_->Branch("lxy"        ,&flxy          ,"lxy/F");
	tree_->Branch("lxyerr"     ,&flxyerr       ,"lxyerr/F");
	// tree_->Branch("lxysig"     ,&flxyerr       ,"lxysig/F");
	tree_->Branch("nvtx"        ,&fnvtx          ,"nvtx/I");
	tree_->Branch("weight"      ,&fweight         ,"weight/F");
	tree_->Branch("channel"     ,&fchannel        ,"channel/I");
	// tree_->Branch("charge"      ,&fcharge         ,"charge/I");
	tree_->Branch("secvtxmass" ,&fsecvtxmass   ,"secvtxmass/F");
	tree_->Branch("secvtxntrk" ,&fsecvtxntrk   ,"secvtxntrk/F");
	// tree_->Branch("pt"         ,&fpt           ,"pt/F");
	// tree_->Branch("eta"        ,&feta          ,"eta/F");
	// tree_->Branch("btag"       ,&fbtag         ,"btag/F");
	// tree_->Branch("puW"        ,&fpuW          ,"puW/F");
	// tree_->Branch("gamma"      ,&fgamma        ,"gamma/F");
	// tree_->Branch("dR"         ,&fdR           ,"dR/F");
	// tree_->Branch("recopt"     ,&frecopt       ,"recopt/F");
	// tree_->Branch("recoeta"    ,&frecoeta      ,"recoeta/F");

}


void TrigAnalysis::clearTreeVariables() {
	flxy        = 99.99;
	flxyerr     = -99.99;
	fnvtx       = -1;
	fweight     = 1.;
	fsecvtxmass = -1.;
	fsecvtxntrk = -1.;
	fchannel    = 0;
	// fcharge = -99;
}


void TrigAnalysis::writeTree(TFile &f) {
	f.cd();
	// tree_->Print();
	tree_->Write();
}
