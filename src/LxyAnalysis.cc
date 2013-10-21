#include "UserCode/llvv_fwk/interface/LxyAnalysis.h"
#include <Math/VectorUtil.h>

using namespace std;
typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >::BetaVector BetaVector;

//
LxyAnalysis::LxyAnalysis(SmartSelectionMonitor &mon,bool runSystematics) : mon_(&mon){
	// Start monitoring histograms of this analysis
	TH1 *h = mon_->addHistogram( new TH2F ("bjetlxybhad",    "; L_{xy}(reco) [cm]; B-hadron; Jets", 50,0,5, 4,0,4) );
	h->GetYaxis()->SetBinLabel(1,"B^{0}");  h->GetXaxis()->SetBinLabel(2,"B^{+}");  h->GetXaxis()->SetBinLabel(3,"B^{0}_{s}");  h->GetXaxis()->SetBinLabel(5,"Others or B^{+}_{c}");
	h = mon_->addHistogram( new TH1F ("lxybhad", "; B-hadron; Jets", 5,0,5) );
	h->GetXaxis()->SetBinLabel(1,"B^{0}");  h->GetXaxis()->SetBinLabel(2,"B^{+}");  h->GetXaxis()->SetBinLabel(3,"B^{0}_{s}");  h->GetXaxis()->SetBinLabel(4,"B^{+}_{c}");  h->GetXaxis()->SetBinLabel(5,"Others");
	for(size_t i=0; i<4; i++){
		TString pf(""); pf+=i;
		mon_->addHistogram( new TH2F ("bjetlxygenb"+pf+"had", "; L_{xy}(reco) [cm]; B-hadron; Jets", 50,0,5,50,0,0.15) );
	}
	mon_->addHistogram( new TH2F ("bjetlxyrespt",   "; L_{xy}(reco)-L_{xy}(B) [cm]; Jet p_{T} [GeV]; Events",  50, -2.5,  2.5, 10, 30., 280.) );
	mon_->addHistogram( new TH2F ("bjetlxyreseta",  "; L_{xy}(reco)-L_{xy}(B) [cm]; Jet #eta; Events",         50, -2.5,  2.5,  4,  0., 2.5 ) );
	mon_->addHistogram( new TH2F ("xbvslxy",        "; p_{T}(B)/p_{T}(b); L_{xy} [cm]; Events",               100,    0,    2, 50,   0,   5 ) );
	mon_->addHistogram( new TH2F ("ptbvslxy",       "; b quark p_{T} [GeV]; L_{xy} [cm]; Events",             100,    0,  200, 50,   0,   5 ) );
	mon_->addHistogram( new TH2F ("topptvslxy",     "; Top p_{T} [GeV]; L_{xy} [cm]; Events",                 100,    0, 1000, 50,   0,   5 ) );
	mon_->addHistogram( new TH1F ("jetlxy",         ";L_{xy} [cm]; Jets",                                      50,   0.,   5    ) );
	mon_->addHistogram( new TH1F ("jetlxyntk",      ";Number of tracks in sec. vertex; Jets",                  10,   0.,   10.  ) );
	mon_->addHistogram( new TH1F ("jetlxysig",      ";Lxy significnace; Jets",                                 40,   0.,   20.  ) );
	mon_->addHistogram( new TH1F ("secvtxmass",     ";sec. vtx. mass [GeV]; Jets",                             40,   0.,   20.  ) );
	mon_->addHistogram( new TH1F ("secvtxntrk",     ";sec. vtx. N tracks; Jets",                               20, -0.5,   19.5 ) );
}

//
void LxyAnalysis::analyze(data::PhysicsObjectCollection_t & leptons, data::PhysicsObjectCollection_t &jets,
	// data::PhysicsObject_t &met,
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > &met, int nvtx, data::PhysicsObjectCollection_t &mctruth, float weight){

	clearTreeVariables();


	//check channel
	TString ch("");
	int lid1(leptons[0].get("id"));
	int lid2=1;
	if     (abs(lid1)==11)  ch="e";
	else if(abs(lid1)==13)  ch="mu";
	if(leptons.size()>1) {
		lid2 = (leptons[1].get("id"));
		if     (abs(lid1)*abs(lid2)==11*11)  ch="ee";
		else if(abs(lid1)*abs(lid2)==11*13)  ch="emu";
		else if(abs(lid1)*abs(lid2)==13*13)  ch="mumu";
	}
	if(ch=="") return;

	// event selection - could move final event selection here...

	//select
	bool passEventSelection(true);
	if(leptons.size()>1) {
		float mll=(leptons[0]+leptons[1]).mass();
		bool isZcand((ch=="ee" || ch=="mumu") && fabs(mll-91)<15);
		bool passMet( ch=="emu" || ((ch=="ee" || ch=="mumu") && met.pt()>40));
		bool isOS( leptons[0].get("id")*leptons[1].get("id") < 0 );
		bool passJet(jets.size()>2);
		if(isZcand || !passMet || !isOS || !passJet) passEventSelection=false;
	}
	if(!passEventSelection) return;

	//get leading Lxy jet
	data::PhysicsObject_t *lxyJet=0;
	float leadingLxy(-1), leadingLxyErr(-1), leadingSecVtxNtrk(-1), leadingSecVtxMass(-1);
	for(size_t ijet=0; ijet<jets.size(); ijet++){
		// Get reconstructed secondary vertex
		const data::PhysicsObject_t &svx = jets[ijet].getObject("svx");
		float lxy = svx.vals.find("lxy")->second;
		if(lxy <= 0) continue;
		if(lxy < leadingLxy) continue;

		//use lxy if it's larger, prefer also more central jets
		if(lxyJet==0 || (fabs(lxyJet->eta())>1.1 && fabs(jets[ijet].eta())<fabs(lxyJet->eta()))){
			lxyJet=&(jets[ijet]);
			leadingLxy        = lxy;
			leadingLxyErr     = svx.vals.find("lxyErr")->second;
			leadingSecVtxNtrk = svx.info.find("ntrk")->second;
			leadingSecVtxMass = svx.M();
		}
	}
	if(leadingLxy<0) return; // Select events with at least one Lxy found
	bool passFiducialCut(fabs(lxyJet->eta())<1.1);

	// MC truth
	// get flavor
	const data::PhysicsObject_t &genJet    = lxyJet->getObject("genJet");
	const data::PhysicsObject_t &genParton = lxyJet->getObject("gen");
	int flavId = genJet.info.find("id")->second;
	int genId  = genParton.info.find("id")->second;

	//match to a b-hadron
	data::PhysicsObject_t *bHad=0,*top=0,*antitop=0;
	for(size_t imc=0; imc<mctruth.size(); imc++){
		int id = mctruth[imc].get("id");
		if     (id==6)  top     = &(mctruth[imc]);
		else if(id==-6) antitop = &(mctruth[imc]);
		else if(abs(id)<500) continue;

		if(deltaR(mctruth[imc],*lxyJet)>0.5) continue;
		bHad=&mctruth[imc];
		break;
	}
	float genLxy(-1);
	if(bHad) genLxy=bHad->getVal("lxy");

	//resolution plots
	if(genLxy>0 && abs(flavId)==5){
		mon_->fillHisto("bjetlxyreseta",ch,leadingLxy-genLxy,fabs(lxyJet->eta()),weight);

		if(passFiducialCut){
			mon_->fillHisto("bjetlxyrespt", ch,leadingLxy-genLxy,lxyJet->pt(),weight);

			int absid=abs(bHad->get("id"));
			int bHadronBin(3),fullBHadronBin(4); // other not so relevant
			if(     absid==511 || absid==10511 || absid==513 || absid==10513 || absid==20513 || absid==515) {bHadronBin=0; fullBHadronBin=0;} // B0 family
			else if(absid==521 || absid==10521 || absid==523 || absid==10523 || absid==20523 || absid==525) {bHadronBin=1; fullBHadronBin=1;} // B+ family
			else if(absid==531 || absid==10531 || absid==533 || absid==10533 || absid==20533 || absid==535) {bHadronBin=2; fullBHadronBin=2;} // B0s family
			else if(absid==541 || absid==10541 || absid==543 || absid==10543 || absid==20543 || absid==545) {bHadronBin=3; fullBHadronBin=3;} // B+c family
			mon_->fillHisto("bjetlxybhad", ch, leadingLxy,  bHadronBin, weight);
			mon_->fillHisto("lxybhad",     ch, fullBHadronBin,          weight);

			if(genParton.pt()>0){
				float xb(bHad->pt()/genParton.pt());
				float bBoost = genParton.BoostToCM().R();
				float genLxyCM=sqrt(1-bBoost*bBoost)*genLxy;

				TString pf(""); pf+=bHadronBin;
				mon_->fillHisto("bjetlxygenb"+pf+"had", ch, leadingLxy,     genLxyCM,   weight);
				mon_->fillHisto("xbvslxy",              ch, xb,             leadingLxy, weight);
				mon_->fillHisto("ptbvslxy",             ch, genParton.pt(), leadingLxy, weight);
			}
		}
	}

	// Lxy distribution
	if(!passFiducialCut) return;
	mon_->fillHisto("jetlxy",     ch, leadingLxy,               weight); // this will be used for the mass measurement later
	mon_->fillHisto("jetlxysig",  ch, leadingLxy/leadingLxyErr, weight);
	mon_->fillHisto("secvtxmass", ch, leadingSecVtxMass,        weight);
	mon_->fillHisto("secvtxntrk", ch, leadingSecVtxNtrk,        weight);

	if(top && antitop && abs(genId) == 5){
		float genpt( genId == -5 ? antitop->pt() : top->pt() );
		mon_->fillHisto("topptvslxy", ch, genpt, leadingLxy, weight );
	}

	// set some tree variables
	flxy        = leadingLxy;
	flxyerr     = leadingLxyErr;
	fchannel    = abs(lid1*lid2);
	fsecvtxmass = leadingSecVtxMass;
	fsecvtxntrk = leadingSecVtxNtrk;
	fnvtx       = nvtx;
	fweight     = weight;

	tree_->Fill();
}

void LxyAnalysis::initializeTreeVariables() {
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


void LxyAnalysis::clearTreeVariables() {
	flxy        = 99.99;
	flxyerr     = -99.99;
	fnvtx       = -1;
	fweight     = 1.;
	fsecvtxmass = -1.;
	fsecvtxntrk = -1.;
	fchannel    = 0;
	// fcharge = -99;
}


void LxyAnalysis::writeTree(TFile &f) {
	f.cd();
	tree_->Print();
	tree_->Write();
}
