#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <vector>

#include "TFile.h"
#include "TTree.h"
#include "TVector3.h"
#include "TH2.h"
#include "TH1.h"
#include "TClonesArray.h"
#include "TRandom3.h"

#include "constants.h"
#include "readhipo_helper.h"

using namespace std;


int main(int argc, char** argv) {
	// check number of arguments
	if( argc < 5 ){
		cerr << "Incorrect number of arugments. Instead use:\n\t./code [outputFile] [MC/DATA] [inputFile] \n\n";
		cerr << "\t\t[outputFile] = ____.root\n";
		cerr << "\t\t[bgFile]  = ____.root\n";
		cerr << "\t\t[incFile] = ____.root\n";
		cerr << "\t\t[sigFile] = ____.root\n\n";
		return -1;
	}

	TRandom3* fRand = new TRandom3(0);

	TFile* bgFile = new TFile(argv[2]);
	TTree* bgTree = (TTree*)bgFile->Get("neutrons");

	TFile* incFile = new TFile(argv[3]);
	TTree* incTree = (TTree*)incFile->Get("electrons");

	TFile* sigFile = new TFile(argv[4]);
	TTree* sigTree = (TTree*)sigFile->Get("tagged");

	// Create output tree
	TFile * outFile = new TFile(argv[1],"RECREATE");
	TTree * outTree = new TTree("tagged","BAND Neutrons and CLAS Electrons");
	//	Event info:
	int Runno		= 11;
	double Ebeam		= 10.2;
	double gated_charge	= 0;
	double livetime		= 0;
	double starttime	= 0;
	double current		= 0;
	bool goodneutron = false;
	int nleadindex = -1;
	double weight		= 0;
	// 	Neutron info:
	int nMult		= 0;
	TClonesArray * nHits = new TClonesArray("bandhit");
	TClonesArray &saveHit = *nHits;
	//	Electron info:
	clashit* eHit = new clashit;
	//	Tagged info:
	TClonesArray * tags = new TClonesArray("taghit");
	TClonesArray &saveTags = *tags;
	// Background info
	int background;
	// 	Event branches:
	outTree->Branch("Runno"		,&Runno			);
	outTree->Branch("Ebeam"		,&Ebeam			);
	outTree->Branch("gated_charge"	,&gated_charge		);
	outTree->Branch("livetime"	,&livetime		);
	outTree->Branch("starttime"	,&starttime		);
	outTree->Branch("current"	,&current		);
	outTree->Branch("weight"	,&weight		);
	//	Neutron branches:
	outTree->Branch("nMult"		,&nMult			);
	outTree->Branch("nHits"		,&nHits			);
	//Branches to store if good Neutron event and leadindex
	outTree->Branch("goodneutron"		,&goodneutron	);
	outTree->Branch("nleadindex"		,&nleadindex			);
	//	Electron branches:
	outTree->Branch("eHit"		,&eHit			);
	//	Tagged branches:
	outTree->Branch("tag"		,&tags			);
	outTree->Branch("background",	&background		);


	bgTree->SetBranchAddress("Runno"		,&Runno			);
	bgTree->SetBranchAddress("Ebeam"		,&Ebeam			);
	bgTree->SetBranchAddress("gated_charge"		,&gated_charge		);
	bgTree->SetBranchAddress("livetime"		,&livetime		);
	bgTree->SetBranchAddress("starttime"		,&starttime		);
	bgTree->SetBranchAddress("current"		,&current		);
	bgTree->SetBranchAddress("weight"		,&weight		);
	
	bgTree->SetBranchAddress("nMult",	&nMult		);
	bgTree->SetBranchAddress("nHits",	&nHits		);
	bgTree->SetBranchAddress("goodneutron",	&goodneutron	);
	bgTree->SetBranchAddress("nleadindex",	&nleadindex	);

	incTree->SetBranchAddress("eHit",	&eHit	);


	int nBG = min(bgTree->GetEntries(), incTree->GetEntries());
	
	// Add background events to output tree
	for(int i = 0; i < nBG; i++) {

		if(i%100000==0){
			cout << "Background event " << i << "/" << nBG << endl;
		}		

		background = 1;

		// Clear all branches
		gated_charge	= 0;
		livetime	= 0;
		starttime 	= 0;
		// Neutron
		nMult		= 0;
		nleadindex = -1;
		goodneutron = false;

		nHits->Clear();
		eHit->Clear();
		tags->Clear();

		bgTree->GetEntry(i);
		incTree->GetEntry(i);

		bandhit these_nHits[maxNeutrons];
		taghit these_tags[maxNeutrons];	

		for(int n = 0; n < nMult; n++) {
			bandhit* nHit_TClones = (bandhit*)nHits->At(n);
			bandhit this_nHit = *nHit_TClones;
			these_nHits[n] = this_nHit;
		}		
	
		TVector3 	beamVec(0,0,Ebeam);
		TVector3	eVec; eVec.SetMagThetaPhi( eHit->getMomentum(), eHit->getTheta(), eHit->getPhi() );
		TVector3	qVec; qVec = beamVec - eVec;

		// Set random TOF
		double min_tof = -50.;
		double max_tof = 100.; 
		double bg_tof = min_tof + fRand->Rndm()*(max_tof - min_tof);
		
		these_nHits[nleadindex].setTof(bg_tof);
		these_nHits[nleadindex].setTofFadc(bg_tof);

		// Loop over all neutrons to combine with the electron
		for( int hit = 0 ; hit < nMult ; hit++ ){

			TVector3	nVec;
			nVec = (these_nHits[hit].getDL()).Unit();


			TVector3 norm_scatter = qVec.Cross( beamVec );
			norm_scatter 	= norm_scatter.Unit();

			TVector3 norm_reaction = qVec.Cross( nVec );
			norm_reaction 	= norm_reaction.Unit();

			double phi_nq 		= norm_scatter.Angle( norm_reaction );
			double theta_nq 	= nVec.Angle( qVec );
			double CosTheta_nq 	= cos(theta_nq);

			TVector3 direction = norm_scatter.Cross(norm_reaction);
			if( direction.Z() > 0 ){ // this means the phi_rq should be between 0 and pi
			}
			else if( direction.Z() < 0 ){ // this means the phi_rq should be between -pi and 0
				phi_nq *= (-1);
			}

			double beta = these_nHits[hit].getDL().Mag() / (these_nHits[hit].getTofFadc()*cAir);
			double p_n = mN / sqrt( 1./pow(beta,2) - 1. );
			nVec.SetMagThetaPhi(p_n,these_nHits[hit].getDL().Theta(),these_nHits[hit].getDL().Phi());

			double E_n 	= sqrt( mN*mN + p_n*p_n );
			double W_primeSq = mD*mD - eHit->getQ2() + mN*mN + 2.*mD*(eHit->getOmega()-E_n) - 2.*eHit->getOmega()*E_n + 2.*eHit->getQ()*p_n*cos(theta_nq);
			double Wp = sqrt(W_primeSq);
			double Xp = eHit->getQ2()/(2.*( eHit->getOmega()*(mD-E_n) + p_n*eHit->getQ()*CosTheta_nq));
			double As = (E_n - p_n*CosTheta_nq)/mN;
			double Xp2 = eHit->getQ2()/(W_primeSq - mN*mN + eHit->getQ2());

			TVector3 Pt;
			TVector3 pN_par_q = nVec.Dot(qVec) / (qVec.Mag2()) * qVec;
			Pt = nVec - pN_par_q;

			these_tags[hit].setMomentumE	(eVec 		);
			these_tags[hit].setMomentumN	(nVec		);
			these_tags[hit].setMomentumQ	(qVec		);
			these_tags[hit].setMomentumB	(beamVec	);

			these_tags[hit].setPhiNQ	(phi_nq		);
			these_tags[hit].setThetaNQ	(theta_nq	);
			these_tags[hit].setWp		(Wp		);
			these_tags[hit].setXp		(Xp		);
			these_tags[hit].setAs		(As		);
			these_tags[hit].setPt		(Pt		);
			these_tags[hit].setXp2		(Xp2		);
		}

		for(int n = 0; n < nMult; n++) {
			new(saveHit[n]) bandhit;
			saveHit[n] = &these_nHits[n];
			
			new(saveTags[n]) taghit;
			saveTags[n] = &these_tags[n];
		}		

		outTree->Fill();

	} 	

	bgFile->Close();	
	incFile->Close();


	// Add signal events to output tree

	sigTree->SetBranchAddress("Runno"		,&Runno			);
	sigTree->SetBranchAddress("Ebeam"		,&Ebeam			);
	sigTree->SetBranchAddress("gated_charge"	,&gated_charge		);
	sigTree->SetBranchAddress("livetime"		,&livetime		);
	sigTree->SetBranchAddress("starttime"		,&starttime		);
	sigTree->SetBranchAddress("current"		,&current		);
	sigTree->SetBranchAddress("weight"		,&weight		);
	
	sigTree->SetBranchAddress("nMult",	&nMult		);
	sigTree->SetBranchAddress("nHits",	&nHits		);
	sigTree->SetBranchAddress("goodneutron",&goodneutron	);
	sigTree->SetBranchAddress("nleadindex",	&nleadindex	);

	sigTree->SetBranchAddress("eHit",	&eHit		);

	sigTree->SetBranchAddress("tag",	&tags		);

	int nSig = sigTree->GetEntries();

	for(int i = 0; i < nSig; i++) {

		if(i%100000==0){
			cout << "Signal event " << i << "/" << nSig << endl;
		}		
		background = 0;
		
		sigTree->GetEntry(i);

		outTree->Fill();

	}

	outFile->cd();
	outTree->Write();
	outFile->Close();

	return 0;
}
