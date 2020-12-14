#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <vector>

#include "TFile.h"
#include "TTree.h"
#include "TVector3.h"

#include "constants.h"

#include "clashit.h"

using namespace std;


int main(int argc, char** argv) {
	// check number of arguments
	if( argc < 3 ){
		cerr << "Incorrect number of arugments. Instead use:\n\t./code [outputFile] [MC/DATA] [inputFile] \n\n";
		cerr << "\t\t[outputFile] = ____.root\n";
		cerr << "\t\t[beamEnergy] = __ GeV\n";
		cerr << "\t\t[inputFile] = ____.root ____.root ____.root ...\n\n";
		return -1;
	}

	double beamE = atof(argv[2]);
	TVector3 k0 = TVector3(0.,0.,sqrt(beamE*beamE - mE*mE)); 

	// Create output tree
	TFile * outFile = new TFile(argv[1],"RECREATE");
	TTree * outTree = new TTree("electrons","CLAS Electrons");
	//	Event info:
	int Runno		= 0;
	double Ebeam		= beamE;
	double gated_charge	= 0;
	double livetime		= 0;
	double starttime	= 0;
	double current		= 0;
	double weight		= 0;
	//	MC info:
	//	Electron info:
	clashit eHit;
	// 	Event branches:
	outTree->Branch("Runno"		,&Runno			);
	outTree->Branch("Ebeam"		,&Ebeam			);
	outTree->Branch("gated_charge"	,&gated_charge		);
	outTree->Branch("livetime"	,&livetime		);
	outTree->Branch("starttime"	,&starttime		);
	outTree->Branch("current"	,&current		);
	outTree->Branch("weight"	,&weight		);
	//	Electron branches:
	outTree->Branch("eHit"		,&eHit			);

	// Load input file
	for( int i = 3 ; i < argc ; i++ ){


		// Read in ROOT file
		TString inputFile = argv[i];
		TFile* inFile = new TFile(inputFile);
		TTree* T = (TTree*)inFile->Get("T");

		int Nevents = T->GetEntries();

		double pe[3];
		T->SetBranchAddress("pe", pe);

		cout << "Working on file " << inputFile << endl;
		
		for(int j = 0; j < Nevents; j++) {

			eHit.Clear();
		
			T->GetEntry(j);
			TVector3 electron = TVector3(pe[0], pe[1], pe[2]);

			TVector3 q = k0 - electron;

			double E_e = sqrt(electron.Mag2() + mE*mE);
			double omega = beamE - E_e;								
			double Q2 = q.Mag2() - omega*omega;
			double xB = Q2/(2*mP*omega);
			double W2 = mP*mP + 2.*omega*mP	- Q2;
		

			eHit.setPID(11);		
			eHit.setCharge(-1);
						
			eHit.setMomentum(electron.Mag());	
			eHit.setTheta(electron.Theta());	
			eHit.setPhi(electron.Phi());
						
			eHit.setQ(q.Mag());		
			eHit.setThetaQ(q.Theta());
			eHit.setPhiQ(q.Phi());
						
			eHit.setOmega(omega);	
			eHit.setQ2(Q2);	
			eHit.setXb(xB);		
			eHit.setW2(W2);	

			outTree->Fill();

		} // end loop over events
	}// end loop over files
	
	outFile->cd();
	outTree->Write();
	outFile->Close();

	return 0;
}



