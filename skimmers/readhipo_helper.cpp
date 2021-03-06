#include "readhipo_helper.h"


void getNeutronInfo( BBand band_hits, hipo::bank band_rawhits, hipo::bank band_adc, hipo::bank band_tdc,
			int& mult, bandhit hits[maxNeutrons],
			double starttime , int thisRun){

	if( band_hits.getRows() > maxNeutrons ) return; // not interested in events with more than max # BAND hits for now
	for( int hit = 0 ; hit < band_hits.getRows() ; hit++ ){
		//if( band_hits.getStatus(hit) != 0 ) continue;	// not interested in an event that has a veto hit

		// Set the hits information
		hits[hit].setSector		(band_hits.getSector		(hit)			);
		hits[hit].setLayer		(band_hits.getLayer		(hit)			);
		hits[hit].setComponent		(band_hits.getComponent		(hit)			);
		hits[hit].setBarID		(band_hits.getBarKey		(hit)			);
		hits[hit].setEdep		(band_hits.getEnergy		(hit)			);
		hits[hit].setTof		(band_hits.getTime		(hit) - starttime	);
		hits[hit].setTofFadc		(band_hits.getTimeFadc		(hit) - starttime	);
		hits[hit].setTdiff		(band_hits.getDifftimeTdc	(hit)			);
		hits[hit].setTdiffFadc		(band_hits.getDifftimeFadc	(hit)			);
		hits[hit].setX			(band_hits.getX			(hit)			);
		hits[hit].setY			(band_hits.getY			(hit)			);
		hits[hit].setZ			(band_hits.getZ			(hit) - VERTEX_OFFSET	);
		hits[hit].setStatus		(band_hits.getStatus		(hit)			);


		// Using the band hit struct, get the raw hit PMT information to use later
		int rawhit_idxL = band_hits.getLpmtindex(hit);
		int rawhit_idxR = band_hits.getRpmtindex(hit);
		// 	Get the raw hit information corresponding to the band hit above
		hits[hit].setRawLtdc		(band_rawhits.getFloat( 7 , rawhit_idxL ) 		);
        	hits[hit].setRawRtdc		(band_rawhits.getFloat( 7 , rawhit_idxR ) 		);
        	hits[hit].setRawLtdccorr	(band_rawhits.getFloat( 9 , rawhit_idxL ) 		);
        	hits[hit].setRawRtdccorr	(band_rawhits.getFloat( 9 , rawhit_idxR ) 		);
        	hits[hit].setRawLtfadc		(band_rawhits.getFloat( 8 , rawhit_idxL ) 		);
        	hits[hit].setRawRtfadc		(band_rawhits.getFloat( 8 , rawhit_idxR ) 		);
        	hits[hit].setRawLamp		(band_rawhits.getFloat( 6 , rawhit_idxL )		);
        	hits[hit].setRawRamp		(band_rawhits.getFloat( 6 , rawhit_idxR )		);
        	hits[hit].setRawLadc		(band_rawhits.getFloat( 5 , rawhit_idxL )		);
        	hits[hit].setRawRadc		(band_rawhits.getFloat( 5 , rawhit_idxR )		);

		// Using the rawhit struct, get the raw PMT information to use later
		int pmtTdcL	= band_rawhits.getInt( 10 , rawhit_idxL );
		int pmtAdcL	= band_rawhits.getInt( 11 , rawhit_idxL );
		int pmtTdcR	= band_rawhits.getInt( 10 , rawhit_idxR );
		int pmtAdcR	= band_rawhits.getInt( 11 , rawhit_idxR );
		//	Get the raw pmt information corresponding to the band hit above
		hits[hit].setPmtLtdc		(band_tdc.getInt( 4 , pmtTdcL )		);
		hits[hit].setPmtRtdc		(band_tdc.getInt( 4 , pmtTdcR )		);
		hits[hit].setPmtLtfadc		(band_adc.getFloat( 6 , pmtAdcL )	);
		hits[hit].setPmtRtfadc		(band_adc.getFloat( 6 , pmtAdcR )	);
		hits[hit].setPmtLamp		(band_adc.getInt( 5 , pmtAdcL )		);
		hits[hit].setPmtRamp		(band_adc.getInt( 5 , pmtAdcR )		);
		hits[hit].setPmtLadc		(band_adc.getInt( 4 , pmtAdcL )		);
		hits[hit].setPmtRadc		(band_adc.getInt( 4 , pmtAdcR )		);
		hits[hit].setPmtLped		(band_adc.getInt( 7 , pmtAdcL )		);
		hits[hit].setPmtRped		(band_adc.getInt( 7 , pmtAdcR )		);

		// Save how many neutron hits we have
		mult++;
	}

}

bool goodNeutronEvent(bandhit hits[maxNeutrons], int nMult, int& leadindex, int mcdataselect){

	double thres = thresBANDhit; //5
	double time_thres = time_thresBANDhit; //300

	bool vetoHit = false;
	bool thresPass = false;
	bool accept = false;
	double adctoMeVee = 0;
	if (mcdataselect == 0) //MC
		adctoMeVee = adctoMeVee_sim;
	else  { //Data
		adctoMeVee = adctoMeVee_data;
	}

	// Loop over all hits in BAND and get the earliest hit in time
	// but filter out hits that have low Edep
	double fastestTime = 1E5;
	int fastestTimeIdx = -1;
	for( int thishit = 0; thishit < nMult ; thishit++){
			bandhit neutron = hits[thishit];
			//if the Edep in this PMT is less than 2MeVee, do not count it
			//MC
			if (mcdataselect == 0 && neutron.getPmtLadc()/adctoMeVee < 2) 	continue;
			//Data
			if( mcdataselect != 0 && neutron.getEdep()/adctoMeVee < 2 ) continue;

			if( neutron.getTof() < fastestTime ){
					fastestTime = neutron.getTof();
					fastestTimeIdx = thishit;
			}
	}

	// Once we have the fastest hit, ask how many other hits there
	// are in the same event and flag how far in time they are,
	// but keep this earliest hit.
	// Accept the event based on if there is another hit > 300ps
	if( fastestTimeIdx != -1 ){
			accept = true;
			for( int thishit = 0; thishit < nMult ; thishit++){

						bandhit neutron = hits[thishit];
						// if we have a veto bar hit and if the Edep is > 0.25MeVee,  need to use getPmtLADC because no Edep
						if( neutron.getLayer() == 6 && neutron.getPmtLadc()/adctoMeVee > 0.25 ) vetoHit = true;

						//if this bar does not have > 2MeVee, do not count the hit
						//MC
						if (mcdataselect == 0 && neutron.getPmtLadc()/adctoMeVee < 2 ) continue;
						//Data
						if( mcdataselect != 0 && neutron.getEdep()/adctoMeVee < 2 ) continue;



						// if this bar does not have > X MeVee, do not count the hit (X = threshold set by user)
						//MC
						if (mcdataselect == 0 && neutron.getPmtLadc()/adctoMeVee > thres ) thresPass = true;
						//Data
						if( mcdataselect != 0 && neutron.getEdep()/adctoMeVee > thres ) thresPass = true;

						double tdiff = neutron.getTof() - fastestTime;
						if( tdiff == 0 ){
							continue;
						}
						bool adjLayer = false;
						bool adjComp = false;
						bool adjSpace = false;
						bandhit lead = hits[thishit];

						int layerDiff = lead.getLayer() - neutron.getLayer();
						if( fabs(layerDiff) <= 1 ) adjLayer = true;
						if( lead.getSector() == neutron.getSector() ){
						// if hits are in the same sector, then just ask for comp diff
							int compDiff = lead.getComponent() - neutron.getComponent();
							if( fabs(compDiff) <= 1 ) adjComp = true;
						}
						else{
						// have to implement some fancy component check due to sec comp differences
							int yDiff = lead.getY() - neutron.getY();
							if( fabs(yDiff) <= 10 ) adjComp = true;
						}
						if( adjComp && adjLayer ) adjSpace = true;

						if( tdiff < time_thres && tdiff != 0 ) accept = false;
						if( tdiff < time_thres && tdiff != 0 && adjSpace ) accept = true;
			}
		}

		leadindex = fastestTimeIdx ;

		return accept && !vetoHit && thresPass;

}




int getRunNumber( string filename ){
	string parsed = filename.substr( filename.find("inc_") );
	string moreparse = parsed.substr(4,6);
	cout << "\t*Intepreted run number from file name: " << stoi(moreparse) << "\n";
        return stoi(moreparse);
}

void getEventInfo( BEvent eventInfo, double &integrated_charge, double &livetime, double &starttime ){
	if( eventInfo.getRows() != 1 ){
		cerr << "getEventInfo::NotImplementedFunction\n";
		exit(-1);
	}
	integrated_charge       = (double)eventInfo.getBCG(0); 	// not calibrated currently
	livetime 		= (double)eventInfo.getLT(0);		// not calibrated currently
	starttime		= (double)eventInfo.getSTTime(0);
	return;
}

void getMcInfo( hipo::bank gen_particles , hipo::bank gen_info , genpart mcParts[maxGens] ,
		double &starttime, double &weight, double &Ebeam , int &genMult ){
	TVector3 	beamVec(0,0,Ebeam);
	TVector3	eVec;
	bool setElectron = false;

	// Grab the weight for the event:
	weight 	= gen_info.getFloat(9,0);
	// Grab the beam energy for this generated file:
	Ebeam 	= gen_info.getFloat(6,0);

	// Loop over all generated particles and find the electron to put that one first
	for( int hit = 0 ; hit < gen_particles.getRows() ; hit++ ){
		int pid	= gen_particles.getInt( 0 , hit );
		if( pid == 11 && setElectron == false ){
			double px = gen_particles.getFloat( 1 , hit );
			double py = gen_particles.getFloat( 2 , hit );
			double pz = gen_particles.getFloat( 3 , hit );
			double vx = gen_particles.getFloat( 4 , hit );
			double vy = gen_particles.getFloat( 5 , hit );
			double vz = gen_particles.getFloat( 6 , hit );
			double vt = gen_particles.getFloat( 7 , hit );
			eVec.SetXYZ( px, py, pz );
			TVector3	vertex; vertex.SetXYZ( vx, vy, vz );
			TVector3	qVec; qVec = beamVec - eVec;
			starttime = vt;

			mcParts[genMult].setPID		( pid 							);
			mcParts[genMult].setMomentum	( eVec.Mag()					);
			mcParts[genMult].setTheta	( eVec.Theta()					);
			mcParts[genMult].setPhi		( eVec.Phi()					);

			mcParts[genMult].setQ		( qVec.Mag()						);
			mcParts[genMult].setThetaQ	( qVec.Theta()						);
			mcParts[genMult].setPhiQ	( qVec.Phi()						);

			mcParts[genMult].setOmega	( Ebeam - sqrt( pow(eVec.Mag(),2) + mE*mE )		);
			mcParts[genMult].setQ2		( qVec.Mag()*qVec.Mag() - pow(mcParts[genMult].getOmega(),2)	);
			mcParts[genMult].setXb		( mcParts[genMult].getQ2()/(2.*mP*mcParts[genMult].getOmega())		);
			mcParts[genMult].setW2		( mP*mP - mcParts[genMult].getQ2() + 2.*mcParts[genMult].getOmega()*mP	);

			mcParts[genMult].setVtx		(vertex.X()						);
			mcParts[genMult].setVty		(vertex.Y()						);
			mcParts[genMult].setVtz		(vertex.Z()						);

			setElectron = true;
			genMult++;
		}
		else if( pid != 11 && setElectron == true ){
			double px = gen_particles.getFloat( 1 , hit );
			double py = gen_particles.getFloat( 2 , hit );
			double pz = gen_particles.getFloat( 3 , hit );
			double vx = gen_particles.getFloat( 4 , hit );
			double vy = gen_particles.getFloat( 5 , hit );
			double vz = gen_particles.getFloat( 6 , hit );
			double vt = gen_particles.getFloat( 7 , hit );

			TVector3	momentum; 	momentum.SetXYZ( px, py, pz );
			TVector3	vertex; 	vertex.SetXYZ( vx, vy, vz );

			mcParts[genMult].setPID		( pid 							);
			mcParts[genMult].setMomentum	( momentum.Mag()					);
			mcParts[genMult].setTheta	( momentum.Theta()					);
			mcParts[genMult].setPhi		( momentum.Phi()					);

			mcParts[genMult].setVtx		(vertex.X()						);
			mcParts[genMult].setVty		(vertex.Y()						);
			mcParts[genMult].setVtz		(vertex.Z()						);

			genMult++;
		}
	}



	return;
}


void getElectronInfo( BParticle particles, BCalorimeter calorimeter, BScintillator scintillator, hipo::bank DC_Track, hipo::bank DC_Traj,

			clashit &electron,
			double starttime , int thisRun , double Ebeam ){

	TVector3	momentum = particles.getV3P(0);
	TVector3	vertex	 = particles.getV3v(0);

	TVector3 	beamVec(0,0,Ebeam);
	TVector3	qVec; qVec = beamVec - momentum;


	electron.setSector		(	calorimeter.getSector(0)				);
	electron.setPID			(	particles.getPid(0)					);
	electron.setCharge		(	particles.getCharge(0)					);
	electron.setStatus		(	particles.getStatus(0)					);

	electron.setTime		(	particles.getVt(0)					);
	electron.setBeta		(	particles.getBeta(0)					);
	electron.setChi2		(	particles.getChi2pid(0)					);
	electron.setEtot		(	calorimeter.getTotE(0)					);
	electron.setEpcal		(	calorimeter.getPcalE(0)					);
	electron.setEecin		(	calorimeter.getECinE(0)					);
	electron.setEecout		(	calorimeter.getECoutE(0)       				);
	electron.setEoP			(	electron.getEtot() / momentum.Mag()			);
	electron.setTimeScint		(	scintillator.getTime(0)-starttime			);
	electron.setPathScint		(	scintillator.getPath(0)					);
	electron.setU			(	calorimeter.getLU(0)					);
	electron.setV			(	calorimeter.getLV(0)					);
	electron.setW			(	calorimeter.getLW(0)					);
	electron.setVtx			(	vertex.X()						);
	electron.setVty			(	vertex.Y()						);
	electron.setVtz			(	vertex.Z()						);

	electron.setMomentum		(	momentum.Mag()						);
	electron.setTheta		(	momentum.Theta()					);
	electron.setPhi			(	momentum.Phi()						);

	electron.setQ			(	qVec.Mag()						);
	electron.setThetaQ		(	qVec.Theta()						);
	electron.setPhiQ		(	qVec.Phi()						);

	electron.setOmega		(	Ebeam - sqrt( pow(momentum.Mag(),2) + mE*mE )		);
	electron.setQ2			(	qVec.Mag()*qVec.Mag() - pow(electron.getOmega(),2)	);
	electron.setXb			(	electron.getQ2()/(2.*mP*electron.getOmega())		);
	electron.setW2			(	mP*mP - electron.getQ2() + 2.*electron.getOmega()*mP	);


        //Adding Tracking bank

	int Ntrack = DC_Track.getRows();
	int index = -1;
	int count =0;
	for(int i =0; i < Ntrack; i++){
	  int pindex   = DC_Track.getInt(1,i);
	  int detector = DC_Track.getInt(2,i);

	  if (pindex ==0 && detector ==6 )
	    {index = i;
	     count ++;
	   }
	}

	//At this point I only use 1 track event
	//The check show around 4% mul_track events

	if (index != -1 && count ==1) {
	  electron.setDC_chi2             (       DC_Track.getFloat(6, index)                                 );
	  electron.setDC_NDF              (       DC_Track.getInt(7, index)                                   );
          electron.setDC_sector           (       DC_Track.getInt(3, index)                                   );

	}


	//Adding Traj bank, Structure is more complicated
	//Also only consider the traj of 1-track event

	int NTraj = DC_Traj.getRows();
	int index_traj = -1;

	for(int i =0; i <NTraj; i++){
          int pindex   = DC_Traj.getInt(0, i);
	  int detector = DC_Traj.getInt(2,i);
	  int layer    = DC_Traj.getInt(3,i);

	  if(pindex ==0 && detector ==6 && layer ==6 && count ==1)

	     {
	       electron.setDC_x1        (  DC_Traj.getFloat(4, i)                                           );
	       electron.setDC_y1        (  DC_Traj.getFloat(5, i)                                           );
	       electron.setDC_z1        (  DC_Traj.getFloat(6, i)                                           );

	     }

	  if(pindex ==0 && detector == 6 && layer ==18 && count ==1)

	    {
	       electron.setDC_x2        (  DC_Traj.getFloat(4, i)                                           );
	       electron.setDC_y2        (  DC_Traj.getFloat(5, i)                                           );
	       electron.setDC_z2        (  DC_Traj.getFloat(6, i)                                           );


	    }

	  if(pindex ==0 && detector == 6 && layer ==36 && count ==1)

	    {
	       electron.setDC_x3        (  DC_Traj.getFloat(4, i)                                           );
	       electron.setDC_y3        (  DC_Traj.getFloat(5, i)                                           );
	       electron.setDC_z3        (  DC_Traj.getFloat(6, i)                                           );

	    }


	}


	//if( pid != 11 || charge != -1 ) return false;
	//if( lV < 2 || lW < 2 ) return false;
	//if( momentum.Mag() < 1 || momentum.Mag() > 4.2 ) return false;
	//if( chi2pid == 0 || chi2pid > 1 ) return false;
	//if( E_tot/momentum.Mag() > 0.4 || E_tot/momentum.Mag() < 0.1 ) return false;
	//if( vertex.X() < -2 || vertex.X() > 2) return false;
	//if( vertex.Y() < -2 || vertex.Y() > 2) return false;
	//if( vertex.Z() < -7 || vertex.Z() > 2) return false;
	//if( time < 15 ) return false;
	//if( momentum.Mag() < 2 || momentum.Mag() > 10.6 ) return false;
	//if( E_tot / momentum.Mag() < 0.15 || E_tot / momentum.Mag() > 0.3 ) return false;

	return;
}

void getTaggedInfo( clashit eHit, bandhit nHit[maxNeutrons], taghit tag[maxNeutrons] ,
		double Ebeam , int nMult ){

	TVector3 	beamVec(0,0,Ebeam);
	TVector3	eVec; eVec.SetMagThetaPhi( eHit.getMomentum(), eHit.getTheta(), eHit.getPhi() );
	TVector3	qVec; qVec = beamVec - eVec;

	// Loop over all neutrons to combine with the electron
	for( int hit = 0 ; hit < nMult ; hit++ ){

		TVector3	nVec;
		nVec = (nHit[hit].getDL()).Unit();


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

		double beta = nHit[hit].getDL().Mag() / (nHit[hit].getTofFadc()*cAir);
		double p_n = mN / sqrt( 1./pow(beta,2) - 1. );
		nVec.SetMagThetaPhi(p_n,nHit[hit].getDL().Theta(),nHit[hit].getDL().Phi());

		double E_n 	= sqrt( mN*mN + p_n*p_n );
		double W_primeSq = mD*mD - eHit.getQ2() + mN*mN + 2.*mD*(eHit.getOmega()-E_n) - 2.*eHit.getOmega()*E_n + 2.*eHit.getQ()*p_n*cos(theta_nq);
		double Wp = sqrt(W_primeSq);
		double Xp = eHit.getQ2()/(2.*( eHit.getOmega()*(mD-E_n) + p_n*eHit.getQ()*CosTheta_nq));
		double As = (E_n - p_n*CosTheta_nq)/mN;
		double Xp2 = eHit.getQ2()/(W_primeSq - mN*mN + eHit.getQ2());

		TVector3 Pt;
		TVector3 pN_par_q = nVec.Dot(qVec) / (qVec.Mag2()) * qVec;
		Pt = nVec - pN_par_q;

		tag[hit].setMomentumE	(eVec 		);
		tag[hit].setMomentumN	(nVec		);
		tag[hit].setMomentumQ	(qVec		);
		tag[hit].setMomentumB	(beamVec	);

		tag[hit].setPhiNQ	(phi_nq		);
		tag[hit].setThetaNQ	(theta_nq	);
		tag[hit].setWp		(Wp		);
		tag[hit].setXp		(Xp		);
		tag[hit].setAs		(As		);
		tag[hit].setPt		(Pt		);
		tag[hit].setXp2		(Xp2		);
	}

	//beta = dL / (ToF*cAir);
	//p_n = mN / sqrt( 1./pow(beta,2) - 1. );
	//nVec.SetMagThetaPhi(p_n,theta_n,phi_n);
	//E_n 	= sqrt( mN*mN + p_n*p_n );
	//double W_primeSq = mD*mD - Q2 + mN*mN + 2.*mD*(nu-E_n) - 2.*nu*E_n + 2.*q*p_n*cos(theta_nq);
	//Wp = sqrt(W_primeSq);
	//Xp = Q2/(2.*( nu*(mD-E_n) + p_n*q*CosTheta_nq));
	//As = (E_n - p_n*CosTheta_nq)/mN;


	return;
}

void shiftsReader::LoadInitBar( string filename ){
	ifstream f;
	int sector, layer, component, barId;
	double pol0, height, mean, sig, temp;

	f.open(filename);
	while(!f.eof()){
		f >> sector;
		f >> layer;
		f >> component;
		barId = 100*sector + 10*layer + component;
		f >> pol0;
		f >> height;
		f >> mean;
		f >> sig;
		InitBar[barId] = mean;
		f >> temp;
		f >> temp;
	}
	f.close();
}
void shiftsReader::LoadInitBarFadc( string filename ){
	ifstream f;
	int sector, layer, component, barId;
	double pol0, height, mean, sig, temp;

	f.open(filename);
	while(!f.eof()){
		f >> sector;
		f >> layer;
		f >> component;
		barId = 100*sector + 10*layer + component;
		f >> pol0;
		f >> height;
		f >> mean;
		f >> sig;
		InitBarFadc[barId] = mean;
		f >> temp;
		f >> temp;
	}
	f.close();
}
void shiftsReader::LoadInitRun( string filename ){
	ifstream f;
	int runno;
	double pol0, height, mean, sig, temp;

	f.open(filename);
	while(!f.eof()){
		f >> runno;
		f >> pol0;
		f >> height;
		f >> mean;
		f >> sig;
		InitRun[runno] = mean;
		f >> temp;
		f >> temp;
	}
	f.close();
}
void shiftsReader::LoadInitRunFadc( string filename ){
	ifstream f;
	int runno;
	double pol0, height, mean, sig, temp;

	f.open(filename);
	while(!f.eof()){
		f >> runno;
		f >> pol0;
		f >> height;
		f >> mean;
		f >> sig;
		InitRunFadc[runno] = mean;
		f >> temp;
		f >> temp;
	}
	f.close();
}

double * shiftsReader::getInitBar(void){
	return InitBar;
}
double * shiftsReader::getInitBarFadc(void){
	return InitBarFadc;
}
double * shiftsReader::getInitRun(void){
	return InitRun;
}
double * shiftsReader::getInitRunFadc(void){
	return InitRunFadc;
}
