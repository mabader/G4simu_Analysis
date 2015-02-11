#include<math.h>
#include <vector>
#include <string>
#include <iostream> 


void SimpleAnalysis(string filename)
{
	// Open file
	TFile *f = TFile::Open(("/home/geant4/G4Projects/G4simu/macros/" + filename + ".root").c_str());

	//step into the subdirectory 'events' to acess the tree
	TDirectory* dir = gFile->GetDirectory("events");	
	TTree *tree = (TTree *) dir->Get("events");

	// Connect tree variables to tree entries
	Float_t E; 			// for total energy deposit ( sum of energy deposit in both detectors)
	vector<int> *vcollID = 0; 	// for collection ID
	vector<float> *ved = 0; 	// for deposited energy
	vector<float> *vxp = 0;		// x pos energy deposit
	vector<float> *vyp = 0;
	vector<float> *vzp = 0;
	
	tree -> SetBranchAddress("etot", &E);
	tree -> SetBranchAddress("collid", &vcollID);  
	tree -> SetBranchAddress("ed", &ved);
	tree -> SetBranchAddress("xp", &vxp);
	tree -> SetBranchAddress("yp", &vyp);
	tree -> SetBranchAddress("zp", &vzp);

	// Geometry parameters from tree:
	TDirectory *dir_geo = gFile ->GetDirectory("detector/geometry/");
	TParameter<double> *pos_NaI = (TParameter<double>*) dir_geo ->Get("NaI_crystal_Position");
	TParameter<double> *pos_Source = (TParameter<double>*) dir_geo ->Get("SourceDisk_Position");
	TParameter<double> *NaI_z = (TParameter<double>*) dir_geo ->Get("NaI_crystal_Z");
	TParameter<double> *Coll_r = (TParameter<double>*) dir_geo ->Get("Collimator_Radius");
	TParameter<double> *Teflon_r = (TParameter<double>*) dir_geo ->Get("Teflon_Rout");
	TParameter<double> *Teflon_thickness = (TParameter<double>*) dir_geo ->Get("Teflon_Thickness");

	// Initialize histograms
	TH1F *h_etot  = new TH1F("h_etot", "Total gamma ray energy deposit; Energy (keV)",500, 0, 1500);
	TH1F *h_Edepos_LXe  = new TH1F("h_Edepos_LXe", "Energy deposit in LXe; Energy [keV]; Rate [Hz]",500, 0, 1500);
	TH1F *h_Edepos_NaI  = new TH1F("h_Edepos_NaI", "Energy deposited in NaI; Energy [keV]; Rate [Hz]",500, 0, 1500);

	TH2F *h_coincidence_xy  = new TH2F("h_coincidence_xy", "Coincidence signal in LXe and NaI xy-plane;  Energy in LXe (keV); Energy in NaI ",40,-25,25,40,-25,25);
	TH2F *h_coincidence_rz  = new TH2F("h_coincidence_rz", "Coincidence signal in LXe and NaI rz-plane;  Energy in LXe (keV); Energy in NaI ",20,-5,630,20,-60,60);

	TH2F *h_hitpattern_xy = new TH2F("h_hitpattern_xy","xy hitpattern in LXe; x [mm] ; y [mm] ",40,-25,25,40,-25,25);
	TH2F *h_hitpattern_rz = new TH2F("h_hitpattern_rz","rz hitpattern in LXe; r**2 [mm**2] ; z [mm] ",20,-5,630,20,-60,60);
	
	TH1F *h_coincidence_LXe_spectrum = new TH1F("h_coincidence_LXe_spectrum", "Na-22 coincidence spectrum in LXe with 511keV trigger in NaI; Energy [keV], Rate [Hz]", 500, 0, 1500 );
	TH1F *h_spread_z = new TH1F("h_spread_z", "Spread of 511 keV gamma rays in detector; z Position [mm]", 100, -50, 50 );


	// Some important quantities 
	const Int_t TotEvent = (Int_t)tree->GetEntries();
	double activity = 6680; 			// Bq = events/s
	double deltaT = TotEvent / activity ; 		// duration of experiment
	double weight = 1/deltaT ;
	
	double dSource_NaI = (pos_NaI-> GetVal() - NaI_z-> GetVal()/2) - pos_Source-> GetVal(); //mm
	double alpha = atan(Coll_r -> GetVal()/dSource_NaI); //mm
	double l = dSource_NaI/cos(alpha); //mm
			
	//double BeamArea = 2*3.14*l*(l-sqrt(l**2-(Coll_r -> GetVal())**2)); //mm2
	//double BeamSphere = 4*3.14*l**2 ; //mm2
	double BeamRadius = Coll_r -> GetVal() * (pos_Source -> GetVal() - ( Teflon_r -> GetVal() - Teflon_thickness -> GetVal()))/ ( pos_NaI -> GetVal() - pos_Source -> GetVal() - NaI_z ->GetVal()/2) ; //mm


	// Some counting parameters
	double w511 = 10;
	int NaI511 = 0;	
	double inside = 0.; 
	double outside = 0.;

	// Read tree and fill histos ---------------------------------------------------------------------------//
	for(Int_t i=0; i<TotEvent; i++){  //all events in a tree
		tree -> GetEntry(i);

		// vcollID -> size() = 0 -> no hit
		// vcollID -> size() = 1 -> hit in LXe (vcollID->at(j) = 0 ) or hit in NaI (vcollID->at(j) =  1)
		// vcollID -> size() = 2 -> coincidence, hit in LXe and NaI

		// Coincidence LXe and NaI, trigger on 511 keV in NaI crystal, fill histo if energy deposit in LXe is unequal to 0
		if ( vcollID -> size() == 2  && TMath::Abs(ved->at(1) - 511) < w511 && ved -> at(0) != 0){
			//Spectrum, weighted to get the rate on the y-axis
			h_coincidence_LXe_spectrum -> Fill(ved->at(0), weight);	

			// Collection of 511 keV gammas that end up in LXe
			if( (ved->at(0) - 511) < w511) {
				h_coincidence_xy -> Fill(vxp -> at(0), vyp -> at(0));
				double A_xy = (vxp -> at(0))**2 + (vyp -> at(0))**2 ;
				h_coincidence_rz -> Fill(A_xy, vzp -> at(0));
			}

			// z-position of 511 keV's in LXe
			//if(TMath::Abs(ved->at(0) - 511) < w511){
			//	if( vzp -> at(0) <= BeamRadius) inside = inside + 1.;
			//	if( vzp -> at(0) > BeamRadius) outside = outside + 1.;
 				h_spread_z -> Fill(vzp -> at(0));
			//}
		}
	

		for(int j = 0; j < vcollID -> size(); j++){
			// Hit in LXe			
			if( vcollID->at(j) == 0 && E != 0){
				// Total energy deposit in LXe, scaled with the rate
				 h_Edepos_LXe -> Fill(ved -> at(j),weight);

				// xy and rz hitpattern LXe
				// histogram entries weighted with the energy deposited
				h_hitpattern_xy -> Fill(vxp -> at(j), vyp -> at(j));
				//distance from center r = sqrt(x**2 + y**2)
				double A_xy = (vxp -> at(j))**2 + (vyp -> at(j))**2 ;
				h_hitpattern_rz -> Fill(A_xy, vzp -> at(j));
			}
			
			// Hit in NaI			
			if( vcollID->at(j) == 1 && E != 0) {
				h_Edepos_NaI -> Fill(ved -> at(j), weight);
				if (TMath::Abs(ved->at(j) - 511) < w511) NaI511 = NaI511+1; // counter	
			}		
		}
	}
	
	// Fit histogramms ----------------------------------------------------------------------------------------//
	TF1 *f1 = new TF1("f1","[0]*exp(-(x-[1])**2/(2*[2]**2))",-50,50);
	f1 -> SetParameters(50, h_spread_z->GetMean(),h_spread_z->GetRMS());
	f1 ->SetParNames("Constant","Mean_value","Sigma");
	h_spread_z->Fit("f1");	
	cout << f1-> GetChisquare()/ f1->GetNDF() << endl;


	// Draw and save stuff ----------------------------------------------------------------------------------------//
	TCanvas *c_spread = new TCanvas ("c_spread", "spread");
	c_spread -> cd();
	h_spread_z -> SetMarkerStyle(20);
	h_spread_z -> Draw("P");
	f1 -> Draw("SAME");
	
	TCanvas *cEdepos =new TCanvas("cEdepos","Energy deposit LXe and NaI");
	cEdepos -> Divide(2,1);
  	cEdepos -> cd(1);
	h_Edepos_LXe -> Draw();
	cEdepos -> cd(2);
	h_Edepos_NaI -> Draw();
	
	TCanvas *c_Spectrum = new TCanvas("c_Spectrum","Na-22 spectrum in LXe from coincidence trigger");
	c_Spectrum -> cd();
	h_coincidence_LXe_spectrum -> Draw();

	TCanvas *c_hitpattern =new TCanvas("c_hitpattern","LXe hitpattern");
	c_hitpattern -> Divide(2,1);
	c_hitpattern -> cd(1);
	h_hitpattern_xy -> Draw("colz");
	h_coincidence_xy -> SetMarkerStyle(4);
	h_coincidence_xy -> Draw("SAME"); 
	c_hitpattern -> cd(2);
	h_hitpattern_rz -> Draw("colz");
	h_coincidence_rz -> SetMarkerStyle(4);
	h_coincidence_rz -> Draw("SAME"); 

	cEdepos -> SaveAs((filename + "_energydeposit.png").c_str());
	c_hitpattern -> SaveAs((filename + "_hitpattern.png").c_str());
	c_Spectrum -> SaveAs((filename + "_spectrum.png").c_str());
	c_spread -> SaveAs((filename + "_spread_z.png").c_str());

	// output some information -------------------------------------------------------------------------------------//
	ofstream outputfile;
	outputfile.open((filename + ".txt").c_str());

	outputfile << " --------------------------------------- Event Info ---------------------------------------" << endl;
	outputfile << "root file: "<< filename << ".root" << endl;
	outputfile << endl;
	outputfile << "Number of generated events: " << TotEvent << endl;
	outputfile << "Source position: " << pos_Source-> GetVal() << " mm" << endl;
	outputfile << "Collimator radius: "<< Coll_r ->GetVal() << " mm"<< endl;
	outputfile << "Beam radius: " << BeamRadius << " mm"<< endl;
	outputfile << "Duration of run: "<< deltaT << " s"<< endl;

	outputfile << endl;
	outputfile << "---- Starting from a 511 keV trigger on NaI -----" << endl;
	outputfile << "Total Number of events in LXe :" << h_coincidence_LXe_spectrum -> GetEntries() << endl;
	outputfile << "Number of events in LXe with 511 kev: "<< h_coincidence_LXe_spectrum -> GetBinContent(h_coincidence_LXe_spectrum -> FindBin(511))*deltaT << endl;
	outputfile << "Total rate of energy deposition in LXe: " << h_coincidence_LXe_spectrum -> GetEntries()/deltaT << " Hz " << endl;
	outputfile << "Rate of energy deposition of 511 keV's in LXe: " <<  h_coincidence_LXe_spectrum -> GetBinContent(h_coincidence_LXe_spectrum -> FindBin(511)) << " Hz"<< endl;

	double activity_improved = 100000; //Bq
	double deltaT_improved = TotEvent/activity_improved ;
	outputfile << endl ;
	outputfile << "With improved activity: " << activity_improved << " Bq " << endl;
	outputfile << "Total rate of energy deposition in LXe: " << h_coincidence_LXe_spectrum -> GetEntries()/deltaT_improved << " Hz " << endl;
	outputfile << "Rate of energy deposition of 511 keV's in LXe: " <<  (h_coincidence_LXe_spectrum -> GetBinContent(h_coincidence_LXe_spectrum -> FindBin(511)) )*deltaT/deltaT_improved<< " Hz"<< endl;
	outputfile.close();


return;

	

    //datarate, , spread in LXe, coincidence plot energy in LXe
  
    
    
}
