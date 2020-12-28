#include <algorithm>
#define FloatVec std::vector<float>
void Drawplot(TString fname){
	TFile *wfile = new TFile("PlotHitmap_"+fname+".root","recreate");
	
	const int Nlayers=14;
	FloatVec* posx_step=NULL;
	FloatVec* posy_step=NULL;
	FloatVec* posz_step=NULL;
	FloatVec* E_step=NULL;

	FloatVec* posx_bar=NULL;
	FloatVec* posy_bar=NULL;
	FloatVec* posz_bar=NULL;
	FloatVec* E_bar=NULL;
	FloatVec* T1_bar=NULL;
	FloatVec* T2_bar=NULL;
	FloatVec* dlayer_bar=NULL;
	FloatVec* slayer_bar=NULL;

	FloatVec* posx_sim=NULL; 
	FloatVec* posy_sim=NULL; 
	FloatVec* posz_sim=NULL; 
	FloatVec* E_sim=NULL; 
	FloatVec* dlayer_sim=NULL; 
	FloatVec* slayer_sim=NULL; 

   FloatVec* posx_rec=NULL;
   FloatVec* posy_rec=NULL;
   FloatVec* posz_rec=NULL; 
   FloatVec* E_rec=NULL;

	TH2D *h_2dxy_step = new TH2D("h_2dxy_step", "", 100,-50,50, 300,0,300);
	TH2D *h_2dxz_step = new TH2D("h_2dxz_step", "", 200,130,330, 300,0,300);
	TH2D *h_2dyz_step = new TH2D("h_2dyz_step", "", 400,50,450, 400,-150,250);
	TH2D *h_2dmap_step[Nlayers];

	TH1D *h_posx_sim = new TH1D("h_posx_sim","",30,0,300);
	TH1D *h_posy_sim = new TH1D("h_posy_sim","",20,-100,100);
	TH1D *h_posz_sim = new TH1D("h_posz_sim","",20,130,330);
	TH2D *h_SimHit = new TH2D("h_SimHit","",40,50,450, 40,-150,250);
   TH2D *h_2dmap_sim_s0[Nlayers];
   TH2D *h_2dmap_sim_s1[Nlayers];
	TH1D *h_posy_simlayer[Nlayers];
	TH1D *h_posz_simlayer[Nlayers];

   TH1D *h_posx_rec = new TH1D("h_posx_rec","",30,0,300);
   TH1D *h_posy_rec = new TH1D("h_posy_rec","",20,-100,100);
   TH1D *h_posz_rec = new TH1D("h_posz_rec","",20,130,330);
	TH2D *h_RecHit = new TH2D("h_RecHit","",40,50,450, 40,-150,250);
	TH2D *h_2dmap_rec[Nlayers];
   TH1D *h_posy_reclayer[Nlayers];
   TH1D *h_posz_reclayer[Nlayers];

//	TH1D *h_Eres1d = new TH1D("h_Eres1d", "", 100,-10,10);
//	TH2D *h_Eres2d[Nlayers];
//	TH2D *h_Eres2d_nowi[Nlayers];
//	TH2D *h_Ediff2d[Nlayers];


	for(int i=0;i<Nlayers;i++){
		TString s_il; s_il.Form("%d",i);
		h_2dmap_step[i] = new TH2D("h_2dmap_step_"+s_il, "", 200,130,330,200,-100,100);
		h_2dmap_sim_s0[i] = new TH2D("h_2dmap_sim_s0_"+s_il, "", 20,130,330,20,-100,100);
		h_2dmap_sim_s1[i] = new TH2D("h_2dmap_sim_s1_"+s_il, "", 20,130,330,20,-100,100);
      h_2dmap_rec[i] = new TH2D("h_2dmap_rec_"+s_il, "", 20,130,330,20,-100,100);
		h_posy_simlayer[i] = new TH1D("h_posy_simlayer"+s_il,"",20,-100,100);
		h_posz_simlayer[i] = new TH1D("h_posz_simlayer"+s_il,"",20, 130,330);
		h_posy_reclayer[i] = new TH1D("h_posy_reclayer"+s_il,"",20,-100,100);
		h_posz_reclayer[i] = new TH1D("h_posz_reclayer"+s_il,"",20, 130,330);
//		h_Eres2d[i] = new TH2D("h_Eres2d"+s_il, "", 20,130,330,20,-100,100);
//		h_Eres2d_nowi[i] = new TH2D("h_Eres2d_nowi"+s_il, "", 20,130,330,20,-100,100);
//		h_Ediff2d[i] = new TH2D("h_Ediff2d"+s_il, "", 20,130,330,20,-100,100);
	}

	TFile *rfile = new TFile("OutTree_"+fname+".root");
	TTree *t_step =(TTree*)rfile->Get("SimStep");
	TTree *t_bar = (TTree*)rfile->Get("SimBarHit");
	TTree *t_sim = (TTree*)rfile->Get("TruthSimHit");
	TTree *t_rec = (TTree*)rfile->Get("RecHit");

	t_step->SetBranchAddress("step_x", &posx_step);
	t_step->SetBranchAddress("step_y", &posy_step);
	t_step->SetBranchAddress("step_z", &posz_step);
	t_step->SetBranchAddress("step_E", &E_step);

	t_sim->SetBranchAddress("simTruth_x", &posx_sim);	
	t_sim->SetBranchAddress("simTruth_y", &posy_sim);	
	t_sim->SetBranchAddress("simTruth_z", &posz_sim);	
	t_sim->SetBranchAddress("simTruth_E", &E_sim);	
	t_sim->SetBranchAddress("simTruth_dlayer", &dlayer_sim);	
	t_sim->SetBranchAddress("simTruth_slayer", &slayer_sim);	
	t_rec->SetBranchAddress("Rec_x", &posx_rec);
	t_rec->SetBranchAddress("Rec_y", &posy_rec);
	t_rec->SetBranchAddress("Rec_z", &posz_rec);
	t_rec->SetBranchAddress("Rec_E", &E_rec);

	for(int ii=0;ii<1;ii++){
		t_step->GetEntry(ii);
		int Nstep = E_step->size();
		float minX = 1800;
		for(int i=0;i<Nstep;i++){
	      for(int il=0;il<Nlayers;il++){
   	   if(posx_step->at(i)>1800+il*20 && posx_step->at(i)<1800+(il+1)*20 ){
      	   h_2dmap_step[il]->Fill(posz_step->at(i), posy_step->at(i), E_step->at(i)); break;
	      }}
			h_2dxy_step->Fill(posy_step->at(i), posx_step->at(i)-minX, E_step->at(i));
			h_2dxz_step->Fill(posz_step->at(i), posx_step->at(i)-minX, E_step->at(i));
			h_2dyz_step->Fill(posz_step->at(i), posy_step->at(i), E_step->at(i));
		}
	}

	for(int ii=0;ii<1;ii++){
		t_sim->GetEntry(ii);
		int Nhit = E_sim->size();
//		float minX = *min_element(posx_sim->begin(),posx_sim->end());
		float minX = 1800;
		for(int i=0;i<Nhit;i++){
		for(int il=0;il<Nlayers;il++){
		if(dlayer_sim->at(i)-1==il){
			if(slayer_sim->at(i)==0) h_2dmap_sim_s0[il]->Fill(posz_sim->at(i), posy_sim->at(i), E_sim->at(i));
			if(slayer_sim->at(i)==1) h_2dmap_sim_s1[il]->Fill(posz_sim->at(i), posy_sim->at(i), E_sim->at(i));
			h_posy_simlayer[il]->Fill(posy_sim->at(i), E_sim->at(i)); 
			h_posz_simlayer[il]->Fill(posz_sim->at(i), E_sim->at(i));
		}}
		h_posx_sim->Fill(posx_sim->at(i)-minX, E_sim->at(i));
		h_posy_sim->Fill(posy_sim->at(i), E_sim->at(i));
		h_posz_sim->Fill(posz_sim->at(i), E_sim->at(i));
		h_SimHit->Fill(posz_sim->at(i), posy_sim->at(i), E_sim->at(i));
		}
	}

   for(int ii=0;ii<1;ii++){
      t_rec->GetEntry(ii);
      int Nhit = E_rec->size();
//      float minX = *min_element(posx_rec->begin(),posx_rec->end());
		float minX = 1800;
      for(int i=0;i<Nhit;i++){
      for(int il=0;il<Nlayers;il++){
      if(posx_rec->at(i)>1800+il*20 && posx_rec->at(i)<1800+(il+1)*20 ){
			h_posy_reclayer[il]->Fill(posy_rec->at(i), E_rec->at(i)); 
			h_posz_reclayer[il]->Fill(posz_rec->at(i), E_rec->at(i));
         h_2dmap_rec[il]->Fill(posz_rec->at(i), posy_rec->at(i), E_rec->at(i)); break;
      }}
      h_posx_rec->Fill(posx_rec->at(i)-minX, E_rec->at(i));
      h_posy_rec->Fill(posy_rec->at(i), E_rec->at(i));
      h_posz_rec->Fill(posz_rec->at(i), E_rec->at(i));
		h_RecHit->Fill(posz_rec->at(i), posy_rec->at(i), E_rec->at(i));
      }
   }

	for(int i=0;i<Nlayers;i++){ cout<<h_2dmap_sim_s0[i]->Integral()+h_2dmap_sim_s1[i]->Integral()<<'\t'<<h_2dmap_rec[i]->Integral()<<endl;}

/*	int Nbinx = h_RecHit->GetNbinsX();
	int Nbiny = h_RecHit->GetNbinsY();
	for(int il=0;il<Nlayers;il++){
	for(int i=0;i<Nbinx;i++){
	for(int j=0;j<Nbiny;j++){
		double Etruth = h_2dmap_sim_s0[il]->GetBinContent(i+1,j+1) + h_2dmap_sim_s1[il]->GetBinContent(i+1,j+1);
		double Erec = h_2dmap_rec[il]->GetBinContent(i+1,j+1);
		double Eres;
		if(Etruth<1e-5){ Etruth=0; Erec=0; Eres=0;}
		else Eres = (Etruth - Erec) / Etruth;

		double weight;
		if(h_2dmap_rec[il]->Integral()>1e-7) weight = Erec/h_2dmap_rec[il]->Integral();
		else weight=0;

		h_Eres1d->Fill(Eres, weight);
		h_Eres2d[il]->SetBinContent(i+1, j+1, Eres*weight);
		h_Eres2d_nowi[il]->SetBinContent(i+1, j+1, Eres);
		h_Ediff2d[il]->SetBinContent(i+1, j+1, Etruth - Erec);
		//cout<<Eres<<'\t'<<weight<<endl;
	}}}
*/
	rfile->Close();
	wfile->Write();
	wfile->Close();

}
