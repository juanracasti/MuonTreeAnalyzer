////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

/////            Muon Studies for CSA14                     ////

/////        A. Calderón (IFCA)   18 / 08 / 2014            ////

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////





#ifndef MUONFR_NEWJETETDEF_H
#define MUONFR_NEWJETETDEF_H 1

#include "PAFAnalysis.h"

#include <TH1F.h>
#include <TMatrix.h>
#include <TH2F.h>
#include "TCounterUI.h"
#include <TLorentzVector.h>
#include "Riostream.h"  
#include "PUWeight.h"


class muonAnalyzer: public PAFAnalysis{

 public:
   muonAnalyzer(TTree *tree=0);
   virtual ~muonAnalyzer() {}

 protected:
   virtual void              Initialise();
   virtual void              InsideLoop();
   virtual void              SetDataMembersAtTermination();
   virtual void              Summary(); 
   void                      GetAllMuons();
   bool                      passTightMuCuts(int, float, float, float, int, int);
   bool                      passHWWMuCuts(int);
   void                      SetGenInfo(TString);
   bool                      MatchGenToReco(int&, int&, TString);
   void                      doEffsGEN(float, TString);
   void                      doEffsRECO(int, float);
   void                      doEffsPtGEN(float, float, float, TString);
   void                      doEffsEtaGEN(float, float, float, TString);
   void                      doEffsNpvGEN(float, float, float, TString);
   void                      doEffsPtRECO(int, float, float, float);
   void                      doEffsEtaRECO(int, float, float, float);
   void                      doEffsNpvRECO(int, float, float, float);
   void                      FillRelEff(string, int, int, double, int);
   void                      Filldz2D(int, int, double, int);
   void                      FillPFIso(int, int, double, string);
   void                      FillTypeMu(double, string);
   void                      FillPtEta(int, int, double, string);
   bool                      passPFIso (int, string);
   double                    getPFRelIso (int, string);	
   void                      GetAllJets(double);
   Int_t                     SelectedVertexIndex();
   Int_t                     SelectedVertexIndex(int);
   Int_t                     SelectedVertexIndex(TString, int);
   Double_t                  get_dz(int, int);
   Double_t                  giveMT(int, double, double);
   Double_t                  getDeltaR (int , int );
   float                     DeltaPhi(float, float);
   float                     projectedMET(int, int);
   float                     deltaPhiJet(int, int);
   float                     ptDilepton(int, int);
   bool                      passesAntiBTagging(int, int);    
   bool                      passesWWSelection(int, int, int);

  
 public:
 
   // My Declarations:
   // Define global variables

   PUWeight* fPUWeight;
 
   ///** GEN INFORMATION
   std::vector<TLorentzVector> G_GEN_PromptMuon_4vec;
   std::vector<TLorentzVector> G_GEN_Muon_4vec;

   //** MUONS
   std::vector<TLorentzVector> G_Muon_4vec;
   std::vector<bool>           G_Muon_GLBPFID;
   std::vector<bool>           G_Muon_HWW_ID;
   std::vector<bool>           G_Muon_dzID;
   std::vector<bool>           G_Muon_fromPVID;
   std::vector<bool>           G_Muon_TightIDbutdz;
   std::vector<bool>           G_Muon_TightIDfromPV;
   std::vector<bool>           G_Muon_TightID;

   std::vector<bool>           G_Muon_ISOR03;
   std::vector<bool>           G_Muon_ISOR04;
   std::vector<bool>           G_Muon_ISOdBetaR03;
   std::vector<bool>           G_Muon_ISOdBetaR04;
   std::vector<bool>           G_Muon_ISOPFWeightsR03;
   std::vector<bool>           G_Muon_ISOPFWeightsR04;
   
   std::vector<bool>           G_Muon_HWW_ISOR03;
   std::vector<bool>           G_Muon_HWW_ISOR04;
   std::vector<bool>           G_Muon_HWW_ISOdBetaR03;
   std::vector<bool>           G_Muon_HWW_ISOdBetaR04;
   std::vector<bool>           G_Muon_HWW_ISOPFWeightsR03;
   std::vector<bool>           G_Muon_HWW_ISOPFWeightsR04;
   
   std::vector<bool>           G_Muon_TightISOR03;
   std::vector<bool>           G_Muon_TightISOR04;
   std::vector<bool>           G_Muon_TightISOdBetaR03;
   std::vector<bool>           G_Muon_TightISOdBetaR04;
   std::vector<bool>           G_Muon_TightISOPFWeightsR03;
   std::vector<bool>           G_Muon_TightISOPFWeightsR04;
   
   //std::vector<double>         G_Muon_PFIsoBeta;
   //std::vector<bool>           G_Muon_MatchW;

   std::vector<bool>           G_Muon_Match_LeadJet;
   std::vector<bool>           G_Muon_Numerator;
   std::vector<bool>           G_Muon_Denominator;
   std::vector<double>         G_Muon_RelCombIsol;
   std::vector<double>         G_Muon_JetEt;
   std::vector<double>         G_Muon_JetEta;
   std::vector<bool>           G_Muon_match_jet;
   std::vector<double>         G_Muon_isoMVA;
  
   //std::vector<double>         G_Muon_ChargeHPFRelISo;
   //std::vector<double>         G_Muon_PFRelISo;
   //std::vector<double>         G_Muon_PFChargedH;
   //std::vector<double>         G_Muon_PFNeutralH;
   //std::vector<double>         G_Muon_PFPhoton;
   //std::vector<double>         G_Muon_PFRho;

   bool        G_isMuMu;
   bool        G_isMuTau;
   bool        G_isTauMu;
   bool        G_isTauTau;
   bool        G_isNonPrompt;


   //** JETS 
   std::vector<TLorentzVector> G_Jet_4vec;   
   int         G_N_Jets;

   //** VERTEX
   int         G_PV_Index;
   


   bool G_Debug_DefineAnalysisVariables;

 
   //* Histograms 

   TH1F       *h_N_PV0_PVLep;
   TH1F       *h_N_dZ_PV0_PVLep;

   TH1F       *h_N_PV;
   TH1F       *h_N_PV2;
   TH1F       *h_N_PV3;
   TH1F       *h_N_Dilep_TypeMu;
   TH1F       *h_N_Dilep_TypeMu_LP;
   TH1F       *h_N_Dilep_TypeMu_HP;
   TH1F       *h_N_WWlevel_TypeMu;
   TH1F       *h_N_WWlevel_TypeMu_LP; 
   TH1F       *h_N_WWlevel_TypeMu_HP;
   TH1F       *h_N_Dilep_TightMuCuts;
   TH1F       *h_N_WWlevel_TightMuCuts;

   // Efficiencies vs pt, eta, and npv

   TH1F       *h_Dilep_Eff_pt_AllRECO[2];
   TH1F       *h_Dilep_Eff_pt_GEN[2];
   TH1F       *h_Dilep_Eff_pt_AllMatched[2];
   TH1F       *h_Dilep_Eff_pt_GLBID[2];
   TH1F       *h_Dilep_Eff_pt_PFID[2];
   TH1F       *h_Dilep_Eff_pt_GLBPFID[2];
   TH1F       *h_Dilep_Eff_pt_dzID[2];
   TH1F       *h_Dilep_Eff_pt_fromPVID[2];
   TH1F       *h_Dilep_Eff_pt_TightIDbutdz[2];
   TH1F       *h_Dilep_Eff_pt_TightIDfromPV[2];
   TH1F       *h_Dilep_Eff_pt_TightIDAndfromPV[2];
   TH1F       *h_Dilep_Eff_pt_TightID[2];
   TH1F       *h_Dilep_Eff_pt_HWWID[2];
   TH1F       *h_Dilep_Eff_pt_TightISOR03[2];
   TH1F       *h_Dilep_Eff_pt_TightISOR04[2];
   TH1F       *h_Dilep_Eff_pt_TightISOdBetaR03[2];
   TH1F       *h_Dilep_Eff_pt_TightISOdBetaR04[2];
   TH1F       *h_Dilep_Eff_pt_TightISOPFWeightsR03[2];
   TH1F       *h_Dilep_Eff_pt_TightISOPFWeightsR04[2];

   TH1F       *h_Dilep_Eff_eta_AllRECO[2];
   TH1F       *h_Dilep_Eff_eta_GEN[2];
   TH1F       *h_Dilep_Eff_eta_AllMatched[2];
   TH1F       *h_Dilep_Eff_eta_GLBID[2];
   TH1F       *h_Dilep_Eff_eta_PFID[2];
   TH1F       *h_Dilep_Eff_eta_GLBPFID[2];
   TH1F       *h_Dilep_Eff_eta_dzID[2];
   TH1F       *h_Dilep_Eff_eta_fromPVID[2];
   TH1F       *h_Dilep_Eff_eta_TightIDbutdz[2];
   TH1F       *h_Dilep_Eff_eta_TightIDfromPV[2];
   TH1F       *h_Dilep_Eff_eta_TightIDAndfromPV[2];
   TH1F       *h_Dilep_Eff_eta_TightID[2];
   TH1F       *h_Dilep_Eff_eta_HWWID[2];
   TH1F       *h_Dilep_Eff_eta_TightISOR03[2];
   TH1F       *h_Dilep_Eff_eta_TightISOR04[2];
   TH1F       *h_Dilep_Eff_eta_TightISOdBetaR03[2];
   TH1F       *h_Dilep_Eff_eta_TightISOdBetaR04[2];
   TH1F       *h_Dilep_Eff_eta_TightISOPFWeightsR03[2];
   TH1F       *h_Dilep_Eff_eta_TightISOPFWeightsR04[2];

   TH1F       *h_Dilep_Eff_npv_AllRECO[2];
   TH1F       *h_Dilep_Eff_npv_GEN[2];
   TH1F       *h_Dilep_Eff_npv_AllMatched[2];
   TH1F       *h_Dilep_Eff_npv_GLBID[2];
   TH1F       *h_Dilep_Eff_npv_PFID[2];
   TH1F       *h_Dilep_Eff_npv_GLBPFID[2];
   TH1F       *h_Dilep_Eff_npv_dzID[2];
   TH1F       *h_Dilep_Eff_npv_fromPVID[2];
   TH1F       *h_Dilep_Eff_npv_TightIDbutdz[2];
   TH1F       *h_Dilep_Eff_npv_TightIDfromPV[2];
   TH1F       *h_Dilep_Eff_npv_TightIDAndfromPV[2];
   TH1F       *h_Dilep_Eff_npv_TightID[2];
   TH1F       *h_Dilep_Eff_npv_HWWID[2];
   TH1F       *h_Dilep_Eff_npv_TightISOR03[2];
   TH1F       *h_Dilep_Eff_npv_TightISOR04[2];
   TH1F       *h_Dilep_Eff_npv_TightISOdBetaR03[2];
   TH1F       *h_Dilep_Eff_npv_TightISOdBetaR04[2];
   TH1F       *h_Dilep_Eff_npv_TightISOPFWeightsR03[2];
   TH1F       *h_Dilep_Eff_npv_TightISOPFWeightsR04[2];
  
   //TH1F       *h_N_Dilep_TightMuCuts_butTkLayers [2]; 
   //TH1F       *h_N_Dilep_TightMuCuts_butSTAHits [2];
   //TH1F       *h_N_WWlevel_TightMuCuts_butTkLayers [2]; 
   //TH1F       *h_N_WWlevel_TightMuCuts_butSTAHits [2];

   //TH1F       *h_N_Dilep_GLBPF_butTkLayers [2]; 
   //TH1F       *h_N_Dilep_GLBPF_butSTAHits [2];
   //TH1F       *h_N_WWlevel_GLBPF_butTkLayers [2]; 
   //TH1F       *h_N_WWlevel_GLBPF_butSTAHits [2];

   //TH1F       *h_Dilep_AllMu_PFIsoBeta_Mu1;
   //TH1F       *h_Dilep_AllMu_PFIsoBeta_Mu2;
   //TH1F       *h_Dilep_HWWMu_PFIsoBeta_Mu1;
   //TH1F       *h_Dilep_HWWMu_PFIsoBeta_Mu2; 
   //TH1F       *h_WWlevel_HWWMu_PFIsoBeta_Mu1;
   //TH1F       *h_WWlevel_HWWMu_PFIsoBeta_Mu2;


   //------>  Plots after Tight Muon ID at dilepton level
   TH1F       *h_Dilep_TightMu_TkLayers [2];
   TH1F       *h_Dilep_TightMu_StaHits [2];
   //TH1F       *h_Dilep_TightMu_PFCH [2];
   //TH1F       *h_Dilep_TightMu_PFNH [2]; 
   //TH1F       *h_Dilep_TightMu_PFPh [2]; 
   //TH1F       *h_Dilep_TightMu_PFRho [2];

   TH1F       *h_Dilep_TightMu_pt [2];
   TH1F       *h_Dilep_TightMu_eta [2];
   //   TH1F       *h_Dilep_TightMu_PFIsoBeta [2];
   TH1F       *h_Dilep_TightMu_RelEff [2];
   TH1F       *h_Dilep_Chi2 [2];
   TH1F       *h_Dilep_StaHits [2];
   TH2F       *h_Dilep_StaHitsEta [2];
   TH2F       *h_Dilep_StaHitsPV [2];
   TH1F       *h_Dilep_StaNStation [2];
   TH2F       *h_Dilep_StaNStationEta [2];
   TH2F       *h_Dilep_StaNStationPV [2];
   TH1F       *h_Dilep_PixelHits [2]; 
   TH1F       *h_Dilep_TkLayers [2]; 
   TH1F       *h_Dilep_dxy [2];
   TH1F       *h_Dilep_dz_GTrack [2];
   TH1F       *h_Dilep_dz_InTrack [2];

   TH1F       *h_Dilep_dz [2];
   TH1F       *h_Dilep_dzMu [2];
   TH1F       *h_Dilep_dz_BestTrack [2];
   TH1F       *h_Dilep_dz_bf_dy [2];
   TH1F       *h_Dilep_dzMu_bf_dy [2];
   TH1F       *h_Dilep_dz_BestTrack_bf_dy [2];
   TH1F       *h_TrueDilep_dz [2];
   TH1F       *h_TrueDilep_dzMu [2];
   TH1F       *h_TrueDilep_dz_BestTrack [2];
   TH1F       *h_TrueDilep_dz_bf_dy [2];
   TH1F       *h_TrueDilep_dzMu_bf_dy [2];
   TH1F       *h_TrueDilep_dz_BestTrack_bf_dy [2];

   TH1F       *h_Dilep_dz_PV_Eq [2];
   TH1F       *h_Dilep_dzMu_PV_Eq [2];
   TH1F       *h_Dilep_dz_PV_NEq [2];
   TH1F       *h_Dilep_dzMu_PV_NEq [2];

   TH1F       *h_2D_Dilep_dz_Eff;
   TH1F       *h_2D_Dilep_dzMu_Eff;
   TH1F       *h_2D_Dilep_dzBT_Eff;
   TH1F       *h_2D_Dilep_dz_bf_dy_Eff;
   TH1F       *h_2D_Dilep_dzMu_bf_dy_Eff;
   TH1F       *h_2D_Dilep_dzBT_bf_dy_Eff;

   TH2F       *h_2D_Dilep_dz;
   TH2F       *h_2D_Dilep_dzMu;
   TH2F       *h_2D_Dilep_dzBT;
   TH2F       *h_2D_Dilep_dz_PV_Eq;
   TH2F       *h_2D_Dilep_dzMu_PV_Eq;
   TH2F       *h_2D_Dilep_dzBT_PV_Eq;
   TH2F       *h_2D_Dilep_dz_PV_NEq;
   TH2F       *h_2D_Dilep_dzMu_PV_NEq;
   TH2F       *h_2D_Dilep_dzBT_PV_NEq;
   TH2F       *h_2D_Dilep_dz_bf_dy;
   TH2F       *h_2D_Dilep_dzMu_bf_dy;
   TH2F       *h_2D_Dilep_dzBT_bf_dy;
   TH2F       *h_2D_Dilep_dz_bf_dy_PV_Eq;
   TH2F       *h_2D_Dilep_dzMu_bf_dy_PV_Eq;
   TH2F       *h_2D_Dilep_dzBT_bf_dy_PV_Eq;
   TH2F       *h_2D_Dilep_dz_bf_dy_PV_NEq;
   TH2F       *h_2D_Dilep_dzMu_bf_dy_PV_NEq;
   TH2F       *h_2D_Dilep_dzBT_bf_dy_PV_NEq;

   TH1F       *h_Dilep_Gen_Muon_MpdgId;
   TH1F       *h_Dilep_Gen_Muon_MpdgId_PV_Eq;
   TH1F       *h_Dilep_Gen_Muon_MpdgId_PV_NEq;
   TH1F       *h_Dilep_Jet_Energy;
   TH1F       *h_Dilep_Jet_Energy_PV_Eq;
   TH1F       *h_Dilep_Jet_Energy_PV_NEq;
   TH1F       *h_Dilep_N_Jets;
   TH1F       *h_Dilep_N_Jets_PV_Eq;
   TH1F       *h_Dilep_N_Jets_PV_NEq;

   //-------- Rel. Iso plots
   TH1F       *h_Dilep_AllMu_PFRelIso_R03[2];
   TH1F       *h_Dilep_AllMu_PFRelIso_R04[2];
   TH1F       *h_Dilep_AllMu_PFRelIso_dBetaR03[2];
   TH1F       *h_Dilep_AllMu_PFRelIso_dBetaR04[2];
   TH1F       *h_Dilep_AllMu_PFRelIso_PFWeightsR03[2];
   TH1F       *h_Dilep_AllMu_PFRelIso_PFWeightsR04[2];
   TH1F       *h_Dilep_HWWMu_PFRelIso_R03[2];
   TH1F       *h_Dilep_HWWMu_PFRelIso_R04[2];
   TH1F       *h_Dilep_HWWMu_PFRelIso_dBetaR03[2];
   TH1F       *h_Dilep_HWWMu_PFRelIso_dBetaR04[2];
   TH1F       *h_Dilep_HWWMu_PFRelIso_PFWeightsR03[2];
   TH1F       *h_Dilep_HWWMu_PFRelIso_PFWeightsR04[2];
   TH1F       *h_Dilep_TightMu_PFRelIso_R03[2];
   TH1F       *h_Dilep_TightMu_PFRelIso_R04[2];
   TH1F       *h_Dilep_TightMu_PFRelIso_dBetaR03[2];
   TH1F       *h_Dilep_TightMu_PFRelIso_dBetaR04[2];
   TH1F       *h_Dilep_TightMu_PFRelIso_PFWeightsR03[2];
   TH1F       *h_Dilep_TightMu_PFRelIso_PFWeightsR04[2];
   TH1F       *h_WWlevel_AllMu_PFRelIso_R03[2];
   TH1F       *h_WWlevel_AllMu_PFRelIso_R04[2];
   TH1F       *h_WWlevel_AllMu_PFRelIso_dBetaR03[2];
   TH1F       *h_WWlevel_AllMu_PFRelIso_dBetaR04[2];
   TH1F       *h_WWlevel_AllMu_PFRelIso_PFWeightsR03[2];
   TH1F       *h_WWlevel_AllMu_PFRelIso_PFWeightsR04[2];
   TH1F       *h_WWlevel_HWWMu_PFRelIso_R03[2];
   TH1F       *h_WWlevel_HWWMu_PFRelIso_R04[2];
   TH1F       *h_WWlevel_HWWMu_PFRelIso_dBetaR03[2];
   TH1F       *h_WWlevel_HWWMu_PFRelIso_dBetaR04[2];
   TH1F       *h_WWlevel_HWWMu_PFRelIso_PFWeightsR03[2];
   TH1F       *h_WWlevel_HWWMu_PFRelIso_PFWeightsR04[2];
   TH1F       *h_WWlevel_TightMu_PFRelIso_R03[2];
   TH1F       *h_WWlevel_TightMu_PFRelIso_R04[2];
   TH1F       *h_WWlevel_TightMu_PFRelIso_dBetaR03[2];
   TH1F       *h_WWlevel_TightMu_PFRelIso_dBetaR04[2];
   TH1F       *h_WWlevel_TightMu_PFRelIso_PFWeightsR03[2];
   TH1F       *h_WWlevel_TightMu_PFRelIso_PFWeightsR04[2];



   //------>  Plots after Tight Muon ID+ISO at dilepton level
   TH1F       *h_Dilep_TightMuISO_pt  [2];
   TH1F       *h_Dilep_TightMuISO_eta [2];

   //------>  Plots after Tight Muon ID at WW level 
   TH1F       *h_WWlevel_TightMu_pt [2];
   TH1F       *h_WWlevel_TightMu_eta [2];
   TH1F       *h_WWlevel_TightMu_PFIsoBeta [2];
   TH1F       *h_WWlevel_Chi2 [2];
   TH1F       *h_WWlevel_StaHits [2];
   TH1F       *h_WWlevel_StaNStation [2];
   TH1F       *h_WWlevel_PixelHits [2]; 
   TH1F       *h_WWlevel_TkLayers [2]; 
   TH1F       *h_WWlevel_dxy [2];
   TH1F       *h_WWlevel_dz [2];

   //------>  Plots after Tight Muon ID+ISO at dilepton level
   TH1F       *h_WWlevel_TightMuISO_pt  [2];
   TH1F       *h_WWlevel_TightMuISO_eta [2];


   // * Input parameters
   TString  Signal;       // Type of Signal
   int      NEvents;      // Total number of events in the sample before skim
   double   Luminosity;   // Total luminosity
   double   XSection;     // Process cross section
   bool     IsDATA; // True if is Data, False in case MC
   int      WhichRun; // 1 in case of RunI samples. 2 In case of RunII samples.;

   ClassDef(muonAnalyzer,0);
};
#endif

