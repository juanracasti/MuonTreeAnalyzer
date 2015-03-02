///////////////////////////////////////////////////////////////////////
//
//    FILE: RunProof.C
// AUTHORS: I. Gonzalez Caballero, A.Y. Rodriguez Marrero
//    DATE: 2010
//
// CONTENT: Main macro to run over MiniTrees or TESCO Trees using PROOF
//          in PROOF-Lite, PROOF-Cluster or Sequential mode
///////////////////////////////////////////////////////////////////////



TProof* proof = 0;

/////////////////////////////////////////////////////////////
/////////////           TOP ANALYSIS            /////////////
/////////////////////////////////////////////////////////////

void RunPROOF_muonAnalyzer(const char* data) {
 
  // This loads all the PROOF Analysis Framework utilities
  gROOT->LoadMacro("$PAFPATH/PAF.C");
  
  double luminosity = 20000;
  TString dataInfo;
  TString mcInfo;

 
  /////////////////////////////////////////////////////////////////////////
  //
  // PROOF SETTINGS
  // ==============
  //
  // (You may inspect scripts/PAFOptions.h to see all the posible settings)
  //
  // Edit the lines below to select tree type, input data sample, output
  // file name ...
  //
  
  ///////////////////////////////
  // PROOF MODE
  // Defines the mode in which you want to run PROOF. Read the documentation
  // for details:
  // * kSequential: No PROOF. Plain sequential code
  // * kLite: PROOF Lite mode
  // * kCluster: PROOF Cluter mode
  // * kPoD: PROOF on Demand mode
  gPAFOptions->proofMode = kCluster;
  //
  // Optional parameters for PROOF Cluster (kCluster):
  //   + The number of slots you would like to use (default is 10)
  //gPAFOptions->NSlots = 10; 
  //   + Proof server and port (default are proof.ifca.es and 1093) 
   gPAFOptions->proofServer = "proof.ifca.es";
   gPAFOptions->proofServerPort = 1093;
  //   + Maximum number of slaves per node (default is 9999, i.e. all)
  // gPAFOptions->maxSlavesPerNode = 9999;
  //
  // Start PROOF
  //
  cout << ">> Starting PROOF..." << endl;
  proof = InitProof(); 
  if (!proof && gPAFOptions->proofMode != kSequential) {
    cerr << "ERROR: I could not initialise a PROOF session!" << endl;
    return;
  }
 

  ///////////////////////////////
  // TREE TYPE
  // Defines the data formats that may be used: MiniTrees (default) or TESCO
  gPAFOptions->SetTreeType(kMiniTrees);
  //gPAFOptions->treeType = kTESCO;


  gPAFOptions->SetTreeDir("demo");
  gPAFOptions->SetTreeName("Tree");   

  ///////////////////////////////
  // INPUT DATA SAMPLE
  //
   //TString dataPath="/gpfs/csic_projects/cms/arodrig/"; //IFCA   (gridui)
  // 1) Set the base path to files
  // TString dataPath="/gpfs/csic_projects/cms/calderon/TreesCSA14"; //IFCA   (gridui)
  TString dataPath="/gpfs/csic_projects/cms/jfernan/TreesPHYS14";
  // TString dataPath="/gpfs/csic_projects/tier3data"; //IFCA   (gridui)
   //        TString dataPath="/hadoop";                      //UniOvi (fanaeui)
  //   TString dataPath="/pool/data1/MiniTrees";        //CERN   (cmsovd02)
  
  // 2) Asign the files to gPAFOptions->dataFiles (a vector<TString>)
  // Ex. MiniTree...
   //  gPAFOptions->dataFiles.push_back(dataPath + "/Data/Data7TeVRun2010A/Tree_Mu_Run2010A_Sep27ReReco_Skim2LPt1010.root");
  //
  // Note: You may consider using DatasetManager for this. See the
  //       documentation in the wiki

  TString Signal = data;

  // *********** DATA

  if (Signal=="DoubleMu_FR") {
   
    //----------  RUN2012 A ---------- // 
    //   gPAFOptions->dataFiles.push_back(dataPath + "/LatinosSkims/MuData8TeVRun2012A/Tree_DoubleMu_558.2_FR.root");  

    bool isdata  = true;
    int nEventsInTheSample = 1; //35537550; 
    double xSection = -999;

    dataInfo = "RunA2012";
  }
  
  else if (Signal=="SingleMu_FR") {

    gPAFOptions->dataFiles.push_back(dataPath + "/LatinosSkims/MuData8TeVRun2012A/Tree_SingleMu_558.2_FR.root");
 
    gPAFOptions->dataFiles.push_back(dataPath + "/LatinosSkims/MuData8TeVRun2012A/Tree_SingleMu_698.4_FR.root");
    gPAFOptions->dataFiles.push_back(dataPath + "/LatinosSkims/MuData8TeVRun2012B/Tree_SingleMu_209.4_FR.root");

    bool isdata  = true;
    int nEventsInTheSample = 1; //35537550; 
    double xSection = -999;

    dataInfo = "RunA2012";
 
  }

  //+++++++++++++++++++++++++++++++++
  // PHYS14 PU20bx25 SAMPLES

  else if (Signal=="MC_GGHWW_PU20bx25") {

    //gPAFOptions->dataFiles.push_back( "/gpfs/csic_projects/cms/calderon/TreesCSA14/Tree_GluGluToHToWWTo2LAndTau2Nu_M-125_13TeV_PU20bx25.root");
    gPAFOptions->dataFiles.push_back("/gpfs/csic_projects/cms/piedra/top/PU20bx25_PHYS14/Skim_2L_Pt17_8/Tree_GluGluToHToWW.root");
    
    bool isdata = false;
    int nEventsInTheSample = 99555; 
    //int nEventsInTheSample = 37497; 
    double xSection =  1.0 ;
    int whichRun = 2; 
 
  }

  else if (Signal=="MC_Wjets_PU20bx25") {

    //gPAFOptions->dataFiles.push_back(dataPath + "/PU20bx25/Tree_WJets_Madgraph.root");
    gPAFOptions->dataFiles.push_back("/gpfs/csic_projects/cms/piedra/top/PU20bx25_PHYS14/Skim_2L_Pt17_8/Tree_WJetsToLNu.root");
      
    bool isdata = false;
    int nEventsInTheSample = 10017462; 
    //int nEventsInTheSample = 428179; 
    double xSection =  20508.9;
    int whichRun = 2;

 }

  else if (Signal=="MC_DY_PU20bx25") {

    gPAFOptions->dataFiles.push_back("/gpfs/csic_projects/cms/piedra/top/PU20bx25_PHYS14/Skim_2L_Pt17_8/Tree_DYJetsToLL.root");
   
    bool isdata = false;
    int nEventsInTheSample = 2829164; 
    //int nEventsInTheSample = 956666; 
    double xSection =  6025.2;
    int whichRun = 2;

 }

  //++++++++++++++++++++++++++++++++++++++++



 
 else if (Signal=="MC_GGHWW_S14") {

    //gPAFOptions->dataFiles.push_back(dataPath + "/Tree_GluGluToHToWWTo2LAndTau2Nu_M-125_13TeV.root");

    gPAFOptions->dataFiles.push_back(dataPath + "/PUS14/Tree_GluGluHWW2L2Nu_NoSkim_PUS14.root");
    

    bool isdata = false;
    int nEventsInTheSample = 99555; 
    double xSection =  1.0 ;
    int whichRun = 2; 
 
  }

 else if (Signal=="MC_Wjets_S14") {

    gPAFOptions->dataFiles.push_back(dataPath + "/PUS10/Tree_WJetsToLNu_13TeV-madgraph-pythia8-tauola_NoSkim_PU_S10_0.root");
    gPAFOptions->dataFiles.push_back(dataPath + "/PUS10/Tree_WJetsToLNu_13TeV-madgraph-pythia8-tauola_NoSkim_PU_S10_1.root");
    gPAFOptions->dataFiles.push_back(dataPath + "/PUS10/Tree_WJetsToLNu_13TeV-madgraph-pythia8-tauola_NoSkim_PU_S10_2.root");
    gPAFOptions->dataFiles.push_back(dataPath + "/PUS10/Tree_WJetsToLNu_13TeV-madgraph-pythia8-tauola_NoSkim_PU_S10_3.root");
    gPAFOptions->dataFiles.push_back(dataPath + "/PUS10/Tree_WJetsToLNu_13TeV-madgraph-pythia8-tauola_NoSkim_PU_S10_4.root");
    gPAFOptions->dataFiles.push_back(dataPath + "/PUS10/Tree_WJetsToLNu_13TeV-madgraph-pythia8-tauola_NoSkim_PU_S10_5.root");
    gPAFOptions->dataFiles.push_back(dataPath + "/PUS10/Tree_WJetsToLNu_13TeV-madgraph-pythia8-tauola_NoSkim_PU_S10_6.root");
    gPAFOptions->dataFiles.push_back(dataPath + "/PUS10/Tree_WJetsToLNu_13TeV-madgraph-pythia8-tauola_NoSkim_PU_S10_7.root");
    gPAFOptions->dataFiles.push_back(dataPath + "/PUS10/Tree_WJetsToLNu_13TeV-madgraph-pythia8-tauola_NoSkim_PU_S10_8.root");
    gPAFOptions->dataFiles.push_back(dataPath + "/PUS10/Tree_WJetsToLNu_13TeV-madgraph-pythia8-tauola_NoSkim_PU_S10_9.root");
    gPAFOptions->dataFiles.push_back(dataPath + "/PUS10/Tree_WJetsToLNu_13TeV-madgraph-pythia8-tauola_NoSkim_PU_S10_10.root");
    gPAFOptions->dataFiles.push_back(dataPath + "/PUS10/Tree_WJetsToLNu_13TeV-madgraph-pythia8-tauola_NoSkim_PU_S10_11.root");
    gPAFOptions->dataFiles.push_back(dataPath + "/PUS10/Tree_WJetsToLNu_13TeV-madgraph-pythia8-tauola_NoSkim_PU_S10_12.root");
    gPAFOptions->dataFiles.push_back(dataPath + "/PUS10/Tree_WJetsToLNu_13TeV-madgraph-pythia8-tauola_NoSkim_PU_S10_13.root");
    gPAFOptions->dataFiles.push_back(dataPath + "/PUS10/Tree_WJetsToLNu_13TeV-madgraph-pythia8-tauola_NoSkim_PU_S10_14.root");
    gPAFOptions->dataFiles.push_back(dataPath + "/PUS10/Tree_WJetsToLNu_13TeV-madgraph-pythia8-tauola_NoSkim_PU_S10_15.root");
    gPAFOptions->dataFiles.push_back(dataPath + "/PUS10/Tree_WJetsToLNu_13TeV-madgraph-pythia8-tauola_NoSkim_PU_S10_16.root");
    gPAFOptions->dataFiles.push_back(dataPath + "/PUS10/Tree_WJetsToLNu_13TeV-madgraph-pythia8-tauola_NoSkim_PU_S10_17.root");
    gPAFOptions->dataFiles.push_back(dataPath + "/PUS10/Tree_WJetsToLNu_13TeV-madgraph-pythia8-tauola_NoSkim_PU_S10_18.root");
    gPAFOptions->dataFiles.push_back(dataPath + "/PUS10/Tree_WJetsToLNu_13TeV-madgraph-pythia8-tauola_NoSkim_PU_S10_19.root");
    gPAFOptions->dataFiles.push_back(dataPath + "/PUS10/Tree_WJetsToLNu_13TeV-madgraph-pythia8-tauola_NoSkim_PU_S10_20.root");
    gPAFOptions->dataFiles.push_back(dataPath + "/PUS10/Tree_WJetsToLNu_13TeV-madgraph-pythia8-tauola_NoSkim_PU_S10_21.root");
    gPAFOptions->dataFiles.push_back(dataPath + "/PUS10/Tree_WJetsToLNu_13TeV-madgraph-pythia8-tauola_NoSkim_PU_S10_22.root");
       

    bool isdata = false;
    int nEventsInTheSample = 41122326; 
    double xSection = 61526.7  ;
    int whichRun = 1; 
  }




  else if (Signal=="MC_Wjets_8TeV") {

    gPAFOptions->dataFiles.push_back("/gpfs/csic_projects/tier3data/MC_Summer12_53X/SUSY/Tree_WJets_Madgraph_0.root");
    gPAFOptions->dataFiles.push_back("/gpfs/csic_projects/tier3data/MC_Summer12_53X/SUSY/Tree_WJets_Madgraph_1.root");
    gPAFOptions->dataFiles.push_back("/gpfs/csic_projects/tier3data/MC_Summer12_53X/SUSY/Tree_WJets_Madgraph_2.root");

    bool isdata = false;
    int nEventsInTheSample = 76102995; 
    double xSection = 37509  ;
    int whichRun = 1; 

 }

  else if (Signal=="MC_GGHWW_8TeV") {
    gPAFOptions->dataFiles.push_back("/gpfs/csic_projects/tier3data/MC_Summer12_53X/SUSY/Tree_HWW125.root");
  
    bool isdata = false;
    int nEventsInTheSample = 299975; 
    double xSection = 0.444 ;
    int whichRun = 1; 
}

  else if (Signal=="MC_DY_8TeV") {
    gPAFOptions->dataFiles.push_back("/gpfs/csic_projects/tier3data/MC_Summer12_53X/SUSY/Tree_ZJets_Madgraph_0.root");
    gPAFOptions->dataFiles.push_back("/gpfs/csic_projects/tier3data/MC_Summer12_53X/SUSY/Tree_ZJets_Madgraph_1.root");
    gPAFOptions->dataFiles.push_back("/gpfs/csic_projects/tier3data/MC_Summer12_53X/SUSY/Tree_ZJets_Madgraph_2.root");
    gPAFOptions->dataFiles.push_back("/gpfs/csic_projects/tier3data/MC_Summer12_53X/SUSY/Tree_ZJets_Madgraph_3.root");
    gPAFOptions->dataFiles.push_back("/gpfs/csic_projects/tier3data/MC_Summer12_53X/SUSY/Tree_ZJets_Madgraph_4.root");
    gPAFOptions->dataFiles.push_back("/gpfs/csic_projects/tier3data/MC_Summer12_53X/SUSY/Tree_ZJets_Madgraph_5.root");
    gPAFOptions->dataFiles.push_back("/gpfs/csic_projects/tier3data/MC_Summer12_53X/SUSY/Tree_ZJets_Madgraph_6.root");
    gPAFOptions->dataFiles.push_back("/gpfs/csic_projects/tier3data/MC_Summer12_53X/SUSY/Tree_ZJets_Madgraph_7.root");
    gPAFOptions->dataFiles.push_back("/gpfs/csic_projects/tier3data/MC_Summer12_53X/SUSY/Tree_ZJets_Madgraph_8.root");
    gPAFOptions->dataFiles.push_back("/gpfs/csic_projects/tier3data/MC_Summer12_53X/SUSY/Tree_ZJets_Madgraph_9.root");
    
  
    bool isdata = false;
    int nEventsInTheSample = 10000000; 
    double xSection = 40000 ;
    int whichRun = 1; 
}

  // *********** DatasetManager
  else {


    gROOT->LoadMacro("/gpfs/csic_users/calderon/UserCode/IGonzalez/DatasetManager/DatasetManager.C+");

    DatasetManager* dm = new DatasetManager("Spring11Latinos");

    dm->RedownloadFiles();

    dm->LoadDataset(Signal);  // Load information about a given dataset

    G_Event_Weight = dm->GetCrossSection() * G_Event_Lumi / dm->GetEventsInTheSample();

    cout << endl;
    cout << "      x-section = " << dm->GetCrossSection()      << endl;
    cout << "     luminosity = " << G_Event_Lumi               << endl;
    cout << "        nevents = " << dm->GetEventsInTheSample() << endl;
    cout << " base file name = " << dm->GetBaseFileName()      << endl;


    int nEventsInTheSample = dm->GetEventsInTheSample();
    double xSection = dm->GetCrossSection();
    bool isdata = false;

    gPAFOptions->dataFiles = dm->GetFiles();

  }

  

  ///////////////////////////////
  // OUTPUT FILE NAME
  // Specify the name of the file where you want your histograms to be saved
  

  std::ostringstream out;
 
  TString outTest = out.str();

  //TString output = TString("csa14_GluGluToHToWWTo2LAndTau2Nu_M-125_13TeV_PUS14_Vertex.root"); 
  //TString output = TString("csa14_W1234JetsToLNu_Tune4C_13TeV_PUS14.root");
  //TString output = TString("csa14_WToMuNu_Tune4C_13TeV_PUS14.root");
  //TString output = TString("csa14_GluGluToHToWWTo2LAndTau2Nu_M-125_13TeV_PU20bx25.root");
  //TString output = TString("csa14_WJetsToLNu_13TeV-madgraph_PU20bx25.root"); 
  //TString output = TString("csa14_WJets_Madgraph_8Tev.root");
  //TString output = TString("csa14_HWW125_8Tev.root");


    TString output = TString("phys14_"+Signal+".root"); 


  //TString output = TString(Signal+".root"); 

  gPAFOptions->outputFile=output;
  


  ///////////////////////////////
  // PARAMETERS FOR THE ANALYSIS
  // This parameters are passed to the analysis class and can be use there.
  // They are stored in a InputParameters object. They are saved to the 
  // output file.
  // See packages/InputParameters/InputParameters.h for information on how
  // to use this class.

 
  gPAFOptions->inputParameters = new InputParameters();

  gPAFOptions->inputParameters->SetNamedBool("IsDATA", isdata);
  gPAFOptions->inputParameters->SetNamedString("Signal", data);
  gPAFOptions->inputParameters->SetNamedDouble("XSection", xSection);
  gPAFOptions->inputParameters->SetNamedDouble("Luminosity", luminosity);
  gPAFOptions->inputParameters->SetNamedInt("NEvents", nEventsInTheSample); // all
  gPAFOptions->inputParameters->SetNamedFloat("luminosityPU", 19468.3);  
  gPAFOptions->inputParameters->SetNamedInt("WhichRun", whichRun);
   

  ////// I.G.
  //Find the total number of entries in the dataset and send it to the input parameters
  /*  TChain* chain = new TChain("Tree", "Tree");
      for (unsigned int i = 0; i < dataFiles.size(); i++) 
      chain->Add(dataFiles[i]);
      gPAFOptions->inputParameters->SetNamedInt("NEventsTotal", chain->GetEntries()); //after skimming
      TString eventsfile(gSystem->pwd());
      eventsfile+="/";
      eventsfile+=Signal;
      eventsfile+="_events.log";
      gPAFOptions->inputParameters->SetNamedString("fFileList", (const char*) eventsfile);
  */

  ///////////////////////////////
  // DYNAMIC HISTOGRAMS
  // Specify the name of the histograms you would like to monitor as they are
  // filled by PROOF
  //
  //  gPAFOptions->dynamicHistograms.push_back("myHistogram");
  //...
  
  /////////////////////////////
  // NUMBER OF EVENTS 
  // Specify the number (Long64_t) of events to process.
  // Set it to -1 to use the full sample.
 
  gPAFOptions->SetNEvents(-1);
  
  //
  /////////////////////////////////////////////////////////////////////////







  /////////////////////////////////////////////////////////////////////////
  //
  // EXTRA proof settings:
  // ====================
  //
  // It is unlikely that you need to edit anything below. At least at the
  // beginning of your PAF experience. However we provide a couple of hooks
  // for extensions.
  //

  ///////////////////////////////
  // NAME OF ANALYSIS CLASS. 
  // If 0 the default name schema will be used, i.e. depending on the value
  // of gPAFOptions->treeType: MyAnalysisTESCO or MyAnalsyisMiniTrees
  //
  
  gPAFOptions->SetAnalysis("muonAnalyzer");


  ///////////////////////////////
  // ADDITIONAL PACKAGES TO BE UPLOADED TO PROOF.
  // The mandatory ones are added automatically in PAFOptions
  //


  gPAFOptions->AddPackage("PUWeight");
  //  gPAFOptions->AddPackage("MuonIsoMVA");


  ///////////////////////////////
  // CONTROL OUTPUT AND CHECKS
  // + If true (default) PAF checks for new version in CVS every time
  // gPAFOptions->checkVersion = true;
  // + If true (default) the output file is reopened so the objects in the
  //   file can be interactively accessed. The object in the output are also
  //   listed
  // gPAFOptions->reopenOutputFile = false;

  //
  /////////////////////////////////////////////////////////////////////////



 
  /////////////////////////////////////////////////////////////////////////
  //
  // RUN THE ANALYSIS
  // ================
  //
  // If something needs to be edited below (or inside), contact the 
  // developers.
  //
  // Run the analysis
  //
  //gPAFOptions->reopenOutputFile = false;
  gPAFOptions->reopenOutputFileRemoved= false;

  if (!RunAnalysis())
    cerr << "ERROR: There was a problem running the analysis!" << endl;
  //
  /////////////////////////////////////////////////////////////////////////


}
