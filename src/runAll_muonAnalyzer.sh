#////////////////////////////////////////////////////////////////
#////////////////////////////////////////////////////////////////

#/////            Muon Studies for CSA14                     ////

#/////        A. Calderón (IFCA)   18 / 08 / 2014            ////

#////////////////////////////////////////////////////////////////
#////////////////////////////////////////////////////////////////


#root -l -b -q 'RunPROOF_muonAnalyzer.C("MC_GGHWW_PU20bx25")';
#root -l -b -q 'RunPROOF_muonAnalyzer.C("MC_Wjets_PU20bx25")';
root -l -b -q 'RunPROOF_muonAnalyzer.C("MC_DY_PU20bx25")';
root -l -b -q 'RunPROOF_muonAnalyzer.C("MC_TTbar_PU20bx25")';
#root -l -b -q 'RunPROOF_muonAnalyzer.C("DY_ISO_PU20bx25")';
#root -l -b -q 'RunPROOF_muonAnalyzer.C("TTbar_ISO_PU20bx25")';
#root -l -b -q 'RunPROOF_muonAnalyzer.C("QCD_ISO_PU20bx25")';
#root -l -b -q 'RunPROOF_muonAnalyzer.C("SingleMu_720")';
#root -l -b -q 'RunPROOF_muonAnalyzer.C("SingleMu_740")';


