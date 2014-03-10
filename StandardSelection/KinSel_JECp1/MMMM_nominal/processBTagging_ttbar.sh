root -l -b -q '../../HbbHbb_bTagging.cc++("Signal", "TTJets_FullLeptMGDecays_8TeV-madgraph_Skim", "MMMM", "nominal", "JECp1")' 
say "b-tagging of fully-leptonic t t-bar Monte Carlo done."
root -l -b -q '../../HbbHbb_bTagging.cc++("Signal", "TTJets_SemiLeptMGDecays_8TeV-madgraph_Skim", "MMMM", "nominal", "JECp1")' 
say "b-tagging of semi-leptonic t t-bar Monte Carlo done."
root -l -b -q '../../HbbHbb_bTagging.cc++("Signal", "TTJets_HadronicMGDecays_8TeV-madgraph_Skim", "MMMM", "nominal", "JECp1")'  
say "b-tagging of hadronic t t-bar Monte Carlo done."
