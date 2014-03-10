root -l -b -q 'HbbHbb_KinematicSelection.cc++("ttbar", "TTJets_FullLeptMGDecays_8TeV-madgraph_Skim")'
say "Kinematic selection of fully-leptonic t t-bar Monte Carlo done."
root -l -b -q 'HbbHbb_KinematicSelection.cc++("ttbar", "TTJets_SemiLeptMGDecays_8TeV-madgraph_Skim")'
say "Kinematic selection of semi-leptonic t t-bar Monte Carlo done."
root -l -b -q 'HbbHbb_KinematicSelection.cc++("ttbar", "TTJets_HadronicMGDecays_8TeV-madgraph_Skim")'
say "Kinematic selection of hadronic t t-bar Monte Carlo done."
