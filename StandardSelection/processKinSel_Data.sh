root -l -b -q 'HbbHbb_KinematicSelection.cc++("Data_V42_V8_Full", "BJetPlusX_Run2012B-13Jul2012-v1_EDMNtuple_V42_ProcV1_Skim")'
say "2012B kinematically selected"

root -l -b -q 'HbbHbb_KinematicSelection.cc++("Data_V42_V8_Full", "BJetPlusX_Run2012C-24Aug2012-v2_EDMNtuple_V42_ProcV1_Skim")'
say "2012C August reconstruction kinematically selected"

root -l -b -q 'HbbHbb_KinematicSelection.cc++("Data_V42_V8_Full", "BJetPlusX_Run2012C-PromptReco-v2_EDMNtuple_V42_ProcV1_Skim")'
say "2012C part 2 kinematically selected"

root -l -b -q 'HbbHbb_KinematicSelection.cc++("Data_V42_V8_Full", "BJetPlusX_Run2012D-PromptReco-v1_EDMNtuple_V42_ProcV1_Skim")'
say "2012D kinematically selected"
