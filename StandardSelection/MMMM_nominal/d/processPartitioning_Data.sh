root -l -b -q '../../HbbHbb_Partitioning.cc++("BJetPlusX_Run2012B-13Jul2012-v1_EDMNtuple_V42_ProcV1_Skim", 107.5)'
root -l -b -q '../../HbbHbb_Partitioning.cc++("BJetPlusX_Run2012C-24Aug2012-v2_EDMNtuple_V42_ProcV1_Skim", 107.5)'
root -l -b -q '../../HbbHbb_Partitioning.cc++("BJetPlusX_Run2012C-PromptReco-v2_EDMNtuple_V42_ProcV1_Skim", 107.5)'
root -l -b -q '../../HbbHbb_Partitioning.cc++("BJetPlusX_Run2012D-PromptReco-v1_EDMNtuple_V42_ProcV1_Skim", 107.5)'
hadd -f Histograms_BJetPlusX_Run2012BCD_Skim.root Histograms_BJetPlusX_Run2012B-13Jul2012-v1_EDMNtuple_V42_ProcV1_Skim.root Histograms_BJetPlusX_Run2012C-24Aug2012-v2_EDMNtuple_V42_ProcV1_Skim.root Histograms_BJetPlusX_Run2012C-PromptReco-v2_EDMNtuple_V42_ProcV1_Skim.root Histograms_BJetPlusX_Run2012D-PromptReco-v1_EDMNtuple_V42_ProcV1_Skim.root
