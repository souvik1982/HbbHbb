import re
import os
for file in os.listdir("."):
  if file.endswith(".txt"):
    filename=file[0:-4]
    os.system("text2workspace.py "+filename+".txt -o "+filename+".root -L ../../GaussExp_cxx.so -L ../../ExpGaussExp_cxx.so")
    os.system("combine -M Asymptotic "+filename+".root -L ../../GaussExp_cxx.so -L ../../ExpGaussExp_cxx.so >& "+filename+".log")
    # os.system("rm "+filename+".root")
    print "Found limit for "+filename

