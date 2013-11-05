#!/bin/bash


for mass in 400 450 500 550 600 650 700 800 900 1000 1100;
do
  for toy in `seq 0 99`;
  do
    echo "mass "$mass", toy "$toy
    sed "s/w_data_$mass.root/w_data_$mass\_toy$toy.root/g" HbbHbb_19invfb_mX$mass.txt > HbbHbb_19invfb_mX$mass\_toy$toy.txt
    mkdir -p mX$mass\_$toy
    text2workspace.py HbbHbb_19invfb_mX$mass\_toy$toy.txt -o HbbHbb_19invfb_mX$mass\_toy$toy.root -L ../../GaussExp_cxx.so -L ../../ExpGaussExp_cxx.so
    echo "done text2workspace"
    combine -M MaxLikelihoodFit HbbHbb_19invfb_mX$mass\_toy$toy.root -v 2 --saveNormalizations --plot --out mX$mass\_$toy -L ../../GaussExp_cxx.so -L ../../ExpGaussExp_cxx.so > HbbHbb_19invfb_mX$mass\_toy$toy\_MaxLikelihood.log
    echo "done maxlikelihoodfit"
  done
done
