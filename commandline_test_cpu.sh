#!/bin/bash
inputdir="/Users/gaok/Downloads/DrugCell-public/datanew/"
gene2idfile=$inputdir"gene2ind.txt"
cell2idfile=$inputdir"cell2ind_new.txt" #c
drug2idfile=$inputdir"drug2ind.txt"
testdatafile=$inputdir"drugcell_test_new.txt" #c

mutationfile=$inputdir"MUT_vec2.txt" #c
drugfile=$inputdir"drug2fingerprint.txt"

modelfile="/Users/gaok/Downloads/DrugCell-public/sample/MUT_Model_Sample_full/model_final.pt" #change
resultdir="/Users/gaok/Downloads/DrugCell-public/datanew/MUT_Result_sample_full" #change
hiddendir="/Users/gaok/Downloads/DrugCell-public/datanew/MUT_Hidden_sample_full" #change


mkdir $resultdir
mkdir $hiddendir

source activate pytorch3drugcellcpu

python -u /Users/gaok/Downloads/DrugCell-public/code/predict_drugcell_cpu.py -gene2id $gene2idfile -cell2id $cell2idfile -drug2id $drug2idfile -genotype $mutationfile -fingerprint $drugfile -hidden $hiddendir -result $resultdir -predict $testdatafile -load $modelfile > test_sample.log
