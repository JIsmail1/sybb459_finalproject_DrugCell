#!/bin/bash
inputdir="/Users/gaok/Downloads/DrugCell-public/datanew/"
gene2idfile=$inputdir"gene2ind.txt"
cell2idfile=$inputdir"cell2ind_new.txt"
drug2idfile=$inputdir"drug2ind.txt"
traindatafile=$inputdir"drugcell_train_new.txt" #change
valdatafile=$inputdir"drugcell_val_new.txt" #change
ontfile=$inputdir"drugcell_ont.txt"

mutationfile=$inputdir"MUT_vec2.txt" #change
drugfile=$inputdir"drug2fingerprint.txt"

cudaid=-1

modeldir=MUT_Model_Sample_full #change

mkdir $modeldir


source activate pytorch3drugcellcpu

python -u /Users/gaok/Downloads/DrugCell-public/code/train_drugcell.py -onto $ontfile -gene2id $gene2idfile -drug2id $drug2idfile -cell2id $cell2idfile -train $traindatafile -test $valdatafile -model $modeldir -cuda $cudaid -genotype $mutationfile -fingerprint $drugfile -genotype_hiddens 6 -drug_hiddens '100,50,6' -final_hiddens 6 -epoch 100 -batchsize 5000 -genotype $mutationfile > train_sample.log

