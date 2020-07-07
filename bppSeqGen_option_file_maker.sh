#!/bin/bash

# this script create an option file for bppSeqGen
# It takes 11 parameters: the file containing trees for 1Mb in newick format, th ename of the output fasta file, and all the substitution model parameters (GTR)
#a, b,c,d,e, theta, theta1, theta2, alpha

a=$3
b=$4
c=$5
d=$6
e=$7
theta=$8
theta1=$9
theta2=${10}
alpha=${11}

echo "alphabet = DNA" > bppSeqGen.options_"${12}"_"${13}"
echo "input.tree.method = multiple" >> bppSeqGen.options_"${12}"_"${13}"
echo "input.tree.format = Newick" >> bppSeqGen.options_"${12}"_"${13}"
echo "input.tree.file = $1" >> bppSeqGen.options_"${12}"_"${13}"
echo "number_of_sites = 1000000" >> bppSeqGen.options_"${12}"_"${13}"
echo "output.sequence.file = $2" >> bppSeqGen.options_"${12}"_"${13}"
echo "output.sequence.format = Fasta" >> bppSeqGen.options_"${12}"_"${13}"
echo "nonhomogeneous = no" >> bppSeqGen.options_"${12}"_"${13}"
echo "model = GTR(a=$3, b=$4, c=$5, d=$6, e=$7, theta=$8, theta1=$9, theta2=$theta2)" >> bppSeqGen.options_"${12}"_"${13}"
echo "rate_distribution = Gamma(n=4, alpha=$alpha)" >> bppSeqGen.options_"${12}"_"${13}"

