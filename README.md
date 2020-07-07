# coalHMM_bias_estimation_workflow

## This gwf workflow evaluates the bias on population parameters (speciation times, ancestral population sizes, recombination rate) estimated by coalHMM.
## It requitres values for each parameters for a given branch (i.e. species quartet) (mutation rate, generation time, starting Ne, GTR  model parameters as well as tau1, tau2, tau3, theta1, theta2, rho as given by coalHMM i.e. rescaled) 
## It frist simulates full ARG with msprime under the form of a series of marginal trees for each recombination segment in newick format. 
## It then generate mutations on those marginal trees following the GTR model using bppseqgen and generates a series of 1Mb fasta alignemnts. 
## Those fasta alignments are combined and turned into a zipped maf file using biopython
## Finally it runs autocoalhmm (rivasiker/autocoalhmm) to allow the comparison of simulated vs estimated parameters and thus an estimation of the bias of coalHMM for each parameter.

## Warnings: in addtion to the content of the repository, the workflow requires 
#     -a tabular file with values for simulated parameters for each branch
#     -bio++ libraries 
#     -the directory ../autocoalhmm-master (path relative to the directory where the present workflow is run) and its content (see rivasiker/autocoalhmm) 
