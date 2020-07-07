import numpy as np
import tskit
import math
import msprime
import os
import pandas as pd
import sys
from gwf import Workflow

gwf = Workflow()

### This workflow simulates 5 times ARG and 10*1Mb fasta alignments for a given quartet of species. 

### recover parameters for each branches in filtered_params.csv and store them in python variables: everything is included in : parameter.csv
# species1
# species2
# species3
# outgroup
# generation time
# mutation rate
# Ne
# tau1
# tau2
# tau3
# theta1
# theta2
# rho
# GTR model parameters

# Define the trio + outgroup list
lst = ['Nycticebus_bengalensis Otolemur_garnettii Lemur_catta Homo_sapiens',
 'Loris_tardigradus Nycticebus_bengalensis Otolemur_garnettii Lemur_catta',
 'Nycticebus_bengalensis Nycticebus_pygmaeus Loris_tardigradus Otolemur_garnettii',
 'Galago_moholi Otolemur_garnettii Nycticebus_bengalensis Lemur_catta',
 'Daubentonia_madagascariensis Lemur_catta Nycticebus_bengalensis Homo_sapiens',
 'Lemur_catta Microcebus_murinus Daubentonia_madagascariensis Homo_sapiens',
 'Lemur_catta Prolemur_simus Microcebus_murinus Daubentonia_madagascariensis',
 'Hylobates_lar Macaca_mulatta Saguinus_midas Lemur_catta',
 'Homo_sapiens Hylobates_lar Macaca_mulatta Saguinus_midas',
 'Homo_sapiens Pongo_abelii Hylobates_lar Macaca_mulatta',
 'Gorilla_gorilla_gorilla Homo_sapiens Pongo_abelii Macaca_mulatta',
 'Homo_sapiens Pan_troglodytes Gorilla_gorilla_gorilla Pongo_abelii',
 'Pan_paniscus Pan_troglodytes Homo_sapiens Pongo_abelii',
 'Pongo_abelii Pongo_pygmaeus Homo_sapiens Macaca_mulatta',
 'Hylobates_lar Nomascus_leucogenys Homo_sapiens Macaca_mulatta',
 'Hylobates_lar Symphalangus_syndactylus Nomascus_leucogenys Homo_sapiens',
 'Hoolock_hoolock Symphalangus_syndactylus Hylobates_lar Homo_sapiens',
 'Hylobates_lar Hylobates_pileatus Symphalangus_syndactylus Homo_sapiens',
 'Colobus_guereza Macaca_mulatta Homo_sapiens Saguinus_midas',
 'Cercopithecus_mona Macaca_mulatta Colobus_guereza Homo_sapiens',
 'Macaca_mulatta Mandrillus_sphinx Cercopithecus_mona Homo_sapiens',
 'Macaca_mulatta Macaca_silenus Mandrillus_sphinx Colobus_guereza',
 'Macaca_nemestrina Macaca_silenus Macaca_mulatta Mandrillus_sphinx',
 'Macaca_assamensis Macaca_mulatta Macaca_silenus Mandrillus_sphinx',
 'Mandrillus_sphinx Papio_hamadryas Macaca_mulatta Colobus_guereza',
 'Cercocebus_atys Mandrillus_sphinx Papio_hamadryas Colobus_guereza',
 'Mandrillus_leucophaeus Mandrillus_sphinx Cercocebus_atys Colobus_guereza',
 'Papio_hamadryas Theropithecus_gelada Mandrillus_sphinx Colobus_guereza',
 'Lophocebus_aterrimus Papio_hamadryas Theropithecus_gelada Colobus_guereza',
 'Papio_anubis Papio_hamadryas Lophocebus_aterrimus Colobus_guereza',
 'Cercopithecus_mona Chlorocebus_sabaeus Macaca_mulatta Colobus_guereza',
 'Chlorocebus_sabaeus Erythrocebus_patas Cercopithecus_mona Colobus_guereza',
 'Chlorocebus_aethiops Chlorocebus_sabaeus Erythrocebus_patas Cercopithecus_mona',
 'Cercopithecus_mitis Cercopithecus_mona Chlorocebus_sabaeus Colobus_guereza',
 'Colobus_guereza Trachypithecus_phayrei Macaca_mulatta Homo_sapiens',
 'Rhinopithecus_roxellana Trachypithecus_phayrei Colobus_guereza Macaca_mulatta',
 'Pygathrix_nemaeus Rhinopithecus_roxellana Trachypithecus_phayrei Colobus_guereza',
 'Rhinopithecus_roxellana Rhinopithecus_strykeri Pygathrix_nemaeus Colobus_guereza',
 'Colobus_guereza Piliocolobus_tephrosceles Trachypithecus_phayrei Macaca_mulatta',
 'Colobus_angolensis_palliatus Colobus_guereza Piliocolobus_tephrosceles Macaca_mulatta',
 'Pithecia_pithecia Saguinus_midas Homo_sapiens Lemur_catta',
 'Ateles_fusciceps Saguinus_midas Pithecia_pithecia Homo_sapiens',
 'Cebus_albifrons Saguinus_midas Ateles_fusciceps Homo_sapiens',
 'Cebus_albifrons Sapajus_apella Saguinus_midas Homo_sapiens',
 'Aotus_nancymaae Saguinus_midas Cebus_albifrons Homo_sapiens',
 'Callithrix_jacchus Saguinus_midas Aotus_nancymaae Homo_sapiens']

# For each branch:
for species in lst:
	#define all variables corresponding to branch specifi parameters:
	df = pd.read_csv('parameters.tsv', sep='\t', header = 0)
	Ne=df[(df['species'] == species)]['Ne'].values[0]
	generation_time=df[(df['species'] == species)]['generation_time'].values[0]
	mutation_rate=df[(df['species'] == species)]['mutation_rate'].values[0]
	tau1=df[(df['species'] == species)]['TRUE_TRUE_tau1'].values[0]
	tau2=df[(df['species'] == species)]['TRUE_TRUE_tau2'].values[0]
	tau3=df[(df['species'] == species)]['TRUE_TRUE_tau3'].values[0]
	theta1=df[(df['species'] == species)]['TRUE_TRUE_theta1'].values[0]
	theta2=df[(df['species'] == species)]['TRUE_TRUE_theta2'].values[0]
	rho=df[(df['species'] == species)]['TRUE_TRUE_rho'].values[0]
	GTR_a=df[(df['species'] == species)]['a_median'].values[0]
	GTR_b=df[(df['species'] == species)]['b_median'].values[0]
	GTR_c=df[(df['species'] == species)]['c_median'].values[0]
	GTR_d=df[(df['species'] == species)]['d_median'].values[0]
	GTR_e=df[(df['species'] == species)]['e_median'].values[0]
	GTR_theta=df[(df['species'] == species)]['GTR_theta_median'].values[0]
	GTR_theta1=df[(df['species'] == species)]['GTR_theta1_median'].values[0]
	GTR_theta2=df[(df['species'] == species)]['GTR_theta2_median'].values[0]
	alpha=df[(df['species'] == species)]['alpha_median'].values[0]
	branch_ID=df[(df['species'] == species)]['sp_short'].values[0] 
	input_list=[]
	
## Some of those parameters are rescaled: tau1, tau2, tau3, theta1, theta2, rho. They are used only in target "msprime". They are "un-"scaled INSIDE msprime_simulations.py.

	# loop on 50 (5*10*1Mb simulated) replicates
	for tree_file_ID in range(1,51):
# simulate tree sequences (10 x 1Mb for each branch ie each set of parameters, all of this x5 for replication)
		gwf.target('msprime_{}_{}'.format(branch_ID, tree_file_ID), 
			inputs=[],
			outputs=['tree_{}_{}.newick'.format(tree_file_ID, branch_ID)],
			cores=5,
			memory='20g',
			walltime= '00:10:00',
			account='Primategenomes') << """
		python msprime_simulations.py {} {} {} {} {} {} {} {} {} {} {}
		""".format(Ne, theta1, theta2, tau1, tau2, tau3, rho, generation_time, mutation_rate, tree_file_ID, branch_ID)

		# transform msprime output in bppGenSeq input -> tree_branch_length_conversion.sh
		gwf.target('tree_conversion_{}_{}'.format(branch_ID, tree_file_ID),
           inputs=['tree_{}_{}.newick'.format(tree_file_ID, branch_ID)],
		   outputs=['modified_tree_{}_{}.newick'.format(tree_file_ID, branch_ID)],
		   cores=4,
    	   	memory='20g',
		   walltime= '00:40:00',
			account='Primategenomes') << """
		python tree_branch_length_conversion.py tree_{}_{}.newick {} {} {} {}
		""".format(tree_file_ID, branch_ID, mutation_rate, generation_time, tree_file_ID, branch_ID)

		#create the bppSeqGen option file: one for each branch -> bppSeqGen_option_file_maker.sh
		gwf.target('bppSeqGen_control_file_{}_{}'.format(branch_ID, tree_file_ID), 
           inputs=['tree_{}_{}.newick'.format(tree_file_ID, branch_ID)], 
		   outputs=['bppSeqGen.options_{}_{}'.format(tree_file_ID, branch_ID)],
		   cores=1,
    	   	memory='2g',
		   walltime= '00:10:00',
			account='Primategenomes') << """
		./bppSeqGen_option_file_maker.sh modified_tree_{}_{}.newick alignment_{}_{}.fasta {} {} {} {} {} {} {} {} {} {} {}
		""".format(tree_file_ID, branch_ID, tree_file_ID, branch_ID, GTR_a, GTR_b, GTR_c, GTR_d, GTR_e, GTR_theta, GTR_theta1, GTR_theta2, alpha, tree_file_ID, branch_ID)

		# Run BppSeqGen
		gwf.target('bppSeqGen_{}_{}'.format(branch_ID, tree_file_ID),
			inputs=['bppSeqGen.options_{}_{}'.format(tree_file_ID, branch_ID), 'modified_tree_{}_{}.newick'.format(tree_file_ID, branch_ID)],
			outputs=['alignment_{}_{}.fasta'.format(tree_file_ID, branch_ID)],
			cores=5,
			memory='20g',
			walltime= '00:60:00',
			account='Primategenomes') << """
		./bppseqgen --noninteractive=yes param=bppSeqGen.options_{}_{} input.tree.file=modified_tree_{}_{}.newick output.sequence.file=alignment_{}_{}.fasta
		""".format(tree_file_ID, branch_ID, tree_file_ID, branch_ID, tree_file_ID, branch_ID)
	
		input_list.append('alignment_{}_{}.fasta'.format(tree_file_ID, branch_ID))

	# out of tree_file_ID loop: from 50 fasta generated by bppGenseq generate 5 maf files with 10 x 1Mb each, separated by more than 100 bp
	# First target to create a list of 50 fasta:
	gwf.target('cat_fasta_{}'.format(branch_ID),
			inputs=input_list,
			outputs=['liste_fasta_{}.txt'.format(branch_ID)],
			cores=1,
			memory='3g',
			walltime= '00:10:00',
			account='Primategenomes') << """
	ls alignment_*_{}.fasta > liste_fasta_{}.txt
	""".format(branch_ID, branch_ID)

	# Second target to turn the 50 fasta into 5 maf files.
	gwf.target('from_fasta_to_maf_{}'.format(branch_ID),
			inputs=['liste_fasta_{}.txt'.format(branch_ID)],
			outputs=['alignment_1_{}.maf.gz'.format(branch_ID), 'alignment_2_{}.maf.gz'.format(branch_ID), 'alignment_3_{}.maf.gz'.format(branch_ID), 'alignment_4_{}.maf.gz'.format(branch_ID), 'alignment_5_{}.maf.gz'.format(branch_ID)],
			cores=1,
			memory='8g',
			walltime= '00:40:00',
			account='Primategenomes') << """
	python from_fasta_to_maf.py liste_fasta_{}.txt {}
	gzip alignment_*_{}.maf
	""".format(branch_ID, branch_ID, branch_ID)

	# run coalHMM workflow for each maf file: 5x47 times normally. Normally the first steps of the worflow should be useless because the provided maf file is already filtered.
	for chrm in range(1,6):
		if not os.path.isdir('./{}'.format(branch_ID)):# Create directory
			os.mkdir('./{}'.format(branch_ID))
		gwf.target('coalHMM_{}_{}'.format(branch_ID, chrm),
	        inputs=['alignment_{}_{}.maf.gz'.format(chrm, branch_ID)], 
            outputs=['./{}/chr_{}/final_table.HDF'.format(branch_ID, chrm)],
            cores=1,
            memory='10g',
            walltime= '00:60:00',
			account='Primategenomes') << """
		mkdir ./{}/chr_{}
		cd ./{}/chr_{}
		python ../../../autocoalhmm-master/autocoalhmm.py species1 species2 species3 species4 species1.chrm{} ../../alignment_{}_{}.maf.gz
		""".format(branch_ID, chrm, branch_ID, chrm, chrm, chrm, branch_ID)
