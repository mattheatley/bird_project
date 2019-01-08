import os, sys

# files for which trees 2b created
indir='../paml_input/' # PAML alignments in here
alle=os.listdir(indir)
# possible species
info=open('/fastdata/bo4mhe/run_files/scripts/4tree/Species_names_lookup_table_V2Joe.csv','rU').readlines() 

treemap={}
# Start process for each file
for a in alle[:]:
	# get species in the alignment
	species=[]
	for x in open(indir+a).readlines()[1:]:
		#specnow=x[:29].rstrip().split('_')
		#specnow2=x[:29].rstrip()
		
		#specnow=x[:29].rstrip().split('_')
		specnow2=x.split(' ')[0]
		
		#if x.find('_')>-1: species.append(specnow[0]+'_'+specnow[1])
		if x.find('_')>-1: species.append(specnow2)
	outwrite=open('newlist.csv','w')
	outwrite.write(info[0])
	counter=0
	
	species.sort()
	print a, counter, len(species)
	print species
	
	
	# extract species from list
	for i in info[1:]:
		#print i.split(',')[5]
		if i.split(',')[0] in species:
			species.remove(i.split(',')[0]) 
			outwrite.write(i)
			counter+=1
	outwrite.close()
	print counter
	print species
	

	
	# rerun script
	os.system('Rscript /fastdata/bo4mhe/run_files/scripts/4tree/getTree.R')
	
	# move output
	os.system('mv MeanRateTree_NaturePaper_unroot.tre trees/'+a+'.tre')


	info=open('trees/'+a+'.tre').read()
	info=info.replace('Anser_cygnoides_domesticus','Anser_cygnoides')
	info=info.replace('Apteryx_australis_mantelli','Apteryx_australis')
	info=info.replace('Setophaga_coronata_coronata','Setophaga_coronata')
	info=info.replace('Aquila_chrysaetos_canadensis','Aquila_chrysaetos')
	info=info.replace('Zosterops_lateralis_melanops',' Zosterops_lateralis')
	info=info.replace('Buceros_rhinoceros_silvestris','Buceros_rhinoceros')
	info=info.replace('Aquila_chrysaetos_canadensis','Aquila_chrysaetos')
	info=info.replace('Tympanuchus_cupido_pinnatus','Tympanuchus_cupido')
	info=info.replace('Lyrurus_tetrix_tetrix','Lyrurus_tetrix')
	info=info.replace('Aquila_chrysaetos_canadensis','Aquila_chrysaetos')
	info=info.replace('Balearica_regulorum_gibbericeps','Balearica_regulorum')
	info=info.replace('Phoenicopterus_ruber_ruber','Phoenicopterus_ruber')
	info=info.replace('Lyrurus_tetrix_tetrix','Lyrurus_tetrix')
	info=info.replace('Struthio_camelus_australis','Struthio_camelus')
	info=info.replace('Corvus_cornix_cornix','Corvus_cornix')
	
	outree=open('final_trees/'+a+'.tre','w')
	outree.write(info)
	outree.close()



