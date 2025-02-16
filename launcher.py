import sys
import os
import json

args = dict()
for a in sys.argv[1:len(sys.argv)]:
	args[a[0:a.find('=')]] = a[a.find('=')+1:len(a)]

print('\nArgs:')
print(args)
print('\n')

installed = os.path.dirname(os.path.realpath(sys.argv[0]))

if not os.path.isfile(installed+"/ref/454/silva.bacteria.fasta"):
	print("THIS SEEMS TO BE THE FIRST TIME YOU ARE USING mbPIPELINE")
	print("mbPipeline is an open source automated workflow for analysing microbiome data")
	print("Developed by Nikhil Gowda, Randall Johnson, Kira Vasquez, Dominique Brown")
	print("Written in Python (utilizing Snakemake)")
	print("Uses mbGraphics package for R developed by Nikhil Gowda")
	print("Uses mothur:")
	print("Schloss, P.D., et al., Introducing mothur: Open-source, platform-independent, community-supported software for describing and comparing microbial communities. Appl Environ Microbiol, 2009. 75(23):7537-41")
	print("")
	print("")
	print("")
	print("Please wait, extracting files...")
	import zipfile
	zip_ref = zipfile.ZipFile(installed+"/ref/454/silva.bacteria.fasta.zip", 'r')
	zip_ref.extractall(installed+"/ref/454")
	zip_ref.close()
	zip_ref = zipfile.ZipFile(installed+"/ref/miseq/silva.bacteria.fasta.zip", 'r')
	zip_ref.extractall(installed+"/ref/miseq")
	zip_ref.close()
	#INSERT OTHER MISEQ EXTRACTIONS HERE
	print("")
	print("Extraction complete")
	print("")

print("Initializing mbPipeline")
print("")
print("Please read any following warnings or messages before launching:")
print("")

######### SETUP DEFAULT VALUES #############
try:
	PATH = args['workdir']
	if len(PATH) > 0:
		os.chdir(PATH)
except KeyError:
	print("Using default working directory.")
except OSError:
	print("Invalid value for workdir, using default directory.")

new = True
if os.path.isfile("run.json"):
	new = False
	print("WARNING: run.json file already in directory. WILL NOT OVERWRITE WITH DEFAULTS OR SPECIFIED ARGUMENTS UNLESS new=true")
	try:
		choice = args['new']
		if choice.lower() == "true":
			new = True
			print("    -> new=true specified, overwriting run.json with defaults and specified arguments!")
			os.system("cp "+installed+"/defaults.json run.json")
	except KeyError:
		print("    -> new!=true or new not specified, will not overwrite!")
		print("    -> However you may manually write in values to run.json and restart.")
else:
	os.system("cp "+installed+"/defaults.json run.json")

if new:
	with open('run.json', 'r+') as f:
		run = json.load(f)
		pipeline = ""
		try:
			pipeline = args['pipeline']
			run["setup"]["pipeline"] = pipeline
		except KeyError:
			raise Exception("Pipeline mode not propvided!")
		if (pipeline != '454') & (pipeline != 'miseq'):
			raise Exception("Proper pipeline name not specified!")

		try:
			metadata = args['metadata']
			run["setup"]["metadata"] = metadata
		except KeyError:
			raise Exception("Metadata file not propvided!")

		try:
			controlsfile = args['controls']
			run["setup"]["controlsfile"] = controlsfile
			run["setup"]["arecontrols"] = "1"
		except KeyError:
			run["setup"]["arecontrols"] = "0"
			print("Warning: No controls file was provided, assuming that there are none.")

		# link reference files to the current working directory
		try:
			REFPATH = args['refpath']
		except KeyError:
			REFPATH = os.system("echo $MBREF")
			if REFPATH == 0 :
				REFPATH = "."
				print("Warning: Reference data path not specified! Will try and use default directory.")

		src_files = os.listdir(installed+ "/ref/" +pipeline+ "/")
		for file_name in src_files:
			full_file_name = os.path.join(installed+ "/ref/" +pipeline+ "/" +file_name)
			if(os.path.isfile(full_file_name)):
				os.system("ln -fs " + full_file_name + " .")

		if os.path.isdir(REFPATH):
			if (pipeline==454) and (not os.path.isfile(REFPATH +"/oligos.txt")):
				raise Exception("Oligos file not in reference directory! Must be named oligos.txt!")
			os.system("ln -fs " + REFPATH + "/oligos.txt .")
			os.system("ln -fs " + REFPATH + "/trainset* .")
			os.system("ln -fs " + REFPATH + "/LookUp_Titanium.pat .")
			os.system("ln -fs " + REFPATH + "/silva.* .")
			if not os.path.isfile(REFPATH + "/" + metadata + ""):
				raise Exception("Matadata file not in reference directory!")
			os.system("ln -fs " + REFPATH + "/" + metadata + " .")
			if run["setup"]["arecontrols"] == "1":
				if not os.path.isfile(REFPATH + "/" + controlsfile + ""):
					raise Exception("Controls file not in reference directory!")
				os.system("ln -fs " + REFPATH + "/" +  controlsfile + " .")
		else:
			raise Exception("Bad value for refpath.")

		# link data files to the current working directory
		try:
			DATAPATH = args['datapath']
			run["setup"]["datapath"] = DATAPATH
		except KeyError:
			DATAPATH = "."
			print("Warning: Data path for sff files not specified! Will try and use deafault directory.")

		if not os.path.isdir(DATAPATH):
			raise Exception("Bad value for datapath.")
		
		try:
			trainset = args['trainset']
			run["setup"][pipeline]["trainset"] = trainset
		except KeyError:
			print("Warning: Using default trainset!")


		# make sure we have the project name ###
		try:
			proj = args['project']
			run["setup"]["proj"] = proj
		except KeyError:
			raise Exception("Missing project name!")

		if pipeline == "454":
			# pdiffs
			try:
				pdiffs = args['pdiffs']
				run["setup"][pipeline]["pdiffs"] = pdiffs
			except KeyError:
				print("Using default pdiffs.")

			# bdiffs
			try:
				bdiffs = args['bdiffs']
				run["setup"][pipeline]["bdiffs"] = bdiffs
			except KeyError:
				print("Using default bdiffs")
		if pipeline == "miseq":

			try:
				seqerror_reference = args['seq_reference']
				run["setup"][pipeline]["seq_reference"] = seqerror_reference
			except KeyError: 
				print("Using default seq_reference.")


		try:
			indvars = args['indvars']
			if not type(indvars) is list:
				print("Warning: indvars is not a list, will use all independent variables in metadata (cat. and cont.) for graphing!")
			else:
				run["setup"]["indvars"] = indvars
		except KeyError:
			print("Will use all independent variables in metadata (cat. and cont.) for graphing!")

		try:
			min_stack_proportion = args['min_stack_proportion']
			run["setup"]["min_stack_proportion"] = min_stack_proportion
		except KeyError:
			print("Minimum stack proportion not provided, will use default!")

		try:
			nprocessors = args['processors']
			globalprocessors = os.system("echo $SLURM_NTASKS")
			if globalprocessors != 0 :
				nprocessors = globalprocessors
				print("Determined number of processors from global variable.")
			run["setup"]["nprocessors"] = nprocessors
		except KeyError:
			print("Number of processors not provided, will use default")
		try:
			classify_cutoff = args['classify_cutoff']
			run["setup"]["classify_cutoff"] = classify_cutoff
		except KeyError:
			print("Using defaults for classify_cutoff")
		try:
			dist_cutoff = args['dist_cutoff']
			run["setup"]["dist_cutoff"] = dist_cutoff
		except KeyError:
			print("Using defaults for dist_cutoff")
		try:
			silva = args['silva']
			run["setup"]["silva"] = silva
		except KeyError:
			print("Using defaults for silva")

		f.seek(0)
		f.write(json.dumps(run))
		f.truncate()

pipeline = ""
with open('run.json') as data_file:
	run = json.load(data_file)
	pipeline = run["setup"]["pipeline"]

if pipeline == "454":
	with open('run.json', 'r+') as f:
		run = json.load(f)
		DATAPATH = run["setup"]["datapath"]
		import subprocess
		sff = subprocess.Popen('find '+DATAPATH+' -name *.sff', shell = True, stdout=subprocess.PIPE).communicate()[0]
		sff = sff.strip()
		sff = sff.decode("UTF-8")
		sff = sff.rsplit("\n")
		sff_file_names = [txt[:-4] for txt in sff]
		print("SFF FILES:")
		print(sff_file_names)

		run["setup"]["data"] = sff_file_names

		f.seek(0)
		f.write(json.dumps(run))
		f.truncate()
elif pipeline == "miseq":
	with open('run.json', 'r+') as f:
		run = json.load(f)
		DATAPATH = run["setup"]["datapath"]
		import subprocess
		fastq = subprocess.Popen('find '+DATAPATH+' -name *.fastq', shell = True, stdout=subprocess.PIPE).communicate()[0]
		fastq = fastq.strip()
		fastq = fastq.rsplit('\n')
		fastq_file_names = [txt[:-4] for txt in fastq]
		print("FASTQ FILES:")
		print(fastq_file_names)

		run["setup"]["data"] = fastq

		f.seek(0)
		f.write(json.dumps(run))
		f.truncate()

	raise Exception("NOT YET SUPPORTED")
	#miseq data setup goes here
else:
	raise Exception("ERROR: PIPELINE SPECIFICATION IN RUN.JSON IS INCORRECT. PLEASE CHECK.")

nprocessors = "1"				
with open('run.json') as data_file:
	run = json.load(data_file)
	rulefiles = run["setup"][pipeline]["rules"]
	nprocessors = run["setup"]["nprocessors"]

stitch = True
from sys import version_info

if os.path.isfile("Snakefile"):
	py3 = version_info[0] > 2 #creates boolean value for test that Python major version > 2

	if py3:
		response = input("WARNING SNAKEFILE EXISTS IN DIRECTORY. Overwrite? (yes or no) : ")
		if not((response == "yes") | (response == "y")) :
			stitch = False

	else:
		response = raw_input("WARNING SNAKEFILE EXISTS IN DIRECTORY. Overwrite? (yes or no) : ")
		if not((response == "yes") | (response == "y")) :
			stitch = False

if stitch:
	with open("Snakefile", "w") as f_snakefile:
		for filename in rulefiles:
			with open(installed+"/rules/"+filename+".rl") as f_rulefile:
				for line in f_rulefile:
					f_snakefile.write(line)

print("")


py3 = version_info[0] > 2 #creates boolean value for test that Python major version > 2

if py3:
	response = input("Setup complete, launch pipeline? (yes or no) : ")
	if not((response == "yes") | (response == "y")) :
		raise Exception("TERMINATING: User did not approve initiation of pipeline.")

else:
	response = raw_input("Setup complete, launch pipeline? (yes or no) : ")
	if not((response == "yes") | (response == "y")) :
		raise Exception("TERMINATING: User did not approve initiation of pipeline.")

print("")
print("LAUNCHING SNAKEMAKE")
if new:
	os.system("snakemake -j "+nprocessors+" --forceall")
else:
	os.system("snakemake -j "+nprocessors+"")

#example launcher start command:
#python /Users/gowdanb/Desktop/mbPipeline/Unified/launcher.py pipeline=454 metadata=nometa.txt refpath=../reference datapath=../data trainset=trainset7_112011.pds project=test1
		 
