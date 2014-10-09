# defaults.py
# set up defaults for the pipeline
# Randall Johnson
# BSP CCR Genetics Core at Frederick National Laboratory
# Leidos Biomedical Research, Inc
# Created January 15, 2014
# Last Modified January 17, 2014

###### variables and defaults ######

# workdir = <current working directory>
# datadir = "/data1/johnsonra/microbiome/refdata/"
# pdiffs = 2
# bdiffs = 1
# diffs = 1 ... used in clustering

###### need to be added ######
# processors = 12

# Switch to correct working directory
try:
      PATH = args['workdir']
      if len(PATH) > 0:
            os.chdir(PATH)
except KeyError:
      print("Using default working directory.")
except OSError:
      print("Invalid value for workdir, using default directory.")

try:
      nprocessors = args['processors']
except KeyError:
      nprocessors = 1
      print("Number of processors not provided, will default to 1!")


# diffs
try:
      diffs = args['diffs']
except KeyError:
      diffs = '2'
      print("Diffs not provided, will default to 2!")

# maxambig
try:
     maxambig = args['maxambig']
except KeyError:
     maxambig = '0'
     print("Maxambig not provided, will default to 0!")

# maxlength
try:
     maxlength = args['maxlength']
except KeyError:
     maxlength = '275'
     print("Maxlength not provided, will default to 275!")

# pcr.seqs_start
try:
     pcrseqs_start = args['pcrseqs_start']
except KeyError:
     pcrseqs_start = '11894'
     print("Pcrseqs_start not provided, will default to 11894")

# pcr.seqs_end
try:
     pcrseqs_end = args['pcrseqs_end']
except KeyError:
     pcrseqs_end = '25319'
     print("Pcrseqs_end not provided, will default to 25319!")
#keepdots

try:
     keepdots = args['keepdots']
except KeyError:
     keepdots = 'F'
     print("Keepdots not provided, will default to F!")
# screen.seqs_start

try:
     screenseqs_start = args['start']
except KeyError:
     screenseqs_start = '1968'
     print("Screenseqs_start not provided, will default to 1968!")
# screen.seqs_end

try:
     screenseqs_end = args['end']
except KeyError:
     screenseqs_end = '11550'
     print("Screenseqs_end not provided, will default to 11550!")

# maxhomop

try:
     maxhomop = args['maxhomop']
except KeyError:
     maxhomop = '8'
     print("Maxhomop not provided, will default to maxhomop!")
#vertical

try:
     vertical = args['vertical']
except KeyError:
     vertical = 'T'
     print("Vertical not provided, will default to T!")
#dereplicate

try:
     dereplicate = args['dereplicate']
except KeyError:
     dereplicate = 't'
     print("Dereplicate not provided, will default to t!")
#Trump

try:
     trump = args['trump']
except KeyError:
     trump = '.'
     print("Trump not provided, will default to . !")
#classify.seqs_cutoff

try:
     classifyseqs_cutoff = args['classifyseqs_cutoff']
except KeyError:
     classifyseqs_cutoff ='80'
     print(" Classify.seqs_cutoff not provided, will default to 80!")
#aligned

try:
     aligned = args['aligned']
except KeyError:
     aligned = 'F'
     print("Aligned not provided, will default to F!")
#dist.seqs_cuttoff
try:
     distseqs_cutoff = args['distseqs_cutoff']
except KeyError:
     distseqs_cutoff = '.20'
     print("Distseqs_cuttoff not provided, will default to .20!")
# taxon
try: 
     taxon = args['taxon']
except KeyError:
     taxon = 'Chloroplast-Mitochondria-unknown-Archaea-Eukaryota'
     print("Taxon not provided, will default to Chloroplast-Mitochondria-unknown-Archaea-Eukaryota")
# groups1
try:
     groups1 = args['groups1']
except KeyError:
     groups1 = 'Mock'
     print("Groups1 not provided, will default to Mock file!")
try:
      REFPATH = args['refpath']
except KeyError:
      REFPATH = os.system("echo $MBREF")
      if REFPATH == 0 :
            REFPATH = "."
            print("Warning: Reference data path not specified! Will try and use default directory.")

if os.path.isdir(REFPATH):
      os.system("ln -fs " + REFPATH + "/silva.bacteria.fasta .")
      os.system("ln -fs " + REFPATH + "/trainset9_032012.pds.fasta .")
      os.system("ln -fs " + REFPATH + "/HMP_MOCK.v35.fasta .")
     # os.system("ln -fs " + REFPATH + "/" + metadata +"")
      #os.system("ln -fs " + REFPATH + "/" +  controlsfile + "")
else:
      raise Exception("Bad value for refpath.")
try:
      DATAPATH = args['datapath']
except KeyError:
      DATAPATH = "."
      print("Warning: Data path for sff files not specified! Will try and use default directory.")

if not os.path.isdir(DATAPATH):
      raise Exception("Bad value for datapath.")
try:
      proj = args['project']
except KeyError:
      raise Exception("Missing project name!")
      
try:
      metadata = args['metadata']
except KeyError:
<<<<<<< Updated upstream
      min_stack_proportion = "0.14"
      print("Minimum stack proportion not provided, will default to 0.14!")
=======
      raise Exception("Metadata file not propvided!")
>>>>>>> Stashed changes

try:
     controlsfile = args['controls']
     arecontrols = True
except KeyError:
     arecontrols = False
     print("Warning: No controls file was provided, assuming that there are none.")
