# defaults.py
# set up defaults for the pipeline
# Dominique Brown
# adpted from Randall Johnson script
# BSP CCR Genetics Core at Frederick National Laboratory
# Leidos Biomedical Research, Inc
# Last Modified October 10, 2014

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
      files = args['files']
except KeyError:
      files = "stability.files"
      print("sff files will default") 

try:
       project = args['project']
except KeyError:
       project = "MiSeq Pipeline"
       print("Project Name: MiSeq Pipeline")
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
     taxonomy = args['taxonomy']
except KeyError:
     taxonomy = "trainset9_032012.pds.tax"
try:
     label = args['label']
except KeyError:
     label = "0.03"
try:
     label2 = args['label2']
except KeyError:
     label2 = "1"
try:
      REFPATH = args['refpath']
except KeyError:
      REFPATH = os.system("echo $MBREF")
      if REFPATH == 0 :
            REFPATH = "."
            print("Warning: Reference data path not specified! Will try and use default directory.")
try:
      pcrseqs_reference = args['pcrseqs_reference']
except KeyError:
      pcrseqs_reference = "silva.bacteria.fasta"
try:
     classifyseqs_reference = args['classifyseqs_reference']
except KeyError:
     classifyseqs_reference = "trainset9_032012.pds.fasta"
try:
     seqerror_reference = args['seqerror_reference']
except KeyError:
     seqerror_reference = "HMP_MOCK.v35.fasta"

if os.path.isdir(REFPATH):
      os.system("ln -fs " + REFPATH + "/+pcr.seqs_reference+ .")
      os.system("ln -fs " + REFPATH + "/+classify.seqs_reference+ .")
      os.system("ln -fs " + REFPATH + "/+seq.error_reference+ .")
else:
      raise Exception("Bad value for refpath.")





















