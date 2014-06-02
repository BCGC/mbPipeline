# defaults.py
# set up defaults for the pipeline
# Randall Johnson
# BSP CCR Genetics Core at Frederick National Laboratory
# Leidos Biomedical Research, Inc
# Created January 15, 2014
# Last Modified January 17, 2014

###### variables and defaults ######
# project       <-- must specify
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
      metadata = args['metadata']
except KeyError:
      raise Exception("Metadata file not propvided!")

# link reference files to the current working directory
try:
      REFPATH = args['refpath']
except KeyError:
      REFPATH = os.system("echo $MBREF")
      if REFPATH == 0 :
            REFPATH = "."
            print("Warning: Reference data path not specified! Will try and use default directory.")

if os.path.isdir(REFPATH):
      os.system("ln -fs " + REFPATH + "/oligos.txt .")
      os.system("ln -fs " + REFPATH + "/trainset* .")
      os.system("ln -fs " + REFPATH + "/LookUp_Titanium.pat .")
      os.system("ln -fs " + REFPATH + "/silva.* .")
      os.system("ln -fs " + REFPATH + "/" + metadata + "")
else:
      raise Exception("Bad value for refpath.")

# link data files to the current working directory
try:
      DATAPATH = args['datapath']
except KeyError:
      DATAPATH = "."
      print("Warning: Data path for sff files not specified! Will try and use deafault directory.")

if not os.path.isdir(DATAPATH):
      raise Exception("Bad value for datapath.")

# make sure we have the project name ###
try:
      proj = args['project']
except KeyError:
      raise Exception("Missing project name!")

# pdiffs
try:
      pdiffs = args['pdiffs']
except KeyError:
      pdiffs = '2'

# bdiffs
try:
      bdiffs = args['bdiffs']
except KeyError:
      bdiffs = '1'

try:
      indvars = args['indvars']
      if not type(indvars) is list:
            print("Warning: indvars is not a list, will use all independent variables in metadata (cat. and cont.) for graphing!")
except KeyError:
      indvars = []
      print("Will use all independent variables in metadata (cat. and cont.) for graphing!")

try:
      min_stack_proportion = args['min_stack_proportion']
except KeyError:
      min_stack_proprtion = 0.14
      print("Minimum stack proportion not provided, will default to 0.14!")
