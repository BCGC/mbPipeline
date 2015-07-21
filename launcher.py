import sys
import os

args = dict()
for a in sys.argv[1:len(sys.argv)]:
    args[a[0:a.find('=')]] = a[a.find('=')+1:len(a)]

print('\nArgs:')
print(args)
print('\n')

installed = os.path.dirname(os.path.realpath(sys.argv[0]))

######### SETUP DEFAULT VALUES #############
try:
    PATH = args['workdir']
    if len(PATH) > 0:
      os.chdir(PATH)
except KeyError:
    print("Using default working directory.")
except OSError:
    print("Invalid value for workdir, using default directory.")

os.system("cp "+installed+"/defaults.json run.json")
    
with open('run.json', 'r+') as f:
        run = json.load(f)
        
        try:
            pipeline = args['pipeline']
            run["setup"]["pipeline"] = pipeline
        except KeyError:
            raise Exception("Pipeline mode not propvided!")
        if pipeline != '454' & pipeline != 'miseq':
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

        if os.path.isdir(REFPATH):
            os.system("ln -fs " + REFPATH + "/oligos.txt .")
            os.system("ln -fs " + REFPATH + "/trainset* .")
            os.system("ln -fs " + REFPATH + "/LookUp_Titanium.pat .")
            os.system("ln -fs " + REFPATH + "/silva.* .")
            os.system("ln -fs " + REFPATH + "/" + metadata + "")
            os.system("ln -fs " + REFPATH + "/" +  controlsfile + "")
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
              maxambig = args['maxambig']
              run["setup"][pipeline]["maxambig"] = maxambig
            except KeyError:
              print("Using default maxambig.")  

            try:
              maxlength = args['maxlength']
              run["setup"][pipeline]["maxlength"] = maxlength
            except KeyError:
              print("Using default maxlength.")

            try:
              pcr_start = args['pcr_start']
              run["setup"][pipeline]["pcr_start"] = pcr_start
            except KeyError:
              print("Using default pcr_start.")

            try:
              pcr_end = args['pcr_end']
              run["setup"][pipeline]["pcr_end"] = pcr_end
            except KeyError:
              print("Using default pcr_end.")

            try:
              screen_start = args['screen_start']
              run["setup"][pipeline]["screen_start"] = screen_start
            except KeyError:
              print("Using default screen_start.")

            try:
              screen_end = args['screen_end']
              run["setup"][pipeline]["screen_end"] = screen_end
            except KeyError:
              print("Using default screen_end.")

            try:
              maxhomop = args['maxhomop']
              run["setup"][pipeline]["maxhomop"] = maxhomop
            except KeyError:
              print("Using default maxhomop.")    

            try:
              taxon = args['taxon']
              run["setup"][pipeline]["taxon"] = taxon
            except KeyError:
              print("Using default taxon.")
   
            try:
              groups1 = args['groups1']
              run["setup"][pipeline]["group1"] = group1
            except KeyError:
              print("Using default group1.")
   
            try:
              seq_reference = args['seq_reference']
              run["setup"][pipeline]["seq_reference"] = seq_reference
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
            nprocessors = args['nprocessors']
            run["setup"]["nprocessors"] = nprocessors
        except KeyError:
            print("Number of processors not provided, will use default")
        try:
            classify_cutoff = args['classify_cutoff']
            run["setup"]["classify_cutoff"]
        except KeyError:
            print("Using defaults for cut_cutoff")
        try:
            dist_cutoff = args['dist_cutoff']
            run["setup"]["dist_cutoff"]
        except KeyError:
            print("Using defaults for dist_cutoff")
        try:
            silva = args['silva']
            run["setup"]["silva"]
        except KeyError:
            print("Using defaults for silva")




        f.seek(0)
        f.write(json.dumps(run))
        f.truncate()

        
with open('run.json') as data_file:
     run = json.load(data_file)
     rulefiles = run["setup"]["+pipeline+"]

with open("Snakefile", w) as f_snakefile:
    for file in rulefiles:
        with open(file) as f_rulefile:
            for line in f_rulefile:
            f_snakefile.write(line)

       
### INSERT CODE TO LAUNCH SNAKEMAKE HERE
     
