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
               try:
                   maxlength = args['maxlength']
                try:
                    pcr_start = args['pcr_start']
                try:
                    pcr_end = args['pcr_end']
                try:
                    screen_start = args['screen_start']
                try:
                    screen_end = args['screen_seqs']
                try:
                    maxhomop = args['maxhomop']
                try:
                    taxon = args['taxon']
                try:
                    groups1 = args['groups1']
                try:
                    seq_reference = args['seq_reference']





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
                run["setup"]["nprocessors"] = nprocessors
            except KeyError:
                print("Number of processors not provided, will use default")


            f.seek(0)
            f.write(json.dumps(run))
            f.truncate()


