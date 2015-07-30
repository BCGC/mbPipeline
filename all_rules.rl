#snakemake rules for microbiome pipeline
#these rules will be seperated into individual files upon completion

import os
import json

with open('run.json') as data_file:
    run = json.load(data_file)
PROJECT = run["setup"]["proj"]
SFF_FILE_NAMES = run["setup"]["sff"]
#SFF_FILE_NAMES = ["one", "two"]
#PROJECT='test'
def sysio_set(cmd, extensions, newprefix):
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True)
    out = p.communicate()[0]
    p.wait()
    for extension in extensions :
        current = out[out[0:out.rfind(extension)].rfind("\n")+1:out[out.rfind(extension):len(out)].find("\n")+out.rfind(extension)]
        new = prefix + extension
        os.system("cp "+current+" "+new+"")

def sysio_get(cmd, extensions):
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True)
    out = p.communicate()[0]
    p.wait()
    outputs = {}
    for extension in extensions :
        current = out[out[0:out.rfind(extension)].rfind("\n")+1:out[out.rfind(extension):len(out)].find("\n")+out.rfind(extension)]
        outputs[extension] = current
    return outputs

rule all:
    input:
        'AlphaDiversity.pdf',
        'BetaDiversity.pdf',
        'NumSequences.pdf',
        'TaxonomicComposition.pdf'

rule graph:
    input:
        graphics_data = '{project}.mb_graphics_data.txt'.format(project=PROJECT),
        beta_data = '{project}.beta_data.out'.format(project=PROJECT),
        taxshared = '{project}.final.tax.shared'.format(project=PROJECT),
        taxconsensus = '{project}.final.tax.consensus'.format(project=PROJECT)
    output:
        'AlphaDiversity.pdf',
        'BetaDiversity.pdf',
        'NumSequences.pdf',
        'TaxonomicComposition.pdf'
    run:
        with open('run.json') as data_file:
            run = json.load(data_file)
        min_stack_proportion = run["setup"]["min_stack_proportion"]
        os.system("Rscript graphall.R "+input.taxconsensus+" "+input.taxshared+" "+min_stack_proportion+"")

rule data_setup:
    input: temp_adiv='{project}.temp.adiv', temp_locs='{project}.temp.locs', temp_nums='{project}.temp.nums'
    output: '{project}.mb_graphics_data.txt',
    run:
        seqs = ["meta", "nseqs"]
        adiv = ["meta", "adiv"]
        barcode = ["meta", "Barcode"]
        with open('run.json') as data_file:
            run = json.load(data_file)
        num_lines = run["storage"]["lines"]
        metadata = run["setup"]["metadata"]
        indvars = run["setup"]["indvars"]

        f = open(input.temp_adiv)
        for i in range(0, num_lines) :
            seqs.append(f.readline())
        f.close()

        f = open(input.temp_adiv)
        for i in range(0, num_lines) :
            adiv.append(f.readline())
        f.close()

        f = open(input.temp_locs)
        for i in range(0, num_lines) :
            barcode.append(f.readline())
        f.close()

        for i in range(2, num_lines+2) :
            barcode[i] = barcode[i][:-2]
            adiv[i] = adiv[i][:-2]
            seqs[i] = seqs[i][:-2]

        num_lines = sum(1 for line in open(metadata))
        f1 = open(metadata)
        lines = f1.readlines()
        f2 = open("final_data.txt", "w")
        #This for loop is terribly overcoded - but hey, it works ;)
        for i in range(0, num_lines) :
            tabs = lines[i].split("\t")
            tabs[len(tabs)-1] = tabs[len(tabs)-1][0:tabs[len(tabs)-1].find('\n')]
            if i==0:
                tabs.append(seqs[i])
                tabs.append(adiv[i])
                f2.write("\t".join(tabs)+"\n")
            if i==1:
                tabs.append(seqs[i])
                tabs.append(adiv[i])
                f2.write("\t".join(tabs)+"\n")
            if i>=2:
                for j in range(2, len(barcode)) :
                    if barcode[j] in tabs: #only continues if barcode is found
                        tabs.append(seqs[j])
                        tabs.append(adiv[j])
                        f2.write("\t".join(tabs)+"\n")
        f1.close()
        f2.close()

        if not len(indvars) == 0 :
            f1 = open("final_data.txt")
            f2 = open("mb_graphics_data.txt", "w")
            lines = f1.readlines()
            numcols = len(lines[0].split("\t"))
            columns_to_ignore = []
            for i in range(0, numcols) :
                if lines[0].split("\t")[i] == "cat" or lines[0].split("\t")[i] == "cont" :
                    if not lines[1].split("\t")[i] in indvars :
                        columns_to_ignore.append(i)
            num_lines=len(lines)
            for i in range(0, num_lines) :
                tabs = lines[i].split("\t")
                tabs[len(tabs)-1] = tabs[len(tabs)-1][0:tabs[len(tabs)-1].find('\n')]
                tabs = [j for k, j in enumerate(tabs) if k not in columns_to_ignore]
                f2.write("\t".join(tabs)+"\n")
            f1.close()
            f2.close()
        else:
            import shutil 
            shutil.copy2("final_data.txt", wildcards.project+"mb_graphics_data.txt")
            
rule beta_data:
    input: '{project}.final.shared'
    output: '{project}.beta_data.out'
    run:
        outputs = sysio_get("mothur \"#summary.shared(shared="+input+", calc=thetayc)\"", [".summary"])
        summary = outputs[".summary"]

        os.system("cut -f2 "+summary+" > .temp_sample1.out")
        num_lines5 = sum(1 for line in open('.temp_sample1.out'))
        sample1 = []
        f = open('.temp_sample1.out')
        for i in range(0, num_lines5):
            sample1.append(f.readline())
        f.close()
        for i in range(0, len(sample1)):
            sample1[i] = sample1[i][:-1]
        sample1[0] = "sample1"

        os.system("cut -f3 "+summary+" > .temp_sample2.out")
        sample2 = []
        f = open('.temp_sample2.out')
        for i in range(0, num_lines5):
            sample2.append(f.readline())
        f.close()
        for i in range(0, len(sample2)):
            sample2[i] = sample2[i][:-1]
        sample2[0] = "sample2"

        os.system("cut -f5 "+summary+" > .temp_bdiv.out")
        bdiv = []
        f = open('.temp_bdiv.out')
        for i in range(0, num_lines5):
            bdiv.append(f.readline())
        f.close()
        for i in range(0, len(bdiv)):
            bdiv[i] = bdiv[i][:-1]
        bdiv[0] = "bdiv"

        os.system("cut -f7 "+summary+" > .temp_cmin.out")
        cmin = []
        f = open('.temp_cmin.out')
        for i in range(0, num_lines5):
            cmin.append(f.readline())
        f.close()
        for i in range(0, len(cmin)):
            cmin[i] = cmin[i][:-1]
        for i in range(1, len(cmin)):
            cmin[i] = 1 - float(cmin[i])
        for i in range(1, len(cmin)):
            cmin[i] = str(cmin[i])
        cmin[0] = "cmin"

        os.system("cut -f6 "+summary+" > "".temp_cmax.out")
        cmax = []
        f = open('.temp_cmax.out')
        for i in range(0, num_lines5):
            cmax.append(f.readline())
        f.close()
        for i in range(0, len(cmax)):
            cmax[i] = cmax[i][:-1]
        for i in range(1, len(cmax)):
            cmax[i] = 1 - float(cmax[i])
        for i in range(1, len(cmax)):
            cmax[i] = str(cmax[i])
        cmax[0] = "cmax"

        with open(wildcards.project+'.beta_data.out', 'w') as f:
            for f1, f2, f3, f4, f5 in zip(sample1, sample2, bdiv, cmin, cmax):
                f.write(f1+"\t"+f2+"\t"+f3+"\t"+f4+"\t"+f5+"\n")
        f.close()

rule process_otu:
    input:
        fasta='{project}.final.fasta',
        names='{project}.final.names',
        taxonomy='{project}.final.taxonomy',
        groups='{project}.final.groups'
    output:
        '{project}.temp.adiv',
        '{project}.final.tax.shared',
        '{project}.final.shared',
        '{project}.final.tax.consensus'
    run:
        with open('run.json') as data_file:
            run = json.load(data_file)
        nprocessors = run["setup"]["nprocessors"] #is there a better way of globally defining this earlier on?
        dist_cutoff = run["setup"]["dist_cutoff"]
        lowest = run["storage"]["lowest"]
        
        ### OTUs ###
        outputs = sysio_get("mothur \"#set.logfile(name=master.logfile, append=T); dist.seqs(fasta="+input.fasta+", cutoff="+dist_cutoff+", processors="+str(nprocessors)+")\"", [".dist"])
        dist = outputs[".dist"]

        outputs = sysio_get("mothur \"#set.logfile(name=master.logfile, append=T); cluster(column="+dist+", name="+input.names+")\"", [".list"])
        anlist = outputs[".list"]

        outputs = sysio_get("mothur \"#set.logfile(name=master.logfile, append=T); make.shared(list="+anlist+", group="+input.groups+", label=0.03)\"", [".shared"])
        shared = outputs[".shared"]
        os.system("cp "+shared+" "+wildcards.project+".final.shared") #CHECK IF THIS IS THE APPOPRIATE WAY TO GET WILDCARD!!!

        outputs = sysio_get("mothur \"#set.logfile(name=master.logfile, append=T); sub.sample(shared="+shared+", size="+str(int(lowest))+")\"", [".subsample.shared"])
        subsample_shared = outputs[".subsample.shared"]

        outputs = sysio_get("mothur \"#set.logfile(name=master.logfile, append=T); classify.otu(list="+anlist+", name="+input.names+", taxonomy="+input.taxonomy+", label=0.03)\"", [])

        outputs = sysio_get("mothur \"#set.logfile(name=master.logfile, append=T); phylotype(taxonomy="+input.taxonomy+", name="+input.names+", label=1)\"", [".tx.list"])
        txlist = outputs[".tx.list"]

        outputs = sysio_get("mothur \"#set.logfile(name=master.logfile, append=T); make.shared(list="+txlist+", group="+groups+", label=1)\"", [".shared"])
        txshared = outputs[".shared"]
        os.system("cp "+txshared+" "+wildcards.project+".final.tax.shared")

        outputs = sysio.get("mothur \"#set.logfile(name=master.logfile, append=T); sub.sample(shared="+txshared+", size="+str(int(lowest))+")\"", [])

        outputs = sysio_get("mothur \"#set.logfile(name=master.logfile, append=T); classify.otu(list="+txlist+", name="+names+", taxonomy="+taxonomy+", label=1)\"", [".cons.taxonomy"])
        txconsensus = outputs[".cons.taxonomy"]
        os.system("cp "+txshared+" "+wildcards.project+".final.tax.consensus")

        ### Alpha Diversity ###

        os.system("mothur \"#set.logfile(name=master.logfile, append=T); collect.single(shared="+shared+", calc=chao-invsimpson, freq=100)\"")

        sample_list = []
        os.system("grep -l '0.03' *.invsimpson > .sample_list.out")
        num_lines3 = sum(1 for line in open('.sample_list.out'))
        f = open('.sample_list.out')
        for i in range(0, num_lines3):
            sample_list.append(f.readline())
            sample_list[i] = sample_list[i][:-1]
        f.close()
        temp1 = []
        summ = 0
        invsimpson = []
        for i in range(0, num_lines3):
            os.system("cut -f2 -s "+sample_list[i]+" | tail -n 5 > .temp_nums.out")
            num_lines4 = sum(1 for line in open('.temp_nums.out'))
            f = open('.temp_nums.out')
            for j in range(0, num_lines4):
                temp1.append(f.readline())
            for z in range(0, num_lines4):
                summ += float(temp1[z])
            temp1 = []
            invsimpson.append(summ/num_lines4)
            summ = 0
            f.close()
        f = open(wildcards.project+'.temp.adiv', 'w')
        for i in range(0, len(invsimpson)):
            f.write(str(invsimpson[i]) + ' \n')
        f.close()
        
rule finalize_sequences:
    input:
        fasta='{project}.remove.fasta',
        names='{project}.remove.names',
        groups='{project}.remove.groups',
    output:
        '{project}.final.fasta',
        '{project}.final.names',
        '{project}.final.taxonomy',
        '{project}.final.groups',
        '{project}.temp.locs',
        '{project}.temp.nums'
    run:
        # final files
        os.system("cp "+input.fasta+" "+wildcards.project+".final.fasta")
        fasta = wildcards.project+'final.fasta'
        os.system("cp "+input.names+" "+wildcards.project+".final.names")
        names = wildcards.project+'final.names'
        os.system("cp "+input.groups+" "+wildcards.project+".final.groups")
        groups = wildcards.project+'final.groups'

        ### get sequence data ###

        os.system("rm .seq_data.out") #in case a prior file by this name existed
        os.system("mothur \"#set.logfile(name=master.logfile, append=T);count.groups(group="+groups+")\" > .seq_data.out")

        ### pull apart data in x.seq_data.out ###

        num_lines = sum(1 for line in open('.seq_data.out'))
        data = []
        f = open('.seq_data.out')
        for i in range(0, num_lines) :
            text = f.readline()
            if 'contains' in text:
                   data.append(text)
        f.close()
        locs = []
        nums = []
        for i in range(0, len(data)):
            data[i] = data[i][:-2]
        for i in range(0, len(data)):
            temp1,_,temp2 = data[i].partition(' contains ')
            locs.append(temp1)
            nums.append(temp2)

        ### print warnings, find optimal sequence size and save ctrl seqs to file ###
        with open('run.json') as data_file:
            run = json.load(data_file)
        arecontrols = run["setup"]["arecontrols"]

        if arecontrols=="1":
            ctrls = []
            num_lines2 = sum(1 for line in open(controlsfile))
            f = open(controlsfile)
            for i in range(0, num_lines2):
                  ctrls.append(f.readline())
            f.close()
            for i in range(0, len(ctrls)):
                  ctrls[i] = ctrls[i][:-1]
            ctrl_nums = []
            ctrl_warn = []
            ctrl_locs = []
            for i in range(0, len(ctrls)):
                  for j in range(0, len(locs)-1):
                    if ctrls[i] == locs[j]:
                     ctrl_locs.append(locs.pop(j))
                     ctrl_nums.append(nums.pop(j))
            for i in range(0, len(ctrl_nums)):
                  if float(ctrl_nums[i]) > 1000:
                     ctrl_warn.append(ctrl_locs[i])

            f = open('.control.seqs', 'w')
            for i in range(0, len(ctrls)):
                  f.write(ctrls[i] + ": " + ctrl_nums[i] + " \n")
            f.close()

            print("")
            print("Warning: the following control samples have an unusually high number of sequences: " + str(ctrl_warn))
        f = open('.temp.numseqs', 'w')
        for i in range(0, len(nums)):
            f.write(str(nums[i]) + " \n")
            f.close()
        with open('run.json', 'r+') as f:
            run = json.load(f)
            run["storage"]["lines"] = ""+sum(1 for line in open('.temp.numseqs'))
            f.seek(0)
            f.write(json.dumps(run))
            f.truncate()

        low_warn = [] #This part grabs all samples with fewer than 3000 sequences
        for i in range(0, len(nums)):
            if float(nums[i]) < 3000:
                   low_warn.append(locs[i])
        print("")
        print("Warning: the following samples have an unusually low number of sequences, they will be thrown out: " + str(low_warn))

        low_seq_nums = []
        for i in range(0, len(low_warn)):
            for j in range(0, len(nums)-1):
                if locs[j] == low_warn[i]:
                     low_seq_nums.append(nums[j])
        print("")
        for i in range(0, len(low_warn)):
            print(low_warn[i] + " has " + low_seq_nums[i] + " sequences.") #Prints those samples and their # of seqs

        #Following loop removes those found low sequences names and numbers from the orig lists
        for i in range(0, len(low_warn)):
            for j in range(0, len(nums)-1):
                if locs[j] == low_warn[i]:
                    locs.pop(j)
                    nums.pop(j)
        highest = 0 #This part finds the sample with the highest number of sequences
        for i in range(0, len(nums)):
            if float(nums[i]) > float(highest):
                highest = float(nums[i])
        lowest = highest
        with open('run.json', 'r+') as f:
            run = json.load(f)
            run["storage"]["lowest"] = lowest
            f.seek(0)
            f.write(json.dumps(run))
            f.truncate()
        
        #The following part finds the sample with the lowest number of sequences (which is consider the ideal lowest)
        for i in range(0, len(nums)):
            if float(nums[i]) < lowest:
                lowest = float(nums[i])
                ideal_loc = locs[i]
        print("")
        print("The lowest number of sequences will be set to " + str(lowest) + " from " + ideal_loc + ".")

        ### remove controls ###

        if arecontrols == "1": #THIS HAS NOT YET BEEN TESTED#######################
            sysio_set("mothur \"#set.logfile(name=master.logfile, append=T);" +
                      "remove.groups(fasta="+fasta+", accnos="+controlsfile+", group="+groups+
                      ", name="+names+".final.names, taxonomy="+taxonomy+")\"", [".fasta",".names",".groups",".taxonomy"], wildcards.project+".final")

            for i in range(0, len(ctrls)):
                for j in range(0, len(nums)-1):
                    if locs[j] == ctrls[i]:
                        locs.pop(j)
                        nums.pop(j)

        f = open(wildcards.project+'.temp.locs', 'w')
        for i in range(0, len(locs)):
            f.write(str(locs[i]) + " \n")
        f.close()
        f = open(wildcards.project+'.temp.nums', 'w')
        for i in range(0, len(nums)):
            f.write(str(nums[i]) + " \n")
        f.close()


rule remove:
    input:
        fasta='{project}.process.fasta',
        names='{project}.process.names',
        groups='{project}.process.groups'
    output:
        '{project}.remove.fasta',
        '{project}.remove.names',
        '{project}.remove.groups',
        '{project}.remove.taxonomy'
    run:
        with open('run.json') as data_file:
            run = json.load(data_file)
        nprocessors = run["setup"]["nprocessors"]
        trainset = run["setup"]["454"]["trainset"]
        
        # identify likely chimeras
        outputs = sysio_get("mothur \"#set.logfile(name=master.logfile, append=T);" +
                            "chimera.uchime(fasta="+input.fasta+", name="+input.names+", group="+input.groups+", processors="+str(nprocessors)+")\"", [".accnos"])

        accnos = outputs[".accnos"]
        tmp = numpy.genfromtxt(accnos, dtype='str')

        # remove identified chimeras, throwing exception if all sequences were flagged as chimeras
        fasta = ""
        names = ""
        groups = ""
        if tmp.shape[0] > 0:
            outputs = sysio_get("mothur \"#set.logfile(name=master.logfile, append=T);" +
                                "remove.seqs(accnos="+accnos+", fasta="+fasta+", name="+names+", " +
                                "group="+groups+", dups=T)\"", [".fasta",".names",".groups"])
            fasta = outputs[".fasta"]
            names = outputs[".names"]
            groups = outputs[".groups"]
        else:
            raise Exception("All sequences flagged as chimeras!")


        # classify sequences using given taxonomy trainset
        #os.system()
        outputs = sysio_get("mothur \"#set.logfile(name=master.logfile, append=T);" +
                            "classify.seqs(fasta="+fasta+", name="+names+", group="+groups+
                            ", template="+trainset+"=.fasta, taxonomy="+trainset+".tax, cutoff=80, processors="+str(nprocessors)+")\"", [".taxonomy"])


        taxonomy = outputs[".taxonomy"]

        # remove contaminant mitochondria/chloroplast sequences
        sysio_set("mothur \"#set.logfile(name=master.logfile, append=T);" + 
                "remove.lineage(fasta="+fasta+", name="+names+", group="+groups+", taxonomy="+taxonomy+
                  ", taxon=Mitochondria-Cyanobacteria_Chloroplast-unknown)\"", [".taxonomy",".names",".groups",",fasta"], wildcards.project+".remove")

rule process_sequences:
    input:
        fasta='{project}.unique.fasta',
        names='{project}.unique.names',
        groups='{project}.load.groups'
    output:
        '{project}.process.fasta',
        '{project}.process.names',
        '{project}.process.groups'
    run:
        with open('run.json') as data_file:
            run = json.load(data_file)
        nprocessors = run["setup"]["nprocessors"]
        silva = run["setup"]["silva"]
        # initial alignment
        outputs = sysio_get("mothur \"#set.logfile(name=master.logfile, append=T);" +
                            "align.seqs(fasta="+input.fasta+", reference="+silva+", flip=F, processors="+str(nprocessors)+")\"", [".align"])

        fastacheck = outputs[".align"]
        p = subprocess.Popen("mothur \"#set.logfile(name=master.logfile, append=T); summary.seqs(fasta="+fastacheck+", name="+input.names+")\"", stdout=subprocess.PIPE, shell=True)
        out = p.communicate()[0]
        p.wait()
        out = out[out.find("97.5%-tile:")+12:len(out)]
        out = out[out.find("\t")+1:len(out)]
        out = out[out.find("\t")+1:len(out)]
        nbasesafter = out[0:out.find("\t")]

        if int(nbasesafter)/int(nbases) <= 0.5 :
            print("Warning: Attempting to flip direction and re-allign sequences.")
            outputs = sysio_get("mothur \"#set.logfile(name=master.logfile, append=T);" +
                                "align.seqs(fasta="+input.fasta+", reference="+silva+", flip=T, processors="+str(nprocessors)+")\"", [".fasta"])
            fastacheck = outputs[".align"]
            p = subprocess.Popen("mothur \"#set.logfile(name=master.logfile, append=T); summary.seqs(fasta="+fastacheck+", name="+input.names+")\"", stdout=subprocess.PIPE, shell=True)
            out = p.communicate()[0]
            p.wait()
            out = out[out.find("97.5%-tile:")+12:len(out)]
            out = out[out.find("\t")+1:len(out)]
            out = out[out.find("\t")+1:len(out)]
            nbasesafter = out[0:out.find("\t")]
            if int(nbasesafter)/int(nbases) <= 0.5 :
                raise Exception("Error in aligning sequences! nbases too low.")
            print("Flipping was successful!")

        fasta = fastacheck
        
        # screen the sequences so we only keep the stuff in the region we are interested in :)
        # 0:seqname 1:start 2:end 3:nbases 4:ambigs 5:polymer 6:numSeqs
        outputs = sysio_get("mothur \"#set.logfile(name=master.logfile, append=T); summary.seqs(fasta="+fasta+", name="+input.names+")\"", [".summary"])
        summary = outputs[".summary"]
        summ = numpy.genfromtxt(summary, skiprows=1, dtype='str')
        end = map(int, summ[:,2])

        if numpy.percentile(end, 25) != numpy.percentile(end, 75):
            warnings.warn("Sequence endings are not consistent. Check to see if they have been flipped.", Warning)
        end = str(int(numpy.percentile(end, 50)))

        outputs = sysio_get("mothur \"#set.logfile(name=master.logfile, append=T);" +
                            "screen.seqs(fasta="+fasta+", name="+input.names+", group="+input.groups+
                            ", end="+end+", optimize=start, criteria=95, processors="+str(nprocessors)+")\"", [".align",".names",".groups"])

        fasta = outputs[".align"]
        names = outputs[".names"]
        groups = outputs[".groups"]
        os.system("cp "+groups+" "+wildcards.project+".process.groups")

        os.system("mothur \"#set.logfile(name=master.logfile, append=T); summary.seqs(fasta="+fasta+", name="+names+")\"")


        # filter the sequences so they all overlap the same region
        outputs = sysio_get("mothur \"#set.logfile(name=master.logfile, append=T);" +
                            "filter.seqs(fasta="+fasta+", vertical=T, trump=., processors="+str(nprocessors)+")\"", [".fasta"])

        fasta = outputs[".fasta"]
        print(fasta)

        # should get some more unique sequences
        outputs = sysio_get("mothur \"#set.logfile(name=master.logfile, append=T); unique.seqs(fasta="+fasta+", name="+names+")\"", [".fasta",".names"])

        fasta = outputs[".fasta"]
        print(fasta)
        names = outputs[".names"]

        os.system("mothur \"#set.logfile(name=master.logfile, append=T); summary.seqs(fasta="+fasta+", name="+names+")\"")


        # precluster to help get rid of sequencing errors - also helps with computational efficiency
        sysio_set("mothur \"#set.logfile(name=master.logfile, append=T);" +
                  "pre.cluster(fasta="+fasta+", name="+names+", group="+groups+", diffs=2)\"", [".fasta",".names"], wildcards.project+".process")

        os.system("mothur \"#set.logfile(name=master.logfile, append=T); summary.seqs(fasta="+fasta+", name="+names+")\"")

#IGNORE THIS RULE
#rule rename_454:
#    input:
#        fasta='{project}.unique.fasta',
#        names='{project}.unique.names'
#        groups='{project}.load.groups'
#    output:
#        '{project}.preprocessed.fasta',
#        '{project}.preprocessed.names',
#        '{project}.preprocessed.groups'
#    run:
#        newfasta = input.fasta[0:input.fasta.find('unique.fasta')] + 'preprocessed.fasta'
#        os.system("cp "+input.fasta+" "+newfasta+"")
#        newnames = input.names[0:input.names.find('unique.names')] + 'preprocessed.names'
#        os.system("cp "+input.names+" "+newnames+"")
#        newgroups = input.groups[0:input.fasta.find('groups')] + 'preprocessed.groups'
#        os.system("cp "+input.groups+" "+newgroups+"")

rule unique_sequences:
    input:
        fasta='{project}.trim.fasta',
        names='{project}.trim.names'
    output:
        '{project}.unique.fasta',
        '{project}.unique.names'
    run:
        sysio_set("mothur \"#unique.seqs(fasta="+input.fasta+", name="+input.names+")\"", [".fasta",".names"], wildcards.project+".unique")

        p = subprocess.Popen("mothur \"#set.logfile(name=master.logfile, append=T); summary.seqs(fasta="+fasta+", name="+names+")\"", stdout=subprocess.PIPE, shell=True)
        out = p.communicate()[0]
        p.wait()
        out = out[out.find("97.5%-tile:")+12:len(out)]
        out = out[out.find("\t")+1:len(out)]
        out = out[out.find("\t")+1:len(out)]
        nbases = out[0:out.find("\t")]
        with open('run.json', 'r+') as f:
            run = json.load(f)
            run["storage"]["nbases"] = nbases
            f.seek(0)
            f.write(json.dumps(run))
            f.truncate()

rule trim_sequences:
    input:
        fasta='{project}.load.fasta',
        names='{project}.load.names'
    output:
        '{project}.trim.fasta',
        '{project}.trim.names'
    run:
        with open('run.json') as data_file:
            run = json.load(data_file)
        pdiffs = run["setup"]["454"]["pdiffs"]
        bdiffs = run["setup"]["454"]["bdiffs"]
        nprocessors = run["setup"]["nprocessors"]
        sysio_set("mothur \"#set.logfile(name=master.logfile, append=T); trim.seqs(fasta="+input.fasta+
                  ", name="+input.names+", oligos=oligos.txt, pdiffs="+pdiffs+", bdiffs="+bdiffs+
                  ", maxhomop=8, minlength=200, flip=T processors="+str(nprocessors)+")\"", [".fasta", ".names"], wildcards.project+".trim")


rule load:
    input: sff=expand('{filename}.sff', filename = SFF_FILE_NAMES)
    output: '{project}.load.fasta', '{project}.load.names', '{project}.load.groups'
    run:
        with open('run.json') as data_file:
            run = json.load(data_file)
        DATAPATH = run["setup"]["datapath"]
        nprocessors = run["setup"]["nprocessors"]
        pdiffs = run["setup"]["454"]["pdiffs"]
        bdiffs = run["setup"]["454"]["bdiffs"]
        project = run["setup"]["proj"]
        
        #sff = subprocess.Popen('find '+DATAPATH+' -name *.sff', shell = True, stdout=subprocess.PIPE).communicate()[0]
        #sff = sff.rsplit('\n')
        sff = input.sff

        # this is the repository for all sff files
        os.system("printf '' > all.flow.files")
        os.system("printf '' > master.logfile")

        for f in sff:
            if os.path.isfile(f):
                x = f[0:f.find('.sff')]
                head, tail = os.path.split(f)
                y = tail[0:tail.find('.sff')]
                os.system("mothur \"#set.logfile(name=master.logfile, append=T); " +
                          "set.dir(output=.); " +
                          "sffinfo(sff="+x+".sff); " +
                          "set.dir(input=.); " +
                          "summary.seqs(fasta="+y+".fasta); " +
                          "trim.flows(flow="+y+".flow, oligos=oligos.txt, pdiffs="+pdiffs+","+"bdiffs="+bdiffs+", processors="+str(nprocessors)+")\"")
                os.system("cat "+y+".flow.files >> all.flow.files")


        flows = 'all.flow.files'
        sysio_set("mothur \"#set.logfile(name=master.logfile, append=T);" +
                  "shhh.flows(file="+flows+", processors="+str(nprocessors)+")\"", [".fasta",".names",".groups"], project+".load.")

        # check our sequences as of right now
        # 0:seqname 1:start 2:end 3:nbases 4:ambigs 5:polymer 6:numSeqs
        outputs = sysio_get("mothur \"#set.logfile(name=master.logfile, append=T); summary.seqs(fasta="+proj+".fasta, name="+proj+".names)\"", [".summary"])
        summary = outputs[".summary"]
        summ = numpy.genfromtxt(summary, skiprows=1, dtype='str')

        tmp = 0
        for i in summ[:,3]:
            if int(i) < 200: # count number of reads that are less than 200 bp long
                tmp += 1

        if tmp / summ.shape[0] > 0.2:
            warnings.warn(str(tmp / summ.shape[0] * 100) +
                          "% of unique reads are shorter than 200 bp.", Warning)

        
