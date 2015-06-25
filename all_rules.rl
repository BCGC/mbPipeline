#snakemake rules for microbiome pipeline
#these rules will be seperated into individual files upon completion

import os
import json

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
        graphics_data = 'mb_graphics_data.txt',
        beta_data = 'beta_data.out',
        taxshared = '{project}.final.tax.shared'
        taxconsensus = '{project}.final.tax.consensus'
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
    input: temp_adiv='.temp.adiv', temp_locs='.temp.locs', temp_nums='.temp.nums'
    output: 'mb_graphics_data.txt',
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
            shutil.copy2("final_data.txt", "mb_graphics_data.txt")
            
rule beta_data:
    input: '{project}.final.shared'
    output: 'beta_data.out'
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

        with open('beta_data.out', 'w') as f:
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
        '.temp.adiv',
        '{project}.final.tax.shared',
        '{project}.final.shared',
        '{project}.final.taxconsensus'
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
        f = open('.temp.adiv', 'w')
        for i in range(0, len(invsimpson)):
            f.write(str(invsimpson[i]) + ' \n')
        f.close()
        
rule finalize_sequences:
    input:
        fasta='{project}.remove.fasta',
        names='{project}.remove.names',
        groups='{project}.remove.groups'
    output:
        '{project}.final.fasta',
        '{project}.final.names',
        '{project}.final.taxonomy',
        '{project}.final.groups',
        '.temp.locs',
        '.temp.nums'
    run:
        #FILLER

rule remove_454:
    input:
        fasta='{project}.processed.fasta',
        names='{project}.processed.names',
        groups='{project}.processed.groups'
    output:
        '{project}.remove.fasta',
        '{project}.remove.names',
        '{project}.remove.groups'
    run:
        #FILLER

rule remove_miseq:
    input:
        fasta='{project}.processed.fasta',
        names='{project}.processed.names',
        groups='{project}.processed.groups'
    output:
        '{project}.remove.fasta',
        '{project}.remove.names',
        '{project}.remove.groups'
    run:
        #FILLER

rule process_sequences:
    input:
        fasta='{project}.preprocess.fasta',
        names='{project}.preprocess.names'
        groups='{project}.preprocess.groups'
    output:
        '{project}.processed.fasta',
        '{project}.processed.names',
        '{project}.processed.groups'
    run:
        #FILLER

rule rename_miseq:
    input:
    output:
    run:

rule rename_454:
    input:
        fasta='{project}.unique.fasta',
        names='{project}.unique.names'
        groups='{project}.groups'
    output:
        '{project}.preprocessed.fasta',
        '{project}.preprocessed.names',
        '{project}.preprocessed.groups'
    run:
        newfasta = input.fasta[0:input.fasta.find('unique.fasta')] + 'preprocessed.fasta'
        os.system("cp "+input.fasta+" "+newfasta+"")
        newnames = input.names[0:input.names.find('unique.names')] + 'preprocessed.names'
        os.system("cp "+input.names+" "+newnames+"")
        newgroups = input.groups[0:input.fasta.find('groups')] + 'preprocessed.groups'
        os.system("cp "+input.groups+" "+newgroups+"")

rule unique_sequences:
    input:
        fasta='{project}.trim.fasta',
        names='{project}.trim.names'
    output:
        '{project}.unique.fasta',
        '{project}.unique.names'
    run:
        #FILLER

rule trim_sequences:
    input:
        fasta='{project}.fasta',
        names='{project}.names'
    output:
        '{project}.trim.fasta',
        '{project}.trim.names'
    run:
        #FILLER

########## MISEQ PREPROCESS #############
rule pcr_sequences:
    input:
    output:
    run:
        os.system("mothur \"#set.logfile(name=master.logfile, append=T); pcr.seqs(fasta="+pcrseqs_reference+", start="+pcrseqs_start+", end="+pcrseqs_end+", keepdots="+keepdots+", processors=8)\"") 

        pcrseqs_reference=pcrseqs_reference[0:pcrseqs_reference.find('fasta')] + 'pcr.fasta'    

        print pcrseqs_reference

        os.system("mothur \"#set.logfile(name=master.logfile, append=T); summary.seqs(fasta="+pcrseqs_reference+")\"")


rule count_sequences:
    input:
    output:
    run:
        os.system("mothur \"#set.logfile(name=master.logfile, append=T); count.seqs(name="+names+", group="+groups+")\"")  

        count=fasta[0:fasta.find('unique.fasta')] + 'count_table'   
        print count

        os.system("mothur \"#set.logfile(name=master.logfile, append=T); summary.seqs(fasta="+count+")\"")


rule unique_sequences:
    input:
        '{project}.good.fasta',
    output:
    run:
        os.system("mothur \"#set.logfile(name=master.logfile, append=T); unique.seqs(fasta="+fasta+")\"")   

        names=fasta[0:fasta.find('fasta')] + 'names'
        fasta=fasta[0:fasta.find('fasta')] + 'unique.fasta'   
        print fasta
        print names


rule screen_sequences:
    input: 
        '{project}.fasta', 
        '{project}.groups'
    output:
        '{project}.good.fasta', '{project}.good.groups'
    run:
        os.system("mothur \"#set.logfile(name=master.logfile, append=T); screen.seqs(fasta="+fasta+", group="+groups+", maxambig="+maxambig+", maxlength="+maxlength+")\"")   

        fasta = fasta[0:fasta.find('fasta')] + 'good.fasta'
        groups = groups[0:groups.find("groups")] + "good.groups"
        print fasta

        os.system("mothur \"#set.logfile(name=master.logfile, append=T); summary.seqs(fasta="+fasta+")\"") 

#####################

rule load_454:
    output: '{project}.fasta', '{project}.names', '{project}.groups'
    run:
        #FILLER

rule load_miseq:
    output: '{project}.fasta', '{project}.groups'
    run:
        os.system("mothur \"#set.logfile(name=master.logfile, append=T);"  "make.       contigs(file="+files+", processors="+str(nprocessors)+       ")\"")  


        fasta = files[0:files.find('files')] + 'trim.contigs.fasta'                   
        groups = files[0:files.find('files')] + 'contigs.groups' 
        print fasta
        print groups

        
