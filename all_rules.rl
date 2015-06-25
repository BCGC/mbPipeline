#snakemake rules for microbiome pipeline
#these rules will be seperated into individual files upon completion

import os
import json

def sysio(cmd, extensions, newprefix):
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True)
    out = p.communicate()[0]
    p.wait()
    for extension in extensions :
        current = out[out[0:out.rfind(extension)].rfind("\n")+1:out[out.rfind(extension):len(out)].find("\n")+out.rfind(extension)]
        new = prefix + extension
        os.system("cp "+current+" "+new+"")

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
        taxconsensus = '{project}.final.taxconsensus'
    output:
        'AlphaDiversity.pdf',
        'BetaDiversity.pdf',
        'NumSequences.pdf',
        'TaxonomicComposition.pdf'
    run:
        min_stack_proportion = #OBTAIN FROM DEFUALTS.JSON
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
    input: '{project}.final.beta.shared'
    output: 'beta_data.out'
    run:
        #FILLER

rule process_otu:
    input:
        fasta='{project}.final.fasta',
        names='{project}.final.names',
        taxonomy='{project}.final.taxonomy',
        groups='{project}.final.groups'
    output:
        '.temp.adiv',
        '{project}.final.tax.shared',
        '{project}.final.beta.shared',
        '{project}.final.taxconsensus'
    run:
        #FILLER
        
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

        
