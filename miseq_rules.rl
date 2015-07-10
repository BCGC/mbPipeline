rule miseq_process_otu:
    input:
        fasta='{project}.final.fasta',
        count='{project}.final.count',
        taxonomy='{project}.final.taxonomy',
        groups='{project}.final.groups'
    output:
        '.temp.adiv',
        '{project}.final.tax.shared',
        '{project}.final.shared',
        '{project}.final.tax.consensus',
        '.temp.locs',
        '.temp.nums'
    run:
        with open('run.json') as data_file:
            run = json.load(data_file)
        nprocessors = run["setup"]["nprocessors"] #is there a better way of globally defining this earlier on?
        dist_cutoff = run["setup"]["distseqs_cutoff"]
        


        ### OTUs ###
        outputs = sysio_get("mothur \"#set.logfile(name=master.logfile, append=T); dist.seqs(fasta="+input.fasta+", cutoff="+distseqs_cutoff+", processors="+str(nprocessors)+")\"", [".dist"])
        dist = outputs[".dist"]
        

        outputs = sysio_get("mothur \"#set.logfile(name=master.logfile, append=T); cluster(column="+input.dist+", count="+input.count+")\"", [".list"])
        anlist = outputs[".list"]

        outputs = sysio_get("mothur \"#set.logfile(name=master.logfile, append=T); make.shared(list="+input.anlist+", count="+input.count+", label=0.03)\"", [".shared"])
        shared = outputs[".shared"]
        os.system("cp "+shared+" "+wildcards.project+".final.shared") #CHECK IF THIS IS THE APPOPRIATE WAY TO GET WILDCARD!!!

        outputs = sysio_get("mothur \"#set.logfile(name=master.logfile, append=T); classify.otu(list="+input.anlist+", count="+input.count+", taxonomy="+input.taxonomy+", label=0.03)\"", [])

        outputs = sysio_get("mothur \"#set.logfile(name=master.logfile, append=T); phylotype(taxonomy="+input.taxonomy+", label=1)\"", [".tx.list"])
        txlist = outputs[".tx.list"]

        outputs = sysio_get("mothur \"#set.logfile(name=master.logfile, append=T); make.shared(list="+input.txlist+", count="+input.count+", label=1)\"", [".shared"])
        txshared = outputs[".shared"]
        os.system("cp "+txshared+" "+wildcards.project+".final.tax.shared")

        outputs = sysio_get("mothur \"#set.logfile(name=master.logfile, append=T); classify.otu(list="+input.txlist+", count="+input.count+", taxonomy="+input.taxonomy2+", label=1)\"", [".cons.taxonomy"])
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
        


rule miseq_finalize_sequences:
    input:
        fasta='{project}.remove.fasta',
        count='{project}.remove.count',
        taxonomy = '{project}.remove.taxonomy'
    output:
        '{project}.final.fasta',
        '{project}.final.count',
        '{project}.final.taxonomy',
        '.temp.locs',
        '.temp.nums'
    run:
        # final files
        os.system("cp "+input.fasta+" "+wildcards.project+".final.fasta")
        fasta = wildcards.project+'final.fasta'
        os.system("cp "+input.count+" "+wildcards.project+".final.count")
        count = wildcards.project+'final.count'
        os.system("cp "+input.taxonomy+" "+wildcards.project+"final.taxonomy")
        taxonomy = wildcards.project+'final.taxonomy'
        run = json.load(data_file)
            groups1 = run["setup"]["miseq"]["groups1"]
            seqerror_reference = run["setup"]["miseq"]["seqerror_reference"]
            aligned = run["setup"]["miseq"]["aligned"]
            distseqs_cutoff = run["setup"]["miseq"]["distseqs_cutoff"]
        ### get sequence data ###

        os.system("rm .seq_data.out") #in case a prior file by this name existed
        os.system("mothur \"#set.logfile(name=master.logfile, append=T);count.groups(count="+input.count+")\" > .seq_data.out")

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

            print ""
            print "Warning: the following control samples have an unusually high number of sequences: " + str(ctrl_warn)
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
        print ""
        print "Warning: the following samples have an unusually low number of sequences, they will be thrown out: " + str(low_warn)

        low_seq_nums = []
        for i in range(0, len(low_warn)):
            for j in range(0, len(nums)-1):
                if locs[j] == low_warn[i]:
                     low_seq_nums.append(nums[j])
        print ""
        for i in range(0, len(low_warn)):
            print low_warn[i] + " has " + low_seq_nums[i] + " sequences." #Prints those samples and their # of seqs

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
        print ""
        print("The lowest number of sequences will be set to " + str(lowest) + " from " + ideal_loc + ".")

        ### remove controls ###

        if arecontrols == "1": #THIS HAS NOT YET BEEN TESTED#######################
            sysio_set("mothur \"#set.logfile(name=master.logfile, append=T);" +
                      "remove.groups(fasta="+inputfasta+", accnos="+controlsfile+", count="+input.count+
                      ", taxonomy="+taxonomy+")\"")

            for i in range(0, len(ctrls)):
                for j in range(0, len(nums)-1):
                    if locs[j] == ctrls[i]:
                        locs.pop(j)
                        nums.pop(j)

        f = open('.temp.locs', 'w')
        for i in range(0, len(locs)):
            f.write(str(locs[i]) + " \n")
        f.close()
        f = open('.temp.nums', 'w')
        for i in range(0, len(nums)):
            f.write(str(nums[i]) + " \n")
        f.close()

        outputs=sysio_get("mothur \"#set.logfile(name=master.logfile, append=T); get.groups(count="+input.count+", fasta="+input.fasta+", groups="+groups1+")\"")
        fasta=outputs[".fasta"]
        count=outputs[".count"]


 

        #measure the error rates
        output=sysio_get("mothur \"#set.logfile(name=master.logfile, append=T); seq.error(fasta="+input.fasta+", reference="+seqerror_reference+", aligned="+aligned+")\"")
        summary=outputs["summary"]


        #can use R to get actual error rates
        #s <- read.table(file="stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.error.summary", header=T)
        #ct <- read.table(file="stability.trim.contigs.good.unique.good.filter.unique.precluster.uchime.pick.pick.pick.count_table", header=T)
        outputs=sysio_get("Rscript seqerror.R  "+input.summary+"  "+input.count+"") 
        #/Users/browndd/Desktop/MiseqPipeline
        #This string of commands will produce a file for you
        outputs=sysio_get("mothur \"#set.logfile(name=master.logfile, append=T); dist.seqs(fasta="+input.fasta+", cutoff="+distseqs_cutoff+")\"")
        column=outputs[".column"]

        outputs=sysio_get("mothur \"#set.logfile(name=master.logfile, append=T); cluster(column="+input.column+", count="+input.count+")\"")     #runs well
        list=outputs[".list"]
        os.system("mothur \"#set.logfile(name=master.logfile, append=T); make.shared(list="+input.list+", count="+input.count+", label=0.03)\"")
        shared=outputs[".shared"]
        outputs=sysio_get("mothur \"#set.logfile(name=master.logfile, append=T); rarefaction.single(shared="+input.shared+")\"")

        #os.system("mothur \"#set.logfile(name=master.logfile, append=T); summary.seqs(fasta = "+shared+")\"") # need fixed




        ###Preparing for analysis

        #we want to remove the Mock sample from our dataset
        outputs=sysio_get("mothur \"#set.logfile(name=master.logfile, append=T); remove.groups(count="+input.count+", fasta="+input.fasta+", taxonomy="+input.taxonomy2+", groups="+groups1+")\""[".fasta",".count",".taxonomy"], wildcards.project+".final")




rule remove_miseq:
    input:
        fasta ='{project}.process.fasta',
        count ='{project}.process.fasta'
        
    output:
        '{project}.remove.fasta',
        '{project}.remove.count',
        '{project}.remove.taxonomy'
    run:
        with open('run.json') as data_file:
            run = json.load(data_file)
        nprocessors = run["setup"]["nprocessors"]
        dereplicate = run["setup"]["miseq"]["dereplicate"]
        classify.seqs_reference = run["setup"]["miseq"]["classify.seqs_reference"]
        
        # identify likely chimeras
        outputs = sysio_get("mothur \"#set.logfile(name=master.logfile, append=T);" +
                            "chimera.uchime(fasta="+input.fasta+", count="+input.count+",dereplicate="+dereplicate+", processors="+str(nprocessors)+")\"", [".accnos",".count"])

        accnos = outputs[".accnos"]
        count = outputs[".count"]
        tmp = numpy.genfromtxt(accnos, dtype='str')

        # remove identified chimeras, throwing exception if all sequences were flagged as chimeras
        fasta = ""
        if tmp.shape[0] > 0:
            outputs = sysio_get("mothur \"#set.logfile(name=master.logfile, append=T);" +
                                "remove.seqs(accnos="+accnos+", fasta="+fasta+")\"", [".fasta"])
            fasta = outputs[".fasta"]
        else:
            raise Exception("All sequences flagged as chimeras!")


        # classify sequences using given taxonomy trainset
        #os.system()
        outputs = sysio_get("mothur \"#set.logfile(name=master.logfile, append=T);" +
                            "classify.seqs(fasta="+input.fasta+", count="+input.count+",
                             reference="+classifyseqs_reference+", taxonomy="+taxonomy+".tax, cutoff=80, processors="+str(nprocessors)+")\"", [".taxonomy"])


        taxonomy2 = outputs[".taxonomy"]

        # remove contaminant mitochondria/chloroplast sequences
        sysio_set("mothur \"#set.logfile(name=master.logfile, append=T);" + 
                "remove.lineage(fasta="+input.fasta+", count="+input.count+", taxonomy="+input.taxonomy2+",
                     taxon="+taxon+")\"", [".taxonomy",".count",",fasta"], wildcards.project+".remove")
        





rule process_sequences:
    input:
        fasta='{project}.preprocess.fasta',
        names='{project}.preprocess.count'
    output:
        '{project}.process.fasta',
        '{project}.process.count'
        
    run:
        with open('run.json') as data_file:
            run = json.load(data_file)
        nprocessors = run["setup"]["nprocessors"]
        pcrseqs_reference = run["setup"]["miseq"]["pcrseqs_reference"]
        pcrseqs_start = run["setup"]["miseq"]["pcrseqs_start"]
        keepdots = run ["setup"]["miseq"]["keepdots"]
        screenseqs_start = run ["setup"]["miseq"]["screenseqs_start"]
        screenseqs_end = run ["setup"]["miseq"]["screenseqs_end"]
        maxhomop = run ["setup"]["miseq"]["maxhomop"]
        vertical = run ["setup"]["miseq"]["vertical"]
        trump = run ["setup"]["miseq"]["trump"]
        

       
        outputs = sysio_get("mothur \"#set.logfile(name=master.logfile, append=T);"+
                            "pcr.seqs(fasta="+pcrseqs_reference+", start="+pcrseqs_start+",end="+pcrseqs_end+", keepdots="+keepdots+", processors=8)\"",[".fasta"])
        #name variable [".fasta"]
        outputs = sysio_get("mothur \"#set.logfile(name=master.logfile, append=T);" +
                            "align.seqs(fasta="+input.fasta+", reference="+pcrseqs_reference+", flip=F, processors="+str(nprocessors)+")\"", [".align"])

        align= outputs[".align"]
        p = subprocess.Popen("mothur \"#set.logfile(name=master.logfile, append=T); summary.seqs(fasta="+align+", count="+input.count+")\"", stdout=subprocess.PIPE, shell=True)
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
            align = outputs[".align"]
            p = subprocess.Popen("mothur \"#set.logfile(name=master.logfile, append=T); summary.seqs(fasta="+input.align+", count="+input.count+")\"")#, stdout=subprocess.PIPE, shell=True)
            out = p.communicate()[0]
            p.wait()
            out = out[out.find("97.5%-tile:")+12:len(out)]
            out = out[out.find("\t")+1:len(out)]
            out = out[out.find("\t")+1:len(out)]
            nbasesafter = out[0:out.find("\t")]
            if int(nbasesafter)/int(nbases) <= 0.5 :
                raise Exception("Error in aligning sequences! nbases too low.")
            print("Flipping was successful!")

        fasta = align
# screen the sequences so we only keep the stuff in the region we are interested in :)
        # 0:seqname 1:start 2:end 3:nbases 4:ambigs 5:polymer 6:numSeqs
        outputs = sysio_get("mothur \"#set.logfile(name=master.logfile, append=T); summary.seqs(fasta="+fasta+", count="+input.count+")\"", [".summary"])
        summary = outputs[".summary"]
        summ = numpy.genfromtxt(summary, skiprows=1, dtype='str')
        end = map(int, summ[:,2])

        if numpy.percentile(end, 25) != numpy.percentile(end, 75):
            warnings.warn("Sequence endings are not consistent. Check to see if they have been flipped.", Warning)
        end = str(int(numpy.percentile(end, 50)))

        outputs = sysio_get("mothur \"#set.logfile(name=master.logfile, append=T);" + 
                            " screen.seqs(fasta="+input.align+", count="+input.count+", summary="+input.summary+", start="+screenseqs_start+", end="+screenseqs_end+", maxhomop="+maxhomop+")\", [".align",".count",".summary"]")

        fasta = outputs[".align"]
        count = outputs[".count"]
        summary = outputs[".summary"]
        

        outputs = sysio_get("mothur \"#set.logfile(name=master.logfile, append=T); summary.seqs(fasta="+input.align+", count="+input.count+")\",[".summary",".align",".count"]")
        summary = outputs[".summary"]
        fasta = outputs[".align"]
        count = outputs[".count"]
        

        # filter the sequences so they all overlap the same region
        outputs = sysio_get("mothur \"#set.logfile(name=master.logfile, append=T);" +
                            "filter.seqs(fasta="+input.align+", vertical=T, trump=., processors="+str(nprocessors)+")\"", [".fasta"])

        fasta = outputs[".fasta"]
        print fasta

        # should get some more unique sequences
        outputs = sysio_get("mothur \"#set.logfile(name=master.logfile, append=T); unique.seqs(fasta="+input.fasta+", count="+input.count+")\"", [".fasta",".count"])

        fasta = outputs[".fasta"]
        fasta = outputs[".count"]
        print fasta
       


        # precluster to help get rid of sequencing errors - also helps with computational efficiency
        sysio_get("mothur \"#set.logfile(name=master.logfile, append=T);" +
                  "pre.cluster(fasta="+input.fasta+", count="+input.count+", diffs=2)\"", [".fasta",".count"], wildcards.project+".process")

    






########## MISEQ PREPROCESS #############
rule pcr_sequences:
    input:
        names = '{project}.unique.names'
        groups = '{project}.screen.groups'
    output:
        '{project}.preprocess.fasta',
        '{project}.preprocess.count'
    run:
        with open('run.json') as data_file:
            run = json.load(data_file)
        pcrseqs_start = run["setup"]["miseq"]["pcrseqs_start"]
        pcrseqs_end = run["setup"]["miseq"]["pcrseqs_end"]
        keepdots= run["setup"]["miseq"]["keepdots"]

        sysio.set("mothur \"#set.logfile(name=master.logfile, append=T); pcr.seqs(fasta="+#create new default  +", start="+pcrseqs_start+", end="+pcrseqs_end+", keepdots="+keepdots+", processors=8)\""[".fasta"], wildcards.project+".preprocess")                                            
         sysio.set("mothur \"#set.logfile(name=master.logfile, append=T); count.seqs(name="+input.names+", group="+input.groups+")\""[".count"], wildcards.project+".preprocess") 


rule unique_sequences:
       input:
         fasta = '{project}.screen.fasta'
    output:
        '{project}.unique.fasta'
        '{project}.unique.names'
    run:
        sysio.set("mothur \"#set.logfile(name=master.logfile, append=T); unique.seqs(fasta="+input.fasta+")\""[".fasta",".names"], wildcards.project+".unique")    
        

rule screen_sequences:
      input: 
        fasta='{project}.fasta' 
        groups='{project}.groups'
    output:
        '{project}.screen.fasta'
        '{project}.screen.groups'
    run:
         with open('run.json') as data_file:
            run = json.load(data_file)
        maxambig= run["setup"]["miseq"]["maxambig"]
        maxlength= run["setup"]["miseq"]["maxlength"]
        sysio.set("mothur \"#set.logfile(name=master.logfile, append=T); screen.seqs(fasta="+input.fasta+", group="+input.groups+", maxambig="+maxambig+", maxlength="+maxlength+")\""[".fasta",".group"], wildcards.project+".screen")  
        

#####################

