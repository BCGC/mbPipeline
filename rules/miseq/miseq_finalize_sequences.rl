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
