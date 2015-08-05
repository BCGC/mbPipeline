rule finalize_sequences:
    input:
        fasta='{project}.remove.fasta',
        names='{project}.remove.names',
        groups='{project}.remove.groups',
        taxonomy='{project}.remove.taxonomy'
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
        fasta = wildcards.project+'.final.fasta'
        os.system("cp "+input.names+" "+wildcards.project+".final.names")
        names = wildcards.project+'.final.names'
        os.system("cp "+input.groups+" "+wildcards.project+".final.groups")
        groups = wildcards.project+'.final.groups'
        os.system("cp "+input.taxonomy+" "+wildcards.project+".final.taxonomy")
        taxonomy = wildcards.project+'.final.taxonomy'

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
            run["storage"]["lines"] = ""+str(sum(1 for line in open('.temp.numseqs')))
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
                #The following part finds the sample with the lowest number of sequences (which is consider the ideal lowest)
        ideal_loc=""
        for i in range(0, len(nums)):
            if float(nums[i]) < lowest:
                lowest = float(nums[i])
                ideal_loc = locs[i]
        print("")
        print("The lowest number of sequences will be set to " + str(lowest) + " from " + ideal_loc + ".")
        with open('run.json', 'r+') as f:
            run = json.load(f)
            run["storage"]["lowest"] = str(lowest)
            f.seek(0)
            f.write(json.dumps(run))
            f.truncate()
        
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
