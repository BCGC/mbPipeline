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
        try:
            lowest = int(lowest)
        except ValueError:
            lowest = int(float(lowest))
        
        ### OTUs ###
        outputs = sysio_get("mothur \"#set.logfile(name=master.logfile, append=T); dist.seqs(fasta="+input.fasta+", cutoff="+dist_cutoff+", processors="+str(nprocessors)+")\"", [".dist"])
        dist = outputs[".dist"]

        outputs = sysio_get("mothur \"#set.logfile(name=master.logfile, append=T); cluster(column="+dist+", name="+input.names+")\"", [".list"])
        anlist = outputs[".list"]

        outputs = sysio_get("mothur \"#set.logfile(name=master.logfile, append=T); make.shared(list="+anlist+", group="+input.groups+", label=0.03)\"", [".shared"])
        shared = outputs[".shared"]
        os.system("cp "+shared+" "+wildcards.project+".final.shared") #CHECK IF THIS IS THE APPOPRIATE WAY TO GET WILDCARD!!!

        outputs = sysio_get("mothur \"#set.logfile(name=master.logfile, append=T); sub.sample(shared="+shared+", size="+str(lowest)+")\"", [".subsample.shared"])
        subsample_shared = outputs[".subsample.shared"]

        outputs = sysio_get("mothur \"#set.logfile(name=master.logfile, append=T); classify.otu(list="+anlist+", name="+input.names+", taxonomy="+input.taxonomy+", label=0.03)\"", [])

        outputs = sysio_get("mothur \"#set.logfile(name=master.logfile, append=T); phylotype(taxonomy="+input.taxonomy+", name="+input.names+", label=1)\"", [".tx.list"])
        txlist = outputs[".tx.list"]

        outputs = sysio_get("mothur \"#set.logfile(name=master.logfile, append=T); make.shared(list="+txlist+", group="+input.groups+", label=1)\"", [".shared"])
        txshared = outputs[".shared"]
        os.system("cp "+txshared+" "+wildcards.project+".final.tax.shared")

        outputs = sysio_get("mothur \"#set.logfile(name=master.logfile, append=T); sub.sample(shared="+txshared+", size="+str(int(lowest))+")\"", [])

        outputs = sysio_get("mothur \"#set.logfile(name=master.logfile, append=T); classify.otu(list="+txlist+", name="+input.names+", taxonomy="+input.taxonomy+", label=1)\"", [".cons.taxonomy"])
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
