rule load:
    input: sff=expand('{filename}.sff', filename = DATA_FILE_NAMES)
    output: '{project}.load.fasta', '{project}.load.names'
    run:
        with open('run.json') as data_file:
            run = json.load(data_file)
        DATAPATH = run["setup"]["datapath"]
        nprocessors = run["setup"]["nprocessors"]
        pdiffs = run["setup"]["454"]["pdiffs"]
        bdiffs = run["setup"]["454"]["bdiffs"]
        project = run["setup"]["proj"]
        
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
                  "shhh.flows(file="+flows+", processors="+str(nprocessors)+")\"", [".fasta",".names"], project+".load")

        # check our sequences as of right now
        # 0:seqname 1:start 2:end 3:nbases 4:ambigs 5:polymer 6:numSeqs
        outputs = sysio_get("mothur \"#set.logfile(name=master.logfile, append=T); summary.seqs(fasta="+project+".load.fasta, name="+project+".load.names)\"", [".summary"])
        summary = outputs[".summary"]
        summ = numpy.genfromtxt(summary, skiprows=1, dtype='str')

        tmp = 0
        for i in summ[:,3]:
            if int(i) < 200: # count number of reads that are less than 200 bp long
                tmp += 1

        if tmp / summ.shape[0] > 0.2:
            warnings.warn(str(tmp / summ.shape[0] * 100) +
                          "% of unique reads are shorter than 200 bp.", Warning)
