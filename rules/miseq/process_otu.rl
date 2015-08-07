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
        silva = run["setup"]["miseq"]["silva"]

            
        
# check how to at on on to silva
#find way to find start and end for pcr and screen
       
        outputs = sysio_get("mothur \"#set.logfile(name=master.logfile, append=T);"+
                            "pcr.seqs(fasta="+silva+", start=11894,end=25319, keepdots=F, processors=8)\"",[".fasta"])
        silva= outputs[".fasta"]
        outputs = sysio_get("mothur \"#set.logfile(name=master.logfile, append=T);" +
                            "align.seqs(fasta="+input.fasta+", reference="+input.silva+", flip=F, processors="+str(nprocessors)+")\"", [".align"])

        align= outputs[".align"]
        p = subprocess.Popen("mothur \"#set.logfile(name=master.logfile, append=T); summary.seqs(fasta="+input.align+", count="+input.count+")\"", stdout=subprocess.PIPE, shell=True)
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
        #summ = numpy.genfromtxt(summary, skiprows=1, dtype='str')
        #end = map(int, summ[:,2])

        #if numpy.percentile(end, 25) != numpy.percentile(end, 75):
         #   warnings.warn("Sequence endings are not consistent. Check to see if they have been flipped.", Warning)
       # end = str(int(numpy.percentile(end, 50)))

        outputs = sysio_get("mothur \"#set.logfile(name=master.logfile, append=T); screen.seqs(fasta="+input.align+", count="+input.count+", summary="+input.summary+", start=1968, end=11550, maxhomop=8)\"", [".align",".count",".summary"])

        fasta = outputs[".align"]
        count = outputs[".count"]
        summary = outputs[".summary"]
        

        outputs = sysio_get("mothur \"#set.logfile(name=master.logfile, append=T); summary.seqs(fasta="+input.align+", count="+input.count+")\"",[".summary",".align",".count"])
        summary = outputs[".summary"]
        fasta = outputs[".align"]
        count = outputs[".count"]
        

        # filter the sequences so they all overlap the same region
        outputs = sysio_get("mothur \"#set.logfile(name=master.logfile, append=T);" +
                            "filter.seqs(fasta="+input.align+", vertical=T, trump=., processors="+str(nprocessors)+")\"", [".fasta"])

        fasta = outputs[".fasta"]
        print (fasta)

        # should get some more unique sequences
        outputs = sysio_get("mothur \"#set.logfile(name=master.logfile, append=T); unique.seqs(fasta="+input.fasta+", count="+input.count+")\"", [".fasta",".count"])

        fasta = outputs[".fasta"]
        fasta = outputs[".count"]

       


        # precluster to help get rid of sequencing errors - also helps with computational efficiency
        outputs = sysio_set("mothur \"#set.logfile(name=master.logfile, append=T); pre.cluster(fasta="+input.fasta+", count="+input.count+", diffs=2)\"", [".fasta",".count"], wildcards.project+".process")

    



    

