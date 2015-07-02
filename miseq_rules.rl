





rule process_sequences:
    input:
        fasta='{project}.preprocess.fasta',
        names='{project}.preprocess.names'
        groups='{project}.preprocess.groups'
    output:
        '{project}.process.fasta',
        '{project}.process.names',
        '{project}.process.groups'
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
        p = subprocess.Popen("mothur \"#set.logfile(name=master.logfile, append=T); summary.seqs(fasta="+align+", count="+input.count+")\"")#, stdout=subprocess.PIPE, shell=True)
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
        os.system("cp "+groups+" "+wildcards.project+".process.groups")

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
       # names = names[0:names.find('names')] + 'filter.names'


        # precluster to help get rid of sequencing errors - also helps with computational efficiency
        sysio_set("mothur \"#set.logfile(name=master.logfile, append=T);" +
                  "pre.cluster(fasta="+fasta+", name="+names+", group="+groups+", diffs=2)\"", [".fasta",".names"], wildcards.project+".process")

        os.system("mothur \"#set.logfile(name=master.logfile, append=T); summary.seqs(fasta="+fasta+", name="+names+")\"")






########## MISEQ PREPROCESS #############
rule pcr_sequences:
    input:
        fasta='silva.bacteria.fasta'
    output:
        '{project}.pcr.fasta' 
    run:
        with open('run.json') as data_file:
            run = json.load(data_file)
        pcrseqs_start = run["setup"]["miseq"]["pcrseqs_start"]
        pcrseqs_end = run["setup"]["miseq"]["pcrseqs_end"]
        keepdots= run["setup"]["miseq"]["keepdots"]
        sysio.set("mothur \"#set.logfile(name=master.logfile, append=T); pcr.seqs(fasta="input.fasta", start="+pcrseqs_start+", end="+pcrseqs_end+", keepdots="+keepdots+", processors=8)\"") 
        print pcrseqs_reference

rule count_sequences:
      input:
        '{project}.screen.fasta'
    output:
        '{project}.unique.fasta'
        '{project}.names'
    run:
        sysio.set("mothur \"#set.logfile(name=master.logfile, append=T); unique.seqs(fasta="+input.fasta+")\"")      
        print fasta
        print names

rule unique_sequences:
       input:
        '{project}.screen.fasta'
    output:
        '{project}.unique.fasta'
        '{project}.names'
    run:
        sysio.set("mothur \"#set.logfile(name=master.logfile, append=T); unique.seqs(fasta="+input.fasta+")\"")      
        print fasta
        print names

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
        sysio.set("mothur \"#set.logfile(name=master.logfile, append=T); screen.seqs(fasta="+input.fasta+", group="+input.groups+", maxambig="+maxambig+", maxlength="+maxlength+")\"")   
        print fasta

#####################

