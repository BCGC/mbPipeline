rule process_sequences:
    input:
        fasta='{project}.unique.fasta',
        names='{project}.unique.names',
        groups='{project}.trim.groups'
    output:
        '{project}.process.fasta',
        '{project}.process.names',
        '{project}.process.groups'
    run:
        with open('run.json') as data_file:
            run = json.load(data_file)
        nprocessors = run["setup"]["nprocessors"]
        silva = run["setup"]["silva"]
        nbases = run["storage"]["nbases"]
        # initial alignment
        outputs = sysio_get("mothur \"#set.logfile(name=master.logfile, append=T);" +
                            "align.seqs(fasta="+input.fasta+", reference="+silva+", flip=F, processors="+str(nprocessors)+")\"", [".align"])

        fastacheck = outputs[".align"]
        p = subprocess.Popen("mothur \"#set.logfile(name=master.logfile, append=T); summary.seqs(fasta="+fastacheck+", name="+input.names+")\"", stdout=subprocess.PIPE, shell=True)
        out = p.communicate()[0].decode("utf-8")
        p.wait()
        out = out[out.find("97.5%-tile:")+12:len(out)]
        out = out[out.find("\t")+1:len(out)]
        out = out[out.find("\t")+1:len(out)]
        nbasesafter = out[0:out.find("\t")]

        if int(nbasesafter)/int(nbases) <= 0.5 :
            print("Warning: Attempting to flip direction and re-allign sequences.")
            outputs = sysio_get("mothur \"#set.logfile(name=master.logfile, append=T);" +
                                "align.seqs(fasta="+input.fasta+", reference="+silva+", flip=T, processors="+str(nprocessors)+")\"", [".align"])
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
        end = list(map(int, summ[:,2]))

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
