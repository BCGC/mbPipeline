rule remove_sequences:
    input:
        fasta='{project}.process.fasta',
        names='{project}.process.names',
        groups='{project}.process.groups'
    output:
        '{project}.remove.fasta',
        '{project}.remove.names',
        '{project}.remove.groups',
        '{project}.remove.taxonomy'
    run:
        with open('run.json') as data_file:
            run = json.load(data_file)
        nprocessors = run["setup"]["nprocessors"]
        trainset = run["setup"]["454"]["trainset"]
        
        # identify likely chimeras
        outputs = sysio_get("mothur \"#set.logfile(name=master.logfile, append=T);" +
                            "chimera.uchime(fasta="+input.fasta+", name="+input.names+", group="+input.groups+", processors="+str(nprocessors)+")\"", [".accnos"])

        accnos = outputs[".accnos"]
        tmp = numpy.genfromtxt(accnos, dtype='str')

        # remove identified chimeras, throwing exception if all sequences were flagged as chimeras
        fasta = input.fasta
        names = input.names
        groups = input.groups
        if tmp.shape[0] > 0:
            outputs = sysio_get("mothur \"#set.logfile(name=master.logfile, append=T);" +
                                "remove.seqs(accnos="+accnos+", fasta="+fasta+", name="+names+", " +
                                "group="+groups+", dups=T)\"", [".fasta",".names",".groups"])
            fasta = outputs[".fasta"]
            names = outputs[".names"]
            groups = outputs[".groups"]
        else:
            raise Exception("All sequences flagged as chimeras!")


        # classify sequences using given taxonomy trainset
        #os.system()
        outputs = sysio_get("mothur \"#set.logfile(name=master.logfile, append=T);" +
                            "classify.seqs(fasta="+fasta+", name="+names+", group="+groups+
                            ", template="+trainset+".fasta, taxonomy="+trainset+".tax, cutoff=80, processors="+str(nprocessors)+")\"", [".taxonomy"])


        taxonomy = outputs[".taxonomy"]

        # remove contaminant mitochondria/chloroplast sequences
        sysio_set("mothur \"#set.logfile(name=master.logfile, append=T);" + 
                "remove.lineage(fasta="+fasta+", name="+names+", group="+groups+", taxonomy="+taxonomy+
                  ", taxon=Mitochondria-Cyanobacteria_Chloroplast-unknown)\"", [".taxonomy",".names",".groups",".fasta"], wildcards.project+".remove")
