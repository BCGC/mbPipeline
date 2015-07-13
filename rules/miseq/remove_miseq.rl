
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
        

