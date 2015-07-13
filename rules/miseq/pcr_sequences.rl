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
        keepdots = run["setup"]["miseq"]["keepdots"]
        silva = run["setup"]["silva"]

        sysio_set("mothur \"#set.logfile(name=master.logfile, append=T); pcr.seqs(fasta="+silva+", start="+pcrseqs_start+", end="+pcrseqs_end+", keepdots="+keepdots+", processors=8)\""[".fasta"], wildcards.project+".preprocess")                                            
         sysio.set("mothur \"#set.logfile(name=master.logfile, append=T); count.seqs(name="+input.names+", group="+input.groups+")\""[".count"], wildcards.project+".preprocess") 
