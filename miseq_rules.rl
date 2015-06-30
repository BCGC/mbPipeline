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

