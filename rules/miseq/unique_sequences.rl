
rule unique_sequences:
       input:
         fasta = '{project}.screen.fasta'
         fasta = '{project}.screen.groups'

    output:
        '{project}.preprocess.fasta',
        '{project}.preprocess.count'
    run:
        outputs = sysio_set("mothur \"#set.logfile(name=master.logfile, append=T); unique.seqs(fasta="+input.fasta+")\""[".fasta, .names"], wildcards.project+".preprocess")
fasta = outputs['.fasta']
names = outputs['.names']

outputs = sysio.set("mothur \"#set.logfile(name=master.logfile, append=T); count.seqs(name="+input.names+", group="+input.groups+")\""[".count"], wildcards.project+".preprocess") 
