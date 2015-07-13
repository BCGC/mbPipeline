
rule unique_sequences:
       input:
         fasta = '{project}.screen.fasta'
    output:
        '{project}.unique.fasta'
        '{project}.unique.names'
    run:
        sysio_set("mothur \"#set.logfile(name=master.logfile, append=T); unique.seqs(fasta="+input.fasta+")\""[".fasta",".names"], wildcards.project+".unique")
