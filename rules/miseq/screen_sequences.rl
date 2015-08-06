rule screen_sequences:
    input: 
        fasta='{project}.fasta' 
        groups='{project}.groups'
    output:
        '{project}.screen.fasta'
        '{project}.screen.groups'
    run:
        outputs = sysio_set("mothur \"#set.logfile(name=master.logfile, append=T); screen.seqs(fasta="+input.fasta+", group="+input.groups+", maxambig=0, maxlength=275)\""[".fasta",".group"], wildcards.project+".screen")  
        

