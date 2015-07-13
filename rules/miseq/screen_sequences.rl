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
        sysio_set("mothur \"#set.logfile(name=master.logfile, append=T); screen.seqs(fasta="+input.fasta+", group="+input.groups+", maxambig="+maxambig+", maxlength="+maxlength+")\""[".fasta",".group"], wildcards.project+".screen")  
        

