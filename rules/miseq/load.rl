rule load_miseq:
    output: '{project}.fasta', '{project}.groups'
    run:
        with open('run.json') as data_file:
            run = json.load(data_file)
        nprocessors= run["setup"]["nprocessors"]
        maxlength= run["setup"]["miseq"]["files"]
      outputs = sysio.set("mothur \"#set.logfile(name=master.logfile, append=T);"  "make.contigs(file="+files+", processors="+str(nprocessors)+")\""[".fasta",".groups"], wildcards.project+)
  
