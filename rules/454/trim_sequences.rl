rule trim_sequences:
    input:
        fasta='{project}.load.fasta',
        names='{project}.load.names'
    output:
        '{project}.trim.fasta',
        '{project}.trim.names',
        '{project}.trim.groups'
    run:
        with open('run.json') as data_file:
            run = json.load(data_file)
        pdiffs = run["setup"]["454"]["pdiffs"]
        bdiffs = run["setup"]["454"]["bdiffs"]
        nprocessors = run["setup"]["nprocessors"]
        sysio_set("mothur \"#set.logfile(name=master.logfile, append=T); trim.seqs(fasta="+input.fasta+
                  ", name="+input.names+", oligos=oligos.txt, pdiffs="+pdiffs+", bdiffs="+bdiffs+
                  ", maxhomop=8, minlength=200, flip=T processors="+str(nprocessors)+")\"", [".fasta", ".names", ".groups"], wildcards.project+".trim")
