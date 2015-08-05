rule unique_sequences:
    input:
        fasta='{project}.trim.fasta',
        names='{project}.trim.names'
    output:
        '{project}.unique.fasta',
        '{project}.unique.names'
    run:
        outputs = sysio_set("mothur \"#unique.seqs(fasta="+input.fasta+", name="+input.names+")\"", [".fasta",".names"], wildcards.project+".unique")
        fasta = outputs[".fasta"]
        names = outputs[".names"]

        p = subprocess.Popen("mothur \"#set.logfile(name=master.logfile, append=T); summary.seqs(fasta="+fasta+", name="+names+")\"", stdout=subprocess.PIPE, shell=True)
        out = p.communicate()[0].decode("utf-8")
        p.wait()
        out = out[out.find("97.5%-tile:")+12:len(out)]
        out = out[out.find("\t")+1:len(out)]
        out = out[out.find("\t")+1:len(out)]
        nbases = out[0:out.find("\t")]
        with open('run.json', 'r+') as f:
            run = json.load(f)
            run["storage"]["nbases"] = nbases
            f.seek(0)
            f.write(json.dumps(run))
            f.truncate()