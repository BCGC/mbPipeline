rule data_setup:
    input: temp_adiv='{project}.temp.adiv', temp_locs='{project}.temp.locs', temp_nums='{project}.temp.nums'
    output: '{project}.mb_graphics_data.txt',
    run:
        seqs = ["meta", "nseqs"]
        adiv = ["meta", "adiv"]
        barcode = ["meta", "Barcode"]
        with open('run.json') as data_file:
            run = json.load(data_file)
        num_lines = int(run["storage"]["lines"])
        metadata = run["setup"]["metadata"]
        indvars = run["setup"]["indvars"]

        f = open(input.temp_adiv)
        for i in range(0, num_lines) :
            seqs.append(f.readline())
        f.close()

        f = open(input.temp_adiv)
        for i in range(0, num_lines) :
            adiv.append(f.readline())
        f.close()

        f = open(input.temp_locs)
        for i in range(0, num_lines) :
            barcode.append(f.readline())
        f.close()

        for i in range(2, num_lines+2) :
            barcode[i] = barcode[i][:-2]
            adiv[i] = adiv[i][:-2]
            seqs[i] = seqs[i][:-2]

        num_lines = sum(1 for line in open(metadata))
        f1 = open(metadata)
        lines = f1.readlines()
        f2 = open("final_data.txt", "w")
        #This for loop is terribly overcoded - but hey, it works ;)
        for i in range(0, num_lines) :
            tabs = lines[i].split("\t")
            tabs[len(tabs)-1] = tabs[len(tabs)-1][0:tabs[len(tabs)-1].find('\n')]
            if i==0:
                tabs.append(seqs[i])
                tabs.append(adiv[i])
                f2.write("\t".join(tabs)+"\n")
            if i==1:
                tabs.append(seqs[i])
                tabs.append(adiv[i])
                f2.write("\t".join(tabs)+"\n")
            if i>=2:
                for j in range(2, len(barcode)) :
                    if barcode[j] in tabs: #only continues if barcode is found
                        tabs.append(seqs[j])
                        tabs.append(adiv[j])
                        f2.write("\t".join(tabs)+"\n")
        f1.close()
        f2.close()

        if not len(indvars) == 0 :
            f1 = open("final_data.txt")
            f2 = open("mb_graphics_data.txt", "w")
            lines = f1.readlines()
            numcols = len(lines[0].split("\t"))
            columns_to_ignore = []
            for i in range(0, numcols) :
                if lines[0].split("\t")[i] == "cat" or lines[0].split("\t")[i] == "cont" :
                    if not lines[1].split("\t")[i] in indvars :
                        columns_to_ignore.append(i)
            num_lines=len(lines)
            for i in range(0, num_lines) :
                tabs = lines[i].split("\t")
                tabs[len(tabs)-1] = tabs[len(tabs)-1][0:tabs[len(tabs)-1].find('\n')]
                tabs = [j for k, j in enumerate(tabs) if k not in columns_to_ignore]
                f2.write("\t".join(tabs)+"\n")
            f1.close()
            f2.close()
        else:
            import shutil 
            shutil.copy2("final_data.txt", wildcards.project+"mb_graphics_data.txt")
