rule beta_data:
    input: shared = '{project}.final.shared'
    output: '{project}.beta_data.out'
    run:
        outputs = sysio_get("mothur \"#summary.shared(shared="+input.shared+", calc=thetayc)\"", [".summary"])
        summary = outputs[".summary"]

        os.system("cut -f2 "+summary+" > .temp_sample1.out")
        num_lines5 = sum(1 for line in open('.temp_sample1.out'))
        sample1 = []
        f = open('.temp_sample1.out')
        for i in range(0, num_lines5):
            sample1.append(f.readline())
        f.close()
        for i in range(0, len(sample1)):
            sample1[i] = sample1[i][:-1]
        sample1[0] = "sample1"

        os.system("cut -f3 "+summary+" > .temp_sample2.out")
        sample2 = []
        f = open('.temp_sample2.out')
        for i in range(0, num_lines5):
            sample2.append(f.readline())
        f.close()
        for i in range(0, len(sample2)):
            sample2[i] = sample2[i][:-1]
        sample2[0] = "sample2"

        os.system("cut -f5 "+summary+" > .temp_bdiv.out")
        bdiv = []
        f = open('.temp_bdiv.out')
        for i in range(0, num_lines5):
            bdiv.append(f.readline())
        f.close()
        for i in range(0, len(bdiv)):
            bdiv[i] = bdiv[i][:-1]
        bdiv[0] = "bdiv"

        os.system("cut -f7 "+summary+" > .temp_cmin.out")
        cmin = []
        f = open('.temp_cmin.out')
        for i in range(0, num_lines5):
            cmin.append(f.readline())
        f.close()
        for i in range(0, len(cmin)):
            cmin[i] = cmin[i][:-1]
        for i in range(1, len(cmin)):
            cmin[i] = 1 - float(cmin[i])
        for i in range(1, len(cmin)):
            cmin[i] = str(cmin[i])
        cmin[0] = "cmin"

        os.system("cut -f6 "+summary+" > "".temp_cmax.out")
        cmax = []
        f = open('.temp_cmax.out')
        for i in range(0, num_lines5):
            cmax.append(f.readline())
        f.close()
        for i in range(0, len(cmax)):
            cmax[i] = cmax[i][:-1]
        for i in range(1, len(cmax)):
            cmax[i] = 1 - float(cmax[i])
        for i in range(1, len(cmax)):
            cmax[i] = str(cmax[i])
        cmax[0] = "cmax"

        with open(wildcards.project+'.beta_data.out', 'w') as f:
            for f1, f2, f3, f4, f5 in zip(sample1, sample2, bdiv, cmin, cmax):
                f.write(f1+"\t"+f2+"\t"+f3+"\t"+f4+"\t"+f5+"\n")
        f.close()
