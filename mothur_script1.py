#!/usr/local/bin/python2.7

# mothur_script1.py
# main driver for Mothur pipeline
# Kira Vasquez and Randall Johnson
# BSP CCR Genetics Core at Frederick National Laboratory
# Leidos Biomedical Research, Inc
# Created September-ish 2013
# Last Modified February 28, 2014

import os
import subprocess
import sys
import numpy
import warnings
import scipy


### set up arguments ###

args = dict()
for a in sys.argv[1:len(sys.argv)]:
      args[a[0:a.find('=')]] = a[a.find('=')+1:len(a)]

print('\nArgs:')
print(args)
print('\n')


### set up the files in my directory ###

pipeline = 'data/'
execfile(pipeline+'defaults.py') # where are these being kept? I should source them using absolute paths!

### process sff files ###

# if not os.path.isfile("all.flow.files"):
# get all sff files in the working current directory and convert to a character array
sff = subprocess.Popen('ls | grep sff', shell = True, stdout=subprocess.PIPE).communicate()[0]
sff = sff.rsplit('\n')

# this is the repository for all sff files
os.system("printf '' > all.flow.files")
os.system("printf '' > master.logfile")

for f in sff:
      if os.path.isfile(f):
            x = f[0:f.find('.sff')]
            os.system("mothur \"#set.logfile(name=master.logfile, append=T); " +
                                "sffinfo(sff="+x+".sff); " +
                                "summary.seqs(fasta="+x+".fasta); " +
                                "trim.flows(flow="+x+".flow, oligos=oligos.txt, pdiffs="+pdiffs+","+"bdiffs="+bdiffs+", processors=12)\"")
            os.system("cat "+x+".flow.files >> all.flow.files")


flows = 'all.flow.files'
os.system("mothur \"#set.logfile(name=master.logfile, append=T);" +
                    "shhh.flows(file="+flows+", processors=12, lookup="+DATAPATH+")\"")

fasta = 'all.shhh.fasta'
names = 'all.shhh.names'
groups = 'all.shhh.groups'


# check our sequences as of right now
# 0:seqname 1:start 2:end 3:nbases 4:ambigs 5:polymer 6:numSeqs
os.system("mothur \"#set.logfile(name=master.logfile, append=T); summary.seqs(fasta="+fasta+", name="+names+")\"")

summ = numpy.genfromtxt(fasta[0:fasta.find('fasta')] + 'summary', skiprows=1, dtype='str')

tmp = 0
for i in summ[:,3]:
    if int(i) < 200: # count number of reads that are less than 200 bp long
        tmp += 1

if tmp / summ.shape[0] > 0.2:
      warnings.warn(str(tmp / summ.shape[0] * 100) +
                    "% of unique reads are shorter than 200 bp.", Warning)


# trim barcodes and primers, make sure everything is xxx bp long
os.system("mothur \"#set.logfile(name=master.logfile, append=T); trim.seqs(fasta="+fasta+
          ", name="+names+", oligos=oligos.txt, pdiffs="+pdiffs+", bdiffs="+bdiffs+
          ", maxhomop=8, minlength=200, flip=T processors=12)\"")

fasta = fasta[0:fasta.find('fasta')] + 'trim.fasta'
names = names[0:names.find('names')] + 'trim.names'
os.system("mothur \"#set.logfile(name=master.logfile, append=T); summary.seqs(fasta="+fasta+", name="+names+")\"")


# save some effort by only analyzing unique sequences
os.system("mothur \"#unique.seqs(fasta="+fasta+", name="+names+")\"")

fasta = fasta[0:fasta.find('fasta')] + 'unique.fasta'
names = names[0:names.find('names')] + 'unique.names'
os.system("mothur \"#set.logfile(name=master.logfile, append=T); summary.seqs(fasta="+fasta+", name="+names+")\"")


# initial alignment
# oops...If you didn't get them flipped in the correct direction - use flip=T
os.system("mothur \"#set.logfile(name=master.logfile, append=T);" +
                    "align.seqs(fasta="+fasta+", reference=silva.bacteria.fasta, flip=T, processors=12)\"")

fasta = fasta[0:fasta.find('fasta')] + 'align'


# screen the sequences so we only keep the stuff in the region we are interested in :)
# 0:seqname 1:start 2:end 3:nbases 4:ambigs 5:polymer 6:numSeqs
os.system("mothur \"#set.logfile(name=master.logfile, append=T); summary.seqs(fasta="+fasta+", name="+names+")\"")

summ = numpy.genfromtxt(fasta[0:fasta.find('align')] + 'summary', skiprows=1, dtype='str')
end = map(int, summ[:,2])

if numpy.percentile(end, 25) != numpy.percentile(end, 75):
    warnings.warn("Sequence endings are not consistent. Check to see if they have been flipped.", Warning)
end = str(int(numpy.percentile(end, 50)))

os.system("mothur \"#set.logfile(name=master.logfile, append=T);" +
                    "screen.seqs(fasta="+fasta+", name="+names+", group="+groups+
                                 ", end="+end+", optimize=start, criteria=95, processors=12)\"")

fasta = fasta[0:fasta.find('align')] + 'good.align'
names = names[0:names.find('names')] + 'good.names'
groups = groups[0:groups.find('groups')] + 'good.groups'

os.system("mothur \"#set.logfile(name=master.logfile, append=T); summary.seqs(fasta="+fasta+", name="+names+")\"")


# filter the sequences so they all overlap the same region
os.system("mothur \"#set.logfile(name=master.logfile, append=T);" +
                    "filter.seqs(fasta="+fasta+", vertical=T, trump=., processors=12)\"")

fasta = fasta[0:fasta.find('align')] + 'filter.fasta' ####################################
print fasta

# should get some more unique sequences
os.system("mothur \"#set.logfile(name=master.logfile, append=T); unique.seqs(fasta="+fasta+", name="+names+")\"")

fasta = fasta[0:fasta.find('fasta')] + 'unique.fasta'
print fasta
names = names[0:names.find('names')] + 'filter.names'

os.system("mothur \"#set.logfile(name=master.logfile, append=T); summary.seqs(fasta="+fasta+", name="+names+")\"")


# precluster to help get rid of sequencing errors - also helps with computational efficiency
os.system("mothur \"#set.logfile(name=master.logfile, append=T);" +
                    "pre.cluster(fasta="+fasta+", name="+names+", group="+groups+", diffs=2)\"")

fasta = fasta[0:fasta.find('fasta')] + 'precluster.fasta'
print fasta
names = names[0:names.find('names')] + 'unique.precluster.names'

os.system("mothur \"#set.logfile(name=master.logfile, append=T); summary.seqs(fasta="+fasta+", name="+names+")\"")


# identify likely chimeras
os.system("mothur \"#set.logfile(name=master.logfile, append=T);" +
                    "chimera.uchime(fasta="+fasta+", name="+names+", group="+groups+", processors=12)\"")

accnos = fasta[0:fasta.find('fasta')] + 'uchime.accnos'
tmp = numpy.genfromtxt(accnos, dtype='str')

# remove identified chimeras, throwing exception if all sequences were flagged as chimeras
if tmp.shape[0] > 0:
    os.system("mothur \"#set.logfile(name=master.logfile, append=T);" +
                        "remove.seqs(accnos="+accnos+", fasta="+fasta+", name="+names+", " +
                                    "group="+groups+", dups=T)\"")
else:
    raise Exception("All sequences flagged as chimeras!")

################# NIKHIL #################

fasta = fasta[0:fasta.find('fasta')] + 'pick.fasta'
print fasta
names = names[0:names.find('names')] + 'pick.names'
groups = groups[0:groups.find('groups')] + 'pick.groups'

# classify sequences using given taxonomy trainset
os.system("mothur \"#set.logfile(name=master.logfile, append=T);" +
          "classify.seqs(fasta="+fasta+", name="+names+", group="+groups+
          ", template=trainset7_112011.pds.fasta, taxonomy=trainset7_112011.pds.tax, cutoff=80, processors=12)\"")

taxonomy = fasta[0:fasta.find('fasta')] + 'pds.taxonomy'
accnos = fasta[0:fasta.find('fasta')] + 'pds.flip.accnos'

# remove contaminant mitochondria/chloroplast sequences
os.system("mothur \"#set.logfile(name=master.logfile, append=T);" + 
          "remove.lineage(fasta="+fasta+", name="+names+", group="+groups+", taxonomy="+taxonomy+
          ", taxon=Mitochondria-Cyanobacteria_Chloroplast-unknown)\"")

taxonomy = taxonomy[0:taxonomy.find('taxonomy')] + 'pick.taxonomy'
names = names[0:names.find('names')] + 'pick.names'
fasta = fasta[0:fasta.find('fasta')] + 'pick.fasta'
groups = groups[0:groups.find('groups')] + 'pick.groups'

# summary??

# final files
os.system("cp "+fasta+" final.fasta")
fasta = 'final.fasta'
os.system("cp "+names+" final.names")
names = 'final.names'
os.system("cp "+groups+" final.groups")
groups = 'final.groups'
os.system("cp "+taxonomy+" final.taxonomy")
taxonomy = 'final.taxonomy'

### get sequence data ###

os.system("mothur \"#set.logfile(name=master.logfile, append=T);" + 
          "count.groups(group=final.groups)\" > .seq_data.out")

### pull apart data in x.seq_data.out ###

num_lines = sum(1 for line in open('.seq_data.out'))
data = []
f = open('.seq_data.out')
for i in range(0, num_lines) :
      text = f.readline()
      if 'contains' in text:
      	   data.append(text)
f.close()
locs = []
nums = []
for i in range(0, len(data)):
      data[i] = data[i][:-2]
for i in range(0, len(data)):
      temp1,_,temp2 = data[i].partition(' contains ')
      locs.append(temp1)
      nums.append(temp2)

### print warnings, find optimal sequence size and save ctrl seqs to file ###

are_controls = raw_input("Do you have controls? Enter 1 for 'yes' or 2 for 'no': ")
are_controls = int(are_controls)
if are_controls == 1:
      ctrls = []
      num_lines2 = sum(1 for line in open('.control.samples'))
      f = open('.control.samples')
      for i in range(0, num_lines2):
      	  ctrls.append(f.readline())
      f.close()
      for i in range(0, len(ctrls)):
      	  ctrls[i] = ctrls[i][:-1]
      ctrl_nums = []
      ctrl_warn = []
      ctrl_locs = []
      for i in range(0, len(ctrls)):
      	  for j in range(0, len(locs)-1):
      	      if ctrls[i] == locs[j]:
	      	 ctrl_locs.append(locs.pop(j))
	      	 ctrl_nums.append(nums.pop(j))
      for i in range(0, len(ctrl_nums)):
      	  if float(ctrl_nums[i]) > 1000:
      	     ctrl_warn.append(ctrl_locs[i])

      f = open('.control.seqs', 'w')
      for i in range(0, len(ctrls)):
      	  f.write(ctrls[i] + ": " + ctrl_nums[i] + " \n")
      f.close()

      print ""
      print "Warning: the following control samples have an unusually high number of sequences: " + str(ctrl_warn)


f = open('.temp.numseqs', 'w')
for i in range(0, len(nums)):
      f.write(str(nums[i]) + " \n")
f.close()

f = open('.temp.locs', 'w')
for i in range(0, len(locs)):
      f.write(str(locs[i]) + " \n")
f.close()

low_warn = []
for i in range(0, len(nums)):
      if float(nums[i]) < 3000:
      	   low_warn.append(locs[i])
print ""
print "Warning: the following samples have an unusually low number of sequences: " + str(low_warn)

### user may choose to keep low-sequence samples ###

low_seq_nums = []
for i in range(0, len(low_warn)):
      for j in range(0, len(nums)-1):
      	   if locs[j] == low_warn[i]:
	      low_seq_nums.append(nums[j])
print ""
for i in range(0, len(low_warn)):
      print low_warn[i] + " has " + low_seq_nums[i] + " sequences."


for i in range(0, len(low_warn)):
      for j in range(0, len(nums)-1):
      	   if locs[j] == low_warn[i]:
	      locs.pop(j)
	      nums.pop(j)
highest = 0
for i in range(0, len(nums)):
      if nums[i] > highest:
      	   highest = nums[i]
lowest = highest

for i in range(0, len(nums)):
      if nums[i] < lowest:
      	   lowest = nums[i]
	   ideal_loc = locs[i]
print ""
lowest = raw_input("We recommend that the lowest number of sequences should be " + lowest + " from " + ideal_loc + ". What would you like to set the lowest allowed number of sequences to? ")

### remove controls ###

if are_controls == 1:
      os.system("mothur \"#set.logfile(name=master.logfile, append=T);" +
                "remove.groups(fasta="+fasta+", accnos="+x+".control.samples, group="+groups+
                ", name="+names+".final.names, taxonomy="+taxonomy+")\"")
      fasta = fasta[0:fasta.find('fasta')] + 'pick.fasta'
      taxonomy = taxonomy[0:taxonomy.find('taxonomy')] + 'pick.taxonomy'
      names = names[0:names.find('names')] + 'pick.names'
      groups = groups[0:groups.find('groups')] + 'pick.groups'

### OTUs ###

os.system("mothur \"#set.logfile(name=master.logfile, append=T);" +
          "dist.seqs(fasta="+fasta+", cutoff=0.15, processors=12)\"")

dist = fasta[0:fasta.find('fasta')] + 'dist'

os.system("mothur \"#set.logfile(name=master.logfile, append=T);" +
          "cluster(column="+dist+", name="+names+")\"")

list = fasta[0:fasta.find('fasta')] + 'an.list'

os.system("mothur \"#set.logfile(name=master.logfile, append=T);" +
          "make.shared(list="+list+", group="+groups+", label=0.03)\"")

shared = list[0:list.find('list')] + 'shared'

os.system("mothur \"#set.logfile(name=master.logfile, append=T);" +
          "sub.sample(shared="+shared+", size="+lowest+")\"")

sharedold = shared #FIGURE OUT WHATS HAPPENING HERE - THIS IS BAD NOMENCLATURE - but works for now ;)
shared = list[0:shared.find('shared')] + '0.03.subsample.shared'

os.system("mothur \"#set.logfile(name=master.logfile, append=T);" +
          "classify.otu(list="+list+", name="+names+", taxonomy="+taxonomy+", label=0.03)\"")

txconsensus = taxonomy[0:taxonomy.find('taxonomy')] + 'an.0.03.cons.taxonomy'

os.system("mothur \"#set.logfile(name=master.logfile, append=T);" +
          "phylotype(taxonomy="+taxonomy+", name="+names+", label=1)\"")

txlist = fasta[0:fasta.find('fasta')] + 'tx.list'

os.system("mothur \"#set.logfile(name=master.logfile, append=T);" +
          "make.shared(list="+txlist+", group="+groups+", label=1)\"")

txshared = txlist[0:txlist.find('list')] + 'shared'

os.system("mothur \"#set.logfile(name=master.logfile, append=T);" +
          "sub.sample(shared="+txshared+", size="+lowest+")\"")

os.system("mothur \"#set.logfile(name=master.logfile, append=T);" +
          "classify.otu(list="+txlist+", name="+names+", taxonomy="+taxonomy+", label=1)\"")
txconsensus = taxonomy[0:taxonomy.find('taxonomy')] + 'tx.1.cons.taxonomy'

### Alpha Diversity ###

os.system("mothur \"#set.logfile(name=master.logfile, append=T);" +
          "collect.single(shared="+shared+", calc=chao-invsimpson, freq=100)\"")

sample_list = []
os.system("grep -l '0.03' *.invsimpson > .sample_list.out")
num_lines3 = sum(1 for line in open('.sample_list.out'))
f = open('.sample_list.out')
for i in range(0, num_lines3):
      sample_list.append(f.readline())
      sample_list[i] = sample_list[i][:-1]
f.close()
temp1 = []
summ = 0
invsimpson = []
for i in range(0, num_lines3):
      os.system("cut -f2 -s "+sample_list[i]+" | tail -n 5 > .temp_nums.out")
      num_lines4 = sum(1 for line in open('.temp_nums.out'))
      f = open('.temp_nums.out')
      for j in range(0, num_lines4):
      	  temp1.append(f.readline())
      for z in range(0, num_lines4):
      	  summ += float(temp1[z])
      temp1 = []
      invsimpson.append(summ/num_lines4)
      summ = 0
      f.close()
f = open('.temp.adiv', 'w')
for i in range(0, len(invsimpson)):
      f.write(str(invsimpson[i]) + ' \n')
f.close()

### Generating Graphics Data File ###
#NEED TO DEVELOP A WAY TO HANDLE METADATA - FOR NOW MANUAL INPUT
#seqs = ["meta", "nseqs"]
#adiv = ["meta", "adiv"]
#barcode = ["meta", "Barcode"]
#variables = []
#num_lines = sum(1 for line in open('.temp.numseqs'))
#print "You must enter at least one set of independent categorical or continuous variables that describe each sample in order to generate plots!"
#cont = "1"
#while cont == "1":
#      newvar = raw_input("Enter the name describing the first variable (eg. gender, age, etc.): ")
#      newvarlist = []
#      success = False
#      while not success:
#            type = raw_input("Enter the type of variable that it is, cat for catergorical or cont for continuous (eg. gender is cat, age is cont): ")
#            if "cat" in type:
#                  newvarlist.append('cat')
#                  success = True
#            if "cont" in type:
#                  newvarlist.append('cont')
#                  success = True
#      newvarlist.append(newvar)
#      f = open('.temp.locs')
#      for i in range(0, num_lines) :
#            barcode = f.readline()
#            value = raw_input("Enter value of " +newvar+ " describing " +barcode+ "(be sure to be consistent!) : ")
#            newvarlist.append(value)
#      f.close()
#      variables.append(newvarlist)
#      print ""
#      print "Entry for variable completed."
#      print ""
#      cont = raw_input("Are there more variables to define and enter? Enter 1 for yes or 2 for no: ")
#
#f = open('.temp.numseqs')
#for i in range(0, num_lines) :
#    seqs.append(f.readline())
#f.close()
#
#f = open('.temp.adiv')
#for i in range(0, num_lines) :
#    adiv.append(f.readline())
#f.close()
#
#f = open('.temp.locs')
#for i in range(0, num_lines) :
#    barcode.append(f.readline())
#f.close()

#for i in range(2, num_lines+2) :
#    barcode[i] = barcode[i][:-2]
#    adiv[i] = adiv[i][:-2]
#    seqs[i] = seqs[i][:-2]
#
#f = open('graphics_data.txt', 'w')
#for i in range(0, num_lines+2):
#      f.write(barcode[i]+"\t"+seqs[i]+"\t"+adiv[i]+"\t")
#      for j in range(0, len(variables)):
#            f.write(variables[j][i]+"\t")
#      f.write("\n")
#f.close()

### Beta Diversity ###

os.system("mothur \"#summary.shared(shared="+sharedold+", calc=thetayc)\"")

summary = sharedold + '.summary'

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

with open('beta_data.out', 'w') as f:
      for f1, f2, f3, f4, f5 in zip(sample1, sample2, bdiv, cmin, cmax):
            f.write(f1+"\t"+f2+"\t"+f3+"\t"+f4+"\t"+f5+"\n")
f.close()


### USING mbGRAPHCIS R PACKAGE TO PRODUCE GRAPHS ###
seqs = ["meta", "nseqs"]
adiv = ["meta", "adiv"]
num_lines = sum(1 for line in open('.temp.numseqs'))

f = open('.temp.numseqs')
for i in range(0, num_lines) :
    seqs.append(f.readline())
f.close()

f = open('.temp.adiv')
for i in range(0, num_lines) :
    adiv.append(f.readline())
f.close()

for i in range(2, num_lines+2) :
    barcode[i] = barcode[i][:-2]
    adiv[i] = adiv[i][:-2]
    seqs[i] = seqs[i][:-2]

num_lines = sum(1 for line in open(metadata))
f1 = open(metadata)
lines = f1.readlines()
f2 = open("mb_graphics_data.txt", "w")
for i in range(0, num_lines) :
      tabs = lines[i].split("\t")
      tabs[len(tabs)-1] = tabs[len(tabs)-1][0:tabs[len(tabs)-1].find('\n')]
      tabs.append(seqs[i])
      tabs.append(adiv[i])
      f2.write("\t".join(tabs)+"\n")
f1.close()
f2.close()

os.system("Rscript graphall.R "+txconsensus+" "+txshared+" 0.14")




#################################### IGNORE VVVVVVVVV



#os.system("mothur \"#summary.seqs(fasta="+x+".shhh.trim.unique.good.filter.unique.precluster.pick.fasta, name="+x+".shhh.trim.unique.good.filter.unique.precluster.pick.names)\"")

#os.system("mothur \"#classify.seqs(fasta="+x+".shhh.trim.unique.good.filter.unique.precluster.pick.fasta, name="+x+".shhh.trim.unique.good.filter.unique.precluster.pick.names, group="+x+".shhh.good.pick.groups, template=trainset7_112011.pds.fasta, taxonomy=trainset7_112011.pds.tax, cutoff=80, processors=12)\"")

#os.system("mothur \"#remove.lineage(fasta="+x+".shhh.trim.unique.good.filter.unique.precluster.pick.fasta, name="+x+".shhh.trim.unique.good.filter.unique.precluster.pick.names, group="+x+".shhh.good.pick.groups, taxonomy="+x+".shhh.trim.unique.good.filter.unique.precluster.pick.pds.taxonomy, taxon=Mitochondria-Cyanobacteria_Chloroplast-unknown)\"")
#os.system("mothur \"#summary.seqs(fasta="+x+".shhh.trim.unique.good.filter.unique.precluster.pick.pick.fasta, name="+x+".shhh.trim.unique.good.filter.unique.precluster.pick.pick.names)\"")

### shortening to final file names ###

#os.system("mothur \"#system(cp "+x+".shhh.trim.unique.good.filter.unique.precluster.pick.pick.fasta "+x+".final.fasta)\"")
#os.system("mothur \"#system(cp "+x+".shhh.trim.unique.good.filter.unique.precluster.pick.pick.names "+x+".final.names)\"")
#os.system("mothur \"#system(cp "+x+".shhh.good.pick.pick.groups "+x+".final.groups)\"")
#os.system("mothur \"#system(cp "+x+".shhh.trim.unique.good.filter.unique.precluster.pick.pds.pick.taxonomy "+x+".final.taxonomy)\"")

# ### get sequence data ###

# #os.system("mothur \"#count.groups(group="+x+".final.groups)\" > "+x+".seq_data.out")

# ### pull apart data in x.seq_data.out ###

# num_lines = sum(1 for line in open(''+x+'.seq_data.out'))
# data = []
# f = open(''+x+'.seq_data.out')
# for i in range(0, num_lines-2) :
#       if i > 28:
#       	   data.append(f.readline())
#       else:
# 	   f.readline()
# f.close()
# locs = []
# nums = []
# for i in range(0, len(data)):
#       data[i] = data[i][:-2]
# for i in range(0, len(data)):
#       temp1,_,temp2 = data[i].partition(' contains ')
#       locs.append(temp1)
#       nums.append(temp2)

# ### print warnings, find optimal sequence size and save ctrl seqs to file ###

# are_controls = raw_input("Do you have controls? Enter 1 for 'yes' or 2 for 'no': ")
# are_controls = int(are_controls)
# if are_controls == 1:
#       ctrls = []
#       num_lines2 = sum(1 for line in open(''+x+'.control.samples'))
#       f = open(''+x+'.control.samples')
#       for i in range(0, num_lines2):
#       	  ctrls.append(f.readline())
#       f.close()
#       for i in range(0, len(ctrls)):
#       	  ctrls[i] = ctrls[i][:-1]
#       ctrl_nums = []
#       ctrl_warn = []
#       ctrl_locs = []
#       for i in range(0, len(ctrls)):
#       	  for j in range(0, len(locs)-1):
#       	      if ctrls[i] == locs[j]:
# 	      	 ctrl_locs.append(locs.pop(j))
# 	      	 ctrl_nums.append(nums.pop(j))
#       for i in range(0, len(ctrl_nums)):
#       	  if float(ctrl_nums[i]) > 1000:
#       	     ctrl_warn.append(ctrl_locs[i])

#       f = open(''+x+'.control.seqs', 'w')
#       for i in range(0, len(ctrls)):
#       	  f.write(ctrls[i] + ": " + ctrl_nums[i] + " \n")
#       f.close()

#       print ""
#       print "Warning: the following control samples have an unusually high number of sequences: " + str(ctrl_warn)


# f = open(''+x+'.temp.numseqs', 'w')
# for i in range(0, len(nums)):
#       f.write(str(nums[i]) + " \n")
# f.close()


# low_warn = []
# for i in range(0, len(nums)):
#       if float(nums[i]) < 3000:
#       	   low_warn.append(locs[i])
# print ""
# print "Warning: the following samples have an unusually low number of sequences: " + str(low_warn)

# ### user may choose to keep low-sequence samples ###

# low_seq_nums = []
# for i in range(0, len(low_warn)):
#       for j in range(0, len(nums)-1):
#       	   if locs[j] == low_warn[i]:
# 	      low_seq_nums.append(nums[j])
# print ""
# for i in range(0, len(low_warn)):
#       print low_warn[i] + " has " + low_seq_nums[i] + " sequences."


# for i in range(0, len(low_warn)):
#       for j in range(0, len(nums)-1):
#       	   if locs[j] == low_warn[i]:
# 	      locs.pop(j)
# 	      nums.pop(j)
# highest = 0
# for i in range(0, len(nums)):
#       if nums[i] > highest:
#       	   highest = nums[i]
# lowest = highest

# for i in range(0, len(nums)):
#       if nums[i] < lowest:
#       	   lowest = nums[i]
# 	   ideal_loc = locs[i]
# print ""
# lowest = raw_input("We recommend that the lowest number of sequences should be " + lowest + " from " + ideal_loc + ". What would you like to set the lowest allowed number of sequences to? ")

# ### remove controls ###

# if are_controls == 1:
#       os.system("mothur \"#remove.groups(fasta="+x+".final.fasta, accnos="+x+".control.samples, group="+x+".final.groups, name="+x+".final.names, taxonomy="+x+".final.taxonomy)\"")

# ### OTU section ###

# #if are_controls == 1:
#       #os.system("mothur \"#dist.seqs(fasta="+x+".final.pick.fasta, cutoff=0.15, processors=12)\"")

#       #os.system("mothur \"#cluster(column="+x+".final.pick.dist, name="+x+".final.pick.names)\"")

#       #os.system("mothur \"#make.shared(list="+x+".final.pick.an.list, group="+x+".final.pick.groups, label=0.03)\"")

#       #os.system("mothur \"#sub.sample(shared="+x+".final.pick.an.shared, size="+lowest+")\"")

#       #os.system("mothur \"#classify.otu(list="+x+".final.pick.an.list, name="+x+".final.pick.names, taxonomy="+x+".final.pick.taxonomy, label=0.03)\"")

#       #os.system("mothur \"#phylotype(taxonomy="+x+".final.pick.taxonomy, name="+x+".final.pick.names, label=1)\"")

#       #os.system("mothur \"#make.shared(list="+x+".final.pick.tx.list, group="+x+".final.pick.groups, label=1)\"")

#       #os.system("mothur \"#sub.sample(shared="+x+".final.pick.tx.shared, size="+lowest+")\"")

#       #os.system("mothur \"#classify.otu(list="+x+".final.pick.tx.list, name="+x+".final.pick.names, taxonomy="+x+".final.pick.taxonomy, label=1)\"")

# #if are_controls != 1:
#       #os.system("mothur \"#dist.seqs(fasta="+x+".final.fasta, cutoff=0.15, processors=12)\"")

#       #os.system("mothur \"#cluster(column="+x+".final.dist, name="+x+".final.names)\"")

#       #os.system("mothur \"#make.shared(list="+x+".final.an.list, group="+x+".final.groups, label=0.03)\"")

#       #os.system("mothur \"#sub.sample(shared="+x+".final.an.shared, size="+lowest+")\"")

#       #os.system("mothur \"#classify.otu(list="+x+".final.an.list, name="+x+".final.names, taxonomy="+x+".final.taxonomy, label=0.03)\"")

#       #os.system("mothur \"#phylotype(taxonomy="+x+".final.taxonomy, name="+x+".final.names, label=1)\"")

#       #os.system("mothur \"#make.shared(list="+x+".final.tx.list, group="+x+".final.groups, label=1)\"")

#       #os.system("mothur \"#sub.sample(shared="+x+".final.tx.shared, size="+lowest+")\"")

#       #os.system("mothur \"#classify.otu(list="+x+".final.tx.list, name="+x+".final.names, taxonomy="+x+".final.taxonomy, label=1)\"")

# ### alpha diversity ###

# #if are_controls == 1:
#       #os.system("mothur \"#collect.single(shared="+x+".final.pick.an.0.03.subsample.shared, calc=chao-invsimpson, freq=100)\"")
# #if are_controls != 1:
#       #os.system("mothur \"#collect.single(shared="+x+".final.an.0.03.subsample.shared, calc=chao-invsimpson, freq=100)\"")

# sample_list = []
# os.system("grep -l '0.03' "+x+"*.invsimpson > "+x+".sample_list.out")
# num_lines3 = sum(1 for line in open(''+x+'.sample_list.out'))
# f = open(''+x+'.sample_list.out')
# for i in range(0, num_lines3):
#       sample_list.append(f.readline())
#       sample_list[i] = sample_list[i][:-1]
# f.close()
# temp1 = []
# summ = 0
# invsimpson = []
# for i in range(0, num_lines3):
#       os.system("cut -f2 -s "+sample_list[i]+" | tail -n 5 > "+x+".temp_nums.out")
#       num_lines4 = sum(1 for line in open(''+x+'.temp_nums.out'))
#       f = open(''+x+'.temp_nums.out')
#       for j in range(0, num_lines4):
#       	  temp1.append(f.readline())
#       for z in range(0, num_lines4):
#       	  summ += float(temp1[z])
#       temp1 = []
#       invsimpson.append(summ/num_lines4)
#       summ = 0
#       f.close()
# f = open(''+x+'.temp.adiv', 'w')
# for i in range(0, len(invsimpson)):
#       f.write(str(invsimpson[i]) + " \n")
# f.close()

# ### beta diversity ###

# if are_controls == 1:
#       os.system("mothur \"#summary.shared(shared="+x+".final.pick.an.shared, calc=thetayc)\"")

#       os.system("cut -f2 "+x+".final.pick.an.shared.summary > "+x+".temp_sample1.out")
#       num_lines5 = sum(1 for line in open(''+x+'.temp_sample1.out'))
#       sample1 = []
#       f = open(''+x+'.temp_sample1.out')
#       for i in range(0, num_lines5):
#       	  sample1.append(f.readline())
#       f.close()
#       for i in range(0, len(sample1)):
#       	  sample1[i] = sample1[i][:-1]
#       sample1[0] = "sample1"

#       os.system("cut -f3 "+x+".final.pick.an.shared.summary > "+x+".temp_sample2.out")
#       sample2 = []
#       f = open(''+x+'.temp_sample2.out')
#       for i in range(0, num_lines5):
#       	  sample2.append(f.readline())
#       f.close()
#       for i in range(0, len(sample2)):
#       	  sample2[i] = sample2[i][:-1]
#       sample2[0] = "sample2"

#       os.system("cut -f5 "+x+".final.pick.an.shared.summary > "+x+".temp_bdiv.out")
#       bdiv = []
#       f = open(''+x+'.temp_bdiv.out')
#       for i in range(0, num_lines5):
#       	  bdiv.append(f.readline())
#       f.close()
#       for i in range(0, len(bdiv)):
#       	  bdiv[i] = bdiv[i][:-1]
#       bdiv[0] = "bdiv"

#       os.system("cut -f6 "+x+".final.pick.an.shared.summary > "+x+".temp_cmin.out")
#       cmin = []
#       f = open(''+x+'.temp_cmin.out')
#       for i in range(0, num_lines5):
#       	  cmin.append(f.readline())
#       f.close()
#       for i in range(0, len(cmin)):
#       	  cmin[i] = cmin[i][:-1]
#       for i in range(1, len(cmin)):
#       	  cmin[i] = 1 - float(cmin[i])
#       for i in range(1, len(cmin)):
#       	  cmin[i] = str(cmin[i])
#       cmin[0] = "cmin"

#       os.system("cut -f7 "+x+".final.pick.an.shared.summary > "+x+".temp_cmax.out")
#       cmax = []
#       f = open(''+x+'.temp_cmax.out')
#       for i in range(0, num_lines5):
#       	  cmax.append(f.readline())
#       f.close()
#       for i in range(0, len(cmax)):
#       	  cmax[i] = cmax[i][:-1]
#       for i in range(1, len(cmax)):
#       	  cmax[i] = 1 - float(cmax[i])
#       for i in range(1, len(cmax)):
#       	  cmax[i] = str(cmax[i])
#       cmax[0] = "cmax"
# if are_controls != 1:
#       os.system("mothur \"#summary.shared(shared="+x+".final.an.shared, calc=thetayc)\"")

#       os.system("cut -f2 "+x+".final.an.shared.summary > "+x+".temp_sample1.out")
#       num_lines5 = sum(1 for line in open(''+x+'.temp_sample1.out'))
#       sample1 = []
#       f = open(''+x+'.temp_sample1.out')
#       for i in range(0, num_lines5):
#       	  sample1.append(f.readline())
#       f.close()
#       for i in range(0, len(sample1)):
#       	  sample1[i] = sample1[i][:-1]
#       sample1[0] = "sample1"

#       os.system("cut -f3 "+x+".final.an.shared.summary > "+x+".temp_sample2.out")
#       sample2 = []
#       f = open(''+x+'.temp_sample2.out')
#       for i in range(0, num_lines5):
#       	  sample2.append(f.readline())
#       f.close()
#       for i in range(0, len(sample2)):
#       	  sample2[i] = sample2[i][:-1]
#       sample2[0] = "sample2"

#       os.system("cut -f5 "+x+".final.an.shared.summary > "+x+".temp_bdiv.out")
#       bdiv = []
#       f = open(''+x+'.temp_bdiv.out')
#       for i in range(0, num_lines5):
#       	  bdiv.append(f.readline())
#       f.close()
#       for i in range(0, len(bdiv)):
#       	  bdiv[i] = bdiv[i][:-1]
#       bdiv[0] = "bdiv"

#       os.system("cut -f6 "+x+".final.an.shared.summary > "+x+".temp_cmin.out")
#       cmin = []
#       f = open(''+x+'.temp_cmin.out')
#       for i in range(0, num_lines5):
#       	  cmin.append(f.readline())
#       f.close()
#       for i in range(0, len(cmin)):
#       	  cmin[i] = cmin[i][:-1]
#       for i in range(1, len(cmin)):
#       	  cmin[i] = 1 - float(cmin[i])
#       for i in range(1, len(cmin)):
#       	  cmin[i] = str(cmin[i])
#       cmin[0] = "cmin"

#       os.system("cut -f7 "+x+".final.an.shared.summary > "+x+".temp_cmax.out")
#       cmax = []
#       f = open(''+x+'.temp_cmax.out')
#       for i in range(0, num_lines5):
#       	  cmax.append(f.readline())
#       f.close()
#       for i in range(0, len(cmax)):
#       	  cmax[i] = cmax[i][:-1]
#       for i in range(1, len(cmax)):
#       	  cmax[i] = 1 - float(cmax[i])
#       for i in range(1, len(cmax)):
#       	  cmax[i] = str(cmax[i])
#       cmax[0] = "cmax"

# with open(''+x+'.beta_data.out', 'w') as f:
#       for f1, f2, f3, f4, f5 in zip(sample1, sample2, bdiv, cmin, cmax):
#       	  f.write(f1+"\t"+f2+"\t"+f3+"\t"+f4+"\t"+f5+"\n")
# f.close()
