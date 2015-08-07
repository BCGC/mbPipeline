#snakemake rules for microbiome pipeline
#these rules will be seperated into individual files upon completion

import os
import json
import subprocess
import numpy
import warnings
from io import StringIO

with open('run.json') as data_file:
    run = json.load(data_file)
PROJECT = run["setup"]["proj"]
DATA_FILE_NAMES = run["setup"]["data"]

def sysio_set(cmd, extensions, newprefix):
    out=""
    with subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True, bufsize=1, universal_newlines=True) as p, StringIO() as buf:
        for line in p.stdout:
            print(line, end='')
            buf.write(line)
        out=buf.getvalue()
    outputs = {}
    for extension in extensions :
        current = out[out[0:out.rfind(extension)].rfind("\n")+1:out[out.rfind(extension):len(out)].find("\n")+out.rfind(extension)]
        wrong_extension = [".scrap"+extension, extension+".report"]
        if any(x in current for x in wrong_extension):
            out1 = out[0:out.rfind(current)+1]
            current = out1[out1[0:out1.rfind(extension)].rfind("\n")+1:out1[out1.rfind(extension):len(out1)].find("\n")+out1.rfind(extension)]
        new = newprefix + extension
        os.system("cp "+current+" "+new+"")
        outputs[extension] = new
    return outputs

def sysio_get(cmd, extensions):
    out=""
    with subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True, bufsize=1, universal_newlines=True) as p, StringIO() as buf:
        for line in p.stdout:
            print(line, end='')
            buf.write(line)
        out=buf.getvalue()
    outputs = {}
    for extension in extensions :
        current = out[out[0:out.rfind(extension)].rfind("\n")+1:out[out.rfind(extension):len(out)].find("\n")+out.rfind(extension)]
        wrong_extension = [".scrap"+extension, extension+".report"]
        if any(x in current for x in wrong_extension):
            out1 = out[0:out.rfind(current)]
            current = out1[out1[0:out1.rfind(extension)].rfind("\n")+1:out1[out1.rfind(extension):len(out1)].find("\n")+out1.rfind(extension)]
        outputs[extension] = current
    return outputs
