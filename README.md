mbPipeline
==========

**Developed by:**
Nikhil Gowda, Randall Johnson, Kira Vasquez, Dominique Brown

==========
mbPipeline is an automated workflow to process microbiome sequence data from both 454 and MiSeq sequencing (MiSeq pipeline not supported yet). Raw output files are generated for the users perusal, but the key feature of this pipeline is the automated graphics produced utilizing user specified metadata.

To use this pipeline please clone it to your repository. You must have Python 3.X, Snakemake, and R installed. In order to generate graphics you must have the package mbGraphics (available at: www.github.com/BCGC/mbGraphics) installed on R. Numpy must be installed on python as well.

The pipeline is launched by calling the launcher script. Here is an example call:
```
python /Users/gowdanb/Desktop/mbPipeline/Unified/launcher.py pipeline=454 metadata=metadata.txt refpath=../reference datapath=../data trainset=trainset7_112011.pds processors=2 project=test1
```


==========
**Arguments:**

```pipeline ```
This specifies the pipeline to run (currently supported options are 454 and miseq).

```refpath, datapath, and workdir ```
This specifies the location of the directories containing reference files, sff files, and the workding directory respectively. If not specified each will default to the current working directory (note specifying workdir will update the working directory before refpath and datapath are set).
Refpath must contain:
metadata.txt, oligos.txt
Refpath can contain the following files:
454: trainset (.fasta and .tax), lookup file, silva bacterial reference file (silva.bacteria.fasta), controls file
Each of these optional reference files will be substituted with default reference files (except the controls file of course)

```trainset ```
This specifies the name of the trainset file (if provided in the reference directory).

```controlsfile ```
This specifies the name of the controlsfile  (if provided in the reference directory).

```silva ```
This specifies the name of the silva file (if provided in the reference directory)

```project ```
This is required and all generated data files will be prefixed by this.

```processors ```
This is the current method of specifying how many processors to use (if running on biowulf the job must be submitted with equal or more processors, and if not specified the launcher will attempt to determine the number of available processors).

```new ```
If a run.json is detected in the working directory, the launcher will not overwrite it and start fresh unless new=true is specified. If this is the case and new=true is not specified then ideally Snakemake will continue where it left off int he workflow.

Other arguments (will be defaulted to value in defaults.json if not specified) :

```min_stack_proportion, classify_cutoff, dist_cutoff ```


454 specific arguments:

```pdiffs, bdiffs ``` (will be defaulted to what is in defaults.json if not specified)

miseq specific arguments:

*pipeline still in development*


==========
**Metadata file format:**

The metadata file is a tab delimited file that is provided by the user that describes each sample. In the future if the metadata file is not provided the pipeline will generate basic graphics using just sample name. For now, if you do not want to use a metadata file please provide a blank metadata file; the pipeline will break at the relevant graphing rules but all raw outputs will still have been generated. 
The metadata text file should be formatted like this:

| meta | meta | cat | cont | cat |
| :--- | :--- | :-- | :--- | :-- |
| **barcode** | **date** | **race** | **age** | **gender** |
| **ABC11** | **01212013** | **eu** | **54** | **m** |
|**ABC12** | **01212013** | **af** | **43** | **f** |

Please remember that this is all tab delimited. any column with the header as meta will be regarded as metadata and will not be used in graphing. Any column with the header as cat will be regarded as a categorical var for graphing, and cont will be regarded as a continuious var (example: age is continuous whereas gender is categorical). *The exact barcode for each sample must be provided under a meta column so that the pipeline will be able to properly match up output data to the appropriate row.


==========
**Notes:**

Feel free to open issues and submit possible fixes. For questions please contact Nikhil Gowda at nikhilbg96@gmail.com. Further documentation is still being provided.
