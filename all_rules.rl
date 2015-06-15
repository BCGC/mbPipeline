#snakemake rules for microbiome pipeline
#these rules will be seperated into individual files upon completion

rule all:
    input:
        'AlphaDiversity.pdf',
        'BetaDiversity.pdf',
        'NumSequences.pdf',
        'TaxonomicComposition.pdf'

rule graph:
    input:
        graphics_data = 'mb_graphics_data.txt',
        beta_data = 'beta_data.out',
        taxonomy = '{project}.final.taxdata'
    output:
        'AlphaDiversity.pdf',
        'BetaDiversity.pdf',
        'NumSequences.pdf',
        'TaxonomicComposition.pdf'
    run:
        #FILLER

rule data_setup:
    input: '.temp.adiv'
    output: 'mb_graphics_data.txt',
    run:
        #FILLER

rule beta_data:
    input: '{project}.final.beta.shared'
    output: 'beta_data.out'
    run:
        #FILLER

rule process_otu:
    input:
        fasta='{project}.final.fasta',
        names='{project}.final.names',
        taxonomy='{project}.final.taxonomy',
        groups='{project}.final.groups'
    output:
        '.temp.adiv', '{project}.final.taxdata', '{project}.final.beta.shared'
    run:
        #FILLER
        
rule finalize_sequences:
    input:
        fasta='{project}.remove.fasta',
        names='{project}.remove.names',
        groups='{project}.remove.groups'
    output:
        '{project}.final.fasta',
        '{project}.final.names',
        '{project}.final.taxonomy',
        '{project}.final.groups'
    run:
        #FILLER

rule remove_454:
    input:
        fasta='{project}.processed.fasta',
        names='{project}.processed.names',
        groups='{project}.processed.groups'
    output:
        '{project}.remove.fasta',
        '{project}.remove.names',
        '{project}.remove.groups'
    run:
        #FILLER

rule remove_miseq:
    input:
        fasta='{project}.processed.fasta',
        names='{project}.processed.names',
        groups='{project}.processed.groups'
    output:
        '{project}.remove.fasta',
        '{project}.remove.names',
        '{project}.remove.groups'
    run:
        #FILLER

rule process_sequences:
    input:
        fasta='{project}.unique.fasta',
        names='{project}.unique.names'
        groups='{project}.groups'
    output:
        '{project}.processed.fasta',
        '{project}.processed.names',
        '{project}.processed.groups'
    run:
        #FILLER

rule unique_sequences:
    input:
        fasta='{project}.trim.fasta',
        names='{project}.trim.names'
    output:
        '{project}.unique.fasta',
        '{project}.unique.names'
    run:
        #FILLER

rule trim_sequences:
    input:
        fasta='{project}.fasta',
        names='{project}.names'
    output:
        '{project}.trim.fasta',
        '{project}.trim.names'
    run:
        #FILLER

rule load_454:
    output: '{project}.fasta', '{project}.names', '{project}.groups'
    run:
        #FILLER

rule load_miseq:
    #FILLER


        
