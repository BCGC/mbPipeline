rule graph:
    input:
        graphics_data = '{project}.mb_graphics_data.txt'.format(project=PROJECT),
        beta_data = '{project}.beta_data.out'.format(project=PROJECT),
        taxshared = '{project}.final.tax.shared'.format(project=PROJECT),
        taxconsensus = '{project}.final.tax.consensus'.format(project=PROJECT)
    output:
        'AlphaDiversity.pdf',
        'BetaDiversity.pdf',
        'NumSequences.pdf',
        'TaxonomicComposition.pdf'
    run:
        with open('run.json') as data_file:
            run = json.load(data_file)
        min_stack_proportion = run["setup"]["min_stack_proportion"]
        os.system("Rscript graphall.R "+input.taxconsensus+" "+input.taxshared+" "+min_stack_proportion+"")
