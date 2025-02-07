#!/usr/bin/env nextflow
params.input = "$baseDir/in"
params.output = "$baseDir/out"
params.gene_result_column = 0
params.gzip = false
params.species = null
params.coverage = null
params.threshold = null
params.acquired = true
params.point = false
params.unknown_mut = false

args = []
if (params.acquired) args.push("-acq")
if (params.point) args.push("-c")
if (params.unknown_mut) args.push("-u")
if (params.threshold) args.push("-t $params.threshold")
if (params.species) args.push("-s $params.species")
if (params.coverage) args.push("-l $params.coverage")

process RESFINDER{
    publishDir params.output, mode: 'copy'

    input:
    path fasta

    output:
    path "*.tsv"

    script:
    def prefix = fasta.getSimpleName()
    def is_compressed = fasta.getName().endsWith(".gz") ? true : false
    def fasta_name = fasta.getName().replace(".gz", "")

    """
    if [ "$is_compressed" == "true" ]; then
        gzip -c -d $fasta > $fasta_name
    fi
    python -m resfinder ${args.join(' ')} -ifa $fasta_name -o '$fasta_name'.result
    mv '$fasta_name'.result/ResFinder_results_tab.txt '$fasta_name'.tsv
    """
}

process CSV{
    publishDir params.output, mode: 'copy'

    input:
    val tables

    output:
    path 'resfinder_results.csv'

    exec:
    gene_list = []
    results = [:]
    tables.each { table ->
        sample_genes = []

        table
            .splitCsv(header:true,sep:"\t")
            .each {row -> sample_genes.push(row["Resistance gene"])}

        sample_genes.unique()
        gene_list += sample_genes
        sample_name = table.getSimpleName()
        results[sample_name] = sample_genes
    }
    result_table = ""
    gene_list.unique().sort()

    results = results.sort()
    results.each{ sample_name, genes ->
        result_row = []
        gene_list.each { gene ->
            if (genes.contains(gene)){
                result_row += 1
            } else{
                result_row += 0
            }
        }
        result_row.push(sample_name)
        result_table += result_row.join(',') + "\n"
    }
    gene_list.sort()
    gene_list.push('Isolate')
    headers = gene_list.join(',') + "\n"
    result_table = headers + result_table

    csv_file = task.workDir.resolve('resfinder_results.csv')
    csv_file.text = result_table
}

process ZIP{
    publishDir params.output, mode: 'copy'

    input:
    path files
    path csv

    output:
    path '*.tar.gz'

    """
    current_date=\$(date +"%Y-%m-%d")
    outfile="resfinder_\${current_date}.tar.gz"
    tar -chzf \${outfile} ${files.join(' ')} $csv
    """
}


workflow {
    input_seqs = Channel
        .fromPath("$params.input/*{fas,gz,fasta,fsa,fsa.gz,fas.gz}")

    RESFINDER(input_seqs)
    CSV(RESFINDER.out.collect())
    if (params.gzip){
        ZIP(RESFINDER.out.collect(),CSV.out)
    }
}
