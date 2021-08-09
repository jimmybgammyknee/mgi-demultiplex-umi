nextflow.enable.dsl=2

params.outdir = '/Users/jimmy.breen/Downloads/results'

// demultiplex mgi
params.demultiplex = '/homes/jimmy.breen/bioinfoCore/MGI/mgi/demultiplexing.py'

// custom refs
params.exons = '/apps/bioinfo/share/bcbio/genomes/Hsapiens/GRCh37/rnaseq/exons.hg19.bed'
params.bitRef = '/apps/bioinfo/share/bcbio/genomes/Hsapiens/hg38/ucsc/hg38.2bit'
params.bwa_idx = ''

include { markDupSambamba } from './modules/sambamba/0.7.1/RunSambambaMarkDup'
include { BwaMapSortedBam } from './modules/bwa/0.7.17/BwaMapSortedBam.nf'

// demultiplex with Ziad's script
// process demulti_umi {

//     tag { "Demulti - ${sample_id}" }
//     publishDir "${outdir}/demultiplex", mode: 'copy'
//     label 'process_low'

//     input:
//       tuple val(lane_id), file(reads)
//       path sample_index
//       val outdir
//     output:
//       path "*_R1.fastq.gz", emit: R1
//       path "*_R2.fastq.gz", emit: R2
//     script:
//       """
//       python ${params.demultiplex} \
//         -R1 reads[] -R2 []  \
//         -M mgi/example/lane3/sample_index -O lane3_split -RC -AM 1
//       """
// }

// some linux voodoo magic to split up the R2
process split_R2 {

    tag { "split_fastq - ${sample_id}" }
    publishDir "${outdir}/demultiplex", mode: 'copy'
    label 'process_low'

    input:
      tuple val(lane_id), file(reads)
    output:
      path "*_R1.fastq.gz", emit: R1
      path "*_R2.fastq.gz", emit: R2
      path "*_I1.fastq.gz", emit: I1
    script:
      """
      
      """
}

process umiadd {

    tag { "Fastp UMI-ADD - ${sample_id}" }
    publishDir "${outdir}/umi_add", mode: 'copy'
    label 'process_low'

    input:
    tuple val(sample_id), file(reads), file(I1)
    val outdir
    val opt_args

    output:
    tuple val(sample_id), file("${sample_id}_U{1,2}.fastq.gz"), emit: reads
    file("${sample_id}.{html,json}")

    script:
    def usr_args = opt_args ?: ''

    """
    fastp \
        ${usr_args} \
        -i ${reads[0]}  \
        -I ${I1}  \
        -o ${sample_id}_U1.fastq.gz \
        -O ${sample_id}_umi1.fastq.gz \
        --json ${sample_id}_U1.json \
        --html ${sample_id}_U1.html \
        --umi --umi_loc=read2 --umi_len=8 \
        -G -Q -A -L -w 1 -u 100 -n 8 -Y 100
    fastp \
        ${usr_args} \
        -i ${reads[1]}  \
        -I ${I1}  \
        -o ${sample_id}_U2.fastq.gz \
        -O ${sample_id}_umi2.fastq.gz \
        --json ${sample_id}_U2.json \
        --html ${sample_id}_U2.html \
        --umi --umi_loc=read2 --umi_len=8 \
        -G -Q -A -L -w 1 -u 100 -n 8 -Y 100
    """
}

workflow {

    ch_reads = Channel.fromFilePairs('data/*_R{1,2}.fastq.gz')

    // Genome files
	ch_exon = file(params.exons)
	ch_genome = file(params.reference)
    ch_bitRef = file(params.bitRef)

    // will needs some samtools flagstats in here

    BwaMapSortedBam(params.bwa_idx, ch_reads)
    markDupSambamba(BwaMapSortedBam.out)
    factera(params.factera, markDupSambamba.out, ch_exon, params.bitRef)

}