
nextflow.enable.dsl=2

params.outdir = '/Users/jimmy.breen/Downloads/results'

// demultiplex mgi
//params.demultiplex = '/homes/jimmy.breen/bioinfoCore/MGI/mgi/demultiplexing.py'

// custom refs

include { fastqc } from './modules/fastqc/0.11.9/fastqc.nf'
include { multiqc } from './modules/multiqc/1.9/multiqc.nf'
include { umiadd } from './modules/fastp/0.20.1/umiadd.nf'
include { starAlign } from './modules/STAR/2.7.3a/StarAlign'
include { UmiDedup } from './modules/umi_tools/1.0.0/UmiDedup.nf'
include { RunGatkBaseRecalibration } from './modules/fgatk/4.1.9/RunGatkBaseRecalibration.nf'
include { RunGatkHC } from './modules/gatk/4.1.9/RunGatkHC.nf'

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
      bash mgi_24bp_r2_split.sh ${reads}
      """
}

// Reordering BAM files
process runReorderSam {
    tag { filename + " - ReorderSam" }

    publishDir "${params.outDir}/${group}/${filename}/STAR/reorderSam", mode: 'copy'

    label 'gatkEnv'

    input:
    file ref from referenceMap.genomeFile
    file dict from referenceMap.genomeDict
    set group, 
        sample, 
        filename, 
        file(bam) from results_markDup

    output:
    set group, 
        sample, 
        filename, 
        file("${filename}.reorder.bam") into results_reorderBam

    script:
    """
    gatk ReorderSam \
    -I ${bam} \
    -O ${filename}.reorder.bam \
    -SD ${dict} \
    --TMP_DIR \${PWD}
    """
}

// SplitNCigarReads
process runSplitNCigarReads {

    tag { filename + " - SplitNCigarReads" }

    publishDir "${params.outDir}/STAR/splitNCigarReads", mode: 'copy'

    label 'gatkEnv'

    input:
    // val tmpDir from params.tmpDir
    file ref from referenceMap.genomeFile
    file idx from referenceMap.genomeIndex
    file dict from referenceMap.genomeDict
    set group, 
        sample, 
        filename, 
        file(bam) from results_reorderBam

    output:
    set group, 
        sample, 
        filename, 
        file("${filename}.split.bam"),
        file("${filename}.split.bai") into results_splitCigar

    script:
    """
    gatk SplitNCigarReads \
    -R ${ref} \
    -I ${bam} \
    -O ${filename}.split.bam \
    --tmp-dir \${PWD}
    """
}

workflow {

    ch_reads = Channel.fromFilePairs('data/*_R{1,2}.fastq.gz')

    // Genome files
	ch_genome_gtf = file(params.genome_gtf)
	ch_STAR_index = file(params.star_index)

    // QC
    fastqc(ch_reads,
           params.outdir,
           params.fastqc_optional_args,
		   params.sampleProject)
    
    fastqc.out.fastqc_output.map {id, f -> return(f)}.collect().set { ch_fastqc_files }

    multiqc(ch_fastqc_files,
            params.outdir,
			params.sampleProject)

    // The fun starts with splitting the weird 24bp R2
    split_R2(ch_reads)
    umiadd(split_R2.out, params.outdir)
    starAlign(ch_reads, ch_STAR_index, params.outdir)

    umiDedup(starAlign.out, params.outdir)

      // Collect files to be single input for feature counts
      umiDedup.out.bam.collect().set { ch_bams }
      umiDedup.out.bai.collect().set { ch_bai }
    
    runReorderSam
    runSplitNCigarReads
    RunGatkBaseRecalibration
    RunGatkHC
}