/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def valid_params = [
    aligners       : ['star_salmon', 'star_rsem', 'hisat2'],
    trimmers       : ['trimgalore', 'fastp'],
    pseudoaligners : ['salmon'],
    rseqc_modules  : ['bam_stat', 'inner_distance', 'infer_experiment', 'junction_annotation', 'junction_saturation', 'read_distribution', 'read_duplication', 'tin']
]

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowRnaseq.initialise(params, log, valid_params)

// Check input path parameters to see if they exist
checkPathParamList = [
    params.input, params.multiqc_config,
    params.fasta, params.transcript_fasta, params.additional_fasta,
    params.gtf, params.gff, params.gene_bed,
    params.ribo_database_manifest, params.splicesites,
    params.star_index, params.hisat2_index, params.rsem_index, params.salmon_index
]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }

// Check rRNA databases for sortmerna
if (params.remove_ribo_rna) {
    ch_ribo_db = file(params.ribo_database_manifest, checkIfExists: true)
    if (ch_ribo_db.isEmpty()) {exit 1, "File provided with --ribo_database_manifest is empty: ${ch_ribo_db.getName()}!"}
}

// Check if file with list of fastas is provided when running BBSplit
if (!params.skip_bbsplit && !params.bbsplit_index && params.bbsplit_fasta_list) {
    ch_bbsplit_fasta_list = file(params.bbsplit_fasta_list, checkIfExists: true)
    if (ch_bbsplit_fasta_list.isEmpty()) {exit 1, "File provided with --bbsplit_fasta_list is empty: ${ch_bbsplit_fasta_list.getName()}!"}
}

// Check alignment parameters
def prepareToolIndices  = []
if (!params.skip_bbsplit) { prepareToolIndices << 'bbsplit' }
if (!params.skip_alignment) { prepareToolIndices << params.aligner }
if (!params.skip_pseudo_alignment && params.pseudo_aligner) { prepareToolIndices << params.pseudo_aligner }

// Get RSeqC modules to run
def rseqc_modules = params.rseqc_modules ? params.rseqc_modules.split(',').collect{ it.trim().toLowerCase() } : []
if (params.bam_csi_index) {
    for (rseqc_module in ['read_distribution', 'inner_distance', 'tin']) {
        if (rseqc_modules.contains(rseqc_module)) {
            rseqc_modules.remove(rseqc_module)
        }
    }
}

// Stage dummy file to be used as an optional input where required
ch_dummy_file = file("$projectDir/assets/dummy_file.txt", checkIfExists: true)

// Check if an AWS iGenome has been provided to use the appropriate version of STAR
def is_aws_igenome = false
if (params.fasta && params.gtf) {
    if ((file(params.fasta).getName() - '.gz' == 'genome.fa') && (file(params.gtf).getName() - '.gz' == 'genes.gtf')) {
        is_aws_igenome = true
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

ch_multiqc_config          = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
ch_multiqc_custom_config   = params.multiqc_config ? Channel.fromPath( params.multiqc_config, checkIfExists: true ) : Channel.empty()
ch_multiqc_logo            = params.multiqc_logo   ? Channel.fromPath( params.multiqc_logo, checkIfExists: true ) : Channel.empty()
ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)

// Header files for MultiQC
ch_pca_header_multiqc        = file("$projectDir/assets/multiqc/deseq2_pca_header.txt", checkIfExists: true)
ch_clustering_header_multiqc = file("$projectDir/assets/multiqc/deseq2_clustering_header.txt", checkIfExists: true)
ch_biotypes_header_multiqc   = file("$projectDir/assets/multiqc/biotypes_header.txt", checkIfExists: true)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Loaded from modules/local/
//
include { BEDTOOLS_GENOMECOV                 } from '../modules/local/bedtools_genomecov'
include { DESEQ2_QC as DESEQ2_QC_STAR_SALMON } from '../modules/local/deseq2_qc'
include { DESEQ2_QC as DESEQ2_QC_RSEM        } from '../modules/local/deseq2_qc'
include { DESEQ2_QC as DESEQ2_QC_SALMON      } from '../modules/local/deseq2_qc'
include { DUPRADAR                           } from '../modules/local/dupradar'
include { MULTIQC                            } from '../modules/local/multiqc'
include { MULTIQC_CUSTOM_BIOTYPE             } from '../modules/local/multiqc_custom_biotype'
include { UMITOOLS_PREPAREFORRSEM as UMITOOLS_PREPAREFORSALMON } from '../modules/local/umitools_prepareforrsem.nf'
//customized modules
// include { TrinityNormalizeReads as TrinityNormalizeReads_SingleEnd } from '../modules/local/TrinityNormalization/trinity_normalization.nf'
// include { TrinityNormalizeReads as TrinityNormalizeReads_DoubleEnd } from '../modules/local/TrinityNormalization/trinity_normalization.nf'
//include { CreateSampleFile } from '../modules/local/Samples_file_for_trinity_normalization.nf'
// include { Staging as Staging_SingleEnd } from '../modules/local/create_samples_file_staging.nf'
// include { Staging as Staging_DoubleEnd } from '../modules/local/create_samples_file_staging.nf'

include { BAMSIFTER  } from '../modules/local/bamsifter.nf'
include { BAMSIFTER as BAMSIFTER_NORMALIZATION_MERGED_BAM} from '../modules/local/bamsifter.nf'
include { SAMTOOLS_MERGE } from '../modules/local/samtools_merge.nf'
include { ASSIGN_STRAND_AFTER_STRINGTIE  } from '../modules/local/assignstrandafterstringtie.nf'
include { AGAT_CONVERTSPGXF2GXF } from '../modules/local/convertspgxf2gxf/main'
include { AGAT_CONVERTSPGXF2GXF as GTF_FINAL_FORMATTING } from '../modules/local/convertspgxf2gxf/main'

// include { TRINITY_NORMALIZATION as TRINITY_NORMALIZATION_PARALLEL_DoubleEnd} from '../modules/local/trinity.nf'
// include { TRINITY_NORMALIZATION as TRINITY_NORMALIZATION_PARALLEL_SingleEnd} from '../modules/local/trinity.nf'

include { RSEM_PREPAREREFERENCE as MAKE_TRANSCRIPTS_FASTA_FROM_NEW_GTF      } from '../modules/nf-core/rsem/preparereference/main'
include { GTF_GENE_FILTER   as   GTF_GENE_FILTER_FROM_NEW_GTF                    } from '../modules/local/gtf_gene_filter.nf'
include { SALMON_INDEX   as   SALMON_INDEX_FROM_NEW_TRANSCRIPT_FASTA                 } from '../modules/nf-core/salmon/index/main'
include { GFFCOMPARE                   } from '../modules/local/gffcompare.nf'
include { GTF_INSERT } from '../modules/local/gtf_insert.nf'

//fastq after trinity normalization
//include { FASTQC as FASTQC_AFTER_TRINITY} from '../modules/nf-core/fastqc/main'

//StringTie merge modules
//include {STRINGTIE_MERGE} from '../modules/local/stringTie_merge/main.nf'


//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
//include { TRINITY_NORMALIZATION  } from '../subworkflows/local/trinity_normalization'
include { INPUT_CHECK    } from '../subworkflows/local/input_check'
include { PREPARE_GENOME } from '../subworkflows/local/prepare_genome'
include { ALIGN_STAR     } from '../subworkflows/local/align_star'
include { QUANTIFY_RSEM  } from '../subworkflows/local/quantify_rsem'
include { QUANTIFY_SALMON as QUANTIFY_STAR_SALMON } from '../subworkflows/local/quantify_salmon'
include { QUANTIFY_SALMON as QUANTIFY_SALMON      } from '../subworkflows/local/quantify_salmon'
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { CAT_FASTQ                   } from '../modules/nf-core/cat/fastq/main'
include { BBMAP_BBSPLIT               } from '../modules/nf-core/bbmap/bbsplit/main'

//!! we add this
include { PRESEQ_LCEXTRAP             } from '../modules/nf-core/preseq/lcextrap/main'
include { QUALIMAP_RNASEQ             } from '../modules/nf-core/qualimap/rnaseq/main'
include { SORTMERNA                   } from '../modules/nf-core/sortmerna/main'
include { STRINGTIE_STRINGTIE         } from '../modules/nf-core/stringtie/stringtie/main'
include { SUBREAD_FEATURECOUNTS       } from '../modules/nf-core/subread/featurecounts/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/custom/dumpsoftwareversions/main'
include { RSEM_PREPAREREFERENCE as MAKE_TRANSCRIPTS_FASTA	} from '../modules/local/preparegenome/main'
include { RSEM_PREPAREREFERENCE as MAKE_TRANSCRIPTS_FASTA_POST       } from '../modules/local/preparegenome/main'
include { SALMON_INDEX as SALMON_INDEX_FINAL	  } from '../modules/nf-core/salmon/index/main'
include { GUNZIP as GUNZIP_FASTA            } from '../modules/nf-core/gunzip/main'
include { SAMTOOLS_INDEX                    } from '../modules/local/samtools/index/main'
include { SAMTOOLS_SORT                     } from '../modules/local/samtools/sort/main'

//
// SUBWORKFLOW: Consisting entirely of nf-core/modules
//
include { FASTQ_SUBSAMPLE_FQ_SALMON        } from '../subworkflows/nf-core/fastq_subsample_fq_salmon/main'
include { FASTQ_FASTQC_UMITOOLS_TRIMGALORE } from '../subworkflows/nf-core/fastq_fastqc_umitools_trimgalore/main'
include { FASTQ_FASTQC_UMITOOLS_TRIMGALORE as FASTQ_FASTQC_UMITOOLS_TRIMGALORE_AFTER_TRINITY } from '../subworkflows/nf-core/fastq_fastqc_umitools_trimgalore/main'
include { FASTQ_FASTQC_UMITOOLS_FASTP      } from '../subworkflows/nf-core/fastq_fastqc_umitools_fastp/main'
include { FASTQ_ALIGN_HISAT2               } from '../subworkflows/nf-core/fastq_align_hisat2/main'
include { BAM_SORT_STATS_SAMTOOLS          } from '../subworkflows/nf-core/bam_sort_stats_samtools/main'
include { BAM_MARKDUPLICATES_PICARD        } from '../subworkflows/nf-core/bam_markduplicates_picard/main'
include { BAM_RSEQC                        } from '../subworkflows/nf-core/bam_rseqc/main'
include { BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS as BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS_GENOME        } from '../subworkflows/nf-core/bam_dedup_stats_samtools_umitools/main'
include { BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS as BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS_TRANSCRIPTOME } from '../subworkflows/nf-core/bam_dedup_stats_samtools_umitools/main'
include { BEDGRAPH_BEDCLIP_BEDGRAPHTOBIGWIG as BEDGRAPH_BEDCLIP_BEDGRAPHTOBIGWIG_FORWARD } from '../subworkflows/nf-core/bedgraph_bedclip_bedgraphtobigwig/main'
include { BEDGRAPH_BEDCLIP_BEDGRAPHTOBIGWIG as BEDGRAPH_BEDCLIP_BEDGRAPHTOBIGWIG_REVERSE } from '../subworkflows/nf-core/bedgraph_bedclip_bedgraphtobigwig/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary
def multiqc_report     = []
def pass_mapped_reads  = [:]
def pass_trimmed_reads = [:]
def pass_strand_check  = [:]

workflow RNASEQ_TRANSCRIPTOME_UPDATE {

    ch_versions = Channel.empty()
    
    //
    // MODULE: Standardize formatting of input gtf
    //
    if (!params.skip_agat) {

        AGAT_CONVERTSPGXF2GXF([[:], params.gtf])
        ch_formatted_gtf = AGAT_CONVERTSPGXF2GXF.out.output_gtf    
    
    } else {
        ch_formatted_gtf = Channel.of(params.gtf)
    }

    //
    // SUBWORKFLOW: Uncompress and prepare reference genome files
    //
    def biotype = params.gencode ? "gene_type" : params.featurecounts_group_type
    
    PREPARE_GENOME (
        params.fasta,
        ch_formatted_gtf.first(),
        params.gff,
        params.additional_fasta,
        params.transcript_fasta,
        params.gene_bed,
        params.splicesites,
        params.bbsplit_fasta_list,
        params.star_index,
        params.rsem_index,
        params.salmon_index,
        params.hisat2_index,
        params.bbsplit_index,
        params.gencode,
        is_aws_igenome,
        biotype,
        prepareToolIndices,
        params.genome_size
    )
    ch_versions = ch_versions.mix(PREPARE_GENOME.out.versions)
    
    // Check if contigs in genome fasta file > 512 Mbp
    if (!params.skip_alignment && !params.bam_csi_index) {
        PREPARE_GENOME
            .out
            .fai
            .map { WorkflowRnaseq.checkMaxContigSize(it, log) }
    }

    //
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //
    INPUT_CHECK (
        ch_input
    )
    .reads
    .map {
        meta, fastq ->
            new_id = meta.id - ~/_T\d+/
            [ meta + [id: new_id], fastq ]
    }
    .groupTuple()
    .branch {
        meta, fastq ->
            single  : fastq.size() == 1
                return [ meta, fastq.flatten() ]
            multiple: fastq.size() > 1
                return [ meta, fastq.flatten() ]
    }
    .set { ch_fastq }
    ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)

    //
    // MODULE: Concatenate FastQ files from same sample if required
    //
    CAT_FASTQ (
        ch_fastq.multiple
    )
    .reads
    .mix(ch_fastq.single)
    .set { ch_cat_fastq }
    ch_versions = ch_versions.mix(CAT_FASTQ.out.versions.first().ifEmpty(null))

    // Branch FastQ channels if 'auto' specified to infer strandedness
    ch_cat_fastq
        .branch {
            meta, fastq ->
                auto_strand : meta.strandedness == 'auto'
                    return [ meta, fastq ]
                known_strand: meta.strandedness != 'auto'
                    return [ meta, fastq ]
        }
        .set { ch_strand_fastq }

    //
    // SUBWORKFLOW: Sub-sample FastQ files and pseudo-align with Salmon to auto-infer strandedness
    //
    // Return empty channel if ch_strand_fastq.auto_strand is empty so salmon index isn't created
    PREPARE_GENOME.out.fasta
        .combine(ch_strand_fastq.auto_strand)
        .map { it.first() }
        .first()
        .set { ch_genome_fasta }

    FASTQ_SUBSAMPLE_FQ_SALMON (
        ch_strand_fastq.auto_strand,
        ch_genome_fasta,
        PREPARE_GENOME.out.transcript_fasta,
        PREPARE_GENOME.out.gtf,
        PREPARE_GENOME.out.salmon_index,
        !params.salmon_index && !('salmon' in prepareToolIndices)
    )
    ch_versions = ch_versions.mix(FASTQ_SUBSAMPLE_FQ_SALMON.out.versions)

    FASTQ_SUBSAMPLE_FQ_SALMON
        .out
        .json_info
        .join(ch_strand_fastq.auto_strand)
        .map { meta, json, reads ->
            return [ meta + [ strandedness: WorkflowRnaseq.getSalmonInferredStrandedness(json) ], reads ]
        }
        .mix(ch_strand_fastq.known_strand)
        .set { ch_strand_inferred_fastq }

    //
    // SUBWORKFLOW: Read QC, extract UMI and trim adapters with TrimGalore!
    // default setting trimmer is params.trimmer = 'trimgalore' in the nextflow.config
    ch_filtered_reads      = Channel.empty()
    ch_fastqc_raw_multiqc  = Channel.empty()
    ch_fastqc_trim_multiqc = Channel.empty()
    ch_trim_log_multiqc    = Channel.empty()
    ch_trim_read_count     = Channel.empty()
    if (params.trimmer == 'trimgalore') {
        FASTQ_FASTQC_UMITOOLS_TRIMGALORE (
            ch_strand_inferred_fastq,
            params.skip_fastqc || params.skip_qc,
            params.with_umi,
            params.skip_umi_extract,
            params.skip_trimming,
            params.umi_discard_read,
            params.min_trimmed_reads
        )
        ch_filtered_reads      = FASTQ_FASTQC_UMITOOLS_TRIMGALORE.out.reads
        ch_fastqc_raw_multiqc  = FASTQ_FASTQC_UMITOOLS_TRIMGALORE.out.fastqc_zip
        ch_fastqc_trim_multiqc = FASTQ_FASTQC_UMITOOLS_TRIMGALORE.out.trim_zip
        ch_trim_log_multiqc    = FASTQ_FASTQC_UMITOOLS_TRIMGALORE.out.trim_log
        ch_trim_read_count     = FASTQ_FASTQC_UMITOOLS_TRIMGALORE.out.trim_read_count
        ch_versions = ch_versions.mix(FASTQ_FASTQC_UMITOOLS_TRIMGALORE.out.versions)
    }

    //
    // SUBWORKFLOW: Read QC, extract UMI and trim adapters with fastp
    //
    if (params.trimmer == 'fastp') {
        FASTQ_FASTQC_UMITOOLS_FASTP (
            ch_strand_inferred_fastq,
            params.skip_fastqc || params.skip_qc,
            params.with_umi,
            params.skip_umi_extract,
            params.umi_discard_read,
            params.skip_trimming,
            [],
            params.save_trimmed,
            params.save_trimmed,
            params.min_trimmed_reads
        )
        ch_filtered_reads      = FASTQ_FASTQC_UMITOOLS_FASTP.out.reads
        ch_fastqc_raw_multiqc  = FASTQ_FASTQC_UMITOOLS_FASTP.out.fastqc_raw_zip
        ch_fastqc_trim_multiqc = FASTQ_FASTQC_UMITOOLS_FASTP.out.fastqc_trim_zip
        ch_trim_log_multiqc    = FASTQ_FASTQC_UMITOOLS_FASTP.out.trim_json
        ch_trim_read_count     = FASTQ_FASTQC_UMITOOLS_FASTP.out.trim_read_count
        ch_versions = ch_versions.mix(FASTQ_FASTQC_UMITOOLS_FASTP.out.versions)
    }

    //
    // Get list of samples that failed trimming threshold for MultiQC report
    //
    ch_trim_read_count
        .map {
            meta, num_reads ->
                pass_trimmed_reads[meta.id] = true
                if (num_reads <= params.min_trimmed_reads.toFloat()) {
                    pass_trimmed_reads[meta.id] = false
                    return [ "$meta.id\t$num_reads" ]
                }
        }
        .collect()
        .map {
            tsv_data ->
                def header = ["Sample", "Reads after trimming"]
                WorkflowRnaseq.multiqcTsvFromList(tsv_data, header)
        }
        .set { ch_fail_trimming_multiqc }

    //
    // MODULE: Remove genome contaminant reads
    //
    if (!params.skip_bbsplit) {
        BBMAP_BBSPLIT (
            ch_filtered_reads,
            PREPARE_GENOME.out.bbsplit_index,
            [],
            [ [], [] ],
            false
        )
        .primary_fastq
        .set { ch_filtered_reads }
        ch_versions = ch_versions.mix(BBMAP_BBSPLIT.out.versions.first())
    }

    //
    // MODULE: Remove ribosomal RNA reads
    //
    ch_sortmerna_multiqc = Channel.empty()
    if (params.remove_ribo_rna) {
        ch_sortmerna_fastas = Channel.from(ch_ribo_db.readLines()).map { row -> file(row, checkIfExists: true) }.collect()
//why we can use collect() to get ch_sortmerna_fastas, and pass ch_sortmerna_fastas into SORTMERNA?
// but do not use collect() for ch_filtered_reads because collect will broke the tuple stucture
//and collect() will not broke the structure for ch_sortmerna_fastas
        SORTMERNA (
            ch_filtered_reads,
            ch_sortmerna_fastas
        )
        .reads
        .view{ "Meta: ${it[0]}, Path: ${it[1]}" }
        .set { ch_filtered_reads }

        ch_sortmerna_multiqc = SORTMERNA.out.log
        ch_versions = ch_versions.mix(SORTMERNA.out.versions.first())
}


if (params.double_end_sample) {

ch_filtered_reads = ch_filtered_reads
    .filter { meta, path ->
        meta.single_end == false || meta.single_end == null  // null is for the case of undefined
    }
}


//In terms of Single end

if (params.single_end_sample) {

    ch_filtered_reads = ch_filtered_reads
        .filter { meta, path ->
            meta.single_end == true
        }

}

   // SUBWORKFLOW: Alignment with STAR and gene/transcript quantification with Salmon

    ch_genome_bam                 = Channel.empty()
    ch_genome_bam_index           = Channel.empty()
    ch_samtools_stats             = Channel.empty()
    ch_samtools_flagstat          = Channel.empty()
    ch_samtools_idxstats          = Channel.empty()
    ch_star_multiqc               = Channel.empty()
    ch_aligner_pca_multiqc        = Channel.empty()
    ch_aligner_clustering_multiqc = Channel.empty()
    if (!params.skip_alignment && params.aligner == 'star_salmon') {
        ALIGN_STAR (
            ch_filtered_reads,
            PREPARE_GENOME.out.star_index.map { [ [:], it ] },
            PREPARE_GENOME.out.gtf.map { [ [:], it ] },
            params.star_ignore_sjdbgtf,
            '',
            params.seq_center ?: '',
            is_aws_igenome,
            PREPARE_GENOME.out.fasta.map { [ [:], it ] }
        )
        //please check subworkflow align_stat, you will found emit： bam            = BAM_SORT_STATS_SAMTOOLS.out.bam
        // and module/nf-core star/align/main.nf

        //we use ch_genome_bam        = ALIGN_STAR.out.bam_sorted
        // to replace original bam            = ALIGN_STAR.out.bam


        ch_genome_bam        = ALIGN_STAR.out.bam
        ch_genome_bam_index  = ALIGN_STAR.out.bai
        ch_transcriptome_bam = ALIGN_STAR.out.bam_transcript
        ch_samtools_stats    = ALIGN_STAR.out.stats
        ch_samtools_flagstat = ALIGN_STAR.out.flagstat
        ch_samtools_idxstats = ALIGN_STAR.out.idxstats
        ch_star_multiqc      = ALIGN_STAR.out.log_final




        if (params.bam_csi_index) {
            ch_genome_bam_index = ALIGN_STAR.out.csi
        }
        ch_versions = ch_versions.mix(ALIGN_STAR.out.versions)



if(!params.skip_bamsifter) {
    ch_bamsifter = BAMSIFTER(ch_genome_bam).normalized_bam.collect()

    ch_bamsifter
        .flatMap()
        .map { item ->
            if(item instanceof Map) {
                item['id'] = params.gene_tx_prefix
                return item
            } else {
                return item
            }
        }
        .collate(2)
        .groupTuple()
        .map { meta, fastq ->
            [ meta, fastq.flatten() ]
        }
        .set{ ch_bamsifter_ready_samtools_merged }
        ch_bamsifter_ready_samtools_merged.view{ "Ready for Samtools_Merge  Meta: ${it[0]}, Path: ${it[1]}" }

    // use BAMSIFTER to normalize the bam files in parallel.
    SAMTOOLS_MERGE(ch_bamsifter_ready_samtools_merged,
        PREPARE_GENOME.out.fasta.map { [ [:], it ] },
        PREPARE_GENOME.out.fai.map { [ [:], it ] }
    )
} else {
    // IF skip BAMSIFTER，use CH_GENOME_BAM chanel
    SAMTOOLS_MERGE(ch_genome_bam,
        PREPARE_GENOME.out.fasta.map { [ [:], it ] },
        PREPARE_GENOME.out.fai.map { [ [:], it ] }
    )
}
ch_versions = ch_versions.mix(SAMTOOLS_MERGE.out.versions.first())



    SAMTOOLS_MERGE.out.bam.view()

    //BAMSIFTER_NORMALIZATION_MERGED_BAM agin after SAMTOOLS_MERGE

    if(!params.skip_bamsifter) {
    BAMSIFTER_NORMALIZATION_MERGED_BAM(SAMTOOLS_MERGE.out.bam)
    BAMSIFTER_NORMALIZATION_MERGED_BAM.out.normalized_bam.set{ch_genome_bam}
    ch_genome_bam.view{ "Ready for StringTie  Meta: ${it[0]}, Path: ${it[1]}" }}
    
    //
    // MODULE PAIR: SAMTOOLS_SORT - SAMTOOLS_INDEX
    //

    SAMTOOLS_SORT(ch_genome_bam)
    SAMTOOLS_INDEX(SAMTOOLS_SORT.out.bam)


    //
    // MODULE: STRINGTIE_STRINGTIE
    //
    if (!params.skip_alignment && !params.skip_stringtie) {
        STRINGTIE_STRINGTIE (
            SAMTOOLS_SORT.out.bam,
            PREPARE_GENOME.out.gtf
        )
        ch_versions = ch_versions.mix(STRINGTIE_STRINGTIE.out.versions.first())
    }
    //
    // Module: ASSIGN_STRAND_AFTER_STRINGTIE
    //
    if(params.stringtie_ignore_gtf){
    ASSIGN_STRAND_AFTER_STRINGTIE(STRINGTIE_STRINGTIE.out.transcript_gtf)}



    //
    // MODULE: GFFCOMPARE
    //

    // combine fasta and fai into one channel，make sure format is tuple val(meta2), path(fasta), path(fai)
ch_combine_fasta_fai = PREPARE_GENOME.out.fasta.combine(PREPARE_GENOME.out.fai)
                          .map { fasta, fai -> [ [:], fasta, fai ] }
    // convert reference_gtf int a format  tuple val(meta3), path(reference_gtf)
ch_reference_gtf = PREPARE_GENOME.out.gtf.map { [ [:], it ] }
// 21th Nov: STRINGTIE_STRINGTIE.out.transcript_gtf will be changed to  ASSIGN_STRAND_AFTER_STRINGTIE.out.processed_gtf (which is a gtf file of deleting unspeficified standness of the stringtie output gtf file without -e option )
    GFFCOMPARE(ASSIGN_STRAND_AFTER_STRINGTIE.out.processed_gtf,ch_combine_fasta_fai,ch_reference_gtf)

    ch_versions = ch_versions.mix(GFFCOMPARE.out.versions.first())


    //December 20th, 2024 we add this gtfinsert subworkflow
    // MODULE: GTFINSERT


    // Call the GTFINSERT subworkflow
    // if you only have one input sample, you can use ch_gtfinsert_input = GFFCOMPARE.out.annotated.gtf otherwise you need to use ch_gtfinsert_input = GFFCOMPARE.out.combined.gtf

    ch_gffcompare_gtf = params.onlyOneInputSample ? GFFCOMPARE.out.annotated_gtf : GFFCOMPARE.out.combined_gtf
    ch_new_gtf = GTF_INSERT(
        ch_gffcompare_gtf,
        GFFCOMPARE.out.tracking,
        ch_reference_gtf,
        GFFCOMPARE.out.loci,
        params.gene_tx_prefix

    )
    
    if (!params.skip_agat) {
        ch_final_gtf = GTF_FINAL_FORMATTING([[:], ch_new_gtf]).out.output_gtf
    } else {
        ch_new_gtf.set { ch_final_gtf }
    }

    if ((params.fasta).endsWith('.gz')) {
        ch_fasta_salmon    = GUNZIP_FASTA ( [ [:], params.fasta ] ).gunzip.map { it[1] }
        ch_versions = ch_versions.mix(GUNZIP_FASTA.out.versions)
    } else {
	ch_fasta_salmon = Channel.value(file(params.fasta))
    }


    MAKE_TRANSCRIPTS_FASTA_POST( ch_fasta_salmon, ch_final_gtf )
    ch_versions         = ch_versions.mix(MAKE_TRANSCRIPTS_FASTA_POST.out.versions)

    SALMON_INDEX_FINAL(ch_fasta_salmon, MAKE_TRANSCRIPTS_FASTA_POST.out.transcript_fasta)
    ch_versions         = ch_versions.mix(SALMON_INDEX_FINAL.out.versions)    

    ch_index  = SALMON_INDEX_FINAL.out.index
    ch_new_tx = MAKE_TRANSCRIPTS_FASTA_POST.out.transcript_fasta

    QUANTIFY_SALMON(
        ch_strand_fastq.auto_strand,
        ch_index.first(),
        ch_new_tx.first(),
        ch_final_gtf.first(),
        params.alignment_mode, 
        params.lib_type
    )    
}
}



/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report, pass_mapped_reads, pass_trimmed_reads, pass_strand_check)
    }

    if (params.hook_url) {
        NfcoreTemplate.IM_notification(workflow, params, summary_params, projectDir, log)
    }

    NfcoreTemplate.summary(workflow, params, log, pass_mapped_reads, pass_trimmed_reads, pass_strand_check)
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

