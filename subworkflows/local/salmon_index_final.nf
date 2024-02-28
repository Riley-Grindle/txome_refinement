
//
// Load in modules required for workflow
//

include { GUNZIP as GUNZIP_FASTA            } from '../../modules/nf-core/gunzip/main'
include { GUNZIP as GUNZIP_GTF              } from '../../modules/nf-core/gunzip/main'
include { RSEM_PREPAREREFERENCE as MAKE_TRANSCRIPTS_FASTA       } from '../../modules/nf-core/rsem/preparereference/main'
include { SALMON_INDEX } from '../../modules/nf-core/salmon/index/main'

//

workflow SALMON_INDEX_FINAL {
    take:
    fasta                //      file: /path/to/genome.fasta
    gtf                  //      file: /path/to/genome.gtf

    main:

    ch_versions = Channel.empty()

    //
    // Uncompress genome fasta file if required
    //
    if (fasta.endsWith('.gz')) {
        ch_fasta    = GUNZIP_FASTA ( [ [:], fasta ] ).gunzip.map { it[1] }
        ch_versions = ch_versions.mix(GUNZIP_FASTA.out.versions)
    } else {
        ch_fasta = Channel.value(file(fasta))
    }     

    ch_transcript_fasta = MAKE_TRANSCRIPTS_FASTA( ch_fasta, gtf ).transcript_fasta
    ch_versions         = ch_versions.mix(MAKE_TRANSCRIPTS_FASTA.out.versions)

    ch_salmon           = SALMON_INDEX(ch_fasta, ch_transcript_fasta)
    ch_versions         = ch_versions.mix(SALMON_INDEX.out.versions)
    
    emit:
    index              = SALMON_INDEX.out.index
    transcript_fasta   = ch_transcript_fasta
    
    versions           = ch_versions
}

















