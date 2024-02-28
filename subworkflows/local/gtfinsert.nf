include { FIND_OVERLAP_TRANSCRIPT } from '../../modules/local/find_overlap_transcript'
include { PARSE_GTF_TO_DICT } from '../../modules/local/parse_gtf_to_dict_by_cmp_ref'
include { EXTRACT_KEY_VALUE } from '../../modules/local/extract_key_value_pairs'
include { CHECK_UNIQUE_KEYS }  from '../../modules/local/check_key_unique'
include { UPDATE_KEYS } from '../../modules/local/update_keys_with_tracking_key_value_pair'
include { PARSE_GTF_TO_DICT_BY_GENEID } from '../../modules/local/parse_gtf_to_dict_by_geneID'
include { INSERT_BY_START_POSITION } from '../../modules/local/insert_by_start_position'
include { JSON_TO_GTF} from '../../modules/local/parse_json_gtf'
include { FIND_NOVEL_TRANSCRIPTS } from '../../modules/local/find_novel_transcripts'


workflow GTFINSERT {
    take:
    combined_gtf   //channel: [ val(meta), [ path ] ]
    tracking_file  //channel: [ val(meta), [ path ] ]
    reference_gtf //channel: [ val(meta), [ path ] ]

    main:
    // Step 1: Find overlapping transcripts
    FIND_OVERLAP_TRANSCRIPT(combined_gtf)

    // Step 2: Parse GTF to dict by cmp_ref
    PARSE_GTF_TO_DICT (FIND_OVERLAP_TRANSCRIPT.out.select_transcript_insert_based_on_gffcompare_class_code)

    // Step 3: Extract key-value pairs
    EXTRACT_KEY_VALUE(tracking_file)

    // Step 4: (Optional) Check if keys are unique
    CHECK_UNIQUE_KEYS(EXTRACT_KEY_VALUE.out.key_value_json)

    // Step 5: Update keys with tracking key-value pair
    UPDATE_KEYS (EXTRACT_KEY_VALUE.out.key_value_json, PARSE_GTF_TO_DICT.out.parsed_gtf_json_by_cmp_ref)

    // Process Reference GTF
    // Step 6: Parse GTF to dict by gene ID
    PARSE_GTF_TO_DICT_BY_GENEID(reference_gtf)

    // Insertion and Sorting
    // Step 7: Insert by start position
    INSERT_BY_START_POSITION(UPDATE_KEYS.out.update_keys_with_tracking_key_value_pair, PARSE_GTF_TO_DICT_BY_GENEID.out.parsed_gtf_json_by_geneID)

    // Finalization
    // Step 8: Parse JSON GTF
    JSON_TO_GTF(INSERT_BY_START_POSITION.out.finished_insert_overlap_transcript_json)

    // Handling Novel Transcripts
    // Step 9: Find novel transcripts
    FIND_NOVEL_TRANSCRIPTS(combined_gtf,JSON_TO_GTF.out.finished_insert_overlap_transcript_gtf)

    emit:
    final_gtf = FIND_NOVEL_TRANSCRIPTS.out.final_annotation_gtf  //select_transcript_append_based_on_gffcompare_class_code_to_generate_final_gtf




}
