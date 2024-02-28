// process CreateSampleFile {
//     publishDir "${params.outdir}/sortmerna", mode: 'copy'
//     input:
//     tuple val(meta), path(reads)
//      // tuple则表示channel传递了一个包含两个元素的元组。

//     output:
//     path "samples.txt"

//     script:
//     """
//     ch_filtered_reads
//     .map { meta, path ->
//         def read1
//         def read2 = ""

//         if (meta.single_end == true) {
//             read1 = path.toString()
//         } else {
//             read1 = path[0].toString()
//             read2 = path[1].toString()
//         }

//         return [meta.id, meta.strandedness, read1, read2].join('\t')
//     }
//     .collectFile(name: 'samples.txt', newLine: true, sort:true)
//     """
// }
