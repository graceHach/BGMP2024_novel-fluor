//I think this just lets you ./ to run:
//#!/usr/bin/env nextflow

////////////////////////////// PARAMETERS //////////////////////////////


params.in_dir_blue = "/projects/bgmp/shared/groups/2024/novel-fluor/shared/upload/BGMP_2024/BLUE/NovaSeq_GC3F_7125/"
params.in_dir_red = "/projects/bgmp/shared/groups/2024/novel-fluor/shared/upload/BGMP_2024/RED/NovaSeq_GC3F_7124/"
params.out_dir_merge_blue = "illu-dat/NovaSeq_merged/blue/"
params.out_dir_merge_red = "illu-dat/NovaSeq_merged/red/"
params.out_dir_trim_blue = "illu-dat/primer_trimmed/blue/"
params.out_dir_trim_red = "illu-dat/primer_trimmed/red/"
params.hts_primers_report = "../reports/hts_primers_reports.txt"
params.forward_primers = "../primers/CVR205stub_FWD.fasta"
params.reverse_primers = "../primers/CVR205stub_REV.fasta"



process merge_reads_red {

    publishDir path: "${params.out_dir_merge_red}", mode: 'copy', overwrite: true
    cpus 8

    input:
    tuple val(key), path(read1_path), path(read2_path)

    output:
    path "${key}_MERGED.fastq"
    
    script:
    """
    conda activate bbmerge
    bbmerge.sh -in1=${read1_path} -in2=${read2_path} -out="${key}_MERGED.fastq" \
    -outu1="${key}_R1_REJCETED.fastq" -outu2="${key}_R2_REJCETED.fastq" 
    """
}

process merge_reads_blue {

    publishDir path: "${params.out_dir_merge_blue}", mode: 'copy', overwrite: true
    cpus 8

    input:
    tuple val(key), path(read1_path), path(read2_path)

    output:
    path "${key}_MERGED.fastq"
    
    script:
    """
    conda activate bbmerge
    bbmerge.sh -in1=${read1_path} -in2=${read2_path} -out="${key}_MERGED.fastq" \
    -outu1="${key}_R1_REJCETED.fastq" -outu2="${key}_R2_REJCETED.fastq" 
    """
}

process trim_reads_red {

    publishDir path: "${params.params.out_dir_trim_red}", mode: 'copy', overwrite: true
    //publishDir '.', saveAs: { it == "*.fastq.gz" ? "${params.params.out_dir_trim_red}${it}" : "${params.hts_primers_report}${it}" }

    input:
    path input_fasta_path
    val forward_primers
    val reverse_primers

    output:
    path "*_TRIMMED_SE.fastq.gz"
    //trimmed_filename=\$(basename "${input_fasta_path}" | sed 's/_MERGED.fastq//')
    
    script:
    """
    conda activate htstream
    hts_Primers -U "${input_fasta_path}" -f "{sed 's/_MERGED.fastq/_TRIMMED/' "${input_fasta_path}"}" \
    -P "${forward_primers}" -Q "${reverse_primers}" -l 5 -x -e 6 -d 6 -F
    """
}

process trim_reads_blue {

    publishDir path: "${params.params.out_dir_trim_blue}", mode: 'copy', overwrite: true
    //publishDir '.', saveAs: { it == "*.fastq.gz" ? "${params.params.out_dir_trim_red}${it}" : "${params.hts_primers_report}${it}" }

    input:
    path input_fasta_path
    val forward_primers
    val reverse_primers

    output:
    path "*_TRIMMED_SE.fastq.gz"
    //trimmed_filename=\$(basename "${input_fasta_path}" | sed 's/_MERGED.fastq//')
    
    script:
    """
    conda activate htstream
    hts_Primers -U "${input_fasta_path}" -f "{sed 's/_MERGED.fastq/_TRIMMED/' "${input_fasta_path}"}" \
    -P "${forward_primers}" -Q "${reverse_primers}" -l 5 -x -e 6 -d 6 -F
    """
}


////////////////////////////// WORKFLOW //////////////////////////////

workflow {

     // Create input channels for blue and red umerged datasets
    blue_files_ch = Channel.fromFilePairs("${params.in_dir_blue}*_{R1,R2}*.fastq.gz", flat: true)
    red_files_ch = Channel.fromFilePairs("${params.in_dir_red}*_{R1,R2}*.fastq.gz", flat: true)

    // Merge
    merge_reads_red(red_files_ch)
    merge_reads_blue(red_files_ch)

    // Create input channels to trim merged reads
    // blue_files_ch = Channel.fromPath("${params.out_dir_merge_blue}*_MERGED.fastq")
    // red_files_ch = Channel.fromPath("${params.out_dir_merge_red}*_MERGED.fastq")

    // VALUES, not paths to preserve directory structure
    forward_primers = Channel.of("${params.forward_primers}")
    reverse_primers = Channel.of("${params.reverse_primers}")

    trim_reads_red(merge_reads_red.out, forward_primers, reverse_primers)
    trim_reads_blue(merge_reads_blue.out, forward_primers, reverse_primers)




}
