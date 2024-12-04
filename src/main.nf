

////////////////////////////// PARAMETERS //////////////////////////////


params.in_dir_blue = "../mini_illum/blue"
params.in_dir_red = "../mini_illum/red"

//params.out_dir_merge_blue = "/projects/bgmp/shared/groups/2024/novel-fluor/malm/illu-dat/NovaSeq_merged/blue"
//params.out_dir_merge_red = "/projects/bgmp/shared/groups/2024/novel-fluor/malm/illu-dat/NovaSeq_merged/red"

//params.out_dir_trim_blue = "/projects/bgmp/shared/groups/2024/novel-fluor/malm/illu-dat/trimmed_blue"
//params.out_dir_trim_red = "/projects/bgmp/shared/groups/2024/novel-fluor/malm/illu-dat/trimmed_red"

//params.out_dir_counts_blue = "/projects/bgmp/shared/groups/2024/novel-fluor/malm/illu-dat/counts_blue"
//params.out_dir_counts_red = "/projects/bgmp/shared/groups/2024/novel-fluor/malm/illu-dat/counts_red"

params.forward_primers = "../primers/CVR205stub_FWD.fasta"
params.reverse_primers = "../primers/CVR205stub_REV.fasta"

params.counts_script = "generate_counts_better.py"


////////////////////////////// WORKFLOW //////////////////////////////

workflow {

    // Input channels
    blue_files_ch = Channel.fromFilePairs("${params.in_dir_blue}/*{R1,R2}*.fastq.gz", flat: true)
    red_files_ch = Channel.fromFilePairs("${params.in_dir_red}/*{R1,R2}*.fastq.gz", flat: true)

    // Merge reads
    merged_blue_ch = merge_reads_blue(blue_files_ch)
    merged_red_ch = merge_reads_red(red_files_ch)

    // Trim reads
    trim_reads_blue(merged_blue_ch) | collect | generate_counts_blue
    trimmed_red_ch = trim_reads_red(merged_red_ch)

   
    // Generate counts
    // counts_blue_input_ch = trimmed_blue_ch.map { it -> "${params.out_dir_trim_blue}/${it}" }.view()
    // counts_red_input_ch = trimmed_red_ch.map { it -> "${params.out_dir_trim_red}/${it}" }.view()

    //counts_blue_ch = generate_counts_blue(trimmed_blue_ch)
    //counts_red_ch = generate_counts_red(trimmed_red_ch)

}

////////////////////////////// PROCESSES //////////////////////////////

// Merge Reads - Blue
process merge_reads_blue {

    //publishDir path: "${params.out_dir_merge_blue}", mode: 'copy', overwrite: true

    input:
    tuple val(key), path(read1_path), path(read2_path)

    output:
    tuple val(key), path("${key}_blue_MERGED.fastq"), emit: merged_blue

    script:
    """
    conda activate bbmerge
    bbmerge.sh -in1=${read1_path} -in2=${read2_path} \
        -out="${key}_blue_MERGED.fastq" \
        -outu1="${key}_blue_R1_REJECTED.fastq" \
        -outu2="${key}_blue_R2_REJECTED.fastq"
    """
}

// Merge Reads - Red
process merge_reads_red {

    //publishDir path: "${params.out_dir_merge_red}", mode: 'copy', overwrite: true

    input:
    tuple val(key), path(read1_path), path(read2_path)

    output:
    tuple val(key), path("${key}_red_MERGED.fastq"), emit: merged_red

    script:
    """
    conda activate bbmerge
    bbmerge.sh -in1=${read1_path} -in2=${read2_path} \
        -out="${key}_red_MERGED.fastq" \
        -outu1="${key}_red_R1_REJECTED.fastq" \
        -outu2="${key}_red_R2_REJECTED.fastq"
    """
}

// Trim Reads - Blue
process trim_reads_blue {

    //publishDir path: "${params.out_dir_trim_blue}", mode: 'copy', overwrite: true

    input:
    tuple val(key), path(merged_blue)

    output:
    // edited to only emit path, rather than tuple with key and path
    path ("*_TRIMMED*"), emit: trimmed_blue_path

    script:
    """
    conda activate htstream 
    out_file="\$(echo ${merged_blue} | sed s'/_MERGED.fastq/_TRIMMED/')"
    hts_Primers -U "${merged_blue}" -f "\${out_file}" \
    -P "$baseDir/${params.forward_primers}" -Q "$baseDir/${params.reverse_primers}" -l 5 -x -e 6 -d 6 -F
    """

}

process trim_reads_red {

    //publishDir path: "${params.out_dir_trim_red}", mode: 'copy', overwrite: true

    input:
    tuple val(key), path(merged_red)

    output:
    // This can't just be * or it will also output the stats.log from hts_Primers
    path ("*_TRIMMED*"), emit: trimmed_red_path

    script:
    """
    conda activate htstream 
    out_file="\$(echo ${merged_red} | sed s'/_MERGED.fastq/_TRIMMED/')"
    hts_Primers -U "${merged_red}" -f "\${out_file}" \
    -P "$baseDir/${params.forward_primers}" -Q "$baseDir/${params.reverse_primers}" -l 5 -x -e 6 -d 6 -F
    """
}

//Generate tsv count files 
process generate_counts_blue {

    //publishDir path: "${params.out_dir_counts_blue}", mode: 'copy', overwrite: true

    input:
    path trimmed_blue_path

    output:
    path "blue_counts.tsv", emit: blue_counts

    script:  
    """
    $baseDir/generate_counts_nf.py \
         --file_list "${trimmed_blue_path}" \
         --out_file "blue_counts.tsv"
    """
}

// process generate_counts_red {

//     publishDir path: "${params.out_dir_counts_red}", mode: 'copy', overwrite: true

//     input:
//     file("${params.out_dir_trim_red}/*.fastq.gz")

//     output:
//     path "red_counts_FINAL.tsv", emit: red_counts

//     script:
//      """
//     /projects/bgmp/shared/groups/2024/novel-fluor/malm/BGMP2024_novel-fluor/src/generate_counts_better.py \
//         --in_dir "${params.out_dir_trim_red}" \
//         --out_file "red_counts_FINAL.tsv"
//     """
// }