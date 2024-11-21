//I think this just lets you ./ to run:
//#!/usr/bin/env nextflow

////////////////////////////// PARAMETERS //////////////////////////////


params.in_dir_blue = "/projects/bgmp/shared/groups/2024/novel-fluor/shared/upload/BGMP_2024/BLUE/NovaSeq_GC3F_7125/"
params.in_dir_red = "/projects/bgmp/shared/groups/2024/novel-fluor/shared/upload/BGMP_2024/RED/NovaSeq_GC3F_7124/"
params.out_dir_merge_blue = "illu-dat/NovaSeq_merged/blue/"
params.out_dir_merge_red = "illu-dat/NovaSeq_merged/red/"
params.out_dir_trim_blue = "/projects/bgmp/shared/groups/2024/novel-fluor/malm/illu-dat/blue_illum/"
params.out_dir_trim_red = "/projects/bgmp/shared/groups/2024/novel-fluor/malm/illu-dat/red_illum/"
params.forward_primers = "/projects/bgmp/shared/groups/2024/novel-fluor/malm/CVR205stub_FWD.fasta"
params.reverse_primers = "/projects/bgmp/shared/groups/2024/novel-fluor/malm/CVR205stub_REV.fasta"





process merge_reads {

    cpus 8

    input:
    tuple val(key), path(read1_path), path(read2_path)
    path path_out

    output:
    path "${path_out}${key}_MERGED.fastq"
    
    script:
    """
    mkdir -p ${path_out}
    conda activate bbmerge
    bbmerge.sh -in1=${read1_path} -in2=${read2_path} -out="${path_out}${key}_MERGED.fastq" \
    """
}


////////////////////////////// WORKFLOW //////////////////////////////

workflow {

    // Channel
    // .fromPath("${params.in_dir_red}*.fastq") 
    // .set { red_input_files }
    blue_files_ch = Channel.fromFilePairs("${params.in_dir_blue}*_{R1,R2}.fastq.gz", flat: true)
    red_merge_out = Channel.fromPath("${params.out_dir_merge_blue}")
        
        //.map { tuple(it.key, it.value[0], it.value[1]) } // (baseName, R1, R2)

    red_files_ch = Channel.fromFilePairs("${params.in_dir_red}*_{R1,R2}*.fastq.gz", flat: true)
    red_merge_out = Channel.fromPath("${params.out_dir_merge_red}")
    
    merge_reads(red_files_ch, red_merge_out)
    merge_reads(blue_files_ch, blue_merge_out)
    //merge_reads(blue_files_ch)
    //help_me()

}
