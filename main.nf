////////////////////////////// PARAMETERS //////////////////////////////

// default use the red dataset 
params.rawfastqs_R12 = '/projects/bgmp/shared/groups/2024/novel-fluor/shared/rawdata/RED/NovaSeq_GC3F_7124/*_R{1,2}_001.fastq.gz'

//params.in_dir_blue = "../mini_illum/blue"
//params.in_dir_red = "../mini_illum/red"

//params.out_dir_merge_blue = "/projects/bgmp/shared/groups/2024/novel-fluor/malm/illu-dat/NovaSeq_merged/blue"
//params.out_dir_merge_red = "/projects/bgmp/shared/groups/2024/novel-fluor/malm/illu-dat/NovaSeq_merged/red"

//params.out_dir_trim_blue = "/projects/bgmp/shared/groups/2024/novel-fluor/malm/illu-dat/trimmed_blue"
//params.out_dir_trim_red = "/projects/bgmp/shared/groups/2024/novel-fluor/malm/illu-dat/trimmed_red"

//params.out_dir_counts_blue = "/projects/bgmp/shared/groups/2024/novel-fluor/malm/illu-dat/counts_blue"
//params.out_dir_counts_red = "/projects/bgmp/shared/groups/2024/novel-fluor/malm/illu-dat/counts_red"

params.forward_primers = "./primers/CVR205stub_FWD.fasta"
params.reverse_primers = "./primers/CVR205stub_REV.fasta"

// params.counts_script = "src/generate_counts_better.py"

////////////////////////////// WORKFLOW //////////////////////////////

workflow {

	// binreads = Channel.fromFilePairs( params.rawfastqs_R12 ) 
	binreads = Channel.fromFilePairs( '/projects/bgmp/shared/groups/2024/novel-fluor/shared/rawdata/RED/NovaSeq_GC3F_7124/*_R{1,2}_001.fastq.gz' )

	merge_reads( binreads ) 

	merge_reads.out.merged.view()

	trim_reads( merge_reads.out.merged ) 

    // Input channels
    // blue_files_ch = Channel.fromFilePairs("${params.in_dir_blue}/*{R1,R2}*.fastq.gz", flat: true)
    // red_files_ch = Channel.fromFilePairs("${params.in_dir_red}/*{R1,R2}*.fastq.gz", flat: true)

    // Merge reads
    // merged_blue_ch = merge_reads_blue(blue_files_ch)
    // merged_red_ch = merge_reads_red(red_files_ch)

    // Trim reads
    // trim_reads_blue(merged_blue_ch) | collect | generate_counts_blue
    // trimmed_red_ch = trim_reads_red(merged_red_ch)

    // Generate counts
    // counts_blue_input_ch = trimmed_blue_ch.map { it -> "${params.out_dir_trim_blue}/${it}" }.view()
    // counts_red_input_ch = trimmed_red_ch.map { it -> "${params.out_dir_trim_red}/${it}" }.view()

    //counts_blue_ch = generate_counts_blue(trimmed_blue_ch)
    //counts_red_ch = generate_counts_red(trimmed_red_ch)

}

////////////////////////////// PROCESSES //////////////////////////////

process merge_reads {
	
	input: 
	tuple val(bin), path(reads) 

	output: 
	path("merged_${bin}.fastq"), emit: merged 

	script: 
	def (r1, r2) = reads 

	"""
    	bbmerge.sh \
		-in1="${r1}" \
		-in2="${r2}" \
        	-out="merged_${bin}.fastq" \
        	-outu1="rejected_${bin}_R1_001.fastq" \
        	-outu2="rejected_${bin}_R2_001.fastq"
	"""
}

process trim_reads {

    //publishDir path: "${params.out_dir_trim_blue}", mode: 'copy', overwrite: true

    input:
    //tuple val(key), path(mergedreads)
	path(mergedreads)

    output:
    // edited to only emit path, rather than tuple with key and path
    //path("trimmed_${mergedreads.baseName}.fastq"), emit: trimmed
    path ("trimmed*"), emit: trimmed

    script:
	def basefilename = mergedreads.baseName
    """
    out_file="\$(echo ${mergedreads} | sed s'/merged/trimmed/')"
    hts_Primers -U $mergedreads -f "\${out_file}" \
    -P $params.forward_primers -Q $params.reverse_primers -l 5 -x -e 6 -d 6 -F
    """
}

//process trim_reads_red {
//
    //publishDir path: "${params.out_dir_trim_red}", mode: 'copy', overwrite: true

    //input:
    //tuple val(key), path(merged_red)

    //output:
    //// This can't just be * or it will also output the stats.log from hts_Primers
    //path ("*_TRIMMED*"), emit: trimmed_red_path

    //script:
    //"""
    //conda activate htstream 
    //out_file="\$(echo ${merged_red} | sed s'/_MERGED.fastq/_TRIMMED/')"
    //hts_Primers -U "${merged_red}" -f "\${out_file}" \
    //-P "${params.forward_primers}" -Q "${params.reverse_primers}" -l 5 -x -e 6 -d 6 -F
    //"""
//}

//Generate tsv count files 
//process generate_counts_blue {

    //publishDir path: "${params.out_dir_counts_blue}", mode: 'copy', overwrite: true

    //input:
    //path trimmed_blue_path

    //output:
    //path "blue_counts.tsv", emit: blue_counts

    //script:  
    //"""
    //$baseDir/generate_counts_nf.py \
         //--file_list "${trimmed_blue_path}" \
         //--out_file "blue_counts.tsv"
    //"""
//}

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
