////////////////////////////// PARAMETERS //////////////////////////////

// default use the red dataset 
params.rawfastqs_R12 = '/projects/bgmp/shared/groups/2024/novel-fluor/shared/rawdata/RED/NovaSeq_GC3F_7124/*_R{1,2}_001.fastq.gz'
// point this to dir with fastq.gz data, both reads 1 and 2. Example: dir/*_R{1,2}.fastq.gz

params.out_dir = "$baseDir/outputs"
params.forward_primers = "$baseDir/primers/CVR205stub_FWD.fasta"
params.reverse_primers = "$baseDir/primers/CVR205stub_REV.fasta"

COUNTS_SCRIPT = "$baseDir/src/generate_counts_better.py"

////////////////////////////// WORKFLOW //////////////////////////////

workflow {

	binreads = Channel.fromFilePairs( params.rawfastqs_R12 ) 

	merge_reads( binreads ) 

	merge_reads.out.merged.view()

	trim_reads( merge_reads.out.merged ) | collect | generate_counts

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

	publishDir path: params.out_dir, mode: 'copy', overwrite: true

	input:
	path(mergedreads)

	output:
    	path ("trimmed*"), emit: trimmed

    	script:
    	"""
    	out_file="\$(echo ${mergedreads} | sed s'/merged/trimmed/')"
    	hts_Primers -U $mergedreads -f "\${out_file}" \
    	-P $params.forward_primers -Q $params.reverse_primers -l 5 -x -e 6 -d 6 -F
    	"""
}

//Generate tsv count files 
process generate_counts {

    publishDir path: "${params.out_dir}", mode: 'copy', overwrite: true

    input:
    path trimmed_path

    output:
    path "counts.tsv", emit: counts

    script:  
    """
    $COUNTS_SCRIPT \
         --file_list "${trimmed_path}" \
         --out_file "counts.tsv"
    """
}
