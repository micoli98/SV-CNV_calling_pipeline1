process Prepare_input {
    cache = true

    input: 
    val make_input
    val sample_info
    val old_sample_info
    val custom_patients

    output:
    path 'ids_table.tsv', emit: sample_ids 
    path 'gridss_input.tsv', emit: gridss_ids
    path 'patients_to_move.tsv', emit: to_move

    script: 
    """
    Rscript $make_input $sample_info $old_sample_info $custom_patients
    """
}

process Move {
    cache = true

    input:
    val pubDir
    path patients_to_move

    output:
    path 'patients_moved.txt', emit: move_confirm

    script:
    """
    # Create or empty the output file
    > patients_moved.txt

    # Read the tsv file and iterate over each patient
    while IFS= read -r patient; do
        # Check if the patient's directory exists in old_path
        if [ -d "/mnt/storageBig8/work/micoli/SCNA_Purple/results/231006/\${patient}" ]; then
            # If it exists, copy the directory to the output_path
            cp -r "/mnt/storageBig8/work/micoli/SCNA_Purple/results/231006/\${patient}" "${pubDir}/\${patient}"
            # Append the patient's name to patients_moved.txt
            echo \$patient >> patients_moved.txt
        else
            # Print a warning if the directory does not exist
            echo "Warning: Directory for patient \${patient} not found."
        fi
    done < <(tail -n +2 ${patients_to_move}) # This will skip the header of the tsv file

    """ 
}

process Gridss {
    cpus 8
    memory '32 GB'
    cache true
    queueSize = 10
    publishDir "$pubDir/${patient}", mode: "copy"

    input:
    val java
    val ref_genome_fa
    path bk_bed
    val pubDir
    tuple val(normalSample), val(patient), val(bams)

    output:
    tuple path("${normalSample}_calls.vcf"), val(normalSample), emit: gridss_output

    script:
    """
    export PATH=/usr/lib/jvm/java-11-openjdk-amd64/bin/:\$PATH

    #GRIDSS process
    /opt/share/gridss-2.13.2/gridss \
    -r $ref_genome_fa \
    -j  /opt/share/gridss-2.13.2/gridss-2.13.2-gridss-jar-with-dependencies.jar \
    -o ${normalSample}_calls.vcf \
    -b $bk_bed \
    ${bams}

    #copy in PoN directory
    ln -sf ${normalSample}_calls.vcf $pubDir/PoN/${normalSample}_calls.vcf

    #create compressed version and index file
    bcftools view ${normalSample}_calls.vcf -Oz -o $pubDir/${patient}/${normalSample}_calls.vcf.gz
    bcftools index $pubDir/${patient}/${normalSample}_calls.vcf.gz -t -o $pubDir/${patient}/${normalSample}_calls.vcf.gz.tbi

    #modify header and create the new version of vcf.gz
    bcftools view -h $pubDir/${patient}/${normalSample}_calls.vcf.gz | sed 's/.bam//g' > header.vcf
    bcftools view -H $pubDir/${patient}/${normalSample}_calls.vcf.gz -o body.vcf
    cat header.vcf body.vcf | bcftools view -Oz -o $pubDir/${patient}/${normalSample}_calls_modified.vcf.gz

    #remove header and body
    rm header.vcf
    rm body.vcf
    """
}

process Pon {
    publishDir "$pubDir", mode: "copy"

    input:  
    val pon_breakend
    val pon_breakpoint
    val pubDir
    val ref_genome_fa
    path pon_input

    output:
    path("gridss_pon_breakpoint.bedpe"), emit: pon_bedpe
    path("gridss_pon_single_breakend.bed"), emit: pon_bed

    script: 
    """
    java -Xmx8g \
	-cp /opt/share/gridss-2.13.2/gridss-2.13.2-gridss-jar-with-dependencies.jar \
	gridss.GeneratePonBedpe \
    \$(ls -1 $pubDir/PoN/*.vcf | awk ' { print "INPUT=" \$0 }' | head -\$n) \
	INPUT_BEDPE= $pon_breakpoint \
	INPUT_BED= $pon_breakend \
	O=gridss_pon_breakpoint.bedpe \
	NORMAL_ORDINAL=0\
	SBO=gridss_pon_single_breakend.bed \
	REFERENCE_SEQUENCE= $ref_genome_fa
    """
}

process Gripss {
    cpus 4
    memory '8 GB'
    cache true
    publishDir "$pubDir/${patient}", mode: "copy"

    input:
    val java
    val pubDir
    val ref_genome_fa
    path pon_breakend
    path pon_breakpoint 
    tuple val(normalSample), val(sample), val(bamFile), val(normalBamFile), val(patient), path("${normalSample}_calls.vcf")

    output:
    tuple val(sample), 
        path("${sample}.gripss.vcf.gz"), 
        path("${sample}.gripss.vcf.gz.tbi"), emit: gripss_result_ch
    tuple val(sample), 
        path("${sample}.gripss.filtered.vcf.gz"), 
        path("${sample}.gripss.filtered.vcf.gz.tbi"), emit: gripss_filtered_ch
    
    script:
    """
    $java -jar /mnt/storageBig8/work/micoli/SCNA_Purple/resources/dependencies/gripss.jar \
        -sample  ${sample}\
        -reference ${normalSample} \
        -ref_genome $ref_genome_fa \
        -pon_sgl_file gridss_pon_single_breakend.bed \
        -pon_sv_file gridss_pon_breakpoint.bedpe \
        -vcf $pubDir/${patient}/${normalSample}_calls.vcf.gz \
        -output_dir .
    """
}

process Cobalt {
    cpus 2
    memory '8 GB'
    cache true
    publishDir "$pubDir/${patient}", mode: "copy"

    input:
    val java
    val ref_genome_fa
    val gc_profile
    val pubDir
    tuple val(sample), val(normalSample), val(bamFile), val(normalBamFile), val(patient)

    output:
    tuple val(sample), path("*"), emit: cobalt_output_ch

    script:
    """
    $java -Xmx8g -cp /mnt/storageBig8/work/micoli/SCNA_Purple/resources/dependencies/cobalt.jar com.hartwig.hmftools.cobalt.CobaltApplication \
    -reference ${normalSample} \
    -reference_bam ${normalBamFile} \
    -tumor ${sample} \
    -tumor_bam ${bamFile} \
    -output_dir . \
    -threads 4 \
    -ref_genome $ref_genome_fa \
    -gc_profile $gc_profile
    """
}

process Amber {
    cpus 2
    memory '32 GB'
    cache true
    publishDir "$pubDir/${patient}", mode: "copy"

    input:
    val java
    val pubDir
    val loci_path
    tuple val(sample), val(normalSample), val(bamFile), val(normalBamFile), val(patient) 

    output:
    tuple val(sample), path("*"), emit: amber_result_ch

    script:
    """
    $java -Xmx32g -cp /mnt/storageBig8/work/micoli/SCNA_Purple/resources/dependencies/amber.jar com.hartwig.hmftools.amber.AmberApplication \
    -reference ${normalSample} \
    -reference_bam ${normalBamFile} \
    -tumor ${sample} \
    -tumor_bam ${bamFile} \
    -output_dir . \
    -ref_genome_version 38 \
    -threads 4 \
    -loci $loci_path/${patient}.vcf.gz
    """
}

process Purple {
    cpus 2
    memory '32 GB'
    publishDir "$pubDir/${patient}", mode: "copy"

    input:
    val java
    val pubDir
    val ref_genome_fa
    val gc_profile
    val somatic_data
    val germline_data
    val ensembl_dir
    val driver_catalog
    val germline_hs
    val somatic_hs
    tuple   val(sample), val(normalSample), val(bamFile), val(normalBamFile), val(patient), /* output from the sample information table */
            path("${sample}.gripss.vcf.gz"), path("${sample}.gripss.vcf.gz.tbi"),  /* soft-filtered SVs from gripss */
            path("${sample}.gripss.filtered.vcf.gz"), path("${sample}.gripss.filtered.vcf.gz.tbi"), /* hard-filtered SVs from gripss */
            path("${patient}/*"), /* output from amber */
            path("${patient}/*") /* output from cobalt */
    
    output:
    path "*" , emit: purple_results_ch

    script:
    """
    $java -Xmx32g -jar /mnt/storageBig8/work/micoli/SCNA_Purple/resources/dependencies/purple_v3.7.2.jar \
    -reference ${normalSample} \
    -tumor ${sample} \
    -amber ${patient} \
    -cobalt ${patient} \
    -gc_profile $gc_profile \
    -ref_genome $ref_genome_fa \
    -ensembl_data_dir $ensembl_dir \
    -structural_vcf ${sample}.gripss.filtered.vcf.gz \
    -sv_recovery_vcf ${sample}.gripss.vcf.gz \
    -germline_vcf $germline_data/${patient}.vcf.gz \
    -somatic_vcf $somatic_data/${patient}.vcf.gz \
    -run_drivers \
    -driver_gene_panel $driver_catalog \
    -somatic_hotspots $somatic_hs \
    -germline_hotspots $germline_hs \
    -ref_genome_version 38 \
    -output_dir . \
    -no_charts
    """
}

process Results {
    publishDir "$pubDir", mode: "move"

    input: 
    val pubDir 
    val sample_info
    val segs_sunrises
    path ("*.purple.{purity|purity.range|somatic.hist|somatic.clonality|cnv.somatic}.tsv") 
    path ("patients_moved.txt")
    
    script: 
    """
    Rscript $segs_sunrises $sample_info $pubDir
    """
}

/* Multiploidy patients are retrieved and some plots are created for them */

process Multiploidy {
    publishDir "$pubDir/multiploidy", mode: "move"

    input:  
    val sample_info
    val multiploidy
    val pubDir
    path ("*.purple.purity.range.tsv")
    path ("patients_moved.txt")

    script: 
    """
    Rscript $multiploidy $sample_info $pubDir
    """
}