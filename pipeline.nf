#!/usr/bin/env nextflow

/*
sCNA analysis of WGS data using Hartwig Medical Foundation tools (HMF)
Author: Giulia Micoli
Partially based on the public pipeline in github: ErasmusMC-Bioinformatics/NextFlow-VC-pipeline

Input (mandatory):
--sample_names: tsv file with tumor and normal labels and paths to bam files
--pubDir: output directory
*/

projectDir = "/mnt/storageBig8/work/micoli/pipeline" /* Modify in favor of the location of the script*/

/* Generation of input channels */
bk_bed = Channel.value(file(params.bk_bed))
pon_breakend = Channel.value(file(params.pon_breakend))
pon_breakpoint = Channel.value(file(params.pon_breakpoint))
ensembl_dir = Channel.value(file(params.ensembl_dir))
somatic_data = Channel.value(file(params.somatic_data))
germline_data = Channel.value(file(params.germline_data))
ref_genome_fa = Channel.value(params.ref_genome_fa)
loci_path = Channel.value(file(params.loci_path))
gc_profile = Channel.value(file(params.gc_profile))
samples = Channel.value(file(params.sample_names))
java = Channel.value("/usr/lib/jvm/java-11-openjdk-amd64/bin/java")
pubDir = Channel.value(params.pubDir)
make_input = Channel.value(file(params.make_input))
multiploidy = Channel.value(file(params.multiploidy))
segs_sunrises = Channel.value(file(params.segs_sunrises))
driver_catalog = Channel.value(file(params.driver_gene_panel))
germline_hs = Channel.value(file(params.germline_hotspots))
somatic_hs = Channel.value(file(params.somatic_hotspots))

log.info """\

         HMF VARIANT CALLING AND CNV ANALYSIS PIPELINE     
         ===================================
         samples      : ${params.sample_names}
         pubdir       : ${params.pubDir}
         """
         .stripIndent()

/* Generation of the input channels for the tools from the sample tsv */
process Prepare_input {
    cache true
    input: 
    val make_input
    val samples
    val pubDir

    output:
    path "ids_table.tsv" into sample_ids 
    path "gridss_input.tsv" into gridss_ids

    script: 
    """
    Rscript $make_input $samples $pubDir 
    """
}

/* The resulting tables are split in rows */
gridss_ids
    .splitCsv(sep:'\t', header: true)
    .map{ row-> tuple(row.normalSample, row.patient, row.bams)}
    .set{ gridss_input }

/* Structural variant caller: GRIDSS*/
process Gridss {
    cpus 8
    memory '32 GB'
    cache true
    queueSize = 10
    publishDir "$pubDir/${patient}", mode: "copy"

    input:
    val java
    val ref_genome_fa
    file bk_bed
    val pubDir
    tuple normalSample, patient, bams from gridss_input

    output:
    tuple path("${normalSample}_calls.vcf"), val(normalSample) into gridss_output

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
    ln -s $pubDir/${patient}/${normalSample}_calls.vcf $pubDir/PoN/${normalSample}_calls.vcf

    #create compressed version and index file
    bcftools view $pubDir/${patient}/${normalSample}_calls.vcf -Oz -o $pubDir/${patient}/${normalSample}_calls.vcf.gz
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

gridss_output.into{
    gridss_pon;
    gridss_rest
}

gridss_pon.map { it.first() }
    .collect()
    .set{ pon_input }

process Pon {
    publishDir "$pubDir", mode: "copy"

    input:  
    val pon_breakend
    val pon_breakpoint
    val pubDir
    val ref_genome_fa
    path("*") from pon_input

    output:
    path("gridss_pon_breakpoint.bedpe") into pon_bedpe
    path("gridss_pon_single_breakend.bed") into pon_bed

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

sample_ids
    .splitCsv(sep:'\t', header: ['sample', 'normalSample', 'bamFile', 'normalBamFile', 'patient'], skip: 1)
    .map{ row-> tuple(row.sample, row.normalSample, row.bamFile, row.normalBamFile, row.patient)}    
    .into{
        sample_info_split_gripss;
        sample_info_split_Amber;
        sample_info_split_Cobalt;
        sample_info_split_Purple
    }

sample_info_split_gripss
    .combine(gridss_rest, by:1)
    .set { gripss_input }

/* Somatic filter: Gripss */
process Gripss {
    cpus 4
    memory '8 GB'
    cache true
    publishDir "$pubDir/${patient}", mode: "copy"

    input:
    val java
    val pubDir
    val ref_genome_fa
    path("gridss_pon_breakpoint.bedpe") from pon_bedpe
    path("gridss_pon_single_breakend.bed") from pon_bed 
    tuple normalSample, sample, bamFile, normalBamFile, patient, path("${normalSample}_calls.vcf") from gripss_input

    output:
    tuple val(sample), 
        path("${sample}.gripss.vcf.gz"), path("${sample}.gripss.vcf.gz.tbi") into gripss_result_ch
    tuple val(sample), 
        path("${sample}.gripss.filtered.vcf.gz"), path("${sample}.gripss.filtered.vcf.gz.tbi") into gripss_filtered_ch
    
    script:
    """
    $java -jar /mnt/storageBig8/work/micoli/pipeline/Scripts_Dependencies/gripss.jar \
        -sample  ${sample}\
        -reference ${normalSample} \
        -ref_genome $ref_genome_fa \
        -pon_sgl_file gridss_pon_single_breakend.bed \
        -pon_sv_file gridss_pon_breakpoint.bedpe \
        -vcf $pubDir/${patient}/${normalSample}_calls.vcf.gz \
        -output_dir .
    """
}

/* Extraction of the read depth and B allele frequency (BAF): Cobalt and Amber */
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
    tuple sample, normalSample, bamFile, normalBamFile, patient from sample_info_split_Cobalt

    output:
    tuple val(sample), path("*") into cobalt_output_ch

    script:
    """
    $java -Xmx8g -cp /mnt/storageBig8/work/micoli/pipeline/Scripts_Dependencies/cobalt.jar com.hartwig.hmftools.cobalt.CobaltApplication \
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
    val pubDir
    val java
    val loci_path
    tuple sample, normalSample, bamFile, normalBamFile, patient from sample_info_split_Amber

    output:
    tuple val(sample), path("*") into amber_result_ch

    script:
    """
    $java -Xmx32g -cp /mnt/storageBig8/work/micoli/pipeline/Scripts_Dependencies/amber.jar com.hartwig.hmftools.amber.AmberApplication \
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

/* Output of the previous processes are joined by the sample_id matching key */
sample_info_split_Purple.join(gripss_result_ch)
    .join(gripss_filtered_ch)
    .join(amber_result_ch)
    .join(cobalt_output_ch)
    .set{ purpleInput }

/* Ploidy and Purity estimation: Purple */
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
    tuple   sample, normalSample, bamFile, normalBamFile, patient, /* output from the sample information table */
            path("${sample}.gripss.vcf.gz"), path("${sample}.gripss.vcf.gz.tbi"),  /* soft-filtered SVs from gripss */
            path("${sample}.gripss.filtered.vcf.gz"), path("${sample}.gripss.filtered.vcf.gz.tbi"), /* hard-filtered SVs from gripss */
            path("${patient}/*"), /* output from amber */
            path("${patient}/*") from purpleInput /* output from cobalt */
    
    output:
    path "*.purple.*" into purple_results_ch

    script:
    """
    $java -jar /mnt/storageBig8/work/micoli/pipeline/Scripts_Dependencies/purple_v3.7.2.jar \
    -reference ${normalSample} \
    -tumor ${sample} \
    -amber ${patient} \
    -cobalt ${patient} \
    -gc_profile $gc_profile \
    -ref_genome $ref_genome_fa \
    -ensembl_data_dir $ensembl_dir \
    -structural_vcf ${sample}.gripss.filtered.vcf.gz \
    -sv_recovery_vcf ${sample}.gripss.vcf.gz \
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

/* Results from purple are collected in three channels so that the final modification of the results is started only after purple is completes */

purple_results_ch.collect()
    .unique()
    .into{
        purple_for_segments;
        purple_for_multiploidy
    }

/* Result process modifies purple results giving a table with all segmentation and sunrise plots */

process Results {
    publishDir "$pubDir", mode: "move"

    input: 
    val pubDir 
    val samples
    val segs_sunrises
    path ("*.purple.{purity|purity.range|somatic.hist|somatic.clonality|cnv.somatic}.tsv") from purple_for_segments 
    
    script: 
    """
    Rscript $segs_sunrises $samples $pubDir
    """
}

/* Multiploidy patients are retrieved and some plots are created for them */

process Multiploidy {
    publishDir "$pubDir/multiploidy", mode: "move"

    input:  
    val samples
    val multiploidy
    val pubDir
    path "*.purple.purity.range.tsv" from purple_for_multiploidy

    script: 
    """
    Rscript $multiploidy $samples $pubDir
    """
}

