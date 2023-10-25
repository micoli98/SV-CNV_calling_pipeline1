#!/usr/bin/env nextflow

projectDir = "/mnt/storageBig8/work/micoli/SCNA_Purple/src/post_processing" /* Modify in favor of the location of the script*/

/* Generation of input channels */
ensembl_dir = Channel.value(file(params.ensembl_dir))
samples_ids = Channel.value(file(params.sample_ids))
java = Channel.value("/usr/lib/jvm/java-11-openjdk-amd64/bin/java")
pubDir = Channel.value(params.pubDir)
// input_script = Channel.value(file("/mnt/storageBig8/work/micoli/CNV_signatures/prepare_input.R"))

/*
process Input {
    input: 
    val input_script
    val samples

    output:
    path "ids.tsv" into linx_ids 

    script: 
    """
    Rscript $input_script $samples  
    """
} */

samples_ids
    .splitCsv(sep:'\t', header: ['sample', 'normalSample', 'bamFile', 'normalBamFile', 'patient'], skip: 1)
    .map{ row-> tuple(row.sample, row.normalSample, row.bamFile, row.normalBamFile, row.patient)}
    .set { linx_input }

process Linx {
    queueSize = 50
    publishDir "$pubDir/${patient}", mode: "copy"

    input:
    val java
    val pubDir
    val ensembl_dir
    tuple sample, normalSample, bamFile, normalBamFile, patient from linx_input

    output:
    tuple sample, patient into linx_out

    script:
    """
    $java -jar /mnt/storageBig8/work/micoli/SCNA_Purple/resources/dependencies/linx_v1.22.jar \
        -sample ${sample} \
        -ref_genome_version 38 \
        -sv_vcf $pubDir/${patient}/${sample}.purple.sv.vcf.gz \
        -purple_dir $pubDir/${patient} \
        -output_dir $pubDir/${patient} \
        -ensembl_data_dir $ensembl_dir \
        -check_fusions \
        -known_fusion_file /mnt/storageBig8/work/micoli/SCNA_Purple/resources/hmf_pipeline_resources.38_v5.31/sv/known_fusion_data.38.csv \
        -fragile_site_file /mnt/storageBig8/work/micoli/SCNA_Purple/resources/hmf_pipeline_resources.38_v5.31/sv/fragile_sites_hmf.38.csv \
        -line_element_file /mnt/storageBig8/work/micoli/SCNA_Purple/resources/hmf_pipeline_resources.38_v5.31/sv/line_elements.38.csv \
        -write_all
    """
} 

/*
process Visualization {
    cpus 8
    publishDir "$pubDir/${patient}", mode: "copy"

    input: 
    val java
    val pubDir
    val ensembl_dir
    tuple sample, patient from linx_out

    output: 
    path "*" into finish

    script:
    """
    mkdir $pubDir/$patient/circos_plots

    $java -cp /mnt/storageBig8/work/micoli/SCNA_Purple/resources/dependencies/linx_v1.22.jar com.hartwig.hmftools.linx.visualiser.SvVisualiser \
        -sample $sample \
        -ensembl_data_dir $ensembl_dir \
        -ref_genome_version 38 \
        -plot_out $pubDir/$patient/circos_plots \
        -data_out $pubDir/$patient/data \
        -vis_file_dir $pubDir/$patient/ \
        -circos /mnt/storageBig8/work/micoli/SCNA_Purple/resources/dependencies/circos-0.69-9/bin/circos \
        -threads 8
    """
} */
