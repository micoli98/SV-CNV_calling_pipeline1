#!/usr/bin/env nextflow

/*
sCNA analysis of WGS data using Hartwig Medical Foundation tools (HMF)
Author: Giulia Micoli
Partially based on the public pipeline in github: ErasmusMC-Bioinformatics/NextFlow-VC-pipeline

Input (mandatory):
--sample_info: tsv file with tumor and normal labels and paths to bam files
--old_sample_info: previous run tsv file 
--pubDir: output directory
*/

projectDir = "/mnt/storageBig8/work/micoli/SCNA_Purple/src" /* Modify in favor of the location of the script*/

nextflow.enable.dsl=2
include { Prepare_input; Move; Gridss; Pon; Gripss; Cobalt; Amber; Purple; Results; Multiploidy } from '/mnt/storageBig8/work/micoli/SCNA_Purple/src/modules.nf' 

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
sample_info = Channel.value(file(params.sample_info))
old_sample_info = Channel.value(file(params.old_sample_info))
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
         sample_info        : ${params.sample_info}
         ols_sample_info    : ${params.old_sample_info}
         pubdir             : ${params.pubDir}
         """
         .stripIndent()

workflow {
    Prepare_input_results = Prepare_input(make_input, sample_info, old_sample_info)

    Move_results = Move(pubDir, Prepare_input_results.to_move)

    gridss_input = Prepare_input_results.gridss_ids
        .splitCsv(sep:'\t', header: true)
        .map{ row-> tuple(row.normalSample, row.patient, row.bams)}
    
    Gridss_results = Gridss(java, ref_genome_fa, bk_bed, pubDir, gridss_input)
    
    pon_input = Gridss_results.gridss_output.map { it.first() }
        .collect()
    
    Pon_results = Pon(pon_breakend, pon_breakpoint, pubDir, ref_genome_fa, pon_input)

    sample_info_split = Prepare_input_results.sample_ids
        .splitCsv(sep:'\t', header: ['sample', 'normalSample', 'bamFile', 'normalBamFile', 'patient'], skip: 1)
        .map{ row-> tuple(row.sample, row.normalSample, row.bamFile, row.normalBamFile, row.patient)}
    
    gripss_input = sample_info_split.combine(Gridss_results.gridss_output, by:1)

    Gripss_results = Gripss(java, pubDir, ref_genome_fa, pon_breakend, pon_breakpoint, gripss_input)

    Cobalt_results = Cobalt(java, ref_genome_fa, gc_profile, pubDir, sample_info_split)

    Amber_results = Amber(java, pubDir, loci_path, sample_info_split)

    purple_input = sample_info_split.join(Gripss_results.gripss_result_ch)
        .join(Gripss_results.gripss_filtered_ch)
        .join(Amber_results.amber_result_ch)
        .join(Cobalt_results.cobalt_output_ch)

    Purple_results = Purple(java, pubDir, ref_genome_fa, gc_profile, somatic_data, germline_data, 
        ensembl_dir, driver_catalog, germline_hs, somatic_hs, purple_input)

    purple_for_analysis = Purple_results.purple_results_ch
        .collect()
        .unique()
    
    Final_results = Results(pubDir, sample_info, segs_sunrises, purple_for_analysis, Move_results)

    Multiploidy_results = Multiploidy(sample_info, multiploidy, pubDir, purple_for_analysis, Move_results)
}

