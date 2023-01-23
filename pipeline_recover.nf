k

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
make_input = Channel.value(file("/mnt/storageBig8/work/micoli/temporary_stuff/prepare_input.R"))
multiploidy = Channel.value(file(params.multiploidy))
segs_sunrises = Channel.value(file(params.segs_sunrises))
// cnGenes = Channel.value(file(params.cnGenes))
// genes_path = Channel.value(file(params.genes))
driver_catalog = Channel.value(file(params.driver_gene_panel))
germline_hs = Channel.value(file(params.germline_hotspots))
somatic_hs = Channel.value(file(params.somatic_hotspots))

pon_bedpe_gridss = Channel.value("/mnt/storageBig8/work/micoli/051122/PoN/gridss_pon_breakpoint.bedpe")
pon_bed_gridss = Channel.value("/mnt/storageBig8/work/micoli/051122/PoN/gridss_pon_single_breakend.bed")


process Prepare_input {
    cache true
    input: 
    val make_input
    val samples
    val pubDir

    output:
    path "ids.tsv" into ids 

    script: 
    """
    Rscript $make_input $samples $pubDir 
    """
} 

/*
samples
    .splitCsv(sep:'\t', header: true)
    .map{ row-> tuple(row.normalSample, row.patient, row.Bams)}
    .set { gridss_input }

/* Structural variant caller: GRIDSS
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
    tuple normalSample, patient, Bams from gridss_input

    output:
    tuple path("${normalSample}_calls.vcf"), val(normalSample) into gridss_output

    script:
    """
    export PATH=/usr/lib/jvm/java-11-openjdk-amd64/bin/:\$PATH

    /opt/share/gridss-2.13.2/gridss \
    -r $ref_genome_fa \
    -j  /opt/share/gridss-2.13.2/gridss-2.13.2-gridss-jar-with-dependencies.jar \
    -o ${normalSample}_calls.vcf \
    -b $bk_bed \
    ${Bams}
    """
} 

/*
process Collect_GRIDSS {
    input:
    tuple sample, patient, bamFile from gridss_input

    output:
    path ("*") into gridss_files 

    script:
    """
    ln -s /mnt/storageBig8/work/micoli/051122/${patient}/${sample}_calls.vcf /mnt/storageBig8/work/micoli/051122/PoN/${sample}_calls.vcf
    echo "${sample}" > ${sample}.txt
    """
}

gridss_files.collect()
    .set{ gridss_collected }

process Pon {
    publishDir "$pubDir/PoN", mode: "copy"

    input:  
    val pon_breakend
    val pon_breakpoint
    val pubDir
    val ref_genome_fa
    path("*") from gridss_collected

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
} */

ids
    .splitCsv(sep:'\t', header: ['sample', 'normalSample', 'patient', 'bamFile', 'normalBamFile'], skip: 1)
    .map{ row-> tuple(row.sample, row.normalSample, row.patient, row.bamFile, row.normalBamFile)}     
    .set{ sample_info_split_Purple}

/*
process Gripss {
    cpus 4
    memory '8 GB'
    cache true
    publishDir "$pubDir/${patient}", mode: "copy"

    input:
    val java
    val pubDir
    val ref_genome_fa
    val pon_bed_gridss
    val pon_bedpe_gridss
    tuple sample, normalSample, patient, bamFile, normalBamFile from sample_info_split_GRIPSS

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
        -pon_sgl_file $pon_bed_gridss \
        -pon_sv_file $pon_bedpe_gridss \
        -vcf /mnt/storageBig8/work/micoli/051122/${patient}/${normalSample}_calls.vcf \
        -output_dir .
    """
}

gripss_result_ch.collect()
    .set{ gripss_collected }


cobalt_ids
    .splitCsv(sep:'\t', header: ['sample', 'normalSample', 'bamFile', 'normalBamFile', 'patient'], skip: 1)
    .map{ row-> tuple(row.sample, row.normalSample, row.bamFile, row.normalBamFile, row.patient)}     
    .set{ sample_info_split_Cobalt} 

/* Extraction of the read depth and B allele frequency (BAF): Cobalt and Amber 
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
    tuple sample, true_normal, patient, bamFile, normalBamFile from cobalt_input

    output:
    tuple val(sample_id), path("*") into cobalt_result_ch

    script:
    """
    $java -Xmx8g -cp /mnt/storageBig8/work/micoli/pipeline/Scripts_Dependencies/cobalt.jar com.hartwig.hmftools.cobalt.CobaltApplication \
    -reference ${true_normal} \
    -reference_bam /mnt/storageBig8/resources/processed_data/HERCULES/WGSbams/${true_normal}.bam \
    -tumor ${sample} \
    -tumor_bam ${bamFile} \
    -output_dir $pubDir/${patient} \
    -threads 4 \
    -ref_genome $ref_genome_fa \
    -gc_profile $gc_profile
    """
} */

/*amber_ids
    .splitCsv(sep:'\t', header: ['sample', 'normalSample', 'bamFile', 'normalBamFile', 'patient'], skip: 1)
    .map{ row-> tuple(row.sample, row.normalSample, row.bamFile, row.normalBamFile, row.patient)}    
    .set{ sample_info_split_Amber}

process Amber {
    cpus 2
    memory '32 GB'
    cache true
    publishDir "$pubDir/${patient}", mode: "copy"

    input:
    val pubDir
    val java
    val loci_path
    tuple sample, true_normal, patient, bamFile, normalBamFile from amber_input

    output:
    tuple val(sample_id), path("*") into amber_result_ch

    script:
    """
    $java -Xmx32g -cp /mnt/storageBig8/work/micoli/pipeline/Scripts_Dependencies/amber.jar com.hartwig.hmftools.amber.AmberApplication \
    -reference ${true_normal} \
    -reference_bam /mnt/storageBig8/resources/processed_data/HERCULES/WGSbams/${true_normal}.bam \
    -tumor ${sample} \
    -tumor_bam ${bamFile} \
    -output_dir $pubDir/${patient} \
    -ref_genome_version 38 \
    -threads 4 \
    -loci $loci_path/${patient}.vcf.gz
    """
} */


process Purple {
    cpus 2
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
    tuple sample, normalSample, patient, bamFile, normalBamFile from sample_info_split_Purple
        
    output:
    path "*" into purple_results_ch

    script:
    """
    $java -jar /mnt/storageBig8/work/micoli/pipeline/Scripts_Dependencies/purple_v3.7.2.jar \
    -reference ${normalSample} \
    -tumor ${sample} \
    -amber $pubDir/${patient} \
    -cobalt $pubDir/${patient} \
    -gc_profile $gc_profile \
    -ref_genome $ref_genome_fa \
    -ensembl_data_dir $ensembl_dir \
    -structural_vcf $pubDir/${patient}/${sample}.gripss.filtered.vcf.gz \
    -sv_recovery_vcf $pubDir/${patient}/${sample}.gripss.vcf.gz \
    -germline_vcf $germline_data/${patient}.vcf.gz \
    -somatic_vcf $somatic_data/${patient}.vcf.gz \
    -run_drivers \
    -driver_gene_panel $driver_catalog \
    -somatic_hotspots $somatic_hs \
    -germline_hotspots $germline_hs \
    -ref_genome_version 38 \
    -output_dir $pubDir/${patient} \
    -no_charts
    """
} 

