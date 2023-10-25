#!/bin/bash

while read patient normalSample samples; do        
        bcftools view -s $samples /mnt/storageBig8/work/micoli/051122/$patient/${normalSample}_calls.vcf > /mnt/storageBig8/work/micoli/051122/$patient/${normalSample}_calls_filtered.vcf
done < /mnt/storageBig8/work/micoli/temporary_stuff/gripss_to_correct.tsv