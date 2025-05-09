# umiQC

QC workflow to assess UMI components. The workflow runs BarcodEX, a tool for extracting Unique Molecular Identifiers (UMIs) from single or paired-end read sequences. bwaMem is used to generate alignments from UMI fastq files with the following deduplication and QC using imported bamQC workflow

## Overview

## Dependencies

* [barcodex-rs 0.1.2](https://github.com/oicr-gsi/barcodex-rs/archive/v0.1.2.tar.gz)
* [rust 1.2](https://www.rust-lang.org/tools/install)
* [umi-tools 1.1.1](https://github.com/CGATOxford/UMI-tools/archive/1.1.1.tar.gz)
* [bwa 0.7.12](https://github.com/lh3/bwa/archive/0.7.12.tar.gz)
* [samtools 1.9](https://github.com/samtools/samtools/archive/0.1.19.tar.gz)
* [cutadapt 1.8.3](https://cutadapt.readthedocs.io/en/v1.8.3/)
* [slicer 0.3.0](https://github.com/OpenGene/slicer/archive/v0.3.0.tar.gz)
* [picard 2.21.2](https://broadinstitute.github.io/picard/command-line-overview.html)
* [python 3.6](https://www.python.org/downloads/)
* [bam-qc-metrics 0.2.5](https://github.com/oicr-gsi/bam-qc-metrics.git)
* [mosdepth 0.2.9](https://github.com/brentp/mosdepth)


## Usage

### Cromwell
```
java -jar cromwell.jar run umiQC.wdl --inputs inputs.json
```

### Inputs

#### Required workflow parameters:
Parameter|Value|Description
---|---|---
`umiList`|String|Reference ile with valid UMIs
`fastq1`|File|Fastq file for read 1
`fastq2`|File|Fastq file for read 2
`pattern1`|String|UMI pattern 1
`pattern2`|String|UMI pattern 2
`reference`|String|the reference genome for input sample
`bwaMem.readGroups`|String|Complete read group header line
`preDedupBamQC.bamQCMetrics_workflowVersion`|String|Workflow version string
`preDedupBamQC.bamQCMetrics_refSizesBed`|String|Path to human genome BED reference with chromosome sizes
`preDedupBamQC.bamQCMetrics_refFasta`|String|Path to human genome FASTA reference
`preDedupBamQC.metadata`|Map[String,String]|JSON file containing metadata
`postDedupBamQC.bamQCMetrics_workflowVersion`|String|Workflow version string
`postDedupBamQC.bamQCMetrics_refSizesBed`|String|Path to human genome BED reference with chromosome sizes
`postDedupBamQC.bamQCMetrics_refFasta`|String|Path to human genome FASTA reference
`postDedupBamQC.metadata`|Map[String,String]|JSON file containing metadata


#### Optional workflow parameters:
Parameter|Value|Default|Description
---|---|---|---
`outputPrefix`|String|"output"|Specifies the start of output files


#### Optional task parameters:
Parameter|Value|Default|Description
---|---|---|---
`extractUMIs.modules`|String|"barcodex-rs/0.1.2 rust/1.45.1"|Required environment modules
`extractUMIs.memory`|Int|24|Memory allocated for this job
`extractUMIs.timeout`|Int|12|Time in hours before task timeout
`bwaMem.adapterTrimmingLog_timeout`|Int|48|Hours before task timeout
`bwaMem.adapterTrimmingLog_jobMemory`|Int|12|Memory allocated indexing job
`bwaMem.indexBam_timeout`|Int|48|Hours before task timeout
`bwaMem.indexBam_modules`|String|"samtools/1.9"|Modules for running indexing job
`bwaMem.indexBam_jobMemory`|Int|12|Memory allocated indexing job
`bwaMem.bamMerge_timeout`|Int|72|Hours before task timeout
`bwaMem.bamMerge_modules`|String|"samtools/1.9"|Required environment modules
`bwaMem.bamMerge_jobMemory`|Int|32|Memory allocated indexing job
`bwaMem.runBwaMem_timeout`|Int|96|Hours before task timeout
`bwaMem.runBwaMem_jobMemory`|Int|32|Memory allocated for this job
`bwaMem.runBwaMem_threads`|Int|8|Requested CPU threads
`bwaMem.runBwaMem_addParam`|String?|None|Additional BWA parameters
`bwaMem.adapterTrimming_timeout`|Int|48|Hours before task timeout
`bwaMem.adapterTrimming_jobMemory`|Int|16|Memory allocated for this job
`bwaMem.adapterTrimming_addParam`|String?|None|Additional cutadapt parameters
`bwaMem.adapterTrimming_modules`|String|"cutadapt/1.8.3"|Required environment modules
`bwaMem.slicerR2_timeout`|Int|48|Hours before task timeout
`bwaMem.slicerR2_jobMemory`|Int|16|Memory allocated for this job
`bwaMem.slicerR2_modules`|String|"slicer/0.3.0"|Required environment modules
`bwaMem.slicerR1_timeout`|Int|48|Hours before task timeout
`bwaMem.slicerR1_jobMemory`|Int|16|Memory allocated for this job
`bwaMem.slicerR1_modules`|String|"slicer/0.3.0"|Required environment modules
`bwaMem.countChunkSize_timeout`|Int|48|Hours before task timeout
`bwaMem.countChunkSize_jobMemory`|Int|16|Memory allocated for this job
`bwaMem.numChunk`|Int|1|number of chunks to split fastq file [1, no splitting]
`bwaMem.doTrim`|Boolean|false|if true, adapters will be trimmed before alignment
`bwaMem.trimMinLength`|Int|1|minimum length of reads to keep [1]
`bwaMem.trimMinQuality`|Int|0|minimum quality of read ends to keep [0]
`bwaMem.adapter1`|String|"AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"|adapter sequence to trim from read 1 [AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC]
`bwaMem.adapter2`|String|"AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"|adapter sequence to trim from read 2 [AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT]
`preDedupBamQC.collateResults_timeout`|Int|1|hours before task timeout
`preDedupBamQC.collateResults_threads`|Int|4|Requested CPU threads
`preDedupBamQC.collateResults_jobMemory`|Int|8|Memory allocated for this job
`preDedupBamQC.collateResults_modules`|String|"python/3.6"|required environment modules
`preDedupBamQC.cumulativeDistToHistogram_timeout`|Int|1|hours before task timeout
`preDedupBamQC.cumulativeDistToHistogram_threads`|Int|4|Requested CPU threads
`preDedupBamQC.cumulativeDistToHistogram_jobMemory`|Int|8|Memory allocated for this job
`preDedupBamQC.cumulativeDistToHistogram_modules`|String|"python/3.6"|required environment modules
`preDedupBamQC.runMosdepth_timeout`|Int|4|hours before task timeout
`preDedupBamQC.runMosdepth_threads`|Int|4|Requested CPU threads
`preDedupBamQC.runMosdepth_jobMemory`|Int|16|Memory allocated for this job
`preDedupBamQC.runMosdepth_modules`|String|"mosdepth/0.2.9"|required environment modules
`preDedupBamQC.bamQCMetrics_timeout`|Int|4|hours before task timeout
`preDedupBamQC.bamQCMetrics_threads`|Int|4|Requested CPU threads
`preDedupBamQC.bamQCMetrics_jobMemory`|Int|16|Memory allocated for this job
`preDedupBamQC.bamQCMetrics_modules`|String|"bam-qc-metrics/0.2.5"|required environment modules
`preDedupBamQC.bamQCMetrics_normalInsertMax`|Int|1500|Maximum of expected insert size range
`preDedupBamQC.markDuplicates_timeout`|Int|4|hours before task timeout
`preDedupBamQC.markDuplicates_threads`|Int|4|Requested CPU threads
`preDedupBamQC.markDuplicates_jobMemory`|Int|16|Memory allocated for this job
`preDedupBamQC.markDuplicates_modules`|String|"picard/2.21.2"|required environment modules
`preDedupBamQC.markDuplicates_picardMaxMemMb`|Int|6000|Memory requirement in MB for running Picard JAR
`preDedupBamQC.markDuplicates_opticalDuplicatePixelDistance`|Int|100|Maximum offset between optical duplicate clusters
`preDedupBamQC.downsampleRegion_timeout`|Int|4|hours before task timeout
`preDedupBamQC.downsampleRegion_threads`|Int|4|Requested CPU threads
`preDedupBamQC.downsampleRegion_jobMemory`|Int|16|Memory allocated for this job
`preDedupBamQC.downsampleRegion_modules`|String|"samtools/1.9"|required environment modules
`preDedupBamQC.downsample_timeout`|Int|4|hours before task timeout
`preDedupBamQC.downsample_threads`|Int|4|Requested CPU threads
`preDedupBamQC.downsample_jobMemory`|Int|16|Memory allocated for this job
`preDedupBamQC.downsample_modules`|String|"samtools/1.9"|required environment modules
`preDedupBamQC.downsample_randomSeed`|Int|42|Random seed for pre-downsampling (if any)
`preDedupBamQC.downsample_downsampleSuffix`|String|"downsampled.bam"|Suffix for output file
`preDedupBamQC.findDownsampleParamsMarkDup_timeout`|Int|4|hours before task timeout
`preDedupBamQC.findDownsampleParamsMarkDup_threads`|Int|4|Requested CPU threads
`preDedupBamQC.findDownsampleParamsMarkDup_jobMemory`|Int|16|Memory allocated for this job
`preDedupBamQC.findDownsampleParamsMarkDup_modules`|String|"python/3.6"|required environment modules
`preDedupBamQC.findDownsampleParamsMarkDup_customRegions`|String|""|Custom downsample regions; overrides chromosome and interval parameters
`preDedupBamQC.findDownsampleParamsMarkDup_intervalStart`|Int|100000|Start of interval in each chromosome, for very large BAMs
`preDedupBamQC.findDownsampleParamsMarkDup_baseInterval`|Int|15000|Base width of interval in each chromosome, for very large BAMs
`preDedupBamQC.findDownsampleParamsMarkDup_chromosomes`|Array[String]|["chr12", "chr13", "chrXII", "chrXIII"]|Array of chromosome identifiers for downsampled subset
`preDedupBamQC.findDownsampleParamsMarkDup_threshold`|Int|10000000|Minimum number of reads to conduct downsampling
`preDedupBamQC.findDownsampleParams_timeout`|Int|4|hours before task timeout
`preDedupBamQC.findDownsampleParams_threads`|Int|4|Requested CPU threads
`preDedupBamQC.findDownsampleParams_jobMemory`|Int|16|Memory allocated for this job
`preDedupBamQC.findDownsampleParams_modules`|String|"python/3.6"|required environment modules
`preDedupBamQC.findDownsampleParams_preDSMultiplier`|Float|1.5|Determines target size for pre-downsampled set (if any). Must have (preDSMultiplier) < (minReadsRelative).
`preDedupBamQC.findDownsampleParams_precision`|Int|8|Number of decimal places in fraction for pre-downsampling
`preDedupBamQC.findDownsampleParams_minReadsRelative`|Int|2|Minimum value of (inputReads)/(targetReads) to allow pre-downsampling
`preDedupBamQC.findDownsampleParams_minReadsAbsolute`|Int|10000|Minimum value of targetReads to allow pre-downsampling
`preDedupBamQC.findDownsampleParams_targetReads`|Int|100000|Desired number of reads in downsampled output
`preDedupBamQC.indexBamFile_timeout`|Int|4|hours before task timeout
`preDedupBamQC.indexBamFile_threads`|Int|4|Requested CPU threads
`preDedupBamQC.indexBamFile_jobMemory`|Int|16|Memory allocated for this job
`preDedupBamQC.indexBamFile_modules`|String|"samtools/1.9"|required environment modules
`preDedupBamQC.countInputReads_timeout`|Int|4|hours before task timeout
`preDedupBamQC.countInputReads_threads`|Int|4|Requested CPU threads
`preDedupBamQC.countInputReads_jobMemory`|Int|16|Memory allocated for this job
`preDedupBamQC.countInputReads_modules`|String|"samtools/1.9"|required environment modules
`preDedupBamQC.updateMetadata_timeout`|Int|4|hours before task timeout
`preDedupBamQC.updateMetadata_threads`|Int|4|Requested CPU threads
`preDedupBamQC.updateMetadata_jobMemory`|Int|16|Memory allocated for this job
`preDedupBamQC.updateMetadata_modules`|String|"python/3.6"|required environment modules
`preDedupBamQC.filter_timeout`|Int|4|hours before task timeout
`preDedupBamQC.filter_threads`|Int|4|Requested CPU threads
`preDedupBamQC.filter_jobMemory`|Int|16|Memory allocated for this job
`preDedupBamQC.filter_modules`|String|"samtools/1.9"|required environment modules
`preDedupBamQC.filter_minQuality`|Int|30|Minimum alignment quality to pass filter
`bamSplitDeduplication.modules`|String|"umi-tools/1.0.0 samtools/1.9"|Required environment modules
`bamSplitDeduplication.memory`|Int|24|Memory allocated for this job
`bamSplitDeduplication.timeout`|Int|6|Time in hours before task timeout
`bamMerge.modules`|String|"samtools/1.9"|Required environment modules
`bamMerge.memory`|Int|24|Memory allocated for indexing job
`bamMerge.timeout`|Int|6|Hours before task timeout
`postDedupBamQC.collateResults_timeout`|Int|1|hours before task timeout
`postDedupBamQC.collateResults_threads`|Int|4|Requested CPU threads
`postDedupBamQC.collateResults_jobMemory`|Int|8|Memory allocated for this job
`postDedupBamQC.collateResults_modules`|String|"python/3.6"|required environment modules
`postDedupBamQC.cumulativeDistToHistogram_timeout`|Int|1|hours before task timeout
`postDedupBamQC.cumulativeDistToHistogram_threads`|Int|4|Requested CPU threads
`postDedupBamQC.cumulativeDistToHistogram_jobMemory`|Int|8|Memory allocated for this job
`postDedupBamQC.cumulativeDistToHistogram_modules`|String|"python/3.6"|required environment modules
`postDedupBamQC.runMosdepth_timeout`|Int|4|hours before task timeout
`postDedupBamQC.runMosdepth_threads`|Int|4|Requested CPU threads
`postDedupBamQC.runMosdepth_jobMemory`|Int|16|Memory allocated for this job
`postDedupBamQC.runMosdepth_modules`|String|"mosdepth/0.2.9"|required environment modules
`postDedupBamQC.bamQCMetrics_timeout`|Int|4|hours before task timeout
`postDedupBamQC.bamQCMetrics_threads`|Int|4|Requested CPU threads
`postDedupBamQC.bamQCMetrics_jobMemory`|Int|16|Memory allocated for this job
`postDedupBamQC.bamQCMetrics_modules`|String|"bam-qc-metrics/0.2.5"|required environment modules
`postDedupBamQC.bamQCMetrics_normalInsertMax`|Int|1500|Maximum of expected insert size range
`postDedupBamQC.markDuplicates_timeout`|Int|4|hours before task timeout
`postDedupBamQC.markDuplicates_threads`|Int|4|Requested CPU threads
`postDedupBamQC.markDuplicates_jobMemory`|Int|16|Memory allocated for this job
`postDedupBamQC.markDuplicates_modules`|String|"picard/2.21.2"|required environment modules
`postDedupBamQC.markDuplicates_picardMaxMemMb`|Int|6000|Memory requirement in MB for running Picard JAR
`postDedupBamQC.markDuplicates_opticalDuplicatePixelDistance`|Int|100|Maximum offset between optical duplicate clusters
`postDedupBamQC.downsampleRegion_timeout`|Int|4|hours before task timeout
`postDedupBamQC.downsampleRegion_threads`|Int|4|Requested CPU threads
`postDedupBamQC.downsampleRegion_jobMemory`|Int|16|Memory allocated for this job
`postDedupBamQC.downsampleRegion_modules`|String|"samtools/1.9"|required environment modules
`postDedupBamQC.downsample_timeout`|Int|4|hours before task timeout
`postDedupBamQC.downsample_threads`|Int|4|Requested CPU threads
`postDedupBamQC.downsample_jobMemory`|Int|16|Memory allocated for this job
`postDedupBamQC.downsample_modules`|String|"samtools/1.9"|required environment modules
`postDedupBamQC.downsample_randomSeed`|Int|42|Random seed for pre-downsampling (if any)
`postDedupBamQC.downsample_downsampleSuffix`|String|"downsampled.bam"|Suffix for output file
`postDedupBamQC.findDownsampleParamsMarkDup_timeout`|Int|4|hours before task timeout
`postDedupBamQC.findDownsampleParamsMarkDup_threads`|Int|4|Requested CPU threads
`postDedupBamQC.findDownsampleParamsMarkDup_jobMemory`|Int|16|Memory allocated for this job
`postDedupBamQC.findDownsampleParamsMarkDup_modules`|String|"python/3.6"|required environment modules
`postDedupBamQC.findDownsampleParamsMarkDup_customRegions`|String|""|Custom downsample regions; overrides chromosome and interval parameters
`postDedupBamQC.findDownsampleParamsMarkDup_intervalStart`|Int|100000|Start of interval in each chromosome, for very large BAMs
`postDedupBamQC.findDownsampleParamsMarkDup_baseInterval`|Int|15000|Base width of interval in each chromosome, for very large BAMs
`postDedupBamQC.findDownsampleParamsMarkDup_chromosomes`|Array[String]|["chr12", "chr13", "chrXII", "chrXIII"]|Array of chromosome identifiers for downsampled subset
`postDedupBamQC.findDownsampleParamsMarkDup_threshold`|Int|10000000|Minimum number of reads to conduct downsampling
`postDedupBamQC.findDownsampleParams_timeout`|Int|4|hours before task timeout
`postDedupBamQC.findDownsampleParams_threads`|Int|4|Requested CPU threads
`postDedupBamQC.findDownsampleParams_jobMemory`|Int|16|Memory allocated for this job
`postDedupBamQC.findDownsampleParams_modules`|String|"python/3.6"|required environment modules
`postDedupBamQC.findDownsampleParams_preDSMultiplier`|Float|1.5|Determines target size for pre-downsampled set (if any). Must have (preDSMultiplier) < (minReadsRelative).
`postDedupBamQC.findDownsampleParams_precision`|Int|8|Number of decimal places in fraction for pre-downsampling
`postDedupBamQC.findDownsampleParams_minReadsRelative`|Int|2|Minimum value of (inputReads)/(targetReads) to allow pre-downsampling
`postDedupBamQC.findDownsampleParams_minReadsAbsolute`|Int|10000|Minimum value of targetReads to allow pre-downsampling
`postDedupBamQC.findDownsampleParams_targetReads`|Int|100000|Desired number of reads in downsampled output
`postDedupBamQC.indexBamFile_timeout`|Int|4|hours before task timeout
`postDedupBamQC.indexBamFile_threads`|Int|4|Requested CPU threads
`postDedupBamQC.indexBamFile_jobMemory`|Int|16|Memory allocated for this job
`postDedupBamQC.indexBamFile_modules`|String|"samtools/1.9"|required environment modules
`postDedupBamQC.countInputReads_timeout`|Int|4|hours before task timeout
`postDedupBamQC.countInputReads_threads`|Int|4|Requested CPU threads
`postDedupBamQC.countInputReads_jobMemory`|Int|16|Memory allocated for this job
`postDedupBamQC.countInputReads_modules`|String|"samtools/1.9"|required environment modules
`postDedupBamQC.updateMetadata_timeout`|Int|4|hours before task timeout
`postDedupBamQC.updateMetadata_threads`|Int|4|Requested CPU threads
`postDedupBamQC.updateMetadata_jobMemory`|Int|16|Memory allocated for this job
`postDedupBamQC.updateMetadata_modules`|String|"python/3.6"|required environment modules
`postDedupBamQC.filter_timeout`|Int|4|hours before task timeout
`postDedupBamQC.filter_threads`|Int|4|Requested CPU threads
`postDedupBamQC.filter_jobMemory`|Int|16|Memory allocated for this job
`postDedupBamQC.filter_modules`|String|"samtools/1.9"|required environment modules
`postDedupBamQC.filter_minQuality`|Int|30|Minimum alignment quality to pass filter
`statsMerge.memory`|Int|16|Memory allocated for this job


### Outputs

Output | Type | Description | Labels
---|---|---|---
`umiCounts`|File|JSON record of UMI counts after extraction|vidarr_label: umiCounts
`extractionMetrics`|File|JSON of metrics relating to extraction process|vidarr_label: extractionMetrics
`umiCountsPerPosition`|File|tsv file tabulates the counts for unique combinations of UMI and position|vidarr_label: umiCountsPerPosition
`preDedupBamMetrics`|File|BamQC JSON report on bam file pre-deduplication|vidarr_label: preDedupBamMetrics
`postDedupBamMetrics`|File|BamQC JSON report on bam file post-deduplication|vidarr_label: postDedupBamMetrics


## Commands
This section lists command(s) run by umiQC workflow
 
* Running umiQC
 
QC workflow to assess UMI components.
 
### Get lengths of paired-end UMIs from kit
 
```
 
       k=($(awk '{ match($1, "([ACTG])+"); print RLENGTH }' ~{umiList} | uniq))
 
       for i in ${k[@]}
       do
           for j in ${k[@]}
           do
               # adding 1 to account for period in UMI
               # 'ATG.ATCG'
               L+=($(($i+$j+1)))
           done
       done
 
       umiLengths=($(tr ' ' '\n' <<< "${L[@]}" | awk '!u[$0]++' | tr ' ' '\n'))
       printf "%s\n" "${umiLengths[@]}"
 
```
 
### Extracting UMIs from FASTQ files and sorting UMI counts json
 
```
             barcodex-rs --umilist ~{umiList} --prefix ~{outputPrefix} --separator "__" inline \
             --pattern1 '~{pattern1}' --r1-in ~{fastq1} \
             --pattern2 '~{pattern2}' --r2-in ~{fastq2} 

             cat ~{outputPrefix}_UMI_counts.json > umiCounts.txt

            tr [,] ',\n' < umiCounts.txt | sed 's/[{}]//' > tmp.txt
            echo "{$(sort -i tmp.txt)}" > new.txt
            tr '\n' ',' < new.txt | sed 's/,$//' > ~{outputPrefix}_UMI_counts.json
```
 
### Splitting and deduplicating BAM files based on UMI lengths
 
```
         samtools view -H ~{bamFile} > ~{outputPrefix}.~{umiLength}.sam
         samtools view ~{bamFile} | grep -P "^.*__\S{~{umiLength}}\t" >> ~{outputPrefix}.~{umiLength}.sam
         samtools view -Sb ~{outputPrefix}.~{umiLength}.sam > ~{outputPrefix}.~{umiLength}.bam
 
         samtools index ~{outputPrefix}.~{umiLength}.bam
 
         umi_tools group -I ~{outputPrefix}.~{umiLength}.bam \
         --group-out=~{outputPrefix}.~{umiLength}.umi_groups.tsv \
         --output-bam > ~{outputPrefix}.~{umiLength}.dedup.bam \
         --log=group.log --paired | samtools view
```
 
### Merge UMI metrics into one TSV file
 
```
         umiMetrics=(~{sep=" " umiMetrics})
         length=${#umiMetrics[@]}
 
         i=0
 
         awk 'NR==1' ${umiMetrics[i]} > mergedUMIMetrics.tsv
         while [ $i -le $length ]
         do
             gawk -i inplace '(NR>1) { match($7, "([ACTG.])+") ;  $9=RLENGTH-1"."$9 ; print}' ${umiMetrics[i]}
     
             cat ${umiMetrics[i]} >> mergedUMIMetrics.tsv
             i=$(( $i+1 ))
         done
         tr -s '[ ,     ]' '\t' < mergedUMIMetrics.tsv > tmp.tsv && mv tmp.tsv mergedUMIMetrics.tsv 
```
### Merge deduplicated BAM files
 
```        
         set -euo pipefail
         samtools merge -c ~{outputPrefix}.dedup.bam ~{sep=" " umiDedupBams}
``` 
## Support

For support, please file an issue on the [Github project](https://github.com/oicr-gsi) or send an email to gsi@oicr.on.ca .
