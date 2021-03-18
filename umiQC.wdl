version 1.0

import "imports/pull_bwaMem.wdl" as bwaMem
import "imports/pull_bamQC.wdl" as bamQC

workflow umiQC {
    input {
        File umiList
        String outputPrefix = "output"
        File fastq1
        File fastq2
    }

    parameter_meta {
        umiList: "File with valid UMIs"
        outputPrefix: "Specifies the start of output files"
        fastq1: "Fastq file for read 1"
        fastq2: "Fastq file for read 2"
    }

    meta {
        author: "Michelle Feng"
        email: "mfeng@oicr.on.ca"
        description: "QC workflow to assess UMI components"
        dependencies: [
            {
                name: "barcodex-rs/0.1.2",
                url: "https://github.com/oicr-gsi/barcodex-rs/archive/v0.1.2.tar.gz"
            },
            {
                name: "rust/1.2",
                url: "https://www.rust-lang.org/tools/install"
            },
            {
                name: "umi-tools/1.1.1",
                url: "https://github.com/CGATOxford/UMI-tools/archive/1.1.1.tar.gz"
            },
            {
            name: "bwa/0.7.12",
            url: "https://github.com/lh3/bwa/archive/0.7.12.tar.gz"
            },
            {
                name: "samtools/1.9",
                url: "https://github.com/samtools/samtools/archive/0.1.19.tar.gz"
            },
            {
                name: "cutadapt/1.8.3",
                url: "https://cutadapt.readthedocs.io/en/v1.8.3/"
            },
            {
                name: "slicer/0.3.0",
                url: "https://github.com/OpenGene/slicer/archive/v0.3.0.tar.gz"
            },
            {
                name: "picard/2.21.2",
                url: "https://broadinstitute.github.io/picard/command-line-overview.html"
            },
            {
                name: "python/3.6",
                url: "https://www.python.org/downloads/"
            },
            {
                name: "bam-qc-metrics/0.2.5",
                url: "https://github.com/oicr-gsi/bam-qc-metrics.git"
            },
            {
                name: "mosdepth/0.2.9",
                url: "https://github.com/brentp/mosdepth"
            }
        ]
        output_meta: {
            umiCounts: "Record of UMI counts after extraction",
            extractionMetrics: "Metrics relating to extraction process",
            preDedupBamMetrics: "BamQC report on bam file pre-deduplication",
            umiMetricsSix: "File mapping read id to read group",
            umiMetricsSeven: "File mapping read id to read group",
            umiMetricsEight: "File mapping read id to read group",
            postDedupBamMetrics: "BamQC report on bam file post-deduplication"
        }
    }

    call extractUMIs { 
        input:
            umiList = umiList,
            outputPrefix = outputPrefix,
            fastq1 = fastq1,
            fastq2 = fastq2
    }

    call bwaMem.bwaMem {
        input:
            fastqR1 = extractUMIs.fastqR1,
            fastqR2 = extractUMIs.fastqR2,
            outputFileNamePrefix = outputPrefix
    }

    call bamQC.bamQC as preDedupBamQC {
        input:
            bamFile = bwaMem.bwaMemBam,
            outputFileNamePrefix = "preDedup"
    }

    call bamSplit {
        input:
            bamFile = bwaMem.bwaMemBam,
            outputPrefix = outputPrefix
    }

    call umiDeduplications {
        input:
            outputPrefix = outputPrefix,
            bamFiles = bamSplit.bamFiles
    }

    call bamQC.bamQC as postDedupBamQC {
        input:
            bamFile = umiDeduplications.umiDedupBam,
            outputFileNamePrefix = "postDedup"
    }

    output {
        # barcodex metrics
        File umiCounts = extractUMIs.umiCounts
        File extractionMetrics = extractUMIs.extractionMetrics

        # pre-collapse bamqc metrics
        File preDedupBamMetrics = preDedupBamQC.result

        # umi-tools metrics
        File umiMetricsSix = umiDeduplications.umiMetricsSix
        File umiMetricsSeven = umiDeduplications.umiMetricsSeven
        File umiMetricsEight = umiDeduplications.umiMetricsEight

        # post-collapse bamqc metrics
        File postDedupBamMetrics = postDedupBamQC.result
    } 
}

task extractUMIs {
    input {
        File umiList
        String outputPrefix
        File fastq1
        File fastq2
        String modules = "barcodex-rs/0.1.2 rust/1.45.1"
        Int memory = 24
        Int timeout = 12
    }

    parameter_meta {
        umiList: "File with valid UMIs"
        outputPrefix: "Specifies the start of the output files"
        fastqR1: "FASTQ file containing read 1"
        fastqR2: "FASTQ file containing read 2"
        modules: "Required environment modules"
        memory: "Memory allocated for this job"
        timeout: "Time in hours before task timeout"
    }

    command <<<
        barcodex-rs --umilist ~{umiList} --prefix ~{outputPrefix} --separator "__" inline \
        --pattern1 "(?P<umi_1>^[ACGT]{3}[ACG])(?P<discard_1>T)|(?P<umi_2>^[ACGT]{3})(?P<discard_2>T)" --r1-in ~{fastq1} \
        --pattern2 "(?P<umi_1>^[ACGT]{3}[ACG])(?P<discard_1>T)|(?P<umi_2>^[ACGT]{3})(?P<discard_2>T)" --r2-in ~{fastq2} 
    >>>

    runtime {
        modules: "~{modules}"
        memory: "~{memory}G"
        timeout: "~{timeout}"
    }

    output {
        File fastqR1 = "~{outputPrefix}_R1.fastq.gz"
        File fastqR2 = "~{outputPrefix}_R2.fastq.gz"
        File discardR1 = "~{outputPrefix}_R1.discarded.fastq.gz"
        File discardR2 = "~{outputPrefix}_R2.discarded.fastq.gz"
        File extractR1 = "~{outputPrefix}_R1.extracted.fastq.gz"
        File extractR2 = "~{outputPrefix}_R2.extracted.fastq.gz"
        File umiCounts = "~{outputPrefix}_UMI_counts.json"
        File extractionMetrics = "~{outputPrefix}_extraction_metrics.json"
    }

    meta {
        output_meta: {
            fastqR1: "Read 1 fastq file with UMIs extracted",
            fastqR2: "Read 2 fastq file with UMIs extracted",
            discardR1: "Reads without a matching UMI pattern in read 1",
            discardR2: "Reads without a matching UMI pattern in read 2",
            extractR1: "Extracted reads (UMIs and any spacer sequences) from read 1",
            extractR2: "Extracted reads (UMIs and any spacer sequences) from read 2",
            umiCounts: "Record of UMI counts after extraction",
            extractionMetrics: "Metrics relating to extraction process"
        }
    }
}

task bamSplit {
    input {
        File bamFile
        String outputPrefix
        Int minLength = 3
        Int maxLength = 4
        String modules = "samtools/1.9"
        Int memory = 24
        Int timeout = 6
    }

    parameter_meta {
        bamFile: "Bam file from bwaMem containing UMIs of varying lengths"
        outputPrefix: "Specifies the start of the output files"
        minLength: "Minimum length of barcode from barcode list"
        maxLength: "Maximum length of barcode from barcode list"
        modules: "Required environment modules"
        memory: "Memory allocated for this job"
        timeout: "Time in hours before task timeout"
    }

    command <<<
        samtools view -H ~{bamFile} > ~{outputPrefix}.~{minLength * 2}.sam
        samtools view ~{bamFile} \
        | grep -P "^.*__\[ACGT]{~{minLength}}\.\[ACGT]{~{minLength}}\t" \
        >> ~{outputPrefix}.~{minLength * 2}.sam
        samtools view -Sb ~{outputPrefix}.~{minLength * 2}.sam \
        > ~{outputPrefix}.~{minLength * 2}.bam

        samtools view -H ~{bamFile} > ~{outputPrefix}.~{minLength + maxLength}.1.sam
        samtools view ~{bamFile} \
        | grep -P "^.*__\[ACGT]{~{minLength}}\.\[ACGT]{~{maxLength}}\t" \
        >> ~{outputPrefix}.~{minLength + maxLength}.1.sam
        samtools view -Sb ~{outputPrefix}.~{minLength + maxLength}.1.sam \
        > ~{outputPrefix}.~{minLength + maxLength}.1.bam

        samtools view -H ~{bamFile} > ~{outputPrefix}.~{minLength + maxLength}.2.sam
        samtools view ~{bamFile} \
        | grep -P "^.*__\[ACGT]{~{maxLength}}\.\[ACGT]{~{minLength}}\t" \
        >> ~{outputPrefix}.~{minLength + maxLength}.2.sam
        samtools view -Sb ~{outputPrefix}.~{minLength + maxLength}.2.sam \
        > ~{outputPrefix}.~{minLength + maxLength}.2.bam

        samtools view -H ~{bamFile} > ~{outputPrefix}.~{maxLength * 2}.sam
        samtools view ~{bamFile} \
        | grep -P "^.*__\[ACGT]{~{maxLength}}\.\[ACGT]{~{maxLength}}\t" \
        >> ~{outputPrefix}.~{maxLength * 2}.sam
        samtools view -Sb ~{outputPrefix}.~{maxLength * 2}.sam \
        > ~{outputPrefix}.~{maxLength * 2}.bam
    >>>

    runtime {
        modules: "~{modules}"
        memory: "~{memory}G"
        timeout: "~{timeout}"
    }

    output {
        File outputSix = "~{outputPrefix}.~{minLength * 2}.bam"
        File outputSevenOne = "~{outputPrefix}.~{minLength + maxLength}.1.bam"
        File outputSevenTwo = "~{outputPrefix}.~{minLength + maxLength}.2.bam"
        File outputEight = "~{outputPrefix}.~{maxLength * 2}.bam"
        Array[File] bamFiles = glob("*.bam")
    }

    meta {
        output_meta: {
            outputSix: "UMIs with total barcode length six",
            outputSevenOne: "UMIs with barcode length three + four",
            outputSevenTwo: "UMIs with barcode length four + three",
            outputEight: "UMIs with total barcode length eight",
            bamFiles: "Array of BAMs with varying lengths of UMIs"
        }
    }
}

task umiDeduplications {
    input {
        Array[File] bamFiles
        Int minLength = 3
        Int maxLength = 4
        String outputPrefix
        String modules = "umi-tools/1.0.0 samtools/1.9"
        Int memory = 24
        Int timeout = 6
    }

    parameter_meta {
        bamFiles: "Array of BAM files with varying lengths of UMIs"
        minLength: "Minimum length of barcode from barcode list"
        maxLength: "Maximum length of barcode from barcode list"
        outputPrefix: "Specifies the start of the output files"
        modules: "Required environment modules"
        memory: "Memory allocated for this job"
        timeout: "Time in hours before task timeout"
    }

    command <<<
        for x in ~{sep=' ' bamFiles}
        do
            samtools index "${x}"
            umi_tools group -I "${x}" \
            --group-out=$(basename "${x}" .bam).umi_groups.tsv \
            --output-bam > $(basename "${x}" .bam).dedup.bam \
            --log=group.log --paired | samtools view
        done;

        samtools merge ~{outputPrefix}.dedup.bam \
        ~{outputPrefix}.~{minLength * 2}.dedup.bam \
        ~{outputPrefix}.~{minLength + maxLength}.dedup.bam \
        ~{outputPrefix}.~{maxLength * 2}.dedup.bam
    >>>

    runtime {
        modules: "~{modules}"
        memory: "~{memory}G"
        timeout: "~{timeout}"
    }

    output {
        File umiDedupBam = "~{outputPrefix}.dedup.bam"
        File umiMetricsSix = "~{outputPrefix}.~{minLength * 2}.umi_groups.tsv"
        File umiMetricsSeven = "~{outputPrefix}.~{minLength + maxLength}.umi_groups.tsv"
        File umiMetricsEight = "~{outputPrefix}.~{maxLength * 2}.umi_groups.tsv"
    }

    meta {
        output_meta: {
            umiDedupBam: "Deduplicated bam file",
            umiMetricsSix: "File mapping read id to read group",
            umiMetricsSeven: "File mapping read id to read group",
            umiMetricsEight: "File mapping read id to read group"
        }
    }
}
