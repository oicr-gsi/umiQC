version 1.0

import "imports/pull_bwaMem.wdl" as bwaMem
import "imports/pull_bamQC.wdl" as bamQC

workflow umiQC {
    input {
    }

    parameter_meta {
        
    }

    meta {
        author: "Michelle Feng"
        email: "mfeng@oicr.on.ca"
        description: ""
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

        }
    }

    call extractUMIs { 
        input:
            umiList = umiList,
            prefix = prefix,
            fastqR1 = fastqR1,
            fastqR2 = fastqR2,
            modules = modules,
            memory = memory,
            timeout = timeout
    }

    call bwaMem.bwaMem {
        input:
            fastqR1 = extractUMIs.fastqR1,
            fastqR2 = extractUMIs.fastqR2,
            outputFileNamePrefix
    }

    call bamQC.bamQC {
        input:
            bamFile = bwaMem.bwaMemBam,
            outputFileNamePrefix = outputFileNamePrefix
    }

    call umiDeduplication {
        input:
            bamFile = bamQC.bamQC.preDedupBamQCOutput,
            prefix = prefix
    }

    call bamQC.bamQC {
        input:
            bamFile = umiDeduplication.umiDedupBam,
            outputFileNamePrefix = outputFileNamePrefix
    }

    output {
        # barcodex metrics
        File umiCounts = extractUMIs.umiCounts
        File extractionMetrics = extractUMIs.extractionMetrics

        # pre-collapse bamqc metrics
        File preDedupBamMetrics = bamQC.preDedupResult

        # umi-tools metrics
        File umiDedupBam = umiDeduplication.umiDedupBam
        File umiGroups = umiDeduplication.umiGroups

        # post-collapse bamqc metrics
        File postDedupBamMetrics = bamQC.postDedupResult
    }

    task extractUMIs {
        input {
            File umiList
            String prefix
            File fastqR1
            File fastqR2
            String modules = "barcodex-rs/0.1.2 rust/1.45.1"
            Int memory = 24
            Int timeout = 12
        }

        parameter_meta {
            umiList: "File with valid UMIs"
            prefix: "Specifies the start of the output files"
            fastqR1: "FASTQ file containing read 1"
            fastqR2: "FASTQ file containing read 2"
            modules: "Required environment modules"
            memory: "Memory allocated for this job"
            timeout: "Time in hours before task timeout"
        }

        command <<<
            barcodex-rs --umilist ~{umiList} --prefix ~{prefix} --separator "__" inline \
            --pattern1 "(?P<umi_1>^[ACGT]{3}[ACG])(?P<discard_1>T)|(?P<umi_2>^[ACGT]{3})(?P<discard_2>T)" --r1-in ~{fastqR1} \
            --pattern2 "(?P<umi_1>^[ACGT]{3}[ACG])(?P<discard_1>T)|(?P<umi_2>^[ACGT]{3})(?P<discard_2>T)" --r2-in ~{fastqR2} 
        >>>

        runtime {
            modules: "~{modules}"
            memory: "~{memory}G"
            timeout: "~{timeout}"
        }

        output {
            File fastqR1 = "~{prefix}_R1.fastq.gz"
            File fastqR2 = "~{prefix}_R2.fastq.gz"
            File discardR1 = "~{prefix}_R1.discarded.fastq.gz"
            File discardR2 = "~{prefix}_R2.discarded.fastq.gz"
            File extractR1 = "~{prefix}_R1.extractedfastq.gz"
            File extractR2 = "~{prefix}_R2.extractedfastq.gz"
            File umiCounts = "~{prefix}_UMI_counts.json"
            File extractionMetrics = "~{prefix}_extraction_metrics.json"
        }

        meta {
            output_meta: {
                fastqR1: "Read 1 fastq file with UMIs extracted"
                fastqR2: "Read 2 fastq file with UMIs extracted"
                discardR1: "Reads without a matching UMI pattern in read 1"
                discardR2: "Reads without a matching UMI pattern in read 2"
                extractR1: "Extracted reads (UMIs and any spacer sequences) from read 1"
                extractR2: "Extracted reads (UMIs and any spacer sequences) from read 2"
                umiCounts: "Record of UMI counts after extraction"
                extractionMetrics: "Metrics relating to extraction process"
            }
        }
    }

    task umiDeduplication {
        input {
            File bamFile
            File prefix
            String modules = "umi-tools/1.0.0"
            Int memory = 24
            Int timeout = 6
        }

        parameter_meta {
            modules: "Required environment modules"
            memory: "Memory allocated for this job"
            timeout: "Time in hours before task timeout"
        }

        command <<<
            umi_tools group -I ~{bamFile} \
            --group-out=~{prefix}.umi_groups.tsv \
            --output-bam > ~{prefix}.dedup.bam
        >>>

        runtime {
            memory: "~{memory}G"
            timeout: "~{timeout}"
        }

        output {
            File umiDedupBam = "~{prefix}.dedup.bam"
            File umiGroups = ~{prefix}.umi_groups.tsv
        }

        meta {
            output_meta: {
                umiDedupBam: "Deduplicated bam file"
                umiGroups: "File mapping read id to read group"
            }
        }
    }
}
