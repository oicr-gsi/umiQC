version 1.0
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
            }
        ]
        output_meta: {

        }
    }

    task extractUMIs {
        input {
            File umiList
            String prefix
            File fastqR1
            File fastqR2
            String separator
            String modules = "barcodex-rs/0.1.2 rust/1.45.1"
            Int memory = 24
            Int timeout = 12
        }

        parameter_meta {
            umiList: "File with valid UMIs"
            prefix: "Specifies the start of the output files"
            fastqR1: "FASTQ file containing read 1"
            fastqR2: "FASTQ file containing read 2"
            separator: "String separating the UMI sequence in the read name"
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
            timeout: "~{timeout}"
            memory: "~{memory}G"
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
                fastqR1: "Read 1 fastq file with UMI extracted"
                fastqR2: "Read 2 fastq file with UMI extraced"
                discardR1: "Portion of read 1 with no UMIs found"
                discardR2: "Portion of read 2 with no UMIs found"
                extractR1: "Portion of read 1 extracted"
                extractR2: "Portion of read 2 extracted"
                umiCounts: ""
                extractionMetrics: ""
            }
        }
    }

    task alignReads {
        input {
            File fastqR1
            File fastqR2
            String outputFileNamePrefix
        }

        parameter_meta {
        }

        command <<<
        >>>

        runtime {
        }

        output {

        }

        meta {
            output_meta: {

            }
        }
    }

    task precollapseBamQC {
        input {
        }

        parameter_meta {
        }

        command <<<
        >>>

        runtime {
        }

        output {

        }

        meta {
            output_meta: {

            }
        }
    }

    task umiDeduplication {
        input {
            String modules = "umi-tools/1.0.0"
        }

        parameter_meta {
            modules: "Required environment modules"
        }

        command <<<
            umi_tools group -I 
            ~{id}.bam --group-out=~{id}

            .umi_groups.tsv --output-bam
            ~{id}

            .dedup.bam
        >>>

        runtime {
        }

        output {

        }

        meta {
            output_meta: {

            }
        }
    }

    task postcollapseBamQC {
        input {
        }

        parameter_meta {
        }

        command <<<
        >>>

        runtime {
        }

        output {

        }

        meta {
            output_meta: {

            }
        }
    }
}
