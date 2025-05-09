version 1.0

import "imports/pull_bwaMem.wdl" as bwaMem
import "imports/pull_bamQC.wdl" as bamQC

struct GenomeResources {
    String bwaRef
    String runBwaMem_modules
}

workflow umiQC {
    input {
        String umiList
        String outputPrefix = "output"
        File fastq1
        File fastq2
        String pattern1 
        String pattern2
        String reference
    }

    Map[String,GenomeResources] resources = {
    "hg19": {
        "bwaRef": "$HG19_BWA_INDEX_ROOT/hg19_random.fa",
        "runBwaMem_modules": "samtools/1.9 bwa/0.7.12 hg19-bwa-index/0.7.12"
    },
    "hg38": {
        "bwaRef": "$HG38_BWA_INDEX_ROOT/hg38_random.fa",
        "runBwaMem_modules": "samtools/1.9 bwa/0.7.12 hg38-bwa-index/0.7.12"
    },
    "mm10":{
        "bwaRef": "$MM10_BWA_INDEX_ROOT/mm10.fa",
        "runBwaMem_modules": "samtools/1.9 bwa/0.7.12 mm10-bwa-index/0.7.12"
    }
    }

    parameter_meta {
        umiList: "Reference ile with valid UMIs"
        outputPrefix: "Specifies the start of output files"
        fastq1: "Fastq file for read 1"
        fastq2: "Fastq file for read 2"
        pattern1: "UMI pattern 1"
        pattern2: "UMI pattern 2"
        reference : "the reference genome for input sample"
    }

    meta {
        author: "Michelle Feng, Murto Hilali and Gavin Peng"
        email: "mfeng@oicr.on.ca, mhilali@oicr.on.ca and gpeng@oicr.on.ca"
        description: "QC workflow to assess UMI components. The workflow runs BarcodEX, a tool for extracting Unique Molecular Identifiers (UMIs) from single or paired-end read sequences. bwaMem is used to generate alignments from UMI fastq files with the following deduplication and QC using imported bamQC workflow"
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
    umiCounts: {
        description: "JSON record of UMI counts after extraction",
        vidarr_label: "umiCounts"
    },
    extractionMetrics: {
        description: "JSON of metrics relating to extraction process",
        vidarr_label: "extractionMetrics"
    },
    preDedupBamMetrics: {
        description: "BamQC JSON report on bam file pre-deduplication",
        vidarr_label: "preDedupBamMetrics"
    },
    postDedupBamMetrics: {
        description: "BamQC JSON report on bam file post-deduplication",
        vidarr_label: "postDedupBamMetrics"
    },
    umiCountsPerPosition: {
        description: "tsv file tabulates the counts for unique combinations of UMI and position",
        vidarr_label: "umiCountsPerPosition"
    }
}
    }

    call getUMILengths {
      input:
          umiList = umiList    

    }

    Array[Int] umiLengths = getUMILengths.umiLengths

    call extractUMIs { 
        input:
            umiList = umiList,
            outputPrefix = outputPrefix,
            fastq1 = fastq1,
            fastq2 = fastq2,
            pattern1 = pattern1,
            pattern2 = pattern2
    }


    call bwaMem.bwaMem {
        input:
            fastqR1 = extractUMIs.fastqR1,
            fastqR2 = extractUMIs.fastqR2,
            outputFileNamePrefix = outputPrefix,
            runBwaMem_bwaRef = resources [ reference ].bwaRef,
            runBwaMem_modules = resources [ reference ].runBwaMem_modules
    }

    call bamQC.bamQC as preDedupBamQC {
        input:
            bamFile = bwaMem.bwaMemBam,
            outputFileNamePrefix = "~{outputPrefix}.preDedup"
    }

    scatter (umiLength in umiLengths) {
        call bamSplitDeduplication {
            input:
                bamFile = bwaMem.bwaMemBam,
                umiLength = umiLength,
                outputPrefix = outputPrefix
        }
    }


    call bamMerge {
        input:
            outputPrefix = outputPrefix,
            umiDedupBams = bamSplitDeduplication.umiDedupBams
    }

    call bamQC.bamQC as postDedupBamQC {
        input:
            bamFile = bamMerge.umiDedupBam,
            outputFileNamePrefix = "~{outputPrefix}.postDedup"
    }

    call statsMerge {
        input:
            umiCountsPerPositions = bamSplitDeduplication.dedupUmiCountsPerPosition,
            outputPrefix = outputPrefix
    }

    output {
        # barcodex metrics
        File umiCounts = extractUMIs.umiCounts
        File extractionMetrics = extractUMIs.extractionMetrics

        # metrix from umi-tools dedup
        File umiCountsPerPosition = statsMerge.umiCountsPerPosition

        # pre-collapse bamqc metrics
        File preDedupBamMetrics = preDedupBamQC.result

        # post-collapse bamqc metrics
        File postDedupBamMetrics = postDedupBamQC.result
    } 
}

task getUMILengths {
  input {
      File umiList
  }

  parameter_meta {
      umiList: "File with valid UMIs"
  }

  command <<<

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

  >>>

  output {

    Array[Int] umiLengths = read_lines(stdout())

  }

  meta {
            output_meta: {
                umiLengths: "An integer array of UMI lengths"
            }
        }

}

task extractUMIs {
        input {
            String umiList
            String outputPrefix
            File fastq1
            File fastq2
            String pattern1
            String pattern2
            String modules = "barcodex-rs/0.1.2 rust/1.45.1"
            Int memory = 24
            Int timeout = 12
        }

        parameter_meta {
            umiList: "Reference file with valid UMIs"
            outputPrefix: "Specifies the start of the output files"
            fastq1: "FASTQ file containing read 1"
            fastq2: "FASTQ file containing read 2"
            pattern1: "UMI RegEx pattern 1"
            pattern2: "UMI RegEx pattern 2"
            modules: "Required environment modules"
            memory: "Memory allocated for this job"
            timeout: "Time in hours before task timeout"
        }

        command <<<

            barcodex-rs --umilist ~{umiList} --prefix ~{outputPrefix} --separator "__" inline \
            --pattern1 '~{pattern1}' --r1-in ~{fastq1} \
            --pattern2 '~{pattern2}' --r2-in ~{fastq2} 

            cat ~{outputPrefix}_UMI_counts.json > umiCounts.txt

            tr [,] ',\n' < umiCounts.txt | sed 's/[{}]//' > tmp.txt
            echo "{$(sort -i tmp.txt)}" > new.txt
            tr '\n' ',' < new.txt | sed 's/,$//' > ~{outputPrefix}_UMI_counts.json
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


task bamSplitDeduplication {
    input {
        Int umiLength
        File bamFile
        String modules = "umi-tools/1.0.0 samtools/1.9"
        String outputPrefix
        Int memory = 24
        Int timeout = 6
    }

    parameter_meta {
        bamFile: "Bam file from bwaMem containing UMIs of varying lengths"
        umiLength: "Specifies the start of the output files"
        outputPrefix: "Specifies the start of the output files"
        modules: "Required environment modules"
        memory: "Memory allocated for this job"
        timeout: "Time in hours before task timeout"
    }

    command <<<
        samtools view -H ~{bamFile} > ~{outputPrefix}.~{umiLength}.sam
        samtools view ~{bamFile} | grep -P "^.*__\S{~{umiLength}}\t" >> ~{outputPrefix}.~{umiLength}.sam
        samtools view -Sb ~{outputPrefix}.~{umiLength}.sam > ~{outputPrefix}.~{umiLength}.bam

        samtools index ~{outputPrefix}.~{umiLength}.bam

        umi_tools dedup -I ~{outputPrefix}.~{umiLength}.bam \
        -S deduplicated.bam \
        --output-stats=deduplicated 
    >>>

    runtime {
        modules: "~{modules}"
        memory: "~{memory}G"
        timeout: "~{timeout}"
    }

    output {
        File umiDedupBams = "deduplicated.bam"
        File dedupUmiCountsPerPosition = "deduplicated_per_umi_per_position.tsv"
    }

    meta {
        output_meta: {
            umiDedupBams: "Bam files with deduplicated UMIs of varying lengths",
            dedupUmiCountsPerPosition: "tsv file tabulates the counts for unique combinations of UMI and position"
        }
    }
}

task bamMerge {
    input {
        Array[File] umiDedupBams
        String modules = "samtools/1.9"
        Int memory = 24
        Int timeout = 6
        String outputPrefix
    }

    parameter_meta {
        umiDedupBams: "Input bam files"
        outputPrefix: "Prefix for output file"
        memory: "Memory allocated for indexing job"
        modules: "Required environment modules"
        timeout: "Hours before task timeout"
    }
    
    String resultMergedBam = "~{outputPrefix}.dedup.bam"

    command <<<        
        set -euo pipefail
        samtools merge -c ~{outputPrefix}.dedup.bam ~{sep=" " umiDedupBams}
    >>>

    runtime {
        modules: "~{modules}"
        memory: "~{memory}G"
        timeout: "~{timeout}"
    }

    output {
        File umiDedupBam = "~{resultMergedBam}"
    }

    meta {
        output_meta: {
            umiDedupBam: "Deduplicated bam file"
        }
    }
}

task statsMerge {
    input{
        Array[File] umiCountsPerPositions
        Int memory = 16
        String outputPrefix
    }

    parameter_meta {
        umiCountsPerPositions: "An array of TSV files with umiCountsPerPosition"
        memory: "Memory allocated for this job"
    }

    command <<<
        set -euo pipefail

        umiCountsPerPositions=(~{sep=" " umiCountsPerPositions})
        length=${#umiCountsPerPositions[@]}
        touch umi.tsv
        i=0
        while [ $i -lt $length ]
        do
            cat umi.tsv > tmp.tsv
            awk  '{s2[$1]+=$2; s3[$1]+=$3} \
            END {for (k in s2) print k"\t"s2[k]"\t"s3[k]}' \
            tmp.tsv <(tail -n +2 ${umiCountsPerPositions[i]}) > umi.tsv
            i=$(( $i+1 ))
        done
        cat <(head -n 1 ${umiCountsPerPositions[0]}) umi.tsv > ~{outputPrefix}.umiCountsPerPosition.tsv
    >>>

    runtime {
    memory: "~{memory}G"
    }

    output {
        File umiCountsPerPosition = "~{outputPrefix}.umiCountsPerPosition.tsv"
    }
}
