version 1.0

import "imports/pull_bwaMem.wdl" as bwaMem
import "imports/pull_bamQC.wdl" as bamQC

workflow umiQC {
    input {
        String umiList
        String outputPrefix = "output"
        String fastq1
        String fastq2
	    String pattern1 
	    String pattern2
    }

    parameter_meta {
        umiList: "File with valid UMIs"
        outputPrefix: "Specifies the start of output files"
        fastq1: "Fastq file for read 1"
        fastq2: "Fastq file for read 2"
        pattern1: "UMI pattern 1"
        pattern2: "UMI pattern 2"
    }

    meta {
        author: "Michelle Feng and Murto Hilali"
        email: "mfeng@oicr.on.ca and mhilali@oicr.on.ca"
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
            mergedUMIMetrics: "TSV of files mapping read id to read group",
            postDedupBamMetrics: "BamQC report on bam file post-deduplication"
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
            outputFileNamePrefix = outputPrefix
    }

    call bamQC.bamQC as preDedupBamQC {
        input:
            bamFile = bwaMem.bwaMemBam,
            outputFileNamePrefix = "preDedup"
    }

    scatter (umiLength in umiLengths) {
        call bamSplitDeduplication {
            input:
                bamFile = bwaMem.bwaMemBam,
                umiLength = umiLength,
                outputPrefix = outputPrefix
        }
    }

    call mergeUMIs {
        input:
            umiMetrics = bamSplitDeduplication.umiMetrics,
    }

    call bamMerge {
        input:
            outputPrefix = outputPrefix,
            umiDedupBams = bamSplitDeduplication.umiDedupBams
    }

    call bamQC.bamQC as postDedupBamQC {
        input:
            bamFile = bamMerge.umiDedupBam,
            outputFileNamePrefix = "postDedup"
    }

    output {
        # barcodex metrics
        File umiCounts = extractUMIs.umiCounts
        File extractionMetrics = extractUMIs.extractionMetrics

        # pre-collapse bamqc metrics
        File preDedupBamMetrics = preDedupBamQC.result

        # umi-tools metrics
        #Array[File] umiMetrics = bamSplitDeduplication.umiMetrics
        
        File mergedUMIMetrics = mergeUMIs.mergedUMIMetrics

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
              L+=($(($i+$j+1)))
          done
      done

      umiLengths=($(tr ' ' '\n' <<< "${L[@]}" | awk '!u[$0]++' | tr ' ' '\n'))
      printf "%s\n" "${umiLengths[@]}"

  >>>

  output {

    Array[Int] umiLengths = read_lines(stdout())

  }

}

task extractUMIs {
        input {
            File umiList
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
            umiList: "File with valid UMIs"
            outputPrefix: "Specifies the start of the output files"
            fastqR1: "FASTQ file containing read 1"
            fastqR2: "FASTQ file containing read 2"
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

        umi_tools group -I ~{outputPrefix}.~{umiLength}.bam \
        --group-out=~{outputPrefix}.~{umiLength}.umi_groups.tsv \
        --output-bam > ~{outputPrefix}.~{umiLength}.dedup.bam \
        --log=group.log --paired | samtools view
    >>>

    runtime {
        modules: "~{modules}"
        memory: "~{memory}G"
        timeout: "~{timeout}"
    }

    output {
        File umiDedupBams = "output.~{umiLength}.dedup.bam"
        File umiMetrics = "output.~{umiLength}.umi_groups.tsv"
    }

    meta {
        output_meta: {
            umiDedupBams: "Bam files with deduplicated UMIs of varying lengths",
            umiMetrics: "File mapping read id to read group"
        }
    }
}

task mergeUMIs {
    input {
        Array[File] umiMetrics
    }

    parameter_meta {
        umiMetrics: "An array of TSV files with UMI metrics"

    }

    command <<<
        
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
        tr -s '[ , 	]' '\t' < mergedUMIMetrics.tsv > tmp.tsv && mv tmp.tsv mergedUMIMetrics.tsv 
    >>>

    output {
        File mergedUMIMetrics = "mergedUMIMetrics.tsv"
    }

    meta {
        output_meta: {
            mergedUMIMetrics: "A TSV of UMI metrics for all UMIs"
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
