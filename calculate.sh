find . \( -type f -size +0 -iname "output.umi_groups.tsv" \) -printf "output umi_groups file exists\n";
find . \( -type f -size +0 -iname "output_extraction_metrics.json" \) -printf "extraction metrics file exists\n";
find . \( -type f -size +0 -iname "preDedup.bamQC_results.json" \) -printf "preDedup bamQC file exists\n";
find . \( -type f -size +0 -iname "postDedup.bamQC_results.json" \) -printf "postDedup bamQC file exists\n";
find . \( -type f -size +0 -iname "output_UMI_counts.json" \) -printf "UMI_counts file exists\n";
