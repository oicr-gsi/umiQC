## 1.2.0 - 2024-06-25
[GRD-797](https://jira.oicr.on.ca/browse/GRD-797) - add vidarr labels to outputs (changes to medata only)
## 1.0.2 - 2022-03-02
- [GP-2621](https://jira.oicr.on.ca/browse/GP-2621) - Overall speed and efficiency improvement
    - Updated UMI pattern regexps to get pulled in from UMI kit
    - Added new task `getUMILengths` to create array for scattered task
    - Combined `bamSplit` and `umiDeduplication` tasks into `bamSplitDeduplication` and applied scatter/gather structure
    - Created new tasks `bamMerge` and `umiMetrics`
    - `umiMetrics` merges all metrics TSVs into one file
    - Updated `calculate.sh` script

## 1.0.1 - 2021-06-21
- Workflow is built and tagged by Michelle Feng
