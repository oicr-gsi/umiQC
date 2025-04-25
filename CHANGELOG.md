# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased] - 2025-04-23
### Added
- [GRD-795](https://jira.oicr.on.ca/browse/GRD-795) - Expanded built-in documentation (metadata changes only).

## [1.2.0] - 2024-06-25
### Added
- [GRD-797](https://jira.oicr.on.ca/browse/GRD-797) - add vidarr labels to outputs (changes to medata only)

## [1.0.2] - 2022-03-02
### Changed
- [GP-2621](https://jira.oicr.on.ca/browse/GP-2621) - Overall speed and efficiency improvement
- Updated UMI pattern regexps to get pulled in from UMI kit
- Combined `bamSplit` and `umiDeduplication` tasks into `bamSplitDeduplication` and applied scatter/gather structure
- `umiMetrics` merges all metrics TSVs into one file
- Updated `calculate.sh` script
### Added
- Added new task `getUMILengths` to create array for scattered task
- Created new tasks `bamMerge` and `umiMetrics`

## [1.0.1] - 2021-06-21
### Added
- Workflow is built and tagged by Michelle Feng
