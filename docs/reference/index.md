# Package index

## Pipeline stages

Core SORTER2 stages run in sequence on target-capture reads.

- [`sorter2_stage1a()`](sorter2_stage1a.md) : Run the SORTER2 Stage1A
  step
- [`sorter2_stage1b()`](sorter2_stage1b.md) : Run the SORTER2 Stage1B
  step
- [`sorter2_stage2()`](sorter2_stage2.md) : Run the SORTER2 Stage2 step
- [`sorter2_stage3()`](sorter2_stage3.md) : Run the SORTER2 Stage3 step

## Post-processing

Downstream analysis tools run after phasing is complete.

- [`sorter2_processor()`](sorter2_processor.md) : Run the SORTER2
  Processor step
- [`sorter2_progenitor_processor()`](sorter2_progenitor_processor.md) :
  Process progenitor sequences from Stage 3 output
- [`sorter2_haplominer()`](sorter2_haplominer.md) : Run the SORTER2
  Hapl-O-Miner step

## Utilities

Helper functions for input preparation and workflow integration.

- [`sorter2_filter_assemblies()`](sorter2_filter_assemblies.md) : Copy
  ingroup assembly folders from Stage 1A output, excluding outgroups
- [`sorter2_infer_hybrids()`](sorter2_infer_hybrids.md) : Identify
  candidate hybrid or polyploid samples from Stage 1B output
- [`sorter2_format_reads()`](sorter2_format_reads.md) : Run the SORTER2
  FormatReads step
- [`sorter2_format_ref()`](sorter2_format_ref.md) : Reformat a reference
  FASTA for SORTER2
- [`sorter2_run()`](sorter2_run.md) : Run a SORTER2 Python script
