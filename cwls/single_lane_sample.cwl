#!/usr/bin/env cwl-runner

class: Workflow

id: "single-lane-sample-workflow"

label: "A CGP workflow to generate mapping stats and gene counts from RNA-seq data using tools in cgpRna"

cwlVersion: v1.0

inputs:
  raw_reads:
    doc: "RAW read input, can be a bam file, or one of a pair of FastQ file. (optionally gzip compressed)."
    type:
      type: array
      items: File

  map_reference:
    type: File
    doc: "The core STAR reference and a GTF file bundled in a tar.gz."

  sample_name:
    type: string
    doc: "Sample name, which will used to prefix output file names and SM tag in the BAM file header."
    default: ''

  stats_reference:
    type: File
    doc: "The reference files bundled in a tar.gz."

  count_reference:
    type: File
    doc: "A reference GTF file."

  bigwig_reference:
    type: File
    doc: "FASTA file of a reference file, which the input BAM file was mapped to."
    secondaryFiles:
    - .fai

  bigwig_threads:
    type: int?
    doc: "Number of threads to use."

  map_threads:
    type: int?
    doc: "Number of threads to use."

  rg_id_tag:
    type: string?
    doc: "Readgroup ID tag value in the output BAM. Default: 1 or taken from the input raw BAM file."

  lb_tag:
    type: string?
    doc: "Sequencing library tag value in the output BAM header. Default: None or taken from the input raw BAM file."

  ds_tag:
    type: string?
    doc: "Description tag value in the output BAM header. Default: None or taken from the input raw BAM file."

  pl_tag:
    type: string?
    doc: "Platform tag value in the output BAM header. Default: None or taken from the input raw BAM file."

  pu_tag:
    type: string?
    doc: "Platform unit tag value in the output BAM header. Default: None or taken from the input raw BAM file."

outputs:
  star_transcriptome_bam:
    type: File
    outputSource: map/star_transcriptome_bam

  dup_marked_bam:
    type: File
    outputSource: map/dup_marked_bam

  dup_marked_bam_dup_met:
    type: File
    outputSource: map/dup_marked_bam_dup_met
    
  dup_marked_bam_md5:
    type: File
    outputSource: map/dup_marked_bam_md5

  rna_bas:
    type: File
    outputSource: stats/rna_bas

  gene_cover_png:
    type: File
    outputSource: stats/gene_cover_png

  gene_body_coverage_rscript:
    type: File
    outputSource: stats/gene_body_coverage_rscript

  gene_body_coverage_txt:
    type: File
    outputSource: stats/gene_body_coverage_txt

  gene_body_coverage_updated_rscript:
    type: File
    outputSource: stats/gene_body_coverage_updated_rscript

  read_dist:
    type: File
    outputSource: stats/read_dist

  out_bw:
    type: File
    outputSource: bigwig/out_bw

  out_count:
    type: File
    outputSource: count/out_count

steps:
  map:
    in:
      raw_reads:
        source: raw_reads
      reference:
        source: map_reference
      sample_name:
        source: sample_name
      threads:
        source: map_threads
      rg_id_tag:
        source: rg_id_tag
      lb_tag:
        source: lb_tag
      ds_tag:
        source: ds_tag
      pl_tag:
        source: pl_tag
      pu_tag:
        source: pu_tag
    out: [star_transcriptome_bam, dup_marked_bam, dup_marked_bam_dup_met, dup_marked_bam_md5]
    run: tools/run-cgprna_star-map.cwl

  stats:
    in:
      sample_bam:
        source: map/dup_marked_bam
      reference:
        source: stats_reference
      transcriptome_bam:
        source: map/star_transcriptome_bam
    out: [rna_bas, gene_cover_png, gene_body_coverage_rscript, gene_body_coverage_txt, gene_body_coverage_updated_rscript, read_dist]
    run: tools/run-cgprna_stats.cwl

  bigwig:
    in:
      sample_bam:
        source: map/dup_marked_bam
      reference:
        source: bigwig_reference
      threads:
        source: bigwig_threads
    out: [out_bw]
    run: tools/run-cgprna_bigwig.cwl

  count:
    in:
      sample_bam:
        source: map/dup_marked_bam
      reference:
        source: count_reference
    out: [out_count]
    run: tools/run-cgprna_htseq-count.cwl

doc: |
  A workflow to generate mapping stats and gene counts from RNA-seq data using cgpRna container. See the [cgpRna](https://github.com/cancerit/cgpRna) website for more information.

$schemas:
  - http://schema.org/docs/schema_org_rdfa.html

$namespaces:
  s: http://schema.org/

s:codeRepository: https://github.com/cancerit/cgpRna
s:license: https://spdx.org/licenses/AGPL-3.0

s:author:
  - class: s:Person
    s:email: mailto:yx2@sanger.ac.uk
    s:name: Yaobo Xu
