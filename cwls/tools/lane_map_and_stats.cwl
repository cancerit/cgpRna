class: Workflow
cwlVersion: v1.0

id: "lane-map-stats"

label: "map raw reads to reference and generate QC stats"

requirements:
  - class: InlineJavascriptRequirement

inputs:
  raw_reads:
    type:
      type: array
      items: File
  map_reference:
    type: File
  sample_name:
    type: string
  stats_reference:
    type: File
  map_threads:
    type: int?
  rg_id_tag:
    type: string?
  lb_tag:
    type: string?
  ds_tag:
    type: string?
  pl_tag:
    type: string?
  pu_tag:
    type: string?

outputs:
  dup_marked_bam:
    type: File
    outputSource: map/dup_marked_bam
  dup_marked_bam_dup_met:
    type: File
    outputSource: map/dup_marked_bam_dup_met
  transcriptome_bam:
    type: File
    outputSource: map/star_transcriptome_bam
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

steps:
  map:
    in:
      raw_reads:
        source: raw_reads
      reference:
        source: map_reference
      sample_name:
        source: sample_name
      output_file_prefix:
        source: sample_name
        valueFrom: |
          ${
            return self + '.lane.' + inputs.raw_reads[0].nameroot;
          }
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
    out: [star_transcriptome_bam, dup_marked_bam, dup_marked_bam_dup_met]
    run: run-cgprna_star-map.cwl
  stats:
    in:
      sample_bam:
        source: map/dup_marked_bam
      reference:
        source: stats_reference
      transcriptome_bam:
        source: map/star_transcriptome_bam
    out: [rna_bas, gene_cover_png, gene_body_coverage_rscript, gene_body_coverage_txt, gene_body_coverage_updated_rscript, read_dist]
    run: run-cgprna_stats.cwl
