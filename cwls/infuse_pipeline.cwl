#!/usr/bin/env cwl-runner

class: Workflow

label: "workflow to detect fusion events"

cwlVersion: v1.0

requirements:
  - class: InlineJavascriptRequirement

inputs:
  in_bam:
    type: File
    doc: "input BAM file."

  sample_name:
    type: string
    doc: "Sample name, which will used to prefix output file names."

  tophat_fusion_reference:
    type: File
    doc: "Tophat fusion reference bundle tar file."

  star_reference:
    type: File
    doc: "STAR reference bundle tar file."

  defuse_reference:
    type: File
    doc: "Defuse reference bundle tar file."

  gtf:
    type: File
    doc: "GTF file which is used to annotate break points."

  vagrent_cache:
    type: File
    doc: "VAGrENT cache file that should be the same reference and gene build as the GTF file."
    secondaryFiles:
      - ".tbi"
      - $(self.basename.replace("cache.gz", "fa"))
      - $(self.basename.replace("cache.gz", "fa.fai"))

  tophat_threads:
    type: int?
    default: 1
    doc: "Number of threads to use for tophat fusion process."

  star_threads:
    type: int?
    default: 1
    doc: "Number of threads to use for tophat fusion process."

  defuse_threads:
    type: int?
    default: 1
    doc: "Number of threads to use for tophat fusion process."

outputs:
  tophat_fusions:
    type: File
    outputSource: tophat_fusion/output

  star_fusions:
    type: File
    outputSource: star_fusion/output

  defuse_fusions:
    type: File
    outputSource: defuse/output

  combined_fusions:
    type: File
    outputSource: combine_fusions/output

steps:
  bam_to_fq:
    in:
      in_bam:
        source: in_bam
    out: [output_fqs]
    run: tools/bam_to_fq.cwl

  tophat_fusion:
    in:
      reads:
        source: bam_to_fq/output_fqs
      reference:
        source: tophat_fusion_reference
      sample_name:
        source: sample_name
      threads:
        source: tophat_threads
    out: [output]
    run: tools/run-cgprna_tophat-fusion.cwl

  star_fusion:
    in:
      reads:
        source: bam_to_fq/output_fqs
      reference:
        source: star_reference
      sample_name:
        source: sample_name
      threads:
        source: star_threads
    out: [output]
    run: tools/run-cgprna_star-fusion.cwl

  defuse:
    in:
      reads:
        source: bam_to_fq/output_fqs
      reference:
        source: defuse_reference
      sample_name:
        source: sample_name
      threads:
        source: defuse_threads
    out: [output]
    run: tools/run-cgprna_defuse.cwl

  combine_fusions:
    in:
      sample_name:
        source: sample_name
      gtf:
        source: gtf
      vagrent_cache:
        source: vagrent_cache
      tophat_fusion_output:
        source: tophat_fusion/output
      star_fusion_output:
        source: star_fusion/output
      defuse_output:
        source: defuse/output
    out: [output]
    run: tools/compare_overlapping_fusions.cwl

doc: |
  A workflow to run cpgRna Infuse pipeline from a single input BAM file (can be aligned or unaligned) of RNA-seq reads using cgpRna container. See the [cgpRna](https://github.com/cancerit/cgpRna) website for more information.

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
