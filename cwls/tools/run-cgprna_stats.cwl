#!/usr/bin/env cwl-runner

class: CommandLineTool

id: "run-cgprna_stats"

label: "cgpRna stats"

cwlVersion: v1.0

doc: |
  ![build_status](https://quay.io/repository/wtsicgp/cgprna/status)
  A Docker container for the cgpRna mapping flow. See the [cgpRna](https://github.com/cancerit/cgpRna) website for more information.

  Please read the relevant [changes](https://github.com/cancerit/cgpRna/blob/dev/CHANGES.md) when upgrading.

  Parameters for a CWL definition are generally described in a json file, but parameters can be provided on the command line.

  To see the parameters descriptions please run: cwltool --tool-help path_to.cwl

requirements:
  - class: DockerRequirement
    dockerPull: "quay.io/wtsicgp/cgprna:2.3.4"

hints:
  - class: ResourceRequirement
    coresMin: 1
    ramMin: 2000

inputs:
  sample_bam:
    type: File
    doc: "Input BAM file, in which reads are mapped to a reference genome (NOT transcriptome)"
    inputBinding:
      prefix: --input
      separate: true
    secondaryFiles:
    - .bai

  reference:
    type: File
    doc: "The reference files bundled in a tar.gz."
    inputBinding:
      prefix: --reference
      separate: true

  transcriptome_bam:
    type: File
    doc: "BAM file, in which reads are mapped to a reference transciptome (NOT genome)."
    inputBinding:
      prefix: --transcriptome-bam
      separate: true
    secondaryFiles:
    - .bai

outputs:
  rna_bas:
    type: File
    outputBinding:
      glob: $(inputs.sample_bam.nameroot).RNA.bas

  gene_cover_png:
    type: File
    outputBinding:
      glob: $(inputs.sample_bam.nameroot).geneBodyCoverage.curves.png

  gene_body_coverage_rscript:
    type: File
    outputBinding:
      glob: $(inputs.sample_bam.nameroot).geneBodyCoverage.r

  gene_body_coverage_txt:
    type: File
    outputBinding:
      glob: $(inputs.sample_bam.nameroot).geneBodyCoverage.txt

  gene_body_coverage_updated_rscript:
    type: File
    outputBinding:
      glob: $(inputs.sample_bam.nameroot).geneBodyCoverage_UPDATED.r

  read_dist:
    type: File
    outputBinding:
      glob: $(inputs.sample_bam.nameroot).read_dist.txt

baseCommand: ["run-cgprna", "stats"]

$schemas:
  - http://schema.org/docs/schema_org_rdfa.html

$namespaces:
  s: http://schema.org/

s:codeRepository: https://github.com/cancerit/cgpRna
s:license: https://spdx.org/licenses/AGPL-3.0-only

s:author:
  - class: s:Person
    s:email: mailto:cgphelp@sanger.ac.uk
    s:name: Yaobo Xu
