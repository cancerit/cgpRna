#!/usr/bin/env cwl-runner

class: CommandLineTool

id: "run-cgprna_comp-fusions"

label: "compare fusions"

cwlVersion: v1.0

doc: |
  ![build_status](https://quay.io/repository/wtsicgp/cgprna/status)
  A Docker container for the cgpRna mapping flow. See the [cgpRna](https://github.com/cancerit/cgpRna) website for more information.

  Please read the relevant [changes](https://github.com/cancerit/cgpRna/blob/dev/CHANGES.md) when upgrading.

  Parameters for a CWL definition are generally described in a json file, but parameters can be provided on the command line.

  To see the parameters descriptions please run: cwltool --tool-help path_to.cwl

requirements:
  - class: DockerRequirement
    dockerPull: "quay.io/wtsicgp/cgprna:2.5.0"
  - class: InlineJavascriptRequirement

hints:
  - class: ResourceRequirement
    ramMin: 8000

inputs:
  sample_name:
    type: string
    doc: "sample name which is used to name the output file."
    inputBinding:
      prefix: -s
      separate: true
      shellQuote: true
      position: 1

  gtf:
    type: File
    doc: "GTF file which is used to annotate break points."
    inputBinding:
      prefix: -g
      separate: true
      shellQuote: true
      position: 2

  vagrent_cache:
    type: File
    doc: "VAGrENT cache file that should be the same reference and gene build as the GTF file."
    inputBinding:
      prefix: -c
      separate: true
      shellQuote: true
      position: 3
    secondaryFiles:
      - ".tbi"
      - $(self.basename.replace("cache.gz", "fa"))
      - $(self.basename.replace("cache.gz", "fa.fai"))

  tophat_fusion_output:
    type: File
    doc: "Output file of Tophat fusion, usually with suffix: '.tophat-fusion.normals.filtered.strand.txt'."
    inputBinding:
      shellQuote: true
      position: 4

  star_fusion_output:
    type: File
    doc: "Output file of Star fusion, usually with suffix: '.star-fusion.normals.filtered.txt'."
    inputBinding:
      shellQuote: true
      position: 5

  defuse_output:
    type: File
    doc: "Output file of tophat fusion, usually with suffix: '.defuse-fusion.normals.ext.filtered.txt'."
    inputBinding:
      shellQuote: true
      position: 6

outputs:
  output:
    type: File
    outputBinding:
      glob: $(inputs.sample_name).infuse.detected.fusions.grass.txt

baseCommand: ["compare_overlapping_fusions.pl", "-o", "."]

$schemas:
  - https://schema.org/version/latest/schema.rdf

$namespaces:
  s: http://schema.org/

s:codeRepository: https://github.com/cancerit/cgpRna
s:license: https://spdx.org/licenses/AGPL-3.0-only

s:author:
  - class: s:Person
    s:email: mailto:cgphelp@sanger.ac.uk
    s:name: Yaobo Xu
