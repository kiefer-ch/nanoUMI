# nanoUMI
Simple tool to remove PCR duplicates from Oxford Nanopore Technologie's PCR-cDNA sequencing runs.

# Algorithm overview

- Reads are assigned to a gene by [featureCounts](https://subread.sourceforge.net/)
- UMIs are extracted from the raw string by edlib
- UMIs are clustered within genes by vsearch
- If the TS tag is '-', the read is replaced by its reverse complement
- The longest read per cluster is reported

# General usage
```
nanoUMI -b {input.bam} -g {input.gtf} -o {output.fa} -c {threads} -l {log}
```

Input files must be basecalled by a recent version of [dorado](https://github.com/nanoporetech/dorado/) and mapped to the genome using e.g., [minimap2](https://github.com/lh3/minimap2).

Dorado writes the raw UMI string to the RX tag and the strand to the TS tag in the BAM output. Keep those by setting -T RX,TS when piping from unmapped bam files to minimap2 using samtools view.