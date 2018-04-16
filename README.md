# rgsam

[![travis-ci](https://travis-ci.org/djhshih/rgsam.svg?branch=master)](https://travis-ci.org/djhshih/rgsam)
[![codecov](https://codecov.io/gh/djhshih/rgsam/branch/master/graph/badge.svg)](https://codecov.io/gh/djhshih/rgsam)

Infer read-group information from read names in SAM or FASTQ file.

# Installation

```{bash}
make
make install
```

# Usage

```{bash}
usage: rgsam [command]

commands:
  collect    collect read-group information from SAM or FASTQ file
  split      split SAM or FASTQ file based on read-group
  tag        tag reads in SAM file with read-group field
  qnames     list supported read name formats
  version    print version
```

Read-group identifier (`ID`) and platform unit (`PU`) are inferred from read
names according to supported read name formats:

```{yaml}
illumina-1.0:
    format: @{flowcell}-{instrument}:{lane}:{tile}:{x}:{y}#{sample}/{pair}
    example: @HWUSI-EAS100R:6:73:941:1973#0/1
illumina-1.8:
    format: @{flowcell}:{run}:{flowcell}:{lane}:{tile}:{x}:{y}
    exaample: @EAS139:136:FC706VJ:2:2104:15343:197393
broad-1.0:
    format: @{flowcell,5}:{barcode}:{lane}:{tile}:{x}:{y}
    example: @H0164ALXX140820:2:1101:10003:23460
```

Platform (`PL`) is assumed to be `illumina`.

Sample (`SM`) and library identifier (`LB`) may be inferred from input file name.

Files with reads from more than one sample or library are *not* supported.

To split BAM or SAM files containing proper `@RG` header lines and reads tagged
with read-group field (e.g. `RG:Z:H1`), use instead:

```{bash}
samtools view -r <rg_id> <in.bam>
```

## Example

Suppose we have a BAM file with no read-group data, then we first infer
the set of read-groups by

```{bash}
samtools view sample.bam | rgsam collect -s sample -o rg.txt
```

Now we can tag the reads with read-group information (any existing read-group
tags will be replaced).

```{bash}
samtools view -h sample.bam | 
  rgsam tag -r rg.txt |
  samtools view -b - > sample.rg.bam
```

Note that we use the `-h` flag of `samtools view` to ensure that other header data
are preserved (any existing `@RG` will be replaced).

