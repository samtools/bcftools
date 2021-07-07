[![Build Status](https://api.cirrus-ci.com/github/samtools/bcftools.svg?branch=develop)](https://api.cirrus-ci.com/github/samtools/bcftools)
[![Build status](https://ci.appveyor.com/api/projects/status/yanx5wnsdsqm4pay?svg=true)](https://ci.appveyor.com/project/samtools/bcftools)
[![Github All Releases](https://img.shields.io/github/downloads/samtools/bcftools/total.svg)](https://github.com/samtools/bcftools/releases/latest)

This is the official development repository for BCFtools. It contains all the vcf* commands
which previously lived in the htslib repository (such as vcfcheck, vcfmerge, vcfisec, etc.)
and the samtools BCF calling from bcftools subdirectory of samtools. 

For a full documentation, see [bcftools GitHub page](http://samtools.github.io/bcftools/). 

Other useful links:
------------------

File format specifications live on [HTS-spec GitHub page](http://samtools.github.io/hts-specs/)
[htslib](https://github.com/samtools/htslib)
[samtools](https://github.com/samtools/samtools)
[tabix](https://github.com/samtools/tabix)

### Citing

Please cite this paper when using BCFtools for your publications. http://samtools.github.io/bcftools/howtos/publications.html

> Twelve years of SAMtools and BCFtools </br>
> Petr Danecek, James K Bonfield, Jennifer Liddle, John Marshall, Valeriu Ohan, Martin O Pollard, Andrew Whitwham, Thomas Keane, Shane A McCarthy, Robert M Davies, Heng Li </br>
> _GigaScience_, Volume 10, Issue 2, February 2021, giab008, https://doi.org/10.1093/gigascience/giab008

```
@article{10.1093/gigascience/giab008,
    author = {Danecek, Petr and Bonfield, James K and Liddle, Jennifer and Marshall, John and Ohan, Valeriu and Pollard, Martin O and Whitwham, Andrew and Keane, Thomas and McCarthy, Shane A and Davies, Robert M and Li, Heng},
    title = "{Twelve years of SAMtools and BCFtools}",
    journal = {GigaScience},
    volume = {10},
    number = {2},
    year = {2021},
    month = {02},
    abstract = "{SAMtools and BCFtools are widely used programs for processing and analysing high-throughput sequencing data. They include tools for file format conversion and manipulation, sorting, querying, statistics, variant calling, and effect analysis amongst other methods.The first version appeared online 12 years ago and has been maintained and further developed ever since, with many new features and improvements added over the years. The SAMtools and BCFtools packages represent a unique collection of tools that have been used in numerous other software projects and countless genomic pipelines.Both SAMtools and BCFtools are freely available on GitHub under the permissive MIT licence, free for both non-commercial and commercial use. Both packages have been installed \\&gt;1 million times via Bioconda. The source code and documentation are available from https://www.htslib.org.}",
    issn = {2047-217X},
    doi = {10.1093/gigascience/giab008},
    url = {https://doi.org/10.1093/gigascience/giab008},
    note = {giab008},
    eprint = {https://academic.oup.com/gigascience/article-pdf/10/2/giab008/36332246/giab008.pdf},
}
```
