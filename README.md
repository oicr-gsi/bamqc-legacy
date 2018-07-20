[![Build Status](https://travis-ci.org/oicr-gsi/bamqc.svg)](https://travis-ci.org/oicr-gsi/bamqc)

# bamqc
Perl scripts for generating quality control stats from BAM files

To install dependencies (requires cpanminus): 

    cpanm --installdeps .

To use the report generating scripts, you will also need to install the [gsi-website](https://github.com/oicr-gsi/gsi-website)  perl module, freely available from Github.

    git clone git@github.com:oicr-gsi/gsi-website.git
    cd gsi-website
    perl Makefile.PL
    make
    make install

## Installation

To install this module:

    perl Makefile.PL && make && make install && make test
    
To run tests alone:

    make test
    


## Documentation

For man pages, see [man/GSI__bamqc](man/GSI__bamqc) and [man/GSI__RunReport](man/GSI__RunReport).

To generate this documentation:

    pod2text lib/GSI/bamqc.pm > man/GSI__bamqc
    pod2text lib/GSI/RunReport.pm > man/GSI__RunReport

For the html pages:

    pod2html --css css/skeleton.css lib/GSI/bamqc.pm --noindex > docs/bamqc.html
    pod2html --css css/skeleton.css lib/GSI/RunReport.pm --noindex > docs/RunReport.html


## Bedtools histogram files

The -H option in bamqc.pl expects a 'bedtools histogram file' for coverage. Here's how to generate this file:

    samtools rmdup MY_BAM.bam - | bedtools coverage -hist -abam stdin -b targets.bed | grep all | sed 's/all/collapsed/g' >> MY_BAM.hist
    bedtools coverage -hist -abam MY_BAM.bam -b targets.bed | grep all | sed 's/all/noncollapsed/g' >> MY_BAM.hist

