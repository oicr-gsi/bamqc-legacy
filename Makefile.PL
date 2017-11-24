

use ExtUtils::MakeMaker;
# See lib/ExtUtils/MakeMaker.pm for details of how to influence
# the contents of the Makefile that is written.
WriteMakefile(
    NAME         => 'bamqc',
    VERSION_FROM => 'lib/GSI/bamqc.pm', # finds \$VERSION
    AUTHOR       => 'Morgan Taschuk (morgan.taschuk@oicr.on.ca)',
    ABSTRACT     => 'Generate quality control statistics from BAM files',
    PREREQ_PM    => {
                     'Test::Simple' => 0.44,
                    },
);
