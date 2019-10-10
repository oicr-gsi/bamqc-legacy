package GSI::util;

use strict;
use warnings;

#################### main pod documentation begin ###################

=begin html

<label for="show-menu" class="show-menu">Show Menu</label>
<input type="checkbox" id="show-menu" role="button">

<h1>Genome Sequence Informatics</h1>

=end html

=head1 util

=head2 NAME

util - Shared utility functions to generate the BamQC run report

=head2 SYNOPSIS

  use GSI::util;

=head2 AUTHOR

L<Genome Sequence Informatics|https://gsi.oicr.on.ca>,
L<Ontario Institute for Cancer Research|https://oicr.on.ca>.
On Github at L<https://github.com/oicr-gsi/bamqc>.

=head2 COPYRIGHT

Copyright (C) 2019 The Ontario Institute for Cancer Research

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 3 of the License, or (at
    your option) any later version.

This program is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
MA 02110-1301, USA.

=head1 Subroutines

=cut

#################### main pod documentation end ###################



BEGIN {
    use Exporter ();
    use vars qw($VERSION @ISA @EXPORT @EXPORT_OK %EXPORT_TAGS);

    $VERSION = 1.00;
    @ISA     = qw(Exporter);
    @EXPORT_OK =
      qw(is_single_read);
}

=for html <hr>

=head2 is_single_read()

Read the JSON hash output by BamQC, and determine if results are for a single read

B<Returns>

1 if single read; 0 otherwise (eg. paired-end)

=cut

sub is_single_read {
    # input: JSON output from BamQC
    my ($jsonHash) = @_;
    my $is_single = 0;
    if (defined($$jsonHash{"paired end"}) && !$$jsonHash{"paired end"}) {
	$is_single = 1; # 3.0+
    } elsif (defined($$jsonHash{"number of ends"}) && $$jsonHash{"number of ends"} eq "single end") {
	$is_single = 1; # 2.x
    }
    return $is_single;
}


1;
