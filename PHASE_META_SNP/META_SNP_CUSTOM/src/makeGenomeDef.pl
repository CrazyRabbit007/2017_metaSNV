#!/usr/bin/env perl
use strict;
use warnings;
use FAlite;

my $f = new FAlite(\*STDIN);
while (my $e=$f->nextEntry) {
    my $def = $e->def; $def =~ s/\s.*//;
    print substr($def, 1)."\t"."1"."\t".length($e->seq)."\n";
}
