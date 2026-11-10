#!/usr/bin/env perl
#
# See the NOTICE file distributed with this work for additional information
# regarding copyright ownership.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#
# Perl-side oracle for `vep bench oracle` (ensembl.vep.bench.oracle): drives
# Bio::EnsEMBL::Variation::Utils::VariationEffect directly against hand-built mock objects, so it
# has no dependency on this repository's Python model. build_cases below must stay aligned, by name
# and order, with CASE_NAMES in ensembl.vep.bench.oracle; a test pins the two together.

use strict;
use warnings;

use Getopt::Long qw(GetOptions);
use JSON::PP qw(encode_json);
use Bio::EnsEMBL::Variation::Utils::VariationEffect ();

package MockTranscript;

sub new {
  my ($class, %args) = @_;
  return bless \%args, $class;
}

sub start { $_[0]->{start} }
sub end { $_[0]->{end} }
sub strand { $_[0]->{strand} }
sub translation { $_[0]->{translation} }
sub coding_region_start { $_[0]->{coding_region_start} }
sub coding_region_end { $_[0]->{coding_region_end} }
sub biotype { $_[0]->{biotype} || 'protein_coding' }

package MockVariationFeature;

sub new {
  my ($class, %args) = @_;
  return bless \%args, $class;
}

sub start { $_[0]->{start} }
sub end { $_[0]->{end} }
sub seq_region_start { $_[0]->{start} }
sub seq_region_end { $_[0]->{end} }

package MockVariationFeatureOverlap;

sub new {
  my ($class, %args) = @_;
  return bless \%args, $class;
}

sub feature { $_[0]->{feature} }
sub base_variation_feature { $_[0]->{base_variation_feature} }

package MockVariationFeatureOverlapAllele;

sub new {
  my ($class, %args) = @_;
  return bless \%args, $class;
}

sub base_variation_feature_overlap { $_[0]->{base_variation_feature_overlap} }
sub base_variation_feature { $_[0]->{base_variation_feature} }

sub _intron_effects {
  my ($self) = @_;
  return {
    intronic => $self->{intron_number} ? 1 : 0,
    within_frameshift_intron => 0,
    splice_region => 0,
  };
}

package main;

sub build_cases {
  return (
    { name => 'coding_exact', start => 10, end => 10, translation => 1, intron_number => 0 },
    { name => 'coding_near_boundary', start => 11, end => 11, translation => 1, intron_number => 0 },
    { name => 'intron_primary', start => 40, end => 40, translation => 0, intron_number => 1 },
    { name => 'intron_secondary', start => 41, end => 41, translation => 0, intron_number => 2 },
    { name => 'upstream_near', start => 0, end => 0, translation => 0, intron_number => 0 },
    { name => 'upstream_far', start => -50, end => -50, translation => 0, intron_number => 0 },
    { name => 'downstream_near', start => 91, end => 91, translation => 0, intron_number => 0 },
    { name => 'downstream_far', start => 140, end => 140, translation => 0, intron_number => 0 },
    { name => 'no_hit', start => 6000, end => 6000, translation => 0, intron_number => 0 },
  );
}

sub build_overlap {
  my ($case) = @_;
  my $transcript = MockTranscript->new(
    start => 1,
    end => 90,
    strand => 1,
    translation => $case->{translation} ? 1 : undef,
    coding_region_start => 1,
    coding_region_end => 90,
    biotype => $case->{translation} ? 'protein_coding' : 'lncRNA',
  );
  my $bvf = MockVariationFeature->new(
    start => $case->{start},
    end => $case->{end},
  );
  my $bvfo = MockVariationFeatureOverlap->new(
    feature => $transcript,
    base_variation_feature => $bvf,
  );
  return MockVariationFeatureOverlapAllele->new(
    base_variation_feature_overlap => $bvfo,
    base_variation_feature => $bvf,
    intron_number => $case->{intron_number},
  );
}

sub consequence_terms {
  my ($overlap) = @_;
  my @terms;

  my $bvfo = $overlap->base_variation_feature_overlap;
  my $feat = $bvfo->feature;
  my $bvf = $overlap->base_variation_feature;

  push @terms, 'coding_sequence_variant'
    if $feat->translation
    && $bvf->start >= $feat->start
    && $bvf->end <= $feat->end;
  push @terms, 'intron_variant' if Bio::EnsEMBL::Variation::Utils::VariationEffect::within_intron($overlap);
  push @terms, 'upstream_gene_variant' if Bio::EnsEMBL::Variation::Utils::VariationEffect::upstream($overlap);
  push @terms, 'downstream_gene_variant' if Bio::EnsEMBL::Variation::Utils::VariationEffect::downstream($overlap);

  return \@terms;
}

sub benchmark {
  my (%args) = @_;
  my $cases = $args{cases};
  my $number = $args{number};
  my $repeats = $args{repeats};

  my $count = 0;
  for (1 .. $number) {
    for (1 .. $repeats) {
      for my $case (@$cases) {
        my $terms = consequence_terms(build_overlap($case));
        $count += @$terms;
      }
    }
  }

  return $count;
}

sub main {
  my $mode = 'bench';
  my $number = 1;
  my $repeats = 1;

  GetOptions(
    'mode=s' => \$mode,
    'number=i' => \$number,
    'repeats=i' => \$repeats,
  ) or die "invalid options\n";

  my @cases = build_cases();

  if ($mode eq 'json') {
    my @terms = map { consequence_terms(build_overlap($_)) } @cases;
    print encode_json(\@terms);
    return 0;
  }

  my $count = benchmark(cases => \@cases, number => $number, repeats => $repeats);
  print $count, "\n";
  return 0;
}

exit main();
