#!/usr/bin/perl -w


=head1 gro_to_openmm_itp.pl
#"""
#Convert a set of SAPT-FF gromacs input itp files into openmm xml force field
#Usage : perl gro_to_openmm_itp.pl [paramfile]
#@author Chang Yun Son <cson@chem.wisc.edu>
#"""
=cut


#==================#
#| Global Imports |#
#==================#

use strict;
use warnings;
use POSIX;
use Math::Trig ':pi';
use File::Copy;
use List::MoreUtils 'pairwise';

my @param = @ARGV;
my $paramfile = shift @param;
my $qposfile = shift @param;
s/^\s+|\s+$//g foreach @param;                #remove white spaces


#Read q position file
my @q;
open (QPOS, "< " . $qposfile) or die "Cannot open input q position coord file " . $qposfile . "!\n";
my $idx = 0;
while(<QPOS>)
{
    my ($x,$y)=split;
    $q[$idx++] = $x;
}
close(QPOS);


#Read aff parameter file
my %hsparams;
open (PARAM, "< " . $paramfile) or die "Cannot open input q position coord file " . $paramfile . "!\n";
<PARAM>;     #Read the comment header line
while(<PARAM>)
{
    my @paramdata = split;
    s/^\s+|\s+$//g foreach @paramdata;                #remove white spaces
    my ($atname, $a0, $b0, $a1, $b1, $a2, $b2, $a3, $b3, $c)=@paramdata;
    my @arra = ($a0, $a1, $a2, $a3);
    my @arrb = ($b0, $b1, $b2, $b3);
    $hsparams{$atname} = [([@arra], [@arrb], $c)];
}
close(PARAM);

my $numpara = @param;
my $qadjust = 0;
if ($numpara > 0) {
    my $writename = shift @param;
    if(@param == 0) {
        if (exists $hsparams{$writename}) {
            my ($refa, $refb, $c) = @{$hsparams{$writename}};
            write_aff($writename, $refa, $refb, $c, 0);
        }
        else { die "There is no AFF param data for atom name $writename\n Please check the aff param file $paramfile\n"; }
    }
    else {
        my $numparams = @param/2;
        if (@param%2 ==0) {
            my (@arra, @arrb);
            my $cval = 0;
            for (my $i=0; $i<$numparams; $i++) {
                my $atname = $param[$i*2];
                if (exists $hsparams{$atname}) {
                    my ($refa, $refb, $c) = @{$hsparams{$atname}};
                    push(@arra, map { $_ * $param[$i*2+1] } @{$refa});
                    push(@arrb,@{$refb});
                    $cval += $c * $param[$i*2+1];
                }
                elsif ($atname eq 'q') {
                    $qadjust = $param[$i*2+1];
                }
                else { die "There is no AFF param data for atom name $atname\n Please check the aff param file $paramfile\n"; }
            }
            write_aff($writename, \@arra, \@arrb, $cval, $qadjust);
            
        }
        else { die "Number of parameters are not properly set\n Usage: write_form_factors.pl [paramfile] [qposfile] aff_write_name {atomi_name atomi_weight}\n"; }
    }
}
else {
    foreach my $atname (keys %hsparams) {
        my ($refa, $refb, $c) = @{$hsparams{$atname}};
        write_aff($atname, $refa, $refb, $c);
    }
}
exit;


sub write_aff {
    my ($atname, $refa, $refb, $c, $q) = @_;
    my @arra = @{$refa};
    my @arrb = @{$refb};
    my $filename = "AFF_$atname";
    my $bakfile = my $prefile = "$filename.out";
    my $idx = 0;
    my $qsum = 0;
    map { $qsum += $_ } @arra;
    $qsum += $c;
    
    if(-e $bakfile) {
     while(-e $bakfile) {
      $bakfile = "#" . $filename . "_$idx";
      $idx++;
     }
     copy($prefile, $bakfile);
    }
    open(AFF, "> " . $prefile);
    my ($qval, $qov4pisq, $fq);
    for (my $i=0; $i<@q; $i++) {
        $qval = $q[$i]/10;          # coefficients are defined in A^-1 but the qpos are defined in nm^-1
        $qov4pisq = ($qval / (4*pi))**2 ;
        $fq = $c;
        $fq += $_ for pairwise { $a * exp(-$b*$qov4pisq)} @arra, @arrb;
        $fq *= ($qsum - $q)/$qsum;
        print AFF "$q[$i], $fq\n";
    }
    close(AFF);
}

