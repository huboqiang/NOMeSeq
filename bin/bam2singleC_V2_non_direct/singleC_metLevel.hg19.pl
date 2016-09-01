#!/usr/bin/perl 
#===============================================================================
#
#         FILE: singleC_metLevel.hg19.pl
#
#        USAGE: ./singleC_metLevel.hg19.pl  
#
#  DESCRIPTION: 
#
#      OPTIONS: ---
# REQUIREMENTS: ---
#         BUGS: ---
#        NOTES: ---
#       AUTHOR: Hu Boqiang (phD Candidate), huboqiang@gmail.com
# ORGANIZATION: BioDynamic Optical Imaging Center
#      VERSION: 1.0
#      CREATED: 10/29/13 13:39:25
#     REVISION: ---
#===============================================================================

use warnings;

my $usage = <<USAGE;

Usage:

      perl $0 /data/Analysis/huboqiang/database/hg19/hg19.fa file1.pileup > file.single5mC

USAGE

die $usage if @ARGV < 1;
open R,"<$ARGV[0]" or die "Cannot open file $ARGV[0]($!) \n";
my $chr_seq;
while (<R>){
   my $chr = $1 if /^>?(\S+)\s/o;
   $/ = '>';
   my $seq = <R>;
   $seq =~ s/\n|>//g;
   $seq = uc $seq;
   $chr_seq{$chr} = $seq;
   $/ = "\n";
}
close R;

open I,"<$ARGV[1]" or die "Cannot open file $ARGV[1]($!) \n";
while (<I>){
   my ($chr,$pos,$ref,$dep,$base,$qual,$oth)  = split;
   my ($str,$tot,$met,$umt,$ratio,$tri,$type) = ();
   $ref = uc $ref;
   next if  ($ref !~ /[CG]/i or $chr =~ /random/ or $chr eq "chrM" or $dep == 0);
   my $beg = $pos - 1;
   if ($ref =~ /C/i){
      $str = "+";
      ($met,$umt) = ($base =~ tr/\./\./, $base =~ tr/T/T/);
      $tot = $met + $umt;
      next if $tot == 0;
      $tri = substr($chr_seq{$chr},$beg,3);
      $tri = uc $tri;
   }
   else {
      $str = "-";
      ($met,$umt) = ($base =~ tr/,/,/,   $base =~ tr/a/a/);
      $tot = $met + $umt;
      next if $tot == 0 or $beg == 0;
      $tri = ( $beg == 1 ) ? ( substr($chr_seq{$chr},0,2) ) : ( substr($chr_seq{$chr},$beg-2,3) );
      $tri = uc $tri;
      $tri =~ tr/ATGC/TACG/;
      $tri = reverse $tri;
   }
   if (length($tri) == 2){
      my $base2 = substr($tri,1,1);
      if ($base2 =~ /G/i){
         $type = "CpG";
      }
   }
   else{
      if($tri=~/C[ATC]G/i){
         $type = "CHG";
      }
      elsif ($tri =~ /C[ATC][ATC]/){
         $type = "CHH";
      }
      elsif ($tri =~ /^CG/i){
         $type = "CpG";
      }
   }
   $ratio = $met / $tot;
   print "$chr\t$pos\t$ref\t$str\t$tot\t$met\t$umt\t$ratio\t$tri\t$type\n" if defined $type;
}
