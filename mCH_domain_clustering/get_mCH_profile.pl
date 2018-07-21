#!/usr/bin/perl
use warnings;
use strict;

my $num_up_bin = 20;
my $num_domain_bin = 10;
my $num_down_bin = 20;

## Reading mCH profile of the whole genome
my $path = "/gale/netapp/home/yupeng/ENCODE/integrative/figures/mCH/";
my $mCH_file = "$path/data/mm10_bin5kb_with_CH_level.tsv.corrected";
my $bin_size = 5000;

open(DATA,"$mCH_file") or die "can't open $mCH_file!\n";
$_ = <DATA>;
chomp;

my @colnames = split(/\t/);

my %organ2ind;

foreach my $ind (3..$#colnames)
{
    my @vals = split(/\_/,$colnames[$ind]);
    # Organ
    my $organ = $vals[$#vals];
    push @{$organ2ind{$organ}},$ind - 3;
    push @{$organ2ind{"all"}},$ind - 3;
}

my %mCH;

while(<DATA>)
{
    chomp;
    next if length($_) == 0;
    my ($chr,$start,$end,@mCH_levels) = split(/\t/);
    ${$mCH{$chr}}{$start} = \@mCH_levels;
}
close(DATA);



my @tissues = (#E11.5
	       "E11_5_CF",
	       "E11_5_FB",
	       "E11_5_HB",
	       "E11_5_HT",
	       "E11_5_LM",
	       "E11_5_MB",
	       "E11_5_NT",
	       "E11_5_LV",
	       #E12.5
	       "E12_5_CF",
	       "E12_5_FB",
	       "E12_5_HB",
	       "E12_5_HT",
	       "E12_5_LM",
	       "E12_5_MB",
	       "E12_5_NT",
	       "E12_5_LV",
	       #E13.5
	       "E13_5_CF",
	       "E13_5_FB",
	       "E13_5_HB",
	       "E13_5_HT",
	       "E13_5_LM",
	       "E13_5_MB",
	       "E13_5_NT",
	       "E13_5_LV",
	       #E14.5
	       "E14_5_CF",
	       "E14_5_FB",
	       "E14_5_HB",
	       "E14_5_HT",
	       "E14_5_IT",
	       "E14_5_KD",
	       "E14_5_LG",
	       "E14_5_LM",
	       "E14_5_MB",
	       "E14_5_NT",
	       "E14_5_ST",
	       "E14_5_LV",
	       #E15.5
	       "E15_5_CF",
	       "E15_5_FB",
	       "E15_5_HB",
	       "E15_5_HT",
	       "E15_5_IT",
	       "E15_5_KD",
	       "E15_5_LG",
	       "E15_5_LM",
	       "E15_5_MB",
	       "E15_5_NT",
	       "E15_5_ST",
	       "E15_5_LV",
	       #E16.5
	       "E16_5_FB",
	       "E16_5_HB",
	       "E16_5_HT",
	       "E16_5_IT",
	       "E16_5_KD",
	       "E16_5_LG",
	       "E16_5_MB",
	       "E16_5_ST",
	       "E16_5_LV",
	       #P0
	       "P0_FB",
	       "P0_HB",
	       "P0_HT",
	       "P0_IT",
	       "P0_KD",
	       "P0_LG",
	       "P0_MB",
	       "P0_ST",
	       "P0_LV");

@tissues = ("merged");

foreach my $tissue (@tissues)
{
    my @vals = split(/\_/,$tissue);
    my $organ = $vals[$#vals];
    $organ = "all" if $tissue eq "merged";

    my $input_file = "$path/final/mCHdomain_$tissue\.bed";

    open(DATA,"$input_file") or die "can't open $input_file!\n";
    open(OUT,">mCHdomain_$tissue\.bed.profile") or die "can't open mCHdomain_$tissue\.bed.profile!\n";
        
    ## Print header
    print OUT join("\t",("chr","start","end"));
    foreach my $ind (@{$organ2ind{$organ}})
    {
	my $tis = $colnames[$ind + 3];
	print OUT "\t",join("\t",("upstream_$tis","domain_$tis","downstream_$tis"));
    }
    print OUT "\n";
        
    while(<DATA>)
    {
	chomp;
	next if length($_) == 0;
	my ($chr,$start,$end) = split(/\t/);
	$chr =~ s/chr//;

	print OUT join("\t",("chr".$chr,$start,$end));
	
	foreach my $ind (@{$organ2ind{$organ}})
	{
	    ## Upstream
	    my @up;
	    foreach my $i (reverse(1...$num_up_bin))
	    {
		my $pos = $start - $i* $bin_size;
		if(defined(${$mCH{$chr}}{$pos}))
		{
		    push @up,${${$mCH{$chr}}{$pos}}[$ind];
		}
		else
		{
		    push @up,"NA";
		}
	    }
	   
	    ## mCH domain
	    my @domain;
	    my $domain_bin_size = int(($end-$start)/$num_domain_bin);
	    foreach my $i (1...$num_domain_bin)
	    {
		my $bin_start = $start + ($i-1)* $domain_bin_size;
		my $bin_end = $start + $i* $domain_bin_size;
		$bin_end = $end if $bin_end > $end;

		my $pos_start = int($bin_start/$bin_size) * $bin_size;
		$pos_start += $bin_size if $pos_start != $bin_start;
		my $pos_end = int($bin_end/$bin_size) * $bin_size;
		
		my $pos = $pos_start;
		my $score = "NA";
		my $effective_len = 0;
		while($pos <= $pos_end)
		{
		    if(defined(${$mCH{$chr}}{$pos}) and 
		       ${${$mCH{$chr}}{$pos}}[$ind] ne "NA")
		    {
			$score = 0 if $score eq "NA";
			$score += ${${$mCH{$chr}}{$pos}}[$ind] * $bin_size;
			$effective_len += $bin_size;
		    }
		    $pos += $bin_size;
		}
		
		## Left boundary
		$pos = int($bin_start/$bin_size)*$bin_size;
		if($pos_start - $bin_start > 0 and
		   ${${$mCH{$chr}}{$pos}}[$ind] ne "NA")
		{
			$score = 0 if $score eq "NA";
			$score += ${${$mCH{$chr}}{$pos}}[$ind] * ($pos_start - $bin_start);
			$effective_len += ($pos_start - $bin_start);
		}

		## Right boundary
		$pos = int($bin_end/$bin_size)*$bin_size + $bin_size;
		if($bin_end - $pos_end > 0 and
		   ${${$mCH{$chr}}{$pos}}[$ind] ne "NA")
		{
		    $score = 0 if $score eq "NA";
		    $score += ${${$mCH{$chr}}{$pos}}[$ind] * ($bin_end - $pos_end);
		    $effective_len += ($bin_end - $pos_end);
		}
		
		## Output
		$score /= $effective_len if $score ne "NA" and $effective_len != 0;
		push @domain,$score;
	    }
	    
	    ## Downstream
	    my @down;
	    foreach my $i (1...$num_down_bin)
	    {
		my $pos = $end - ($i-1) * $bin_size;
		if(defined(${$mCH{$chr}}{$pos}))
		{
		    push @down,${${$mCH{$chr}}{$pos}}[$ind];
		}
		else
		{
		    push @down,"NA";
		}
	    }
	    print OUT "\t",join("\,",@up),"\t",join("\,",@domain),"\t",join("\,",@down);
	}
	print OUT "\n";
    }
    close(DATA);
    close(OUT);
	
}
