#!/usr/bin/perl
use warnings;
use strict;

my $max_iter = 5;

my @control_samples = (#E10.5
		       "E10_5_FB_1","E10_5_FB_2",
		       "E10_5_MB_1","E10_5_MB_2",
		       "E10_5_HB_1","E10_5_HB_2",
		       "E10_5_CF_1","E10_5_CF_2",
		       "E10_5_LM_1","E10_5_LM_2",
		       #"E10_5_HT_1","E10_5_HT_2",
		       #E11.5
		       "E11_5_CF_1","E11_5_CF_2",
		       "E11_5_FB_1","E11_5_FB_2",
		       "E11_5_HB_1","E11_5_HB_2",
		       #"E11_5_HT_1","E11_5_HT_2",
		       "E11_5_LM_1","E11_5_LM_2",
		       "E11_5_MB_1","E11_5_MB_2",
		       "E11_5_NT_1","E11_5_NT_2",
		       "E11_5_LV_1","E11_5_LV_2");

`mkdir -p tmp/ final/`;
`rm -f tmp/excl_mCHdomain.bed`;

print "Start mCH domain calling\n";
my $t = gmtime();
print $t,"\n\n";
    
my $iter = 1;
while(1)
{
    print "Iteration $iter: ";
    
    `Rscript call_changepoint.R Control`;
    
    my @excl_domain_files;
    foreach my $sample (@control_samples)
    {
	
	my $input_file = "tmp/mCH_segment_$sample\.bed";
	my $output_file = "tmp/mCHdomain_$sample\.bed";
	push @excl_domain_files, $output_file;
	open(DATA,"$input_file") or die "can't open $input_file!\n";
	my %cand;
	while(<DATA>)
	{
	    chomp;
	    next if length($_) == 0;
	    my ($chr) = split(/\t/,$_);
	    push @{$cand{$chr}},$_;
	}
	close(DATA);
	
	open(OUT,">$output_file") or die "can't create $output_file!\n";
	## get mCH domain
	foreach my $chr (keys(%cand))
	{
	    my @cand_chr = @{$cand{$chr}};
	    if(@cand_chr == 1)
	    {
		next;
	    }
	    elsif(@cand_chr == 2)
	    {
		next;
	    }
	    elsif(@cand_chr > 2)
	    {
		my @mCH_array;
		foreach my $i (0...$#cand_chr)
		{
		    my ($chr,$start,$end,$mCH) = split(/\t/,$cand_chr[$i]);
		    push @mCH_array,$mCH;
		}
		
		foreach my $i (1...($#cand_chr-1))
		{		
		    #if($mCH_array[$i] > $mCH_array[$i-1] and
		    #   $mCH_array[$i] > $mCH_array[$i+1])
		    if( ($mCH_array[$i] + 0.001) / ($mCH_array[$i-1] + 0.001) >= 1.5 and
			($mCH_array[$i] + 0.001) / ($mCH_array[$i+1] + 0.001) >= 1.5)
		    {
			print OUT "chr".$cand_chr[$i],"\t";
			print OUT join(",",($mCH_array[$i-1],$mCH_array[$i],$mCH_array[$i+1])),"\n";
		    }
		}	    
	    }
	}
	close(OUT);
    }


    my $cmd = "cat " . join(" ",@excl_domain_files) . "|sort -k 1,1 -k 2,2g";
    $cmd .= "|mergeBed -i -|sed s/chr// >> tmp/excl_mCHdomain.bed";
    `$cmd`;

    ## Check if all fake mCH domains are excluded in this iterative process

    $cmd = "cat " . join(" ",@excl_domain_files) . "|wc -l";
    my $count = `$cmd`;
    chomp($count);
    print $count," mCH domains in control samples.\n";

    ## Print out the time
    $t = gmtime();
    print $t,"\n\n";

    ## Break out of the circle?
    last if($count == 0 or $iter > $max_iter );
    $iter++;
}
print "Iteration finished!\n";

$t = gmtime();
print $t,"\n\n";

`cp tmp/excl_mCHdomain.bed final/`;

## Call mCH domains
my @samples = (#E10.5
	       "E10_5_FB_1","E10_5_FB_2",
	       "E10_5_MB_1","E10_5_MB_2",
	       "E10_5_HB_1","E10_5_HB_2",
	       "E10_5_CF_1","E10_5_CF_2",
	       "E10_5_LM_1","E10_5_LM_2",
	       "E10_5_HT_1","E10_5_HT_2",
	       #E11.5
	       "E11_5_CF_1","E11_5_CF_2",
	       "E11_5_FB_1","E11_5_FB_2",
	       "E11_5_HB_1","E11_5_HB_2",
	       "E11_5_HT_1","E11_5_HT_2",
	       "E11_5_LM_1","E11_5_LM_2",
	       "E11_5_MB_1","E11_5_MB_2",
	       "E11_5_NT_1","E11_5_NT_2",
	       "E11_5_LV_1","E11_5_LV_2",
	       #E12.5
	       "E12_5_CF_1","E12_5_CF_2",
	       "E12_5_FB_1","E12_5_FB_2",
	       "E12_5_HB_1","E12_5_HB_2",
	       "E12_5_HT_1","E12_5_HT_2",
	       "E12_5_LM_1","E12_5_LM_2",
	       "E12_5_LV_1","E12_5_LV_2",
	       "E12_5_MB_1","E12_5_MB_2",
	       "E12_5_NT_1","E12_5_NT_2",
	       #E13.5
	       "E13_5_CF_1","E13_5_CF_2",
	       "E13_5_FB_1","E13_5_FB_2",
	       "E13_5_HB_1","E13_5_HB_2",
	       "E13_5_HT_1","E13_5_HT_2",
	       "E13_5_LM_1","E13_5_LM_2",
	       "E13_5_MB_1","E13_5_MB_2",
	       "E13_5_NT_1","E13_5_NT_2",
	       "E13_5_LV_1","E13_5_LV_2",
	       #E14.5
	       "E14_5_CF_1","E14_5_CF_2",
	       "E14_5_FB_1","E14_5_FB_2",
	       "E14_5_HB_1","E14_5_HB_2",
	       "E14_5_HT_1","E14_5_HT_4",
	       "E14_5_IT_1","E14_5_IT_2",
	       "E14_5_KD_1","E14_5_KD_4",
	       "E14_5_LG_1","E14_5_LG_2",
	       "E14_5_LM_1","E14_5_LM_2",
	       "E14_5_MB_1","E14_5_MB_2",
	       "E14_5_NT_1","E14_5_NT_2",
	       "E14_5_ST_1","E14_5_ST_4",
	       "E14_5_LV_2","E14_5_LV_3",
	       #E15.5
	       "E15_5_CF_1","E15_5_CF_2",
	       "E15_5_HT_1","E15_5_HT_2",
	       "E15_5_LG_1","E15_5_LG_2",
	       "E15_5_KD_1","E15_5_KD_2",
	       "E15_5_IT_1","E15_5_IT_2",
	       "E15_5_ST_1","E15_5_ST_2",
	       "E15_5_LM_1","E15_5_LM_2",
	       "E15_5_FB_1","E15_5_FB_2",
	       "E15_5_MB_1","E15_5_MB_2",
	       "E15_5_HB_1","E15_5_HB_2",
	       "E15_5_NT_1","E15_5_NT_2",
	       "E15_5_LV_1","E15_5_LV_2",
	       #E16.5
	       "E16_5_HB_1","E16_5_HB_3",#"E16_5_HB_2",
	       "E16_5_MB_1","E16_5_MB_3",#"E16_5_MB_2",
	       "E16_5_FB_1","E16_5_FB_2",
	       "E16_5_ST_1","E16_5_ST_2",
	       "E16_5_HT_1","E16_5_HT_2",
	       "E16_5_IT_1","E16_5_IT_2",
	       "E16_5_KD_1","E16_5_KD_2",
	       "E16_5_LG_1","E16_5_LG_2",
	       "E16_5_LV_1","E16_5_LV_2",	   
	       #P0
	       "P0_FB_1","P0_FB_2",
	       "P0_HB_1","P0_HB_2",
	       "P0_HT_1","P0_HT_2",
	       "P0_IT_1","P0_IT_2",
	       "P0_KD_1","P0_KD_2",
	       "P0_LG_1","P0_LG_2",
	       "P0_MB_1","P0_MB_2",
	       "P0_ST_1","P0_ST_2",
	       "P0_LV_1","P0_LV_2");

`mkdir -p final/`;

print "Calling mCH domain for all samples: ";

`Rscript call_changepoint.R`;
foreach my $sample (@samples)
{
    my $input_file = "tmp/mCH_segment_$sample\.bed";
    open(DATA,"$input_file") or die "can't open $input_file!\n";
    my %cand;
    while(<DATA>)
    {
	chomp;
	next if length($_) == 0;
	my ($chr) = split(/\t/,$_);
	push @{$cand{$chr}},$_;
    }
    close(DATA);

    open(OUT,">final/mCHdomain_$sample\.bed") or die "can't create output for $sample!\n";
    ## get mCH domain
    foreach my $chr (keys(%cand))
    {
	my @cand_chr = @{$cand{$chr}};
	if(@cand_chr == 1)
	{
	    next;
	}
	elsif(@cand_chr == 2)
	{
	    next;
	}
	elsif(@cand_chr > 2)
	{
	    my @mCH_array;
	    foreach my $i (0...$#cand_chr)
	    {
		my ($chr,$start,$end,$mCH) = split(/\t/,$cand_chr[$i]);
		push @mCH_array,$mCH;
	    }

	    foreach my $i (1...($#cand_chr-1))
	    {		
		#if($mCH_array[$i] > $mCH_array[$i-1] and
		#   $mCH_array[$i] > $mCH_array[$i+1])
		if( ($mCH_array[$i] + 0.001) / ($mCH_array[$i-1] + 0.001) >= 1.5 and
		    ($mCH_array[$i] + 0.001) / ($mCH_array[$i+1] + 0.001) >= 1.5)
		{
		    print OUT "chr".$cand_chr[$i],"\t";
		    print OUT join(",",($mCH_array[$i-1],$mCH_array[$i],$mCH_array[$i+1])),"\n";
		}
	    }	    
	}
    }
    close(OUT);
}

print "done!\n";
$t = gmtime();
print $t,"\n\n";


my @all_files;
foreach my $ind (0...($#samples-1)/2)
{
    my $sample1 = $samples[$ind*2];
    my $sample2 = $samples[$ind*2+1];
    my @vals = split(/\_/,$sample1);
    my $tissue = join("\_",@vals[0...($#vals-1)]);
	
    my $cmd = "intersectBed -a final/mCHdomain_$sample1\.bed ";
    $cmd .= "-b final/mCHdomain_$sample2\.bed|cut -f 1-3";
    $cmd .= "awk \'(\$3-\$2)>10000\' > final/mCHdomain_$tissue\.bed";
    `$cmd`;
    
    push @all_files,"final/mCHdomain_$tissue\.bed";
}

my $cmd = "cat " . join("\ ",@all_files) . " |sort -k 1,1 -k 2,2g|mergeBed -i - > final/mCHdomain_merged.bed";
`$cmd`;
