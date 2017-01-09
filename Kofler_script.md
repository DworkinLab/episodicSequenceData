#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use File::Basename;
use File::Path;
use FindBin qw/$RealBin/;
use lib "$RealBin/../Modules";
use MaxCoverage;
use Synchronized;
use SynchronizeUtility;
use MajorAlleles; # get the two major allele
use Test;


# Author: Robert Kofler


# Define the variables
my $input;
my $output="";
my $usertimeseries;
my $help=0;
my $test=0;
my $verbose=1;

my $mincount=2;
my $mincoverage=4;
my $usermaxcoverage;
my $minlogpvalue=0.0;
my $removetemp=0;

# --input /Users/robertkofler/pub/PoPoolation2/Walkthrough/demo-data/cmh/small-test.sync --output /Users/robertkofler/pub/PoPoolation2/Walkthrough/demo-data/cmh/small-test.cmh --population 1,2,3,4 --min-count 2 --min-coverage 4 --max-coverage 200

GetOptions(
    "input=s"	    =>\$input,
    "output=s"	    =>\$output,
    "min-count=s"   =>\$mincount,
    "min-coverage=i"=>\$mincoverage,
    "max-coverage=s"=>\$usermaxcoverage,
    "time-series=s"  =>\$usertimeseries,
    "min-logpvalue=f"  =>\$minlogpvalue,
    "remove-temp"   =>\$removetemp,
    "test"          =>\$test,
    "help"	    =>\$help
) or pod2usage(-msg=>"Wrong options",-verbose=>1);

pod2usage(-verbose=>2) if $help;
die "Tests currently not supported" if $test;
pod2usage(-msg=>"A input file has to be provided\n",-verbose=>1) unless -e $input;
pod2usage(-msg=>"A output file has to be provided\n",-verbose=>1) unless $output;
#pod2usage(-msg=>"Minimum coverage must be equal or larger than minimum count",-verbose=>1) unless $mincoverage>= $mincount;
pod2usage(-msg=>"Maximum coverage has to be provided",-verbose=>1) unless $usermaxcoverage;
pod2usage(-msg=>"The time series has to be provided (--time-series)",-verbose=>1) unless $usertimeseries;


################# write param file

my $paramfile=$output.".params";
open my $pfh, ">",$paramfile or die "Could not open $paramfile\n";
print $pfh "Using input\t$input\n";
print $pfh "Using output\t$output\n";
print $pfh "Using min-count\t$mincount\n";
print $pfh "Using min-coverage\t$mincoverage\n";
print $pfh "Using max-coverage\t$usermaxcoverage\n";
print $pfh "Using time-series\t$usertimeseries\n";
print $pfh "Using min-logpvalue\t$minlogpvalue\n";
print $pfh "Remove temporary files\t$removetemp\n";
print $pfh "Using test\t$test\n";
print $pfh "Using help\t$help\n";
close $pfh;




my $maxcoverage=get_max_coverage($input,$usermaxcoverage);
my $timeseries=CMHUtil::resolve_timeseries($usertimeseries);
my $syncparser=get_sumsnp_synparser($mincount,$mincoverage,$maxcoverage);

#my $syncparser=get_sumsnp_synparser($mincount,$mincoverage,$maxcoverage);

my $rinput=$output.".rin";
my $routput=$output.".rout";

print "Reading sync file and writing temporary R output file\n";
CMHUtil::write_Rinput($input,$rinput,$syncparser,$timeseries);

print "Calling R, to calculate the gls test statistic\n";
system("R --vanilla --slave <$rinput >$routput");

print "Parsing R-output and writing output file\n";
CMHUtil::write_output($routput,$output,$minlogpvalue);

if($removetemp)
{
	print "Removing temporary files\n";
	unlink($rinput);
	unlink($routput);
}
print "Done\n";

exit(0);



{
	package CMHUtil;
	use strict;
	use warnings;
	use List::Util qw[min max];
	use FindBin qw/$RealBin/;
	use lib "$RealBin/Modules";
	use MaxCoverage;
	use Synchronized;
	use SynchronizeUtility;
	use MajorAlleles; # get the two major allele
	
	sub write_output
	{
		my $routput=shift;
		my $output=shift;
		my $minlogpvalue=shift;
		
		open my $ifh,"<", $routput or die "Could not open input file\n";
		open my $ofh,">",$output or die "Could not open output file\n";
		
		while(1)
		{
			#[1] "2R\t2296\tN\t90:10:0:0:0:0\t100:0:0:0:0:0\t100:0:0:0:0:0\t100:0:0:0:0:0"
			#[1] 0.003583457
			my $line=<$ifh>;
			last unless $line;
			my $pvalue=<$ifh>;
			chomp $line; chomp $pvalue;
			$line=~s/^\S+\s//;
			$line=~s/^"//;
			$line=~s/"$//;
			$line=~s/\\t/\t/g;
			$pvalue=~s/^\S+\s//;
			#$pvalue="1.0" if $pvalue eq "NaN"; 	# stupid mantelhaenszeltest prodcues NaN for example mantelhaen.test(array(c(100,100,0,0,100,100,0,0,100,100,0,0),dim=c(2,2,3)),alternative=c("two.sided"))
								# this is clearly no differentiation thus 1.0 (necessary as it fucks up sorting by significance)
			next if $pvalue <  $minlogpvalue;
			print $ofh $line."\t".$pvalue."\n";
		}
		close $ofh;
		close $ifh;
	}
	
	sub resolve_timeseries
	{
		my $usertimeseries=shift;
		die "At least two time-series need to be specified, e.g.: 1-2,3-4" unless $usertimeseries=~m/,/;
		die "At least two time points are necessary in each time-series, e.g.: 1-2,3-4" unless $usertimeseries=~m/-/;
		
		
		my $timeseries=[];
		my @temp=split /,/,$usertimeseries;
		my $lengcomp=undef;
		foreach my $t (@temp)
		{
			my $temp=[];
			my @a=split /-/,$t;
			foreach my $asplit (@a)
			{
				push @$temp,$asplit;
			}
			$lengcomp=length(@$temp) unless $lengcomp;
			die "All time-series need to have same length" unless $lengcomp == length(@$temp);
			push @$timeseries,$temp;
		}
		
		## $timeseries->[0] 1-2-3-4-5-6
		##  $timerseries->[1] 7-8-9-10-11-12-13-14
		return $timeseries;
	}
	

	
	sub _write_common
	{
		my $fh=shift;
		my $timeseriescount=shift;
		my $timeserieslength=shift;
print $fh <<PERLSUCKS;
library(nlme)
tscount=$timeseriescount
tslength=$timeserieslength
time <- factor(c(0:(tslength-1)))
# e.g.: factor(c(0:6)) = Levels: 0 1 2 3 4 5 6 
PERLSUCKS

	}
	
	sub write_Rinput
	{
		my $syncfile=shift;
		my $rinput=shift;
		my $syncparser=shift;
		my $timeseries=shift;
		
		
		
		my $timeseriescount=@$timeseries;
		my $timeserieslength=@{$timeseries->[0]};
		my $flattimeseries=_flaten_time_series($timeseries);
		
		open my $ifh, "<", $syncfile or die "Could not open input file";
		open my $ofh, ">", $rinput or die "Could not open routput file";
		_write_common($ofh,$timeseriescount,$timeserieslength);
		while(my $line=<$ifh>)
		{
			chomp $line;
			my $e=$syncparser->($line);
			next unless $e->{ispuresnp};
			_write_single_gls($ofh,$line,$e,$flattimeseries,$timeseriescount,$timeserieslength,$timeseries);
		}
		close $ofh;
		close $ifh;
	}
	
	
	sub _update_freqar_old
	{
		# old version which looks for identical timeseries
		my $freqar=shift;
		my $timeseries=shift;
		
		
		my $isidentical=1;
		my $timeseriescount=@$timeseries;
		my $timeserieslength=@{$timeseries->[0]};
		for(my $i=0; $i<$timeserieslength; $i++)
		{
			my $val=undef;
			for(my $k=0; $k<$timeseriescount; $k++)
			{
				my $index=$i+$timeserieslength*$k; 
				my $activeval=$freqar->[$index];
				$val=$activeval unless(defined($val));
				if($val != $activeval)
				{
					$isidentical=0;
				}

			}
		}
		
		if($isidentical)
		{
			$freqar->[0] += 0.001
		}
		return $freqar;

	}
	
	
	sub _is_fixed
	{
		# check if last element in timeseries is always fixed
		my $freqar=shift;
		my $timeseries=shift;
		
		
		my $isfixed=1;
		my $timeseriescount=@$timeseries;
		my $timeserieslength=@{$timeseries->[0]};
		

		
		for(my $k=0; $k<$timeseriescount; $k++)
		{
			my $index=($timeserieslength-1)+$timeserieslength*$k; 
			my $activeval=$freqar->[$index];
			$isfixed=0 if($activeval!=1)
		}
		return $isfixed;
	}
	
	
	
	sub _update_freqar
	{
		# new version which just updates it
		my $freqar=shift;
		my $timeseries=shift;
		
		
		my $timeseriescount=@$timeseries;
		my $timeserieslength=@{$timeseries->[0]};
	
		for(my $k=0; $k<$timeseriescount; $k++)
		{
			my $index=$timeserieslength*$k;
			my $mod=$index % 3;
			if($mod==0){
				$freqar->[$index]+=0.001;
			}
			if($mod==2)
			{
				$freqar->[$index]-=0.001;
			}
		}
		
		return $freqar;
	}
	
	
	sub _flaten_time_series
	{
		my $timeseries=shift;
		my $flattimeseries=[];
		foreach my $temp (@$timeseries)
		{
			foreach my $t (@$temp)
			{
				push @$flattimeseries,$t;
			}
		}
		return $flattimeseries;
	}
	
	
	sub _write_single_gls
	{
		my $ofh=shift;
		my $line=shift;
		my $entry=shift;
		my $flattimeseries=shift;
		my $timeseriescount=shift;
		my $timeserieslength=shift;
		my $truedimtimeseries=shift;
		
		
		my ($freqar,$covar) = _get_freq_cov_ar($entry->{samples},$flattimeseries);
		#  $freqar=_update_freqar($freqar,$tdtimeseries);
		my $isfixed=_is_fixed($freqar,$truedimtimeseries);
		
		if($isfixed)
		{
			# ROUT
	# [1] "2L\t89\tA\t1942:0:58:0:0:0\t1978:0:22:0:0:0\t1973:0:27:0:0:0\t1972:0:28:0:0:0\t1981:0:19:0:0:0\t1998:0:2:0:0:0\t2000:0:0:0:0:0\t1942:0:58:0:0:0\t1891:0:109:0:0:0\t1907:0:93:0:0:0\t1898:0:102:0:0:0\t1911:0:89:0:0:0\t1908:0:92:0:0:0\t1904:0:96:0:0:0\t1942:0:58:0:0:0\t1971:0:29:0:0:0\t1989:0:11:0:0:0\t1991:0:9:0:0:0\t1988:0:12:0:0:0\t2000:0:0:0:0:0\t2000:0:0:0:0:0"
	# [1] 0.06554403
	print $ofh <<PERLSUCKS;
print("$line")
print(0.0)
PERLSUCKS
		}
		else
		{
			$freqar=_update_freqar($freqar,$truedimtimeseries);
			my ($freqstr,$covstr) = _get_freq_cov_string($freqar,$covar);
			print $ofh <<PERLSUCKS;
print("$line")
data.all <- data.frame(props=$freqstr, ns=$covstr, times=rep(time,tscount), time.n=rep(c(0:(tslength-1)),tscount), reps=factor(rep(1:tscount,each=tslength)) )
res0a<- gls(props~1,correlation=corAR1(form=~time.n|reps),weights=varFixed(~(1/ns)),data=data.all,method="ML")
res1a<- gls(props~times,correlation=corAR1(form=~time.n|reps),weights=varFixed(~(1/ns)),data=data.all,method="ML")
print(-log(anova(res0a,res1a)\$"p-value"[2])/log(10)) # -log10()! don't forget
PERLSUCKS
		}
		
		
		
		

		
		



	}
	
	sub _get_freq_cov_ar
	{
		my $samples=shift;
		my $timeseries=shift;
		my ($major,$minor) = MajorAlleles::get_major_minor_alleles($samples);
		
		my $freqar=[];
		my $covar=[];
		foreach my $popnr (@$timeseries)
		{
			my $majcount=$samples->[$popnr-1]{$major};
			my $mincount=$samples->[$popnr-1]{$minor};
			my $cov=$majcount+$mincount;
			my $freq=$majcount/$cov;
			
			push @$freqar,$freq;
			push @$covar,$cov;
		}
		return ($freqar,$covar);
	}
	
	
	
	sub _get_freq_cov_string
	{
		my $freqar=shift;
		my $covar=shift;
		
		my $str_freq = join(",",@$freqar);
		my $str_cov = join(",",@$covar);
		
		my $freqstring="c($str_freq)";
		my $covstring="c($str_cov)";
		return ($freqstring,$covstring);
	}
	
	
	
    
}








=head1 NAME

cmh-test.pl - This script calculates the Cochran-Mantel-Haenszel test for each SNP  

=head1 SYNOPSIS

 perl cmh-test.pl --input input.sync --output output.cmh --min-count 2 --min-coverage 4 --max-coverage 1000 --population 1-2,3-4,5-6 --remove-temp --select-population 1,2,3,4,5,6

=head1 OPTIONS

=over 4

=item B<--input>

The input file has to be synchronized pileup file. Mandatory parameter

=item B<--output>

The output file. Mandatory parameter

=item B<--min-count>

the minimum count of the minor allele; used for SNP identification. SNPs will be identified considering all populations simultanously. default=2
The minimum count may be provided as one of the following two ways:

 '2' The minimum count of the minor allele should be at least 2 or higher by considering all populations OR selected populations given in parameter --population
 '2%' The minimum minor allele frequency should be at least 2% or higher by considering all populations OR selected populations given in parameter --population
If user uses the --select-population parameter then minimum allele count will be checked for selected populations only. If not using --select-population the minimum allele count will be checked for all populations.
Please note that when you are using min-count in percent then the run time will increase because for each locus program will calculate frequency of 4 alleles.

=item B<--min-coverage>

the minimum coverage; used for SNP identification, the coverage in ALL selected populations has to be higher or equal to this threshold, otherwise no SNP will be called. default=4

=item B<--max-coverage>

The maximum coverage; All selected populations are required to have coverages lower or equal than the maximum coverage; Mandatory
The maximum coverage may be provided as one of the following:

 '500' a maximum coverage of 500 will be used for all populations
 '300,400,500' a maximum coverage of 300 will be used for the first population, a maximum coverage of 400 for the second population and so on
 '2%' the 2% highest coverages will be ignored, this value is independently estimated for every population
  
=item B<--min-logpvalue>

the minimum -log10(p-value) cut off  to filter all snp with > min-pvalue cutoff; default=0.0 [Optional parameter]
For example if the user provides 5 than all pvalues > 0.00001 will be discarded;  

=item B<--time-series>

the time-series to be tested.


Two time series have to be seperated by a C<,> and the populations within a time series by a C<->.
As an example: 1-2-3-4-5-6-7,8-9-10-11-12-13-14
At least two time series have to be provided; All time-series have to have the same length
Mandatory parameter


=item B<--remove-temp>

flag; remove the temporary files at the end; default=off

=item B<--test>

Run the unit tests for this script. 

=item B<--help>

Display help for this script

=back

=head1 Details

=head2 Input

Input is a single tab delimited file which contains a lightwight representation of every pileup file.
Every pileup file represents a population and will be parsed into a list of A-count:T-count:C-count:G-count:N-count:*-count

 2L	5002	G	0:0:0:17:0:0	0:0:0:28:0:0	0:0:0:31:0:0	0:0:0:35:0:0	0:1:0:33:0:0	0:3:0:31:0:0
 2L	5009	A	16:0:0:0:0:0	26:0:0:0:0:0	29:0:1:0:0:0	36:0:0:0:0:0	34:0:0:0:0:0	32:0:1:0:0:0
 2L	5233	G	0:0:5:46:0:0	0:0:0:43:0:0	0:0:0:60:0:0	0:0:3:61:0:0	0:0:0:56:0:0	0:0:0:48:0:0
 
 col 1: reference chromosome
 col 2: position in the reference chromosome
 col 3: reference genome base
 col 4: population 1
 col 5: population 2
 col n: population n-3
 
 population data are in the form
 A:T:C:G:N:*
 A: count of character A
 T: count of character T
 C: count of character C
 G: count of character G
 N: count of character N
 *: deletion, count of deletion
 
=head2 Output

 2L	5002	G	0:0:0:17:0:0	0:0:0:28:0:0	0:0:0:31:0:0	0:0:0:35:0:0	0:1:0:33:0:0	0:3:0:31:0:0	0.609
 2L	5009	A	16:0:0:0:0:0	26:0:0:0:0:0	29:0:1:0:0:0	36:0:0:0:0:0	34:0:0:0:0:0	32:0:1:0:0:0	0.957
 2L	5233	G	0:0:5:46:0:0	0:0:0:43:0:0	0:0:0:60:0:0	0:0:3:61:0:0	0:0:0:56:0:0	0:0:0:48:0:0	0.8088


 col 1: reference chromosome
 col 2: position in the reference chromosome
 col 3: reference genome base
 col 4: population 1
 col 5: population 2
 col n: population n-3
 col n+1: cmh -log10(p-value)
 Note: If user gives --population  1-13,2-6,3-7 and --select-population 1,13,2,6,3,7 then SNP calling and p-value will be calculated only for selected populations but still all population will be printed in output file just to keep all sync file information.

=head1 Technical details

This script identifies the two major alleles for every SNPs and than runs the run Cochran mental haenszel test.
The script creates two temporary output files for R, having the extensions C<.rin> and C<.rout> which may be automatically removed using the option C<--remove-temp>.
Also note that the CMH test sometimes produces the pvalue 'NaN' (eg: mantelhaen.test(array(c(100,100,0,0,100,100,0,0,100,100,0,0),dim=c(2,2,3)),alternative=c("two.sided")))
This NaN will be replaced by 1.0


=head2 Test statistic

You use the Cochran–Mantel–Haenszel test (which is sometimes called the Mantel–Haenszel test) for repeated tests of independence.
There are three nominal variables; you want to know whether two of the variables are independent of each other, and the third variable identifies the repeats.
For more details: http://udel.edu/~mcdonald/statcmh.html

=head1 AUTHORS

 Ram Vinay Pandey
 Robert Kofler
 Viola Nolte
 Christian Schloetterer

=cut
