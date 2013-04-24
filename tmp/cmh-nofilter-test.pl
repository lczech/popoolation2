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


# Author: Ram Vinay Pandey
# Author: Robert Kofler


# Define the variables
my $input;
my $output="";
my $userpopulation;
my $selectpopulation="";
my $help=0;
my $test=0;
my $verbose=1;

my $mincount=0;
my $mincoverage=0;
my $usermaxcoverage=100000000;
my $minlogpvalue=0.0;
my $removetemp=0;

# --input /Users/robertkofler/pub/PoPoolation2/Walkthrough/demo-data/cmh/small-test.sync --output /Users/robertkofler/pub/PoPoolation2/Walkthrough/demo-data/cmh/small-test.cmh --population 1,2,3,4 --min-count 2 --min-coverage 4 --max-coverage 200

GetOptions(
    "input=s"	    =>\$input,
    "output=s"	    =>\$output,
    "population=s"  =>\$userpopulation,
    "select-population=s"  =>\$selectpopulation,
    "remove-temp"   =>\$removetemp,
    "test"          =>\$test,
    "help"	    =>\$help
) or pod2usage(-msg=>"Wrong options",-verbose=>1);

pod2usage(-verbose=>2) if $help;
CMHTest::runTests() && exit if $test;
pod2usage(-msg=>"A input file has to be provided\n",-verbose=>1) unless -e $input;
pod2usage(-msg=>"A output file has to be provided\n",-verbose=>1) unless $output;
pod2usage(-msg=>"The pairwise comparisions have to be provided (--population)",-verbose=>1) unless $userpopulation;


################# write param file

my $paramfile=$output.".params";
open my $pfh, ">",$paramfile or die "Could not open $paramfile\n";
print $pfh "Using input\t$input\n";
print $pfh "Using output\t$output\n";
print $pfh "Using min-count\t$mincount\n";
print $pfh "Using min-coverage\t$mincoverage\n";
print $pfh "Using max-coverage\t$usermaxcoverage\n";
print $pfh "Using population\t$userpopulation\n";

if ($selectpopulation) {
	print $pfh "Using select-population\t$selectpopulation\n";
}
print $pfh "Using min-logpvalue\t$minlogpvalue\n";
print $pfh "Remove temporary files\t$removetemp\n";
print $pfh "Using test\t$test\n";
print $pfh "Using help\t$help\n";
close $pfh;

my $selected_population = [];
if ($selectpopulation) {
	$selected_population = resolve_selected_populations($input,$selectpopulation);
}


my $maxcoverage=get_max_coverage($input,$usermaxcoverage);
my $populations=CMHUtil::resolve_population($userpopulation);
CMHUtil::resolve_selected_user_population($populations,$selected_population);

my $syncparser;
if ($selectpopulation) {
	$syncparser=get_sumsnp_synparser_selected_pop($mincount,$mincoverage,$maxcoverage,$selected_population);
}
else {
	$syncparser=get_sumsnp_synparser($mincount,$mincoverage,$maxcoverage);
}

#my $syncparser=get_sumsnp_synparser($mincount,$mincoverage,$maxcoverage);

my $rinput=$output.".rin";
my $routput=$output.".rout";

print "Reading sync file and writing temporary R output file\n";
CMHUtil::write_Rinput($input,$rinput,$syncparser,$populations);

print "Calling R, to calculate the Cochran-Mantel-Haenszel test statistic\n";
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
			 # next if $pvalue <  $minlogpvalue; NO FILTERING
			print $ofh $line."\t".$pvalue."\n";
		}
		close $ofh;
		close $ifh;
	}
	
	sub resolve_population
	{
		my $userpopulation=shift;
		die "At least two pairwise comparisions need to be specified, e.g.: 1-2,3-4" unless $userpopulation=~m/,/;
		die "At least two pairwise comparisions need to be specified, e.g.: 1-2,3-4" unless $userpopulation=~m/-/;
		
		
		my $populations=[];
		my @temp=split /,/,$userpopulation;
		foreach my $t (@temp)
		{
			my @a=split /-/,$t;
			die "At least two pairwise comparisions need to be specified" if scalar(@a) !=2;
			push @$populations,$a[0];
			push @$populations,$a[1];
		}
		
		die "Pairwise comparisions must be an even number (user provided $userpopulation)" if(scalar(@$populations) % 2);
		return $populations;
	}
	
	sub resolve_selected_user_population
	{
		my $user_populations=shift;
		my $selected_populations=shift;
		
		if (scalar(@$selected_populations)>0) {
		
			die "--select-population parameter should have exactly same populations as given in --population parameter (user provided --select-population $selectpopulation)" if(scalar(@$selected_populations) % 2);
			
			
			for my $pop (@$user_populations) {
				my $popct=0;
				for my $pop1 (@$selected_populations) {
					if ($pop==$pop1) {
						$popct++;
					}
				}
				
				die "--select-population parameter and --population parameter does not have same populations (user provided --population $userpopulation and --select-population $selectpopulation)" if($popct < 1 || $popct>1);
				$popct=0;
			}
		
		}

	}
	
	sub write_mantelro
	{
		my $fh=shift;
print $fh <<PERLSUCKS;
mantelro<-function (x, y = NULL, z = NULL, alternative = c("two.sided", 
                                                 "less", "greater"), correct = TRUE, exact = FALSE, conf.level = 0.95) 
{
  DNAME <- deparse(substitute(x))
  if (is.array(x)) {
    if (length(dim(x)) == 3L) {
      if (any(is.na(x))) 
        stop("NAs are not allowed")
      if (any(dim(x) < 2L)) 
        stop("each dimension in table must be >= 2")
    }
    else stop("'x' must be a 3-dimensional array")
  }
  else {
    if (is.null(y)) 
      stop("if 'x' is not an array, 'y' must be given")
    if (is.null(z)) 
      stop("if 'x' is not an array, 'z' must be given")
    if (any(diff(c(length(x), length(y), length(z))) != 0L)) 
      stop("'x', 'y', and 'z' must have the same length")
    DNAME <- paste(DNAME, "and", deparse(substitute(y)), 
                   "and", deparse(substitute(z)))
    OK <- complete.cases(x, y, z)
    x <- factor(x[OK])
    y <- factor(y[OK])
    if ((nlevels(x) < 2L) || (nlevels(y) < 2L)) 
      stop("'x' and 'y' must have at least 2 levels")
    else x <- table(x, y, z[OK])
  }
  if (any(apply(x, 3L, sum) < 2)) 
    stop("sample size in each stratum must be > 1")
  I <- dim(x)[1L]
  J <- dim(x)[2L]
  K <- dim(x)[3L]
  if ((I == 2) && (J == 2)) {
    alternative <- match.arg(alternative)
    if (!missing(conf.level) && (length(conf.level) != 1 || 
                                   !is.finite(conf.level) || conf.level < 0 || conf.level > 
                                   1)) 
      stop("'conf.level' must be a single number between 0 and 1")
    NVAL <- c(`common odds ratio` = 1)
    if (!exact) {
      s.x <- apply(x, c(1L, 3L), sum)
      s.y <- apply(x, c(2L, 3L), sum)
      n <- as.double(apply(x, 3L, sum))
      DELTA <- sum(x[1, 1, ] - s.x[1, ] * s.y[1, ]/n)
      YATES <- if (correct && (abs(DELTA) >= 0.5)) 
        0.5
      else 0
      STATISTIC <- ((abs(DELTA) - YATES)^2/sum(apply(rbind(s.x, 
                                                           s.y), 2L, prod)/(n^2 * (n - 1))))
      PARAMETER <- 1
      if (alternative == "two.sided") 
        PVAL <- -pchisq(STATISTIC, PARAMETER, lower.tail = FALSE,log.p=TRUE)*log10(exp(1))
      else {
        z <- sign(DELTA) * sqrt(STATISTIC)
        PVAL <- -pnorm(z, log.p=TRUE,lower.tail = (alternative == 
                                         "less"))*log10(exp(1))
      }
      names(STATISTIC) <- "Mantel-Haenszel X-squared"
      names(PARAMETER) <- "df"
      METHOD <- paste("Mantel-Haenszel chi-squared test", 
                      if (YATES) 
                        "with"
                      else "without", "continuity correction")
      s.diag <- sum(x[1L, 1L, ] * x[2L, 2L, ]/n)
      s.offd <- sum(x[1L, 2L, ] * x[2L, 1L, ]/n)
      ESTIMATE <- s.diag/s.offd
      sd <- sqrt(sum((x[1L, 1L, ] + x[2L, 2L, ]) * x[1L, 
                                                     1L, ] * x[2L, 2L, ]/n^2)/(2 * s.diag^2) + sum(((x[1L, 
                                                                                                       1L, ] + x[2L, 2L, ]) * x[1L, 2L, ] * x[2L, 1L, 
                                                                                                                                              ] + (x[1L, 2L, ] + x[2L, 1L, ]) * x[1L, 1L, ] * 
                                                                                                      x[2L, 2L, ])/n^2)/(2 * s.diag * s.offd) + sum((x[1L, 
                                                                                                                                                       2L, ] + x[2L, 1L, ]) * x[1L, 2L, ] * x[2L, 1L, 
                                                                                                                                                                                              ]/n^2)/(2 * s.offd^2))
      CINT <- switch(alternative, less = c(0, ESTIMATE * 
                                             exp(qnorm(conf.level) * sd)), greater = c(ESTIMATE * 
                                                                                         exp(qnorm(conf.level, lower.tail = FALSE) * sd), 
                                                                                       Inf), two.sided = {
                                                                                         ESTIMATE * exp(c(1, -1) * qnorm((1 - conf.level)/2) * 
                                                                                                          sd)
                                                                                       })
      RVAL <- list(statistic = STATISTIC, parameter = PARAMETER, 
                   p.valuelog = PVAL)
    }
    else {
      METHOD <- paste("Exact conditional test of independence", 
                      "in 2 x 2 x k tables")
      mn <- apply(x, c(2L, 3L), sum)
      m <- mn[1L, ]
      n <- mn[2L, ]
      t <- apply(x, c(1L, 3L), sum)[1L, ]
      s <- sum(x[1L, 1L, ])
      lo <- sum(pmax(0, t - n))
      hi <- sum(pmin(m, t))
      support <- lo:hi
      dc <- .C(C_d2x2xk, as.integer(K), as.double(m), as.double(n), 
               as.double(t), d = double(hi - lo + 1))\$d
      logdc <- log(dc)
      dn2x2xk <- function(ncp) {
        if (ncp == 1) 
          return(dc)
        d <- logdc + log(ncp) * support
        d <- exp(d - max(d))
        d/sum(d)
      }
      mn2x2xk <- function(ncp) {
        if (ncp == 0) 
          return(lo)
        if (ncp == Inf) 
          return(hi)
        sum(support * dn2x2xk(ncp))
      }
      pn2x2xk <- function(q, ncp = 1, upper.tail = FALSE) {
        if (ncp == 0) {
          if (upper.tail) 
            return(as.numeric(q <= lo))
          else return(as.numeric(q >= lo))
        }
        if (ncp == Inf) {
          if (upper.tail) 
            return(as.numeric(q <= hi))
          else return(as.numeric(q >= hi))
        }
        d <- dn2x2xk(ncp)
        if (upper.tail) 
          sum(d[support >= q])
        else sum(d[support <= q])
      }
      PVAL <- switch(alternative, less = pn2x2xk(s, 1), 
                     greater = pn2x2xk(s, 1, upper.tail = TRUE), two.sided = {
                       relErr <- 1 + 10^(-7)
                       d <- dc
                       sum(d[d <= d[s - lo + 1] * relErr])
                     })
      mle <- function(x) {
        if (x == lo) 
          return(0)
        if (x == hi) 
          return(Inf)
        mu <- mn2x2xk(1)
        if (mu > x) 
          uniroot(function(t) mn2x2xk(t) - x, c(0, 1))\$root
        else if (mu < x) 
          1/uniroot(function(t) mn2x2xk(1/t) - x, c(.Machine\$double.eps, 
                                                    1))\$root
        else 1
      }
      ESTIMATE <- mle(s)
      ncp.U <- function(x, alpha) {
        if (x == hi) 
          return(Inf)
        p <- pn2x2xk(x, 1)
        if (p < alpha) 
          uniroot(function(t) pn2x2xk(x, t) - alpha, 
                  c(0, 1))\$root
        else if (p > alpha) 
        {
          1/uniroot(function(t) pn2x2xk(x, 1/t) - alpha, 
                    c(.Machine\$double.eps, 1))\$root
        }
        else 1
      }
      ncp.L <- function(x, alpha) {
        if (x == lo) 
          return(0)
        p <- pn2x2xk(x, 1, upper.tail = TRUE)
        if (p > alpha) 
          uniroot(function(t) pn2x2xk(x, t, upper.tail = TRUE) - 
                    alpha, c(0, 1))\$root
        else if (p < alpha) 
          1/uniroot(function(t) pn2x2xk(x, 1/t, upper.tail = TRUE) - 
                      alpha, c(.Machine\$double.eps, 1))\$root
        else 1
      }
      CINT <- switch(alternative, less = c(0, ncp.U(s, 
                                                    1 - conf.level)), greater = c(ncp.L(s, 1 - conf.level), 
                                                                                  Inf), two.sided = {
                                                                                    alpha <- (1 - conf.level)/2
                                                                                    c(ncp.L(s, alpha), ncp.U(s, alpha))
                                                                                  })
      STATISTIC <- c(S = s)
      RVAL <- list(statistic = STATISTIC, p.value = PVAL)
    }
    names(ESTIMATE) <- names(NVAL)
    attr(CINT, "conf.level") <- conf.level
    RVAL <- c(RVAL, list(conf.int = CINT, estimate = ESTIMATE, 
                         null.value = NVAL, alternative = alternative))
  }
  else {
    df <- (I - 1) * (J - 1)
    n <- m <- double(length = df)
    V <- matrix(0, nrow = df, ncol = df)
    for (k in 1:K) {
      f <- x[, , k]
      ntot <- sum(f)
      rowsums <- apply(f, 1L, sum)[-I]
      colsums <- apply(f, 2L, sum)[-J]
      n <- n + c(f[-I, -J])
      m <- m + c(outer(rowsums, colsums, "*"))/ntot
      V <- V + (kronecker(diag(ntot * colsums, nrow = J - 
                                 1) - outer(colsums, colsums), diag(ntot * rowsums, 
                                                                    nrow = I - 1) - outer(rowsums, rowsums))/(ntot^2 * 
                                                                                                                (ntot - 1)))
    }
    n <- n - m
    STATISTIC <- c(crossprod(n, qr.solve(V, n)))
    PARAMETER <- df
    PVAL <- -pchisq(STATISTIC, PARAMETER, lower.tail = FALSE,log.p=TRUE)*log10(exp(1))
    names(STATISTIC) <- "Cochran-Mantel-Haenszel M^2"
    names(PARAMETER) <- "df"
    METHOD <- "Cochran-Mantel-Haenszel test"
    RVAL <- list(statistic = STATISTIC, parameter = PARAMETER, 
                 p.valuelog = PVAL)
  }
  RVAL <- c(RVAL, list(method = METHOD, data.name = DNAME))
  class(RVAL) <- "htest"
  return(RVAL)
}
PERLSUCKS

	}
	
	sub write_Rinput
	{
		my $syncfile=shift;
		my $rinput=shift;
		my $syncparser=shift;
		my $populations=shift;
		
		my $third_dim =int(scalar(@$populations)/2);
		my $dim_str="c(2,2,$third_dim)";
		
		open my $ifh, "<", $syncfile or die "Could not open input file";
		open my $ofh, ">", $rinput or die "Could not open routput file";
		write_mantelro($ofh);
		while(my $line=<$ifh>)
		{
			chomp $line;
			my $e=$syncparser->($line);
			# next unless $e->{ispuresnp}; # NO FIlTERING for any kind of SNPs
			
			my $pop_str=_get_populationstring($e->{samples},$populations);
			my $ar_str="array($pop_str,dim=$dim_str)";
			my $mantel_str="mantelro($ar_str,alternative=c(\"two.sided\"))\$p.valuelog";
			print $ofh "print(\"$line\")\n";
			print $ofh $mantel_str."\n";
		}
		close $ofh;
		close $ifh;
	}
	
	
	sub _get_populationstring
	{
		my $samples=shift;
		my $populations=shift;
		
		my ($major,$minor) = MajorAlleles::get_major_minor_alleles($samples);
		my @ar=();
		for(my $i=1; $i<@$populations; $i+=2)
		{
			my $basenr=$populations->[$i-1];
			my $derivednr=$populations->[$i];
			$basenr--; $derivednr--;
			my $base=$samples->[$basenr];
			my $derived=$samples->[$derivednr];
			
			# set default to 1; otherwise cmh-test crashes..
			$base->{$major}=1 if $base->{$major}==0;
			$base->{$minor}=1 if $base->{$minor}==0;
			$derived->{$major}=1 if $derived->{$major}==0;
			$derived->{$minor}=1 if $derived->{$minor}==0;
			
			push @ar,$base->{$major};
			push @ar,$derived->{$major};
			push @ar,$base->{$minor};
			push @ar,$derived->{$minor};
		}
		my $string_all_allele = join(",",@ar);
		my $popstring="c($string_all_allele)";
		return $popstring;
	}
    
	
    
}




{
    package CMHTest;
    use FindBin qw/$RealBin/;
    use lib "$RealBin/Modules";
    use Test;
    use Pileup;
    use Test::TSynchronized;
    use Test::TMaxCoverage;
    use MajorAlleles;
    

    sub runTests
    {
        run_MaxCoverageTests();
        run_SynchronizedTests();
	run_majorminor();

        exit;
    }
    
    sub run_majorminor
    {
	my($maj,$min);
	
	($maj,$min)=get_major_minor_alleles([{A=>11,T=>0,C=>0,G=>0},{A=>0,T=>10,C=>0,G=>0},{A=>0,T=>0,C=>9,G=>0},{A=>0,T=>0,C=>0,G=>9}]);
	is($maj,"A","identification of major allele; correct allele");
	is($min,"T","identification of minor allele; correct allele");
	
	($maj,$min)=get_major_minor_alleles([{A=>3,T=>0,C=>0,G=>0},{A=>3,T=>10,C=>0,G=>0},{A=>3,T=>0,C=>9,G=>0},{A=>2,T=>0,C=>0,G=>9}]);
	is($maj,"A","identification of major allele; correct allele");
	is($min,"T","identification of minor allele; correct allele");	
	
	($maj,$min)=get_major_minor_alleles([{A=>3,T=>0,C=>0,G=>0},{A=>3,T=>10,C=>0,G=>0},{A=>3,T=>2,C=>9,G=>0},{A=>2,T=>0,C=>0,G=>9}]);
	is($maj,"T","identification of major allele; correct allele");
	is($min,"A","identification of minor allele; correct allele");		
	
	($maj,$min)=get_major_minor_alleles([{A=>3,T=>0,C=>4,G=>0},{A=>3,T=>10,C=>0,G=>0},{A=>3,T=>2,C=>9,G=>0},{A=>2,T=>0,C=>0,G=>9}]);
	is($maj,"C","identification of major allele; correct allele");
	is($min,"T","identification of minor allele; correct allele");	
	    
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

=item B<--population>

the pairwise comparsions which will be used for the cmh-test.
Pairwise comparisions have to be separated by a C<,> and the two populations which will be compared by a C<->. For example when the user provides 1-3,2-4 the script will compare population 1 with population 3 and population 2 with population 4;
Note also comparisions involving one population for several times are possible (e.g.: --population 1-7,3-9); Mandatory parameter

=item B<--select-population>

A comma seperated list of populations. Optional parameter. This parameter should have same populations as in --population parameter.
If user user --select-population parameter then only selected populations will be used for the minimum allele count check and p-value calculations. . If not using --select-population parameter the minimum allele count will be checked for all populations and p-value will be calculated for the population pair given with --population parameter.
(e.g.: --select-population 1,7,3,9);


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
