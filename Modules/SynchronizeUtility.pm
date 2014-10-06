{
    package SynchronizeUtility;
    use strict;
    use warnings;
    use FindBin qw/$Bin/;
    use lib "$Bin";
    use Synchronized;
    use Test;
    use IO::Uncompress::Gunzip;
    
    require Exporter;
    our @ISA = qw(Exporter);
    our @EXPORT =qw(format_parsed_pileup get_popcount_forsyncfile format_synchronized syncsample2string resolve_selected_populations);
    
    
    sub syncsample2string
    {
        my $sync=shift;
        
        
        my $toret="A" x $sync->{A} . "T" x $sync->{T} . "C" x $sync->{C} . "G" x $sync->{G} . "-" x $sync->{del} . "N" x $sync->{N};
        return $toret; 
    }
    
    
    sub format_synchronized
    {
        my $sync=shift;
        
        my @temp=();
        # chr, pos, refchar, minsampcov, samples (A T C G N del)
        push @temp,$sync->{chr};
        push @temp,$sync->{pos};
        push @temp,$sync->{refchar};
        foreach my $s (@{$sync->{samples}})
        {
            
            push @temp,format_parsed_pileup($s);
        }
        
        my $toret=join("\t",@temp);
        return $toret;
    }
    
    sub format_parsed_pileup
    {
        my $pp=shift;
        return "$pp->{A}:$pp->{T}:$pp->{C}:$pp->{G}:$pp->{N}:$pp->{del}" if $pp;
        return "-";
    }
    
    sub get_popcount_forsyncfile
    {
        my $syncfile=shift;
	
	my $ifh = undef;
	if($syncfile=~/\.gz$/i) {
		$ifh = new IO::Uncompress::Gunzip $syncfile or die "Could not open file gzipped file $syncfile  $!";
	}
	else {
		open $ifh, "<", $syncfile  or die "Could not open file handle, $!";
	}
	
        my $pp=get_basic_syncparser();
        my $firstline=<$ifh>;
        chomp $firstline;
        my $pl=$pp->($firstline);
        my $count=scalar(@{$pl->{samples}});
        close $ifh;
        return $count;
        
        
    }
    
    sub resolve_selected_populations
    {
        my $syncfile=shift;
        my $selected_pop=shift;
        
        my $popCpunt = get_popcount_forsyncfile($syncfile);
        
        my $populations=[];
	my @temp=split /,/,$selected_pop;
		foreach my $t (@temp)
		{
                    chomp($t);
		    push @$populations,$t;
		}
        @$populations = sort { $a <=> $b } @$populations;
        
        my $first_pop=0;
        my $last_pop=0;
        
        $first_pop=$populations->[0];
        $last_pop=$populations->[-1];
        if ($first_pop< 1) {
            die "Any selected population in --select-population parameter can not be smallar than 1 (user provided --select-population $selected_pop)\n";
        }
        
        if ($last_pop> $popCpunt) {
            die "Any selected population in --select-population parameter can not be greater than available populations. (user provided --select-population $selected_pop)\n There are only $popCpunt populations available in input file:$syncfile\n";
        }
 
        return $populations;

    }


}

1;