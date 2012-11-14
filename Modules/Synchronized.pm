{
    package Synchronized;
    use strict;
    use warnings;
    use FindBin qw/$Bin/;
    use lib "$Bin";
    use Test;
    use List::Util qw[min max];
    require Exporter;
    our @ISA = qw(Exporter);
    our @EXPORT     =qw(get_basic_syncparser get_sumsnp_synparser get_sumsnp_synparser_selected_pop);
    
    
    #my $sp=get_basic_syncparser(),
    #my $line
    #my $pl=$sp->($line);
    # $pl->{chr}
    # $pl->{samples}[0]{A}
    
    
    sub get_basic_syncparser
    {
        return sub
        {
            my $line=shift;
            chomp $line;
            my @a=split /\s+/,$line;
            my $chr=shift @a;
            my $pos=shift @a;
            my $rc=shift @a;
            $rc=uc($rc);
            
            
            # Parsing
            my @samp;
            for(my $i=0; $i<@a; $i++)
            {
                my $col=$a[$i];
                my $e;
                
                if($col=~/-/)
                {
                    $e={index=>$i,eucov=>0,totcov=>0,A=>0,T=>0,C=>0,G=>0,N=>0,del=>0};
                }
                else
                {
                    my @parts=split /:/,$col;
                    die "failed parsing $col; does not have the correct number of entries" unless @parts ==6;
                    $e = {
                         A=>$parts[0],
                         T=>$parts[1],
                         C=>$parts[2],
                         G=>$parts[3],
                         N=>$parts[4],
                         del=>$parts[5]
                        };
                    
                    $e->{eucov}  = ($e->{A}+$e->{T}+$e->{C}+$e->{G});
                    $e->{totcov} = ($e->{eucov}+$e->{N}+$e->{del});
                    $e->{index}  = $i;
                }
                push @samp,$e;
            }
            
            my $minsampcov=min(map {$_->{eucov}} @samp);
            
            my $en={
            chr=>$chr,
            pos=>$pos,
            refchar=>$rc,
            minsampcov=>$minsampcov,
            samples=>\@samp
            };
            return $en;
        # chr, pos, refchar, minsampcov, samples (A T C G N del)
        }
    }
    
    sub get_sumsnp_synparser #mincount, mincoverage, maxcoverage
    {
        my $mincount=shift;
        my $mincoverage=shift;
        my $maxcoverage=shift;
        
        my $bp=get_basic_syncparser();
        
        return sub
        {
            my $line=shift;
            my $p=$bp->($line);
            
            
            # check if the position is a snp
            my @samp=@{$p->{samples}};
            my $issnp=0;
            my $taintedsnp=0;
            my ($ca,$ct,$cc,$cg,$cn,$cdel)=(0,0,0,0,0,0);
            my $is_suficient_covered=1;
            my $allelecount=0;
            
            for(my $i=0; $i<@samp; $i++)
            {
                my $s=$samp[$i];
                my $maxcov=$maxcoverage->[$i];
                my $iscov=1;
                $iscov=0 if($s->{eucov}<$mincoverage || $s->{eucov}>$maxcov);
                $s->{iscov}=$iscov;
                $ca+= $s->{A};
                $ct+= $s->{T};
                $cc+= $s->{C};
                $cg+= $s->{G};
                $cn+= $s->{N};
                $cdel+= $s->{del};
                #print "$s->{T}\n";
                # also reset the general appropriate coverage unless every entry has the correct coverage;
                $is_suficient_covered=0 unless $iscov;
                
            }
            
            if($mincount=~/%$/) {
        	my $mincount1=$mincount;
                
                $mincount1=~s/%$//;
                chomp($mincount1);
                my $total_cov = $ca+$ct+$cc+$cg;
                my ($a_percent,$t_percent,$c_percent,$g_percent,$del_percent)=(0,0,0,0,0);

                if ($total_cov>0) {
                    $a_percent = ($ca/$total_cov)*100;
                    $t_percent = ($ct/$total_cov)*100;
                    $c_percent = ($cc/$total_cov)*100;
                    $g_percent = ($cg/$total_cov)*100;
                    $del_percent = ($cdel/$total_cov)*100;
                }
                
                
                $allelecount++ if $a_percent >= $mincount1;
                $allelecount++ if $t_percent >= $mincount1;
                $allelecount++ if $c_percent >= $mincount1;
                $allelecount++ if $g_percent >= $mincount1;
                
                # the SNP is tainted if there are any deletions at the given position
                $taintedsnp=($del_percent >= $mincount1)? 1:0;
               
                
            }
            else {
                $allelecount++ if $ca >= $mincount;
                $allelecount++ if $ct >= $mincount;
                $allelecount++ if $cc >= $mincount;
                $allelecount++ if $cg >= $mincount;
                
                # the SNP is tainted if there are any deletions at the given position
                $taintedsnp=($cdel >= $mincount)? 1:0;
            }

            $issnp=1 if $allelecount > 1;          
        
        # unset the snp if not sufficiently covered
        $issnp=0 unless $is_suficient_covered; # no SNP will ever be allowed at a position which is not sufficiently covered in all data files!!

        
        $p->{ignore}=$taintedsnp;
        $p->{iscov}=$is_suficient_covered;
        $p->{issnp}=$issnp;
        $p->{ispuresnp}=$p->{issnp};
        if($taintedsnp)
        {
            $p->{iscov}=0;
            $p->{ispuresnp}=0;
            $p->{issnp}=0;
        }
        
        # iscov, issnp, ispuresnps
        # chr, pos, refchar, minsampcov, samples (A T C G N del)
        return $p;
        };
    }
    
    sub get_sumsnp_synparser_selected_pop #mincount, mincoverage, maxcoverage, requested populations 
    {
        my $mincount=shift;
        my $mincoverage=shift;
        my $maxcoverage=shift;
        my $populations=shift;
        
        
        my $bp=get_basic_syncparser();
        
        return sub
        {
            my $line=shift;
            my $p=$bp->($line);
            
            
            foreach (keys %$p) {
            	#print "$_\n";
            	
            }
            # check if the position is a snp
            my @samp=@{$p->{samples}};
            my $issnp=0;
            my $taintedsnp=0;
            my ($ca,$ct,$cc,$cg,$cn,$cdel)=(0,0,0,0,0,0);
            my $is_suficient_covered=1;
            my $allelecount=0;
            
            ## made change on 24-10-2012 because un-used populations are also tested for min and max cov
            foreach my $popindex (@$populations) {
            	
            	my $i=$popindex-1;
            	my $s=$samp[$i];
                my $maxcov=$maxcoverage->[$i];
                my $iscov=1;
                $iscov=0 if($s->{totcov}<$mincoverage || $s->{totcov}>=$maxcov);
                $s->{iscov}=$iscov;
                $ca+= $s->{A};
                $ct+= $s->{T};
                $cc+= $s->{C};
                $cg+= $s->{G};
                $cn+= $s->{N};
                $cdel+= $s->{del};
                # also reset the general appropriate coverage unless every entry has the correct coverage;
                $is_suficient_covered=0 unless $iscov;
            	
            }
            
            #for(my $i=0; $i<@samp; $i++)
            #{
            #	my $s=$samp[$i];
            #	my $iscov=1;
            #	$s->{iscov}=$iscov;
            #    $ca+= $s->{A};
            #    $ct+= $s->{T};
            #    $cc+= $s->{C};
            #    $cg+= $s->{G};
            #    $cn+= $s->{N};
            #    $cdel+= $s->{del};
            #}
	
            
            
            
            if($mincount=~/%$/) {
        	my $mincount1=$mincount;
                
                $mincount1=~s/%$//;
                chomp($mincount1);
                my $total_cov = $ca+$ct+$cc+$cg;
                my ($a_percent,$t_percent,$c_percent,$g_percent,$del_percent)=(0,0,0,0,0);

                if ($total_cov>0) {
                    $a_percent = ($ca/$total_cov)*100;
                    $t_percent = ($ct/$total_cov)*100;
                    $c_percent = ($cc/$total_cov)*100;
                    $g_percent = ($cg/$total_cov)*100;
                    $del_percent = ($cdel/$total_cov)*100;
                }
                
                
                $allelecount++ if $a_percent >= $mincount1;
                $allelecount++ if $t_percent >= $mincount1;
                $allelecount++ if $c_percent >= $mincount1;
                $allelecount++ if $g_percent >= $mincount1;
                
                # the SNP is tainted if there are any deletions at the given position
                $taintedsnp=($del_percent >= $mincount1)? 1:0;              
                
            }
            else {
                $allelecount++ if $ca >= $mincount;
                $allelecount++ if $ct >= $mincount;
                $allelecount++ if $cc >= $mincount;
                $allelecount++ if $cg >= $mincount;
                
                # the SNP is tainted if there are any deletions at the given position
                $taintedsnp=($cdel >= $mincount)? 1:0;
            }
            $issnp=1 if $allelecount > 1;
            
            #print "$issnp\t$allelecount\n";
        # unset the snp if not sufficiently covered
        #$issnp=0 unless $is_suficient_covered; # no SNP will ever be allowed at a position which is not sufficiently covered in all data files!!
        if (($issnp>0) and ($is_suficient_covered>0)) {
        	$issnp=1;
        }
        else {
        	$issnp=0;	
        }
        # the SNP is tainted if there are any deletions at the given position
        #my $taintedsnp=($cdel >= $mincount)? 1:0;
        
        $p->{iscov}=$is_suficient_covered;
        $p->{issnp}=$issnp;
        $p->{ispuresnp}=($p->{issnp} and not $taintedsnp)?1:0;
        
        # iscov, issnp, ispuresnps
        # chr, pos, refchar, minsampcov, samples (A T C G N del)
        return $p;
        };
    }
    
}
1;