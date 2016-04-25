package HmmTools;
require Exporter;
use strict;
use warnings;
use Carp;
use Data::Dumper;
our @ISA    = qw (Exporter);
our @EXPORT =
  qw(read_hmmer3_output read_hmmer2_output print_htab hmm_database_info get_cutoffs_for_hmm_accession build_alignment read_hmmer3_output2);
our @EXPORT_OK = qw ();

sub read_hmmer3_output {
    my ($path) = @_;
    my $retval = {};

    my $in_result = 0;
    my $in_hit_scores = 0;
    my $in_domain_scores = 0;

    if ( $path ne '' ) { 
	chomp $path;
	my @statd = stat $path;
	$retval->{'info'}->{ 'search_date' } = ( ( localtime( $statd[ 9 ] ) )[ 3 ] ) . "-" 
	    . ( ( localtime( $statd[ 9 ] ) )[ 4 ] + 1 ) . "-" 
	    . ( ( localtime( $statd[ 9 ] ) )[ 5 ] + 1900 );
    } else {
	$retval->{'info'}->{ 'search_date' } = ( ( localtime )[ 3 ] ) . "-" 
	    . ( ( localtime )[ 4 ] + 1 ) . "-" 
	    . ( ( localtime )[ 5 ] + 1900 );
    }   

    open(my $fh, "< $path") or die("Unable to open $path: $!");
    while( my $line = <$fh> ) {
	next if( $line =~ /^\s*$/ || $line =~ /inclusion threshold/ );

	# The program line
	if ( $line =~ /^\#\s*((hmmscan)\.*)/ ) {
	    $retval->{'info'}->{'program'} = $1;

	 # The version
	 } elsif( $line =~ /^\#\s*HMMER ([\d\.]+)\s+\(([^\)]+)\)/ ) {
	     ( $retval->{'info'}->{'version'}, $retval->{'info'}->{'release'} ) = ( $1, $2 );

	 # The hmm database searched
	 } elsif( $line =~ /^\#\s+target HMM database:\s+(\S+)/ ) {
	     $retval->{'info'}->{'hmm_file'} = $1;

	 # The query file
	 } elsif( $line =~ /^\#\squery sequence (file|database):\s+(.+)/ ) {
	     $retval->{'info'}->{'sequence_file'} = $2;

	 # This indicates we're parsing a hit.
	 } elsif( $line =~ /^Query\:\s+(\S+)/ ) {
	     my $data = &_parse_hmmpfam3_hit( $fh );
	     $retval->{'queries'}->{$1} = $data;

	 }
    }

    close($fh);
    return $retval;    
}

sub _parse_hmmpfam3_hit {
    my ($fh) = @_;
    my $data = {};
    
    my $in_hit_scores = 0;
    my $in_domain_scores = 0;
    my $hit_acc;
    
    while( my $line = <$fh> ) {
	chomp($line);
	last if( $line =~ m|//| );
	next if( $line =~ /^\s*$/ );

	if( $line =~ /^Domain annotation/ ) {
	    $in_hit_scores = 0;

	} elsif( $in_hit_scores && $line !~ /^\s+--/ && $line !~ /inclusion_threshold/ ) {
	    if( $line =~ /No hits detected/ ) {
		$data->{'hits'} = {};
		last;
	    }
	    
	    my @c = split(/\s+/, $line);
	    my $t_hit_acc = $c[9];
	    $data->{'hits'}->{$t_hit_acc}->{total_evalue} = $c[1];
	    $data->{'hits'}->{$t_hit_acc}->{total_score} = $c[2];
	    $data->{'hits'}->{$t_hit_acc}->{accession} = $c[9];
	    $data->{'hits'}->{$t_hit_acc}->{hit_description} = join(" ",@c[10..(@c-1)]);
	    $data->{'hits'}->{$t_hit_acc}->{domain_count} = $c[8];
	    $data->{'hits'}->{$t_hit_acc}->{frame} = "";
	} elsif( $line =~ /^>>\s*(\S+)/ ) {
	    $in_domain_scores = 1;
	    $hit_acc = $1;
	    
	} elsif( $in_domain_scores && $line =~ /Alignments for each/ ) {
	    undef $hit_acc;
	    $in_domain_scores = 0;

	} elsif( $in_domain_scores && $line !~ /^\s*[>\#-]/ ) {
	    die("Didn't parse hit accession from header line before getting to domain table")
		unless( $hit_acc );
	    die("The hit accession [$hit_acc] didn't exist in lookup") unless( exists( $data->{'hits'}->{$hit_acc} ) );

	    if( $line =~ /\[No individual domains that/ ) {
		$data->{'hits'}->{$hit_acc}->{'domains'} = {};
		$in_domain_scores = 0;
		next;
	    }

	    my @c = split(/\s+/, $line);

	    if( $c[1] eq 'targets' ||  $c[1] eq 'reported' || $c[1] eq 'Fwd') {
		print Dumper( $data->{'hits'}->{$hit_acc} );
		print "LINE: $line\n";
		print Dumper( @c );
		die("Issue parsing");
	    }
	    
	    my $a = $c[1];
	    my $b = $a + 0;
	    if( $b ne $a ) {
		print Dumper( $data->{'hits'}->{$hit_acc} );
		print "LINE: $line\n";
		print Dumper( \@c );
		die("c[1] not numeric");
	    }
		
	    
	    $data->{'hits'}->{$hit_acc}->{'domains'}->{$c[1]}->{'seq_f'} = $c[10];
	    $data->{'hits'}->{$hit_acc}->{'domains'}->{$c[1]}->{'seq_t'} = $c[11];
	    $data->{'hits'}->{$hit_acc}->{'domains'}->{$c[1]}->{'hmm_f'} = $c[7];
	    $data->{'hits'}->{$hit_acc}->{'domains'}->{$c[1]}->{'hmm_t'} = $c[8];
	    $data->{'hits'}->{$hit_acc}->{'domains'}->{$c[1]}->{'domain_score'} = $c[3];
	    $data->{'hits'}->{$hit_acc}->{'domains'}->{$c[1]}->{'domain_evalue'} = $c[6];
	} elsif( $line =~ /^\s+E-value\s+score/ ) {
	    $in_hit_scores = 1;
	}
    }
    
    return $data;
}

## Subroutine to parse hmmscan (Hmmer3.0) output
sub read_hmmer3_output_old {
    my $path = shift;
    my $data = {};
    my @lines;
    if ( $path ne '' ) { 
	chomp $path;
	my @statd = stat $path;
	$data->{ 'search_date' } = ( ( localtime( $statd[ 9 ] ) )[ 3 ] ) . "-" 
	    . ( ( localtime( $statd[ 9 ] ) )[ 4 ] + 1 ) . "-" 
	    . ( ( localtime( $statd[ 9 ] ) )[ 5 ] + 1900 );
	open( FH, "$path" ) || die "Can't open $path for reading: $!\n";
	chomp( @lines = <FH> );
	close FH; 
    } else {
	chomp( @lines = <STDIN> );
	$data->{ 'search_date' } = ( ( localtime )[ 3 ] ) . "-" 
	    . ( ( localtime )[ 4 ] + 1 ) . "-" 
	    . ( ( localtime )[ 5 ] + 1900 );
    }   
    if ( !@lines ) { 
	carp "No data read from input $path";
	return undef;
    }
    my $i = 0;
    while ($i < @lines) {
	if ( $lines[ $i ] =~ /^#\s*((hmmscan)\.*)/ ) {
	     $data->{ 'program' } = $1;
	     my $version_line = $lines[ ++$i ];
	     $version_line =~ /^#\s*HMMER (\d+)\.(\S+)/;
	     ( $data->{ 'version' }, $data->{ 'release' } ) = ( $1, $2 );
	     $i +=2;
	     last;
	 }
	++$i;
    }
    until ( $lines[ $i ] =~ /^\s*$/ ) {
	if ( $lines[ $i ] =~ /^#\s*target HMM database:\s+(\S+)/ ) {
	     $data->{ 'hmm_file' } = $1;
	 } elsif ( $lines[ $i ] =~ /^#\s*query sequence (file|database):\s+(.+)/ ) {
		   $data->{ 'sequence_file' } = $2;
	       }
	$i++;
	die "Failure to parse" if $i > @lines;
    }
    $i++;
    if ( $lines[ $i ] =~ /^Query:\s+(\S+)/ ) {
	$data->{ 'query' } = $1;
	print "Found query: $1\n";
	$i++;
    }   	
    until ( $lines[$i] =~ /^\s*$/ ) {
	if ( $lines[ $i ] =~ /^Scores for/ ) {
	    $i++; # this skips the separator row
	    my $headers = $lines[ ++$i ];
	    $i++;        # this skips the separator row
	} elsif ( $lines[ $i ] !~ /No hits detected that satisfy reporting thresholds/i ) {
	    $i++ if ($lines[$i] =~ /inclusion threshold/g);
	    my @c = split /\s+/, $lines[ $i ],11;
	    my $hit_index = $c[9];
	    $data->{hit}{$hit_index}{total_evalue} = $c[1];
	    $data->{hit}{$hit_index}{total_score} = $c[2];
	    $data->{hit}{$hit_index}{accession} = $c[9];
	    $data->{hit}{$hit_index}{hit_description} = $c[10];
	    $data->{hit}{$hit_index}{domain_count} = $c[8];
	    $data->{hit}{$hit_index}{frame} = "";
	} else {
	    return $data;
	}
	$i++;
	die "Failure to parse" if $i > @lines;
    }
    $i++;
    if ( $lines[ $i ] =~ /^Domain annotation for each model/ ) { 
	$i++;
    }
    until ($lines[$i] =~ /^\/\/$/) {
	if ( $lines[ $i ] !~ /No targets detected that satisfy reporting thresholds/ ) { 
	    if($lines[$i] =~ />>/) {
		my @c = split /\s+/, $lines[ $i ]; 
		my $hit_index = $c[ 1 ];
		$i += 3;

		if ( !defined $data->{ 'hit' }->{ $hit_index } ) { 
		    warn "Why doesn't '$hit_index' match an existing identifier?";
		} else {
		    until ($lines[$i] =~ /^\s*$/) {
			my @res = split /\s+/, $lines[ $i ];
			$data->{hit}{$hit_index}{domain}{$res[1]}{seq_f} = $res[10];
			$data->{hit}{$hit_index}{domain}{$res[1]}{seq_t} = $res[11];
			$data->{hit}{$hit_index}{domain}{$res[1]}{hmm_f} = $res[7];
			$data->{hit}{$hit_index}{domain}{$res[1]}{hmm_t} = $res[8];
			$data->{hit}{$hit_index}{domain}{$res[1]}{domain_score} = $res[3];
			$data->{hit}{$hit_index}{domain}{$res[1]}{domain_evalue} = $res[6];
			$i++;
		    }
		}
	    }
	} else {
	    return $data;
	}
	$i++;
    }
    return $data;   
}
## Use with caution. Or just don't use this method.
sub read_hmmer2_output {
	my $path = shift;
	my $data = {};
# die "Bad Parser under reconstruction.\n\n";
# drink in data
	my @lines;

# drink in the output from file or stdin
	if ( $path ne '' ) {
		chomp $path;
		my @statd = stat $path;
		$data->{ 'search_date' } =
			( ( localtime( $statd[ 9 ] ) )[ 3 ] ) . "-"
			. ( ( localtime( $statd[ 9 ] ) )[ 4 ] + 1 ) . "-"
			. ( ( localtime( $statd[ 9 ] ) )[ 5 ] + 1900 );
		open( FH, "$path" )
			|| die "Can't open $path for reading: $!\n";
		chomp( @lines = <FH> );
		close FH;
	}
	else {
		chomp( @lines = <STDIN> );
		$data->{ 'search_date' } =
			( ( localtime )[ 3 ] ) . "-"
			. ( ( localtime )[ 4 ] + 1 ) . "-"
			. ( ( localtime )[ 5 ] + 1900 );
	}
	if ( !@lines ) {
		carp "No data read from input $path";
		return undef;
	}
	my $i = 0;

# first line grouping is company, package and license info
# warn "Parsing License. Line $i\n";
# amahurkar:1/15/08 Seems like the current output does not have licesning info, so commenting this
# so commenting out these lines 
#until ( $lines[ $i ] eq "" ) 
#    $data->{ 'header' } .= $lines[ $i ] . "\n";
#    $i++;

#$i++;

# next group is program and version
# warn "Parsing Program and version. Line $i\n";
# amahurkar:1/15/08 the format has changed and now there is no blank space
# after program name, so we are using '- -' as tha match param
# to stop parsing for program name
#until ( $lines[ $i ] eq "" ) 
#until ( $lines[ $i ] =~ m/^-\s-/)  # It doesn't work with LDhmmpfam v1.5.4
#until ( $lines[ $i ] =~ m/^-{8}\s+-/) 
	while ($i < @lines) {
		if ( $lines[ $i ] =~ /^((hmmpfam|hmmsearch)\.*)/ ) {
			$data->{ 'program' } = $1;
			my $version_line = $lines[ ++$i ];
			$version_line =~ /^HMMER (\d+)\.(\S+)/;
			( $data->{ 'version' }, $data->{ 'release' } ) = ( $1, $2 );
			$i +=2;
			last;
		}
		++$i;
	}

# next group is program parameters
# warn "Parsing Parameters. Line $i\n";
# amahurkar:1/15/08 the format has changed and now there is no blank space
# after program name, so we are using '- -' as tha match param
# to stop parsing for program parameters
	until ( $lines[ $i ] =~ /^\s*$/ ) {
		if ( $lines[ $i ] =~ /^HMM file:\s+(\S+)/ ) {
			$data->{ 'hmm_file' } = $1;
		} elsif ( $lines[ $i ] =~ /^Sequence (file|database):\s+(.+)/ ) {
			$data->{ 'sequence_file' } = $2;
		} elsif ( $lines[ $i ] =~ /^per-sequence score cutoff:\s+(.+)/ ) {
			$data->{ 'total_score_cutoff' } = $1;
		} elsif ( $lines[ $i ] =~ /^per-domain score cutoff:\s+(.+)/ ) {
			$data->{ 'domain_score_cutoff' } = $1;
		} elsif ( $lines[ $i ] =~ /^per-sequence E-value cutoff:\s+(.+)/ ) {
			$data->{ 'total_evalue_cutoff' } = $1;
		} elsif ( $lines[ $i ] =~ /^per-domain E-value cutoff:\s+(.+)/ ) {
			$data->{ 'domain_evalue_cutoff' } = $1;
		}
		$i++;
		die "Failure to parse" if $i > @lines;
	}
	$i++;

# get query info
# warn "Parsing Query Info. Line $i\n";
#until ( $lines[ $i ] eq "" ) 
	until ( $lines[$i] =~ /^\s*$/ ) {
		if ( $lines[ $i ] =~ /^Query (?:HMM|sequence):\s+(.+)/ ) {
			$data->{ 'query' } = $1;
		}
		elsif ( $lines[ $i ] =~ /^Accession:\s+(.+)/ ) {
			$data->{ 'query_accession' } = $1;
		}
		elsif ( $lines[ $i ] =~ /^Description:\s+(.+)/ ) {
			$data->{ 'query_description' } = $1;
		}
		$i++;
		die "Failure to parse" if $i > @lines;
	}
	$i++;

# next section is global search results
# warn "Parsing Global Search Results. Line $i\n";
	my $find_frame = 0;  # is datbase nucleotide sequence?
	my $hit_index;
	until ( $lines[$i] =~ /^\s*$/ ) {
		if ( $lines[ $i ] =~ /^Scores for/ ) {
			my $headers = $lines[ ++$i ];
			if ( $headers =~ /\bFr\b/ ) {
				$data->{ 'program' } .= "-frames";
				$find_frame = 1;
			}
			$i++;        # this skips the separator row
		} elsif ( $lines[ $i ] !~ /no hits above thresholds/i ) {
			my @c = split /\s+/, $lines[ $i ];
			if ( $find_frame ) {
				$hit_index = $c[ 0 ] . $c[ $#c ];
				$data->{ 'hit' }->{ $hit_index }->{ 'frame' } = pop @c;
			}
			else {
# $hit_index = $c[ 0 ]; # AP 20090807
				($hit_index = $c[0]) =~ s/\.\d+$//;
			}
			$data->{hit}{$hit_index}{accession}       = shift @c;
			$data->{hit}{$hit_index}{domain_count}    = pop @c;
			$data->{hit}{$hit_index}{total_evalue}    = pop @c;
			$data->{hit}{$hit_index}{total_score}     = pop @c;
			$data->{hit}{$hit_index}{hit_description} = join " ", @c;
		}
		else {
			return $data;
		}
		$i++;
		die "Failure to parse" if $i > @lines;
	}
	$i++;

# next section is domain breakdown
# warn "Parsing Domain Breakdown. Line $i\n";
#until ( $lines[ $i ] eq "" ) 
	until ($lines[$i] =~ /^\s*$/) {
		if ( $lines[ $i ] =~ /^Parsed for domains/ ) {
			$i += 2;  # to skip header and separator
		} elsif ( $lines[ $i ] !~ /no hits above thresholds/ ) {

# $lines[$i] =~ /^(\w+)\s+(\d+)\/\d+\s+(\d+)\s+(\d+)\s[\[\.\]]{2}\s+(\d+)\s+(\d+)\s[\[\.\]]{2}\s+(-?[\.\d]+)\s+([\.\-e\d]+)\s*(\-?\d)?/;
			my @c = split /\s+/, $lines[ $i ];
			if ( $find_frame ) {
				$hit_index = $c[ 0 ] . $c[ $#c ];
			}
			else {
#$hit_index = $c[ 0 ]; # AP 20090807
				($hit_index = $c[0]) =~ s/\.\d+$//;
			}
			if ( !defined $data->{ 'hit' }->{ $hit_index } ) {
				warn
					"Why doesn't '$hit_index' match an existing identifier?";
			}

# $data->{'hit'}->{$hit_index}->{'domain'}->{$2}->{'seq_f'} = $3;
			my ( $d, $t ) = split /\//, $c[ 1 ];
			$data->{hit}{$hit_index}{domain}{$d}{seq_f} = $c[2];
			$data->{hit}{$hit_index}{domain}{$d}{seq_t} = $c[3];
			$data->{hit}{$hit_index}{domain}{$d}{hmm_f} = $c[5];
			$data->{hit}{$hit_index}{domain}{$d}{hmm_t} = $c[6];
			$data->{hit}{$hit_index}{domain}{$d}{domain_score} = $c[8];
			$data->{hit}{$hit_index}{domain}{$d}{domain_evalue} = $c[9];
		}
		$i++;
		die "Failure to parse" if $i > @lines;
	}
	$i++;
	if ( $data->{ 'program' } =~ /hmmsearch/ ) {

# next section is alignments
# warn "Parsing Alignments. Line $i\n";
	    if ( $lines[ $i ] =~ /^Alignments of top-scoring domains/ ) {
		$i++;
		## kgalens [11/13/2012]
		## Not sure how this ever worked? $domain isn't in scope here.
		## I've added use strict; and this module wouldn't. So adding the next line
		## just so it will run.
		my $domain;
		my $hit;
		my $hit_index;
		## END kgalens [11/13/2012]

		until ( $lines[$i] =~ /^\s*$/ ) {
		    if ( $lines[ $i ] =~ /^(\S+): domain (\d+)/ ) {
			($hit_index, $hit, $domain ) = ( $1, $1, $2 );
			$hit =~ s/(.{10}).*/$1/;
# warn "  Parsing hit $hit_index. Line $i\n";
			if ( $find_frame ) {
			    if (   $lines[ $i ] =~ /Fr = ([\-\d]+)/ || $lines[ $i ] =~ /\. frame ([\-\d]+)/ ) {
				$hit_index .= $1;
			    } else {
				warn "ERROR: Couldn't find frame from:\n '$lines[$i]'";
			    }
			}
			if ( !defined $data->{ 'hit' }->{ $hit_index } ) {
			    warn "Why doesn't '$hit_index' match an existing identifier?";
			}		
			$i++;
		    }
		    if ( $lines[ $i ] =~ /\bRF\b/ ) { ## <rar> WHY??
			$i++;
		    }

# capture aligned hmm consensus
		    my $hmm_seq = $lines[ $i ];
		    $hmm_seq =~ s/\s+//g;
		    $hmm_seq =~ s/[\*\-\>\<]//g;

		    $data->{ 'hit' }->{ $hit_index }->{ 'domain' }->{ $domain }->{hmm_seq} .= $hmm_seq;
		    until ( $lines[ $i ] =~ /^\s+\Q$hit\E/) { ## changed from /\b\Q$hit\E\b/ ) 
			$i++;
		    }
		    my $prot_seq = $lines[ $i ];
		    if ( $prot_seq =~ /\w+\s+(\d+|\-)\s+(\S+)\s+(\d+|\-)/ ) {
			$data->{ 'hit' }->{ $hit_index }->{ 'domain' }->{ $domain }->{prot_seq} .= $2;
		    }
		    $i += 2;  # skip the blank line and move on to the next
		    die "Failure to parse" if $i > @lines;
		}	
	    }
	    $i++;

# next (last) section is statistics
# warn "Parsing Statistics. Line $i\n";
	    my @data;
	    while ( $i < @data ) {
		if ( $lines[ $i ] =~ /^\s+mu =\s+(-?\d+)/ ) {
		    $data->{ 'mu' } = $1;
		} elsif ( $lines[ $i ] =~ /^\s+lambda =\s(-?\d+)/ ) {
		    $data->{ 'lambda' } = $1;
		} elsif ( $lines[ $i ] =~ /chi-sq statistic =\s(\d+)/ ) {
		    $data->{ 'chisq' } = $1;
		} elsif ( $lines[ $i ] =~ /Total sequences searched:\s*(\d+)/ ) {
		    $data->{ 'tot_seq_searched' } = $1;
		} elsif ( $lines[ $i ] =~ /Whole sequence top hits/ ) {
		    $lines[ ++$i ] =~ /(\d+)/;
		    $data->{ 'total_hits' } = $1;
		    $lines[ ++$i ] =~ /(\d+)/;
		    $data->{ 'total_hits_above_evalue_cutoff' } = $1;
		} elsif ( $lines[ $i ] =~ /Domain top hits/ ) {
		    $lines[ ++$i ] =~ /(\d+)/;
		    $data->{ 'domain_hits' } = $1;
		    $lines[ ++$i ] =~ /(\d+)/;
		    $data->{ 'domain_hits_above_evalue_cutoff' } = $1;
		}
		$i++;
		die "Failure to parse" if $i > @lines;
	    }
	}
	return $data;
}

sub hmm_database_info {
    my $dbh   = shift;
    my $hmm_q =
      "SELECT hmm_acc, hmm_len, trusted_cutoff, noise_cutoff, hmm_com_name,"
      . " trusted_cutoff2, noise_cutoff2, gathering_cutoff, gathering_cutoff2"
      . " FROM hmm2"
      . " WHERE is_current = 1";
    my $HMM = $dbh->selectall_hashref( $hmm_q, 'hmm_acc' );
    return $HMM;
}

sub print_htab {
#NOTE:  This will produce results if the data hash was created through the 'read_hmmer3_output' subroutine
#		This is because of a change in the naming of some of the property keys (ex 'hit' to 'hits')
    my $data   = shift;
    my $HMM    = shift;
    my $output = shift;
    foreach my $qry_id ( keys %{$data->{'queries'}} ) {
    	foreach my $hit (
        	sort {
                    $data->{'queries'}->{$qry_id}->{'hits'}->{ $b }->{ 'total_score' } <=> $data->{'queries'}->{$qry_id}->
                    	{'hits'}->{ $a }->{ 'total_score' }
        	} keys %{ $data->{'queries'}->{$qry_id}->{'hits'} } )
	{
            my $h = $data->{'queries'}->{$qry_id}->{'hits'}->{ $hit };
            next if (scalar keys %{$h->{'domains'}} == 0);	#skip model hits that have no domain hits
            foreach my $domain ( sort { $a <=> $b } keys %{ $h->{ 'domains' } } )
            {
            	# for convenience
            	my $dh = $h->{ 'domains' }->{ $domain };
            	if ( $data->{'info'}->{ 'program' } =~ /hmmsearch/ ) {	#hmmer2 is currently deprecated so this will probably error
                	my $hmm_com_name =
                	    $HMM->{ $data->{ 'query' } }->{ 'hmm_com_name' }
               	   		? $HMM->{ $data->{ 'query' } }->{ 'hmm_com_name' }
                  		: $data->{ 'query_description' };
                	print $output "$data->{query}"
                  	    . "\t$data->{search_date}"
              		    . "\t$HMM->{$data->{query}}->{hmm_len}"
                	    . "\t$data->{program}"
                 	    . "\t$data->{sequence_file}"
             		    . "\t$h->{accession}"
             		    . "\t$dh->{hmm_f}"
             		    . "\t$dh->{hmm_t}"
              		    . "\t$dh->{seq_f}"
             		    . "\t$dh->{seq_t}"
                	    . "\t$h->{frame}"
               		    . "\t$dh->{domain_score}"
               		    . "\t$h->{total_score}"
                	    . "\t$domain"
                	    . "\t$h->{domain_count}"
               		    . "\t$hmm_com_name"
               		    . "\t$h->{hit_description}"
                	    . "\t$HMM->{$data->{query}}->{trusted_cutoff}"
                	    . "\t$HMM->{$data->{query}}->{noise_cutoff}"
                	    . "\t$h->{total_evalue}"
                	    . "\t$dh->{domain_evalue}"
                	    . "\t$HMM->{$data->{query}}->{trusted_cutoff2}"
                	    . "\t$HMM->{$data->{query}}->{noise_cutoff2}"
                	    . "\t$HMM->{$data->{query}}->{gathering_cutoff}"
                	    . "\t$HMM->{$data->{query}}->{gathering_cutoff2}" . "\n";
            	}
            	elsif ( $data->{'info'}->{ 'program' } =~ /hmmscan|hmmpfam/ ) {
                    my $hmm_com_name =
                	    $HMM->{ $hit }->{ 'hmm_com_name' }
                	  ? $HMM->{ $hit }->{ 'hmm_com_name' }
                	  : $h->{ 'hit_description' };
                    print $output "$h->{accession}"
                	  . "\t$data->{'info'}->{search_date}"
                	  . "\t$HMM->{$hit}->{hmm_len}"
                	  . "\t$data->{'info'}->{program}"
                	  . "\t$data->{'info'}->{hmm_file}"
                 	  . "\t$qry_id"
                	  . "\t$dh->{hmm_f}"
                 	  . "\t$dh->{hmm_t}"
                  	  . "\t$dh->{seq_f}"
                  	  . "\t$dh->{seq_t}"
               		  . "\t$h->{frame}"
                	  . "\t$dh->{domain_score}"
                	  . "\t$h->{total_score}"
                	  . "\t$domain"
                	  . "\t$h->{domain_count}"
                	  . "\t$hmm_com_name"
                	  . "\t$h->{hit_description}"
                	  . "\t$HMM->{$h->{accession}}->{trusted_cutoff}"
                	  . "\t$HMM->{$h->{accession}}->{noise_cutoff}"
                	  . "\t$h->{total_evalue}"
                	  . "\t$dh->{domain_evalue}"
                	  . "\t$HMM->{$h->{accession}}->{trusted_cutoff2}"
                	  . "\t$HMM->{$h->{accession}}->{noise_cutoff2}"
                	  . "\t$HMM->{$h->{accession}}->{gathering_cutoff}"
                	  . "\t$HMM->{$h->{accession}}->{gathering_cutoff2}" . "\n";
            	}
            }
        }
    }
}

sub build_alignment {
    my $data         = shift;
    my $instructions = shift;

    # build output file name
    my $output_file;
    $output_file =
      $instructions->{file_prefix} . "." . $instructions->{file_format};
    open my $OUT, ">$output_file"
      or croak "Can't open '$output_file' as output file: $!\n";
    select $OUT;

    # retrieve aligned sequences
    my %screened;
    foreach my $hit ( keys %{ $data->{ 'hit' } } ) {

        # screen for total score cutoffs
        if ( $data->{hit}->{ $hit }->{total_score} >=
               $instructions->{total_bit_cutoff}
            && $data->{hit}->{ $hit }->{total_evalue} <=
            $instructions->{total_evalue_cutoff} )
        {
            foreach my $domain ( keys %{ $data->{hit}->{ $hit }->{domain} } )
            {
                if ( $data->{hit}->{ $hit }->{domain}->{ $domain }
                    ->{domain_score} >= $instructions->{domain_bit_cutoff}
                    && $data->{hit}->{ $hit }->{domain}->{ $domain }
                    ->{domain_evalue} <=
                    $instructions->{domain_evalue_cutoff} )
                {
                    $screened{ $hit } = $domain;
                }
            }
        }
    }

    # Now that we have sequences aligned to hmm sequence, we have to translate
    # this into a multiple alignment. Assign each position in each alignment to
    # a position on the hmm 'sequence', and keep track of gaps in the hmm alignment
    my %DIST;
    my $ref;
    foreach my $hit ( keys %screened ) {
        $ref =
          $data->{ 'hit' }->{ $hit }->{ 'domain' }->{ $screened{ $hit } };

        # split aligned hmm seq and aligned protein seq into arrays
        my @hmma  = split / */, $ref->{hmm_seq};
        my @prota = split / */, $ref->{prot_seq};

        # these should be the same length. If not, there's an error.
        if ( @hmma != @prota ) {
            croak "Length of hmm alignment (" . @hmma . ")"
              . " is not equal to protein alignment ("
              . @prota . ")"
              . ": $data->{query}/$data->{hit}->{$hit}->{accession}\n"
              . "$ref->{hmm_seq}\n@hmma\n$ref->{prot_seq}\n@prota\n";
        }

       # assign each position in the protein alignment to its hmm alignment position,
       # if one exists
        my $hmm_pos = $ref->{hmm_f};
        my $gap     = 0;
        for ( my $i = 0 ; $i < @hmma ; $i++ ) {
            if ( $hmma[ $i ] ne "." ) {

                # assign position in the protein alignment to its hmm alignment position,
                $prota[ $i ] = $hmm_pos;

                # record max gap distance between hmm alignment positions.
                $DIST{ $hmm_pos } = $gap
                  if ( $gap >= $DIST{ $hmm_pos } );
                $gap = 0;
                $hmm_pos++;
            }
            else {
                $gap++;
            }
        }
        $ref->{aln_map} = \@prota;
    }

    # Now go back through (now that we've fully expanded the hmm alignment
    # to include any and all gaps) and make aligned protein sequence
    foreach my $hit ( keys %screened ) {
        $ref =
          $data->{ 'hit' }->{ $hit }->{ 'domain' }->{ $screened{ $hit } };
        my @prot_seq = split / */, $ref->{prot_seq};

        # start our aligned protein with any gap resulting from a partial HMM hit.
        my $aln_prot = "." x ( $ref->{hmm_f} - 1 );
        my $insert = 0;
        for ( my $i = 0 ; $i < @prot_seq ; $i++ ) {

            # grab the hmm alignment position for each protein alignment position
            my $pos = $ref->{aln_map}->[ $i ];
            if ( $pos =~ /\d+/ ) {

                # if it maps to a position, first insert any gap from the hmm alignment
                $aln_prot .= "." x ( $DIST{ $pos } - $insert );

                # then add the aa.
                $aln_prot .= $prot_seq[ $i ];
                $insert = 0;
            }

            # if it is an insertion (ie the hmm alignment shows gap), insert a gap
            else {
                $aln_prot .= $prot_seq[ $i ];
                $insert++;
            }
        }
        $aln_prot =~ s/\-/\./g;
        $ref->{aln_prot} = $aln_prot;
    }

    # Now print out in selected format
    # Stockholm format
    if ( $instructions->{file_format} eq "mul" ) {
        print "# STOCKHOLM 1.0\n";
        foreach my $hit (
            sort {
                $data->{ 'hit' }->{ $b }
                  ->{ 'total_score' } <=> $data->{ 'hit' }->{ $a }
                  ->{ 'total_score' }
            } keys %screened
          )
        {
            my $domain     = $screened{ $hit };
            my $hit_ref    = $data->{hit}->{ $hit };
            my $domain_ref = $hit_ref->{domain}->{ $domain };

            # each line should look like:
            # prot_acc/coord-coord sequence
            printf "%-40s%s\n",
              (
                "$hit_ref->{accession}/$domain_ref->{seq_f}-$domain_ref->{seq_t}",
                $domain_ref->{aln_prot}
              );
        }
    }

    # FASTA format
    elsif ($instructions->{file_format} eq "fasta"
        || $instructions->{file_format} eq "fa" )
    {
        use TIGR::FASTArecord;
        foreach my $hit (
            sort {
                $data->{ 'hit' }->{ $b }
                  ->{ 'total_score' } <=> $data->{ 'hit' }->{ $a }
                  ->{ 'total_score' }
            } keys %screened
          )
        {
            my $domain     = $screened{ $hit };
            my $hit_ref    = $data->{hit}->{ $hit };
            my $domain_ref = $hit_ref->{domain}->{ $domain };
            my $aln_prot = $domain_ref->{aln_prot};
            $aln_prot =~ s/(.{60})/$1\n/g;
            $aln_prot =~ s/\n+/\n/g;
            chomp $aln_prot;
            my $header =
              ">$hit_ref->{accession}/$domain_ref->{seq_f}-$domain_ref->{seq_t}\n";
            print $header . $aln_prot . "\n";
        }
    }

    # MSF format
    elsif ( $instructions->{file_format} eq "msf" ) {
        my $head_len = 40;
        my ( %new_acc, %tmp_seq );
        my $header_line = 0;
        my @alignment;
        foreach my $hit (
            sort {
                $data->{ 'hit' }->{ $b }
                  ->{ 'total_score' } <=> $data->{ 'hit' }->{ $a }
                  ->{ 'total_score' }
            } keys %screened
          )
        {
            my $domain     = $screened{ $hit };
            my $hit_ref    = $data->{hit}->{ $hit };
            my $domain_ref = $hit_ref->{domain}->{ $domain };
            my $len        = length( $domain_ref->{aln_prot} );
            if ( $header_line == 0 ) {

                #added DUMMY checksum for HMMER 2.2g hmmbuild. DHH.
                print
                  "PileUp\n\n   MSF:  $len  Type: P  Check: 1111  ..\n\n";
                $header_line = 1;
            }

            # print accession list at top
            printf " Name: %-40s Len:  $len  Check:   0  Weight:  1.0\n",
              "$hit_ref->{accession}/$domain_ref->{seq_f}-$domain_ref->{seq_t}";

            # prepare alignment for bottom
            my @tmp_pep = split //, $domain_ref->{aln_prot};
            for ( my $i = 0 ; $i < ( $len / 50 ) ; $i++ ) {
                $alignment[ $i ] .= sprintf "%-40s",
                  "$hit_ref->{accession}/$domain_ref->{seq_f}-$domain_ref->{seq_t}";
                for ( my $b = 0 ; $b < 5 ; $b++ ) {
                    for ( my $a = 0 ; $a < 10 ; $a++ ) {
                        $alignment[ $i ] .=
                          $tmp_pep[ $a + ( $b * 10 ) + ( 50 * $i ) ];
                    }
                    $alignment[ $i ] .= " ";
                }
                $alignment[ $i ] .= "\n";
            }
        }
        print "\n//\n\n\n";
        foreach my $block ( @alignment ) {
            print "$block\n\n";
        }
    }
    else {
        croak
          "Don't recognize alignment file format $instructions->{file_format}\n";
    }
    select STDOUT;
}

sub get_cutoffs_for_hmm_accession {
    my $dbh       = shift;
    my $accession = shift;
    my $hmm_q     =
      "select trusted_cutoff, trusted_cutoff2, noise_cutoff from egad..hmm2 where hmm_acc = '$accession'";
    return $dbh->selectrow_hashref( $hmm_q );
}
1;
