
=head1 DESCRIPTION

Too much code is written over and over.  My own useful utilities are
tucked away here.  Take them as you want, edit your copies, feel free
to let me know if you have generic useful things I'm missing.

Joshua Orvis
jorvis@gmail.com

=cut


use strict;


my  %genetic_code = (
    # Alanine
    'GCA' => 'A', 'GCC' => 'A', 'GCG' => 'A', 'GCT' => 'A',
    # Arginine
    'CGA' => 'R', 'CGC' => 'R', 'CGG' => 'R', 'CGT' => 'R', 'AGA' => 'R', 'AGG' => 'R',
    # Asparagine
    'AAC' => 'N', 'AAT' => 'N',
    # Aspartic Acid
    'GAC' => 'D', 'GAT' => 'D', 
    # Cysteine
    'TGC' => 'C', 'TGT' => 'C',
    # Glutamic acid
    'GAA' => 'E', 'GAG' => 'E',
    # Glutamine
    'CAA' => 'Q', 'CAG' => 'Q',
    # Glycine
    'GGA' => 'G', 'GGC' => 'G', 'GGG' => 'G', 'GGT' => 'G',
    # Histidine
    'CAC' => 'H', 'CAT' => 'H',
    # Isoleucine
    'ATA' => 'I', 'ATC' => 'I', 'ATT' => 'I',
    # Leucine
    'TTA' => 'L', 'TTG' => 'L', 'CTA' => 'L', 'CTC' => 'L', 'CTG' => 'L','CTT' => 'L',
    # Lysine
    'AAA' => 'K', 'AAG' => 'K',
    # Methionine
    'ATG' => 'M',
    # Phenylalanine
    'TTC' => 'F', 'TTT' => 'F',
    # Proline
    'CCA' => 'P', 'CCC' => 'P', 'CCG' => 'P', 'CCT' => 'P',
    # Serine
    'TCA' => 'S', 'TCC' => 'S', 'TCG' => 'S', 'TCT' => 'S', 'AGC' => 'S', 'AGT' => 'S',
    # Threonine
    'ACA' => 'T', 'ACC' => 'T', 'ACG' => 'T', 'ACT' => 'T',
    # Tryptophan
    'TGG' => 'W',
    # Tyrosine
    'TAC' => 'Y', 'TAT' => 'Y',
    # Valine
    'GTA' => 'V', 'GTC' => 'V', 'GTG' => 'V', 'GTT' => 'V',
    # STOP
    'TAA' => '*', 'TAG' => '*', 'TGA' => '*',
);


sub read_list_file {
    my ($file, $regex) = @_;
    
    my $files = [];
    
    open(my $ifh, $file) || die "failed to read input list file ($file): $!";
    
    while (my $line = <$ifh>) {
        chomp $line;
        next unless $line =~ /\S/;
        
        if ( defined $regex ) {
            if ( $line =~ /$regex/ ) {
                push @$files, $line;
            }
        } else {
            push @$files, $line;
        }
    }
    
    return $files;
}




sub characters_per_line {
    my ($string, $char_count) = @_;
    
    $string =~ s/\s//g;
    
    ## set a default if needed
    if (! defined $char_count ) {
        $char_count = 60;   ## accepted FASTA default
    }
    
    my @new_string_parts  = ();
    
    while ( $string =~ /(.{1,$char_count})/g ) {
        push @new_string_parts, $1;
    }

    return join("\n", @new_string_parts);
}



sub gff3_get_column_9_value {
    my ($col9_str, $key) = @_;
    my $value = undef;

    ## have to check for versions with and without a closing semi-colon
    if ( $col9_str =~ /${key}\=(.+?)\;/ ) {
        $value = $1;
    } elsif ( $col9_str =~ /${key}\=(.+)/ ) {
        $value = $1;
    }

    return $value;
}

sub translate_sequence {
    my ($dna, $label) = @_;
    
    my $polypeptide_seq = '';
    
    for(my $i=0; $i<(length($dna) - 2); $i+=3) {
        $polypeptide_seq .= codon2aa( substr($dna,$i,3) );
    }
    
    ## report if this had an internal stop
    if ( $polypeptide_seq =~ /\*.*?.$/ ) {
        _log("ERROR: internal stops detected in feature with feature_id: $label: ($polypeptide_seq)");
    }
    
    return $polypeptide_seq;
}


## taken from Tisdall's book
sub codon2aa {
    my($codon) = @_;

    $codon = uc $codon;

    if (exists $genetic_code{$codon}) {
        return $genetic_code{$codon};
    } else {
        return 'X';
    }
}







## speak the truth, even if your voice shakes.
1==1;
