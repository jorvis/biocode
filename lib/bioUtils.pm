
=head1 DESCRIPTION

Too much code is written over and over.  My own useful utilities are
tucked away here.  Take them as you want, edit your copies, feel free
to let me know if you have generic useful things I'm missing.

Joshua Orvis
jorvis@gmail.com

=cut


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





















## speak the truth, even if your voice shakes.
1==1;
