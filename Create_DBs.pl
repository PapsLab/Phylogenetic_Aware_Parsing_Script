#!/usr/bin/perl -w
#Jordi Paps, Oxford 2015
#Beginning Perl for Bioinformatics
use strict;
use AnyDBM_File;
use Fcntl;

#Define file variables
my $MCL_out = "Input/01_Genomes2014_BLAST_MCL2.out"; 
my $MCL_columns_parsed = "Input/02_Genomes2014_BLAST_MCL2_gene_numbers_parsed.out";
my $FASTA_file = "Input/00_Genomes2014.fasta";

#Check if files exists, open them
check_file ($MCL_out);
check_file ($MCL_columns_parsed);
check_file ($FASTA_file);

#Create DBM files from the MCL output file for fast lookups.
tie (my %MCL_out,"AnyDBM_File", "$MCL_out.db", O_CREAT|O_RDWR, 0777) or die + "Can't open database: $!\n";
%MCL_out = MCL_output_to_hash ($MCL_out);
print "Database $MCL_out.db created\n";

#Create DBM files from the MCL parsed file for fast lookups.
tie (my %MCL_columns_parsed,"AnyDBM_File", "$MCL_columns_parsed.db", O_CREAT|O_RDWR, 0777) or die + "Can't open database: $!\n";
%MCL_columns_parsed = MCL_columns_parsed_to_hash ($MCL_columns_parsed);
print "Database $MCL_columns_parsed.db created\n";

#Create DBM files from the FASTA file for fast lookups
my @FASTA_file = open_file ($FASTA_file);
tie (my %fasta,"AnyDBM_File", "$FASTA_file.db", O_CREAT|O_RDWR, 0666) or die + "Can't open database: $!\n";
%fasta = fasta_to_hash (@FASTA_file);
print "Database $FASTA_file.db created\n";

#Untie hash files
untie %MCL_out;
untie %MCL_columns_parsed;
untie %fasta;

#End of program
print "\nend\n\n";
exit;
################################################################################
#                                Subroutines                                   #
################################################################################

#OPEN_FILE
sub open_file {
    my ($file) = @_;
    unless (open (FILE, $file)) {print "Can't open file $file\n" && die;}
    my @file = <FILE>;
    close FILE;
    #print "\nShowing first 10 lines of file $file:\n", @file [0..9],"\n";
    return @file;
}

#CHECK_FILE: checks file properties, checks if file exists, is flat, is empty, or can be opened
#-e exists, -f flat file, -s empty
sub check_file {
    my ($filename) = @_;
    unless (-e $filename) {print "File $filename does NOT exists\n" and exit;}
    unless (-f $filename) {print "File $filename is NOT a flat file\n" and exit;}
    unless (-s $filename) {print "File $filename is empty\n" and exit;}
    unless (open (FH, $filename)) {print "File $filename can not be opened\n" and exit;}
    close FH;
    return;
}

#MCL_COLUMNS_PARSED_TO_HASH: subrout to parse a MCL output file and store it in a hash, in which the
#keys are the group homology and the values the genes found in that group
sub MCL_columns_parsed_to_hash {
	my ($filename) = @_;
	my %hash;
	my @line = '';
	my $line_counter = 0;
	
	open(FH, $filename);
	my @file = <FH>;
	foreach my $line (@file) {	#Parse file one line at a time
		$line_counter++;
		chomp $line;
		$line =~ s/^\d*\s//;	#To remove first digits, which indicate the group of homology number/ID
		#$line =~ s/\s/\t/g;
		#$line =~ s/\t$//g;
		my $key = "HG_".$line_counter;
		$hash{$key} = "$line";	
	}
	close FH;
	return %hash;
}

#MCL_TO_HASH: subrout to parse a MCL output file and store it in a hash, in which the
#keys are the group homology and the values the genes found in that group
sub MCL_output_to_hash {
	my ($filename) = @_;
	my %hash;
	my $line_counter = 0;
	
	open(FH, $filename);
	my @file = <FH>;
	foreach my $line (@file) {   			#Parse file one line at a time
		$line_counter++;
		chomp $line;
		#$line =~ s/\s/\t/g;
		#$line =~ s/\t$//g;
		my $key = "HG_".$line_counter;
		$hash{$key} = "$line";	
	}
	close FH;
	return %hash;
}

#FASTA_TO_HASH
#Takes a fasta file, and dumps seq names in an array, sequence data in another, then
#populates a hash with seq names as keys and sequence data as values.
sub fasta_to_hash {
    my (@FASTA_file) = @_;
    my $joint_sequence = '';
    my @sequence = ();
    my @seq_names = ();
    
    foreach my $line (@FASTA_file) {
        #Discard blank and comment lines
        if ($line =~ /^\s*$/ || $line =~ /^\s*#/) { next;
        #Dump fasta header in @seq_name and sequence data in @sequence (joining lines)
        } elsif ($line =~ /^>/) {
            #If you find a header, dump the data sequence from previous iteration in @sequences
            unless ($joint_sequence eq '') {
                push (@sequence, $joint_sequence);}
            #Empty $joint_sequence for next sequence round of datat dumping
            $joint_sequence = '';
            #Get header line, clean it
            $line =~ s/^>//;
            $line =~ s/\n$//; $line =~ s/\r$//;            
            push (@seq_names, $line);
        } else {
            $line =~ s/\s//g;   #remove non-sequence data (in this case, whitespace)
            $line =~ s/\n$//; $line =~ s/\r$//;
            $joint_sequence .= $line;
        }
    }       
    push (@sequence, $joint_sequence); 
    #Populate the hash with 2 arrays
    @fasta{@seq_names}= @sequence;
    return %fasta;
}