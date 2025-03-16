#!/usr/bin/env perl

use strict;
use warnings;
use 5.012;

use Cwd qw/abs_path/;
use File::Basename qw/basename/;
use File::Copy qw/move/;
use File::Fetch;
use File::Find;
use File::Temp;
use Getopt::Long;
use IO::Compress::Gzip qw/$GzipError/;
use IPC::Cmd qw/can_run run/;
use JSON qw/from_json/;

use constant ERROR => 'ERROR';
use constant WARNING => 'WARNING';
use constant INFO => 'INFO';

my $TAR = can_run('tar')
    // plog( ERROR, 'setup', "'tar' is required but not found\n" );
my $BLDBCMD = can_run('blastdbcmd') // plog(
    ERROR,
    'setup',
    "'blastdbcmd' is required but not found"
    . " (perhaps you forgot to install BLAST?)"
);

my $base_url = 'https://ftp.ncbi.nih.gov/blast/db';

my $db;
my $fo_fasta; 
my $dir_blast;
my $fo_taxmap;
my $quiet = 0;
my $gzip = 0;
my $tmpdir;

GetOptions(
    'db=s' => \$db,
    'dir_out=s' => \$dir_blast,
    'tax_map=s' => \$fo_taxmap,
    'fasta=s' => \$fo_fasta,
    'tmpdir=s' => \$tmpdir,
    'gzip' => \$gzip,
    'quiet' => \$quiet,
);

$fo_fasta = abs_path($fo_fasta)
    if (defined $fo_fasta);
$dir_blast = abs_path($dir_blast)
    if (defined $dir_blast);
$fo_taxmap = abs_path($fo_taxmap)
    if (defined $fo_taxmap);

$tmpdir //= $ENV{TMPDIR} // '/tmp';
my $meta = fetch_meta(
    $db
);
plog( ERROR, 'setup', "DB name mismatch in fetched metadata\n" )
    if ($db ne $meta->{dbname});
my $expected_size = $meta->{'bytes-total'};

my $staging = File::Temp->newdir(DIR => $tmpdir, CLEANUP => 1);
my $n_files = scalar @{ $meta->{files} };
my $n_downloaded = 0;
for my $fn (@{ $meta->{files} }) {
    ++$n_downloaded;
    plog( INFO, 'download', "Downloading $n_downloaded/$n_files: $fn" );
    download($fn, $staging);
}
my $fetched_size = 0;
find({wanted => \&count_size, no_chdir => 1}, $staging);
if (defined $fo_taxmap) {
    plog( INFO, 'output', "Writing tax map to $fo_taxmap" );
    my $n_seqs = write_taxmap("$staging/$db", $fo_taxmap);
    if ($n_seqs < $meta->{'number-of-sequences'}) {
        # tax map count can be (and often is) greater than the sequence count
        # because of redundant sequences with multiple IDs. Just check that we
        # have *at least* as many entries as sequences
        plog( ERROR, 'output', sprintf(
            "Conversion to tax map returned too few mappings"
            . " (expected %s, got %s)\n",
            $meta->{'number-of-sequences'},
            $n_seqs,
        ));
    }
}
if (defined $fo_fasta) {
    plog( INFO, 'output', "Writing FASTA to $fo_fasta" );
    my $n_seqs = write_fasta("$staging/$db", $fo_fasta);
    if ($n_seqs != $meta->{'number-of-sequences'}) {
        plog( ERROR, 'output', sprintf(
            "Conversion to FASTA returned wrong sequence count"
            . " (expected %s, got %s)\n",
            $meta->{'number-of-sequences'},
            $n_seqs,
        ));
    }
}

if (defined $dir_blast) {

    find({wanted => \&final_move, no_chdir => 1}, $staging);
    $fetched_size = 0;
    find({wanted => \&count_size, no_chdir => 1}, $dir_blast);
    plog( ERROR, 'output', sprintf(
        "Size mismatch of final database: expected %s, got %s\n",
        $expected_size,
        $fetched_size,
    )) if ($expected_size != $fetched_size);

}

exit;

sub write_taxmap {

    my ($db, $fn_out) = @_;

    my @cmd = (
        'blastdbcmd',
        '-db' => "$db",
        '-entry' => 'all',
        '-outfmt' =>  '%a %T %t',
    );
    my $n_seqs = 0;
    my $out;
    if ($gzip) {
        $fn_out .= '.gz' if ($fn_out !~ /\.gz$/);
        $out = IO::Compress::Gzip->new("$fn_out")
            or plog( ERROR, 'output',
                "Failed to open $fn_out for writing: $GzipError" );
    }
    else {
        open $out, '>', $fn_out;
    }
    open my $stream, '-|', @cmd;
    while (my $line = <$stream>) {
        my @f = split ' ', $line;
        say {$out} join "\t",
            (shift @f),
            (shift @f),
            join(' ', @f),
        ;
        ++$n_seqs;
    }
    close $out;
    close $stream
        or plog( ERROR, 'output', "Error writing tax map: $@\n" );
    return $n_seqs;

}

sub write_fasta {

    my ($db, $fn_out) = @_;

    my @cmd = (
        'blastdbcmd',
        '-db' => "$db",
        '-entry' => 'all',
        '-outfmt' =>  '%f',
    );
    my $out;
    if ($gzip) {
        $fn_out .= '.gz' if ($fn_out !~ /\.gz$/);
        $out = IO::Compress::Gzip->new("$fn_out")
            or plog( ERROR, 'output',
                "Failed to open $fn_out for writing: $GzipError" );
    }
    else {
        open $out, '>', $fn_out;
    }
    open my $stream, '-|', @cmd;
    my $n_seqs = 0;
    while (my $line = <$stream>) {
        ++ $n_seqs if ($line =~ /^>/);
        print {$out} $line;
    }
    close $out;
    close $stream
        or plog( ERROR, 'output',  "Error writing FASTA: $@\n");
    return $n_seqs;

}

sub count_size {

    return if (! -f $_);
    return if (basename($_) !~ /^$db\./);
    $fetched_size += -s $_;

}

sub final_move {

    return if (! -f $_);
    move $_, $dir_blast;

}

sub download {
    
    my ($fn_in, $dir_out) = @_;

    $fn_in =~ s/^ftp/rsync/
        if can_run('rsync');

    my $fn_out = join '/', $dir_out, basename($fn_in);
    my $ua = File::Fetch->new(uri => $fn_in);
    my $ret;
    for (0..2) {
        $ret = $ua->fetch(to => $dir_out);
        if (! $ret) {
            unlink $fn_out if (-e $fn_out);
            next;
        }
    }
    plog( ERROR, 'download',
        sprintf("Failed to download %s: %s\n", $fn_in, $ua->error)
    ) if (! $ret);
    my @cmd = (
        $TAR,
        '-C' => $dir_out,
        '-xf' => $fn_out,
    );
    my ($ok, $err) = run(command => \@cmd);
    plog( ERROR, 'download', "Failed to decompress $fn_out: $err\n" )
        if (! $ok);
    unlink $fn_out;

}

sub fetch_meta {

    my ($db) = @_;

    plog( INFO, 'setup', "Fetching metadata for database: $db" );
    my $fn = "$db-nucl-metadata.json";
    my $scratch = File::Temp->newdir(DIR => $tmpdir, CLEANUP => 1);
    my $ua = File::Fetch->new(uri => "$base_url/$fn");
    my $where = $ua->fetch(to => $scratch)
        or plog(
            ERROR,
            'setup',
            "Error downloading metadata for $db: is this a valid database?"
        );
    local $/ = undef;
    open my $in, '<', $where;
    my $data = <$in>;
    close $in;
    my $meta = from_json($data)
        or plog( ERROR, 'setup', "Failed to parse metadata file: $@" );
    return $meta;

}

sub plog {

    my ($lvl, $unit, $msg) = @_;
    if ($lvl ne ERROR) {
        say STDERR "[dl_blast::$unit] $lvl $msg"
            if (! $quiet);
    }
    else {
        die "[dl_blast::$unit] $lvl $msg\n";
    }

}
