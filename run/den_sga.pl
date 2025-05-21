#!/usr/bin/perl
use 5.8.1; use strict; use warnings;

my ($cmd, $path_exe, $path_in, $path_out, $path_genome, $genome, $foreground, $background, $errm);
my ($motif_base, $motif_ext, $n_motifs, $dist_ext, $bin_ext, $sga_log);
my ($i, $foreground_base, $foreground_ext, $background_ext, $genome_prom, $binary);

if(scalar(@ARGV)==0){ die "Wrong arguments!";}

$path_exe=           $ARGV[0]; # path to executable
$path_in=            $ARGV[1]; # input path, motifs
$path_out=           $ARGV[2]; # output path
$path_genome=        $ARGV[3]; # path to executable
$motif_base=         $ARGV[4]; # motif base name
$motif_ext=          $ARGV[5]; # motif extention name
$n_motifs=           $ARGV[5]; # number of motifs
$errm=               $ARGV[7]; # ERRmax threshold
$foreground_base=    $ARGV[8]; # forground fasta without ext
$foreground_ext=     $ARGV[9]; # forground fasta
$background_ext=     $ARGV[10]; # background fasta
$genome_prom=        $ARGV[11]; # genome promoters fasta

$dist_ext = ".dist";
$bin_ext = ".binary";
$sga_log ="sga.log";

if ( -d "$path_out"){
    print "Directory already exist.\n";
}
else{
   mkdir($path_out)
   or die("Can't create directory \"$path_out\": $!\n");
}

$binary = $path_out . $motif_base;
$binary = $binary . $bin_ext;
unlink $binary;

for($i=1;$i<=$n_motifs;$i++)
{

$cmd= "$path_exe/sitega_thr_dist_mat.exe ${path_genome} ${path_in}${motif_base}${i}${motif_ext} ${genome_prom} ${path_out}${motif_base}${i}${dist_ext} ${binary} 0.002 0.0000005 ${sga_log} ab";
print "$cmd\n";
system $cmd;
}

$cmd= "$path_exe/pauc_forback_2motif0s.exe ${path_in}${foreground_base}${foreground_ext} ${path_in}${foreground_base}${background_ext} ${binary} ${n_motifs} ${errm} ${path_out}${foreground_base}.auc_mat ${path_out}${foreground_base}.auc_list ${path_out}${foreground_base}.auc_log1 ${path_out}${foreground_base}.auc_log2 ${path_out}${foreground_base}.roc";
print "$cmd\n";
system $cmd;

