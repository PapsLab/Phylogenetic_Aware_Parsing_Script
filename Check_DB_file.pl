#!/usr/bin/perl

   use strict;
   use AnyDBM_File;
   use Fcntl;

   if (@ARGV < 1) {
      print STDERR "\nUsage $0 dbfile\n\n"; exit;
   }

   my ($dbfile) = @ARGV;
   my (%DB, $key, $val);

   if (!dbmopen(%DB, $dbfile, 0444) ) {
      print STDERR "\nError opening file $dbfile\n\n"; exit;
   }

   while ( ($key, $val) = each %DB) {
      print "$key,$val\n";
   }

   dbmclose(%DB);