#!/usr/bin/perl -w

use Statistics::Basic qw(:all nofill);
$args=$ARGV[0];
#$block_size=$ARGV[1];
#$iter=$ARGV[2];
$taskset=$ARGV[1];
$sched=$ARGV[2];
$total=$ARGV[3];
$fast=$ARGV[4];
$reps=10;

$strict = $ARGV[6];
$steal = $ARGV[5];

$ENV{NX_SCHEDULE} = $sched;
$ENV{NX_HP_FROM} = 1;
$ENV{NX_HP_TO} = $fast;
$min=10000000000;
$max=0;
my $vectorValues = vector( )->set_size($reps);
if($sched eq "bf" || $sched eq "heft") { 
for (my $i=0; $i < ${reps}; $i++) {
  print "Iteration: $i Running with $sched and $taskset\n";
  $OUTPUT=`taskset -c ${taskset} ../cholesky ${args}`;
  $sum+=$OUTPUT;
  $vectorValues->insert($OUTPUT);
  if($OUTPUT > $max) {
     $max=$OUTPUT;
  }
  if($OUTPUT < $min) {
     $min=$OUTPUT;
  }
}
my $std = stddev($vectorValues);
my $mean = mean($vectorValues);
#print "STDDEV = $std";
$res=$sum/$reps;
open(my $fd, ">>${sched}_$args.txt");
print $fd "$total $fast $res $std $min $max $mean\n";
#print $res;
}
if($sched eq "botlev" || $sched eq "cpath") {
#   for (my $strict=0; $strict <=1; $strict++) {
      $ENV{NX_STRICTB} = $strict; 
#      for (my $steal=0; $steal<=1; $steal++) {
         $ENV{NX_STEALB} = $steal;
         for (my $i=0; $i < ${reps}; $i++) {
            print "Iteration: $i Running with $sched and $taskset STRICT: $strict STEAL: $steal\n";
            $OUTPUT=`taskset -c ${taskset} ../cholesky ${args}`;
            $sum+=$OUTPUT;
            $vectorValues->insert($OUTPUT);
            if($OUTPUT > $max) {
               $max=$OUTPUT;
            }
            if($OUTPUT < $min) {
               $min=$OUTPUT;
            }
         }
         my $std = stddev($vectorValues);
         my $mean = mean($vectorValues);
         $res=$sum/$reps;
         open(my $fd, ">>ss_${steal}_strict${strict}/${sched}_$args.txt");
         print $fd "$total $fast $res $std $min $max $mean\n";

#      }
      $sum = 0;
      $min=10000000000;
      $max=0;
#   }
}
else {
   print "sched equals to ${sched} and is not recognized\n";
}
