use strict;

my $jobID=$ARGV[0];
my $jobDir=$ARGV[1];
my $mode=$ARGV[2];
my $FrstFile="tertiary_frustration.dat";

my @CA_coords=qx(awk -F "" '{if(\$0~/ CA /) print \$31\$32\$33\$34\$35\$36\$37\$38, \$39\$40\$41\$42\$43\$44\$45\$46, \$47\$48\$49\$50\$51\$52\$53\$54 }' $jobDir/$jobID);
my @Conts_coords=qx(awk '{print \$1, \$2, \$5, \$6, \$7, \$8, \$9, \$10, \$19}' $jobDir/$FrstFile);
my @Positions=qx(awk '{print \$3 }' $jobDir/$jobID\_equivalences.txt);
my @ResChain=qx(awk '{print \$1 }' $jobDir/$jobID\_equivalences.txt);

chomp @CA_coords;
chomp @Conts_coords;
chomp @Positions;
chomp @ResChain;

open(VPs, ">$jobDir/$jobID.vps");
for (my $i=0; $i<@Conts_coords; $i++)
{
  my @splitted=split " ", $Conts_coords[$i];
  my $x=($splitted[5]+$splitted[2])/2;
  my $y=($splitted[6]+$splitted[3])/2;
  my $z=($splitted[7]+$splitted[4])/2;
  print VPs "$splitted[0] $splitted[1] $x $y $z $splitted[8]\n";
}

close VPs;

open (FiveADENS, ">$jobDir/$jobID\_$mode\_5adens");

for (my $i=0; $i<@CA_coords; $i++)
{
  my @Coords=split " ", $CA_coords[$i];
  my $total_density=qx( awk '{ if ( (sqrt( ($Coords[0] - \$3)^2 + ($Coords[1] - \$4)^2 + ($Coords[2] - \$5)^2 )) < 5 ) print}' $jobDir/$jobID.vps | wc -l);
  my $highly_frustrated=qx( awk '{ if ( (sqrt( ($Coords[0] - \$3)^2 + ($Coords[1] - \$4)^2 + ($Coords[2] - \$5)^2 )) < 5 && \$6 < -1  ) print}' $jobDir/$jobID.vps | wc -l);
  my $neutral_frustrated=qx( awk '{ if ( (sqrt( ($Coords[0] - \$3)^2 + ($Coords[1] - \$4)^2 + ($Coords[2] - \$5)^2 )) < 5 && \$6 >= -1 && \$6 < 0.68) print}' $jobDir/$jobID.vps | wc -l);
  my $minimally_frustrated=qx( awk '{ if ( (sqrt( ($Coords[0] - \$3)^2 + ($Coords[1] - \$4)^2 + ($Coords[2] - \$5)^2 )) < 5 && \$6 >= 0.68  ) print}' $jobDir/$jobID.vps | wc -l);
  my $relHighlyFrustratedDensity=0;
  my $relNeutralFrustratedDensity=0;
  my $relMinimallyFrustratedDensity=0;
  chomp $total_density; chomp $highly_frustrated; chomp $neutral_frustrated; chomp $minimally_frustrated;
  
  if($total_density>0)
  {
    $relHighlyFrustratedDensity=$highly_frustrated/$total_density;
    $relNeutralFrustratedDensity=$neutral_frustrated/$total_density;
    $relMinimallyFrustratedDensity=$minimally_frustrated/$total_density;
  }

  print FiveADENS "$Positions[$i] $ResChain[$i] $total_density $highly_frustrated $neutral_frustrated $minimally_frustrated $relHighlyFrustratedDensity $relNeutralFrustratedDensity $relMinimallyFrustratedDensity\n";
}

close FiveADENS;
