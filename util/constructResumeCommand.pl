#!/usr/bin/perl
use warnings;
use strict;
use Class::Struct;
use File::Basename;
use Zach::Util::File qw(OpenFileHandle);

sub ParseLogFile($);
sub ParsePriorFile($);
sub ConcatExecPath($$);

my %FlagSet = ('enable-gradient-descent' => 1);

if(@ARGV < 2){
    die "Usage: ".basename($0). " runBasename run_abc_path [execPath=.] > cmd\n".
        "\trunBasename is the extensionless path to the output of the ABC2run\n".
        "\trun_abc_path is the path to the run_abc executable from here\n".
        "\texecPath is where the previous ABC2 run was called from\n";
}

sub main {
    my($baseName, $runPath, $execPath) = @_;
    $execPath = "." unless(defined $execPath);
    my $resFile = "$baseName.res";
    my $logFile = "$baseName.log";
    my %commandArgs = ParseLogFile($logFile);
    my $priorPath = ConcatExecPath($execPath,$commandArgs{prior});
    my @paramNames = ParsePriorFile("$execPath/$commandArgs{prior}");
    my $rScaleDict = $commandArgs{'updated-proposal-scale'};
    my @initScaleList = split(",",$commandArgs{'initial-proposal-scale'});
    my @scaleList;
    for(my $i = 0; $i < @paramNames; $i++){
        my $param = $paramNames[$i];
        my $scale = $initScaleList[$i % scalar(@initScaleList)];
        if(exists $rScaleDict->{$param}){
            $scale = $rScaleDict->{$param};
        }
        push(@scaleList,$scale);
    }
    $commandArgs{'initial-proposal-scale'} = join(",",@scaleList);
    delete $commandArgs{'updated-proposal-scale'};
    my @cmdParts = ($runPath, "--resume-from $resFile");
    foreach my $key (keys %FlagSet){
        next unless($commandArgs{$key} eq "FALSE");
        delete $commandArgs{$key}
    }
    foreach my $option (sort keys %commandArgs){
        my $value = (exists $FlagSet{$option} ) ? "" : $commandArgs{$option};
        push(@cmdParts,"--$option $value");
    }
    push(@cmdParts,(">| ${baseName}_resumed.res","2>| ${baseName}_resumed.res"));
    print join(" ",@cmdParts),"\n";
} main(@ARGV);

#Parses the log file in reverse to get the most recent values for proposal scaling and
#alpha
#Inputs - a path to an ABC2 log file
#Output - a flag keyed hash of flag values
sub ParseLogFile($){
    my $file = shift;
    my $fh = OpenFileHandle("tac $file |", "logFile","ERROR");
    my %commandArgs;
    $commandArgs{'updated-proposal-scale'} = {};
    while(my $line = <$fh>){
        chomp($line);
        if($line =~ m/updated/){ #Proposal 
            #Only care a bout values on the first chain
            next unless($line =~ m/Chain 0\)/);
            if($line =~ m/scale for ([a-zA-z_\-0-9]+) proposals updated to (-?[0-9.]+)/){
                my $param = $1;
                my $value = $2;
                #Only want to lates update
                next if(exists $commandArgs{'updated-proposal-scale'}->{$param});
                $commandArgs{'updated-proposal-scale'}->{$param} = exp($value);
            } elsif($line =~ m/variance handling updated to (-?[0-9.]+)/){
                my $value = $1;
                #Only want to lates update
                next if(exists $commandArgs{'initial-sim-variance-alpha'});
                $commandArgs{'initial-sim-variance-alpha'} = 10**$value;
            } elsif($line =~ m/increment updated to (-?[0-9.]+)/){
                my $value = $1;
                next if(exists $commandArgs{'initial-temperature-increment'});
                $commandArgs{'initial-temperature-increment'} = $value;
            }
        }
        if($line =~ m/\t(.+):\t(.+)$/){
            my $option = $1;
            my $value = $2;
            next if(exists $commandArgs{$option});
            $commandArgs{$option} = $value;
        }
    }
    close($fh);
    return %commandArgs;
}


#Parses the prior file (specified in the log) to determine the complete set of parameters
#   this is done to ensure the proposal-scaling-values are properly ordered
#Inputs - a prior file
#Output - an  array of parameter names in alphabetical order
sub ParsePriorFile($){
    my $fh = OpenFileHandle(shift,"prior","ERROR");
    my $priorNameLine = <$fh>;
    my @params;
    while(my $line = <$fh>){
        chomp($line);
        next if($line =~ m/^Init[AL]/);
        my($param) = split(/\t/,$line);
        push(@params,$param);
    }
    return sort @params;
}

#Given an executable path and a path to a file prepends the executable path unless the
#later is an abosulte path
#Inputs - an exectuable path
#       - another path
#Output - a string representing the location of the second file
sub ConcatExecPath($$){
    my ($preFix,$suffix) = @_;
    if($suffix =~ m/^\//){
        return $suffix;
    }
    return "$preFix/$suffix";
}
