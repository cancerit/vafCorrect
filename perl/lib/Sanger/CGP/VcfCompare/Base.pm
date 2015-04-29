package Sanger::CGP::VcfCompare::Base;

use strict;
use warnings FATAL => 'all';

use Data::Dumper;
use Carp;
use autodie qw(:all);
use English qw( -no_match_vars );
use List::Util qw(first);
use Const::Fast qw(const);
use Capture::Tiny ':all';
use Try::Tiny qw(try catch finally);
use File::Copy;

use Log::Log4perl;
my $log = Log::Log4perl->get_logger(__PACKAGE__);

use Cwd 'abs_path';
my $prog_path = abs_path($PROGRAM_NAME);

use Sanger::CGP::Database::Conn;

sub new {
    my ($class) = @_;
    my $self = { };
    bless $self, $class;
    return $self;
}

sub connection {
  my ($self, $conn) = @_;
  if(! defined( $self->{'conn'})) {
    if($conn) {
      $self->{'conn'} = $conn;
    }
    else {
      $log->debug("Attempting to setup connection.");
      try {
        $self->{'conn'} = Sanger::CGP::Database::Conn->new($self->db_type);
      }catch{
        $log->logcroak($_) if($_);
      };
      $log->debug("Connection succeeded.");
    }
  }
  return $self->{'conn'};
}

sub db_type {
  my ($self, $db_type) = @_;
  if($db_type) {
    if($self->{'db_type'}) {
      $log->logcroak('DB_TYPE is already defined as: '.$db_type);
    }
    if($db_type ne 'test' && $db_type ne 'live') {
      $log->logcroak('DB_TYPE can only be test|live, '.$db_type.' is not valid');
    }
    $log->debug("Setting DB_TYPE to $db_type.");
    $self->{'db_type'} = $db_type;
  }
  return $self->{'db_type'};
}

sub commit {
  return shift->connection->commit();
}

sub rollback {
  return shift->connection->rollback();
}

sub clear_connection {
  my ($self) = @_;
  undef $self->{'conn'};
  return;
}

sub rsync {
  my ($from, $to, $full) = @_;

  #Attempt to make path for local copys
  if($to !~ m/:/ && !-e $to) {
    if($to =~ m/\/$/) {
      make_path($to);
    }
  }

  if(-e $to || $from =~ m/:/ || $to =~ m/:/) {
    my $size_only = q{-I };
    $size_only = '--size-only ' if(!$full);
    my $command = 'rsync -L '.$size_only.$from." ".$to;
    _run_external($command, 'rsync', undef, undef, 'no_data');
  }
  else {
    copy($from, $to) or $log->logcroak("Failed to copy $from -> $to (no previous file and local copy so not rsync)");
  }
  return;
}

sub _run_external{
  my ($command) = @_;
  my @args = ();
  $log->warn( ">>>>> $command\n");
  #Use capture tiny to run an external command
  my ($stdout, $stderr, $exit) = capture {
    system( $command, @args );
  };
  if($exit != 0){ #Non zero error code
    $log->logcroak("Error in external command: '$command'\n",$stderr);
  }
  $log->warn( "<<<<<\n");
  chomp($stdout);
  $log->debug("job submission STDOUT: ".$stdout);
  $log->debug("job submission STDERR: ".$stderr);
  my @results = split(/\n/,$stdout);
  return \@results;
}

# recursively cleanups folder and underlying substructure  
sub cleanup_results {
    my ($self,$dir)=@_;
    $log->logcroak("Unable to find cleanup dir:$dir") if(! -d $dir);
    remove(\1,"$dir");
    $log->logcroak("Unable to remove cleanup dir:$dir") if(-d $dir);
    $log->debug("Dir: $dir cleaned successfully");
    return;
}

=head3 clear_path

=over

Implemented as File::Path remove_tree (rmtree) is naff
the idea is to remove the content of all of the directories
and then use remove_tree to remove the directories
needed to handle very big structures from devil

=back

=cut
sub clear_path {
  my ($self, $root_path) = @_;
  my ($dir_count, $file_count) = (0,0);
  my @dirs = ($root_path);
  $dir_count++;
  if((scalar @dirs) > 0) {
    my $curr_path = shift @dirs;
    opendir my $CLEAN, $curr_path or $log->logcroak("Path does not exist $curr_path");
    while(my $thing = readdir $CLEAN) {
      next if($thing =~ m/^\.{1,2}$/);
      my $full_path = $curr_path.'/'.$thing;
      if(-d $full_path) {
        push @dirs, $full_path;
        $dir_count++;
      }
      else {
        unlink $full_path or $log->logcroak("Unable to delete $full_path");
        $file_count++;
      }
    }
    closedir $CLEAN;
  }
  remove_tree($root_path);
  return ($dir_count, $file_count);
}





1;
