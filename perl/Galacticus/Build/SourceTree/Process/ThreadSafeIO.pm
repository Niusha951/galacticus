# Contains a Perl module which implements locking Fortran IO to make it thread-safe.

package Galacticus::Build::SourceTree::Process::ThreadSafeIO;
use strict;
use warnings;
use Cwd;
use lib $ENV{'GALACTICUS_EXEC_PATH'}."/perl";
use Fortran::Utils;

# Insert hooks for our functions.
$Galacticus::Build::SourceTree::Hooks::processHooks       {'threadSafeIO'} = \&Lock_IO;
$Galacticus::Build::SourceTree::Hooks::processDependencies{'threadSafeIO'} = [ "metaPropertyDatabase", "stateStore", "functionClass", "stateStorable" ];

sub Lock_IO {
    ## <workaround type="gfortran" PR="92836" url="https:&#x2F;&#x2F;gcc.gnu.org&#x2F;bugzilla&#x2F;show_bug.cgi=92836">
    ##  <description>Internal file I/O in gfortran can be non-thread safe.</description>
    ## </workaround>
    
    # Get the tree.
    my $tree = shift();
    # Check if debugging is required.
    my $lock = 0;
    if ( exists($ENV{'GALACTICUS_FCFLAGS'}) ) {
	$lock = 1
	    if ( grep {$_ eq "-DTHREADSAFEIO"} split(" ",$ENV{'GALACTICUS_FCFLAGS'}) );
    }
    return
	unless ( $lock );
    # Do not apply thread locking inside error reporting functions. These can be triggered from anywhere in the code - including
    # within another I/O critical section which would then lead to deadlocking. Since at this point we are in an error condition
    # anyway, we can just take our chances with thread race conditions.
    my $depth  = 0;
    my $module = $tree;
    while ( $module ) {
	return
	    if ( $module->{'type'} eq "module" && $module->{'name'} eq "Error" );
	$module = &Galacticus::Build::SourceTree::Walk_Tree($module,\$depth);
    }
    # Walk the tree, looking for IO statements.
    my @ioStatements = ( 'open', 'close', 'read', 'write' );
    my $node         = $tree;
    $depth           = 0;
    while ( $node ) {
	if ( $node->{'type'} eq "code" ) {
	    my $newContent    ;
	    my $inIO       = 0;
	    my $inCritical = 0;
	    open(my $content,"<",\$node->{'content'});
	    while ( ! eof($content) ) {
		&Fortran::Utils::Get_Fortran_Line($content,my $rawLine, my $processedLine, my $bufferedComments);
		# Detect entry and exit to existing critical sections - we don't want to add locks inside such a section as that
		# would result in deadlocks.
		$inCritical = 1
		    if ( $processedLine =~ m/^\s*!\$omp\s+critical\s*\(gfortranInternalIO\)\s*$/ );
		$inCritical = 0
		    if ( $processedLine =~ m/^\s*!\$omp\s+end\s+critical\s*\(gfortranInternalIO\)\s*$/ );
		# Convert "FoX_DOM_Access" critical sections to "gfortranInternalIO". This still prevents multi-threaded access to
		# the FoX library, but also restricts access to gfortran internal IO within the FoX library.
		$rawLine =~ s/FoX_DOM_Access/gfortranInternalIO/
		    if ( $rawLine =~ m/^\s*!\$omp\s+(end\s+)??critical\s*\(FoX_DOM_Access\)\s*$/ );
		# Detect IO statements and lock if not already in a critical section.
		my $isIO = 0;
		foreach my $ioStatement ( @ioStatements ) {
		    if ( $processedLine =~ m/^\s*(!\$)??\s*$ioStatement\s*\(/ ) {
			$isIO = 1;
			last;
		    }
		}
		$isIO = 1
		    if ( $processedLine =~ m/^\s*(!\$)??\s*call\s+flush\s*\(/i );
		if ( $isIO && ! $inIO && ! $inCritical ) {
		    $inIO        = 1;
		    $newContent .= "    !\$omp critical(gfortranInternalIO)\n";
		}
		if ( ! $isIO && $inIO && ! $inCritical ) {
		    $inIO        = 0;
		    $newContent .= "    !\$omp end critical(gfortranInternalIO)\n";
		}
		$newContent .= $rawLine;
	    }
	    close($content);
	    $newContent .= "    !\$omp end critical(gfortranInternalIO)\n"
		if ( $inIO && ! $inCritical );
	    $node->{'content'} = $newContent;
	}
	$node = &Galacticus::Build::SourceTree::Walk_Tree($node,\$depth);
    }
}

1;
