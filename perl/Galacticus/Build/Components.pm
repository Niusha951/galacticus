# Contains a Perl module which implements processing of "component" directives in the Galacticus build system.

package Component;
my $galacticusPath;
if ( exists($ENV{"GALACTICUS_ROOT_V091"}) ) {
 $galacticusPath = $ENV{"GALACTICUS_ROOT_V091"};
 $galacticusPath .= "/" unless ( $galacticusPath =~ m/\/$/ );
} else {
 $galacticusPath = "./";
}
unshift(@INC, $galacticusPath."perl"); 
use strict;
use warnings;
use utf8;
use DateTime;
use Switch;
use Data::Dumper;
use Text::Table;
use Sort::Topological qw(toposort);
use Scalar::Util 'reftype';
use Carp 'verbose';
$SIG{ __DIE__ } = sub { Carp::confess( @_ ) };
require Galacticus::Build::Hooks;

# Insert hooks for our functions.
%Hooks::moduleHooks = 
    (
     %Hooks::moduleHooks,
     component => {parse => \&Components_Parse_Directive, generate => \&Components_Generate_Output}
    );

# Include debugging code.
my $debugging                           = 0;

# Global verbosity level.
my $verbosityLevel                      = 1;

# Switch to control gfortran workarounds.
my $workaround                          = 1;

# Records of the longest component and method names.
my $classNameLengthMax                  = 0;
my $implementationNameLengthMax         = 0;
my $fullyQualifiedNameLengthMax         = 0;
my $propertyNameLengthMax               = 0;
my $linkedDataNameLengthMax             = 0;
my $implementationPropertyNameLengthMax = 0;

sub Components_Parse_Directive {
    # Parse content for a "component" directive.
    my $buildData = shift;

    # Assert that we have a prefix, currentDocument and directive.
    die("Galacticus::Build::Components::Components_Parse_Directive: no currentDocument present")
	unless ( exists($buildData->{'currentDocument'}           ) );
    die("Galacticus::Build::Components::Components_Parse_Directive: no name present"           )
	unless ( exists($buildData->{'currentDocument'}->{'name' }) );
    die("Galacticus::Build::Components::Components_Parse_Directive: no class present"          )
	unless ( exists($buildData->{'currentDocument'}->{'class'}) );

    # Construct an ID for this component.
    my $componentID = ucfirst($buildData->{'currentDocument'}->{'class'}).ucfirst($buildData->{'currentDocument'}->{'name'});
		
    # Store a copy of the component's defining document.
    $buildData->{'components'}->{$componentID} = $buildData->{'currentDocument'};

}

sub Components_Generate_Output {
    # Generate output for a "component" directive.
    my $buildData = shift;

    # Construct a list of all component names.
    @{$buildData->{'componentIdList'}} = keys(%{$buildData->{'components'}});
    
    # Iterate over all functions, calling them with the build data object.
    &{$_}($buildData)
	foreach (
	    # Construct null implementations for all component classes.
	    \&Construct_Null_Components                              ,
	    # Construct component class list and membership lists for all classes.
	    \&Construct_Class_Membership                             ,
	    # Distribute class default values to all members of a class.
	    \&Distribute_Class_Defaults                              ,
	    # Set defaults for unspecified attributes.
	    \&Set_Default_Attributes                                 ,
	    # Construct linked data for component methods.
	    \&Construct_Linked_Data                                  ,
	    # Generate the nodeComponent type.
	    \&Generate_Node_Component_Type                           ,
	    # Generate component class types.
	    \&Generate_Component_Classes                             ,
	    # Sort component implementations such that they will be defined after any component of which they are an extension.
	    \&Sort_Implementations                                   ,
	    # Generate implementation types.
	    \&Generate_Implementations                               ,
	    # Generate the treeNode object.
	    \&Generate_Tree_Node_Object                              ,
	    # Create an initialization method.
	    \&Generate_Initialization_Function                       ,
	    # Generate a finalization method.
	    \&Generate_Finalization_Function                         ,
	    # Generate functions to map other functions over components.
	    \&Generate_Map_Functions                                 ,
	    # Generate functions to dump nodes.
	    \&Generate_Node_Dump_Function                            ,
	    # Generate functions to output nodes.
	    \&Generate_Node_Output_Functions                         ,
	    # Generate functions to serialize/deserialize nodes to/from arrays.
	    \&Generate_Node_Serialization_Functions                  ,
	    # Generate functions to get property names from a supplied index.
	    \&Generate_Node_Property_Name_From_Index_Function        ,
	    # Generate functions to serialize/deserialize nodes to/from arrays.
	    \&Generate_Node_ODE_Initialization_Functions             ,
	    # Generate function to copy one node to another.
	    \&Generate_Node_Copy_Function                            ,
	    # Generate function to move one node to another.
	    \&Generate_Node_Move_Function                            ,
	    # Generate a tree node creation function.
	    \&Generate_Tree_Node_Creation_Function                   ,
	    # Generate a tree node destruction function.
	    \&Generate_Tree_Node_Destruction_Function                ,
	    # Generate type name functions.
	    \&Generate_Type_Name_Functions                           ,
	    # Generate component assignment function.
	    \&Generate_Component_Assignment_Function                 ,
	    # Generate component class destruction functions.
	    \&Generate_Component_Class_Destruction_Functions         ,
	    # Generate component class removal functions.
	    \&Generate_Component_Class_Removal_Functions             ,
	    # Generate component class move functions.
	    \&Generate_Component_Class_Move_Functions                ,
	    # Generate dump functions for each component class.
	    \&Generate_Component_Class_Dump_Functions                ,
	    # Generate output functions for each component class.
	    \&Generate_Component_Class_Output_Functions              ,
	    # Generate functions to determine if an implementation is active.
	    \&Generate_Is_Active_Functions                           ,
	    # Generate component implementation destruction functions.
	    \&Generate_Component_Implementation_Destruction_Functions,
	    # Generate ODE solver initialization functions.
	    \&Generate_ODE_Initialization_Functions                  ,
	    # Generate dump functions for each implementation.
	    \&Generate_Implementation_Dump_Functions                 ,
	    # Generate output functions for each implementation.
	    \&Generate_Implementation_Output_Functions               ,
	    # Generate functions to get the name of properties from an index.
	    \&Generate_Implementation_Name_From_Index_Functions      ,
	    # Generate serialization/deserialization functions for each implementation.
	    \&Generate_Implementation_Serialization_Functions        ,
	    # Generate component count methods.
	    \&Generate_Component_Count_Functions                     ,
	    # Generate component get methods.
	    \&Generate_Component_Get_Functions                       ,
	    # Generate component destruction functions.
	    \&Generate_Component_Destruction_Functions               ,
	    # Generate component creation functions.
	    \&Generate_Component_Creation_Functions                  ,
	    # Generate component class default value functions.
	    \&Generate_Component_Class_Default_Value_Functions       ,
	    # Generate functions for getting/setting/rating value via a deferred function.
	    \&Generate_Deferred_GSR_Function                         ,
	    # Generate functions for getting/setting/rating value directly.
	    \&Generate_GSR_Functions                                 ,
	    # Insert code for type-definitions.
	    \&Insert_Type_Definitions                                ,
	    # Generate objects that record which type of component will be used by default.
	    \&Generate_Default_Component_Sources                     ,
	    # Generate records of which component implementations are selected.
	    \&Generate_Active_Implementation_Records                 ,
	    # Generate deferred procedure pointers.
	    \&Generate_Deferred_Procedure_Pointers                   ,
	    # Generate required null binding functions.
	    \&Generate_Null_Binding_Functions                        ,
	    # Insert the "contains" line.
	    \&Insert_Contains
	);

    # Insert all functions into content.
    $buildData->{'content'} .= join("\n",@{$buildData->{'code'}->{'functions'}})."\n";
    
    # Insert include statements to bring in all functions associated with components.
    my @includeDependencies;
    foreach my $component ( @{$buildData->{'componentIdList'}} ) {
     	if ( exists($buildData->{'components'}->{$component}->{'functions'}) ) {
     	    $buildData->{'content'} .= "  include \"".$buildData->{'components'}->{$component}->{'functions'}."\"\n";
     	    push(@includeDependencies,$buildData->{'components'}->{$component}->{'functions'});
     	}
    }

    # Create a Makefile to specify dependencies on these include files.
    open(makeFile,">./work/build/Makefile_Component_Includes");
    print makeFile "./work/build/objects.nodes.components.Inc:".join("",map {" ./work/build/".$_} @includeDependencies)
	if ( scalar(@includeDependencies) > 0 );
    close(makeFile);

}

sub Get_Type {
    # Returns the type of a method of pipe.
    my $buildData->{'currentDocument'} = shift;
    # Assume scalar type by default
    my $type = "scalar";
    # If a type is specified, then return it instead.
    if ( exists($buildData->{'currentDocument'}->{'type'}) ) {$type = $buildData->{'currentDocument'}->{'type'}};
    return $type;
}

sub Get_Suffix {
    # Returns the suffix for a method.
    my $methodType = shift;
    # Determine the suffix to use.
    my $suffix = "";
    switch ( $methodType ) {
	case ( "scalar"      ) {
	    $suffix = "";
	}
	case ( "array"       ) {
	    $suffix = "_Array";
	}
	case ( "history"     ) {
	    $suffix = "_History";
	}
	case ( "abundances"  ) {
	    $suffix = "_Abundances";
	}
	case ( "chemicalAbundances"  ) {
	    $suffix = "_Chemical_Abundances";
	}
	case ( "keplerOrbit" ) {
	    $suffix = "_Kepler_Orbit";
	}
	else {
	    die("Build_Include_File.pl: unrecognized method type");
	}
    }
    return $suffix;
}

sub dataObjectName {
    # Construct and return the name of the object to use for data of given type and rank.
    my $dataObject = shift;
    # Variable to store the object name.
    my $name = "nodeData";
    # Extract the type.
    if ( exists($dataObject->{'type'}) ) {
	switch ( $dataObject->{'type'} ) {
	    case ( "integer"     ) {$name .= "Integer"            }
	    case ( "logical"     ) {$name .= "Logical"            }
	    case ( "real"        ) {$name .= "Double"             }
	    case ( "chemicals"   ) {$name .= "ChemicalAbundances" }
	    case ( "abundances"  ) {$name .= "Abundances"         }
	    case ( "history"     ) {$name .= "History"            }
	    case ( "keplerOrbit" ) {$name .= "KeplerOrbit"        }
	    else {die "Build_Include_File.pl::dataObjectName: 'type' specifier is unknown"}
	}
    } else {
	die "Build_Include_File.pl::dataObjectName: no 'type' specifier present";
    }
    if ( exists($dataObject->{'rank'}) ) {
	switch ( $dataObject->{'rank'} ) {
	    case ( 0 ) {$name .= "Scalar"}
	    case ( 1 ) {$name .= "1d"}
	    else {die "Build_Include_File.pl::dataObjectName: 'rank' specifier is unknown"}
	}
    } else {
	die "Build_Include_File.pl::dataObjectName: no 'rank' specifier present";
    }
    $name .= "Evolvable" if ( $dataObject->{'isEvolvable'} eq "true" );
    return $name;
}

sub dataObjectPrimitiveName {
    # Construct and return the name and attributes of the primitive data class to use for data of given type and rank.
    my $dataObject = shift;
    my %options;
    if ( $#_ >= 1 ) {(%options) = @_};

    # Variables to store the object name and attributes.
    my $name;
    my $type;
    my @attributes;
    # Extract the type.
    if ( exists($dataObject->{'type'}) ) {
	switch ( $dataObject->{'type'} ) {
	    case ( "integer"     ) {
		$name = "integer"                 ;
		$type = "Integer"                 ;
	    }
	    case ( "logical"     ) {
		$name = "logical"                 ;
		$type = "Logical"                 ;
	    }
	    case ( "real"        ) {
		$name = "double precision"        ;
		$type = "Double"                  ;
	    }
	    case ( "abundances"  ) {
		$name = "type(abundances        )";
		$type = "Abundances"              ;
	    }
	    case ( "history"     ) {
		$name = "type(history)"           ;
		$type = "History"                 ;
	    }
	    case ( "chemicals"   ) {
		$name = "type(chemicalAbundances)";
		$type = "ChemicalAbundances"      ;
	    }
	    case ( "keplerOrbit" ) {
		$name = "type(keplerOrbit)"       ;
		$type = "KeplerOrbit"             ;
	    }
	    else {die "Build_Include_File.pl::dataObjectPrimitiveName: 'type' specifier [".$dataObject->{'type'}."] is unknown"}
	}
    } else {
	die "Build_Include_File.pl::dataObjectPrimitiveName: no 'type' specifier present";
    }
    if ( exists($dataObject->{'rank'}) ) {
	switch ( $dataObject->{'rank'} ) {
	    case ( 0 ) {} # Nothing to do.
	    case ( 1 ) {
		push(@attributes,"dimension(:)");
		push(@attributes,"allocatable" ) unless ( exists($options{'matchOnly'}) && $options{'matchOnly'} == 1 );
	    }
	    else {die "Build_Include_File.pl::dataObjectPrimitiveName: 'rank' specifier is unknown"}
	}
    } else {
	die "Build_Include_File.pl::dataObjectPrimitveName: no 'rank' specifier present";
    }
    my $attributeList = "";
    $attributeList = ", ".join(", ",@attributes) if ( scalar(@attributes) > 0 );
    return ($name,$type,$attributeList);
}

sub Data_Object_Definition {
    # Construct and return the name and attributes of the primitive data class to use for data of given type and rank.
    my $dataObject = shift;
    my %options;
    if ( $#_ >= 1 ) {(%options) = @_};

    # Variables to store the object name and attributes.
    my $intrinsicName;
    my $type         ;
    my $label        ;
    my @attributes   ;
    # Extract the type.
    if ( exists($dataObject->{'type'}) ) {
	switch ( $dataObject->{'type'} ) {
	    case ( "integer"     ) {
		$intrinsicName = "integer"           ;
		$label         = "Integer"           ;
	    }
	    case ( "logical"     ) {
		$intrinsicName = "logical"           ;
		$label         = "Logical"           ;
	    }
	    case ( "real"        ) {
		$intrinsicName = "double precision"  ;
		$label         = "Double"            ;
	    }
	    case ( "abundances"  ) {
		$intrinsicName = "type"              ;
		$type          = "abundances"        ;
		$label         = "Abundances"        ;
	    }
	    case ( "history"     ) {
		$intrinsicName = "type"              ;
		$type          = "history"           ;
		$label         = "History"           ;
	    }
	    case ( "chemicals"   ) {
		$intrinsicName = "type"              ;
		$type          = "chemicalAbundances";
		$label         = "ChemicalAbundances";
	    }
	    case ( "keplerOrbit" ) {
		$intrinsicName = "type"              ;
		$type          = "keplerOrbit"       ;
		$label         = "KeplerOrbit"       ;
	    }
	    else {die "Build_Include_File.pl::dataObjectPrimitiveName: 'type' specifier [".$dataObject->{'type'}."] is unknown"}
	}
    } else {
	die "Build_Include_File.pl::dataObjectPrimitiveName: no 'type' specifier present";
    }
    if ( exists($dataObject->{'rank'}) ) {
	switch ( $dataObject->{'rank'} ) {
	    case ( 0 ) {} # Nothing to do.
	    case ( 1 ) {
		push(@attributes,"dimension(:)");
		push(@attributes,"allocatable" )
		    unless ( exists($options{'matchOnly'}) && $options{'matchOnly'} == 1 );
	    }
	    else {die "Build_Include_File.pl::dataObjectPrimitiveName: 'rank' specifier is unknown"}
	}
    } else {
	die "Build_Include_File.pl::dataObjectPrimitveName: no 'rank' specifier present";
    }
    # Construct the definitions.
    my $dataDefinition;
    $dataDefinition  ->{'intrinsic' }  = $intrinsicName;
    $dataDefinition  ->{'type'      }  = $type
	if ( defined($type) );
    @{$dataDefinition->{'attributes'}} = @attributes
	if ( @attributes    );
    # Return the data definition and label.
    return ($dataDefinition,$label);
}

sub pad {
    # Pad a string to give nicely aligned formatting in the output code.
    die("pad() requires two arguments")
	unless (scalar(@_) == 2);
    my $text       = shift;
    my $padLength  = $_[0];
    my $paddedText = $text." " x ($padLength-length($text));
    return $paddedText;
}

sub padImplementationProperty {
    # Pad a string to give nicely aligned formatting in the output code.
    my $text     = shift;
    my @extraPad = @{$_[0]};
    my $padLength = $implementationPropertyNameLengthMax+$extraPad[0];
    $padLength = $extraPad[1] if ($extraPad[1] > $padLength);
    my $paddedText = $text." " x ($padLength-length($text));
    return $paddedText;
}

sub padComponentClass {
    # Pad a string to give nicely aligned formatting in the output code.
    my $text     = shift;
    my @extraPad = @{$_[0]};
    my $padLength = $classNameLengthMax+$extraPad[0];
    $padLength = $extraPad[1] if ($extraPad[1] > $padLength);
    my $paddedText = $text." " x ($padLength-length($text));
    return $paddedText;
}

sub padImplementation {
    # Pad a string to give nicely aligned formatting in the output code.
    my $text       = shift;
    my @extraPad   = @{$_[0]};
    my $padLength  = $implementationNameLengthMax+$extraPad[0];
    $padLength     = $extraPad[1] if ($extraPad[1] > $padLength);
    my $paddedText = $text." " x ($padLength-length($text));
    return $paddedText;
}

sub padFullyQualified {
    # Pad a string to give nicely aligned formatting in the output code.
    my $text       = shift;
    my @extraPad   = @{$_[0]};
    my $padLength  = $fullyQualifiedNameLengthMax+$extraPad[0];
    $padLength     = $extraPad[1] if ($extraPad[1] > $padLength);
    my $paddedText = $text." " x ($padLength-length($text));
    return $paddedText;
}

sub padProperty {
    # Pad a string to give nicely aligned formatting in the output code.
    my $text     = shift;
    my @extraPad = @{$_[0]};
    my $padLength = $propertyNameLengthMax+$extraPad[0];
    $padLength = $extraPad[1] if ($extraPad[1] > $padLength);
    my $paddedText = $text." " x ($padLength-length($text));
    return $paddedText;
}

sub padLinkedData {
    # Pad a string to give nicely aligned formatting in the output code.
    my $text     = shift;
    my @extraPad = @{$_[0]};
    my $padLength = $linkedDataNameLengthMax+$extraPad[0];
    $padLength = $extraPad[1] if ($extraPad[1] > $padLength);
    my $paddedText = $text." " x ($padLength-length($text));
    return $paddedText;
}

sub Construct_Class_Membership {
    # Generates a null implementation for each component class and makes it the default if no default is specified.
    my $buildData = shift;

    # Iterate over component IDs.
    foreach my $componentID ( @{$buildData->{'componentIdList'}} ) {
	# Get the component.
	my $component               = $buildData->{'components'}->{$componentID};
	# Get the name of the component class.
	my $componentClass          = $component->{'class'};
	# Get the name of the implementation.
	my $componentImplementation = $component->{'name'};
	# Append the component name to the list of members for its class.
	push(
	    @{$buildData->{'componentClasses'}->{$componentClass}->{'members'}},
	    $component->{'name'}
	    );
	# Record the longest class name, implementation name and component ID.
	$classNameLengthMax          = length($componentClass         )
	    if (length($componentClass         ) > $classNameLengthMax         );
	$fullyQualifiedNameLengthMax = length($componentID            )
	    if (length($componentID            ) > $fullyQualifiedNameLengthMax);
	$implementationNameLengthMax = length($componentImplementation)
	    if (length($componentImplementation) > $implementationNameLengthMax);
    }

    # Construct a list of component classes.
    @{$buildData->{'componentClassList'}} = keys(%{$buildData->{'componentClasses'}});
}

sub Distribute_Class_Defaults {
    # Distribute class defaults to all members of a class.
    my $buildData = shift;

    # Iterate over component classes.
    foreach my $componentClassName ( @{$buildData->{'componentClassList'}} ) {
	# Initialize hash to hold defaults.
	my %classDefaults;
	# Iterate over class members.
	foreach my $componentName ( @{$buildData->{'componentClasses'}->{$componentClassName}->{'members'}} ) {
	    # Get the component.
	    my $componentID = ucfirst($componentClassName).ucfirst($componentName);
	    my $component   = $buildData->{'components'}->{$componentID};
	    # Iterate over the properties of this implementation.
	    foreach my $methodName ( keys(%{$component->{'methods'}->{'method'}}) ) {
		# Get the method.
		my $method = $component->{'methods'}->{'method'}->{$methodName};
		# Check for class defaults.
		if ( exists($method->{'classDefault'}) ) {
		    my $code;
		    if ( ref($method->{'classDefault'}) && exists($method->{'classDefault'}->{'content'}) ) {
			$code = $method->{'classDefault'}->{'content'};
		    } else {
			$code = $method->{'classDefault'};
		    }
		    if ( exists($classDefaults{$componentID.$methodName}) ) {
			die("Distribute_Class_Defaults: inconsistent class defaults for ".$componentID." ".$method)
			    unless ($code eq $classDefaults{$componentID.$methodName}->{'code'} );
		    } else {
			$classDefaults{$componentID.$methodName}->{'code'} = $code;
		    }
		    if ( ref($method->{'classDefault'}) && exists($method->{'classDefault'}->{'modules'}) ) {
			my @requiredModules = split(/\s*,\s*/,$method->{'classDefault'}->{'modules'});
			push(
			    @{$classDefaults{$componentID.$methodName}->{'modules'}},
			    @requiredModules
			    );
		    }
		    if ( ref($method->{'classDefault'}) && exists($method->{'classDefault'}->{'count'}) ) {
			if ( exists($classDefaults{$componentID.$methodName}->{'count'}) ) {
			    die("Distribute_Class_Defaults: inconsistent class default counts for ".$componentID." ".$method)
				unless ($method->{'classDefault'}->{'count'} eq $classDefaults{$componentID.$methodName}->{'count'} );
			} else {
			    $classDefaults{$componentID.$methodName}->{'count'} = $method->{'classDefault'}->{'count'};
			}
		    }
		}
	    }
	}
	# Iterate over class members.
	foreach my $componentName ( @{$buildData->{'componentClasses'}->{$componentClassName}->{'members'}} ) {
	    # Get the component.
	    my $componentID = ucfirst($componentClassName).ucfirst($componentName);
	    my $component   = $buildData->{'components'}->{$componentID};
	    # Iterate over the properties of this implementation.
	    foreach my $methodName ( keys(%{$component->{'methods'}->{'method'}}) ) {
		# Get the method.
		my $method = $component->{'methods'}->{'method'}->{$methodName};
		# Set class default if available.
		$method->{'classDefault'} = $classDefaults{$componentID.$methodName}
		    if ( exists($classDefaults{$componentID.$methodName}) );
	    }
	}
    }

    # Construct a list of component classes.
    @{$buildData->{'componentClassList'}} = keys(%{$buildData->{'componentClasses'}});
}

sub Construct_Linked_Data {
    # Generates linked data for component methods.
    my $buildData = shift;

    # Iterate over all component implementations.
    foreach my $componentID ( @{$buildData->{'componentIdList'}} ) {
	# Iterate over all methods belonging to this component.	
	if ( exists($buildData->{'components'}->{$componentID}->{'methods'}) ) {
	    foreach my $methodName ( keys(%{$buildData->{'components'}->{$componentID}->{'methods'}->{'method'}}) ) {
		my $method = $buildData->{'components'}->{$componentID}->{'methods'}->{'method'}->{$methodName};
		# Record the longest method name.
		$propertyNameLengthMax = length($methodName) if (length($methodName) > $propertyNameLengthMax);
		# Check for a pre-defined linkedData element.
		my $linkedDataName;
		if ( exists($method->{'linkedData'}) ) {		    
		    # A linkedData element has been explicitly declared, write an informational message.
		    print " -> INFO: linkedData can be created automatically for ".$methodName." method of the ".lcfirst($componentID)." component\n";
		    # Get the name of the linked data and the data itself.
		    $linkedDataName = $method->{'linkedData'};
		    my $linkedData  = $buildData->{'components'}->{$componentID}->{'content'}->{'data'}->{$linkedDataName};
		    # Set the isEvolvable flag on the linked data.
		    $linkedData->{'isEvolvable'} = $method->{'attributes'}->{'isEvolvable'};
		    # Create a copy of the linked data.
		    $method->{'data'} = $linkedData;
		    $method->{'type'} = $linkedData->{'type'};
		    $method->{'rank'} = $linkedData->{'rank'};
		} else {
		    # No linkedData element is explicitly declared. Therefore, we must have type, and rank specified.
		    foreach my $requiredElement ( "type", "rank" ) {
			die("No ".$requiredElement." was specified for ".$methodName." method of the ".lcfirst($componentID)." component")
			    unless ( exists($method->{$requiredElement}) );
		    }
		    # If no isVirtual element is present, assume "no" by default.
		    $method->{'isVirtual'} = "no"
			unless ( exists($method->{'isVirtual'}) );
		    # Copy the attributes to the data element.
		    $method->{'data'} = 
		    {
			type        => $method->{'type'      }                 ,
			rank        => $method->{'rank'      }                 ,
			isEvolvable => $method->{'attributes'}->{'isEvolvable'}
		    };
		    # Unless this method is virtual, create a linked data object for it.
		    unless ( $method->{'isVirtual'} eq "yes" ) {
			# Write a message.
			print " -> Creating linked data object for ".$methodName." method of the ".lcfirst($componentID)." component\n";
			# Create the linked data name.
			$linkedDataName = $methodName."Data";
			# Create the linked data object.
			$method->{'linkedData'} = $linkedDataName;
			$buildData->{'components'}->{$componentID}->{'content'}->{'data'}->{$linkedDataName} = $method->{'data'};
		    }
		}
		# Record the longest linked data name.
		if ( defined($linkedDataName) ) {
		    $linkedDataNameLengthMax = length($linkedDataName) 
			if (length($linkedDataName) > $linkedDataNameLengthMax);
		    # Record the longest possible implementation plus property name length.
		    my $implementationPropertyName = $buildData->{'components'}->{$componentID}->{'name'}.$linkedDataName;
		    $implementationPropertyNameLengthMax = length($implementationPropertyName) 
			if (length($implementationPropertyName) > $implementationPropertyNameLengthMax);
		}
	    }
	}
    }
}

sub Construct_Null_Components {
    # Generates a null implementation for each component class and makes it the default if no default is specified.
    my $buildData = shift;

    # Iterate over components to determine which classes need a null case building.
    my %classes;
    foreach my $componentID ( @{$buildData->{'componentIdList'}} ) {
	# Get the component object.
	my $component = $buildData->{'components'}->{$componentID};
	# Initialize this class if it hasn't been seen before.
	unless ( exists($classes{$component->{'class'}}) ) {
	    $classes{$component->{'class'}}->{'hasNull'   } = 0;
	    $classes{$component->{'class'}}->{'hasDefault'} = 0;
	}
	# Record if a null component already exists.
	$classes{$component->{'class'}}->{'hasNull'   } = 1
	    if ( $component->{'name'} eq "null" );
	# Record if a default is already specified.
	$classes{$component->{'class'}}->{'hasDefault'} = 1
	    if ( $component->{'isDefault'} eq "yes" );
    }

    # Iterate over classes, creating null components as necessary.
    foreach my $class ( keys(%classes) ) {       
	# Test for pre-existing null component.
	if ( $classes{$class}->{'hasNull'} == 0 ) {
	    # No pre-existing null component is present, so simply insert one into the build data.
	    my $componentID = ucfirst($class)."Null";
	    my $isDefault   = "no";
	    $isDefault      = "yes"
		if ( $classes{$class}->{'hasDefault'} == 0 );
	    $buildData->{'components'}->{$componentID}->{'class'    } = $class;
	    $buildData->{'components'}->{$componentID}->{'name'     } = "null";
	    $buildData->{'components'}->{$componentID}->{'isDefault'} = $isDefault;
	    # Append this new component ID to the component ID list.
	    push(@{$buildData->{'componentIdList'}},$componentID);
	    # Display a message.
	    if ( $verbosityLevel >= 1 ) {
		print " -> Adding null implementation ";
		print "as default "
		    if ( $classes{$class}->{'hasDefault'} == 0 );
		print "for ".$class." class\n";
	    }
	} elsif ( $verbosityLevel >= 1 ) {
	    # Advise that null components don't need to be explicitly specified.
	    print " -> INFO: a pre-existing null component exists for ".$class." class,\n";
	    print " ->       but would be built automatically.\n";
	}
    }
}

sub Set_Default_Attributes{
    # Set any missing attributes to default values.
    my $buildData = shift;

    # Iterate over component implementations.
    foreach my $componentID ( @{$buildData->{'componentIdList'}} ) {
	# Get the component.
	my $component = $buildData->{'components'}->{$componentID};
	# Add a fully-qualified name to the component.
	$component->{'fullyQualifiedName'} = $componentID;
	# If a create function is specified, set it to be non-deferred by default.
	$component->{'createFunction'}->{'isDeferred'} = "false"
	    if ( exists($component->{'createFunction'}) && ! exists($component->{'createFunction'}->{'isDeferred'}) );
	# Iterate over methods.
	foreach my $methodName ( keys(%{$component->{'methods'}->{'method'}}) ) {
	    # Get the method.
	    my $method = $component->{'methods'}->{'method'}->{$methodName};
	    # Add the method's name.
	    $method->{'name'} = $methodName;
	    # Binding.
	    $method->{'attributes'}->{'bindsTo'} = "component"
		unless ( exists($method->{'attributes'}->{'bindsTo'}) );
	    # Auto-creation.
	    $method->{'attributes'}->{'createIfNeeded'} = "false"
		unless ( exists($method->{'attributes'}->{'createIfNeeded'}) );
	    # Deferred status.
	    $method->{'attributes'}->{'isDeferred'} = ""
		unless ( exists($method->{'attributes'}->{'isDeferred'}) );
	    # Generic status.
	    $method->{'attributes'}->{'makeGeneric'} = "false"
		unless ( exists($method->{'attributes'}->{'makeGeneric'}) );
	    # isEvolable synonym.
	    $method->{'attributes'}->{'isRatetable'} = $method->{'attributes'}->{'isEvolvable'};
	    # Rate function.
	    $method->{'rateFunction'} = $componentID.ucfirst($methodName)."Rate"
		unless ( exists($method->{'rateFunction'}) );
	    # Get function.
	    if ( exists($method->{'getFunction'}) ) {
		# A getFunction element was specified.
		if ( defined(reftype($method->{'getFunction'})) ) {
		    # The getFunction element contains structure, so simple set its bindsTo element if not already defined.
		    $method->{'getFunction'}->{'bindsTo'} = "component"
			unless ( exists($method->{'getFunction'}->{'bindsTo'}) );
		} else {
		    # The getFunction element is simply the function name. Replace with a structure with default binding.
		    $method->{'getFunction'} = 
		    {
			content => $method->{'getFunction'},
			bindsTo => "component"
		    };
		}
		# Since getFunction was specified, we will not need to build a get function.
		$method->{'getFunction'}->{'build'} = "false";
	    } else {
		# No getFunction element was specified, assign a default function and record that a get function must be built.
		$method->{'getFunction'} = 
		{
		    content => lcfirst($componentID).ucfirst($methodName)."Get",
		    bindsTo => "component"                               ,
		    build   => "true"
		};
	    }
	    # Set function.
	    if ( exists($method->{'setFunction'}) ) {
		# A setFunction element was specified.
		if ( defined(reftype($method->{'setFunction'})) ) {
		    # The setFunction element contains structure, so simple set its bindsTo element if not already defined.
		    $method->{'setFunction'}->{'bindsTo'} = "component"
			unless ( exists($method->{'setFunction'}->{'bindsTo'}) );
		} else {
		    # The setFunction element is simply the function name. Replace with a structure with default binding.
		    $method->{'setFunction'} = 
		    {
			content => $method->{'setFunction'},
			bindsTo => "component"
		    };
		}
		# Since setFunction was specified, we will not need to build a set function.
		$method->{'setFunction'}->{'build'} = "false";
	    } else {
		# No setFunction element was specified, assign a default function and record that a set function must be built.
		$method->{'setFunction'} = 
		{
		    content => lcfirst($componentID).ucfirst($methodName)."Set",
		    bindsTo => "component"                               ,
		    build   => "true"
		};
	    }
	}
    }
}

sub Generate_Node_Component_Type{
    # Generate the top-level object in the class hierachy: nodeComponent.
    my $buildData = shift;
    # Define type-bound functions.
    my @typeBoundFunctions = 
	(
	 {type => "procedure", name => "type"                   , function => "Node_Component_Generic_Type"            },
	 {type => "procedure", name => "host"                   , function => "Node_Component_Host_Node"               },
	 {type => "procedure", name => "destroy"                , function => "Node_Component_Generic_Destroy"         },
	 {type => "procedure", name => "serializeCount"         , function => "Node_Component_Serialize_Count_Zero"    },
	 {type => "procedure", name => "serializeValues"        , function => "Node_Component_Serialize_Null"          },
	 {type => "procedure", name => "serializeRates"         , function => "Node_Component_Serialize_Null"          },
	 {type => "procedure", name => "serializeScales"        , function => "Node_Component_Serialize_Null"          },
	 {type => "procedure", name => "deserializeValues"      , function => "Node_Component_Deserialize_Null"        },
	 {type => "procedure", name => "deserializeRates"       , function => "Node_Component_Deserialize_Null"        },
	 {type => "procedure", name => "deserializeScales"      , function => "Node_Component_Deserialize_Null"        },
	 {type => "procedure", name => "odeStepRatesInitialize" , function => "Node_Component_ODE_Step_Initialize_Null"},
	 {type => "procedure", name => "odeStepScalesInitialize", function => "Node_Component_ODE_Step_Initialize_Null"},
	 {type => "procedure", name => "dump"                   , function => "Node_Component_Dump_Null"               },
	 {type => "procedure", name => "outputCount"            , function => "Node_Component_Output_Count_Null"       },
	 {type => "procedure", name => "outputNames"            , function => "Node_Component_Output_Names_Null"       },
	 {type => "procedure", name => "output"                 , function => "Node_Component_Output_Null"             }
	);
    # Specify the data content.
    my @dataContent =
	(
	 {
	     intrinsic  => "type",
	     type       => "treeNode",
	     attributes => [ "pointer", "public" ],
	     variables  => [ "hostNode" ]
	 }
	);
    # Create the nodeComponent class.
    $buildData->{'types'}->{'nodeComponent'} = {
	name           => "nodeComponent",
	comment        => "A class for components in \\glspl{node}.",
	isPublic       => "true",
	boundFunctions => \@typeBoundFunctions,
	dataContent    => \@dataContent
    };
    push(@{$buildData->{'typesOrder'}},'nodeComponent');


    # Insert an interface for assignment.
    if ( $workaround == 1 ) {
	$buildData->{'content'} .= "  interface assignment(=)\n";
	$buildData->{'content'} .= "    module procedure Node_Component_Assign\n";
	$buildData->{'content'} .= "  end interface assignment(=)\n";
    }

}

sub Generate_Component_Classes{
    # Generate object types for each component class.
    my $buildData = shift;
 
    # Iterate over all component classes.
    my %classGetDefaults;
    foreach my $componentClass ( @{$buildData->{'componentClassList'}} ) {
	# Define a hash to record which methods have already been created.
	my %methodsCreated;

	# Create a list for type-bound functions.
	my @typeBoundFunctions;

  	# Insert definitions for each method associated with a component implementation of this component class.
    	foreach my $implementationName ( @{$buildData->{'componentClasses'}->{$componentClass}->{'members'}} ) {
	    # Construct a fully-qualified name for this implementation.
	    my $componentName = ucfirst($componentClass).ucfirst($implementationName);
	    # Iterate over methods beloning to this implementation.
	    foreach my $methodName ( keys(%{$buildData->{'components'}->{$componentName}->{'methods'}->{'method'}}) ) {
		# Get the method.
		my $method = $buildData->{'components'}->{$componentName}->{'methods'}->{'method'}->{$methodName};
		# Create functions to set/get/evolve each method as necessary.
		if ( 
		       $method->{'attributes'}->{'isGettable' } eq "true"
		    || $method->{'attributes'}->{'isSettable' } eq "true"
		    || $method->{'attributes'}->{'isEvolvable'} eq "true"
		    )
		{
		    # Name of function being processed.
		    my $functionName;
		    # Get a fully-qualified type identfier for this method.
		    (my $intrinsic,my $type,my $attributes) = &dataObjectPrimitiveName($method);
		    $type .= $method->{'rank'}."InOut";
		    # Record the null bindings needed.
		    $buildData->{'nullMethods'}->{$componentClass}->{"Integer0In"} =
		    {
			type   => "integer",
			rank   => 0        ,
			intent => "in"
		    };
		    $buildData->{'nullMethods'}->{$componentClass}->{$type       } = 
		    {
			type   => $method->{'type'},
			rank   => $method->{'rank'},
			intent => "inout"
		    };
		    # Create the "isSettable" function.
		    $functionName = $methodName."IsSettable";
		    unless ( exists($methodsCreated{$functionName}) ) {
			push(
			    @typeBoundFunctions,
			    {type => "procedure", pass => "nopass", name => $functionName, function => "Boolean_False"}
			    );
			$methodsCreated{$functionName} = 1;
		    }
		    # Handle set functions and related functions.
		    unless ( $method->{'attributes'}->{'isSettable'} eq "false" ) {
			# Create a "set" function if one does not already exist.
			$functionName = $methodName."Set";
			unless ( exists($methodsCreated{$functionName}) ) {
			    my $boundTo;
			    if ( $method->{'setFunction'}->{'bindsTo'} eq "componentClass" )
			    {
				# A setFunction was specified that binds to the component class, so bind to it here.
				$boundTo = $method->{'setFunction'}->{'content'};
			    } else {
				# Create a binding to a null function here. 
				$boundTo = $componentClass."NullBindingSet".$type;
			    }
			    push(
				@typeBoundFunctions,
				{type => "procedure", name => $functionName, function => $boundTo}
				);
			    $methodsCreated{$functionName} = 1;
			}
		    }

		    # Handle evolve functions.
		    unless ( $method->{'attributes'}->{'isEvolvable'} eq "false" ) {
			# Create the "count" function.
			$functionName = $methodName."Count";
			unless ( exists($methodsCreated{$functionName}) ) {
			    push(
				@typeBoundFunctions,
				{type => "procedure", name => $functionName, function => $componentClass."NullBindingInteger0In"}
				);
			    $methodsCreated{$functionName} = 1;
			}
			# Create the "rate" function.
			$functionName = $methodName."Rate";
			unless (  exists($methodsCreated{$functionName}) ) {
			    unless ( 
				   $method->{'attributes'}->{'isDeferred' } =~ m/rate/
				&& $method->{'attributes'}->{'bindsTo'    } eq "top"
				) {
				my $boundTo;
				if ( $method->{'attributes'}->{'createIfNeeded'} eq "true" ) 
				{
				    $boundTo = $componentClass.ucfirst($methodName)."Rate";
				} else {
				    $boundTo = $componentClass."NullBindingRate".$type;
				}
				push(
				    @typeBoundFunctions,
				    {type => "procedure", name => $functionName, function => $boundTo}
				    );
			    }
			    # Create a "scale" function unless this is a virtual method.
			    push(
				@typeBoundFunctions,
				{type => "procedure", name => $methodName."Scale", function => $componentClass."NullBindingSet".$type}
				)
				unless ( $method->{'isVirtual'} eq "true" );
			    $methodsCreated{$functionName} = 1;
			}
		    }
		    # Add any bindings which bind at the component class level.
		    if ( exists($buildData->{'components'}->{$componentName}->{'bindings'}) ) {
			foreach ( @{$buildData->{'components'}->{$componentName}->{'bindings'}->{'binding'}} ) {
			    push(
				@typeBoundFunctions,
				{type => "procedure", name => $_->{'method'}, function => $_->{'nullFunction'}}
				)
				if ( $_->{'bindsTo'} eq "componentClass" );
			}
		    }
		}
	    }
	}
	# Create the type.
	$buildData->{'types'}->{'nodeComponent'.ucfirst($componentClass)} = {
	    name           => "nodeComponent".ucfirst($componentClass),
	    comment        => "Type for the {\\tt ".$componentClass."} component class.",
	    isPublic       => "true",
	    extends        => "nodeComponent",
	    boundFunctions => \@typeBoundFunctions,
	};
	push(@{$buildData->{'typesOrder'}},'nodeComponent'.ucfirst($componentClass));
    }
}

sub Sort_Implementations {
    # Sort component implementations such that they will be defined after any component of which they are an extension.
    my $buildData = shift;
    # Construct depenencies for type extension.
    my %dependencies;
    foreach my $componentID ( @{$buildData->{'componentIdList'}} ) {
	# Get the component.
	my $component = $buildData->{'components'}->{$componentID};
	# Get the parent class.
	my $componentClassName = $component->{'class'};
    	# By default this component will be an extension of the base "nodeComponent" class.
    	my $extensionOf = "nodeComponent".ucfirst($componentClassName);
    	# If it specifies a particular component that it should extend, use that instead.
    	$extensionOf = ucfirst($component->{'extends'}->{'class'}).ucfirst($component->{'extends'}->{'name'})
	    if (exists($component->{'extends'}));
	push(@{$dependencies{$extensionOf}},$componentID);
    }
    # Perform a dependency sort on the implementation list.
    @{$buildData->{'componentIdList'}} = toposort(sub { @{$dependencies{$_[0]} || []}; }, \@{$buildData->{'componentIdList'}});

}

sub Generate_Implementations {
    # Generate a type for each component implementation.
    my $buildData = shift;
    # Create classes for each specific implementation.
    foreach my $componentID ( @{$buildData->{'componentIdList'}} ) {
	# Get the implementation.
	my $component          = $buildData->{'components'}->{$componentID};
	# Get the parent class.
	my $componentClassName = $component->{'class'};
    	# By default this component will be an extension of the base "nodeComponent" class.
    	my $extensionOf = "nodeComponent".ucfirst($componentClassName);
    	# If it specifies a particular component that it should extend, use that instead.
    	$extensionOf = "nodeComponent".ucfirst($component->{'extends'}->{'class'}).ucfirst($component->{'extends'}->{'name'})
	    if (exists($component->{'extends'}));
     	# Create data objects to store all of the linked data for this component.
	my @dataContent;
    	foreach ( keys(%{$component->{'content'}->{'data'}}) ) {
    	    my $type = &dataObjectName($component->{'content'}->{'data'}->{$_});
	    push(
		@dataContent,
		{
		    intrinsic  => "type",
		    type       => $type,
		    variables  => [ $_ ]
		}
		);
    	}
	# Create a list for type-bound functions.
	my @typeBoundFunctions;
	# Add binding for deferred create function set function.
	push(
	    @typeBoundFunctions,
	    {type => "procedure", pass => "nopass", name => "createFunctionSet", function => $componentID."CreateFunctionSet"}
	    )
	    if ( 
	    	   exists($component->{'createFunction'})
		&&        $component->{'createFunction'}->{'isDeferred'} eq "true" 
	    );
     	# If this component has bindings defined, scan through them and create an appropriate method.
    	if ( exists($component->{'bindings'}) ) {
    	    foreach ( @{$component->{'bindings'}->{'binding'}} ) {
		push(
		    @typeBoundFunctions,
		    {type => "procedure", name => $_->{'method'}, function => $_->{'function'}},
		    );
	    }
	}
	# Iterate over methods.
	foreach my $methodName ( keys(%{$component->{'methods'}->{'method'}}) ) {
	    # Get the method.
	    my $method = $component->{'methods'}->{'method'}->{$methodName};
	    push(
		@typeBoundFunctions,
		{type => "procedure", pass => "nopass", name => $methodName."IsGettable", function => "Boolean_".ucfirst($method->{'attributes'}->{'isGettable'})},
		{type => "procedure", pass => "nopass", name => $methodName."IsSettable", function => "Boolean_".ucfirst($method->{'attributes'}->{'isSettable'})}
		);
	}
	# Create the type.
	$buildData->{'types'}->{'nodeComponent'.ucfirst($componentID)} = {
	    name           => "nodeComponent".ucfirst($componentID),
	    comment        => "Class for the ".$component->{'name'}." implementation of the ".$componentClassName." component.",
	    isPublic       => "true",
	    extends        => $extensionOf,
	    boundFunctions => \@typeBoundFunctions,
	    dataContent    => \@dataContent
	};
	push(@{$buildData->{'typesOrder'}},'nodeComponent'.ucfirst($componentID));
    }
}

sub Generate_Active_Implementation_Records{
    # Generate records of which component implementations are selected.
    my $buildData = shift;
    # Create a table.
    my $recordTable = Text::Table->new(
	{
	    is_sep => 1,
	    body   => "  logical :: "
	},
	{
	    align  => "left"
	},
	{
	    is_sep => 1,
	    body   => "=.false."
	}
    );
    # Iterate over all component implementations.
    foreach ( @{$buildData->{'componentIdList'}} ) {
	$recordTable->add("nodeComponent".$_."IsActive");
    }
    # Insert into the document.
    $buildData->{'content'} .= "  ! Records of which component implementations are active.\n";
    $buildData->{'content'} .= $recordTable->table()."\n";
}

sub Generate_Deferred_Procedure_Pointers {
    # Generate deferred procedure pointers.
    my $buildData = shift;
    # Initialize record of pointers which have been created.
    my %createdPointers;
    # Insert comment.
    $buildData->{'content'} .= "  ! Procedure pointers for deferred custom functions.\n";
    # Initialize data content.
    my @dataContent;
    # Iterate over component implementations.
    foreach my $componentID ( @{$buildData->{'componentIdList'}} ) {
	# Get the component.
	my $component = $buildData->{'components'}->{$componentID};
	# Get the component class name.
	my $componentClassName = $component->{'class'};
	# Create pointer for deferred create functions.
	push(
	    @dataContent,
	    {
		intrinsic  => "procedure",
		type       => "",
		attributes => [ "pointer" ],
		variables  => [ $componentID."CreateFunction" ]
	    }
	    )
	    if (
		   exists($component->{'createFunction'})
		&&        $component->{'createFunction'}->{'isDeferred'} eq "true"
	    );
	# Iterate over methods.
	foreach my $methodName ( keys(%{$component->{'methods'}->{'method'}}) ) {
	    my $method = $component->{'methods'}->{'method'}->{$methodName};
	    unless ( $method->{'attributes' }->{'isDeferred'} eq "" ) {
		my $selfType = "generic";
		$selfType = $component->{'class'}
		    unless ( $method->{'attributes'}->{'bindsTo'} eq "top" );
		(my $dataObject, my $label) = &Data_Object_Definition($method);
		my $dataType = $label.$method->{'rank'};
		# Iterate over attributes.
		foreach ( "get", "set", "rate" ) {
		    # Determine if this attribute is deferred and has not yet had a procedure pointer created.
		    if (
			   $method->{'attributes' }->{'isDeferred'} =~ m/$_/ 
			&& $method->{'attributes' }->{'is'.ucfirst($_).'table'} eq "true"
			&& ! exists($createdPointers{$componentClassName.ucfirst($methodName).$_})
			) {
			# Construct the template function.
			my $template = $selfType."NullBinding".ucfirst($_).$dataType."InOut";
			$template = $componentID.ucfirst($methodName).ucfirst($_)
			    if ( $_ eq "get" );
			# Generate the procedure pointer and a boolean to indicate if is has been attached.
			push(
			    @dataContent,
			    {
				intrinsic  => "procedure",
				type       => $template,
				attributes => [ "pointer" ],
				variables  => [ $componentClassName.ucfirst($methodName).ucfirst($_)."Deferred" ]
			    },
			    {
				intrinsic  => "logical",
				variables  => [ $componentClassName.ucfirst($methodName).ucfirst($_)."IsAttachedValue=.false." ]
			    },
			    );
			# Add the required null method to the list.
			$buildData->{'nullMethods'}->{$selfType}->{$dataType."InOut"} =
			{
			    type   => $method->{'type'},
			    rank   => $method->{'rank'},
			    intent => "inout"
			};
			# Record that this procedure pointer has been created.
			$createdPointers{$componentClassName.ucfirst($methodName).$_} = 1;
		    }		    
		}
	    }
	}
    }
    # Insert data content.
    $buildData->{'content'} .= &Format_Variable_Defintions(\@dataContent, indent => 2)."\n";
}

sub Generate_Default_Component_Sources{
    # Generate records of which component implementations are selected.
    my $buildData = shift;
    # Create a table.
    my $recordTable = Text::Table->new(
	{
	    is_sep => 1,
	    body   => "  class(nodeComponent"
	},
	{
	    align  => "left"
	},
	{
	    is_sep => 1,
	    body   => "), allocatable, public :: default"
	},
	{
	    align  => "left"
	}
    );
    # Iterate over all component implementations.
    foreach ( @{$buildData->{'componentClassList'}} ) {
	$recordTable->add(ucfirst($_),ucfirst($_)."Component");
    }
    # Insert into the document.
    $buildData->{'content'} .= "  ! Objects that will record which type of each component is to be used by default.\n";
    $buildData->{'content'} .= $recordTable->table()."\n";
    # Create a table for the class types.
    $recordTable = Text::Table->new(
	{
	    is_sep => 1,
	    body   => "  type(nodeComponent"
	},
	{
	    align  => "left"
	},
	{
	    is_sep => 1,
	    body   => ") :: "
	},
	{
	    align  => "left"
	}
    );
    # Iterate over all component implementations.
    foreach ( @{$buildData->{'componentClassList'}} ) {
	$recordTable->add(ucfirst($_),ucfirst($_)."Class");
    }
    # Insert into the document.
    $buildData->{'content'} .= "  ! Objects that will record which type of each component is to be used by default.\n";
    $buildData->{'content'} .= $recordTable->table()."\n";
}

sub Generate_Tree_Node_Object{
    # Generate the treeNode object.
    my $buildData = shift;

    # Define bound functions.
    my @typeBoundFunctions = 
	(
	 {type => "procedure", name => "type"                     , function => "Tree_Node_Type"                          },
	 {type => "procedure", name => "index"                    , function => "Tree_Node_Index"                         },
	 {type => "procedure", name => "indexSet"                 , function => "Tree_Node_Index_Set"                     },
	 {type => "procedure", name => "uniqueID"                 , function => "Tree_Node_Unique_ID"                     },
	 {type => "procedure", name => "initialize"               , function => "treeNodeInitialize"                      },
	 {type => "procedure", name => "destroy"                  , function => "treeNodeDestroy"                         },
	 {type => "procedure", name => "removeFromHost"           , function => "Tree_Node_Remove_From_Host"              },
	 {type => "procedure", name => "removeFromMergee"         , function => "Tree_Node_Remove_From_Mergee"            },
	 {type => "procedure", name => "isPrimaryProgenitor"      , function => "Tree_Node_Is_Primary_Progenitor"         },
	 {type => "procedure",                                    , function => "Tree_Node_Is_Primary_Progenitor_Of_Index"},
	 {type => "procedure",                                    , function => "Tree_Node_Is_Primary_Progenitor_Of_Node" },
	 {type => "generic"  , name => "isPrimaryProgenitorOf"    , function => ["Tree_Node_Is_Primary_Progenitor_Of_Index","Tree_Node_Is_Primary_Progenitor_Of_Node"]},
	 {type => "procedure", name => "isOnMainBranch"           , function => "Tree_Node_Is_On_Main_Branch"             },
	 {type => "procedure", name => "isSatellite"              , function => "Tree_Node_Is_Satellite"                  },
	 {type => "procedure", name => "lastSatellite"            , function => "Tree_Node_Get_Last_Satellite"            },
	 {type => "procedure", name => "earliestProgenitor"       , function => "Tree_Node_Get_Earliest_Progenitor"       },
	 {type => "procedure", name => "mergesWith"               , function => "Tree_Node_Merges_With_Node"              },
	 {type => "procedure", name => "walkBranch"               , function => "Tree_Node_Walk_Branch"                   },
	 {type => "procedure", name => "walkBranchWithSatellites" , function => "Tree_Node_Walk_Branch_With_Satellites"   },
	 {type => "procedure", name => "walkTree"                 , function => "Tree_Node_Walk_Tree"                     },
	 {type => "procedure", name => "walkTreeUnderConstruction", function => "Tree_Node_Walk_Tree_Under_Construction"  },
	 {type => "procedure", name => "walkTreeWithSatellites"   , function => "Tree_Node_Walk_Tree_With_Satellites"     }
	);
    # Add data content.
    my @dataContent =
	(
	 {
	     intrinsic  => "integer",
	     type       => "kind=kind_int8",
	     variables  => [ "indexValue", "uniqueIdValue" ]
	 },
	 {
	     intrinsic  => "type",
	     type       => "treeNode",
	     attributes => [ "pointer", "public" ],
	     variables  => [ "parent", "firstChild", "sibling", "firstSatellite", "mergeTarget", "firstMergee", "siblingMergee", "formationNode" ]
	 },
	 {
	     intrinsic  => "logical",
	     attributes => [ "public" ],
	     variables  => [ "isPhysicallyPlausible" ]
	 }
	);
    foreach ( @{$buildData->{'componentClassList'}} ) {
	push(
	    @dataContent,
	    {
	     intrinsic  => "class",
	     type       => "nodeComponent".ucfirst($_),
	     attributes => [ "allocatable", "dimension(:)" ],
	     variables  => [ "component".ucfirst($_) ],
	     comment    => "A generic ".$_." object."
	    }
	    );
    }
    # Create the tree node class.
    $buildData->{'types'}->{'treeNode'} = {
	name           => "treeNode",
	comment        => "A class for \\glspl{node} in merger trees.",
	isPublic       => "true",
	boundFunctions => \@typeBoundFunctions,
	dataContent    => \@dataContent
    };
    push(@{$buildData->{'typesOrder'}},"treeNode");
}

sub Generate_Initialization_Function {
    # Generate an initialization function.
    my $buildData = shift;
    # Generate the function code.
    my $functionCode;
    $functionCode .= "  subroutine Galacticus_Nodes_Initialize()\n";
    $functionCode .= "    !% Initialize the \\glc\\ object system.\n";
    $functionCode .= "    use Input_Parameters\n";
    $functionCode .= "    use ISO_Varying_String\n";
    $functionCode .= "    use Memory_Management\n";
    $functionCode .= "    implicit none\n";
    # Generate variables.
    my @dataContent =
	(
	 {
	     intrinsic  => "type",
	     type       => "varying_string",
	     variables  => [ "methodSelection", "message" ]
	 }
	);
    foreach my $componentClass ( @{$buildData->{'componentClassList'}} ) {
    	foreach my $implementationName ( @{$buildData->{'componentClasses'}->{$componentClass}->{'members'}} ) {
	    my $fullName = ucfirst($componentClass).ucfirst($implementationName);
	    push(
		@dataContent,
		{
		    intrinsic  => "type",
		    type       => "nodeComponent".$fullName,
		    variables  => [ "default".$fullName."Component" ]
		}
		);
    	}
    }
    $functionCode .= &Format_Variable_Defintions(\@dataContent)."\n";
    # Iterate over all component classes.
    $buildData->{'content'} .= "  ! Parameters controlling output.\n\n";
    foreach my $componentClass ( @{$buildData->{'componentClassList'}} ) {
	# Identify the default implementation.
    	my $defaultMethod;
    	foreach my $implementationName ( @{$buildData->{'componentClasses'}->{$componentClass}->{'members'}} ) {
    	    my $fullName = ucfirst($componentClass).ucfirst($implementationName);
    	    $defaultMethod = $implementationName
		if ( $buildData->{'components'}->{$fullName}->{'isDefault'} eq "yes" );
    	}
	die("No default method was found for ".$componentClass." class")
	    unless ( defined($defaultMethod) );
	# Insert a function call to get the parameter controlling the choice of implementation for this class.
    	$functionCode .= "    call Get_Input_Parameter('treeNodeMethod".padComponentClass(ucfirst($componentClass)."'",[1,0]).",methodSelection,defaultValue='".padImplementation($defaultMethod."'",[1,0]).")\n";
    	foreach my $implementationName ( @{$buildData->{'componentClasses'}->{$componentClass}->{'members'}} ) {
    	    my $fullName  = ucfirst($componentClass).ucfirst($implementationName);
	    my $component = $buildData->{'components'}->{$fullName};
    	    $functionCode .= "    if (methodSelection == '".padImplementation($implementationName."'",[1,0]).") then\n";
	    $functionCode .= "       allocate(default".padComponentClass(ucfirst($componentClass)."Component",[9,0]).",source=default".padFullyQualified($fullName."Component",[9,0]).")\n";
	    $functionCode .= "       nodeComponent".padFullyQualified($fullName."IsActive",[8,0])."=.true.\n";
	    until ( $fullName eq "" ) {
		if ( exists($buildData->{'components'}->{$fullName}->{'extends'}) ) {
		    $fullName = ucfirst($buildData->{'components'}->{$fullName}->{'extends'}->{'class'}).ucfirst($buildData->{'components'}->{$fullName}->{'extends'}->{'name'});
		    $functionCode .= "       nodeComponent".padFullyQualified($fullName."IsActive",[8,0])."=.true.\n";
		} else {
		    $fullName = "";
		}
	    }
	    $functionCode .= "    end if\n";
	    # Insert code to read and parameters controlling outputs.
	    foreach my $methodName ( keys(%{$component->{'methods'}->{'method'}}) ) {
		my $method = $component->{'methods'}->{'method'}->{$methodName};
		# Check for output and output condition.
		if (
		    exists($method->{'output'}               )                       &&
		    exists($method->{'output'}->{'condition'})                       &&
		           $method->{'output'}->{'condition'} =~ m/\[\[([^\]]+)\]\]/
		    ) 
		{
		    my $parameterName = $1;
		    $functionCode .= "call Get_Input_Parameter('".$parameterName."',".$parameterName.",defaultValue=.false.)\n";
		    $buildData->{'content'} .= "  logical :: ".$parameterName."\n";
		}
	    }

    	}
    	$functionCode .= "    if (.not.allocated(default".padComponentClass(ucfirst($componentClass)."Component",[9,0]).")) then\n";
    	$functionCode .= "       message='unrecognized method \"'//methodSelection//'\" for \"".$componentClass."\" component'\n";
    	$functionCode .= "       call Galacticus_Error_Report('Galacticus_Nodes_Initialize',message)\n";
    	$functionCode .= "    end if\n";
    }
    $buildData->{'content'} .= "\n";
    $functionCode .= "    return\n";
    $functionCode .= "  end subroutine Galacticus_Nodes_Initialize\n";	
    # Insert into the function list.
    push(
	@{$buildData->{'code'}->{'functions'}},
	$functionCode
	);
}

sub Generate_Finalization_Function {
    # Generate a finalization function.
    my $buildData = shift;
    # Create a table for the deallocation code.
    my $table = Text::Table->new(
	{
	    is_sep => 1,
	    body   => "    deallocate(default"
	},
	{
	    align  => "left"
	},
	{
	    is_sep => 1,
	    body   => ")"
	},
	);
    # Populate the table.
    foreach ( @{$buildData->{'componentClassList'}} ) {
	$table->add(ucfirst($_)."Component");
    }
    # Generate the function code.
    my $functionCode;
    $functionCode .= "  subroutine Galacticus_Nodes_Finalize()\n";
    $functionCode .= "    !% Finialize the \\glc\\ object system.\n";
    $functionCode .= "    implicit none\n\n";
    $functionCode .= $table->table();
    $functionCode .= "    return\n";
    $functionCode .= "  end subroutine Galacticus_Nodes_Finalize\n";
    # Insert into the function list.
    push(
	@{$buildData->{'code'}->{'functions'}},
	$functionCode
	);
}

sub Generate_Map_Functions {
    # Generate functions to map other functions over components.
    my $buildData = shift;

    # Function for mapping a void function.
    # Generate variables.
    my @dataContent =
	(
	 {
	     intrinsic  => "class",
	     type       => "treeNode",
	     attributes => [ "intent(inout)" ],
	     variables  => [ "self" ]
	 },
	 {
	     intrinsic  => "procedure",
	     type       => "",
	     attributes => [ "pointer" ],
	     variables  => [ "mapFunction" ]
	 }
	);
    my $functionCode;
    $functionCode .= "  subroutine mapComponentsVoid(self,mapFunction)\n";
    $functionCode .= "    !% Map a void function over components.\n";
    $functionCode .= "    implicit none\n";
    $functionCode .= &Format_Variable_Defintions(\@dataContent)."\n";
    foreach ( @{$buildData->{'componentClassList'}} ) {	    
     	$functionCode .= "    if (allocated(self%component".padComponentClass(ucfirst($_),[19,0]).")) call mapFunction(self%component".padComponentClass(ucfirst($_),[19,0]).")\n";
    }
    $functionCode .= "    return\n";
    $functionCode .= "  end subroutine mapComponentsVoid\n\n";
    # Insert into the function list.
    push(
	@{$buildData->{'code'}->{'functions'}},
	$functionCode
	);
    # Insert a type-binding for this function into the treeNode type.
    push(
	@{$buildData->{'types'}->{'treeNode'}->{'boundFunctions'}},
	{type => "procedure", name => "mapVoid", function => "mapComponentsVoid"}
	);

    # Function for mapping a scalar double function.
    @dataContent =
	(
	 {
	     intrinsic  => "class",
	     type       => "treeNode",
	     attributes => [ "intent(inout)" ],
	     variables  => [ "self" ]
	 },
	 {
	     intrinsic  => "procedure",
	     type       => "double precision",
	     attributes => [ "pointer" ],
	     variables  => [ "mapFunction" ]
	 },
	 {
	     intrinsic  => "integer",
	     attributes => [ "intent(in   )" ],
	     variables  => [ "reduction" ]
	 },
	 {
	     intrinsic  => "double precision",
	     variables  => [ "componentValue" ]
	 }
	);
    $functionCode  = "  double precision function mapComponentsDouble0(self,mapFunction,reduction)\n";
    $functionCode .= "    !% Map a scalar double function over components with a specified {\\tt reduction}.\n";
    $functionCode .= "    implicit none\n";
    $functionCode .= &Format_Variable_Defintions(\@dataContent)."\n";
    $functionCode .= "    select case (reduction)\n";
    $functionCode .= "    case (reductionSummation)\n";
    $functionCode .= "      mapComponentsDouble0=0.0d0\n";
    $functionCode .= "    case (reductionProduct  )\n";
    $functionCode .= "      mapComponentsDouble0=1.0d0\n";
    $functionCode .= "    end select\n";
    foreach ( @{$buildData->{'componentClassList'}} ) {	    
     	$functionCode .= "    if (allocated(self%component".padComponentClass(ucfirst($_),[0,0]).")) then\n";
     	$functionCode .= "      componentValue=mapFunction(self%component".padComponentClass(ucfirst($_),[0,0]).")\n";
     	$functionCode .= "      select case (reduction)\n";
     	$functionCode .= "      case (reductionSummation)\n";
     	$functionCode .= "        mapComponentsDouble0=mapComponentsDouble0+componentValue\n";
     	$functionCode .= "      case (reductionProduct  )\n";
     	$functionCode .= "        mapComponentsDouble0=mapComponentsDouble0*componentValue\n";
     	$functionCode .= "      end select\n";
     	$functionCode .= "    end if\n";
    }
    $functionCode .= "    return\n";
    $functionCode .= "  end function mapComponentsDouble0\n\n";
    # Insert into the function list.
    push(
	@{$buildData->{'code'}->{'functions'}},
	$functionCode
	);
    # Insert a type-binding for this function into the treeNode type.
    push(
	@{$buildData->{'types'}->{'treeNode'}->{'boundFunctions'}},
	{type => "procedure", name => "mapDouble0", function => "mapComponentsDouble0"}
	);
}

sub Generate_Node_Dump_Function {
    # Generate function to dump node properties.
    my $buildData = shift;

    # Create the function.
    my @dataContent =
	(
	 {
	     intrinsic  => "class",
	     type       => "treeNode",
	     attributes => [ "intent(in   )" ],
	     variables  => [ "self" ]
	 },
	 {
	     intrinsic  => "integer",
	     variables  => [ "i" ]
	 },
	 {
	     intrinsic  => "type",
	     type       => "varying_string",
	     variables  => [ "message" ]
	 }
	);
    my $functionCode;
    $functionCode  = "  subroutine Node_Dump(self)\n";
    $functionCode .= "    !% Dump node content.\n";
    $functionCode .= "    use ISO_Varying_String\n";
    $functionCode .= "    use Galacticus_Display\n";
    $functionCode .= "    use String_Handling\n";
    $functionCode .= "    implicit none\n";
    $functionCode .= &Format_Variable_Defintions(\@dataContent)."\n";
    # Iterate over pointers.
    $functionCode .= "    message='Dumping node '\n";
    $functionCode .= "    message=message//self%index()\n";
    $functionCode .= "    call Galacticus_Display_Indent(message)\n";
    $functionCode .= "    call Galacticus_Display_Indent('pointers')\n";
    foreach my $pointer ( "parent", "firstChild", "sibling", "firstSatellite", "mergeTarget", "firstMergee", "siblingMergee", "formationNode" ) {
	$functionCode .= "   message='".(" " x (14-length($pointer))).$pointer.": '\n";
	$functionCode .= "   message=message//self%".pad($pointer,14)."%index()\n";
	$functionCode .= "   call Galacticus_Display_Message(message)\n";
    }
    $functionCode .= "   call Galacticus_Display_Unindent('done')\n";
    # Iterate over all component classes
    foreach ( @{$buildData->{'componentClassList'}} ) {	    
	$functionCode .= "    if (allocated(self%component".padComponentClass(ucfirst($_),[0,0]).")) then\n";
	$functionCode .= "      do i=1,size(self%component".padComponentClass(ucfirst($_),[0,0]).")\n";
	$functionCode .= "        call self%component".padComponentClass(ucfirst($_),[0,0])."(i)%dump()\n";
	$functionCode .= "      end do\n";
	$functionCode .= "    end if\n";
    }
    $functionCode .= "    call Galacticus_Display_Unindent('done')\n";
    $functionCode .= "    return\n";
    $functionCode .= "  end subroutine Node_Dump\n";
    # Insert into the function list.
    push(
	@{$buildData->{'code'}->{'functions'}},
	$functionCode
	);
    # Insert a type-binding for this function into the treeNode type.
    push(
	@{$buildData->{'types'}->{'treeNode'}->{'boundFunctions'}},
	{type => "procedure", name => "dump", function => "Node_Dump"}
	);
}

sub Generate_Node_Output_Functions {
    # Generate functions to output node properties.
    my $buildData = shift;

    # Create an output count function.
    my @dataContent =
	(
	 {
	     intrinsic  => "class",
	     type       => "treeNode",
	     attributes => [ "intent(inout)" ],
	     variables  => [ "self" ]
	 },
	 {
	     intrinsic  => "integer",
	     attributes => [ "intent(inout)" ],
	     variables  => [ "integerPropertyCount", "doublePropertyCount" ]
	 },
	 {
	     intrinsic  => "double precision",
	     attributes => [ "intent(in   )" ],
	     variables  => [ "time" ]
	 },
	 {
	     intrinsic  => "integer",
	     variables  => [ "i" ]
	 }
	);
    my $functionCode;
    $functionCode  = "  subroutine Node_Output_Count(self,integerPropertyCount,doublePropertyCount,time)\n";
    $functionCode .= "    !% Increment the count of properties to output for this node.\n";
    $functionCode .= "    implicit none\n";
    $functionCode .= &Format_Variable_Defintions(\@dataContent)."\n";
    # Iterate over all component classes
    foreach ( @{$buildData->{'componentClassList'}} ) {	    
	$functionCode .= "    if (allocated(self%component".padComponentClass(ucfirst($_),[0,0]).")) then\n";
	$functionCode .= "      do i=1,size(self%component".padComponentClass(ucfirst($_),[0,0]).")\n";
	$functionCode .= "        call self%component".padComponentClass(ucfirst($_),[0,0])."(i)%outputCount(integerPropertyCount,doublePropertyCount,time,instance=i)\n";
	$functionCode .= "      end do\n";
	$functionCode .= "    end if\n";
    }
    $functionCode .= "    return\n";
    $functionCode .= "  end subroutine Node_Output_Count\n";
    # Insert into the function list.
    push(
	@{$buildData->{'code'}->{'functions'}},
	$functionCode
	);
    # Insert a type-binding for this function into the treeNode type.
    push(
	@{$buildData->{'types'}->{'treeNode'}->{'boundFunctions'}},
	{type => "procedure", name => "outputCount", function => "Node_Output_Count"}
	);
    # Create an output property names function.
    @dataContent =
	(
	 {
	     intrinsic  => "class",
	     type       => "treeNode",
	     attributes => [ "intent(inout)" ],
	     variables  => [ "self" ]
	 },
	 {
	     intrinsic  => "double precision",
	     attributes => [ "intent(in   )" ], 
	     variables  => [ "time" ]
	 },
	 {
	     intrinsic  => "integer", 
	     attributes => [ "intent(inout)" ], 
	     variables  => [ "integerProperty", "doubleProperty" ]
	 },
	 {
	     intrinsic  => "character",
	     type       => "len=*",
	     attributes => [ "intent(inout)", "dimension(:)" ], 
	     variables  => [ "integerPropertyNames", "integerPropertyComments", "doublePropertyNames", "doublePropertyComments" ]
	 },
	 {
	     intrinsic  => "double precision",
	     attributes => [ "intent(inout)", "dimension(:)" ],
	     variables  => [ "integerPropertyUnitsSI", "doublePropertyUnitsSI" ]
	 },
	 {
	     intrinsic  => "integer",
	     variables  => [ "i" ]
	 }
	);
    undef($functionCode);
    $functionCode  = "  subroutine Node_Output_Names(self,integerProperty,integerPropertyNames,integerPropertyComments,integerPropertyUnitsSI,doubleProperty,doublePropertyNames,doublePropertyComments,doublePropertyUnitsSI,time)\n";
    $functionCode .= "    !% Establish the names of properties to output for this node.\n";
    $functionCode .= "    implicit none\n";
    $functionCode .= &Format_Variable_Defintions(\@dataContent)."\n";
    # Iterate over all component classes
    foreach ( @{$buildData->{'componentClassList'}} ) {	    
	$functionCode .= "    if (allocated(self%component".padComponentClass(ucfirst($_),[0,0]).")) then\n";
	$functionCode .= "      do i=1,size(self%component".padComponentClass(ucfirst($_),[0,0]).")\n";
	$functionCode .= "        call self%component".padComponentClass(ucfirst($_),[0,0])."(i)%outputNames(integerProperty,integerPropertyNames,integerPropertyComments,integerPropertyUnitsSI,doubleProperty,doublePropertyNames,doublePropertyComments,doublePropertyUnitsSI,time,instance=i)\n";
	$functionCode .= "      end do\n";
	$functionCode .= "    end if\n";
    }
    $functionCode .= "    return\n";
    $functionCode .= "  end subroutine Node_Output_Names\n";
    # Insert into the function list.
    push(
	@{$buildData->{'code'}->{'functions'}},
	$functionCode
	);
    # Insert a type-binding for this function into the treeNode type.
    push(
	@{$buildData->{'types'}->{'treeNode'}->{'boundFunctions'}},
	{type => "procedure", name => "outputNames", function => "Node_Output_Names"}
	);
  # Create an output function.
    @dataContent =
	(
	 {
	     intrinsic  => "class",
	     type       => "treeNode",
	     attributes => [ "intent(inout)" ],
	     variables  => [ "self" ]
	 },
	 {
	     intrinsic  => "double precision",
	     attributes => [ "intent(in   )" ], 
	     variables  => [ "time" ]
	 },
	 {
	     intrinsic  => "integer", 
	     attributes => [ "intent(inout)" ], 
	     variables  => [ "integerProperty", "integerBufferCount", "doubleProperty", "doubleBufferCount" ]
	 },
	 {
	     intrinsic  => "integer",
	     type       => "kind=kind_int8",
	     attributes => [ "intent(inout)", "dimension(:,:)" ],
	     variables  => [ "integerBuffer" ]
	 },
	 {
	     intrinsic  => "double precision",
	     attributes => [ "intent(inout)", "dimension(:,:)" ],
	     variables  => [ "doubleBuffer" ]
	 },
	 {
	     intrinsic  => "integer",
	     variables  => [ "i" ]
	 }
	);
    undef($functionCode);
    $functionCode  = "  subroutine Node_Output(self,integerProperty,integerBufferCount,integerBuffer,doubleProperty,doubleBufferCount,doubleBuffer,time)\n";
    $functionCode .= "    ! Output properties for this node.\n";
    $functionCode .= "    implicit none\n";
    $functionCode .= &Format_Variable_Defintions(\@dataContent)."\n";
    # Iterate over all component classes
    foreach ( @{$buildData->{'componentClassList'}} ) {	    
	$functionCode .= "    if (allocated(self%component".padComponentClass(ucfirst($_),[0,0]).")) then\n";
	$functionCode .= "      do i=1,size(self%component".padComponentClass(ucfirst($_),[0,0]).")\n";
	$functionCode .= "        call self%component".padComponentClass(ucfirst($_),[0,0])."(i)%output(integerProperty,integerBufferCount,integerBuffer,doubleProperty&
       &,doubleBufferCount,doubleBuffer,time,instance=i)\n";
	$functionCode .= "      end do\n";
	$functionCode .= "    end if\n";
    }
    $functionCode .= "    return\n";
    $functionCode .= "  end subroutine Node_Output\n";
    # Insert into the function list.
    push(
	@{$buildData->{'code'}->{'functions'}},
	$functionCode
	);
    # Insert a type-binding for this function into the treeNode type.
    push(
	@{$buildData->{'types'}->{'treeNode'}->{'boundFunctions'}},
	{type => "procedure", name => "output", function => "Node_Output"}
	);
}

sub Generate_Node_Property_Name_From_Index_Function {
    # Generate function to get the name of a property given an index.
    my $buildData = shift;

    # Define variables.
    my @dataContent =
	(
	 {
	     intrinsic  => "class",
	     type       => "treeNode",
	     attributes => [ "intent(in   )" ],
	     variables  => [ "self" ]
	 },
	 {
	     intrinsic  => "integer",
	     attributes => [ "intent(in   )" ],
	     variables  => [ "index" ]
	 },
	 {
	     intrinsic  => "type",
	     type       => "varying_string",
	     variables  => [ "name" ]
	 },
	 {
	     intrinsic  => "integer",
	     variables  => [ "count", "i" ]
	 }
	);
    my $functionCode;
    $functionCode .= "  function Node_Property_Name_From_Index(self,index) result (name)\n";
    $functionCode .= "    !% Return the name of a property given its index.\n";
    $functionCode .= "    use ISO_Varying_String\n";
    $functionCode .= "    implicit none\n";
    $functionCode .= &Format_Variable_Defintions(\@dataContent)."\n";

    # Loop over all component classes
    $functionCode .= "  name='unknown'\n";
    $functionCode .= "  count=index\n";
    foreach ( @{$buildData->{'componentClassList'}} ) {	    
     	$functionCode .= "    if (allocated(self%component".padComponentClass(ucfirst($_),[0,0]).")) then\n";
	$functionCode .= "      do i=1,size(self%component".padComponentClass(ucfirst($_),[0,0]).")\n";
	if ( $workaround == 1 ) { # Workaround for http://gcc.gnu.org/bugzilla/show_bug.cgi?id=53876
	    $functionCode .= "        select type (component => self%component".padComponentClass(ucfirst($_),[0,0])."(i))\n";
	    foreach my $implementationName ( @{$buildData->{'componentClasses'}->{$_}->{'members'}} ) {
		$functionCode .= "      type is (nodeComponent".ucfirst($_).ucfirst($implementationName).")\n";
		$functionCode .= "        call Node_Component_".ucfirst($_).ucfirst($implementationName)."_Name_From_Index(component,count,name)\n";
	    }
	    $functionCode .= "        end select\n";
	} else {
	    $functionCode .= "        call self%component".padComponentClass(ucfirst($_),[0,0])."(i)%nameFromIndex(count,name)\n";
	}
	$functionCode .= "        if (count <= 0) return\n";
	$functionCode .= "      end do\n";
	$functionCode .= "    end if\n";
    }
    $functionCode .= "    if (name == 'unknown') call Galacticus_Error_Report('Node_Property_Name_From_Index','property index out of range')\n";
    $functionCode .= "    return\n";
    $functionCode .= "  end function Node_Property_Name_From_Index\n";
    # Insert into the function list.
    push(
	@{$buildData->{'code'}->{'functions'}},
	$functionCode
	);
    # Insert a type-binding for this function into the treeNode type.
    push(
	@{$buildData->{'types'}->{'treeNode'}->{'boundFunctions'}},
	{type => "procedure", name => "nameFromIndex", function => "Node_Property_Name_From_Index"}
	);
}

sub Generate_Node_Serialization_Functions {
    # Generate functions to serialize/deserialize nodes to/from arrays.
    my $buildData = shift;

    # Function computing a count of the serialization length.
    my @dataContent =
	(
	 {
	     intrinsic  => "class",
	     type       => "treeNode",
	     attributes => [ "intent(in   )" ],
	     variables  => [ "self" ]
	 },
	 {
	     intrinsic  => "integer",
	     variables  => [ "count", "i" ]
	 }
	);
    my $functionCode;
    $functionCode .= "  function SerializeToArrayCount(self) result (count)\n";
    $functionCode .= "    !% Return a count of the size of the serialized {\\tt treeNode} object.\n";
    $functionCode .= "    implicit none\n";
    $functionCode .= &Format_Variable_Defintions(\@dataContent)."\n";
    $functionCode .= "    count=0\n";
    # Loop over all component classes
    foreach ( @{$buildData->{'componentClassList'}} ) {	    
     	$functionCode .= "    if (allocated(self%component".padComponentClass(ucfirst($_),[0,0]).")) then\n";
	$functionCode .= "      do i=1,size(self%component".padComponentClass(ucfirst($_),[0,0]).")\n";
	if ( $workaround == 1 ) { # Workaround for http://gcc.gnu.org/bugzilla/show_bug.cgi?id=53876
	    $functionCode .= "        select type (component => self%component".padComponentClass(ucfirst($_),[0,0])."(i))\n";
	    foreach my $implementationName ( @{$buildData->{'componentClasses'}->{$_}->{'members'}} ) {
		$functionCode .= "      type is (nodeComponent".ucfirst($_).ucfirst($implementationName).")\n";
		$functionCode .= "      write (0,*) 'DEBUG -> SerializeToArrayCount -> nodeComponent".ucfirst($_).ucfirst($implementationName)."',i,Node_Component_".ucfirst($_).ucfirst($implementationName)."_Count(component)\n"
		    if ( $debugging == 1 );
		$functionCode .= "        count=count+Node_Component_".ucfirst($_).ucfirst($implementationName)."_Count(component)\n";
	    }
	    $functionCode .= "        end select\n";
	} else {
	    $functionCode .= "        count=count+self%component".padComponentClass(ucfirst($_),[0,0])."%serializeCount()\n";
	}
	$functionCode .= "      end do\n";
	$functionCode .= "    end if\n";
    }
    $functionCode .= "    return\n";
    $functionCode .= "  end function SerializeToArrayCount\n\n";
    # Insert into the function list.
    push(
	@{$buildData->{'code'}->{'functions'}},
	$functionCode
	);
    # Insert a type-binding for this function into the treeNode type.
    push(
	@{$buildData->{'types'}->{'treeNode'}->{'boundFunctions'}},
	{type => "procedure", name => "serializeCount", function => "serializeToArrayCount"}
	);

    # Iterate over all property-associated data for which we need serialization/deserialization functions.
    foreach my $content ( "value", "scale", "rate" ) {
	# Create the serialization function.
	@dataContent =
	    (
	     {
		 intrinsic  => "class",
		 type       => "treeNode",
		 attributes => [ "intent(in   )" ],
		 variables  => [ "self" ]
	     },
	     {
		 intrinsic  => "integer",
		 variables  => [ "count", "offset", "i" ],
	     },
	     {
		 intrinsic  => "double precision",
		 attributes => [ "dimension(:)", "intent(  out)" ],
		 variables  => [ "array" ]
	     }
	    );
	$functionCode  = "  subroutine SerializeToArray".pad(ucfirst($content)."s",5)."(self,array)\n";
	$functionCode .= "    !% Serialize ".$content."s to array.\n";
	$functionCode .= "    use Memory_Management\n";
	$functionCode .= "    implicit none\n";
	$functionCode .= &Format_Variable_Defintions(\@dataContent)."\n";
	$functionCode .= "    offset=1\n";
	# Loop over all component classes
	foreach ( @{$buildData->{'componentClassList'}} ) {	    
	    $functionCode .= "    if (allocated(self%component".padComponentClass(ucfirst($_),[0,0]).")) then\n";
	    $functionCode .= "      do i=1,size(self%component".padComponentClass(ucfirst($_),[0,0]).")\n";
	    if ( $workaround == 1 ) { # Workaround for http://gcc.gnu.org/bugzilla/show_bug.cgi?id=53876
		$functionCode .= "        count=0\n";
		$functionCode .= "        select type (component => self%component".padComponentClass(ucfirst($_),[0,0])."(i))\n";
		foreach my $implementationName ( @{$buildData->{'componentClasses'}->{$_}->{'members'}} ) {
		    $functionCode .= "      type is (nodeComponent".ucfirst($_).ucfirst($implementationName).")\n";
		    $functionCode .= "        count=Node_Component_".ucfirst($_).ucfirst($implementationName)."_Count(component)\n";
		    $functionCode .= "        write (0,*) 'DEBUG -> SerializeToArray".ucfirst($content)."s -> nodeComponent".ucfirst($_).ucfirst($implementationName)."',i,count,offset,size(array)\n"
			if ( $debugging == 1 );
		    $functionCode .= "        if (count > 0) call Node_Component_".ucfirst($_).ucfirst($implementationName)."_Serialize_".pad(ucfirst($content)."s",5)."(component,array(offset:))\n";
		    $functionCode .= "        if (count > 0 .and. any(array(offset:offset+count-1) <= 0.0d0)) write (0,*) 'DEBUG -> SerializeToArray".ucfirst($content)."s: non-positive scale found for ".$_." ".$implementationName."'\n"
			if ( $content eq "scale" && $debugging == 1 );
		}
		$functionCode .= "        end select\n";
	    } else {
		$functionCode .= "        count=self%component".padComponentClass(ucfirst($_),[0,0])."(i)%serializeCount()\n";
		$functionCode .= "        if (count > 0) call self%component".padComponentClass(ucfirst($_),[0,0])."(i)%serialize".pad(ucfirst($content)."s",5)."(array(offset:))\n";
	    }
	    $functionCode .= "        offset=offset+count\n";
	    $functionCode .= "      end do\n";
	    $functionCode .= "    end if\n";
	}
	$functionCode .= "    return\n";
	$functionCode .= "  end subroutine SerializeToArray".ucfirst($content)."s\n\n";
	# Insert into the function list.
	push(
	    @{$buildData->{'code'}->{'functions'}},
	    $functionCode
	    );
	# Insert a type-binding for this function into the treeNode type.
	push(
	    @{$buildData->{'types'}->{'treeNode'}->{'boundFunctions'}},
	    {type => "procedure", name => "serialize".ucfirst($content)."s", function => "serializeToArray".ucfirst($content)."s"}
	    );
	# Create the deserialization function.
	@dataContent =
	    (
	     {
		 intrinsic  => "class",
		 type       => "treeNode",
		 attributes => [ "intent(inout)" ],
		 variables  => [ "self" ]
	     },
	     {
		 intrinsic  => "integer",
		 variables  => [ "count", "offset", "i" ],
	     },
	     {
		 intrinsic  => "double precision",
		 attributes => [ "dimension(:)", "intent(in   )" ],
		 variables  => [ "array" ]
	     }
	    );
	$functionCode  = "  subroutine DeserializeFromArray".pad(ucfirst($content)."s",5)."(self,array)\n";
	$functionCode .= "    !% Deserialize ".$content."s from {\\tt array}.\n";
	$functionCode .= "    use Memory_Management\n";
	$functionCode .= "    implicit none\n";
	$functionCode .= &Format_Variable_Defintions(\@dataContent)."\n";
	$functionCode .= "    offset=1\n";
	# Loop over all component classes
	foreach ( @{$buildData->{'componentClassList'}} ) {	    
	    $functionCode .= "    if (allocated(self%component".padComponentClass(ucfirst($_),[0,0]).")) then\n";
	    $functionCode .= "      do i=1,size(self%component".padComponentClass(ucfirst($_),[0,0]).")\n";
	    if ( $workaround == 1 ) { # Workaround for http://gcc.gnu.org/bugzilla/show_bug.cgi?id=53876
		$functionCode .= "        count=0\n";
		$functionCode .= "        select type (component => self%component".padComponentClass(ucfirst($_),[0,0])."(i))\n";
		foreach my $implementationName ( @{$buildData->{'componentClasses'}->{$_}->{'members'}} ) {
		    $functionCode .= "      type is (nodeComponent".ucfirst($_).ucfirst($implementationName).")\n";
		    $functionCode .= "        count=Node_Component_".ucfirst($_).ucfirst($implementationName)."_Count(component)\n";
		    $functionCode .= "        if (count > 0) call Node_Component_".ucfirst($_).ucfirst($implementationName)."_Deserialize_".pad(ucfirst($content)."s",5)."(component,array(offset:))\n";
		}
		$functionCode .= "        end select\n";
	    } else {		
		$functionCode .= "        count=self%component".padComponentClass(ucfirst($_),[0,0])."(i)%serializeCount()\n";
		$functionCode .= "        if (count > 0) call self%component".padComponentClass(ucfirst($_),[0,0])."(i)%deserialize".pad(ucfirst($content)."s",5)."(array(offset:))\n";
	    }
	    $functionCode .= "        offset=offset+count\n";
	    $functionCode .= "      end do\n";
	    $functionCode .= "    end if\n";
	}
	$functionCode .= "    return\n";
	$functionCode .= "  end subroutine DeserializeFromArray".ucfirst($content)."s\n\n";
	# Insert into the function list.
	push(
	    @{$buildData->{'code'}->{'functions'}},
	    $functionCode
	    );
	# Insert a type-binding for this function into the treeNode type.
	push(
	    @{$buildData->{'types'}->{'treeNode'}->{'boundFunctions'}},
	    {type => "procedure", name => "deserialize".ucfirst($content)."s", function => "deserializeFromArray".ucfirst($content)."s"}
	    );
    }
}

sub Generate_Node_ODE_Initialization_Functions {
    # Generate functions initialize a node for an ODE step.
    my $buildData = shift;
    # Create functions to initialize property rates for an ODE step.
    my @dataContent =
	(
	 {
	     intrinsic  => "class",
	     type       => "treeNode",
	     attributes => [ "intent(inout)" ],
	     variables  => [ "self" ]
	 },
	 {
	     intrinsic  => "integer",
	     variables  => [ "i" ]
	 }
	);
    my $functionCode;
    $functionCode  = "  subroutine Tree_Node_ODE_Step_Rates_Initialize(self)\n";
    $functionCode .= "    !% Initialize the rates in components of tree node {\\tt self} in preparation for an ODE solver step.\n";
    $functionCode .= "    implicit none\n";
    $functionCode .= &Format_Variable_Defintions(\@dataContent)."\n";
    # Loop over all component classes
    foreach ( @{$buildData->{'componentClassList'}} ) {	    
     	$functionCode .= "    if (allocated(self%component".padComponentClass(ucfirst($_),[0,0]).")) then\n";
	$functionCode .= "      do i=1,size(self%component".padComponentClass(ucfirst($_),[0,0]).")\n";
	if ( $workaround == 1 ) { # Workaround for http://gcc.gnu.org/bugzilla/show_bug.cgi?id=53876
	    $functionCode .= "        select type (component => self%component".padComponentClass(ucfirst($_),[0,0])."(i))\n";
	    foreach my $implementationName ( @{$buildData->{'componentClasses'}->{$_}->{'members'}} ) {
		$functionCode .= "      type is (nodeComponent".ucfirst($_).ucfirst($implementationName).")\n";
		$functionCode .= "        call Node_Component_".ucfirst($_).ucfirst($implementationName)."_ODE_Step_Rates_Init(component)\n";
	    }
	    $functionCode .= "    end select\n";
	} else {
	    $functionCode .= "        call self%component".padComponentClass(ucfirst($_),[0,0])."(i)%odeStepRatesInitialize()\n";
	}
	$functionCode .= "      end do\n";
	$functionCode .= "    end if\n";	
    }
    $functionCode .= "    return\n";
    $functionCode .= "  end subroutine Tree_Node_ODE_Step_Rates_Initialize\n\n";
    # Insert into the function list.
    push(
	@{$buildData->{'code'}->{'functions'}},
	$functionCode
	);
    # Insert a type-binding for this function into the treeNode type.
    push(
	@{$buildData->{'types'}->{'treeNode'}->{'boundFunctions'}},
	{type => "procedure", name => "odeStepRatesInitialize", function => "Tree_Node_ODE_Step_Rates_Initialize"},
	);    
    # Create functions to initialize property scales for an ODE step.
    $functionCode  = "  subroutine Tree_Node_ODE_Step_Scales_Initialize(self)\n";
    $functionCode .= "    !% Initialize the scales in components of tree node {\\tt self} in preparation for an ODE solver step.\n";
    $functionCode .= "    implicit none\n";
    $functionCode .= &Format_Variable_Defintions(\@dataContent)."\n";
    # Loop over all component classes
    foreach ( @{$buildData->{'componentClassList'}} ) {	    
     	$functionCode .= "    if (allocated(self%component".padComponentClass(ucfirst($_),[0,0]).")) then\n";
	$functionCode .= "      do i=1,size(self%component".padComponentClass(ucfirst($_),[0,0]).")\n";
	if ( $workaround == 1 ) { # Workaround for http://gcc.gnu.org/bugzilla/show_bug.cgi?id=53876
	    $functionCode .= "        select type (component => self%component".padComponentClass(ucfirst($_),[0,0])."(i))\n";
	    foreach my $implementationName ( @{$buildData->{'componentClasses'}->{$_}->{'members'}} ) {
		$functionCode .= "      type is (nodeComponent".ucfirst($_).ucfirst($implementationName).")\n";
		$functionCode .= "        call Node_Component_".ucfirst($_).ucfirst($implementationName)."_ODE_Step_Scales_Init(component)\n";
	    }
	    $functionCode .= "    end select\n";
	} else {
	    $functionCode .= "        call self%component".padComponentClass(ucfirst($_),[0,0])."(i)%odeStepScalesInitialize()\n";
	}
	$functionCode .= "      end do\n";
	$functionCode .= "    end if\n";
    }
    $functionCode .= "    return\n";
    $functionCode .= "  end subroutine Tree_Node_ODE_Step_Scales_Initialize\n\n";
    # Insert into the function list.
    push(
	@{$buildData->{'code'}->{'functions'}},
	$functionCode
	);
    # Insert a type-binding for this function into the treeNode type.
    push(
	@{$buildData->{'types'}->{'treeNode'}->{'boundFunctions'}},
	{type => "procedure", name => "odeStepScalesInitialize" , function => "Tree_Node_ODE_Step_Scales_Initialize"}
	);
}

sub Generate_Implementation_Dump_Functions {
    # Generate dump for each component implementation.
    my $buildData = shift;
    # Iterate over component implementations.
    foreach my $componentID ( @{$buildData->{'componentIdList'}} ) {
	# Get the component.
	my $component = $buildData->{'components'}->{$componentID};
	# Initialize function code.
	my $functionCode;
	# Initialize data content.
	my @dataContent =
	    (
	     {
		 intrinsic  => "class",
		 type       => "nodeComponent".ucfirst($componentID),
		 attributes => [ "intent(in   )" ],
		 variables  => [ "self" ]
	     },
	     {
		 intrinsic  => "type",
		 type       => "varying_string",
		 variables  => [ "message" ]
	     },
	     {
		 intrinsic  => "character",
		 type       => "len=12",
		 variables  => [ "label" ]
	     },
	     {
		 intrinsic  => "integer",
		 variables  => [ "i" ]
	     }
	    );
	# Format labels for different data types.
	my %formatLabel = 
	    (
	     real    => "'(e12.6)'",
	     integer => "'(i)'",
	     logical => "*"
	    );
	# Generate dump function.
	$functionCode  = "  subroutine Node_Component_".ucfirst($componentID)."_Dump(self)\n";
	$functionCode .= "    !% Dump the contents of a ".$component->{'name'}." implementation of the ".$component->{'class'}." component.\n";
	$functionCode .= "    use Galacticus_Display\n";
	$functionCode .= "    use ISO_Varying_String\n";
	$functionCode .= "    use String_Handling\n";
	$functionCode .= "    implicit none\n";
	$functionCode .= &Format_Variable_Defintions(\@dataContent)."\n";
	unless ( $component->{'name'} eq "null" ) {
	    # Dump the parent type if necessary.
	    $functionCode .= "    call self%nodeComponent".ucfirst($component->{'extends'}->{'class'}).ucfirst($component->{'extends'}->{'name'})."%dump()\n"
		if ( exists($component->{'extends'}) );
	    $functionCode .= "    call Galacticus_Display_Indent('".$component->{'class'}.": ".(" " x ($fullyQualifiedNameLengthMax-length($component->{'class'}))).$component->{'name'}."')\n";
	    foreach my $methodName ( keys(%{$component->{'methods'}->{'method'}}) ) {
		my $method = $component->{'methods'}->{'method'}->{$methodName};
		# Check if this method has any linked data in this component.
		if ( exists($method->{'linkedData'}) ) {
		    my $linkedDataName = $method->{'linkedData'};
		    my $linkedData     = $component->{'content'}->{'data'}->{$linkedDataName};
		    if ( $linkedData->{'rank'} == 0 ) {
			switch ( $linkedData->{'type'} ) {
			    case ( [ "real", "integer", "logical" ] ) {
				$functionCode .= "    write (label,".$formatLabel{$linkedData->{'type'}}.") self%".padLinkedData($linkedDataName,[0,0])."%value\n";
				$functionCode .= "    message='".$methodName.": ".(" " x ($implementationPropertyNameLengthMax-length($methodName)))."'//label\n";
				$functionCode .= "    call Galacticus_Display_Message(message)\n";
			    }
			    else {
				$functionCode .= "    message='".$methodName.":'\n";
				$functionCode .= "    call Galacticus_Display_Indent(message)\n";
				$functionCode .= "    call self%".padLinkedData($linkedDataName,[0,0])."%value%dump()\n";
				$functionCode .= "    call Galacticus_Display_Unindent('end')\n";
				}
			}
		    } elsif ( $linkedData->{'rank'} == 1 ) {
			switch ( $linkedData->{'type'} ) {
			    case ( [ "real", "integer", "logical" ] ) {
				$functionCode .= "    do i=1,size(self%".$linkedDataName."%value)\n";
				$functionCode .= "       write (label,'(i3)') i\n";
				$functionCode .= "       message='".$methodName.": ".(" " x ($implementationPropertyNameLengthMax-length($methodName)))." '//trim(label)\n";
				$functionCode .= "       write (label,".$formatLabel{$linkedData->{'type'}}.") self%".$linkedDataName."%value(i)\n";
				$functionCode .= "       message=message//label\n";
				$functionCode .= "       call Galacticus_Display_Message(message)\n";
				$functionCode .= "    end do\n";
			    }
			    else {
				$functionCode .= "    do i=1,size(self%".$linkedDataName."%value)\n";
				$functionCode .= "       write (label,'(i3)') i\n";
				$functionCode .= "       message='".$methodName.": ".(" " x ($implementationPropertyNameLengthMax-length($methodName)))." '//trim(label)\n";
				$functionCode .= "       call Galacticus_Display_Indent(message)\n";
				$functionCode .= "       call self%".$linkedDataName."%value(i)\n";
				$functionCode .= "       call Galacticus_Display_Unindent('end')\n";
				$functionCode .= "    end do\n";
			    }
			}			
		    }
		}
	    }
	    $functionCode .= "    call Galacticus_Display_Unindent('done')\n";
	}
	$functionCode .= "    return\n";
	$functionCode .= "  end subroutine Node_Component_".ucfirst($componentID)."_Dump\n";
	# Insert into the function list.
	push(
	    @{$buildData->{'code'}->{'functions'}},
	    $functionCode
	    );
	# Insert a type-binding for this function into the implementation type.
	push(
	    @{$buildData->{'types'}->{'nodeComponent'.ucfirst($componentID)}->{'boundFunctions'}},
	    {type => "procedure", name => "dump", function => "Node_Component_".ucfirst($componentID)."_Dump"},
	    );
    }
}

sub Generate_Implementation_Output_Functions {
    # Generate output functions for each component implementation.
    my $buildData = shift;
    # Iterate over component implementations.
    foreach my $componentID ( @{$buildData->{'componentIdList'}} ) {
	# Get the component.
	my $component = $buildData->{'components'}->{$componentID};
	# Find modules required.
	my %modulesRequired;
	foreach my $methodName ( keys(%{$component->{'methods'}->{'method'}}) ) {
	    my $method = $component->{'methods'}->{'method'}->{$methodName};
	    # Check if this method is to be output.
	    if ( exists($method->{'output'}) ) {
		if ( exists($method->{'output'}->{'modules'}) ) {
		    my $moduleList = $method->{'output'}->{'modules'};
		    $moduleList =~ s/^\s*//;
		    $moduleList =~ s/\s*$//;
		    my @modules = split(/\s*,\s*/,$moduleList);
		    foreach ( @modules ) {
			$modulesRequired{$_} = 1;
		    }
		}
	    }
	}
	# Initialize function code.
	my $functionCode;
	# Create property count function.
	my @dataContent =
	    (
	     {
		 intrinsic  => "class",
		 type       => "nodeComponent".ucfirst($componentID),
		 attributes => [ "intent(inout)" ],
		 variables  => [ "self" ]
	     },
	     {
		 intrinsic  => "integer",
		 attributes => ["intent(inout)" ],
		 variables  => [ "integerPropertyCount", "doublePropertyCount" ]
	     },
	     {
		 intrinsic  => "double precision",
		 attributes => [ "intent(in   )" ],
		 variables  => [ "time" ]
	     },
	     {
		 intrinsic  => "integer",
		 attributes => ["intent(in   )" ],
		 variables  => [ "instance" ]
	     },
	     {
		 intrinsic  => "integer",
		 variables  => [ "i" ]
	     }
	    );
	# Find all derived types to be output.
	my %outputTypes;
	foreach my $methodName ( keys(%{$component->{'methods'}->{'method'}}) ) {
	    my $method = $component->{'methods'}->{'method'}->{$methodName};
	    # Check if property is to be output.
	    if ( exists($method->{'output'}) ) {
		# Get the type of this component.
		my $type;
		if ( exists($method->{'linkedData'}) ) {
		    my $linkedDataName = $method->{'linkedData'};
		    my $linkedData     = $component->{'content'}->{'data'}->{$linkedDataName};
		    $type   = $linkedData->{'type'};
		} else {
		    $type = $method->{'type'};		
		}
		$outputTypes{$type} = 1
		    unless ( $type eq "real" || $type eq "integer" );
	    }
	}
	my @outputTypes;
	foreach ( keys(%outputTypes) ){
	    push(
		@outputTypes,
		{
		    intrinsic => "type",
		    type      => $_,
		    variables => [ "output".ucfirst($_) ]
		}
		);
	}
	push(@dataContent,@outputTypes);
	# Generate output count function.
	undef($functionCode);
	$functionCode  = "  subroutine Node_Component_".ucfirst($componentID)."_Output_Count(self,integerPropertyCount,doublePropertyCount,time,instance)\n";
	$functionCode .= "    !% Increment output property count for a ".$component->{'name'}." implementation of the ".$component->{'class'}." component.\n";
	$functionCode .= "    use ".$_."\n"
	    foreach ( keys(%modulesRequired) );
	$functionCode .= "    implicit none\n";
	$functionCode .= &Format_Variable_Defintions(\@dataContent)."\n";
	unless ( $component->{'name'} eq "null" ) {
	    # Act on the parent type if necessary.
	    $functionCode .= "    call self%nodeComponent".ucfirst($component->{'extends'}->{'class'}).ucfirst($component->{'extends'}->{'name'})."%outputCount(integerPropertyCount,doublePropertyCount,time,instance)\n"
		if ( exists($component->{'extends'}) );
	    # Check that this instance is to be output.
	    my $checkAdded = 0;
	    if (
		exists($component->{'output'}               )           &&
		exists($component->{'output'}->{'instances'})           &&
		       $component->{'output'}->{'instances'} eq "first"
		) {
		$functionCode .= "    if (instance == 1) then\n";
		$checkAdded = 1;
	    }
	    # Initialize counts.
	    my %typeCount =
		(
		 real    => 0,
		 integer => 0
		);
	    my %typeMap =
		(
		 real    => "double",
		 integer => "integer"
		);
	    foreach my $methodName ( keys(%{$component->{'methods'}->{'method'}}) ) {
		my $method = $component->{'methods'}->{'method'}->{$methodName};
		# Check if this method is to be output.
		if ( exists($method->{'output'}) ) {
		    # Define rank, type and value.
		    my $rank;
		    my $type;
		    # Check if this method has any linked data in this component.
		    if ( exists($method->{'linkedData'}) ) {
			my $linkedDataName = $method->{'linkedData'};
			my $linkedData     = $component->{'content'}->{'data'}->{$linkedDataName};
			$rank   = $linkedData->{'rank'};
			$type   = $linkedData->{'type'};
		    } elsif ( $method->{'isVirtual'} eq "yes" && $method->{'attributes'}->{'isGettable'} eq "true" ) {
			$rank = $method->{'rank'};
			$type = $method->{'type'};
		    } else {
			die("Generate_Implementation_Output_Functions(): can not output [".$methodName."]");
		    }
		    # Determine count.
		    my $count;
		    if ( $rank == 0 ) {
			$count = 1;
		    } elsif ( $rank ==1 ) {
			die("Generate_Implementation_Output_Functions(): output of rank>0 objects requires a labels attribute")
			    unless ( exists($method->{'output'}->{'labels'}) );	
			if ( $method->{'output'}->{'labels'} =~ m/^\[(.*)\]$/ ) {
			    my $labelText = $1;
			    $labelText    =~ s/^\s*//;
			    $labelText    =~ s/\s*$//;
			    my @labels    = split(/\s*,\s*/,$labelText);
			    $count = scalar(@labels);
			} elsif ( exists($method->{'output'}->{'count'}) ) {
			    $count = $method->{'output'}->{'count'};
			} else {
			    die('Generate_Implementation_Output_Functions(): no method to determine output count of property');
			}
		    } else {
			die("Generate_Implementation_Output_Functions(): output of rank>1 arrays not supported");
		    }
		    # Increment the counters.
		    switch ( $type ) {
			case ( [ "real", "integer" ] ) {
			    if ( $rank == 0 ) {
				if ( exists($method->{'output'}->{'condition'}) ) {
				    my $condition = $method->{'output'}->{'condition'};
				    $condition =~ s/\[\[([^\]]+)\]\]/$1/g;
				    $functionCode .= "    if (".$condition.") ".$typeMap{$type}."PropertyCount=".$typeMap{$type}."PropertyCount+".$count."\n";
				} elsif ( $count =~ m/^\d/ ) {
				    $typeCount{$type} += $count;
				} else {
				    $functionCode .= "    ".$typeMap{$type}."PropertyCount=".$typeMap{$type}."PropertyCount+".$count."\n";
				}
			    } elsif ( $rank == 1 ) {
				if ( exists($method->{'output'}->{'condition'}) ) {
				    my $condition = $method->{'output'}->{'condition'};
				    $condition =~ s/\[\[([^\]]+)\]\]/$1/g;
				    $condition =~ s/\{i\}/i/g;
				    $functionCode .= "    do i=1,".$count."\n";
				    $functionCode .= "    if (".$condition.") ".$typeMap{$type}."PropertyCount=".$typeMap{$type}."PropertyCount+1\n";
				    $functionCode .= "    end do\n";
				} else {
				    if ( $count =~ m/^\d/ ) {
					$typeCount{$type} += $count;
				    } else {
					$functionCode .= "    ".$typeMap{$type}."PropertyCount=".$typeMap{$type}."PropertyCount+".$count."\n";
				    }  
				}
			    }
			}
			else {
			    $functionCode .= "    output".ucfirst($type)."=self%".$methodName."()\n";
			    $functionCode .= "    call output".ucfirst($type)."%outputCount(integerPropertyCount,doublePropertyCount,time)\n";
			}
		    }
		}
	    }
	    $functionCode .= "    doublePropertyCount =doublePropertyCount +".$typeCount{'real'   }."\n"
		unless ( $typeCount{'real'   } == 0 );
	    $functionCode .= "    integerPropertyCount=integerPropertyCount+".$typeCount{'integer'}."\n"
		unless ( $typeCount{'integer'} == 0 );
	    $functionCode .= "    end if\n"
		if ( $checkAdded == 1 );
	}
	$functionCode .= "    return\n";
	$functionCode .= "  end subroutine Node_Component_".ucfirst($componentID)."_Output_Count\n";
	# Insert into the function list.
	push(
	    @{$buildData->{'code'}->{'functions'}},
	    $functionCode
	    );
	# Insert a type-binding for this function into the implementation type.
	push(
	    @{$buildData->{'types'}->{'nodeComponent'.ucfirst($componentID)}->{'boundFunctions'}},
	    {type => "procedure", name => "outputCount", function => "Node_Component_".ucfirst($componentID)."_Output_Count"},
	    );
	# Create property names function.
	@dataContent =
	    (
	     {
		 intrinsic  => "class",
		 type       => "nodeComponent".ucfirst($componentID),
		 attributes => [ "intent(inout)" ],
		 variables  => [ "self" ]
	     },
	     {
		 intrinsic  => "double precision",
		 attributes => [ "intent(in   )" ], 
		 variables  => [ "time" ]
	     },
	     {
		 intrinsic  => "integer", 
		 attributes => [ "intent(inout)" ], 
		 variables  => [ "integerProperty", "doubleProperty" ]
	     },
	     {
		 intrinsic  => "character",
		 type       => "len=*",
		 attributes => [ "intent(inout)", "dimension(:)" ], 
		 variables  => [ "integerPropertyNames", "integerPropertyComments", "doublePropertyNames", "doublePropertyComments" ]
	     },
	     {
		 intrinsic  => "double precision",
		 attributes => [ "intent(inout)", "dimension(:)" ],
		 variables  => [ "integerPropertyUnitsSI", "doublePropertyUnitsSI" ]
	     },
	     {
		 intrinsic  => "integer",
		 attributes => ["intent(in   )" ],
		 variables  => [ "instance" ]
	     },
	     {
		 intrinsic  => "integer",
		 variables  => [ "i" ]
	     }
	    );
	push(@dataContent,@outputTypes);
	undef($functionCode);
	$functionCode  = "  subroutine Node_Component_".ucfirst($componentID)."_Output_Names(self,integerProperty,integerPropertyNames,integerPropertyComments,integerPropertyUnitsSI,doubleProperty,doublePropertyNames,doublePropertyComments,doublePropertyUnitsSI,time,instance)\n";
	$functionCode .= "    !% Establish property names for a ".$component->{'name'}." implementation of the ".$component->{'class'}." component.\n";
	$functionCode .= "    use ".$_."\n"
	    foreach ( keys(%modulesRequired) );
	$functionCode .= "    implicit none\n";
	$functionCode .= &Format_Variable_Defintions(\@dataContent)."\n";
	unless ( $component->{'name'} eq "null" ) {
	    # Act on the parent type if necessary.
	    $functionCode .= "    call self%nodeComponent".ucfirst($component->{'extends'}->{'class'}).ucfirst($component->{'extends'}->{'name'})."%outputNames(integerProperty,integerPropertyNames,integerPropertyComments,integerPropertyUnitsSI,doubleProperty,doublePropertyNames,doublePropertyComments,doublePropertyUnitsSI,time,instance)\n"
		if ( exists($component->{'extends'}) );
	    # Check that this instance is to be output.
	    my $checkAdded = 0;
	    if (
		exists($component->{'output'}               )           &&
		exists($component->{'output'}->{'instances'})           &&
		       $component->{'output'}->{'instances'} eq "first"
		) {
		$functionCode .= "    if (instance == 1) then\n";
		$checkAdded = 1;
	    }
	    my %typeMap =
		(
		 real    => "double",
		 integer => "integer"
		);
	    foreach my $methodName ( keys(%{$component->{'methods'}->{'method'}}) ) {
		my $method = $component->{'methods'}->{'method'}->{$methodName};
		# Check if this method is to be output.
		if ( exists($method->{'output'}) ) {
		    # Define rank, type and value.
		    my $rank;
		    my $type;
		    my $object;
		    # Check if this method has any linked data in this component.
		    if ( exists($method->{'linkedData'}) ) {
			my $linkedDataName = $method->{'linkedData'};
			my $linkedData     = $component->{'content'}->{'data'}->{$linkedDataName};
			$rank   = $linkedData->{'rank'};
			$type   = $linkedData->{'type'};
			$object = "self%".$linkedDataName;
		    } elsif ( $method->{'isVirtual'} eq "yes" && $method->{'attributes'}->{'isGettable'} eq "true" ) {
			$rank = $method->{'rank'};
			$type = $method->{'type'};
		    } else {
			die("Generate_Implementation_Output_Functions(): can not output [".$methodName."]");
		    }		   
		    # Increment the counters.
		    switch ( $type ) {
			case ( [ "real", "integer" ] ) {
			    if ( $rank == 0 ) {
				if ( exists($method->{'output'}->{'condition'}) ) {
				    my $condition = $method->{'output'}->{'condition'};
				    $condition =~ s/\[\[([^\]]+)\]\]/$1/g;
				    $functionCode .= "    if (".$condition.") then\n";
				}
				$functionCode .= "       ".$typeMap{$type}."Property=".$typeMap{$type}."Property+1\n";
				$functionCode .= "       ".$typeMap{$type}."PropertyNames   (".$typeMap{$type}."Property)='".$component->{'class'}.ucfirst($methodName)."'\n";
				$functionCode .= "       ".$typeMap{$type}."PropertyComments(".$typeMap{$type}."Property)='".$method->{'output'}->{'comment'  }."'\n";
				$functionCode .= "       ".$typeMap{$type}."PropertyUnitsSI (".$typeMap{$type}."Property)=".$method->{'output'}->{'unitsInSI'}."\n";
				$functionCode .= "    end if\n"
				    if ( exists($method->{'output'}->{'condition'}) );
			    } elsif ( $rank == 1 ) {
				die("Generate_Implementation_Output_Functions(): output of rank>0 objects requires a labels attribute")
				    unless ( exists($method->{'output'}->{'labels'}) );

				if ( $method->{'output'}->{'labels'} =~ m/^\[(.*)\]$/ ) {
				    if ( exists($method->{'output'}->{'condition'}) ) {
					my $condition = $method->{'output'}->{'condition'};
					$condition =~ s/\[\[([^\]]+)\]\]/$1/g;
					$functionCode .= "    if (".$condition.") then\n";
				    }
				    my $labelText = $1;
				    $labelText =~ s/^\s*//;
				    $labelText =~ s/\s*$//;
				    my @labels    = split(/\s*,\s*/,$labelText);
				    foreach my $label ( @labels ) {
					$functionCode .= "       ".$typeMap{$type}."Property=".$typeMap{$type}."Property+1\n";
					$functionCode .= "       ".$typeMap{$type}."PropertyNames   (".$typeMap{$type}."Property)='".$component->{'class'}.ucfirst($methodName).$label."'\n";
					$functionCode .= "       ".$typeMap{$type}."PropertyComments(".$typeMap{$type}."Property)='".$method->{'output'}->{'comment'  }." [".$label."]'\n";
					$functionCode .= "       ".$typeMap{$type}."PropertyUnitsSI (".$typeMap{$type}."Property)=".$method->{'output'}->{'unitsInSI'}."\n";	
				    }
				    $functionCode .= "    end if\n"
					if ( exists($method->{'output'}->{'condition'}) );
				} elsif ( exists($method->{'output'}->{'count'}) ) {
				    my $count = $method->{'output'}->{'count' };
				    my $label = $method->{'output'}->{'labels'};
				    $label =~ s/\{i\}/i/g;
				    $functionCode .= "       do i=1,".$count."\n";
				    if ( exists($method->{'output'}->{'condition'}) ) {
					my $condition = $method->{'output'}->{'condition'};
					$condition =~ s/\[\[([^\]]+)\]\]/$1/g;
					$condition =~ s/\{i\}/i/g;
					$functionCode .= "    if (".$condition.") then\n";
				    }
				    $functionCode .= "         ".$typeMap{$type}."Property=".$typeMap{$type}."Property+1\n";
				    $functionCode .= "         ".$typeMap{$type}."PropertyNames   (".$typeMap{$type}."Property)='".$component->{'class'}.ucfirst($methodName)."'//".$label."\n";
				    $functionCode .= "         ".$typeMap{$type}."PropertyComments(".$typeMap{$type}."Property)='".$method->{'output'}->{'comment'  }." ['//".$label."//']'\n";
				    $functionCode .= "         ".$typeMap{$type}."PropertyUnitsSI (".$typeMap{$type}."Property)=".$method->{'output'}->{'unitsInSI'}."\n";	
				    $functionCode .= "    end if\n"
					if ( exists($method->{'output'}->{'condition'}) );
				    $functionCode .= "end do\n";
				} else {
				    die('Generate_Implementation_Output_Functions(): no method to determine output count of property');
				}
			    } else {
				die("Generate_Implementation_Output_Functions(): can not output rank>1 properties");
			    }
			}
			else {
			    $functionCode .= "    output".ucfirst($type)."=self%".$methodName."()\n";
			    $functionCode .= "    call output".ucfirst($type)."%outputNames(integerProperty,integerPropertyNames,integerPropertyComments,integerPropertyUnitsSI,doubleProperty,doublePropertyNames,doublePropertyComments,doublePropertyUnitsSI,time,'".$component->{'class'}.ucfirst($methodName)."','".$method->{'output'}->{'comment'}."',".$method->{'output'}->{'unitsInSI'}.")\n";
			}
		    }
		}
	    }
	    $functionCode .= "    end if\n"
		if ( $checkAdded == 1 );
	}
	$functionCode .= "    return\n";
	$functionCode .= "  end subroutine Node_Component_".ucfirst($componentID)."_Output_Names\n";
	# Insert into the function list.
	push(
	    @{$buildData->{'code'}->{'functions'}},
	    $functionCode
	    );
	# Insert a type-binding for this function into the implementation type.
	push(
	    @{$buildData->{'types'}->{'nodeComponent'.ucfirst($componentID)}->{'boundFunctions'}},
	    {type => "procedure", name => "outputNames", function => "Node_Component_".ucfirst($componentID)."_Output_Names"},
	    );
    }
}

sub Generate_Implementation_Name_From_Index_Functions {
    # Generate serialization/deserialization functions for each component implementation.
    my $buildData = shift;
    # Initialize function code.
    my $functionCode;
    # Initialize data content.
    my @dataContent;
    # Iterate over component implementations.
    foreach my $componentID ( @{$buildData->{'componentIdList'}} ) {
	# Get the component.
	my $component = $buildData->{'components'}->{$componentID};
	# Generate data content.
	@dataContent =
	    (
	     {
		 intrinsic  => "class",
		 type       => "nodeComponent".ucfirst($componentID),
		 attributes => [ "intent(in   )" ],
		 variables  => [ "self" ]
	     },
	     {
		 intrinsic  => "integer",
		 attributes => [ "intent(inout)" ],
		 variables  => [ "count" ]
	     },
	     {
		 intrinsic  => "type",
		 type       => "varying_string",	
		 attributes => [ "intent(inout)" ],
		 variables  => [ "name" ]
	     }
	    );
	# Generate the function.
  	$functionCode  = "  subroutine Node_Component_".ucfirst($componentID)."_Name_From_Index(self,count,name)\n";
	$functionCode .= "    !% Return the name of the property of given index for a ".$component->{'name'}." implementation of the ".$component->{'class'}." component.\n";
	$functionCode .= "    use ISO_Varying_String\n";
	$functionCode .= "    implicit none\n";
	$functionCode .= &Format_Variable_Defintions(\@dataContent)."\n";
	# Iterate over properties.
	foreach my $methodName ( keys(%{$component->{'methods'}->{'method'}}) ) {
	    my $method = $component->{'methods'}->{'method'}->{$methodName};
   	    # Check if this method has any linked data in this component.
	    if ( exists($method->{'linkedData'}) ) {
		# For each linked datum count if necessary.
		my $linkedDataName = $method->{'linkedData'};
		my $linkedData     = $component->{'content'}->{'data'}->{$linkedDataName};
		if ( $linkedData->{'isEvolvable'} eq "true" ) {
		    if ( $linkedData->{'rank'} == 0 ) {
			switch ( $linkedData->{'type'} ) {
			    case ( "real" ) {
				$functionCode .= "    count=count-1\n";
			    }
			    else {
				$functionCode .= "    count=count-self%".padLinkedData($linkedDataName,[0,0])."%value%serializeCount()\n";
			    }
			}
		    } else {
			$functionCode .= "    if (allocated(self%".padLinkedData($linkedDataName,[0,0])."%value)) count=count-size(self%".padLinkedData($linkedDataName,[0,0])."%value)\n";
		    }
		    $functionCode .= "    if (count <= 0) then\n";
		    $functionCode .= "      name='".$component->{'class'}.":".$component->{'name'}.":".$methodName."'\n";
		    $functionCode .= "      return\n";
		    $functionCode .= "    end if\n";
		}
	    }
	}
	$functionCode .= "    return\n";
	$functionCode .= "  end subroutine Node_Component_".ucfirst($componentID)."_Name_From_Index\n\n";
	# Insert into the function list.
	push(
	    @{$buildData->{'code'}->{'functions'}},
	    $functionCode
	    );
	# Insert a type-binding for this function into the implementation type.
	push(
	    @{$buildData->{'types'}->{'nodeComponent'.ucfirst($componentID)}->{'boundFunctions'}},
	    {type => "procedure", name => "nameFromIndex", function => "Node_Component_".ucfirst($componentID)."_Name_From_Index"}
	    );
    }
    # Generate data content.
    @dataContent =
	(
	 {
	     intrinsic  => "class",
	     type       => "nodeComponent",
	     attributes => [ "intent(in   )" ],
	     variables  => [ "self" ]
	 },
	 {
	     intrinsic  => "integer",
	     attributes => [ "intent(inout)" ],
	     variables  => [ "count" ]
	 },
	 {
	     intrinsic  => "type",
	     type       => "varying_string",	
	     attributes => [ "intent(inout)" ],
	     variables  => [ "name" ]
	 }
	);
    # Generate the function.
    $functionCode  = "  subroutine Node_Component_Name_From_Index(self,count,name)\n";
    $functionCode .= "    !% Return the name of the property of given index.\n";
    $functionCode .= "    use ISO_Varying_String\n";
    $functionCode .= "    implicit none\n";
    $functionCode .= &Format_Variable_Defintions(\@dataContent)."\n";
    $functionCode .= "    return\n";
    $functionCode .= "  end subroutine Node_Component_Name_From_Index\n\n";
    # Insert into the function list.
    push(
	@{$buildData->{'code'}->{'functions'}},
	$functionCode
	);
    # Insert a type-binding for this function into the implementation type.
    push(
	@{$buildData->{'types'}->{'nodeComponent'}->{'boundFunctions'}},
	{type => "procedure", name => "nameFromIndex", function => "Node_Component_Name_From_Index"}
	);
}

sub Generate_Implementation_Serialization_Functions {
    # Generate serialization/deserialization functions for each component implementation.
    my $buildData = shift;
    # Initialize function code.
    my $functionCode;
    # Initialize data content.
    my @dataContent;
    # Iterate over component implementations.
    foreach my $componentID ( @{$buildData->{'componentIdList'}} ) {
	# Get the component.
	my $component = $buildData->{'components'}->{$componentID};
	# Generate data content.
	@dataContent =
	    (
	     {
		 intrinsic  => "class",
		 type       => "nodeComponent".ucfirst($componentID),
		 attributes => [ "intent(in   )" ],
		 variables  => [ "self" ]
	     }
	    );
	# Generate a count function.
  	$functionCode  = "  integer function Node_Component_".ucfirst($componentID)."_Count(self)\n";
	$functionCode .= "    !% Return a count of the serialization of a ".$component->{'name'}." implementation of the ".$component->{'class'}." component.\n";
	$functionCode .= "    implicit none\n";
	$functionCode .= &Format_Variable_Defintions(\@dataContent)."\n";
	$functionCode .= "    Node_Component_".ucfirst($componentID)."_Count=0\n";
	# Initialize a count of scalar properties.
	my $scalarPropertyCount = 0;
	# Iterate over properties.
	foreach my $methodName ( keys(%{$component->{'methods'}->{'method'}}) ) {
	    my $method = $component->{'methods'}->{'method'}->{$methodName};
   	    # Check if this method has any linked data in this component.
	    if ( exists($method->{'linkedData'}) ) {
		# For each linked datum count if necessary.
		my $linkedDataName = $method->{'linkedData'};
		my $linkedData     = $component->{'content'}->{'data'}->{$linkedDataName};
		if ( $linkedData->{'isEvolvable'} eq "true" ) {
		    if ( $linkedData->{'rank'} == 0 ) {
			switch ( $linkedData->{'type'} ) {
			    case ( "real" ) {
				++$scalarPropertyCount;
			    }
			    else {
				$functionCode .= "    Node_Component_".ucfirst($componentID)."_Count=Node_Component_".ucfirst($componentID)."_Count+self%".padLinkedData($linkedDataName,[0,0])."%value%serializeCount()\n";
			    }
			}
		    } else {
			$functionCode .= "    if (allocated(self%".padLinkedData($linkedDataName,[0,0])."%value)) Node_Component_".ucfirst($componentID)."_Count=Node_Component_".ucfirst($componentID)."_Count+size(self%".padLinkedData($linkedDataName,[0,0])."%value)\n";
		    }
		}
	    }
	}
	# Insert the final count of scalar properties.
	$functionCode .= "    Node_Component_".ucfirst($componentID)."_Count=Node_Component_".ucfirst($componentID)."_Count+".$scalarPropertyCount."\n"
	    if ($scalarPropertyCount > 0);
	$functionCode .= "    return\n";
	$functionCode .= "  end function Node_Component_".ucfirst($componentID)."_Count\n\n";
	# Insert into the function list.
	push(
	    @{$buildData->{'code'}->{'functions'}},
	    $functionCode
	    );
	# Insert a type-binding for this function into the treeNode type.
	push(
	    @{$buildData->{'types'}->{'nodeComponent'.ucfirst($componentID)}->{'boundFunctions'}},
	    {type => "procedure", name => "serializeCount", function => "Node_Component_".ucfirst($componentID)."_Count"}
	);
	# Iterate over content types.
	foreach my $content ( "value", "scale", "rate" ) {
	    # Specify data content for serialization functions.
	    @dataContent =
	    	(
	    	 {
	    	     intrinsic  => "class",
	    	     type       => "nodeComponent".ucfirst($componentID),
	    	     attributes => [ "intent(in   )" ],
	    	     variables  => [ "self" ]
	    	 },
	    	 {
	    	     intrinsic  => "double precision",
	    	     attributes => [ "intent(  out)", "dimension(:)" ],
	    	     variables  => [ "array" ]
	    	 }
	    	);
	    # Generate serialization function.
	    $functionCode  = "  subroutine Node_Component_".ucfirst($componentID)."_Serialize_".ucfirst($content)."s(self,array)\n";
	    $functionCode .= "    !% Serialize ".$content."s of a ".$component->{'name'}." implementation of the ".$component->{'class'}." component.\n";
	    $functionCode .= "    implicit none\n";
	    my $serializationCode;
	    my $needCount = 0;
	    foreach my $methodName ( keys(%{$component->{'methods'}->{'method'}}) ) {
		my $method = $component->{'methods'}->{'method'}->{$methodName};
	    	# Check if this method has any linked data in this component.
	    	if ( exists($method->{'linkedData'}) ) {
	    	    my $linkedDataName = $method->{'linkedData'};
	    	    my $linkedData     = $component->{'content'}->{'data'}->{$linkedDataName};
	    	    if ( $linkedData->{'isEvolvable'} eq "true" ) {
	    		if ( $linkedData->{'rank'} == 0 ) {
	    		    switch ( $linkedData->{'type'} ) {
	    			case ( "real" ) {
				    $serializationCode .= "    write (0,*) 'DEBUG -> Node_Component_".ucfirst($componentID)."_Serialize_".ucfirst($content)."s -> ".$linkedDataName."',offset,size(array)\n"
					if ( $debugging == 1 );
	    			    $serializationCode .= "    array(offset)=self%".padLinkedData($linkedDataName,[0,0])."%".$content."\n";
				    $serializationCode .= "    if (array(offset) <= 0.0d0) write (0,*) 'DEBUG -> Node_Component_".ucfirst($componentID)."_Serialize_".ucfirst($content)."s: non-positive scale found for ".$linkedDataName."'\n"
					if ( $content eq "scale" && $debugging == 1 );
	    			    $serializationCode .= "    offset=offset+1\n";
	    			}
	    			else {
	    			    $serializationCode .= "    count=self%".padLinkedData($linkedDataName,[0,0])."%".pad($content,5)."%serializeCount(                            )\n";
				    $serializationCode .= "    write (0,*) 'DEBUG -> Node_Component_".ucfirst($componentID)."_Serialize_".ucfirst($content)."s -> ".$linkedDataName."',offset,count,size(array)\n"
					if ( $debugging == 1 );
	    			    $serializationCode .= "    if (count > 0) call  self%".padLinkedData($linkedDataName,[0,0])."%".pad($content,5)."%serialize     (array(offset:offset+count-1))\n";
				    $serializationCode .= "    if (count > 0 .and. any(array(offset:offset+count-1) <= 0.0d0)) write (0,*) 'DEBUG -> Node_Component_".ucfirst($componentID)."_Serialize_".ucfirst($content)."s: non-positive scale found for ".$linkedDataName."'\n"
					if ( $content eq "scale" && $debugging == 1 );
	    			    $serializationCode .= "    offset=offset+count\n";
				    $needCount = 1;
	    			}
	    		    }
	    		} else {
	    		    $serializationCode .= "    if (allocated(self%".padLinkedData($linkedDataName,[0,0])."%".pad($content,5).")) then\n";
	    		    $serializationCode .= "       count=size(self%".padLinkedData($linkedDataName,[0,0])."%".pad($content,5).")\n";
			    $serializationCode .= "    write (0,*) 'DEBUG -> Node_Component_".ucfirst($componentID)."_Serialize_".ucfirst($content)."s -> ".$linkedDataName."',offset,count,size(array)\n"
				if ( $debugging == 1 );
	    		    $serializationCode .= "       array(offset:offset+count-1)=reshape(self%".padLinkedData($linkedDataName,[0,0])."%".$content.",[count])\n";
			    $serializationCode .= "       if (any(array(offset:offset+count-1) <= 0.0d0)) write (0,*) 'DEBUG -> Node_Component_".ucfirst($componentID)."_Serialize_".ucfirst($content)."s: non-positive scale found for ".$linkedDataName."'\n"
					if ( $content eq "scale" && $debugging == 1 );
	    		    $serializationCode .= "       offset=offset+count\n";
	    		    $serializationCode .= "    end if\n";
			    $needCount = 1;
	    		}
	    	    }
	    	}
	    }
	    if ( defined($serializationCode) ) {
		my @variables = ( "offset" );
		push(@variables,"count") 
		    if ( $needCount == 1 );
		push(
		    @dataContent,
		    {
			intrinsic  => "integer",
			variables  => \@variables
		    }    
		    );
		$serializationCode = "    offset=1\n".$serializationCode;
	    }
	    $functionCode .= &Format_Variable_Defintions(\@dataContent)."\n";
	    $functionCode .= $serializationCode
		if ( defined($serializationCode) );
	    $functionCode .= "    return\n";
	    $functionCode .= "  end subroutine Node_Component_".ucfirst($componentID)."_Serialize_".ucfirst($content)."s\n\n";
	    # Insert into the function list.
	    push(
	    	@{$buildData->{'code'}->{'functions'}},
	    	$functionCode
	    	);
	    # Insert a type-binding for this function into the implementation type.
	    push(
		@{$buildData->{'types'}->{'nodeComponent'.ucfirst($componentID)}->{'boundFunctions'}},
		{type => "procedure", name => "serialize".ucfirst($content)."s", function => "Node_Component_".ucfirst($componentID)."_Serialize_"  .ucfirst($content)."s"},
		);
	    # Specify data content for deserialization functions.
	    @dataContent =
	    	(
	    	 {
	    	     intrinsic  => "class",
	    	     type       => "nodeComponent".ucfirst($componentID),
	    	     attributes => [ "intent(inout)" ],
	    	     variables  => [ "self" ]
	    	 },
	    	 {
	    	     intrinsic  => "double precision",
	    	     attributes => [ "intent(in   )", "dimension(:)" ],
	    	     variables  => [ "array" ]
	    	 }
	    	);
	    # Generate deserialization function.
	    $functionCode  = "  subroutine Node_Component_".ucfirst($componentID)."_Deserialize_".ucfirst($content)."s(self,array)\n";
	    $functionCode .= "    !% Serialize ".$content."s of a ".$component->{'name'}." implementation of the ".$component->{'class'}." component.\n";
	    $functionCode .= "    implicit none\n";
	    my $deserializationCode;
	    $needCount = 0;
	    foreach my $methodName ( keys(%{$component->{'methods'}->{'method'}}) ) {
	    	my $method = $component->{'methods'}->{'method'}->{$methodName};
	    	# Check if this method has any linked data in this component.
	    	if ( exists($method->{'linkedData'}) ) {
	    	    # For each linked datum count if necessary.
	    	    my $linkedDataName = $method->{'linkedData'};
	    	    my $linkedData     = $component->{'content'}->{'data'}->{$linkedDataName};
	    	    if ( $linkedData->{'isEvolvable'} eq "true" ) {
	    		if ( $linkedData->{'rank'} == 0 ) {
	    		    switch ( $linkedData->{'type'} ) {
	    			case ( "real" ) {
	    			    $deserializationCode .= "    self%".padLinkedData($linkedDataName,[0,0])."%".pad($content,5)."=array(offset)\n";
	    			    $deserializationCode .= "    offset=offset+1\n";
	    			}
	    			else {
	    			    $deserializationCode .= "    count=self%".padLinkedData($linkedDataName,[0,0])."%".pad($content,5)."%serializeCount(                            )\n";
	    			    $deserializationCode .= "    call  self%".padLinkedData($linkedDataName,[0,0])."%".pad($content,5)."%deserialize   (array(offset:offset+count-1))\n";
	    			    $deserializationCode .= "    offset=offset+count\n";
				    $needCount = 1;
	    			}
	    		    }
	    		} else {
	    		    $deserializationCode .= "    if (allocated(self%".padLinkedData($linkedDataName,[0,0])."%".pad($content,5).")) then\n";
	    		    $deserializationCode .= "       count=size(self%".padLinkedData($linkedDataName,[0,0])."%".pad($content,5).")\n";
	    		    $deserializationCode .= "       self%".padLinkedData($linkedDataName,[0,0])."%".pad($content,5)."=reshape(array(offset:offset+count-1),shape(self%".padLinkedData($linkedDataName,[0,0])."%".pad($content,5)."))\n";
	    		    $deserializationCode .= "       offset=offset+count\n";
	    		    $deserializationCode .= "    end if\n";
			    $needCount = 1;
	    		}
	    	    }
	    	}
	    }
	    if ( defined($deserializationCode) ) {
		my @variables = ( "offset" );
		push(@variables,"count") 
		    if ( $needCount == 1 );
		push(
		    @dataContent,
		    {
			intrinsic  => "integer",
			variables  => \@variables
		    }    
		    );
		$deserializationCode = "    offset=1\n".$deserializationCode;
	    }
	    $functionCode .= &Format_Variable_Defintions(\@dataContent)."\n";
	    $functionCode .= $deserializationCode
		if ( defined($deserializationCode) );
	    $functionCode .= "    return\n";
	    $functionCode .= "  end subroutine Node_Component_".ucfirst($componentID)."_Deserialize_".ucfirst($content)."s\n\n";
	    # Insert into the function list.
	    push(
	    	@{$buildData->{'code'}->{'functions'}},
	    	$functionCode
	    	);
	    # Insert a type-binding for this function into the implementation type.
	    push(
		@{$buildData->{'types'}->{'nodeComponent'.ucfirst($componentID)}->{'boundFunctions'}},
		{type => "procedure", name => "deserialize".ucfirst($content)."s", function => "Node_Component_".ucfirst($componentID)."_Deserialize_"  .ucfirst($content)."s"},
		);
	}
    }
}

sub Generate_Component_Count_Functions {
    # Generate component count functions.
    my $buildData = shift;
    # Initialize function code.
    my $functionCode;
    # Initialize data content.
    my @dataContent;
    # Create methods to get components.
    foreach my $componentClassName ( @{$buildData->{'componentClassList'}} ) {
	# Specify data content for component get function.
	@dataContent =
	    (
	     {
		 intrinsic  => "class",
		 type       => "treeNode",
		 attributes => [ "intent(inout)", "target" ],
		 variables  => [ "self" ]
	     }
	    );
	# Generate code for the count function.
    	$functionCode  = "  integer function ".$componentClassName."CountLinked(self)\n";
	$functionCode .= "    !% Returns the number of {\\tt ".$componentClassName."} components in {\\tt self}.\n";
    	$functionCode .= "    implicit none\n";
	$functionCode .= &Format_Variable_Defintions(\@dataContent)."\n";
    	$functionCode .= "    select type (self)\n";
    	$functionCode .= "    class is (treeNode)\n";
	$functionCode .= "     if (allocated(self%component".ucfirst($componentClassName).")) then\n";
	$functionCode .= "       select type (component => self%component".ucfirst($componentClassName)."(1))\n";
	$functionCode .= "       type is (nodeComponent".ucfirst($componentClassName).")\n";
	$functionCode .= "         ".$componentClassName."CountLinked=0\n";
	$functionCode .= "       class default\n";
	$functionCode .= "         ".$componentClassName."CountLinked=size(self%component".ucfirst($componentClassName).")\n";
	$functionCode .= "       end select\n";
	$functionCode .= "     else\n";
	$functionCode .= "        ".$componentClassName."CountLinked=0\n";
	$functionCode .= "     end if\n";
    	$functionCode .= "    end select\n";
    	$functionCode .= "    return\n";
    	$functionCode .= "  end function ".$componentClassName."CountLinked\n";
	# Insert into the function list.
	push(
	    @{$buildData->{'code'}->{'functions'}},
	    $functionCode
	    );
	# Bind this function to the treeNode type.
	push(
	    @{$buildData->{'types'}->{'treeNode'}->{'boundFunctions'}},
	    {type => "procedure", name => $componentClassName."Count", function => $componentClassName."CountLinked"}
	    );
    }
}

sub Generate_Component_Get_Functions {
    # Generate component get methods.
    my $buildData = shift;
    # Initialize function code.
    my $functionCode;
    # Initialize data content.
    my @dataContent;
    # Create methods to get components.
    foreach my $componentClassName ( @{$buildData->{'componentClassList'}} ) {
	# Specify data content for component get function.
	@dataContent =
	    (
	     {
		 intrinsic  => "class",
		 type       => "nodeComponent".ucfirst($componentClassName),
		 attributes => [ "pointer" ],
		 variables  => [ $componentClassName."Get" ]
	     },
	     {
		 intrinsic  => "class",
		 type       => "treeNode",
		 attributes => [ "intent(inout)", "target" ],
		 variables  => [ "self" ]
	     },
	     {
		 intrinsic  => "integer",
		 attributes => [ "intent(in   )", "optional" ],
		 variables  => [ "instance" ]
	     },
	     {
		 intrinsic  => "logical",
		 attributes => [ "intent(in   )", "optional" ],
		 variables  => [ "autoCreate" ]
	     },
	     {
		 intrinsic  => "integer",
		 variables  => [ "instanceActual" ]
	     },
	     {
		 intrinsic  => "logical",
		 variables  => [ "autoCreateActual" ]
	     }
	    );
	# Generate code for the get function.
    	$functionCode  = "  recursive function ".$componentClassName."Get(self,instance,autoCreate)\n";
	$functionCode .= "    !% Returns the {\\tt ".$componentClassName."} component of {\\tt self}.\n";
	$functionCode .= "    use Galacticus_Error\n";   
 	$functionCode .= "    implicit none\n";
	$functionCode .= &Format_Variable_Defintions(\@dataContent)."\n";
    	$functionCode .= "    select type (self)\n";
    	$functionCode .= "    class is (treeNode)\n";
    	$functionCode .= "       instanceActual=1\n";
    	$functionCode .= "       if (present(instance)) instanceActual=instance\n";
    	$functionCode .= "       autoCreateActual=.false.\n";
    	$functionCode .= "       if (present(autoCreate)) autoCreateActual=autoCreate\n";
	$functionCode .= "       if (allocated(self%component".ucfirst($componentClassName).")) then\n";
	# If we are allowed to autocreate the component and it still has generic type then deallocate it to
	# force it to be created later.
	$functionCode .= "         if (autoCreateActual.and.same_type_as(self%component".ucfirst($componentClassName)."(1),".ucfirst($componentClassName)."Class)) deallocate(self%component".ucfirst($componentClassName).")\n";
	$functionCode .= "       end if\n";
    	$functionCode .= "       if (.not.allocated(self%component".ucfirst($componentClassName).")) then\n";
    	$functionCode .= "         if (autoCreateActual) then\n";
	$functionCode .= "            call self%".lc($componentClassName)."Create"."()\n";
	$functionCode .= "         else\n";
	$functionCode .= "            call Galacticus_Error_Report('".$componentClassName."Get','component is not allocated')\n";
	$functionCode .= "         end if\n";
    	$functionCode .= "       end if\n";
    	$functionCode .= "       ".$componentClassName."Get => self%component".ucfirst($componentClassName)."(instanceActual)\n";
    	$functionCode .= "    end select\n";
    	$functionCode .= "    return\n";
    	$functionCode .= "  end function ".$componentClassName."Get\n";
	# Insert into the function list.
	push(
	    @{$buildData->{'code'}->{'functions'}},
	    $functionCode
	    );
	# Bind this function to the treeNode type.
	push(
	    @{$buildData->{'types'}->{'treeNode'}->{'boundFunctions'}},
	    {type => "procedure", name => $componentClassName, function => $componentClassName."Get"}
	    );
	# Specify data content for create-by-interrupt function.
	@dataContent =
	    (
	     {
		 intrinsic  => "type",
		 type       => "treeNode",
		 attributes => [ "pointer", "intent(inout)" ],
		 variables  => [ "self" ]
	     },
	     {
		 intrinsic  => "class",
		 type       => "nodeComponent".ucfirst($componentClassName),
		 attributes => [ "pointer" ],
		 variables  => [ $componentClassName ]
	     }
	    );
	# Generate function to create component via an interrupt.
   	$functionCode  = "  subroutine ".$componentClassName."CreateByInterrupt(self)\n";
	$functionCode .= "    !% Create the {\\tt ".$componentClassName."} component of {\\tt self} via an interrupt.\n";
    	$functionCode .= "    implicit none\n";
	$functionCode .= &Format_Variable_Defintions(\@dataContent)."\n";
	$functionCode .= "    ".$componentClassName." => self%".$componentClassName."(autoCreate=.true.)\n";
	# Loop over instances of this class, and call custom create routines if necessary.
	my $foundCreateFunctions = 0;
    	foreach my $componentName ( @{$buildData->{'componentClasses'}->{$componentClassName}->{'members'}} ) {
	    my $componentID = ucfirst($componentClassName).ucfirst($componentName);
	    my $component = $buildData->{'components'}->{$componentID};
	    if ( exists($component->{'createFunction'}) ) {
		if ( $foundCreateFunctions == 0 ) {
		    $functionCode .= "    select type (".$componentClassName.")\n";
		    $foundCreateFunctions = 1;
		}
		$functionCode .= "    type is (nodeComponent".padFullyQualified(ucfirst($componentID),[0,0]).")\n";
		my $createFunction = $component->{'createFunction'};
		$createFunction = $component->{'createFunction'}->{'content'}
		    if ( exists($component->{'createFunction'}->{'content'}) );
		$createFunction = $componentID."CreateFunction"
		    if (
			exists($component->{'createFunction'}->{'isDeferred'}          ) &&
		      	       $component->{'createFunction'}->{'isDeferred'} eq "true"
		    );
		$functionCode .= "       call ".$createFunction."(".$componentClassName.")\n";
	    }
	}
	$functionCode .= "    end select\n"
	    unless ( $foundCreateFunctions == 0 );
    	$functionCode .= "    return\n";
    	$functionCode .= "  end subroutine ".$componentClassName."CreateByInterrupt\n";
	# Insert into the function list.
	push(
	    @{$buildData->{'code'}->{'functions'}},
	    $functionCode
	    );
	# If any create function is deferred, create a function to set it at runt time.
    	foreach my $componentName ( @{$buildData->{'componentClasses'}->{$componentClassName}->{'members'}} ) {
	    my $componentID = ucfirst($componentClassName).ucfirst($componentName);
	    my $component = $buildData->{'components'}->{$componentID};
	    if (
		exists($component->{'createFunction'}                           ) && 
		exists($component->{'createFunction'}->{'isDeferred'}           ) &&
		       $component->{'createFunction'}->{'isDeferred'} eq "true"
		) {
		$functionCode  = "   subroutine ".$componentID."CreateFunctionSet(createFunction)\n";
		$functionCode .= "     !% Set the create function for the {\\tt ".$componentID."} component.\n";
		$functionCode .= "     implicit none\n";
		$functionCode .= "     external createFunction\n";
		my $createFunction = $componentID."CreateFunction"; 
		$functionCode .= "     ".$createFunction." => createFunction\n";
		$functionCode .= "     return\n";
		$functionCode .= "   end subroutine ".$componentID."CreateFunctionSet\n";
		# Insert into the function list.
		push(
		    @{$buildData->{'code'}->{'functions'}},
		    $functionCode
		    );
	    }
	}
    }
}

sub Generate_Component_Destruction_Functions {
    # Generate component destruction functions.
    my $buildData = shift;
    # Iterate over component classes.
    foreach my $componentClassName ( @{$buildData->{'componentClassList'}} ) {
	my @dataContent =
	    (
	     {
		 intrinsic  => "class",
		 type       => "treeNode",
		 attributes => [ "intent(inout)" ],
		 variables  => [ "self" ]
	     },
	     {
		 intrinsic  => "integer",
		 variables  => [ "i" ]
	     }
	    );
    	my $functionCode = "  subroutine ".$componentClassName."DestroyLinked(self)\n";
	$functionCode   .= "    !% Destroy the {\\tt ".$componentClassName."} component of {\\tt self}.\n";
    	$functionCode   .= "    implicit none\n";
	$functionCode   .= &Format_Variable_Defintions(\@dataContent)."\n";
    	$functionCode   .= "    if (allocated(self%component".ucfirst($componentClassName).")) then\n";
	$functionCode   .= "      do i=1,size(self%component".ucfirst($componentClassName).")\n";
	$functionCode   .= "        call        self%component".ucfirst($componentClassName)."(i)%destroy()\n";
	$functionCode   .= "      end do\n";
	$functionCode   .= "      deallocate (self%component".ucfirst($componentClassName).")\n";
	$functionCode   .= "    end if\n";
    	$functionCode   .= "    return\n";
    	$functionCode   .= "  end subroutine ".$componentClassName."DestroyLinked\n\n";
	# Insert into the function list.
	push(
	    @{$buildData->{'code'}->{'functions'}},
	    $functionCode
	    );
	# Bind this function to the treeNode type.
	push(
	    @{$buildData->{'types'}->{'treeNode'}->{'boundFunctions'}},
	    {type => "procedure", name => $componentClassName."Destroy" , function => $componentClassName."DestroyLinked" }
	    );
    }
}

sub Generate_Component_Creation_Functions {
    # Generate component creation functions.
    my $buildData = shift;
    # Iterate over component classes.
    foreach my $componentClassName ( @{$buildData->{'componentClassList'}} ) {
	# Specify data content.
	my @dataContent =
	    (
	     {
		 intrinsic  => "class",
		 type       => "treeNode",
		 attributes => [ "target", "intent(inout)" ],
		 variables  => [ "self" ]
	     },
	     {
		 intrinsic  => "class",
		 type       => "nodeComponent".ucfirst($componentClassName),
		 attributes => [ "intent(in   )", "optional" ],
		 variables  => [ "template" ]
	     },
	     {
		 intrinsic  => "type",
		 type       => "varying_string",
		 variables  => [ "message" ]
	     },
	     {
		 intrinsic  => "integer",
		 variables  => [ "i" ]
	     }
	    );
	# Generate value initialization code.
	my $initializeCode;
	my %requiredComponents;
    	foreach my $componentName ( @{$buildData->{'componentClasses'}->{$componentClassName}->{'members'}} ) {
	    my $componentID = ucfirst($componentClassName).ucfirst($componentName);
	    my $component = $buildData->{'components'}->{$componentID};
	    $initializeCode .= "    do i=1,size(self%component".ucfirst($componentClassName).")\n";
	    $initializeCode .= "      select type (component => self%component".ucfirst($componentClassName)."(i))\n";
	    $initializeCode .= "    class is (nodeComponent".$componentID.")\n";
 	    foreach my $methodName ( keys(%{$component->{'methods'}->{'method'}}) ) {
		my $method = $component->{'methods'}->{'method'}->{$methodName};
		if ( exists($method->{'linkedData'}) ) {
		    my $linkedDataName = $method->{'linkedData'};
		    my $linkedData     = $component->{'content'}->{'data'}->{$linkedDataName};
		    # Set to a class default value if available.
		    if ( exists($method->{'classDefault'}) ) {
			my $default = $method->{'classDefault'}->{'code'};
			while ( $default =~ m/self([a-zA-Z]+)Component\s*%/ ) {
			    $requiredComponents{$1} = 1;
			    $default =~ s/self([a-zA-Z]+)Component\s*%//;
			}
			$default = $method->{'classDefault'}->{'code'};
			$default =~ s/self\s*%/component%/g;
			if ( exists($method->{'classDefault'}->{'count'}) ) {
			    my @gsr = ( "value" );
			    push
				(
				 @gsr,
				 "rate",
				 "scale"
				)
				if ( $method->{'attributes'}->{'isEvolvable'} eq "true" );
			    foreach ( @gsr ) {
				$initializeCode .= "           call Alloc_Array(component%".padLinkedData($linkedDataName,[0,0])."%".$_.",[".$method->{'classDefault'}->{'count'}."])\n";
			    }
			}
			$initializeCode .= "            component%".padLinkedData($linkedDataName,[0,0])."%value=".$default."\n";
		    } else {
			# Set to null.
			switch ( $linkedData->{'type'} ) {
			    case ( "real"    ) {
				if ( $linkedData->{'rank'} == 0 ) {
				    $initializeCode .= "            component%".padLinkedData($linkedDataName,[0,0])."%value=0.0d0\n";
				} else {
				    $initializeCode .= "            call Alloc_Array(component%".padLinkedData($linkedDataName,[0,0])."%value,[".join(",","0" x $linkedData->{'rank'})."])\n";
				}
			    }
			    case ( "integer" ) {
				$initializeCode .= "            component%".padLinkedData($linkedDataName,[0,0])."%value=0\n";
			    }
			    case ( "logical" ) {
				$initializeCode .= "            component%".padLinkedData($linkedDataName,[0,0])."%value=.false.\n";
			    }
			    else {
				$initializeCode .= "       call component%".padLinkedData($linkedDataName,[0,0])."%value%reset()\n";			    
			    }
			}
		    }
		}
	    }
	    $initializeCode .= "    end select\n";
	    $initializeCode .= "    end do\n";
	}
	# Add pointers for each required component.
	push(
	    @dataContent,
	    {
		intrinsic  => "class",
		type       => "nodeComponent".ucfirst($_),
		attributes => [ "pointer" ],
		variables  => [ "self".ucfirst($_)."Component" ]
	    }
	    )
	    foreach ( keys(%requiredComponents) );
	# Generate function code.
	my $functionCode;
    	$functionCode .= "  subroutine ".$componentClassName."CreateLinked(self,template)\n";
	$functionCode .= "    !% Create the {\\tt ".$componentClassName."} component of {\\tt self}.\n";
	$functionCode .= "    use ISO_Varying_String\n";
	$functionCode .= "    use Galacticus_Display\n";
	$functionCode .= "    use String_Handling\n";
	# Insert any required modules.
	my %requiredModules;
    	foreach my $componentName ( @{$buildData->{'componentClasses'}->{$componentClassName}->{'members'}} ) {
	    my $componentID = ucfirst($componentClassName).ucfirst($componentName);
	    my $component = $buildData->{'components'}->{$componentID};
 	    foreach my $methodName ( keys(%{$component->{'methods'}->{'method'}}) ) {
		my $method = $component->{'methods'}->{'method'}->{$methodName};
		if ( exists($method->{'classDefault'}) && exists($method->{'classDefault'}->{'modules'}) ) {
		    foreach ( @{$method->{'classDefault'}->{'modules'}} ) {
			$requiredModules{$_} = 1;
		    }
		}
	    }
	}
	foreach ( keys(%requiredModules) ) {
	    $functionCode .= "    use ".$_."\n";
	}
    	$functionCode .= "    implicit none\n";
 	$functionCode .= &Format_Variable_Defintions(\@dataContent)."\n";
	$functionCode .= "    if (Galacticus_Verbosity_Level() >= verbosityInfo) then\n";
	$functionCode .= "      message='Creating ".$componentClassName." in node '\n";
	$functionCode .= "      message=message//self%index()\n";
	$functionCode .= "      call Galacticus_Display_Message(message,verbosityInfo)\n";
	$functionCode .= "    end if\n";
    	$functionCode .= "    if (present(template)) then\n";
	$functionCode .= "       allocate(self%component".ucfirst($componentClassName)."(1),source=template)\n";
    	$functionCode .= "    else\n"; 
	$functionCode .= "       allocate(self%component".ucfirst($componentClassName)."(1),source="."default".ucfirst($componentClassName)."Component)\n";
   	$functionCode .= "    end if\n";
     	$functionCode .= "    select type (self)\n";
	$functionCode .= "    type is (treeNode)\n";
	$functionCode .= "      do i=1,size(self%component".ucfirst($componentClassName).")\n";
	$functionCode .= "        self%component".ucfirst($componentClassName)."(i)%hostNode => self\n";
	$functionCode .= "      end do\n";
    	$functionCode .= "    end select\n";
	foreach my $requiredComponent ( keys(%requiredComponents) ) {
	    $functionCode .= "     self".$requiredComponent."Component => self%".lc($requiredComponent)."()\n";
	}
	$functionCode .= $initializeCode;
    	$functionCode .= "    return\n";
    	$functionCode .= "  end subroutine ".$componentClassName."CreateLinked\n\n";
	# Insert into the function list.
	push(
	    @{$buildData->{'code'}->{'functions'}},
	    $functionCode
	    );
	# Bind this function to the treeNode type.
	push(
	    @{$buildData->{'types'}->{'treeNode'}->{'boundFunctions'}},
	    {type => "procedure", name => $componentClassName."Create" , function => $componentClassName."CreateLinked" }
	    );
    }
}

sub Generate_Node_Copy_Function {
    # Generate function to copy one node to another.
    my $buildData = shift;
    # Specify variables.
    my @dataContent =
	(
	 {
	     intrinsic  => "class",
	     type       => "treeNode",
	     attributes => [ "intent(in   )" ],
	     variables  => [ "self" ]
	 },
	 {
	     intrinsic  => "class",
	     type       => "treeNode",
	     attributes => [ "intent(inout)" ],
	     variables  => [ "targetNode" ]
	 },
	 {
	     intrinsic  => "logical",
	     attributes => [ "intent(in   )", "optional" ],
	     variables  => [ "skipFormationNode" ]
	 },
	 {
	     intrinsic  => "logical",
	     variables  => [ "skipFormationNodeActual" ]
	 }
	);
    push(
	@dataContent,
	{
	    intrinsic  => "integer",
	    variables  => [ "i" ]
	}
	)
	if ( $workaround == 1 );
    # Generate the code.
    my $functionCode;
    $functionCode .= "  subroutine Tree_Node_Copy_Node_To(self,targetNode,skipFormationNode)\n";
    $functionCode .= "    !% Make a copy of {\\tt self} in {\\tt targetNode}.\n";
    $functionCode .= "    implicit none\n";
    $functionCode .= &Format_Variable_Defintions(\@dataContent)."\n";
    $functionCode .= "    skipFormationNodeActual=.false.\n";
    $functionCode .= "    if (present(skipFormationNode)) skipFormationNodeActual=skipFormationNode\n";
    $functionCode .= "    targetNode%".padComponentClass($_,[8,14])." =  self%".$_."\n"
	foreach ( "uniqueIdValue", "indexValue" );
    $functionCode .= "    targetNode%".padComponentClass($_,[8,14])." => self%".$_."\n"
	foreach ( "parent", "firstChild", "sibling", "firstSatellite", "mergeTarget", "firstMergee", "siblingMergee");
    $functionCode .= "    if (.not.skipFormationNodeActual) targetNode%formationNode => self%formationNode\n";
    # Loop over all component classes
    if ( $workaround == 1 ) {
	foreach my $componentClassName ( @{$buildData->{'componentClassList'}} ) {
	    $functionCode .= "    if (allocated(targetNode%component".padComponentClass(ucfirst($componentClassName),[0,0]).")) deallocate(targetNode%component".padComponentClass(ucfirst($componentClassName),[0,0]).")\n";
	    $functionCode .= "    allocate(targetNode%component".padComponentClass(ucfirst($componentClassName),[0,0])."(size(self%component".padComponentClass(ucfirst($componentClassName),[0,0]).")),source=self%component".padComponentClass(ucfirst($componentClassName),[0,0]).")\n";
	    $functionCode .= "    do i=1,size(self%component".padComponentClass(ucfirst($componentClassName),[0,0]).")\n";
	    foreach my $implementationName ( @{$buildData->{'componentClasses'}->{$componentClassName}->{'members'}} ) {
		$functionCode .= "      select type (from => self%component".padComponentClass(ucfirst($componentClassName),[0,0]).")\n";
		$functionCode .= "      type is (nodeComponent".padFullyQualified(ucfirst($componentClassName).ucfirst($implementationName),[0,0]).")\n";
		$functionCode .= "        select type (to => targetNode%component".padComponentClass(ucfirst($componentClassName),[0,0]).")\n";
		$functionCode .= "        type is (nodeComponent".padFullyQualified(ucfirst($componentClassName).ucfirst($implementationName),[0,0]).")\n";
		$functionCode .= "          to=from\n";
		$functionCode .= "        end select\n";
		$functionCode .= "      end select\n";
	    }
	    $functionCode .= "    end do\n";
	}
    } else {
	foreach ( @{$buildData->{'componentClassList'}} ) {
	    $functionCode .= "    targetNode%component".padComponentClass(ucfirst($_),[0,14])."=  self%component".ucfirst($_)."\n";
	}
    }
    # Update target node pointers.
    $functionCode .= "    select type (targetNode)\n";
    $functionCode .= "    type is (treeNode)\n";
    foreach ( @{$buildData->{'componentClassList'}} ) {
	$functionCode .= "      do i=1,size(self%component".padComponentClass(ucfirst($_),[0,0]).")\n";
	
	$functionCode .= "        targetNode%component".padComponentClass(ucfirst($_),[0,14])."(i)%hostNode =>  targetNode\n";
	$functionCode .= "      end do\n";
    }
    $functionCode .= "    end select\n";
    $functionCode .= "    return\n";
    $functionCode .= "  end subroutine Tree_Node_Copy_Node_To\n\n";
    # Insert into the function list.
    push(
	@{$buildData->{'code'}->{'functions'}},
	$functionCode
	);
    # Insert a type-binding for this function into the treeNode type.
    push(
	@{$buildData->{'types'}->{'treeNode'}->{'boundFunctions'}},
	{type => "procedure", name => "copyNodeTo", function => "Tree_Node_Copy_Node_To"}
	);
}

sub Generate_Node_Move_Function {
    # Generate function to move one node to another.
    my $buildData = shift;
    # Specify variables.
    my @dataContent =
	(
	 {
	     intrinsic  => "class",
	     type       => "treeNode",
	     attributes => [ "intent(inout)" ],
	     variables  => [ "self" ]
	 },
	 {
	     intrinsic  => "type",
	     type       => "treeNode",
	     attributes => [ "intent(inout)", "target" ],
	     variables  => [ "targetNode" ]
	 },
	 {
	     intrinsic  => "integer",
	     variables  => [ "i" ]
	 }
	);
    # Generate the code.
    my $functionCode;
    # Create functions for moving node components.
    $functionCode .= "  subroutine Tree_Node_Move_Components(self,targetNode)\n";
    $functionCode .= "    !% Move components from {\\tt self} to {\\tt targetNode}.\n";
    $functionCode .= "    implicit none\n";
    $functionCode .= &Format_Variable_Defintions(\@dataContent)."\n";
     # Loop over all component classes
    foreach ( @{$buildData->{'componentClassList'}} ) {	    
	$functionCode .= "    if (allocated(targetNode%component".padComponentClass(ucfirst($_),[0,0]).")) then\n";
	$functionCode .= "      do i=1,size(targetNode%component".padComponentClass(ucfirst($_),[0,0]).")\n";
	$functionCode .= "        call targetNode%component".padComponentClass(ucfirst($_),[0,0])."(i)%destroy()\n";
	$functionCode .= "      end do\n";
	$functionCode .= "    end if\n";
	$functionCode .= "    if (allocated(self      %component".padComponentClass(ucfirst($_),[0,0]).")) then\n";
	$functionCode .= "       call Move_Alloc(self%component".padComponentClass(ucfirst($_),[0,0]).",targetNode%component".padComponentClass(ucfirst($_),[0,0]).")\n";
	$functionCode .= "       do i=1,size(targetNode%component".padComponentClass(ucfirst($_),[0,0]).")\n";
	$functionCode .= "         targetNode%component".padComponentClass(ucfirst($_),[0,0])."(i)%hostNode => targetNode\n";
	$functionCode .= "       end do\n";
	$functionCode .= "    end if\n";
    }
    $functionCode .= "    return\n";
    $functionCode .= "  end subroutine Tree_Node_Move_Components\n\n";
    # Insert into the function list.
    push(
	@{$buildData->{'code'}->{'functions'}},
	$functionCode
	);
    # Insert a type-binding for this function into the treeNode type.
    push(
	@{$buildData->{'types'}->{'treeNode'}->{'boundFunctions'}},
	{type => "procedure", name => "moveComponentsTo", function => "Tree_Node_Move_Components"}
	);
}

sub Generate_Deferred_Function_Attacher {
    # Generate functions to attach a function to a deferred method and to query the attachment state.
    my $component = shift;
    my $method    = shift;
    my $buildData = shift;
    my $gsr       = shift;
    my $gsrSuffix = "";
    $gsrSuffix = ucfirst($gsr)
	unless ( $gsr eq "get" );
    # Get the component fully-qualified, class and metho names.
    my $componentClassName = $component->{'class'             };
    my $componentName      = $component->{'fullyQualifiedName'};
    my $methodName         = $method->{'name'};
    # Skip if this function was already created.
    my $recordName = $componentClassName.$methodName.$gsr;
    unless ( exists($buildData->{'deferredFunctionComponentClassMethodsMade'}->{$recordName}) ) {
	# Define the function name.
	my $functionName = $componentClassName.ucfirst($methodName).ucfirst($gsr)."Function";
	# Define the data content.
	my $type = "";
	$type = $componentName.ucfirst($methodName).ucfirst($gsr)
	    if ( $gsr eq "get" );
	my @dataContent =
	    (
	     {
		 intrinsic  => "procedure",
		 type       => $type,
		 variables  => [ "deferredFunction" ]
	     },
	    );
	# Construct the function code.
	my $functionCode;
	$functionCode  = "  subroutine ".$functionName."(deferredFunction)\n";
	$functionCode .= "    !% Set the function to be used for ".$gsr." of the {\\tt ".$methodName."} property of the {\\tt ".$componentClassName."} component class.\n";
	$functionCode .= "    implicit none\n";
	$functionCode .= &Format_Variable_Defintions(\@dataContent)."\n";
	$functionCode .= "    ".$componentClassName.ucfirst($methodName).ucfirst($gsr)."Deferred       => deferredFunction\n";
	$functionCode .= "    ".$componentClassName.ucfirst($methodName).ucfirst($gsr)."IsAttachedValue=  .true.\n";
	$functionCode .= "    return\n";
	$functionCode .= "  end subroutine ".$functionName."\n";
	# Insert into the function list.
	push(
	    @{$buildData->{'code'}->{'functions'}},
	    $functionCode
	    );
	# Bind this function to the relevant type.
	if ( 
	    ( $method->{'attributes'}->{'bindsTo'} ne "top" && ( $gsr eq "get" || $gsr eq "set" ) ) ||
	    (                                                    $gsr eq "rate"                   )
	    ) {
	    push(
		@{$buildData->{'types'}->{"nodeComponent".ucfirst($componentClassName)}->{'boundFunctions'}},
		{type => "procedure", pass => "nopass", name => $methodName.$gsrSuffix."Function", function => $functionName}
		);
	}
	# Also create a function to return whether or not the deferred function has been attached.
	$functionCode  = "  logical function ".$componentClassName.ucfirst($methodName).ucfirst($gsr)."IsAttached()\n";
	$functionCode .= "    !% Return true if the deferred function used to ".$gsr." the {\\tt ".$methodName."} property of the {\\tt ".$componentClassName."} component class has been attached.\n";
	$functionCode .= "    implicit none\n";
	$functionCode .= "    ".$componentClassName.ucfirst($methodName).ucfirst($gsr)."IsAttached=".$componentClassName.ucfirst($methodName).ucfirst($gsr)."IsAttachedValue\n";
	$functionCode .= "    return\n";
	$functionCode .= "  end function ".$componentClassName.ucfirst($methodName).ucfirst($gsr)."IsAttached\n";
	# Insert into the function list.
	push(
	    @{$buildData->{'code'}->{'functions'}},
	    $functionCode
	    );
	# Bind this function to the relevant type.
	if ( 
	    ( $method->{'attributes'}->{'bindsTo'} ne "top" && ( $gsr eq "get" || $gsr eq "set" ) ) ||
	    (                                                    $gsr eq "rate"                   )
	    ) {
	    push(
		@{$buildData->{'types'}->{"nodeComponent".ucfirst($componentClassName)}->{'boundFunctions'}},
		{type => "procedure", pass => "nopass", name => $methodName.$gsrSuffix."IsAttached", function => $componentClassName.ucfirst($methodName).ucfirst($gsr)."IsAttached"}
		);
	}
	# Record that these functions have now been created.
	$buildData->{'deferredFunctionComponentClassMethodsMade'}->{$recordName} = 1;
    }
}

sub Generate_Deferred_GSR_Function {
    # Generate function to get/set/rate the value of a method via a deferred function.
    my $buildData = shift;
    # Record bindings already made.
    my %bindings;
    # Iterate over component implementations
    foreach my $componentID ( @{$buildData->{'componentIdList'}} ) {
	# Get the component.
	my $component = $buildData->{'components'}->{$componentID};
	# Iterate over methods.
	foreach my $methodName ( keys(%{$component->{'methods'}->{'method'}}) ) {
	    # Get the method.
	    my $method = $component->{'methods'}->{'method'}->{$methodName};
	    # Get the component fully-qualified and class names.
	    my $componentClassName = $component->{'class'             };
	    my $componentName      = $component->{'fullyQualifiedName'};
	    # Ignore non-deferred functions.
	    unless ( $method->{'attributes'}->{'isDeferred'} eq "" ) {
		# Function code data.
		my $functionCode;
		# Get the name of the method.
		my $methodName = $method->{'name'};
		# Get properties of the data type needed.
		(my $dataDefinition, my $label) = &Data_Object_Definition($method,matchOnly => 1);
		# Identify methods with a deferred get function to be built.
		if (
		    $method->{'attributes' }->{'isDeferred'} =~ m/get/ &&
		    $method->{'attributes' }->{'isGettable'} eq "true" &&
		    $method->{'getFunction'}->{'build'     } eq "true"
		    )
		{
		    # Define data content of this function.
		    @{$dataDefinition->{'variables'}} = ( $componentName.ucfirst($methodName)."Get" );
		    my @dataContent =
			(
			 $dataDefinition,
			 {
			     intrinsic  => "class",
			     type       => "nodeComponent".ucfirst($componentName),
			     attributes => [ "intent(inout)" ],
			     variables  => [ "self" ]
			 },
			);
		    # Construct the function code.
		    $functionCode  = "  function ".$componentName.ucfirst($methodName)."Get(self)\n";
		    $functionCode .= "    !% Get the value of the {\\tt ".$methodName."} property of the {\\tt ".$componentName."} component using a deferred function.\n";
		    $functionCode .= "    implicit none\n";
		    $functionCode .= &Format_Variable_Defintions(\@dataContent)."\n";
		    $functionCode .= "    ".$componentName.ucfirst($methodName)."Get=".$componentClassName.ucfirst($methodName)."GetDeferred(self)\n";
		    $functionCode .= "    return\n";
		    $functionCode .= "  end function ".$componentName.ucfirst($methodName)."Get\n";
		    # Insert into the function list.
		    push(
			@{$buildData->{'code'}->{'functions'}},
			$functionCode
			);
		    # Bind this function to the relevant type.
		    push(
			@{$buildData->{'types'}->{"nodeComponent".ucfirst($componentName)}->{'boundFunctions'}},
			{type => "procedure", name => $methodName, function => $componentName.ucfirst($methodName)."Get"}
			);
		    # Generate an attacher function.
		    &Generate_Deferred_Function_Attacher($component,$method,$buildData,"get");
		}
		# Add an "intent(in)" attribute to the data definition for set and rate functions.
		push(@{$dataDefinition->{'attributes'}},"intent(in   )");
		# Identify methods with a deferred set function to be built.
		if (
		    $method->{'attributes' }->{'isDeferred'} =~ m/set/
		    && $method->{'attributes' }->{'isSettable'} eq "true" 
		    && $method->{'setFunction'}->{'build'     } eq "true"
		    )
		{
		    @{$dataDefinition->{'variables'}} = ( "setValue" );
		    my @dataContent =
			(
			 $dataDefinition,
			 {
			     intrinsic  => "class",
			     type       => "nodeComponent".ucfirst($componentName),
			     attributes => [ "intent(inout)" ],
			     variables  => [ "self" ]
			 },
			);
		    $functionCode  = "  subroutine ".$componentName.ucfirst($methodName)."Set(self,setValue)\n";
		    $functionCode .= "    !% Set the value of the {\\tt ".$methodName."} property of the {\\tt ".$componentName."} component using a deferred function.\n";
		    $functionCode .= "    implicit none\n";
		    $functionCode .= &Format_Variable_Defintions(\@dataContent)."\n";
		    $functionCode .= "    call ".$componentClassName.ucfirst($methodName)."SetDeferred(self,setValue)\n";
		    $functionCode .= "    return\n";
		    $functionCode .= "  end subroutine ".$componentName.ucfirst($methodName)."Set\n\n";
		    # Insert into the function list.
		    push(
			@{$buildData->{'code'}->{'functions'}},
			$functionCode
			);
		    # Bind this function to the relevant type.
		    push(
			@{$buildData->{'types'}->{"nodeComponent".ucfirst($componentName)}->{'boundFunctions'}},
			{type => "procedure", name => $methodName."Set", function => $componentName.ucfirst($methodName)."Set"}
			);
		    # Generate an attacher function.
		    &Generate_Deferred_Function_Attacher($component,$method,$buildData,"set");		  
		}
		# Identify methods with a deferred rate function to be built.
		if (
		    $method->{'attributes' }->{'isDeferred' } =~ m/rate/
		    && $method->{'attributes' }->{'isEvolvable'} eq "true" 
		    )
		{
		    # Define data content of this function.
		    @{$dataDefinition->{'variables'}} = ( "setValue" );
		    my $type = "nodeComponent";
		    $type .= ucfirst($componentName)
			unless ( $method->{'attributes' }->{'bindsTo'} eq "top" );
		    my @dataContent =
			(
			 $dataDefinition,
			 {
			     intrinsic  => "class",
			     type       => $type,
			     attributes => [ "intent(inout)" ],
			     variables  => [ "self" ]
			 },
			 {
			     intrinsic  => "logical",
			     attributes => [ "intent(inout)", "optional" ],
			     variables  => [ "interrupt" ]
			 },
			 {
			     intrinsic  => "procedure",
			     type       => "",
			     attributes => [ "pointer", "optional", "intent(inout)" ],
			     variables  => [ "interruptProcedure" ]
			 }
			);
		    $functionCode  = "  subroutine ".$componentName.ucfirst($methodName)."Rate(self,setValue,interrupt,interruptProcedure)\n";
		    $functionCode .= "    !% Set the rate of the {\\tt ".$methodName."} property of the {\\tt ".$componentName."} component using a deferred function.\n";
		    $functionCode .= "    implicit none\n";
		    $functionCode .= &Format_Variable_Defintions(\@dataContent)."\n";
		    $functionCode .= "    call ".$componentClassName.ucfirst($methodName)."RateDeferred(self,setValue,interrupt,interruptProcedure)\n";
		    $functionCode .= "    return\n";
		    $functionCode .= "  end subroutine ".$componentName.ucfirst($methodName)."Rate\n\n";
		    # Insert into the function list.
		    push(
			@{$buildData->{'code'}->{'functions'}},
			$functionCode
			);
		    # Bind this function to the relevant type.
		    my $bindingName = $type.$methodName."Rate";
		    unless ( exists($bindings{$bindingName}) && $method->{'attributes' }->{'bindsTo'} eq "top" ) {
			push(
			    @{$buildData->{'types'}->{$type}->{'boundFunctions'}},
			    {type => "procedure", name => $methodName."Rate", function => $componentName.ucfirst($methodName)."Rate"}
			    );
			$bindings{$bindingName} = 1;
		    }
		    # Generate an attacher function.
		    &Generate_Deferred_Function_Attacher($component,$method,$buildData,"rate");
		}
	    }
	}
    }
}

sub Generate_GSR_Functions {
    # Generate functions to get/set/rate the value of a method directly.
    my $buildData = shift;
    # Initialize function code.
    my $functionCode;
    # Initialize data content.
    my @dataContent;
    # Initialize records of functions created.
    my %classRatesCreated;
    my %deferredFunctionComponentClassMethodsMade;
    # Iterate over component implementations.
    foreach my $componentID ( @{$buildData->{'componentIdList'}} ) {
	my $component = $buildData->{'components'}->{$componentID};
	# Get the parent class.
	my $componentClassName = $component->{'class'};
	# Iterate over methods.
	foreach my $methodName ( keys(%{$component->{'methods'}->{'method'}}) ) {
	    my $method = $component->{'methods'}->{'method'}->{$methodName};
	    # Handle cases where a get function is explicitly specified for a non-deferred virtual method.
	    if (
		$method->{'attributes' }->{'isGettable'} eq "true"      &&
		$method->{'getFunction'}->{'build'     } eq "false"     &&
		$method->{'getFunction'}->{'bindsTo'   } eq "component" &&
		$method->{'attributes' }->{'isDeferred'} !~ m/get/
		) {
		# No need to build the function - just insert a type-binding into the implementation type.
		push(
		    @{$buildData->{'types'}->{'nodeComponent'.ucfirst($componentID)}->{'boundFunctions'}},
		    {type => "procedure", name => $methodName, function => $method->{'getFunction'}->{'content'}}
		    );
	    }
	    # Handle cases where a set function is explicitly specified for a non-deferred virtual method.
	    if (
		$method->{'attributes' }->{'isSettable'} eq "true"      &&
		$method->{'setFunction'}->{'build'     } eq "false"     &&
		$method->{'setFunction'}->{'bindsTo'   } eq "component" &&
		$method->{'attributes' }->{'isDeferred'} !~ m/set/
		) {
		# No need to build the function - just insert a type-binding into the implementation type.
		push(
		    @{$buildData->{'types'}->{'nodeComponent'.ucfirst($componentID)}->{'boundFunctions'}},
		    {type => "procedure", name => $methodName."Set", function => $method->{'setFunction'}->{'content'}}
		    );	
	    }
   	    # Check if this method has any linked data in this component.
	    if ( exists($method->{'linkedData'}) ) {
		# Get the linked data.
		my $linkedDataName = $method->{'linkedData'};
		my $linkedData     = $component->{'content'}->{'data'}->{$linkedDataName};
		# Create a "get" function if the property is gettable.
		if ( $method->{'attributes'}->{'isGettable'} eq "true" ) {
		    # Skip get function creation if a custom function which binds at the component level has been specified.
		    unless (
			$method->{'getFunction'}->{'build'  } eq "false"     &&
			$method->{'getFunction'}->{'bindsTo'} eq "component"
			)
		    {
			# Determine the suffix for this function.
			my $suffix = "";
			$suffix = "Value"
			    if ( $method->{'attributes' }->{'isDeferred'} =~ m/get/ );
			# Specify the data content.
			(my $dataDefinition,my $label) = &Data_Object_Definition($linkedData);
			push(@{$dataDefinition->{'variables'}},$componentID.ucfirst($methodName)."Get".$suffix);
			@dataContent = (
			    $dataDefinition,
			    {
				intrinsic  => "class",
				type       => "nodeComponent".ucfirst($componentID),
				attributes => [ "intent(inout)" ],
				variables  => [ "self" ]
			    }
			    );
			# Generate the code.
			$functionCode  = "  function ".$componentID.ucfirst($methodName)."Get".$suffix."(self)\n";
			$functionCode .= "    !% Return the {\\tt ".$methodName."} property of the {\\tt ".$componentID."} component implementation.\n";
			$functionCode .= "    implicit none\n";
			$functionCode .= &Format_Variable_Defintions(\@dataContent)."\n";
			$functionCode .= "    ".$componentID.$methodName."Get".$suffix."=self%".$linkedDataName."%value\n";
			$functionCode .= "    return\n";
			$functionCode .= "  end function ".$componentID.ucfirst($methodName)."Get".$suffix."\n\n";
			# Insert into the function list.
			push(
			    @{$buildData->{'code'}->{'functions'}},
			    $functionCode
			    );
			# Insert a type-binding for this function into the implementation type.
			push(
			    @{$buildData->{'types'}->{'nodeComponent'.ucfirst($componentID)}->{'boundFunctions'}},
			    {type => "procedure", name => $methodName.$suffix, function => $componentID.ucfirst($methodName)."Get".$suffix}
			    );
		    }
		}
		# Create a "set" method unless the property is not settable or a custom set function has been specified.
		if ( $method->{'attributes' }->{'isSettable'} eq "true" ) {
		    if ( $method->{'setFunction'}->{'build'} eq "true" ) {
			# Determine the suffix for this function.
			my $suffix = "";
			$suffix = "Value"
			    if ( $method->{'attributes' }->{'isDeferred'} =~ m/set/ );
			# Specify the data content.
			(my $dataDefinition,my $label) = &Data_Object_Definition($linkedData,matchOnly => 1);
			push(@{$dataDefinition->{'attributes'}},"intent(in   )");
			push(@{$dataDefinition->{'variables' }},"setValue"     );
			@dataContent = (
			    $dataDefinition,
			    {
				intrinsic  => "class",
				type       => "nodeComponent".ucfirst($componentID),
				attributes => [ "intent(inout)" ],
				variables  => [ "self" ]
			    }
			    );
			# Generate the function code.
			$functionCode  = "  subroutine ".$componentID.ucfirst($methodName)."Set".$suffix."(self,setValue)\n";
			$functionCode .= "    !% Set the {\\tt ".$methodName."} property of the {\\tt ".$componentID."} component implementation.\n";
			$functionCode .= "    use Memory_Management\n";
			$functionCode .= "    implicit none\n";
			$functionCode .= &Format_Variable_Defintions(\@dataContent)."\n";
			# For non-real properties we also set the rate and scale content. This ensures that they get reallocated to
			# the correct size.
			my @setContent = ( "value" );
			push(@setContent,"rate","scale")
			    if ( $method->{'attributes' }->{'isEvolvable'} eq "true" );
			switch ( $linkedData->{'rank'} ) {
			    case ( 0 ) {
				$functionCode .= "    self%".$linkedDataName."%".pad($_,5)."=setValue\n"
				    foreach ( @setContent );
			    }
			    case ( 1 ) {
				foreach ( @setContent ) {
				    $functionCode .= "    if (.not.allocated(self%".$linkedDataName."%".pad($_,5).")) then\n";
				    $functionCode .= "       call    Alloc_Array  (self%".$linkedDataName."%".pad($_,5).",shape(setValue))\n";
				    $functionCode .= "    else\n";
				    $functionCode .= "       if (size(self%".$linkedDataName."%".pad($_,5).") /= size(setValue)) then\n";
				    $functionCode .= "          call Dealloc_Array(self%".$linkedDataName."%".pad($_,5)."                )\n";
				    $functionCode .= "          call Alloc_Array  (self%".$linkedDataName."%".pad($_,5).",shape(setValue))\n";
				    $functionCode .= "       end if\n";
				    $functionCode .= "    end if\n";
				}
				$functionCode .= "    self%".$linkedDataName."%".pad($_,5)."=setValue\n"
				    foreach ( @setContent );
			    }
			}
			$functionCode .= "    return\n";
			$functionCode .= "  end subroutine ".$componentID.ucfirst($methodName)."Set".$suffix."\n\n";
			# Insert into the function list.
			push(
			    @{$buildData->{'code'}->{'functions'}},
			    $functionCode
			    );
			# Insert a type-binding for this function into the implementation type.
			push(
			    @{$buildData->{'types'}->{'nodeComponent'.ucfirst($componentID)}->{'boundFunctions'}},
			    {type => "procedure", name => $methodName."Set".$suffix, function => $componentID.ucfirst($methodName)."Set".$suffix}
			    );
		    }
		}
		# Create "count", "rate" and "scale" functions if the property is evolvable.
		if ( $method->{'attributes'}->{'isEvolvable'} eq "true" ) {
		    # Specify the "count" function data content.
		    @dataContent = (
			{
			    intrinsic  => "class",
			    type       => "nodeComponent".ucfirst($componentID),
			    attributes => [ "intent(in   )" ],
			    variables  => [ "self" ]
			}
			);
		    # Generate the "count" function code.
		    $functionCode  = "  integer function ".$componentID.ucfirst($methodName)."Count(self)\n";
		    $functionCode .= "    !% Return a count of the number of scalar properties in the {\\tt ".$methodName."} property of the {\\tt ".lcfirst($componentID)."} component implementation.\n";
		    $functionCode .= "    implicit none\n";
		    $functionCode .= &Format_Variable_Defintions(\@dataContent)."\n";
		    switch ( $linkedData->{'rank'} ) {
			case ( 0 ) {
			    $functionCode .= "    ".$componentID.$methodName."Count=1\n";
			}
			case ( 1 ) {
			    $functionCode .= "    if (allocated(self%".$linkedDataName."%value)) then\n";
			    $functionCode .= "    ".$componentID.$methodName."Count=size(self%".$linkedDataName."%value)\n";
			    $functionCode .= "    else\n";
			    $functionCode .= "    ".$componentID.$methodName."Count=0\n";
			    $functionCode .= "    end if\n";
			}
		    }
		    $functionCode .= "    return\n";
		    $functionCode .= "  end function ".$componentID.ucfirst($methodName)."Count\n\n";
		    # Insert into the function list.
		    push(
			@{$buildData->{'code'}->{'functions'}},
			$functionCode
			);
		    # Insert a type-binding for this function into the implementation type.
		    push(
			@{$buildData->{'types'}->{'nodeComponent'.ucfirst($componentID)}->{'boundFunctions'}},
			{type => "procedure", name => $methodName."Count", function => $componentID.ucfirst($methodName)."Count"}
			);
		    # Get the data content for remaining functions.
		    (my $dataDefinition,my $label) = &Data_Object_Definition($linkedData,matchOnly => 1);
		    push(@{$dataDefinition->{'variables' }},"setValue"     );
		    push(@{$dataDefinition->{'attributes'}},"intent(in   )");
		    # Skip rate function creation if the rate function is deferred.
		    unless ( $method->{'attributes'}->{'isDeferred'} =~ m/rate/ ) {
			# Specify the "rate" function data content.
			@dataContent = (
			    $dataDefinition,
			    {
				intrinsic  => "class",
				type       => "nodeComponent".ucfirst($componentID),
				attributes => [ "intent(inout)" ],
				variables  => [ "self" ]
			    },
			    {
				intrinsic  => "logical",
				attributes => [ "optional", "intent(inout)" ],
				variables  => [ "interrupt" ]
			    },
			    {
				intrinsic  => "procedure",
				type       => "",
				attributes => [ "optional", "intent(inout)", "pointer" ],
				variables  => [ "interruptProcedure" ]
			    }
			    );
			# Generate the rate function code.
			$functionCode  = "  subroutine ".$componentID.ucfirst($methodName)."Rate(self,setValue,interrupt,interruptProcedure)\n";
			$functionCode .= "    !% Accumulate to the {\\tt ".$methodName."} property rate of change of the {\\tt ".$componentID."} component implementation.\n";
			$functionCode .= "    implicit none\n";
			$functionCode .= &Format_Variable_Defintions(\@dataContent)."\n";
			if ( $linkedData->{'type'} eq "real" ) {
			    $functionCode .= "    self%".$linkedDataName."%rate=self%".$linkedDataName."%rate+setValue\n";
			} else {
			    $functionCode .= "    call self%".$linkedDataName."%rate%increment(setValue)\n";
			}
			$functionCode .= "    return\n";
			$functionCode .= "  end subroutine ".$componentID.ucfirst($methodName)."Rate\n\n";
			# Insert into the function list.
			push(
			    @{$buildData->{'code'}->{'functions'}},
			    $functionCode
			    );
			# Insert a type-binding for this function into the implementation type.
			push(
			    @{$buildData->{'types'}->{'nodeComponent'.ucfirst($componentID)}->{'boundFunctions'}},
			    {type => "procedure", name => $methodName."Rate", function => $componentID.ucfirst($methodName)."Rate"}
			    );
		    }
		    if ( $method->{'attributes' }->{'makeGeneric'} eq "true" ) {
			# Create a version of this rate function which binds to the top-level class, and so is suitable for
			# attaching to inter-component pipes.
			# Specify the data content.		
			@dataContent = (
			    $dataDefinition,
			    {
				intrinsic  => "class",
				type       => "nodeComponent",
				attributes => [ "intent(inout)" ],
				variables  => [ "self" ]
			    },
			    {
				intrinsic  => "logical",
				attributes => [ "optional", "intent(inout)" ],
				variables  => [ "interrupt" ]
			    },
			    {
				intrinsic  => "procedure",
				type       => "",
				attributes => [ "optional", "intent(inout)", "pointer" ],
				variables  => [ "interruptProcedure" ]
			    },
			    {
				intrinsic  => "class",
				type       => "nodeComponent".ucfirst($componentClassName),
				attributes => [ "pointer" ],
				variables  => [ "this".ucfirst($componentClassName) ]
			    },
			    {
				intrinsic  => "type",
				type       => "treeNode",
				attributes => [ "pointer" ],
				variables  => [ "thisNode" ]
			    }
			    );
			# Generate the function code.
			$functionCode  = "  subroutine ".$componentID.ucfirst($methodName)."RateGeneric(self,setValue,interrupt,interruptProcedure)\n";
			$functionCode .= "    !% Set the rate of the {\\tt ".$methodName."} property of the {\\tt ".$componentID."} component via a generic {\\tt nodeComponent}.\n";
			$functionCode .= "    use Galacticus_Error\n"
			    if ( $method->{'attributes'}->{'createIfNeeded'} eq "true" );
			$functionCode .= "    implicit none\n";
			$functionCode .= &Format_Variable_Defintions(\@dataContent)."\n";
			$functionCode .= "    thisNode => self%host()\n";
			$functionCode .= "    this".ucfirst($componentClassName)." => thisNode%".$componentClassName."()\n";
			if ( $method->{'attributes'}->{'createIfNeeded'} eq "true" ) {
			    $functionCode .= "    select type (this".ucfirst($componentClassName).")\n";
			    $functionCode .= "    type is (nodeComponent".ucfirst($componentClassName).")\n";
			    $functionCode .= "      ! No specific component exists, we must interrupt and create one.\n";
			    if ( $linkedData->{'rank'} == 0 ) {
				switch ( $linkedData->{'type'} ) {
				    case ( "real" ) {
					$functionCode .= "   if (setValue == 0.0d0) return\n";
				    }
				    case ( "integer" ) {
					die('integer should not be evolvable')
					}
				    case ( "logical" ) {
					die('logical should not be evolvable')
					}
				    else {
					$functionCode .= "   if (setValue%isZero()) return\n";
				    }
				}
			    } else {
				switch ( $linkedData->{'type'} ) {
				    case ( "real" ) {
					$functionCode .= "   if (all(setValue == 0.0d0)) return\n";
				    }
				    case ( "integer" ) {
					die('integer should not be evolvable')
					}
				    case ( "logical" ) {
					die('logical should not be evolvable')
					}
				    else {
					die('auto-create of rank>0 objects not supported');
				    }
				}
			    }
			    $functionCode .= "      if (.not.(present(interrupt).and.present(interruptProcedure))) call Galacticus_Error_Report('".$componentID.ucfirst($methodName)."RateGeneric','interrupt required, but optional arguments missing')\n";
			    $functionCode .= "      interrupt=.true.\n";
			    $functionCode .= "      interruptProcedure => ".$componentClassName."CreateByInterrupt\n";
			    $functionCode .= "      return\n";
			    $functionCode .= "    end select\n";
			}
			$functionCode .= "    call this".ucfirst($componentClassName)."%".$methodName."Rate(setValue,interrupt,interruptProcedure)\n";
			$functionCode .= "    return\n";
			$functionCode .= "  end subroutine ".$componentID.ucfirst($methodName)."RateGeneric\n\n";
			# Insert into the function list.
			push(
			    @{$buildData->{'code'}->{'functions'}},
			    $functionCode
			    );
		    }
		    if ( $method->{'attributes' }->{'createIfNeeded'} eq "true" ) {
			# Create a version of this rate function which binds to the component class, and so can auto-create the component as needed.
			my $label = $componentClassName.ucfirst($methodName);
			unless ( exists($classRatesCreated{$label}) ) {
			    $classRatesCreated{$label} = 1;
			    # Specify the data content.		
			    @dataContent = (
				$dataDefinition,
				{
				    intrinsic  => "class",
				    type       => "nodeComponent".ucfirst($componentClassName),
				    attributes => [ "intent(inout)" ],
				    variables  => [ "self" ]
				},
				{
				    intrinsic  => "logical",
				    attributes => [ "optional", "intent(inout)" ],
				    variables  => [ "interrupt" ]
				},
				{
				    intrinsic  => "procedure",
				    type       => "",
				    attributes => [ "optional", "intent(inout)", "pointer" ],
				    variables  => [ "interruptProcedure" ]
				}
				);
			    # Generate the function code.
			    $functionCode  = "  subroutine ".$componentClassName.ucfirst($methodName)."Rate(self,setValue,interrupt,interruptProcedure)\n";
			    $functionCode .= "    !% Accept a rate set for the {\\tt ".$methodName."} property of the {\\tt ".$componentClassName."} component class. Trigger an interrupt to create the component.\n";
			    $functionCode .= "    use Galacticus_Error\n";
			    $functionCode .= "    implicit none\n";
			    $functionCode .= &Format_Variable_Defintions(\@dataContent)."\n";
			    $functionCode .= "    ! No specific component exists, so we must interrupt and create one unless the rate is zero.\n";
			    if ( $linkedData->{'rank'} == 0 ) {
				switch ( $linkedData->{'type'} ) {
				    case ( "real" ) {
					$functionCode .= "   if (setValue == 0.0d0) return\n";
				    }
				    case ( "integer" ) {
					die('integer should not be evolvable')
				    }
				    case ( "logical" ) {
					die('logical should not be evolvable')
				    }
				    else {
					$functionCode .= "   if (setValue%isZero()) return\n";
				    }
				}
			    } else {
				switch ( $linkedData->{'type'} ) {
				    case ( "real" ) {
					$functionCode .= "   if (all(setValue == 0.0d0)) return\n";
				    }
				    case ( "integer" ) {
					die('integer should not be evolvable')
				    }
				    case ( "logical" ) {
					die('logical should not be evolvable')
				    }
				    else {
					die('auto-create of rank>0 objects not supported');
				    }
				}
			    }
			    $functionCode .= "    if (.not.(present(interrupt).and.present(interruptProcedure))) call Galacticus_Error_Report('".$componentClassName.ucfirst($methodName)."Rate','interrupt required, but optional arguments missing')\n";
			    $functionCode .= "    interrupt=.true.\n";
			    $functionCode .= "    interruptProcedure => ".$componentClassName."CreateByInterrupt\n";
			    $functionCode .= "    return\n";
			    $functionCode .= "  end subroutine ".$componentClassName.ucfirst($methodName)."Rate\n\n";
			    # Insert into the function list.
			    push(
				@{$buildData->{'code'}->{'functions'}},
				$functionCode
				);
			}
		    }
		    # Specify the data content for the "scale" function.		
		    @dataContent = (
			$dataDefinition,
			{
			    intrinsic  => "class",
			    type       => "nodeComponent".ucfirst($componentID),
			    attributes => [ "intent(inout)" ],
			    variables  => [ "self" ]
			}
			);
		    # Generate a function to set the "scale".
		    $functionCode  = "  subroutine ".$componentID.ucfirst($methodName)."Scale(self,setValue)\n";
		    $functionCode .= "    !% Set the {\\tt ".$methodName."} property scale of the {\\tt ".$componentID."} component implementation.\n";
		    $functionCode .= "    implicit none\n";
		    $functionCode .= &Format_Variable_Defintions(\@dataContent)."\n";
 		    $functionCode .= "    self%".$linkedDataName."%scale=setValue\n";
		    $functionCode .= "    return\n";
		    $functionCode .= "  end subroutine ".$componentID.ucfirst($methodName)."Scale\n\n";
		    # Insert into the function list.
		    push(
			@{$buildData->{'code'}->{'functions'}},
			$functionCode
			);
		    # Insert a type-binding for this function into the implementation type.
		    push(
			@{$buildData->{'types'}->{'nodeComponent'.ucfirst($componentID)}->{'boundFunctions'}},
			{type => "procedure", name => $methodName."Scale", function => $componentID.ucfirst($methodName)."Scale"}
			);
		}
	    }
	}
    }
}

sub Generate_Tree_Node_Creation_Function {
    # Generate a tree node creation function.
    my $buildData = shift;
    # Specify data content.
    my @dataContent =
	(
	 {
	     intrinsic  => "class",
	     type       => "treeNode",
	     attributes => [ "intent(inout)", "target" ],
	     variables  => [ "self" ]
	 },
	 {
	     intrinsic  => "integer",
	     type       => "kind=kind_int8",
	     attributes => [ "intent(in   )", "optional" ],
	     variables  => [ "index" ]
	 }
	);
    # Create the function code.
    my $functionCode;
    $functionCode .= "  subroutine treeNodeInitialize(self,index)\n";
    $functionCode .= "    !% Initialize a {\\tt treeNode} object.\n";
    $functionCode .= "    use Galacticus_Error\n";
    $functionCode .= "    implicit none\n";
    $functionCode .= &Format_Variable_Defintions(\@dataContent)."\n";
    $functionCode .= "    ! Ensure pointers are nullified.\n";
    $functionCode .= "    nullify (self%".padComponentClass($_,[9,14]).")\n"
	foreach ( "parent", "firstChild", "sibling", "firstSatellite", "mergeTarget", "firstMergee", "siblingMergee", "formationNode" );
    foreach ( @{$buildData->{'componentClassList'}} ) {
    	$functionCode .= "    allocate(self%".padComponentClass("component".ucfirst($_),[9,14])."(1))\n";
    }
    $functionCode .= "    select type (self)\n";
    $functionCode .= "    type is (treeNode)\n";
    foreach ( @{$buildData->{'componentClassList'}} ) {
	$functionCode .= "       self%component".padComponentClass(ucfirst($_),[0,0])."(1)%hostNode => self\n";
    }
    $functionCode .= "    end select\n";
    $functionCode .= "    ! Assign index if supplied.\n";
    $functionCode .= "    if (present(index)) call self%indexSet(index)\n";
    $functionCode .= "    ! Assign a unique ID.\n";
    $functionCode .= "    !\$omp atomic\n";
    $functionCode .= "    uniqueIDCount=uniqueIDCount+1\n";
    $functionCode .= "    if (uniqueIDCount <= 0) call Galacticus_Error_Report('treeNodeInitialize','ran out of unique ID numbers')\n";
    $functionCode .= "    self%uniqueIdValue=uniqueIDCount\n";
    $functionCode .= "    return\n";
    $functionCode .= "  end subroutine treeNodeInitialize\n";	
    # Insert into the function list.
    push(
	@{$buildData->{'code'}->{'functions'}},
	$functionCode
	);
}

sub Generate_Tree_Node_Destruction_Function {
    # Generate a tree node destruction function.
    my $buildData = shift;
    # Specify data content.
    my @dataContent =
	(
	 {
	     intrinsic  => "class",
	     type       => "treeNode",
	     attributes => [ "intent(inout)" ],
	     variables  => [ "self" ]
	 }
	);
    # Create the function code.
    my $functionCode;
    $functionCode .= "  subroutine treeNodeDestroy(self)\n";
    $functionCode .= "    !% Destroy a {\\tt treeNode} object.\n";
    $functionCode .= "    implicit none\n";
    $functionCode .= &Format_Variable_Defintions(\@dataContent)."\n";
    foreach my $componentClass ( @{$buildData->{'componentClassList'}} ) {
     	$functionCode .= "    call self%".padComponentClass(lc($componentClass)."Destroy",[7,0])."()\n";
    }
    $functionCode .= "    return\n";
    $functionCode .= "  end subroutine treeNodeDestroy\n";	
    # Insert into the function list.
    push(
	@{$buildData->{'code'}->{'functions'}},
	$functionCode
	);
}

sub Generate_Type_Name_Functions {
    # Generate a type name functions.
    my $buildData = shift;
    # Initialize data content.
    my @dataContent;
    # Initialize the function code.
    my $functionCode;
    # Iterate over component classes.
    foreach ( @{$buildData->{'componentClassList'}} ) {
	# Specify data content.
	@dataContent =
	    (
	     {
		 intrinsic  => "class",
		 type       => "nodeComponent".ucfirst($_),
		 attributes => [ "intent(in   )" ],
		 variables  => [ "self" ]
	     },
	     {
		 intrinsic  => "type",
		 type       => "varying_string",
		 variables  => [ "Node_Component_".ucfirst($_)."_Type" ]
	     }
	    );
	# Create the function code.
	$functionCode  = "  function Node_Component_".ucfirst($_)."_Type(self)\n";
	$functionCode .= "     !% Returns the type for the ".$_." component.\n";
	$functionCode .= "     implicit none\n";
	$functionCode .= &Format_Variable_Defintions(\@dataContent)."\n";
	$functionCode .= "     ".padComponentClass("Node_Component_".ucfirst($_)."_Type",[20,0])."='nodeComponent:".$_."'\n";
	$functionCode .= "     return\n";
	$functionCode .= "  end function Node_Component_".ucfirst($_)."_Type\n\n";
	# Insert into the function list.
	push(
	    @{$buildData->{'code'}->{'functions'}},
	    $functionCode
	    );

	# Bind this function to the relevant type.
	push(
	    @{$buildData->{'types'}->{'nodeComponent'.ucfirst($_)}->{'boundFunctions'}},
	    {type => "procedure", name => "type", function => "Node_Component_".ucfirst($_)."_Type"}
	    );
    }
    # Iterate over implementations.
    foreach my $componentName ( @{$buildData->{'componentIdList'}} ) {
	my $component = $buildData->{'components'}->{$componentName};
	# Specify data content.
	@dataContent =
	    (
	     {
		 intrinsic  => "class",
		 type       => "nodeComponent".ucfirst($componentName),
		 attributes => [ "intent(in   )" ],
		 variables  => [ "self" ]
	     },
	     {
		 intrinsic  => "type",
		 type       => "varying_string",
		 variables  => [ "Node_Component_".ucfirst($componentName)."_Type" ]
	     }
	    );
	# Create the function code.
  	$functionCode  = "  function Node_Component_".ucfirst($componentName)."_Type(self)\n";
	$functionCode .= "    !% Returns the type for the ".$component->{'name'}." implementation of the ".$component->{'class'}." component.\n";
	$functionCode .= "    implicit none\n";
	$functionCode .= &Format_Variable_Defintions(\@dataContent)."\n";
	$functionCode .= "    ".padImplementationProperty("Node_Component_".ucfirst($componentName)."_Type",[20,0])."='nodeComponent:".$component->{'class'}.":".$component->{'name'}."'\n";
	$functionCode .= "    return\n";
	$functionCode .= "  end function Node_Component_".ucfirst($componentName)."_Type\n\n";
	# Insert into the function list.
	push(
	    @{$buildData->{'code'}->{'functions'}},
	    $functionCode
	    );	
	# Bind this function to the relevant type.
	push(
	    @{$buildData->{'types'}->{'nodeComponent'.ucfirst($componentName)}->{'boundFunctions'}},
	    {type => "procedure", name => "type", function => "Node_Component_".ucfirst($componentName)."_Type"}
	    );
    }
}

sub Generate_Component_Assignment_Function {
    # Generate a type name functions.
    my $buildData = shift;
    # Specify data content.
    my @dataContent =
	(
	 {
	     intrinsic  => "class",
	     type       => "nodeComponent",
	     attributes => [ "intent(  out)" ],
	     variables  => [ "to" ]
	 },
	 {
	     intrinsic  => "class",
	     type       => "nodeComponent",
	     attributes => [ "intent(in   )" ],
	     variables  => [ "from" ]
	 }
	);
    # Generate the function code.
    my $functionCode;
    # Component assignment functions.
    $functionCode  = "  subroutine Node_Component_Assign(to,from)\n";
    $functionCode .= "    !% Assign a node component to another node component.\n";
    $functionCode .= "    implicit none\n";
    $functionCode .= &Format_Variable_Defintions(\@dataContent)."\n";
    $functionCode .= "    to%hostNode => from%hostNode\n";
    $functionCode .= "    select type (to)\n";
    foreach my $componentName ( @{$buildData->{'componentIdList'}} ) {
	my $component = $buildData->{'components'}->{$componentName};
	$functionCode .= "    type is (nodeComponent".padFullyQualified($componentName,[0,0]).")\n";
	$functionCode .= "       select type (from)\n";
	$functionCode .= "       type is (nodeComponent".padFullyQualified($componentName,[0,0]).")\n";
	foreach my $methodName ( keys(%{$component->{'methods'}->{'method'}}) ) {
	    my $method = $component->{'methods'}->{'method'}->{$methodName};
	    if ( exists($method->{'linkedData'}) ) {
		my $linkedDataName = $method->{'linkedData'};		
		$functionCode .= "          to%".padLinkedData($linkedDataName,[0,0])."%value=from%".padLinkedData($linkedDataName,[0,0])."%value\n";
	    }
	}
	$functionCode .= "       end select\n";
    }
    $functionCode .= "    end select\n";
    $functionCode .= "    return\n";
    $functionCode .= "  end subroutine Node_Component_Assign\n\n";
    # Insert into the function list.
    push(
	@{$buildData->{'code'}->{'functions'}},
	$functionCode
	);	
    # Bind this function to the treeNode type.
    push(
	@{$buildData->{'types'}->{'nodeComponent'}->{'boundFunctions'}},
	{type => "procedure", name => "assign"       , function => "Node_Component_Assign"},
	{type => "generic"  , name => "assignment(=)", function => "assign"               }
	)
	unless ( $workaround == 1 );
}

sub Generate_Component_Class_Destruction_Functions {
    # Generate class destruction functions.
    my $buildData = shift;
    # Initialize data content.
    my @dataContent;
    # Generate the function code.
    my $functionCode;
    foreach my $componentClassName ( @{$buildData->{'componentClassList'}} ) {
	# Specify data content.
	@dataContent =
	    (
	     {
		 intrinsic  => "class",
		 type       => "nodeComponent".ucfirst($componentClassName),
		 attributes => [ "intent(inout)" ],
		 variables  => [ "self" ]
	     }
	    );
	# Generate the function code.
	$functionCode  = "  subroutine Node_Component_".ucfirst($componentClassName)."_Destroy(self)\n";
	$functionCode .= "    !% Destroys a ".$componentClassName." component.\n";
	$functionCode .= "    implicit none\n";
	$functionCode .= &Format_Variable_Defintions(\@dataContent)."\n";
	$functionCode .= "    ! Do nothing.\n";
	$functionCode .= "    return\n";
	$functionCode .= "  end subroutine Node_Component_".ucfirst($componentClassName)."_Destroy\n\n";
	# Insert into the function list.
	push(
	    @{$buildData->{'code'}->{'functions'}},
	    $functionCode
	    );
	# Bind this function to the treeNode type.
	push(
	    @{$buildData->{'types'}->{'nodeComponent'.ucfirst($componentClassName)}->{'boundFunctions'}},
	    {type => "procedure", name => "destroy", function => "Node_Component_".ucfirst($componentClassName)."_Destroy"}
	    );
    }
}

sub Generate_Component_Class_Removal_Functions {
    # Generate class removal functions.
    my $buildData = shift;

    # Initialize data content.
    my @dataContent;
    # Generate the function code.
    my $functionCode;
    foreach my $componentClassName ( @{$buildData->{'componentClassList'}} ) {
	# Specify data content.
	@dataContent =
	    (
	     {
		 intrinsic  => "class",
		 type       => "treeNode",
		 attributes => [ "intent(inout)" ],
		 variables  => [ "self" ]
	     },
	     {
		 intrinsic  => "integer",
		 attributes => [ "intent(in   )" ],
		 variables  => [ "instance" ]
	     },
	     {
		 intrinsic  => "integer",
		 variables  => [ "instanceCount" ]
	     },
	     {
		 intrinsic  => "class",
		 type       => "nodeComponent".ucfirst($componentClassName),
		 attributes => [ "allocatable, dimension(:)" ],
		 variables  => [ "instancesTemporary" ]
	     }
	    );
	# Generate the function code.
	$functionCode  = "  subroutine Node_Component_".ucfirst($componentClassName)."_Remove(self,instance)\n";
	$functionCode .= "    !% Removes an instance of the ".$componentClassName." component, shifting other instances to keep the array contiguous.\n";
	$functionCode .= "    use Galacticus_Error\n";
	$functionCode .= "    implicit none\n";
	$functionCode .= &Format_Variable_Defintions(\@dataContent)."\n";
	$functionCode .= "    instanceCount=self%".$componentClassName."count()\n";
	$functionCode .= "    if (instance < 1 .or. instance > instanceCount) call Galacticus_Error_Report('Node_Component_".ucfirst($componentClassName)."_Remove','instance out of range')\n";
	$functionCode .= "    call self%component".ucfirst($componentClassName)."(instance)%destroy()\n";
	$functionCode .= "    if (instanceCount == 1) then\n";
	$functionCode .= "      ! Only one instance of this component. Deallocate it and reallocate with generic type.\n";
	$functionCode .= "      deallocate(self%component".ucfirst($componentClassName).")\n";
	$functionCode .= "      allocate(self%component".ucfirst($componentClassName)."(1))\n";
	$functionCode .= "    else\n";
	$functionCode .= "      ! Multiple instances, so remove the specified instance.\n";
	$functionCode .= "      allocate(instancesTemporary(instanceCount-1),source=self%component".ucfirst($componentClassName).")\n";
	if ( $workaround == 1 ) {
	    foreach my $implementationName ( @{$buildData->{'componentClasses'}->{$componentClassName}->{'members'}} ) {
		$functionCode .= "      select type (from => self%component".ucfirst($componentClassName).")\n";
		$functionCode .= "      type is (nodeComponent".padFullyQualified(ucfirst($componentClassName).ucfirst($implementationName),[0,0]).")\n";
		$functionCode .= "        select type (to => instancesTemporary)\n";
		$functionCode .= "        type is (nodeComponent".padFullyQualified(ucfirst($componentClassName).ucfirst($implementationName),[0,0]).")\n";
		$functionCode .= "          if (instance >             1) to(       1:instance     -1)=from(         1:instance     -1)\n";
		$functionCode .= "          if (instance < instanceCount) to(instance:instanceCount-1)=from(instance+1:instanceCount  )\n";
		$functionCode .= "        end select\n";
		$functionCode .= "      end select\n";
	    }
	} else {
	    $functionCode .= "      if (instance >             1) instancesTemporary(       1:instance     -1)=self%component".ucfirst($componentClassName)."(         1:instance     -1)\n";
	    $functionCode .= "      if (instance < instanceCount) instancesTemporary(instance:instanceCount-1)=self%component".ucfirst($componentClassName)."(instance+1:instanceCount  )\n";
	}
	$functionCode .= "      deallocate(self%component".ucfirst($componentClassName).")\n";
	$functionCode .= "      call Move_Alloc(instancesTemporary,self%component".ucfirst($componentClassName).")\n";
	
	$functionCode .= "    end if\n";
	$functionCode .= "    return\n";
	$functionCode .= "  end subroutine Node_Component_".ucfirst($componentClassName)."_Remove\n";
	# Insert into the function list.
	push(
	    @{$buildData->{'code'}->{'functions'}},
	    $functionCode
	    );
	# Bind this function to the treeNode type.
	push(
	    @{$buildData->{'types'}->{'treeNode'}->{'boundFunctions'}},
	    {type => "procedure", name => $componentClassName."Remove", function => "Node_Component_".ucfirst($componentClassName)."_Remove"}
	    );
    }
}

sub Generate_Component_Class_Move_Functions {
    # Generate class move functions.
    my $buildData = shift;
    # Initialize data content.
    my @dataContent;
    # Generate the function code.
    my $functionCode;
    foreach my $componentClassName ( @{$buildData->{'componentClassList'}} ) {
	# Specify data content.
	@dataContent =
	    (
	     {
		 intrinsic  => "class",
		 type       => "treeNode",
		 attributes => [ "intent(inout)" ],
		 variables  => [ "self" ]
	     },
	     {
		 intrinsic  => "type",
		 type       => "treeNode",
		 attributes => [ "intent(inout)", "target" ],
		 variables  => [ "targetNode" ]
	     },
	     {
		 intrinsic  => "integer",
		 variables  => [ "instanceCount", "targetCount", "i" ]
	     },
	     {
		 intrinsic  => "class",
		 type       => "nodeComponent".ucfirst($componentClassName),
		 attributes => [ "allocatable, dimension(:)" ],
		 variables  => [ "instancesTemporary" ]
	     }
	    );
	# Generate the function code.
	$functionCode  = "  subroutine Node_Component_".ucfirst($componentClassName)."_Move(self,targetNode)\n";
	$functionCode .= "    !% Move instances of the ".$componentClassName." component, from one node to another.\n";
	$functionCode .= "    use Galacticus_Error\n";
	$functionCode .= "    implicit none\n";
	$functionCode .= &Format_Variable_Defintions(\@dataContent)."\n";
	$functionCode .= "    instanceCount=self      %".$componentClassName."count()\n";
	$functionCode .= "    targetCount  =targetNode%".$componentClassName."count()\n";
	$functionCode .= "    if (instanceCount == 0) return\n";
	$functionCode .= "    if (targetCount == 0) then\n";
	$functionCode .= "      deallocate(targetNode%component".ucfirst($componentClassName).")\n";
	$functionCode .= "      call Move_Alloc(self%component".ucfirst($componentClassName).",targetNode%component".ucfirst($componentClassName).")\n";
	$functionCode .= "    else\n";
	$functionCode .= "      ! Multiple instances, so remove the specified instance.\n";
	$functionCode .= "      allocate(instancesTemporary(instanceCount+targetCount),source=self%component".ucfirst($componentClassName)."(1))\n";
	if ( $workaround == 1 ) {
	    foreach my $implementationName ( @{$buildData->{'componentClasses'}->{$componentClassName}->{'members'}} ) {
		$functionCode .= "      select type (from => targetNode%component".ucfirst($componentClassName).")\n";
		$functionCode .= "      type is (nodeComponent".padFullyQualified(ucfirst($componentClassName).ucfirst($implementationName),[0,0]).")\n";
		$functionCode .= "        select type (to => instancesTemporary)\n";
		$functionCode .= "        type is (nodeComponent".padFullyQualified(ucfirst($componentClassName).ucfirst($implementationName),[0,0]).")\n";
		$functionCode .= "          to(1:targetCount)=from\n";
		$functionCode .= "        end select\n";
		$functionCode .= "      end select\n";
	    }
	    foreach my $implementationName ( @{$buildData->{'componentClasses'}->{$componentClassName}->{'members'}} ) {
		$functionCode .= "      select type (from => self%component".ucfirst($componentClassName).")\n";
		$functionCode .= "      type is (nodeComponent".padFullyQualified(ucfirst($componentClassName).ucfirst($implementationName),[0,0]).")\n";
		$functionCode .= "        select type (to => instancesTemporary)\n";
		$functionCode .= "        type is (nodeComponent".padFullyQualified(ucfirst($componentClassName).ucfirst($implementationName),[0,0]).")\n";
		$functionCode .= "          to(targetCount+1:targetCount+instanceCount)=from\n";
		$functionCode .= "        end select\n";
		$functionCode .= "      end select\n";
	    }
	} else {
	    $functionCode .= "      instancesTemporary(            1:targetCount              )=targetNode%component".ucfirst($componentClassName)."\n";
	    $functionCode .= "      instancesTemporary(targetCount+1:targetCount+instanceCount)=self      %component".ucfirst($componentClassName)."\n";
	}
	$functionCode .= "      call targetNode%".$componentClassName."Destroy()\n";
	$functionCode .= "      call self      %".$componentClassName."Destroy()\n";
	$functionCode .= "      call Move_Alloc(instancesTemporary,targetNode%component".ucfirst($componentClassName).")\n";
	$functionCode .= "      allocate(self%component".ucfirst($componentClassName)."(1))\n";
	$functionCode .= "    end if\n";
	$functionCode .= "    do i=1,size(targetNode%component".ucfirst($componentClassName).")\n";
	$functionCode .= "       targetNode%component".ucfirst($componentClassName)."(i)%hostNode => targetNode\n";
	$functionCode .= "    end do\n";
	$functionCode .= "    return\n";
	$functionCode .= "  end subroutine Node_Component_".ucfirst($componentClassName)."_Move\n";
	# Insert into the function list.
	push(
	    @{$buildData->{'code'}->{'functions'}},
	    $functionCode
	    );
	# Bind this function to the treeNode type.
	push(
	    @{$buildData->{'types'}->{'treeNode'}->{'boundFunctions'}},
	    {type => "procedure", name => $componentClassName."Move", function => "Node_Component_".ucfirst($componentClassName)."_Move"}
	    );
    }
}

sub Generate_Component_Class_Dump_Functions {
    # Generate dump for each component class.
    my $buildData = shift;
    # Iterate over component classes.
    foreach my $componentClassName ( @{$buildData->{'componentClassList'}} ) {
	# Initialize function code.
	my $functionCode;
	# Initialize data content.
	my @dataContent =
	    (
	     {
		 intrinsic  => "class",
		 type       => "nodeComponent".ucfirst($componentClassName),
		 attributes => [ "intent(in   )" ],
		 variables  => [ "self" ]
	     }
	    );
	# Generate dump function.
	$functionCode  = "  subroutine Node_Component_".ucfirst($componentClassName)."_Dump(self)\n";
	$functionCode .= "    !% Dump the contents of a generic ".$componentClassName." component.\n";
	$functionCode .= "    use Galacticus_Display\n";
	$functionCode .= "    use ISO_Varying_String\n";
	$functionCode .= "    implicit none\n";
	$functionCode .= &Format_Variable_Defintions(\@dataContent)."\n";
	$functionCode .= "    call Galacticus_Display_Indent('".$componentClassName.": ".(" " x ($fullyQualifiedNameLengthMax-length($componentClassName)))."generic')\n";
	$functionCode .= "    call Galacticus_Display_Unindent('done')\n";
	$functionCode .= "    return\n";
	$functionCode .= "  end subroutine Node_Component_".ucfirst($componentClassName)."_Dump\n";
	# Insert into the function list.
	push(
	    @{$buildData->{'code'}->{'functions'}},
	    $functionCode
	    );
	# Insert a type-binding for this function into the implementation type.
	push(
	    @{$buildData->{'types'}->{'nodeComponent'.ucfirst($componentClassName)}->{'boundFunctions'}},
	    {type => "procedure", name => "dump", function => "Node_Component_".ucfirst($componentClassName)."_Dump"},
	    );
    }
}

sub Generate_Component_Class_Output_Functions {
    # Generate output for each component class.
    my $buildData = shift;
    # Iterate over component classes.
    foreach my $componentClassName ( @{$buildData->{'componentClassList'}} ) {
	# Initialize function code.
	my $functionCode;
	# Create property count function.
	my @dataContent =
	    (
	     {
		 intrinsic  => "class",
		 type       => "nodeComponent".ucfirst($componentClassName),
		 attributes => [ "intent(inout)" ],
		 variables  => [ "self" ]
	     },
	     {
		 intrinsic  => "integer",
		 attributes => [ "intent(inout)" ],
		 variables  => [ "integerPropertyCount", "doublePropertyCount" ]
	     },
	     {
		 intrinsic  => "double precision",
		 attributes => [ "intent(in   )" ],
		 variables  => [ "time" ]
	     },
	     {
		 intrinsic  => "integer",
		 attributes => ["intent(in   )" ],
		 variables  => [ "instance" ]
	     }
	    );
	undef($functionCode);
	$functionCode  = "  subroutine Node_Component_".ucfirst($componentClassName)."_Output_Count(self,integerPropertyCount,doublePropertyCount,time,instance)\n";
	$functionCode .= "    !% Increment the count of properties to output for a generic ".$componentClassName." component.\n";
	$functionCode .= "    implicit none\n";
	$functionCode .= &Format_Variable_Defintions(\@dataContent)."\n";
	$functionCode .= "    call default".ucfirst($componentClassName)."Component%outputCount(integerPropertyCount,doublePropertyCount,time,instance)\n";
	$functionCode .= "    return\n";
	$functionCode .= "  end subroutine Node_Component_".ucfirst($componentClassName)."_Output_Count\n";
	# Insert into the function list.
	push(
	    @{$buildData->{'code'}->{'functions'}},
	    $functionCode
	    );
	# Insert a type-binding for this function into the implementation type.
	push(
	    @{$buildData->{'types'}->{'nodeComponent'.ucfirst($componentClassName)}->{'boundFunctions'}},
	    {type => "procedure", name => "outputCount", function => "Node_Component_".ucfirst($componentClassName)."_Output_Count"},
	    );
	# Create property names function.
	@dataContent =
	    (
	     {
		 intrinsic  => "class",
		 type       => "nodeComponent".ucfirst($componentClassName),
		 attributes => [ "intent(inout)" ],
		 variables  => [ "self" ]
	     },
	     {
		 intrinsic  => "double precision",
		 attributes => [ "intent(in   )" ], 
		 variables  => [ "time" ]
	     },
	     {
		 intrinsic  => "integer", 
		 attributes => [ "intent(inout)" ], 
		 variables  => [ "integerProperty", "doubleProperty" ]
	     },
	     {
		 intrinsic  => "character",
		 type       => "len=*",
		 attributes => [ "intent(inout)", "dimension(:)" ], 
		 variables  => [ "integerPropertyNames", "integerPropertyComments", "doublePropertyNames", "doublePropertyComments" ]
	     },
	     {
		 intrinsic  => "double precision",
		 attributes => [ "intent(inout)", "dimension(:)" ],
		 variables  => [ "integerPropertyUnitsSI", "doublePropertyUnitsSI" ]
	     },
	     {
		 intrinsic  => "integer",
		 attributes => ["intent(in   )" ],
		 variables  => [ "instance" ]
	     }
	    );
	undef($functionCode);
	$functionCode  = "  subroutine Node_Component_".ucfirst($componentClassName)."_Output_Names(self,integerProperty,integerPropertyNames,integerPropertyComments,integerPropertyUnitsSI,doubleProperty,doublePropertyNames,doublePropertyComments,doublePropertyUnitsSI,time,instance)\n";
	$functionCode .= "    !% Establish property names for a generic ".$componentClassName." component.\n";
	$functionCode .= "    implicit none\n";
	$functionCode .= &Format_Variable_Defintions(\@dataContent)."\n";
	$functionCode .= "    call default".ucfirst($componentClassName)."Component%outputNames(integerProperty,integerPropertyNames,integerPropertyComments,integerPropertyUnitsSI,doubleProperty,doublePropertyNames,doublePropertyComments,doublePropertyUnitsSI,time,instance)\n";
	$functionCode .= "    return\n";
	$functionCode .= "  end subroutine Node_Component_".ucfirst($componentClassName)."_Output_Names\n";
	# Insert into the function list.
	push(
	    @{$buildData->{'code'}->{'functions'}},
	    $functionCode
	    );
	# Insert a type-binding for this function into the implementation type.
	push(
	    @{$buildData->{'types'}->{'nodeComponent'.ucfirst($componentClassName)}->{'boundFunctions'}},
	    {type => "procedure", name => "outputNames", function => "Node_Component_".ucfirst($componentClassName)."_Output_Names"},
	    );
	# Create output function.
	@dataContent =
	    (
	     {
		 intrinsic  => "class",
		 type       => "nodeComponent".ucfirst($componentClassName),
		 attributes => [ "intent(inout)" ],
		 variables  => [ "self" ]
	     },
	     {
		 intrinsic  => "double precision",
		 attributes => [ "intent(in   )" ], 
		 variables  => [ "time" ]
	     },
	     {
		 intrinsic  => "integer", 
		 attributes => [ "intent(inout)" ], 
		 variables  => [ "integerProperty", "integerBufferCount", "doubleProperty", "doubleBufferCount" ]
	     },
	     {
		 intrinsic  => "integer",
		 type       => "kind=kind_int8",
		 attributes => [ "intent(inout)", "dimension(:,:)" ],
		 variables  => [ "integerBuffer" ]
	     },
	     {
		 intrinsic  => "double precision",
		 attributes => [ "intent(inout)", "dimension(:,:)" ],
		 variables  => [ "doubleBuffer" ]
	     },
	     {
		 intrinsic  => "integer",
		 type       => "kind=kind_int8",
		 attributes => [ "allocatable", "dimension(:)" ],
		 variables  => [ "outputRank1Integer" ]
	     },
	     {
		 intrinsic  => "double precision",
		 attributes => [ "allocatable", "dimension(:)" ],
		 variables  => [ "outputRank1Double" ]
	     },
	     {
		 intrinsic  => "integer",
		 attributes => ["intent(in   )" ],
		 variables  => [ "instance" ]
	     },
	     {
		 intrinsic  => "integer",
		 variables  => [ "i" ]
	     }
	    );
	# Find all derived types to be output.
	my %outputTypes;
	foreach my $componentName ( @{$buildData->{'componentClasses'}->{$componentClassName}->{'members'}} ) {
	    # Get the component.
	    my $componentID  = ucfirst($componentClassName).ucfirst($componentName);
	    my $component    = $buildData->{'components'}->{$componentID};
	    foreach my $methodName ( keys(%{$component->{'methods'}->{'method'}}) ) {
		my $method = $component->{'methods'}->{'method'}->{$methodName};
		# Check if property is to be output.
		if ( exists($method->{'output'}) ) {
		    # Get the type of this component.
		    my $type;
		    if ( exists($method->{'linkedData'}) ) {
			my $linkedDataName = $method->{'linkedData'};
			my $linkedData     = $component->{'content'}->{'data'}->{$linkedDataName};
			$type   = $linkedData->{'type'};
		    } else {
			$type = $method->{'type'};		
		    }
		    $outputTypes{$type} = 1
			unless ( $type eq "real" || $type eq "integer" );
		}
	    }
	}
	my @outputTypes;
	foreach ( keys(%outputTypes) ){
	    push(
		@outputTypes,
		{
		    intrinsic => "type",
		    type      => $_,
		    variables => [ "output".ucfirst($_) ]
		}
		);
	}
	push(@dataContent,@outputTypes);
	undef($functionCode);
	# Find modules required.
	my %modulesRequired;
	foreach my $componentName ( @{$buildData->{'componentClasses'}->{$componentClassName}->{'members'}} ) {
	    # Get the component.
	    my $componentID  = ucfirst($componentClassName).ucfirst($componentName);
	    my $component    = $buildData->{'components'}->{$componentID};
	    foreach my $methodName ( keys(%{$component->{'methods'}->{'method'}}) ) {
		my $method = $component->{'methods'}->{'method'}->{$methodName};
		# Check if this method is to be output.
		if ( exists($method->{'output'}) ) {
		    if ( exists($method->{'output'}->{'modules'}) ) {
			my $moduleList = $method->{'output'}->{'modules'};
			$moduleList =~ s/^\s*//;
			$moduleList =~ s/\s*$//;
			my @modules = split(/\s*,\s*/,$moduleList);
			foreach ( @modules ) {
			    $modulesRequired{$_} = 1;
			}
		    }
		}
	    }
	}
	$functionCode  = "  subroutine Node_Component_".ucfirst($componentClassName)."_Output(self,integerProperty,integerBufferCount,integerBuffer,doubleProperty,doubleBufferCount,doubleBuffer,time,instance)\n";
	$functionCode .= "    !% Output ptoperties for a ".$componentClassName." component.\n";
	$functionCode .= "    use ".$_."\n" 
	    foreach ( keys(%modulesRequired) );
	$functionCode .= "    implicit none\n";
	$functionCode .= &Format_Variable_Defintions(\@dataContent)."\n";
	my %typeMap =
	    (
	     real    => "double",
	     integer => "integer"
	    );
	foreach my $componentName ( @{$buildData->{'componentClasses'}->{$componentClassName}->{'members'}} ) {
	    # Get the component.
	    my $componentID  = ucfirst($componentClassName).ucfirst($componentName);
	    my $component    = $buildData->{'components'}->{$componentID};
	    my $activeCheck  = "    if (default".ucfirst($componentClassName)."Component%".$componentName."IsActive()";
	    $activeCheck .= ".and.instance == 1"
		if (
		    exists($component->{'output'}               )           &&
		    exists($component->{'output'}->{'instances'})           &&
		           $component->{'output'}->{'instances'} eq "first"
		);
	    $activeCheck .= ") then\n";
	    my $outputsFound = 0;
	    foreach my $methodName ( keys(%{$component->{'methods'}->{'method'}}) ) {
		my $method = $component->{'methods'}->{'method'}->{$methodName};
		# Check if this method is to be output.
		if ( exists($method->{'output'}) ) {
		    # Add conditional statement if necessary.
		    if ( $outputsFound == 0 ) {
			$functionCode .= $activeCheck;
			$outputsFound  = 1;
		    }
		    # Define rank, type and value.
		    my $rank;
		    my $type;
		    # Check if this method has any linked data in this component.
		    if ( exists($method->{'linkedData'}) ) {
			my $linkedDataName = $method->{'linkedData'};
			my $linkedData     = $component->{'content'}->{'data'}->{$linkedDataName};
			$rank   = $linkedData->{'rank'};
			$type   = $linkedData->{'type'};
		    } elsif ( $method->{'isVirtual'} eq "yes" && $method->{'attributes'}->{'isGettable'} eq "true" ) {
			$rank = $method->{'rank'};
			$type = $method->{'type'};
		    } else {
			die("Generate_Component_Class_Output_Functions(): can not output [".$methodName."]");
		    }		   
		    # Determine count.
		    my $count;
		    if ( $rank == 0 ) {
			$count = 1;
		    } elsif ( $rank ==1 ) {
			die("Generate_Component_Class_Output_Functions(): output of rank>0 objects requires a labels attribute")
			    unless ( exists($method->{'output'}->{'labels'}) );	
			if ( $method->{'output'}->{'labels'} =~ m/^\[(.*)\]$/ ) {
			    my $labelText = $1;
			    $labelText    =~ s/^\s*//;
			    $labelText    =~ s/\s*$//;
			    my @labels    = split(/\s*,\s*/,$labelText);
			    $count = scalar(@labels);
			} elsif ( exists($method->{'output'}->{'count'}) ) {
			    $count = $method->{'output'}->{'count'};
			} else {
			    die('Generate_Component_Class_Output_Functions(): no method to determine output size for rank-1 property');
			}
		    } else {
			die("Generate_Component_Class_Output_Functions(): output of rank>1 arrays not supported");
		    }
		    # Increment the counters.
		    switch ( $type ) {
			case ( [ "real", "integer" ] ) {
			    if ( $rank == 0 ) {
				if ( exists($method->{'output'}->{'condition'}) ) {
				    my $condition = $method->{'output'}->{'condition'};
				    $condition =~ s/\[\[([^\]]+)\]\]/$1/g;
				    $functionCode .= "    if (".$condition.") then\n";
				}
				$functionCode .= "       ".$typeMap{$type}."Property=".$typeMap{$type}."Property+1\n";
				$functionCode .= "       ".$typeMap{$type}."Buffer(".$typeMap{$type}."BufferCount,".$typeMap{$type}."Property)=self%".$methodName."()\n";
				$functionCode .= "    end if\n"
				    if ( exists($method->{'output'}->{'condition'}) );
			    } else {
				if ( exists($method->{'output'}->{'condition'}) ) {
				    die("Generate_Component_Class_Output_Functions(): conditions for rank>1 properties not supported")
					unless ( $rank == 1 );
				    my $condition = $method->{'output'}->{'condition'};
				    $condition =~ s/\[\[([^\]]+)\]\]/$1/g;
				    $condition =~ s/\{i\}/i/g;
				    $functionCode .= "    outputRank1".ucfirst($typeMap{$type})."=self%".$methodName."()\n";
				    $functionCode .= "    do i=1,".$method->{'output'}->{'count'}."\n";
				    $functionCode .= "      if (".$condition.") then\n";
				    $functionCode .= "        ".$typeMap{$type}."Property=".$typeMap{$type}."Property+1\n";
				    $functionCode .= "        ".$typeMap{$type}."Buffer(".$typeMap{$type}."BufferCount,".$typeMap{$type}."Property)=outputRank1".ucfirst($typeMap{$type})."(i)\n";
				    $functionCode .= "      end if\n";
				    $functionCode .= "    end do\n";
				    $functionCode .= "    deallocate(outputRank1".ucfirst($typeMap{$type}).")\n";
				} else {
				    $functionCode .= "       ".$typeMap{$type}."Buffer(".$typeMap{$type}."BufferCount,".$typeMap{$type}."Property+1:".$typeMap{$type}."Property+".$count.")=reshape(self%".$methodName."(),[".$count."])\n";
				    $functionCode .= "       ".$typeMap{$type}."Property=".$typeMap{$type}."Property+".$count."\n";
				}
			    }
			}
			else {
			    $functionCode .= "      output".ucfirst($type)."=self%".$methodName."()\n";
			    $functionCode .= "      call output".ucfirst($type)."%output(integerProperty,integerBufferCount,integerBuffer,doubleProperty,doubleBufferCount,doubleBuffer,time)\n";
			}
		    }
		}
	    }
	    $functionCode .= "    end if\n"
		if ( $outputsFound == 1 );
	}
	$functionCode .= "    return\n";
	$functionCode .= "  end subroutine Node_Component_".ucfirst($componentClassName)."_Output\n";
	# Insert into the function list.
	push(
	    @{$buildData->{'code'}->{'functions'}},
	    $functionCode
	    );
	# Insert a type-binding for this function into the implementation type.
	push(
	    @{$buildData->{'types'}->{'nodeComponent'.ucfirst($componentClassName)}->{'boundFunctions'}},
	    {type => "procedure", name => "output", function => "Node_Component_".ucfirst($componentClassName)."_Output"},
	    );
    }
}

sub Generate_Is_Active_Functions {
    # Generate "isActive" functions.
    my $buildData = shift;
    # Initialize function code.
    my $functionCode;
    # Iterate over component implementations.
    foreach my $componentID ( @{$buildData->{'componentIdList'}} ) {
	my $component = $buildData->{'components'}->{$componentID};
	$functionCode  = "  logical function Node_Component_".ucfirst($componentID)."_Is_Active()\n";
	$functionCode .= "    !% Return true if the ".$component->{'name'}." implementation of the ".$component->{'class'}." component is the active choice.\n";
	$functionCode .= "    implicit none\n\n";
	$functionCode .= "    Node_Component_".ucfirst($componentID)."_Is_Active=nodeComponent".ucfirst($componentID)."IsActive\n";
	$functionCode .= "    return\n";
	$functionCode .= "  end function Node_Component_".ucfirst($componentID)."_Is_Active\n";
	# Insert into the function list.
	push(
	    @{$buildData->{'code'}->{'functions'}},
	    $functionCode
	    );
	# Bind this function to the implementation type.
	push(
	    @{$buildData->{'types'}->{'nodeComponent'.ucfirst($component->{'class'})}->{'boundFunctions'}},
	    {type => "procedure", pass => "nopass", name => lcfirst($component->{'name'})."IsActive", function => "Node_Component_".ucfirst($componentID)."_Is_Active"}
	    );
    }
}

sub Generate_Component_Implementation_Destruction_Functions {
    # Generate component implementation destruction functions.
    my $buildData = shift;
    # Initialize function code.
    my $functionCode;
    foreach my $componentID ( @{$buildData->{'componentIdList'}} ) {
	my $component = $buildData->{'components'}->{$componentID};
	# Specify data content.
	my @dataContent =
	    (
	     {
		 intrinsic  => "class",
		 type       => "nodeComponent".ucfirst($componentID),
		 attributes => [ "intent(inout)" ],
		 variables  => [ "self" ]
	     }
	    );
	# Generate function code.
  	$functionCode  = "  subroutine Node_Component_".ucfirst($componentID)."_Destroy(self)\n";
	$functionCode .= "    !% Destroy a ".$component->{'name'}." implementation of the ".$component->{'class'}." component.\n";
	$functionCode .= "    use Memory_Management\n";
	$functionCode .= "    implicit none\n";
	$functionCode .= &Format_Variable_Defintions(\@dataContent)."\n";
	# Iterate over methods.
	foreach my $methodName ( keys(%{$component->{'methods'}->{'method'}}) ) {
	    my $method = $component->{'methods'}->{'method'}->{$methodName};
   	    # Check if this method has any linked data in this component.
	    if ( exists($method->{'linkedData'}) ) {
		# For each linked datum deallocate if necessary.
		my $linkedDataName = $method->{'linkedData'};
		my $linkedData     = $component->{'content'}->{'data'}->{$linkedDataName};
		switch ( $linkedData->{'type'} ) {
		    case ( "real"    ) {
			# Nothing to do in this case.
		    }
		    case ( "integer" ) {
			# Nothing to do in this case.
		    }
		    case ( "logical" ) {
			# Nothing to do in this case.
		    }
		    else {
			my @contents = ( "value" );
			push(@contents,"rate ","scale")
			    if ( $method->{'attributes'}->{'isEvolvable'} eq "true" );
			$functionCode .= "    call self%".padLinkedData($linkedDataName,[0,0])."%".$_."%destroy()\n"
			    foreach ( @contents );
		    }
		}
		if ( $linkedData->{'rank'} > 0 ) {
		    my @contents = ( "value" );
		    push(@contents,"rate ","scale")
			if ( $method->{'attributes'}->{'isEvolvable'} eq "true" );
		    $functionCode .= "    if (allocated(self%".padLinkedData($linkedDataName,[0,0])."%".$_.")) call Dealloc_Array(self%".padLinkedData($linkedDataName,[0,0])."%".$_.")\n"
			foreach ( @contents );
		}
	    }
	}
	$functionCode .= "    return\n";
	$functionCode .= "  end subroutine Node_Component_".ucfirst($componentID)."_Destroy\n\n";
	# Insert into the function list.
	push(
	    @{$buildData->{'code'}->{'functions'}},
	    $functionCode
	    );
	# Bind this function to the implementation type.
	push(
	    @{$buildData->{'types'}->{'nodeComponent'.ucfirst($componentID)}->{'boundFunctions'}},
	    {type => "procedure", name => "destroy", function => "Node_Component_".ucfirst($componentID)."_Destroy"}
	    );    
    }
}

sub Generate_ODE_Initialization_Functions {
    # Generate ODE solver initialization functions.
    my $buildData = shift;
    # Initialize function code.
    my $functionCode;
    # Initialize data content.
    my @dataContent;
    # Iterate over component implementations.
    foreach my $componentID ( @{$buildData->{'componentIdList'}} ) {
	my $component = $buildData->{'components'}->{$componentID};
	# Specify data content.
	@dataContent =
	    (
	     {
		 intrinsic  => "class",
		 type       => "nodeComponent".ucfirst($componentID),
		 attributes => [ "intent(inout)" ],
		 variables  => [ "self" ]
	     }
	    );
	# Generate rate initialization function code.
  	$functionCode  = "  subroutine Node_Component_".ucfirst($componentID)."_ODE_Step_Rates_Init(self)\n";
	$functionCode .= "    !% Initialize rates in a ".$component->{'name'}." implementation of the ".$component->{'class'}." component for an ODE solver step.\n";
	$functionCode .= "    implicit none\n";
	$functionCode .= &Format_Variable_Defintions(\@dataContent)."\n";
	foreach my $methodName ( keys(%{$component->{'methods'}->{'method'}}) ) {
	    my $method = $component->{'methods'}->{'method'}->{$methodName};	    
	    if ( exists($method->{'linkedData'}) ) {
		my $linkedDataName = $method->{'linkedData'};
		my $linkedData     = $component->{'content'}->{'data'}->{$linkedDataName};
		if ( $component->{'content'}->{'data'}->{$linkedDataName}->{'isEvolvable'} eq "true" ) {
		    switch ( $linkedData->{'type'} ) {
			case ( "real"    ) {
			    if ( $linkedData->{'rank'} == 0 ) {
				$functionCode .= "         self%".$linkedDataName."%rate =0.0d0\n";
			    } else {
				$functionCode .= "         if (allocated(self%".$linkedDataName."%rate)) self%".$linkedDataName."%rate =0.0d0\n";
			    }
			}
			case ( "integer" ) {
			    die "Build_Include_File.pl: integer data type should not be evolvable";
			}
			case ( "logical" ) {
			    die "Build_Include_File.pl: logical data type should not be evolvable";
			}
			else {
			    $functionCode .= "    call self%".$linkedDataName."%rate %reset()\n";
			}
		    }
		}
	    }
	}
	$functionCode .= "    return\n";
	$functionCode .= "  end subroutine Node_Component_".ucfirst($componentID)."_ODE_Step_Rates_Init\n\n";
	# Insert into the function list.
	push(
	    @{$buildData->{'code'}->{'functions'}},
	    $functionCode
	    );
	# Bind this function to the implementation type.
	push(
	    @{$buildData->{'types'}->{'nodeComponent'.ucfirst($componentID)}->{'boundFunctions'}},
	    {type => "procedure", name => "odeStepRatesInitialize" , function => "Node_Component_".ucfirst($componentID)."_ODE_Step_Rates_Init"}
	    );    
	# Generate scale initialization code.
  	$functionCode  = "  subroutine Node_Component_".ucfirst($componentID)."_ODE_Step_Scales_Init(self)\n";
	$functionCode .= "    !% Initialize scales in a ".$component->{'name'}." implementation of the ".$component->{'class'}." component for an ODE solver step.\n";
	$functionCode .= "    implicit none\n";
	$functionCode .= &Format_Variable_Defintions(\@dataContent)."\n";
	foreach my $methodName ( keys(%{$component->{'methods'}->{'method'}}) ) {
	    my $method = $component->{'methods'}->{'method'}->{$methodName};
	    if ( exists($method->{'linkedData'}) ) {
		my $linkedDataName = $method->{'linkedData'};
		my $linkedData     = $component->{'content'}->{'data'}->{$linkedDataName};
		if ( $component->{'content'}->{'data'}->{$linkedDataName}->{'isEvolvable'} eq "true" ) {
		    switch ( $linkedData->{'type'} ) {
			case ( "real"    ) {
			    if ( $linkedData->{'rank'} == 0 ) {
				$functionCode .= "         self%".$linkedDataName."%scale=1.0d0\n";
			    } else {
				$functionCode .= "         if (allocated(self%".$linkedDataName."%scale)) self%".$linkedDataName."%scale=1.0d0\n";
			    }
			}
			case ( "integer" ) {
			    die "Integer data type should not be evolvable";
			}
			case ( "logical" ) {
			    die "Logical data type should not be evolvable";
			}
			else {
			    $functionCode .= "    call self%".$linkedDataName."%scale%setToUnity()\n";
			}
		    }
		}
	    }
	}
	$functionCode .= "    end subroutine Node_Component_".ucfirst($componentID)."_ODE_Step_Scales_Init\n\n";
	# Insert into the function list.
	push(
	    @{$buildData->{'code'}->{'functions'}},
	    $functionCode
	    );
	# Bind this function to the implementation type.
	push(
	    @{$buildData->{'types'}->{'nodeComponent'.ucfirst($componentID)}->{'boundFunctions'}},
	    {type => "procedure", name => "odeStepScalesInitialize", function => "Node_Component_".ucfirst($componentID)."_ODE_Step_Scales_Init"}
	    );    
    }
}

sub Generate_Null_Binding_Functions {
    # Generate null binding functions.
    my $buildData = shift;
    # Iterate over component classes.
    foreach my $componentClassName ( keys(%{$buildData->{'nullMethods'}}) ) {
	# Get the null functions required for this component class.
	my $componentClass = $buildData->{'nullMethods'}->{$componentClassName};
	# Iterate over required null functions for this component class.
	foreach my $nullFunctionName ( keys(%{$componentClass}) ) {
	    # Get the null function definition.
	    my $nullFunction = $componentClass->{$nullFunctionName};
	    # Construct a datatype for this null function.
	    (my $dataDefinition, my $label) = &Data_Object_Definition($nullFunction,matchOnly => 1);
	    # Build a label describing the intrinsic type of the data.
	    my $intrinsicType = $dataDefinition->{'intrinsic'};
	    $intrinsicType .= $dataDefinition->{'type'}
	        if ( exists($dataDefinition->{'type'}) );
	    # Append rank to the label for this data.
	    $label .= $nullFunction->{'rank'};
	    # Extract the intent of the function.
	    my $intent = $nullFunction->{'intent'};
	    # Construct the type of "self" for this function.
	    my $selfType = "nodeComponent";
	    $selfType .= ucfirst($componentClassName)
		unless ( $componentClassName eq "generic" );
	    # Add an intent to the attributes of the datatype and specify the name of the input datatype.
	    push(@{$dataDefinition->{'attributes'}},"intent(in   )");
	    @{$dataDefinition->{'variables'}} = "setValue";
	    # Build code for the null set function.
	    my @dataContent =
		(
		 $dataDefinition,
		 {
		     intrinsic  => "class",
		     type       => $selfType,
		     attributes => [ "intent(".$intent.")" ],
		     variables  => [ "self" ]
		 }
		);
	    my $functionCode;
	    $functionCode  = "  subroutine ".$componentClassName."NullBindingSet".$label.$intent."(self,setValue)\n";
	    $functionCode .= "    !% A null set function for rank ".$nullFunction->{'rank'}." ".lc($intrinsicType)."s.\n";
	    $functionCode .= "    implicit none\n";
	    $functionCode .= &Format_Variable_Defintions(\@dataContent)."\n";
	    $functionCode .= "    return\n";
	    $functionCode .= "  end subroutine ".$componentClassName."NullBindingSet".$label.$intent."\n";
	    # Insert into the function list.
	    push(
		@{$buildData->{'code'}->{'functions'}},
		$functionCode
		);
	    # Build code for the null rate function.
	    @dataContent =
		(
		 $dataDefinition,
		 {
		     intrinsic  => "class",
		     type       => $selfType,
		     attributes => [ "intent(".$intent.")" ],
		     variables  => [ "self" ]
		 },
		 {
		     intrinsic  => "logical",
		     attributes => [ "intent(inout)", "optional" ],
		     variables  => [ "interrupt" ]
		 },
		 {
		     intrinsic  => "procedure",
		     type       => "", 
		     attributes => [ "intent(inout)", "optional", "pointer" ],
		     variables  => [ "interruptProcedure" ]
		 }
		);
	    $functionCode  = "  subroutine ".$componentClassName."NullBindingRate".$label.$intent."(self,setValue,interrupt,interruptProcedure)\n";
	    $functionCode .= "    !% A null rate function for rank ".$nullFunction->{'rank'}." ".lc($intrinsicType)."s.\n";
	    $functionCode .= "    implicit none\n";
	    $functionCode .= &Format_Variable_Defintions(\@dataContent)."\n";
	    $functionCode .= "    return\n";
	    $functionCode .= "  end subroutine ".$componentClassName."NullBindingRate".$label.$intent."\n";
	    # Insert into the function list.
	    push(
		@{$buildData->{'code'}->{'functions'}},
		$functionCode
		);
	    # Build code for the null get function.
	    pop(@{$dataDefinition->{'attributes'}});
	    push(@{$dataDefinition->{'attributes'}},"allocatable")
		if ( $nullFunction->{'rank'} > 0 );
	    @{$dataDefinition->{'variables'}} = $componentClassName."NullBinding".$label.$intent;
	     @dataContent =
		(
		 $dataDefinition,
		 {
		     intrinsic  => "class",
		     type       => $selfType,
		     attributes => [ "intent(".$intent.")" ],
		     variables  => [ "self" ]
		 },
	
		);
	    $functionCode  = "  function ".$componentClassName."NullBinding".$label.$intent."(self)\n";
	    $functionCode .= "    !% A null get function for rank ".$nullFunction->{'rank'}." ".lc($intrinsicType)."s.\n";
	    $functionCode .= "    implicit none\n";
	    $functionCode .= &Format_Variable_Defintions(\@dataContent)."\n";
	    $functionCode .= "    return\n";
	    $functionCode .= "  end function ".$componentClassName."NullBinding".$label.$intent."\n\n";
	    # Insert into the function list.
	    push(
		@{$buildData->{'code'}->{'functions'}},
		$functionCode
		);
	}
    }
}

sub Generate_Component_Class_Default_Value_Functions {
    # Generate component class default value functions.
    my $buildData = shift;
    # Iterate over component classes.
    foreach my $componentClassName ( @{$buildData->{'componentClassList'}} ) {
	# Initialize hash to track which method have been created already.
	my %methodsCreated;
	# Iterate over implementations in this class.
    	foreach my $componentName ( @{$buildData->{'componentClasses'}->{$componentClassName}->{'members'}} ) {
	    # Get the component.
	    my $componentID = ucfirst($componentClassName).ucfirst($componentName);
	    my $component   = $buildData->{'components'}->{$componentID};
	    # Iterate over the properties of this implementation.
	    foreach my $methodName ( keys(%{$component->{'methods'}->{'method'}}) ) {
		# Get the method.
		my $method = $component->{'methods'}->{'method'}->{$methodName};
		# Get the linked data.
		my $linkedData;
		if ( exists($method->{'linkedData'}) ) {
		    my $linkedDataName = $method->{'linkedData'};
		    $linkedData = $component->{'content'}->{'data'}->{$linkedDataName};
		} else {
		    $linkedData = $method;
		}
		# Specify required data content.
		(my $dataDefinition, my $label ) = &Data_Object_Definition($linkedData);
		push(@{$dataDefinition->{'variables' }},ucfirst($componentClassName).ucfirst($methodName));
		my @dataContent = (
		    $dataDefinition,
		    {
			intrinsic  => "class",
			type       => "nodeComponent".ucfirst($componentClassName),
			attributes => [ "intent(inout)" ],
			variables   => [ "self" ]
		    }
		    );
		# Skip if this method has already been created.
		unless ( exists($methodsCreated{$methodName}) ) {
		    # Generate code for "isGettable" function.
		    my $functionCode;
		    $functionCode  = "   logical function ".ucfirst($componentClassName).ucfirst($methodName)."IsGettable()\n";
		    $functionCode .= "     !% Returns true if the {\\tt ".$methodName."} method is gettable for the {\\tt ".$componentClassName."} component class.\n\n"; 
		    $functionCode .= "     implicit none\n";
		    $functionCode .= "     ".ucfirst($componentClassName).ucfirst($methodName)."IsGettable=.false.\n";
		    foreach my $componentName2 ( @{$buildData->{'componentClasses'}->{$componentClassName}->{'members'}} ) {
			my $component2ID = ucfirst($componentClassName).ucfirst($componentName2);
			my $component2   = $buildData->{'components'}->{$component2ID};
			$functionCode .= "     if (nodeComponent".ucfirst($component2ID)."IsActive) ".ucfirst($componentClassName).ucfirst($methodName)."IsGettable=.true.\n"
			    if (
				exists($component2->{'methods'}->{'method'}->{$methodName}                  ) && 
				exists($component2->{'methods'}->{'method'}->{$methodName}->{'classDefault'})
			    );
		    }
		    $functionCode .= "     return\n";
		    $functionCode .= "   end function ".ucfirst($componentClassName).ucfirst($methodName)."IsGettable\n";
		    # Insert into the function list.
		    push(
			@{$buildData->{'code'}->{'functions'}},
			$functionCode
			);
		    # Bind this function to the implementation type.
		    push(
			@{$buildData->{'types'}->{'nodeComponent'.ucfirst($componentClassName)}->{'boundFunctions'}},
			{type => "procedure", pass => "nopass", name => $methodName."IsGettable", function => ucfirst($componentClassName).ucfirst($methodName)."IsGettable"}
			);
		    # Generate code for default value function.
		    $functionCode  = "  function ".ucfirst($componentClassName).ucfirst($methodName)."(self)\n";
		    $functionCode .= "    !% Returns the default value for the {\\tt ".$methodName."} method for the {\\tt ".$componentClassName."} component class.\n";
		    # Insert any required modules.
		    if ( exists($method->{'classDefault'}) && exists($method->{'classDefault'}->{'modules'}) ) {
			foreach ( @{$method->{'classDefault'}->{'modules'}} ) {
			    $functionCode .= "    use ".$_."\n";
			}
		    }
		    $functionCode .= "    implicit none\n";
		    # Build default value code, and accumulate which additional components are needed.
		    my $defaultLines = "";
		    my %requiredComponents;
		    foreach my $componentName2 ( @{$buildData->{'componentClasses'}->{$componentClassName}->{'members'}} ) {
			my $component2ID = ucfirst($componentClassName).ucfirst($componentName2);
			my $component2   = $buildData->{'components'}->{$component2ID};
			if ( exists($component2->{'methods'}->{'method'}->{$methodName}) ) {
			    my $method2 = $component2->{'methods'}->{'method'}->{$methodName};
			    if ( exists($method2->{'classDefault'}) ) {
				$defaultLines .= "     if (nodeComponent".ucfirst($component2ID)."IsActive) then\n";
				my %selfComponents;
				my $default = $method2->{'classDefault'}->{'code'};
				while ( $default =~ m/self([a-zA-Z]+)Component\s*%/ ) {
				    $selfComponents{$1} = 1;
				    $requiredComponents{$1} = 1;
				    $default =~ s/self([a-zA-Z]+)Component\s*%//;
				}
				$defaultLines .= "    selfNode => self%host()\n" if ( scalar(keys(%selfComponents)) > 0 );
				foreach my $selfComponent ( keys(%selfComponents) ) {
				    $defaultLines .= "     self".$selfComponent."Component => selfNode%".lc($selfComponent)."()\n";
				}
				$defaultLines .= "       call Alloc_Array(".ucfirst($componentClassName).ucfirst($methodName).",[".$method2->{'classDefault'}->{'count'}."])\n"
				    if ( exists($method2->{'classDefault'}->{'count'}) );
				$defaultLines .= "       ".ucfirst($componentClassName).ucfirst($methodName)."=".$method2->{'classDefault'}->{'code'}."\n";
				$defaultLines .= "       return\n";
				$defaultLines .= "     end if\n";
			    }
			}
		    }
		    # Add a self node pointer if other components are required.
		    push(
			@dataContent,
			{
			    intrinsic  => "type",
			    type       => "treeNode",
			    attributes => [ "pointer" ],
			    variables  => [ "selfNode" ]
			}
			) if ( scalar(keys(%requiredComponents)) > 0 );
		    # Add pointers for each required component.
		    push(
			@dataContent,
			{
			    intrinsic  => "class",
			    type       => "nodeComponent".ucfirst($_),
			    attributes => [ "pointer" ],
			    variables  => [ "self".ucfirst($_)."Component" ]
			}
			)
			foreach ( keys(%requiredComponents) );
		    # Insert data content.
		    $functionCode .= &Format_Variable_Defintions(\@dataContent)."\n";
		    # Insert code to set required default.
		    $functionCode .= $defaultLines;
		    # Insert code to return zero values by default.
		    if ( $linkedData->{'rank'} == 0 ) {
			switch ( $linkedData->{'type'} ) {
			    case ( "real"    ) {
				$functionCode .= "    ".ucfirst($componentClassName).ucfirst($methodName)."=0.0d0\n";
			    }
			    case ( "integer" ) {
				$functionCode .= "    ".ucfirst($componentClassName).ucfirst($methodName)."=0\n";
			    }
			    case ( "logical" ) {
				$functionCode .= "    ".ucfirst($componentClassName).ucfirst($methodName)."=.false.\n";
			    }
			    else {
				$functionCode .= "     call ".ucfirst($componentClassName).ucfirst($methodName)."%reset()\n";
			    }
			}
		    } else {
			$functionCode .= "    ".ucfirst($componentClassName).ucfirst($methodName)."=null".$label.$linkedData->{'rank'}."d\n";
		    }
		    # Close the function.
		    $functionCode .= "    return\n";
		    $functionCode .= "  end function ".ucfirst($componentClassName).ucfirst($methodName)."\n";
		    # Insert into the function list.
		    push(
			@{$buildData->{'code'}->{'functions'}},
			$functionCode
			);
		    # Bind this function to the implementation type.
		    push(
			@{$buildData->{'types'}->{'nodeComponent'.ucfirst($componentClassName)}->{'boundFunctions'}},
			{type => "procedure", name => $methodName, function => ucfirst($componentClassName).ucfirst($methodName)}
			);
		    # Record that this method has been created.
		    $methodsCreated{$methodName} = 1;
		}
	    }
	}
    }
}

sub Bound_Function_Table {
    # Get the list of type-bound functions.
    my @typeBoundFunctions = @{$_[0]};
    # Create a text table object suitable for type-bound function definitions.
    my $table =  Text::Table->new(
	{
	    is_sep => 1,
	    body   => "     "
	},
	{
	    align  => "left"
	},
	{
	    align  => "left"
	},
	{
	    align  => "left"
	},
	{
	    align  => "left"
	},
 	{
	    align  => "left"
	},
        {
	    align  => "left"
	},
	{
	    align  => "left"
	}
       );
    # Iterate over type-bound functions and insert them into the table.
    foreach ( @typeBoundFunctions ) {
	# Determine pass status.
	my $pass = "";
	$pass = ", ".$_->{'pass'}
	    if ( exists($_->{'pass'}) );
	# Determine the connector to use.
	my $connector = "";
	$connector = " => "
	    if ( exists($_->{'name'}) );
	# Determine the name to use.
	my $name = "";
	$name = $_->{'name'}
	    if ( exists($_->{'name'}) );
	# Add a row to the table.
	if ( defined(reftype($_->{'function'})) && reftype($_->{'function'}) eq "ARRAY" ) {
	    # Multiple functions specified. List them, one per row.
	    my $i = 0;
	    foreach my $function ( @{$_->{'function'}} ) {
		++$i;
		# Determine a suitable suffix for this line.
		my $suffix = ", &";
		$suffix = ""
		    if ( $i == scalar(@{$_->{'function'}}) );
		# Add the line to the table.
		if ( $i == 1 ) {
		    $table->add($_->{'type'},$pass," :: ",$name,$connector,$function,$suffix);
		} else {
		    $table->add("     &"    ,""   ,""    ,""   ,""        ,$function,$suffix);
		}
	    }
	} else {
	    # Single function specified. Simply add to the table.
	    $table->add($_->{'type'},$pass," :: ",$name,$connector,$_->{'function'},"");
	}
    }
    # Return the table.
    return $table;
}

sub Insert_Type_Definitions {
    # Generate and insert code for all type definitions.
    my $buildData = shift;
    # Iterate over types.
    foreach ( @{$buildData->{'typesOrder'}} ) {
	# Get the type.
	my $type = $buildData->{'types'}->{$_};
	# Insert the type opening.
	$buildData->{'content'} .= "  type";
	$buildData->{'content'} .= ", public"
	    if ( exists($type->{'isPublic'}) && $type->{'isPublic'} eq "true" );
	$buildData->{'content'} .= ", extends(".$type->{'extends'}.")"
	    if ( exists($type->{'extends'}) );
	$buildData->{'content'} .= " :: ".$type->{'name'}."\n";
	# Insert any comment.
	$buildData->{'content'} .= "  !% ".$type->{'comment'}."\n"
	    if ( exists($type->{'comment'}) );
	# Declare contents private.
	$buildData->{'content'} .= "    private\n";
	# Process any data content.
	$buildData->{'content'} .= &Format_Variable_Defintions($type->{'dataContent'})
	    if ( exists($type->{'dataContent'}) );
	# Generate and insert a type-bound function table.
	if ( exists($type->{'boundFunctions'}) ) {
	    $buildData->{'content'} .= "   contains\n";
	    my $boundFunctionTable = &Bound_Function_Table($type->{'boundFunctions'});   
	    $buildData->{'content'} .= $boundFunctionTable->table();
	}
	# Insert the type closing.
	$buildData->{'content'} .= "  end type ".$type->{'name'}."\n\n";
    }
}

sub Insert_Contains {
    # Insert the "contains" line.
    my $buildData = shift;
    $buildData->{'content'} .= "contains\n\n";
}

sub Format_Variable_Defintions {
    # Generate formatted variable definitions.
    my $variables = shift;
    my %options;
    if ( $#_ >= 1 ) {(%options) = @_};
    # Scan data content searching for repeated attributes.
    my %attributes;
    foreach my $datum ( @{$variables} ) {
	if ( exists($datum->{'attributes'}) ) {
	    foreach ( @{$datum->{'attributes'}} ) {
		(my $attributeName = $_) =~ s/^([^\(]+).*/$1/;
		++$attributes{$attributeName}->{'count'};
		$attributes{$attributeName}->{'column'} = -1;
	    }
	}
    }
    # Find column for aligned attributes.
    my $columnCountMaximum = -1;
    foreach my $datum ( @{$variables} ) {
	if ( exists($datum->{'attributes'}) ) {
	    my $columnCount  = -1;
	    foreach ( sort(@{$datum->{'attributes'}}) ) {
		(my $attributeName = $_) =~ s/^([^\(]+).*/$1/;
		++$columnCount;
		if ( $attributes{$attributeName}->{'count'} > 1 ) {
		    if ( $columnCount > $attributes{$attributeName}->{'column'} ) {
			foreach my $otherAttribute ( sort(keys(%attributes)) ) {
			    ++$attributes{$otherAttribute}->{'column'}
			       if (
				   $attributes{$otherAttribute}->{'column'} >= $columnCount &&
				   $attributes{$otherAttribute}->{'count' } > 1             &&
				   $otherAttribute ne $attributeName
				   );
			}
			$attributes{$attributeName}->{'column'} = $columnCount;
		    }
		    $columnCount = $attributes{$attributeName}->{'column'};
		}
	    }
	    $columnCountMaximum = $columnCount+1
		if ( $columnCount+1 > $columnCountMaximum );
	}
    }
    foreach ( keys(%attributes) ) {
	$columnCountMaximum = $attributes{$_}->{'column'}
	    if ( $attributes{$_}->{'column'} > $columnCountMaximum);
    }
    ++$columnCountMaximum;
    my @attributeColumns;
    push @attributeColumns, {is_sep => 1, body => ""},{align => "left"} foreach (1..$columnCountMaximum);
    # Construct indentation.
    my $indent = "    ";
    $indent = " " x $options{'indent'}
        if ( exists($options{'indent'}) );
    # Create a table for data content.
    my $dataTable = Text::Table->new
	(
	 {
	     is_sep => 1,
	     body   => $indent
	 },
	 {
	     align  => "left"
	 },
	 {
	     is_sep => 1,
	     body   => ""
	 },
	 {
	     align  => "left"
	 },
	 {
	     is_sep => 1,
	     body   => ""
	 },
	 {
	     align  => "left"
	 },
	 {
	     is_sep => 1,
	     body   => ""
	 },
	 {
	     align  => "left"
	 },
	 @attributeColumns,
	 {
	     is_sep => 1,
	     body   => " :: "
	 },
	 {
	     align  => "left"
	 },
	 {
	     is_sep => 1,
	     body   => ""
	 },
	 {
	     align  => "left"
	 }
	);
    # Iterate over all data content.
    foreach ( @{$variables} ) {
	# Construct the type definition.
	my @typeDefinition = ( "", "", "" );
	@typeDefinition = ( "(", $_->{'type'}, ")" )
	    if ( exists($_->{'type'}) );
	# Add attributes.
	my @attributeList;
	if ( exists($_->{'attributes'}) ) {
	    foreach ( sort(@{$_->{'attributes'}}) ) {
		(my $attributeName = $_) =~ s/^([^\(]+).*/$1/;
		if ( $attributes{$attributeName}->{'column'} >= 0 ) {
		    push(@attributeList,"")
			while ( scalar(@attributeList) < $attributes{$attributeName}->{'column'} );
		}
		push(@attributeList,", ".$_);
	    }
	    push(@attributeList,"")
		while ( scalar(@attributeList) < $columnCountMaximum);
	} else {
	    @attributeList = ("") x $columnCountMaximum;
	}
	# Construct any comment.
	my $comment = "";
	$comment = " ! ".$_->{'comment'}
	if ( exists($_->{'comment'}) );
	# Add a row to the table.
	$dataTable->add(
	    $_->{'intrinsic'},
	    @typeDefinition,
	    @attributeList,
	    join(",",@{$_->{'variables'}}),
	    $comment
	    );
    }
    # Return the data content table.
    return $dataTable->table();
}

1;
