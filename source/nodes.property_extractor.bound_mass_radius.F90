!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023, 2024
!!    Andrew Benson <abenson@carnegiescience.edu>
!!
!! This file is part of Galacticus.
!!
!!    Galacticus is free software: you can redistribute it and/or modify
!!    it under the terms of the GNU General Public License as published by
!!    the Free Software Foundation, either version 3 of the License, or
!!    (at your option) any later version.
!!
!!    Galacticus is distributed in the hope that it will be useful,
!!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!!    GNU General Public License for more details.
!!
!!    You should have received a copy of the GNU General Public License
!!    along with Galacticus.  If not, see <http://www.gnu.org/licenses/>.

!!{
Contains a module which implements a property extractor class that extracts the radius enclosing the current bound mass.
!!}

  use :: Galactic_Structure, only : galacticStructureClass

  !![
  <nodePropertyExtractor name="nodePropertyExtractorRadiusBoundMass">
   <description>
    A property extractor class that extracts the radius enclosing the current bound mass.
   </description>
  </nodePropertyExtractor>
  !!]
  type, extends(nodePropertyExtractorScalar) :: nodePropertyExtractorRadiusBoundMass
     !!{
     A property extractor class that extracts the radius enclosing the current bound mass.
     !!}
     private
     class(galacticStructureClass), pointer :: galacticStructure_ => null()
   contains
     final     ::                radiusBoundMassDestructor
     procedure :: extract     => radiusBoundMassExtract
     procedure :: name        => radiusBoundMassName
     procedure :: description => radiusBoundMassDescription
     procedure :: unitsInSI   => radiusBoundMassUnitsInSI
  end type nodePropertyExtractorRadiusBoundMass

  interface nodePropertyExtractorRadiusBoundMass
     !!{
     Constructors for the ``radiusBoundMass'' output analysis class.
     !!}
     module procedure radiusBoundMassConstructorParameters
     module procedure radiusBoundMassConstructorInternal
  end interface nodePropertyExtractorRadiusBoundMass

contains

  function radiusBoundMassConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily radiusBoundMass} property extractor class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type (nodePropertyExtractorRadiusBoundMass)                :: self
    type (inputParameters                     ), intent(inout) :: parameters
    class(galacticStructureClass              ), pointer       :: galacticStructure_
    
    !![
    <objectBuilder class="galacticStructure" name="galacticStructure_" source="parameters"/>
    !!]
    self=nodePropertyExtractorRadiusBoundMass(galacticStructure_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="galacticStructure_"/>
    !!]
    return
  end function radiusBoundMassConstructorParameters

  function radiusBoundMassConstructorInternal(galacticStructure_) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily radiusBoundMass} property extractor class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type (nodePropertyExtractorRadiusBoundMass)                        :: self
    class(galacticStructureClass              ), intent(in   ), target :: galacticStructure_
    !![
    <constructorAssign variables="*galacticStructure_"/>
    !!]

    return
  end function radiusBoundMassConstructorInternal
  
  subroutine radiusBoundMassDestructor(self)
    !!{
    Destructor for the {\normalfont \ttfamily radiusBoundMass} property extractor class which takes a parameter set as input.
    !!}
    implicit none
    type(nodePropertyExtractorRadiusBoundMass), intent(inout) :: self
    
    !![
    <objectDestructor name="self%galacticStructure_"/>
    !!]
    return
  end subroutine radiusBoundMassDestructor

  double precision function radiusBoundMassExtract(self,node,instance)
    !!{
    Implement a bound mass radius property extractor.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentSatellite
    implicit none
    class           (nodePropertyExtractorRadiusBoundMass), intent(inout), target   :: self
    type            (treeNode                            ), intent(inout), target   :: node
    type            (multiCounter                        ), intent(inout), optional :: instance
    class           (nodeComponentSatellite              )               , pointer  :: satellite
    double precision                                                                :: massTotal
    !$GLC attributes unused :: instance

    satellite              => node%satellite()
    massTotal              =  self%galacticStructure_%massEnclosed       (node                                          )
    radiusBoundMassExtract =  self%galacticStructure_%radiusEnclosingMass(node,mass=min(satellite%boundMass(),massTotal))
    return
  end function radiusBoundMassExtract

  function radiusBoundMassName(self)
    !!{
    Return the name of the bound mass radius property.
    !!}
    implicit none
    type (varying_string                      )                :: radiusBoundMassName
    class(nodePropertyExtractorRadiusBoundMass), intent(inout) :: self
    !$GLC attributes unused :: self

    radiusBoundMassName=var_str('satelliteRadiusBoundMass')
    return
  end function radiusBoundMassName

  function radiusBoundMassDescription(self)
    !!{
    Return a description of the bound mass radius property.
    !!}
    implicit none
    type (varying_string                      )                :: radiusBoundMassDescription
    class(nodePropertyExtractorRadiusBoundMass), intent(inout) :: self
    !$GLC attributes unused :: self

    radiusBoundMassDescription=var_str('Radius enclosing the bound mass of the halo [Mpc].')
    return
  end function radiusBoundMassDescription

  double precision function radiusBoundMassUnitsInSI(self)
    !!{
    Return the units of the bound mass radius property in the SI system.
    !!}
    use :: Numerical_Constants_Astronomical, only : megaParsec
    implicit none
    class(nodePropertyExtractorRadiusBoundMass), intent(inout) :: self
    !$GLC attributes unused :: self

    radiusBoundMassUnitsInSI=megaParsec
    return
  end function radiusBoundMassUnitsInSI


