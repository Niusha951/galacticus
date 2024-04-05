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
An implementation of the hot halo temperature class which uses a fixed temperature.
!!}

  use :: Dark_Matter_Halo_Scales, only : darkMatterHaloScaleClass

  !![
  <hotHaloTemperatureProfile name="hotHaloTemperatureProfileFixed">
   <description>
    A hot halo temperature profile class which assumes a temperature equal to a fixed amount given in the parameter file.
   </description>
  </hotHaloTemperatureProfile>
  !!]
  type, extends(hotHaloTemperatureProfileClass) :: hotHaloTemperatureProfileFixed
     !!{
     An implementation of the hot halo temperature profile class which uses a fixed temperature.
     !!}
     private
     double precision                         :: Temp 
     class(darkMatterHaloScaleClass), pointer :: darkMatterHaloScale_ => null()
   contains
     final     ::                        fixedDestructor
     procedure :: temperature         => fixedTemperature
     procedure :: temperatureLogSlope => fixedTemperatureLogSlope
  end type hotHaloTemperatureProfileFixed

  interface hotHaloTemperatureProfileFixed
     !!{
     Constructors for the fixed hot halo temperature profile class.
     !!}
     module procedure fixedConstructorParameters
     module procedure fixedConstructorInternal
  end interface hotHaloTemperatureProfileFixed

contains

  function fixedConstructorParameters(parameters) result(self)
    !!{
    Constructor for the fixed class which builds the object from a parameter set.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type (hotHaloTemperatureProfileFixed )                :: self
    type (inputParameters                ), intent(inout) :: parameters
    class(darkMatterHaloScaleClass       ), pointer       :: darkMatterHaloScale_
    double precision                                      :: Temp 

    !![
    <inputParameter>
      <name>Temp</name>
      <defaultValue>1.0d3</defaultValue>
      <description>The amount of fixed temperature given in the parameter file.</description>
      <source>parameters</source>
    </inputParameter>
    <objectBuilder class="darkMatterHaloScale" name="darkMatterHaloScale_" source="parameters"/>
    !!]
    self=hotHaloTemperatureProfileFixed(Temp,darkMatterHaloScale_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="darkMatterHaloScale_"/>
    !!]
    return
  end function fixedConstructorParameters

  function fixedConstructorInternal(Temp,darkMatterHaloScale_) result(self)
    !!{
    Internal constructor for the fixed class.
    !!}
    implicit none
    type (hotHaloTemperatureProfileFixed )                        :: self
    double precision                      , intent(in   )         :: Temp 
    class(darkMatterHaloScaleClass       ), intent(in   ), target :: darkMatterHaloScale_
    !![
    <constructorAssign variables="Temp, *darkMatterHaloScale_"/>
    !!]

    return
  end function fixedConstructorInternal

  subroutine fixedDestructor(self)
    !!{
    Destructor for the {\normalfont \ttfamily fixed} hot halo temperature profile class.
    !!}
    implicit none
    type(hotHaloTemperatureProfileFixed), intent(inout) :: self

    !![
    <objectDestructor name="self%darkMatterHaloScale_"/>
    !!]
    return
  end subroutine fixedDestructor

  double precision function fixedTemperature(self,node,radius)
    !!{
    Return the density in a {\normalfont \ttfamily fixed} hot halo mass distribution.
    !!}
    implicit none
    class           (hotHaloTemperatureProfileFixed ), intent(inout) :: self
    type            (treeNode                       ), intent(inout) :: node
    double precision                                 , intent(in   ) :: radius
    !$GLC attributes unused :: radius

    fixedTemperature=self%Temp
    return
  end function fixedTemperature

  double precision function fixedTemperatureLogSlope(self,node,radius)
    !!{
    Return the logarithmic slope of the density profile in a {\normalfont \ttfamily fixed} hot halo mass
    distribution.
    !!}
    implicit none
    class           (hotHaloTemperatureProfileFixed), intent(inout) :: self
    type            (treeNode                       ), intent(inout) :: node
    double precision                                 , intent(in   ) :: radius
    !$GLC attributes unused :: self, node, radius

    fixedTemperatureLogSlope=0.0d0
    return
  end function fixedTemperatureLogSlope

