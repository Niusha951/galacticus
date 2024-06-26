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
  Implements a model of the ram pressure stripping force from hot halos based on orbital position within the host halo.
  !!}

  use :: Hot_Halo_Mass_Distributions, only : hotHaloMassDistributionClass
  use :: Galactic_Structure         , only : galacticStructureClass

  !![
  <hotHaloRamPressureForce name="hotHaloRamPressureForceRelativePosition">
   <description>
    A hot halo ram pressure force class which computes the force based on the current position relative to the host
    halo. Specifically, the ram pressure force is    
    \begin{equation}
    \mathcal{F}_\mathrm{ram, hot, host} = \rho_\mathrm{hot, host}(r) v^2(r),
    \end{equation}
    where $\rho_\mathrm{hot, host}(r)$ is the hot halo density profile of the nodes host halo, $v(r)$ is the orbital velocity
    of the node in that host, and $r$ is the instantaneous distance to the host halo.
   </description>
  </hotHaloRamPressureForce>
  !!]
  type, extends(hotHaloRamPressureForceClass) :: hotHaloRamPressureForceRelativePosition
     !!{
     Implementation of a hot halo ram pressure force class based on orbital position within the host halo.
     !!}
     private
     class(hotHaloMassDistributionClass), pointer :: hotHaloMassDistribution_ => null()
     class(galacticStructureClass      ), pointer :: galacticStructure_       => null()
   contains
     final     ::          relativePositionDestructor
     procedure :: force => relativePositionForce
  end type hotHaloRamPressureForceRelativePosition

  interface hotHaloRamPressureForceRelativePosition
     !!{
     Constructors for the {\normalfont \ttfamily relativePosition} hot halo ram pressure force class.
     !!}
     module procedure relativePositionConstructorParameters
     module procedure relativePositionConstructorInternal
  end interface hotHaloRamPressureForceRelativePosition

contains

  function relativePositionConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily relativePosition} hot halo ram pressure force class which builds the object from a parameter set.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type (hotHaloRamPressureForceRelativePosition)                :: self
    type (inputParameters                        ), intent(inout) :: parameters
    class(hotHaloMassDistributionClass           ), pointer       :: hotHaloMassDistribution_
    class(galacticStructureClass                 ), pointer       :: galacticStructure_

    !![
    <objectBuilder class="hotHaloMassDistribution" name="hotHaloMassDistribution_" source="parameters"/>
    <objectBuilder class="galacticStructure"       name="galacticStructure_"       source="parameters"/>
    !!]
    self=hotHaloRamPressureForceRelativePosition(hotHaloMassDistribution_,galacticStructure_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="hotHaloMassDistribution_"/>
    <objectDestructor name="galacticStructure_"      />
    !!]
    return
  end function relativePositionConstructorParameters

  function relativePositionConstructorInternal(hotHaloMassDistribution_,galacticStructure_) result(self)
    !!{
    Internal constructor for the {\normalfont \ttfamily relativePosition} hot halo ram pressure force class.
    !!}
    use :: Array_Utilities , only : operator(.intersection.)
    use :: Error           , only : Error_Report             , Component_List
    use :: Galacticus_Nodes, only : defaultPositionComponent
    implicit none
    type (hotHaloRamPressureForceRelativePosition)                        :: self
    class(hotHaloMassDistributionClass           ), intent(in   ), target :: hotHaloMassDistribution_
    class(galacticStructureClass                 ), intent(in   ), target :: galacticStructure_
    !![
    <constructorAssign variables="*hotHaloMassDistribution_, *galacticStructure_"/>
    !!]

    ! Ensure that required methods are supported.
    if     (                                                                                                                          &
         &  .not.                                                                                                                     &
         &       (                                                                                                                    &
         &        defaultPositionComponent%positionIsGettable().and.                                                                  &
         &        defaultPositionComponent%velocityIsGettable()                                                                       &
         &  )                                                                                                                         &
         & ) call Error_Report                                                                                                        &
         &        (                                                                                                                   &
         &         'this method requires that position, and velocity properties must all be gettable for the `position` component.'// &
         &         Component_List(                                                                                                    &
         &                        'position'                                                                                       ,  &
         &                        defaultPositionComponent%positionAttributeMatch(requireGettable=.true.).intersection.               &
         &                        defaultPositionComponent%velocityAttributeMatch(requireGettable=.true.)                             &
         &                       )                                                                                                 // &
         &         {introspection:location}                                                                                           &
         &        )
    
    return
  end function relativePositionConstructorInternal

  subroutine relativePositionDestructor(self)
    !!{
    Destructor for the {\normalfont \ttfamily relativePosition} hot halo ram pressure force class.
    !!}
    implicit none
    type(hotHaloRamPressureForceRelativePosition), intent(inout) :: self

    !![
    <objectDestructor name="self%hotHaloMassDistribution_"/>
    <objectDestructor name="self%galacticStructure_"      />
    !!]
    return
  end subroutine relativePositionDestructor

  double precision function relativePositionForce(self,node) result(force)
    !!{
    Return a ram pressure force due to the hot halo based on orbital position within the host halo.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentPosition, nodeComponentBasic
    use :: Vectors         , only : Vector_Magnitude
    implicit none
    class           (hotHaloRamPressureForceRelativePosition), intent(inout) :: self
    type            (treeNode                               ), intent(inout) :: node
    class           (nodeComponentPosition                  ), pointer       :: position       , positionHost
    class           (nodeComponentBasic                     ), pointer       :: basic          , basicPrevious   , &
         &                                                                      basicCurrent   , basicHost
    type            (treeNode                               ), pointer       :: nodeHost       , nodeHostPrevious, &
         &                                                                      nodeHostCurrent
    double precision                                                         :: radiusRelative , velocityRelative

    ! Find the host node. Seek the descendant of the node closest in time to our satellite node. This is necessary as satellites
    ! can evolve ahead of their hosts.
    basic            => node%basic ()
    nodeHostPrevious => node%parent
    nodeHostCurrent  => node%parent
    basicPrevious => nodeHostPrevious%basic()
    do while (associated(nodeHostCurrent))
       basicCurrent => nodeHostCurrent%basic()
       if (basicCurrent%time() > basic%time()) exit
       basicPrevious    => basicCurrent
       nodeHostPrevious => nodeHostCurrent
       nodeHostCurrent  => nodeHostCurrent%parent
    end do
    if     (                                        &
         &   abs(basicPrevious%time()-basic%time()) &
         &  <                                       &
         &   abs(basicCurrent %time()-basic%time()) &
         & ) then
       nodeHost  => nodeHostPrevious
       basicHost => basicPrevious
    else
       nodeHost  => nodeHostCurrent
       basicHost => basicCurrent
    end if
    ! Get the position components.
    position         =>  node    %position()
    positionHost     =>  nodeHost%position()
    ! Compute orbital position and velocity.
    radiusRelative   =  +Vector_Magnitude(position%position()-positionHost%position())
    velocityRelative =  +Vector_Magnitude(position%velocity()-positionHost%velocity())
    ! Find the ram pressure force at this orbital radius.
    force                =  +self%hotHaloMassDistribution_%density        (nodeHost,radiusRelative)    &
         &                  *                              velocityRelative                        **2
    return
  end function relativePositionForce
