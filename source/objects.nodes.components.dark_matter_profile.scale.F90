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
Contains a module which implements a dark matter profile method that provides a scale radius.
!!}

module Node_Component_Dark_Matter_Profile_Scale
  !!{
  Implements a dark matter profile method that provides a scale radius.
  !!}
  use :: Dark_Matter_Halo_Scales, only : darkMatterHaloScaleClass
  implicit none
  private
  public :: Node_Component_Dark_Matter_Profile_Scale_Scale_Set        , Node_Component_Dark_Matter_Profile_Scale_Plausibility       , &
       &    Node_Component_Dark_Matter_Profile_Scale_Thread_Initialize, Node_Component_Dark_Matter_Profile_Scale_Thread_Uninitialize, &
       &    Node_Component_Dark_Matter_Profile_Scale_State_Store      , Node_Component_Dark_Matter_Profile_Scale_State_Restore

  !![
  <component>
   <class>darkMatterProfile</class>
   <name>scale</name>
   <isDefault>true</isDefault>
   <properties>
    <property>
      <name>scale</name>
      <type>double</type>
      <rank>0</rank>
      <attributes isSettable="true" isGettable="true" isEvolvable="true" />
      <output unitsInSI="megaParsec" comment="Scale radius of the dark matter profile [Mpc]."/>
      <classDefault>-1.0d0</classDefault>
    </property>
   </properties>
  </component>
  !!]

  ! Objects used by this component.
  class(darkMatterHaloScaleClass), pointer :: darkMatterHaloScale_
  !$omp threadprivate(darkMatterHaloScale_)
  
contains

  !![
  <nodeComponentThreadInitializationTask>
   <unitName>Node_Component_Dark_Matter_Profile_Scale_Thread_Initialize</unitName>
  </nodeComponentThreadInitializationTask>
  !!]
  subroutine Node_Component_Dark_Matter_Profile_Scale_Thread_Initialize(parameters_)
    !!{
    Initializes the tree node scale dark matter profile module.
    !!}
    use :: Galacticus_Nodes, only : defaultDarkMatterProfileComponent
    use :: Input_Parameters, only : inputParameter                   , inputParameters
    implicit none
    type(inputParameters), intent(inout) :: parameters_
    !$GLC attributes unused :: parameters_

    if (defaultDarkMatterProfileComponent%scaleIsActive()) then
       !![
       <objectBuilder class="darkMatterHaloScale" name="darkMatterHaloScale_" source="parameters_"/>
       !!]
     end if
     return
  end subroutine Node_Component_Dark_Matter_Profile_Scale_Thread_Initialize

  !![
  <nodeComponentThreadUninitializationTask>
   <unitName>Node_Component_Dark_Matter_Profile_Scale_Thread_Uninitialize</unitName>
  </nodeComponentThreadUninitializationTask>
  !!]
  subroutine Node_Component_Dark_Matter_Profile_Scale_Thread_Uninitialize()
    !!{
    Uninitializes the tree node scale dark matter profile module.
    !!}
    use :: Galacticus_Nodes, only : defaultDarkMatterProfileComponent
    implicit none

    if (defaultDarkMatterProfileComponent%scaleIsActive()) then
       !![
       <objectDestructor name="darkMatterHaloScale_"/>
       !!]
    end if
    return
  end subroutine Node_Component_Dark_Matter_Profile_Scale_Thread_Uninitialize

  !![
  <radiusSolverPlausibility>
   <unitName>Node_Component_Dark_Matter_Profile_Scale_Plausibility</unitName>
   <after>Node_Component_Basic_Standard_Plausibility</after>
  </radiusSolverPlausibility>
  !!]
  subroutine Node_Component_Dark_Matter_Profile_Scale_Plausibility(node)
    !!{
    Determines whether the dark matter profile is physically plausible for radius solving tasks.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentDarkMatterProfile, nodeComponentDarkMatterProfileScale, treeNode
    implicit none
    type   (treeNode                      ), intent(inout) :: node
    class  (nodeComponentDarkMatterProfile), pointer       :: darkMatterProfile

    ! Return immediately if already non-plausible.
    if (.not.node%isPhysicallyPlausible.and..not.node%isSolvable) return
    ! Get the dark matter profile component.
    darkMatterProfile => node%darkMatterProfile()
    ! Ensure that it is of the scale class.
    select type (darkMatterProfile)
    class is (nodeComponentDarkMatterProfileScale)
       if     (                                                                      &
            &   darkMatterProfile%scale() <= 0.0d0                                   &
            &  .or.                                                                  &
            &   darkMatterProfile%scale() >  darkMatterHaloScale_%radiusVirial(node) &
            & ) then
          node%isPhysicallyPlausible=.false.
          node%isSolvable           =.false.
       end if
    end select
    return
  end subroutine Node_Component_Dark_Matter_Profile_Scale_Plausibility

  !![
  <scaleSetTask>
   <unitName>Node_Component_Dark_Matter_Profile_Scale_Scale_Set</unitName>
  </scaleSetTask>
  !!]
  subroutine Node_Component_Dark_Matter_Profile_Scale_Scale_Set(node)
    !!{
    Set scales for properties of {\normalfont \ttfamily node}.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentDarkMatterProfile, nodeComponentDarkMatterProfileScale, treeNode
    implicit none
    type            (treeNode                      ), intent(inout), pointer :: node
    class           (nodeComponentDarkMatterProfile)               , pointer :: darkMatterProfile
    double precision                                , parameter              :: radiusScaleMinimum=1.0d-6

    ! Get the dark matter profile component.
    darkMatterProfile => node%darkMatterProfile()
    ! Ensure it is of the scale class.
    select type (darkMatterProfile)
    class is (nodeComponentDarkMatterProfileScale)
       ! Set scale for the scale radius.
       call darkMatterProfile%scaleScale(max(darkMatterProfile%scale(),radiusScaleMinimum))
    end select
    return
  end subroutine Node_Component_Dark_Matter_Profile_Scale_Scale_Set

  !![
  <stateStoreTask>
   <unitName>Node_Component_Dark_Matter_Profile_Scale_State_Store</unitName>
  </stateStoreTask>
  !!]
  subroutine Node_Component_Dark_Matter_Profile_Scale_State_Store(stateFile,gslStateFile,stateOperationID)
    !!{
    Store object state,
    !!}
    use            :: Display      , only : displayMessage, verbosityLevelInfo
    use, intrinsic :: ISO_C_Binding, only : c_ptr         , c_size_t
    implicit none
    integer          , intent(in   ) :: stateFile
    integer(c_size_t), intent(in   ) :: stateOperationID
    type   (c_ptr   ), intent(in   ) :: gslStateFile

    call displayMessage('Storing state for: componentDarkMatterProfile -> scale',verbosity=verbosityLevelInfo)
    !![
    <stateStore variables="darkMatterHaloScale_"/>
    !!]
    return
  end subroutine Node_Component_Dark_Matter_Profile_Scale_State_Store

  !![
  <stateRetrieveTask>
   <unitName>Node_Component_Dark_Matter_Profile_Scale_State_Restore</unitName>
  </stateRetrieveTask>
  !!]
  subroutine Node_Component_Dark_Matter_Profile_Scale_State_Restore(stateFile,gslStateFile,stateOperationID)
    !!{
    Retrieve object state.
    !!}
    use            :: Display      , only : displayMessage, verbosityLevelInfo
    use, intrinsic :: ISO_C_Binding, only : c_ptr         , c_size_t
    implicit none
    integer          , intent(in   ) :: stateFile
    integer(c_size_t), intent(in   ) :: stateOperationID
    type   (c_ptr   ), intent(in   ) :: gslStateFile

    call displayMessage('Retrieving state for: componentDarkMatterProfile -> scale',verbosity=verbosityLevelInfo)
    !![
    <stateRestore variables="darkMatterHaloScale_"/>
    !!]
    return
  end subroutine Node_Component_Dark_Matter_Profile_Scale_State_Restore

end module Node_Component_Dark_Matter_Profile_Scale
