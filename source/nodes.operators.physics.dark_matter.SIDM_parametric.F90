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
  Implements a node operator class that maps the CDM solution to SIDM based on the parametric model introduced in Yang et al. 2023: arXiv:2305.16176
  !!}
  
  use :: Dark_Matter_Halo_Mass_Accretion_Histories, only : darkMatterHaloMassAccretionHistory  , darkMatterHaloMassAccretionHistoryClass

  !![
  <nodeOperator name="nodeOperatorSIDMParametric">
   <description>
     A node operator class that maps the CDM solution to SIDM based on the parametric model introduced in Yang et al. 2023: arXiv:2305.16176
   </description>
  </nodeOperator>
  !!]
  type, extends(nodeOperatorClass) :: nodeOperatorSIDMParametric
     !!{
     A node operator class that maps the CDM solution to SIDM based on the parametric model introduced in Yang et al. 2023: arXiv:2305.16176
       !!}
     private
     class(darkMatterParticleClass                ), pointer :: darkMatterParticle_        => null()
     class(darkMatterHaloMassAccretionHistoryClass), pointer :: darkMatterHaloMassAccretionHistory_ => null()
     class(darkMatterProfileDMOClass              ), pointer :: darkMatterProfileDMO_
 
     integer                                                 :: tauID, VmaxSIDMID, RmaxSIDMID, nodeFormationTimeSIDMID
   contains
     !![
     <methods>
       <method method="calculateTau" description="Computes tau, VmaxSIDM and RmaxSIDM for each node."/>
     </methods>
     !!]
     final     ::                                        SIDMParametricDestructor
     procedure :: nodeTreeInitialize                  => SIDMParametricNodeTreeInitialize
     procedure :: nodePromote                         => SIDMParametricNodePromote
     procedure :: differentialEvolutionScales         => SIDMParametriCalculateTauDifferentialEvolutionScale
     procedure :: differentialEvolution               => SIDMParametriCalculateTauDifferentialEvolution
!     procedure :: differentialEvolutionAnalytics      => SIDMParametriDifferentialVmaxAnalytics
!     procedure :: differentialEvolutionSolveAnalytics => SIDMParametriDifferentialVmaxSolveAnalytics
     procedure :: getTauID                            => getTauID_SIDMParametric
     procedure :: getVmaxSIDMID                       => getVmaxSIDMID_SIDMParametric
     procedure :: getRmaxSIDMID                       => getRmaxSIDMID_SIDMParametric
  end type nodeOperatorSIDMParametric
  
  interface nodeOperatorSIDMParametric
     !!{
     Constructors for the {\normalfont \ttfamily SIDMParametric} node operator class.
     !!}
     module procedure SIDMParametricConstructorParameters
     module procedure SIDMParametricConstructorInternal
  end interface nodeOperatorSIDMParametric
  
contains
  
  function SIDMParametricConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily SIDMParametric} node operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type (nodeOperatorSIDMParametric)                  :: self
    type (inputParameters             ), intent(inout) :: parameters
    class(darkMatterParticleClass                ), pointer :: darkMatterParticle_
    class(darkMatterHaloMassAccretionHistoryClass), pointer :: darkMatterHaloMassAccretionHistory_
    class(darkMatterProfileDMOClass              ), pointer :: darkMatterProfileDMO_

    !![
    <objectBuilder class="darkMatterParticle" name="darkMatterParticle_" source="parameters"/>
    <objectBuilder class="darkMatterHaloMassAccretionHistory" name="darkMatterHaloMassAccretionHistory_" source="parameters"/>
    <objectBuilder class="darkMatterProfileDMO" name="darkMatterProfileDMO_" source="parameters"/>
    !!]
    self=nodeOperatorSIDMParametric(darkMatterParticle_,darkMatterHaloMassAccretionHistory_, darkMatterProfileDMO_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="darkMatterParticle_"/>
    <objectDestructor name="darkMatterHaloMassAccretionHistory_"/>
    <objectDestructor name="darkMatterProfileDMO_"/>
    !!]
    return
  end function SIDMParametricConstructorParameters

  function SIDMParametricConstructorInternal(darkMatterParticle_, darkMatterHaloMassAccretionHistory_, darkMatterProfileDMO_) result(self)
    !!{
    Internal constructor for the {\normalfont \ttfamily SIDMParametric} node operator class.
    !!}
    implicit none
    type (nodeOperatorSIDMParametric)                        :: self
    class(darkMatterParticleClass   ), intent(in   ), target :: darkMatterParticle_
    class(darkMatterHaloMassAccretionHistoryClass), intent(in   ),target  :: darkMatterHaloMassAccretionHistory_
    class(darkMatterProfileDMOClass), intent(in   ), target :: darkMatterProfileDMO_

    !![
    <constructorAssign variables="*darkMatterParticle_, *darkMatterHaloMassAccretionHistory_, *darkMatterProfileDMO_"/>
    !!]
    
    !![
    <addMetaProperty component="darkMatterProfile" name="tau"      id="self%tauID"      isEvolvable="yes" isCreator="yes"/>
    <addMetaProperty component="darkMatterProfile" name="VmaxSIDM" id="self%VmaxSIDMID" isEvolvable="yes"  isCreator="yes"/>
    <addMetaProperty component="darkMatterProfile" name="RmaxSIDM" id="self%RmaxSIDMID" isEvolvable="yes"  isCreator="yes"/>
    <addMetaProperty component="basic" name="nodeFormationTimeSIDM" id="self%nodeFormationTimeSIDMID" isEvolvable="no" isCreator="yes"/>
    !!]
    return
  end function SIDMParametricConstructorInternal

  subroutine SIDMParametricDestructor(self)
    !!{
    Destructor for the {\normalfont \ttfamily SIDMParametric} node operator class.
    !!}
    implicit none
    type(nodeOperatorSIDMParametric), intent(inout) :: self

    !![
    <objectDestructor name="self%darkMatterParticle_"/>
    <objectDestructor name="self%darkMatterHaloMassAccretionHistory_"/>
    <objectDestructor name="self%darkMatterProfileDMO_"/>
    !!]
    return
  end subroutine SIDMParametricDestructor

  subroutine SIDMParametricNodeTreeInitialize(self,node)
    !!{
    Initialize the SIDMParametric of all nodes in the tree.    
    !!}

    use :: Galacticus_Nodes  , only : nodeComponentBasic
    use :: Dark_Matter_Halo_Formation_Times, only : Dark_Matter_Halo_Formation_Time

    implicit none
    class(nodeOperatorSIDMParametric), intent(inout), target  :: self
    type (treeNode                  ), intent(inout), target  :: node
    type (treeNode                  ),                pointer :: nodeParent, nodeFormation
    class    (nodeComponentBasic        ), pointer       :: basic, basicParent, basicNode
    
    double precision :: formationMassFraction = 0.5d0
    double precision :: timeFormation

    print *, 'Test inside SIDM parametric NodeTreeInitialize ...'
    call self%nodeInitialize(node)

!    basicNode => node%basic()
!    print *, 'node time: ', basicNode%time()

    !calculating the formation time
    nodeParent    => node
    basic => node%basic()
    do while (associated(nodeParent))
       basicParent => nodeParent%basic()
       print *, 'what is parent time: ', basicParent%time()
       print *, 'parent index: ', nodeParent%index()
       timeFormation =  Dark_Matter_Halo_Formation_Time(node=nodeParent, formationMassFraction=formationMassFraction, darkMatterHaloMassAccretionHistory_=self%darkMatterHaloMassAccretionHistory_)
       if (nodeParent%isPrimaryProgenitor()) then
          nodeParent => nodeParent%parent
       else
          nodeParent => null()
       end if
    end do
    call basic%floatRank0MetaPropertySet(self%nodeFormationTimeSIDMID,timeFormation)
    print *, 'Formation time: ', basic%floatRank0MetaPropertyGet(self%nodeFormationTimeSIDMID)
    print *, 'node index: ', node%index()

!    !calculating the formation time
!    if (.not.associated(node%firstChild)) then
!       print *, 'Node that has no first child, has the time: ', basicNode%time()
!       nodeParent    => node
!       nodeFormation => node
!       do while (associated(nodeParent))
!!          basic         => nodeParent%basic()
!          basic => nodeFormation%basic()
!          basicParent => nodeParent%basic()
!          print *, 'if it has no child, now what is its parent time: ', basicParent%time()
!          timeFormation =  Dark_Matter_Halo_Formation_Time(node=nodeParent, formationMassFraction=formationMassFraction, darkMatterHaloMassAccretionHistory_=self%darkMatterHaloMassAccretionHistory_)
!          if (nodeParent%isPrimaryProgenitor()) then
!             nodeParent => nodeParent%parent
!          else
!             nodeParent => null()
!          end if 
!       end do
!       call basic%floatRank0MetaPropertySet(self%nodeFormationTimeID,timeFormation)
!       print *, 'Formation time: ', basic%floatRank0MetaPropertyGet(self%nodeFormationTimeID)
!    end if

    print *, 'Test end of initialization ...'

!    call self%calculateTau(node, basic%floatRank0MetaPropertyGet(self%nodeFormationTimeID))
    return
  end subroutine SIDMParametricNodeTreeInitialize

  subroutine SIDMParametricNodePromote(self,node)
    !!{
    Ensure that {\normalfont \ttfamily node} is ready for promotion to its parent.
    !!}

    use :: Error             , only : Error_Report
    use :: Galacticus_Nodes  , only : nodeComponentBasic, nodeComponentDarkMatterProfile
    use :: ISO_Varying_String, only : var_str           , varying_string, operator(//)
    use :: String_Handling   , only : operator(//)

    implicit none
    class    (nodeOperatorSIDMParametric), intent(inout) :: self
    type     (treeNode                  ), intent(inout) :: node
    type     (treeNode                  ), pointer       :: nodeParent
    class    (nodeComponentBasic        ), pointer       :: basicParent, basic
    class(nodeComponentDarkMatterProfile), pointer       :: darkMatterProfile, darkMatterProfileParent
    type     (varying_string            )                :: message
    character(len=12                    )                :: label
    !$GLC attributes unused :: self
    
    nodeParent  => node      %parent
    basic       => node      %basic ()
    basicParent => nodeParent%basic ()
    darkMatterProfile => node%darkMatterProfile()    
    darkMatterProfileParent => nodeParent%darkMatterProfile()
    ! Ensure the two halos exist at the same time.
    if (basic%time() /= basicParent%time()) then
       message=var_str("node [")//node%index()//"] has not been evolved to its parent ["//nodeParent%index()//"]"//char(10)
       write (label,'(f12.6)') basic%time()
       message=message//"    node is at time: "//label//" Gyr"//char(10)
       write (label,'(f12.6)') basicParent%time()
       message=message//"  parent is at time: "//label//" Gyr"
       call Error_Report(message//{introspection:location})
    end if


    print *, 'Test inside nodePromote ...'

!    call basic%floatRank0MetaPropertySet(self%nodeFormationTimeID,basicParent%floatRank0MetaPropertyGet(self%nodeFormationTimeID))
!    call darkMatterProfile%floatRank0MetaPropertySet(self%tauID,darkMatterProfileParent%floatRank0MetaPropertyGet(self%tauID))
!    call darkMatterProfile%floatRank0MetaPropertySet(self%VmaxSIDMID,darkMatterProfile%floatRank0MetaPropertyGet(self%VmaxSIDMID))    
!    call darkMatterProfile%floatRank0MetaPropertySet(self%RmaxSIDMID,darkMatterProfileParent%floatRank0MetaPropertyGet(self%RmaxSIDMID))
!    call basic%floatRank0MetaPropertySet(self%nodeFormationTimeID, basicParent%floatRank0MetaPropertyGet(self%nodeFormationTimeID))
    return
  end subroutine SIDMParametricNodePromote
 

  subroutine SIDMParametriCalculateTauDifferentialEvolutionScale(self, node)
    use Galacticus_Nodes, only: treeNode, nodeComponentBasic, nodeComponentDarkMatterProfile

    type(treeNode), intent(inout) :: node
    class(nodeOperatorSIDMParametric), intent(inout) :: self
    class(nodeComponentDarkMatterProfile), pointer :: darkMatterProfile

    print *, 'Test inside evolutionScale ...'

    darkMatterProfile => node%darkMatterProfile()
    call darkMatterProfile%floatRank0MetaPropertyScale(self%tauID, 1.0d0)
!    call darkMatterProfile%floatRank0MetaPropertyScale(self%VmaxSIDMID, 1.0d0)
    return
  end subroutine SIDMParametriCalculateTauDifferentialEvolutionScale

  subroutine SIDMParametriCalculateTauDifferentialEvolution(self, node,interrupt,functionInterrupt,propertyType)
    use Galacticus_Nodes, only: treeNode, nodeComponentBasic, nodeComponentDarkMatterProfile
!    use Dark_Matter_Halo_Mass_Accretion_Histories, only: darkMatterHaloMassAccretionHistoryClass
!    use Dark_Matter_Profiles_DMO , only : darkMatterProfileDMOClass
    use :: Dark_Matter_Halo_Formation_Times, only : Dark_Matter_Halo_Formation_Time

    type(treeNode), intent(inout), target :: node
    type(treeNode), pointer :: nodeFinal
    type(treeNode), pointer :: nodeWork
    type(treeNode), pointer :: nodeStart
    type(treeNode), pointer :: nodeParent
    !    type(darkMatterProfileDMO), pointer :: darkMatterProfileDMO_

    class(nodeOperatorSIDMParametric), intent(inout), target :: self
!    class(darkMatterHaloMassAccretionHistoryClass), target :: darkMatterHaloMassAccretionHistory_
!    class(darkMatterHaloMassAccretionHistoryClass), pointer :: darkMatterHaloMassAccretionHistory_
    class(nodeComponentBasic), pointer :: basic, basicParent
    class(nodeComponentDarkMatterProfile), pointer :: darkMatterProfile
!    class(darkMatterProfileDMOClass), pointer :: darkMatterProfileDMO_

    double precision :: timeFormation
    double precision :: formationMassFraction = 0.5d0
    double precision :: tau, time, timePrevious, tc
    double precision :: RmaxNFW0, VmaxNFW0, r_sNFW0, rho_sNFW0, rho_s, r_s, r_c
    double precision :: VmaxSIDM=0.0d0, RmaxSIDM=0.0d0, VmaxSIDMPrevious, RmaxSIDMPrevious, VmaxCDM, RmaxCDM, VmaxCDMPrevious, RmaxCDMPrevious
    double precision :: dtr
    double precision :: tau, VmaxSIDM, tc, dvdt

    logical, intent(inout) :: interrupt
    procedure(interruptTask), intent(inout), pointer :: functionInterrupt
    integer, intent(in   ) :: propertyType

    !$GLC attributes unused :: interrupt, functionInterrupt, propertyType


    print *, 'Test inside DiffEvolution ...'

    print *, 'node index: ', node%index()
    !nodeParent => node%parent
    print *, 'parent index: ', node%Parent%index()

    basic => node%basic()
    time = basic%time()
!!    timeFormation = basic%floatRank0MetaPropertyGet(self%nodeFormationTimeSIDMID)
    timeFormation = 2.1903d0
!    timeFormation = 1.7887838371779103d0

!    !calculating the formation time
!    nodeParent => node
!!    basic => node%basic()
!    do while (associated(nodeParent))
!       basicParent => nodeParent%basic()
!       print *, 'what is parent time: ', basicParent%time()
!       timeFormation =  Dark_Matter_Halo_Formation_Time(node=nodeParent, formationMassFraction=formationMassFraction, darkMatterHaloMassAccretionHistory_=self%darkMatterHaloMassAccretionHistory_)
!       if (nodeParent%isPrimaryProgenitor()) then
!          nodeParent => nodeParent%parent
!       else
!          nodeParent => null()
!       end if
!    end do
!    print *, 'Formation time: ', timeFormation

    darkMatterProfile => node%darkMatterProfile()

    print *, 'Time: ', time, timeFormation
!    print *, 'nodeFormationTimeID: ', self%nodeFormationTimeID

!    timeFormation = 2.1903d0
    ! Open the output file
!    open(unit=20, file='test_output_data.txt', status='replace', action='write')
!    write(20, '(A)') 'time timeFormation M_vir R_vir R_max V_max tau V_maxSIDM' ! Write header line

    if (time > timeFormation) then
!            call darkMatterProfile%floatRank0MetaPropertySet(self%VmaxSIDMID,darkMatterProfile%floatRank0MetaPropertyGet(self%VmaxSIDMID)+self%darkMatterProfileDMO_%circularVelocityMaximum(node))

            VmaxSIDM = darkMatterProfile%floatRank0MetaPropertyGet(self%VmaxSIDMID)+self%darkMatterProfileDMO_%circularVelocityMaximum(node)
            tc = get_tc(self, node, self%darkMatterProfileDMO_%circularVelocityMaximum(node), self%darkMatterProfileDMO_%radiusCircularVelocityMaximum(node), VmaxSIDM)
            call darkMatterProfile%floatRank0MetaPropertyRate(self%tauID, 1.0d0/tc)

            tau = darkMatterProfile%floatRank0MetaPropertyGet(self%tauID)
            dvdt = dvmaxt(tau, self%darkMatterProfileDMO_%circularVelocityMaximum(node)) * (1.0d0) / tc
!            call self%differentialVmaxSolveAnalytics(node)
            call darkMatterProfile%floatRank0MetaPropertyRate(self%VmaxSIDMID, dvdt)

            print *, 'Read Here!'
            print *, time, timeFormation, tc, tau, self%darkMatterProfileDMO_%circularVelocityMaximum(node), darkMatterProfile%floatRank0MetaPropertyGet(self%VmaxSIDMID)+self%darkMatterProfileDMO_%circularVelocityMaximum(node)

    else
!            print *, 'test1'
            call darkMatterProfile%floatRank0MetaPropertyRate(self%tauID, 0.0d0)
!            print *, 'test2'
            call darkMatterProfile%floatRank0MetaPropertyRate(self%VmaxSIDMID, 0.0d0)
    end if
    
    !call darkMatterProfile%floatRank0MetaPropertySet(self%VmaxSIDMID,darkMatterProfile%floatRank0MetaPropertyGet(self%VmaxSIDMID)+self%darkMatterProfileDMO_%circularVelocityMaximum(node))

    return
  end subroutine SIDMParametriCalculateTauDifferentialEvolution

!  subroutine SIDMParametriDifferentialVmaxAnalytics(self, node)
!    !!{
!    Mark analytically-solvable properties.
!    !!}
!    use :: Galacticus_Nodes, only : nodeComponentBasic, nodeComponentDarkMatterProfile, treeNode
!    class(nodeOperatorSIDMParametric), intent(inout) :: self
!    class(nodeComponentDarkMatterProfile), pointer :: darkMatterProfile
!    type(treeNode), intent(inout), target :: node

!    darkMatterProfile => node%darkMatterProfile()
!    call darkMatterProfile%floatRank0MetaPropertyAnalytic(self%VmaxSIDMID)
!    return
!  end subroutine SIDMParametriDifferentialVmaxAnalytics


!  subroutine SIDMParametriDifferentialVmaxSolveAnalytics(self, node)

!    use :: Galacticus_Nodes, only : nodeComponentBasic, nodeComponentDarkMatterProfile, treeNode
!    class(nodeOperatorSIDMParametric), intent(inout) :: self
!    class(nodeComponentDarkMatterProfile), pointer :: darkMatterProfile
!    type(treeNode), intent(inout), target :: node
!    double precision :: tau, VmaxSIDM, tc, dvdt

!    tau = darkMatterProfile%floatRank0MetaPropertyGet(self%tauID)
!    VmaxSIDM = darkMatterProfile%floatRank0MetaPropertyGet(self%VmaxSIDMID)+self%darkMatterProfileDMO_%circularVelocityMaximum(node)
!    tc = get_tc(self, node, self%darkMatterProfileDMO_%circularVelocityMaximum(node), self%darkMatterProfileDMO_%radiusCircularVelocityMaximum(node), VmaxSIDM)
!    dvdt = dvmaxt(tau, self%darkMatterProfileDMO_%circularVelocityMaximum(node)) * (1.0d0) / tc

!    call darkMatterProfile%floatRank0MetaPropertyRate(self%VmaxSIDMID, dvdt)
!    return
!  end subroutine SIDMParametriDifferentialVmaxSolveAnalytics 


  double precision function get_tc(self, node, Vmax, Rvmax, VmaxSIDM)
    !!{
    Evaluating tc based on Eq. 2.2 from Yang et al. 2024:https://arxiv.org/pdf/2305.16176
    !!}
    use Numerical_Constants_Math, only: Pi
!    use Numerical_Constants_Physical, only: gravitationalConstant
    use :: Numerical_Constants_astronomical, only : gravitationalConstantGalacticus
    use :: Error, only : Error_Report
    use :: Dark_Matter_Particles, only : darkMatterParticleSelfInteractingDarkMatter

    class(nodeOperatorSIDMParametric), intent(inout) :: self
    type(treeNode), intent(in) :: node
    double precision, intent(in) :: Vmax, Rvmax, VmaxSIDM

    double precision :: sigmaeff, reff, rhoeff, C = 0.75, gravitationalConstant !GG = 4.30073e-6, C = 0.75, pi = 3.1415926535897932384626433832795d0
    !type(effectiveCrossSection) :: sigmaeff

    print *, 'Vmax: ', Vmax
    print *, 'VmaxSIDM: ', VmaxSIDM

    select type (darkMatterParticle_ => self%darkMatterParticle_)
    class is (darkMatterParticleSelfInteractingDarkMatter)
       print *, 'darkMatterParticle_ is of type SIDM'
       sigmaeff = darkMatterParticle_%effectiveCrossSection(VmaxSIDM)
!       sigmaeff = darkMatterParticle_%effectiveCrossSection(Vmax)
    class default
       call Error_Report('unexpected class'//{introspection:location})
    end select

    print *, 'right after sigma effective calculations...', sigmaeff
    gravitationalConstant = gravitationalConstantGalacticus*1e3
    print *, 'gravitaional constant: ', gravitationalConstant
    print *, 'Rmax: ', Rvmax
    print *, 'Vmax: ', Vmax

!    write(20, '(F20.6, 2X, F20.6)') 0.0d0, sigmaeff 

    !we need the conversion of Rvmax to kpc rather than Mpc
    reff = Rvmax*1e3 / 2.1626 
    rhoeff = (Vmax / (1.648 * reff))**2 / gravitationalConstant
!    sigmaeff = self%darkMatterParticle_%effectiveCrossSection(Vmax)

    print *, 'reff: ', reff
    print *, 'rhoeff: ', rhoeff

!    get_tc = (150.0d0 / C) * (1.0d0 / (sigmaeff * rhoeff * reff)) * (4.0d0 * Pi * gravitationalConstant * rhoeff) ** (-0.5)
    !the constant 2.09e-10 is multiplied for unit conversion of sigma
    get_tc = (150.0d0 / C) * (1.0d0 / (sigmaeff * 2.09e-10 * rhoeff * reff)) * (1.0d0 / sqrt(4.0d0 * Pi * gravitationalConstant * rhoeff))

    print *, 'get_tc value: ', get_tc

  end function get_tc

  double precision function dvmaxt(tau, Vmaxt)
    double precision, intent(in   ) :: tau, Vmaxt
    double precision :: tau_local

    ! Create a local copy of tau
    tau_local = tau

    if (tau_local > 1.1d0) then
      tau_local = 1.1d0
    end if

    dvmaxt = 0.17774902d0 - 13.195824689999998d0 * tau_local ** 2 + 66.62092676d0 * tau_local ** 3 - 94.33706049999999d0 * tau_local ** 4 + 63.53766111d0 * tau_local ** 6 - 21.925108889999997d0 * tau_local ** 8
    dvmaxt = dvmaxt * Vmaxt
  end function dvmaxt

  double precision function drmaxt(tau, Rmaxt)
    double precision, intent(in   ) :: tau, Rmaxt
    double precision :: tau_local

    tau_local = tau

    if (tau_local > 1.1d0) then
      tau_local = 1.1d0
    end if

    drmaxt = 0.00762288d0 - 1.43996392d0 * tau_local + 1.01282643d0 * tau_local ** 2 - 0.55015288d0 * tau_local ** 3
    drmaxt = drmaxt * Rmaxt 
  end function drmaxt

  double precision function Rmax_NFW(RmaxSIDM, tau)
    double precision, intent(in) :: RmaxSIDM, tau

    Rmax_NFW = RmaxSIDM / (1 + 0.007623d0 * tau - 0.7200d0 * tau ** 2 + 0.3376d0 * tau ** 3 - 0.1375d0 * tau ** 4)
  end function Rmax_NFW

  double precision function Vmax_NFW(VmaxSIDM, tau)
    double precision, intent(in) :: VmaxSIDM, tau

    Vmax_NFW = VmaxSIDM / (1 + 0.1777d0 * tau - 4.399d0 * tau ** 3 + 16.66d0 * tau ** 4 - 18.87d0 * tau ** 5 + 9.077d0 * tau ** 7 - 2.436d0 * tau ** 9)
  end function Vmax_NFW

  double precision function r_s0(Rmax)
    double precision, intent(in) :: Rmax

    r_s0 = Rmax / 2.163d0
  end function r_s0

  double precision function rho_s0(Rs, Vmax)
    use Numerical_Constants_Math, only: Pi
    use Numerical_Constants_Physical, only: gravitationalConstant
    double precision, intent(in) :: Rs, Vmax

    rho_s0 = Vmax ** 2 / (0.465d0 ** 2 * 4.0d0 * Pi * gravitationalConstant * Rs ** 2)
  end function rho_s0

  double precision function get_rho_s(rho_s0, tau)
    double precision, intent(in) :: rho_s0, tau

    get_rho_s = rho_s0 * (2.033d0 + 0.7381d0 * tau + 7.264d0 * tau ** 5 - 12.73d0 * tau ** 7 + 9.915d0 * tau ** 9 + (1.0d0 + 2.033d0) * log(tau + 0.001d0) / log(0.001d0))
  end function get_rho_s

  double precision function get_r_s(r_s0, tau)
    double precision, intent(in) :: r_s0, tau

    get_r_s = r_s0 * (0.7178d0 - 0.1026d0 * tau + 0.2474d0 * tau ** 2 - 0.4079d0 * tau ** 3 + (1.0d0 - 0.7178d0) * log(tau + 0.001d0) / log(0.001d0))
  end function get_r_s

  double precision function get_r_c(r_s0, tau)
    double precision, intent(in) :: r_s0, tau

    get_r_c = r_s0 * (2.555d0 * sqrt(tau) - 3.632d0 * tau + 2.131d0 * tau ** 2 - 1.415d0 * tau ** 3 + 0.4683d0 * tau ** 4)
  end function get_r_c

  double precision function SIDMDensityProfileIsolated(rho_s, r_s, r_c, r)
    double precision, intent(in) :: rho_s, r_s, r_c, r
    double precision :: beta = 4.0d0, term1, term2

    term1 = ((r ** beta + r_c ** beta) ** (1.0d0 / beta) / r_s)
    term2 = (1.0d0 + r / r_s) ** 2

    SIDMDensityProfileIsolated = rho_s / (term1 * term2)
  end function SIDMDensityProfileIsolated

  function getTauID_SIDMParametric(self) result(tauID)
    implicit none
    class(nodeOperatorSIDMParametric), intent(in) :: self
    integer :: tauID
    tauID = self%tauID
  end function getTauID_SIDMParametric

  function getVmaxSIDMID_SIDMParametric(self) result(VmaxSIDMID)
    implicit none
    class(nodeOperatorSIDMParametric), intent(in) :: self
    integer :: VmaxSIDMID
    VmaxSIDMID = self%VmaxSIDMID
  end function getVmaxSIDMID_SIDMParametric

  function getRmaxSIDMID_SIDMParametric(self) result(RmaxSIDMID)
    implicit none
    class(nodeOperatorSIDMParametric), intent(in) :: self
    integer :: RmaxSIDMID
    RmaxSIDMID = self%RmaxSIDMID
  end function getRmaxSIDMID_SIDMParametric


 
