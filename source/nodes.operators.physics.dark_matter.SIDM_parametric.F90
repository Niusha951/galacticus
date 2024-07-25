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
 
     integer                                                 :: tauID, VmaxSIDMID, RmaxSIDMID
   contains
     !![
     <methods>
       <method method="calculateTau" description="Computes tau, VmaxSIDM and RmaxSIDM for each node."/>
     </methods>
     !!]
     final     ::                                        SIDMParametricDestructor
     procedure :: nodeTreeInitialize                  => SIDMParametricNodeTreeInitialize
     procedure :: calculateTau                        => SIDMParametriCalculateTau
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
    <addMetaProperty component="darkMatterProfile" name="tau"      id="self%tauID"      isEvolvable="no" isCreator="yes"/>
    <addMetaProperty component="darkMatterProfile" name="VmaxSIDM" id="self%VmaxSIDMID" isEvolvable="no"  isCreator="yes"/>
    <addMetaProperty component="darkMatterProfile" name="RmaxSIDM" id="self%RmaxSIDMID" isEvolvable="no"  isCreator="yes"/>
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
    implicit none
    class(nodeOperatorSIDMParametric), intent(inout), target  :: self
    type (treeNode                  ), intent(inout), target  :: node
    
    call self%nodeInitialize(node)
    return
  end subroutine SIDMParametricNodeTreeInitialize
 
  subroutine SIDMParametriCalculateTau(self, node)
    use Galacticus_Nodes, only: treeNode, nodeComponentBasic, nodeComponentDarkMatterProfile
!    use Dark_Matter_Halo_Mass_Accretion_Histories, only: darkMatterHaloMassAccretionHistoryClass
!    use Dark_Matter_Profiles_DMO , only : darkMatterProfileDMOClass
    use :: Dark_Matter_Halo_Formation_Times, only : Dark_Matter_Halo_Formation_Time

    type(treeNode), intent(inout), target :: node
    type(treeNode), pointer :: nodeFinal
    type(treeNode), pointer :: nodeWork
!    type(darkMatterProfileDMO), pointer :: darkMatterProfileDMO_

    class(nodeOperatorSIDMParametric), intent(inout) :: self
!    class(darkMatterHaloMassAccretionHistoryClass), target :: darkMatterHaloMassAccretionHistory_
!    class(darkMatterHaloMassAccretionHistoryClass), pointer :: darkMatterHaloMassAccretionHistory_
    class(nodeComponentBasic), pointer :: basic
    class(nodeComponentDarkMatterProfile), pointer :: darkMatterProfile
!    class(darkMatterProfileDMOClass), pointer :: darkMatterProfileDMO_

    double precision :: timeFormation
    double precision :: formationMassFraction = 0.5d0
    double precision :: tau, time, timePrevious, tc
    double precision :: RmaxNFW0, VmaxNFW0, r_sNFW0, rho_sNFW0, rho_s, r_s, r_c
    double precision :: VmaxSIDM=0.0d0, RmaxSIDM=0.0d0, VmaxSIDMPrevious, RmaxSIDMPrevious, VmaxCDM, RmaxCDM, VmaxCDMPrevious, RmaxCDMPrevious
    double precision :: dtr

!    print *, 'Does this node have children?'
!    if (associated(node%firstChild)) then 
!       print *, 'Yes!'
!    else 
!       print *, 'No!'
!    end if

    basic => node%basic()
!    print *, 'time associated to this node: ', basic%time()
!    print *, 'mass associated to this node: ', basic%mass()

    ! Check if the node has children, if so, return
    if (associated(node%firstChild)) return

    nodeFinal => node
    do while (nodeFinal%isPrimaryProgenitor())
!      print *, 'repeat is num. of loop!'
      nodeFinal => nodeFinal%parent
    end do

    basic => nodefinal%basic()
!    print *, 'time associated to nodefinal: ', basic%time()
!    print *, 'mass associated to nodefinal: ', basic%mass()
!    print *, 'End of cheching for children of nodes!'

    ! Calculate the formation time
!    timeFormation = Dark_Matter_Halo_Formation_Time(node, formationMassFraction, self%darkMatterHaloMassAccretionHistory_)

    !Just for test!
    timeFormation = 2.1903d0

    print *, 'Formation time: ', timeFormation

    basic => node%basic()
    print *, 'time associated to node: ', basic%time()
    print *, 'mass associated to node: ', basic%mass()

    tau = 0.0d0
    timePrevious = 0.0d0
    VmaxSIDMPrevious = 0.0d0
    RmaxSIDMPrevious = 0.0d0
!    nodeWork => node
    nodework => nodefinal

    print *, 'Starting the while loop.'

    ! Open the output file
    open(unit=20, file='/home/nahvazi/Galacticus_SIDM_parametric/gal/galacticus/SIDM_output_data.txt', status='replace', action='write')
    write(20, '(A)') 'time tau VmaxSIDM RmaxSIDM' ! Write header line
!     open(unit=20, file='/home/nahvazi/Galacticus_SIDM_parametric/gal/galacticus/SIDM_output_data_sigmaEff.txt', status='replace', action='write')
!     write(20, '(A)') 'time sigma_effective'

    do while (associated(nodeWork))
      basic => nodeWork%basic()
      time = basic%time()

      print *, 'Initiate basic, and read time for nodeWork: ', time

      if (time > timeFormation) then
        print *, 'Test1: inside if loop.'
        ! Compute tc and increment ?~D.
        tc = get_tc(self, nodeWork, self%darkMatterProfileDMO_%circularVelocityMaximum(nodeWork), self%darkMatterProfileDMO_%radiusCircularVelocityMaximum(nodeWork), VmaxSIDMPrevious)
        print *, 'after tc: ', tc

!        tau = tau + (time - timePrevious) / tc
        tau = (time - timeFormation) / tc

        print *, 'after tau: ', tau

        ! Check if dtr is greater than 0.05 and print a message if true
        dtr = (time - timePrevious) / tc
        if (dtr > 0.05d0) then
          print *, 'Gravothermal evolution too fast, dtr: ', dtr
        end if

        VmaxCDM = self%darkMatterProfileDMO_%circularVelocityMaximum(nodeWork)
        RmaxCDM = self%darkMatterProfileDMO_%radiusCircularVelocityMaximum(nodeWork)

        print *, 'VmaxCDM: ', VmaxCDM
        print *, 'VmaxCDMPrevious: ', VmaxCDMPrevious
        print *, 'subtract: ', VmaxCDM - VmaxCDMPrevious

        VmaxSIDM = (VmaxCDM - VmaxCDMPrevious) + VmaxSIDMPrevious + dvmaxt(tau, self%darkMatterProfileDMO_%circularVelocityMaximum(nodeWork)) * (time - timePrevious) / tc
!        VmaxSIDM = VmaxCDM + dvmaxt(tau, self%darkMatterProfileDMO_%circularVelocityMaximum(nodeWork)) * (time - timePrevious) / tc
        RmaxSIDM = (RmaxCDM - RmaxCDMPrevious) + RmaxSIDMPrevious + drmaxt(tau, self%darkMatterProfileDMO_%radiusCircularVelocityMaximum(nodeWork)) * (time - timePrevious) / tc
!        RmaxSIDM = RmaxCDM + drmaxt(tau, self%darkMatterProfileDMO_%radiusCircularVelocityMaximum(nodeWork)) * (time - timePrevious) / tc

        print *, 'VmaxSIDM: ', VmaxSIDM
        print *, 'RmaxSIDM: ', RmaxSIDM


        ! Store the value of tau, VmaxSIDM, and RmaxSIDM in the dark matter
        ! profile component of this node.
        darkMatterProfile => node%darkMatterProfile()
        call darkMatterProfile%floatRank0MetaPropertySet(self%tauID, tau)
        call darkMatterProfile%floatRank0MetaPropertySet(self%VmaxSIDMID, VmaxSIDM)
        call darkMatterProfile%floatRank0MetaPropertySet(self%RmaxSIDMID, RmaxSIDM)

      end if

      ! Write the data to the output file
      write(20, '(F20.6, 2X, F20.6, 2X, F20.6, 2X, F20.6)') time, tau, VmaxSIDM, RmaxSIDM
!      if (time<timeFormation) then
!         write(20, '(F20.6, 2X, F20.6)') time, 0.0d0 
!      end if

      RmaxNFW0 = Rmax_NFW(RmaxSIDM, tau)
      VmaxNFW0 = Vmax_NFW(VmaxSIDM, tau)


      r_sNFW0 = r_s0(RmaxNFW0)
      rho_sNFW0 = rho_s0(r_sNFW0, VmaxNFW0)

      rho_s = get_rho_s(rho_sNFW0, tau)
      r_s = get_r_s(r_sNFW0, tau)
      r_c = get_r_c(r_sNFW0, tau)

!      print *, 'Testing VmaxSIDMPrevious: ', VmaxSIDMPrevious

      if (time > timeFormation) then
        print *, 'check if inside the correct if condition: time > timeFormation'
        VmaxSIDMPrevious = VmaxSIDM
        RmaxSIDMPrevious = RmaxSIDM
        VmaxCDMPrevious = VmaxCDM
        RmaxCDMPrevious = RmaxCDM
      else
        print *, 'Test if it goes through the correct section of if/else'
        VmaxSIDMPrevious = self%darkMatterProfileDMO_%circularVelocityMaximum(nodeWork)
!        print *, 'Testing VmaxSIDMPrevious after association: ', VmaxSIDMPrevious
        RmaxSIDMPrevious = self%darkMatterProfileDMO_%radiusCircularVelocityMaximum(nodeWork)
        VmaxCDMPrevious = self%darkMatterProfileDMO_%circularVelocityMaximum(nodeWork)
        RmaxCDMPrevious = self%darkMatterProfileDMO_%radiusCircularVelocityMaximum(nodeWork)

      end if

      print *, 'End of setting previous Vmax and Rmax', VmaxSIDMPrevious, RmaxSIDMPrevious

      ! Update the timePrevious variable so that we can compute the time
      ! interval on the next pass through this loop.
      timePrevious = time
      !VmaxCDMPrevious = VmaxCDM
      !RmaxCDMPrevious = RmaxCDM

      ! Move to the next node.
!      if (associated(nodeWork, nodeFinal)) then
      if (associated(nodeWork, node)) then
        ! We've reached the final node in the branch - set our worker node to
        ! null
        ! so we exit the loop.
        nodeWork => null()
      else
        ! Move up the branch.
!        nodeWork => nodeWork%parent
        nodeWork => nodeWork%firstChild
      end if
    end do

    return
  end subroutine SIDMParametriCalculateTau

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


 
