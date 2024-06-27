!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023, 2024
!!    Andrew Benson <abenson@carnegiescience.edu>
!!
!! This file is part of Galacticus.

<tauCalculation name="tauCalculation">
  <description>
  </description>
</tauCalculation>

module tauCalculationClassModule
  use Galacticus_Nodes, only: treeNode
  use Dark_Matter_Halo_Mass_Accretion_Histories, only: darkMatterHaloMassAccretionHistoryClass
  use dark_matter_particle_self_interacting, only: effectiveCrossSection
  !use :: Dark_Matter_Particles, only : darkMatterParticleSelfInteractingDarkMatterConstant


  implicit none

  ! Define the new class
  type :: tauCalculation
    double precision :: VmaxSIDM
    double precision :: RmaxSIDM
  end type tauCalculation

contains

  subroutine tauCalculationClass(self, node)

    type(tauCalculation), intent(out) :: self
    type(treeNode), intent(inout) :: node
    type(treeNode), pointer :: nodeFinal   
    type(treeNode), pointer :: nodeWork

    class(darkMatterHaloMassAccretionHistoryClass), intent(in), target :: darkMatterHaloMassAccretionHistory_
    class(nodeComponentBasic), pointer :: basic
    class(nodeComponentDarkMatterProfile), pointer :: darkMatterProfile
    class(darkMatterProfileDMO), pointer :: darkMatterProfileDMO_ 

    double precision :: timeFormation
    double precision :: formationMassFraction = 50   
    double precision :: tau, time, timePrevious, tc
    double precision :: RmaxNFW0, VmaxNFW0, r_sNFW0, rho_sNFW0, rho_s, r_s, r_c
    double precision :: VmaxSIDM, RmaxSIDM, VmaxSIDMPrevious, RmaxSIDMPrevious
    double precision :: dtr 

    ! Removed intent declarations that were misplaced
    !integer :: i
    !double precision :: temp1, temp2, temp3, temp4, temp5, temp6, temp7, temp8

    ! Check if the node has children, if so, return
    if (associated(node%firstChild)) return

    nodeFinal => node
    do while (nodeFinal%isPrimaryProgenitor())
      nodeFinal => nodeFinal%parent
    end do

    ! Calculate the formation time
    timeFormation = Dark_Matter_Halo_Formation_Time(nodeFinal, formationMassFraction, self%darkMatterHaloMassAccretionHistory_)

    tau = 0.0d0
    timePrevious = 0.0d0
    nodeWork => node

    do while (associated(nodeWork))
      basic => nodeWork%basic()
      time = basic%time()

      if (time > timeFormation) then
        ! Compute tc and increment τ.
        tc = get_tc(nodeWork, darkMatterProfileDMO_%circularVelocityMaximum(nodeWork), darkMatterProfileDMO_%radiusCircularVelocityMaximum(nodeWork))
        tau = tau + (time - timePrevious) / tc

        ! Check if dtr is greater than 0.05 and print a message if true
        dtr = (time - timePrevious) / tc
        if (dtr > 0.05d0) then
          print *, 'Gravothermal evolution too fast, dtr: ', dtr
        end if

        VmaxSIDM = VmaxSIDMPrevious + dvmaxt(tau, darkMatterProfileDMO_%circularVelocityMaximum(nodeWork)) * (time - timePrevious) / tc
        RmaxSIDM = RmaxSIDMPrevious + drmaxt(tau, darkMatterProfileDMO_%radiusCircularVelocityMaximum(nodeWork)) * (time - timePrevious) / tc
      end if

      ! Store the value of τ, VmaxSIDM, and RmaxSIDM in the dark matter profile component of this node.
      darkMatterProfile => node%darkMatterProfile()
      call darkMatterProfile%floatRank0MetaPropertySet(self%tauID, tau)
      call darkMatterProfile%floatRank0MetaPropertySet(self%VmaxSIDMID, VmaxSIDM)
      call darkMatterProfile%floatRank0MetaPropertySet(self%RmaxSIDMID, RmaxSIDM)

      RmaxNFW0 = Rmax_NFW(RmaxSIDM, tau)
      VmaxNFW0 = Vmax_NFW(VmaxSIDM, tau)

      r_sNFW0 = r_s0(RmaxNFW0)
      rho_sNFW0 = rho_s0(r_s0, VmaxNFW0)

      rho_s = get_rho_s(rho_sNFW0, tau)
      r_s = get_r_s(r_sNFW0, tau)
      r_c = get_r_c(r_sNFW0, tau)

      if (time > timeFormation) then
        VmaxSIDMPrevious = VmaxSIDM
        RmaxSIDMPrevious = RmaxSIDM
      else
        VmaxSIDMPrevious = darkMatterProfileDMO_%circularVelocityMaximum(nodeWork)
        RmaxSIDMPrevious = darkMatterProfileDMO_%radiusCircularVelocityMaximum(nodeWork)
      end if

      ! Update the timePrevious variable so that we can compute the time interval on the next pass through this loop.
      timePrevious = time

      ! Move to the next node.
      if (associated(nodeWork, nodeFinal)) then
        ! We've reached the final node in the branch - set our worker node to null
        ! so we exit the loop.
        nodeWork => null()
      else
        ! Move up the branch.
        nodeWork => nodeWork%parent
      end if
    end do

    return
  end subroutine tauCalculationClass

  double precision function get_tc(node, Vmax, Rvmax)
    !!{
    Evaluating tc based on Eq. 2.2 from Yang et al. 2024:https://arxiv.org/pdf/2305.16176
    !!}
    use Dark_Matter_Particle_Self_Interacting_Dark_Matter, only: effectiveCrossSection

    type(treeNode), intent(in) :: node
    double precision, intent(in) :: Vmax, Rvmax

    double precision :: reff, rhoeff, GG = 4.30073e-6, C = 0.75, pi = 3.1415926535897932384626433832795d0
    type(effectiveCrossSection) :: sigmaeff


    reff = Rvmax / 2.1626    
    rhoeff = (Vmax / 1.648 / reff) ** 2 / GG
    sigmaeff = effectiveCrossSection(Vmax)

    get_tc = (150.0d0 / C) * (1.0d0 / (sigmaeff * rhoeff * reff)) * (4.0d0 * pi * GG * rhoeff) ** (-0.5)
  end function get_tc

  double precision function dvmaxt(tau, Vmaxt)
    double precision, intent(in) :: tau, Vmaxt

    if (tau > 1.1d0) then
      tau = 1.1d0
    end if

    dvmaxt = 0.17774902d0 - 13.195824689999998d0 * tau ** 2 + 66.62092676d0 * tau ** 3 - 94.33706049999999d0 * tau ** 4 + 63.53766111d0 * tau ** 6 - 21.925108889999997d0 * tau ** 8
    dvmaxt = dvmaxt * Vmaxt
  end function dvmaxt

  double precision function drmaxt(tau, Rmaxt)
    double precision, intent(in) :: tau, Rmaxt

    if (tau > 1.1d0) then
      tau = 1.1d0
    end if

    drmaxt = 0.00762288d0 - 1.43996392d0 * tau + 1.01282643d0 * tau ** 2 - 0.55015288d0 * tau ** 3
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

end module tauCalculationClassModule

