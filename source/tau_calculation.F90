!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023, 2024
!!    Andrew Benson <abenson@carnegiescience.edu>
!!
!! This file is part of Galacticus.

  <tau_calculation name="tau_calculation">
    <description>
    </description>
  </tau_calculation>

  module MyNewClassModule
    implicit none
  
    ! Import the pre-written classes

    use :: PreWrittenClassA, only : PreWrittenClassAType
    use :: PreWrittenClassB, only : PreWrittenClassBType

  
    ! Define the new class
    type :: tau_calculation
      ! Components specific to the new class
      integer :: myInteger
      real    :: myReal
    
      ! Components from pre-written classes
      type(PreWrittenClassAType) :: classAInstance
      type(PreWrittenClassBType) :: classBInstance
    end type MyNewClassType
  
  contains

    ! Constructor for MyNewClassType



    subroutine MyNewClass(self, node, integerInput, realInput)

      use :: Galacticus_Nodes, only : treeNode
      use :: Dark_Matter_Halo_Mass_Accretion_Histories, only : darkMatterHaloMassAccretionHistoryClass

      type(tau_calculation), intent(  out) :: self
      type(treeNode),        intent(inout) :: node
      type(treeNode),        pointer       :: nodeFinal   
      type(treeNode),        pointer       :: nodeWork

      class(darkMatterHaloMassAccretionHistoryClass), intent(in   ), target :: darkMatterHaloMassAccretionHistory_
      class(nodeComponentBasic),                                    pointer :: basic
      class(nodeComponentDarkMatterProfile),                        pointer :: darkMatterProfile
      class(darkMatterProfileDMO),                                  pointer :: darkMatterProfileDMO_ 

      double precision                     :: timeFormation
      double precision                     :: formationMassFraction = 50   
      double precision                     :: tau, time, timePrevious


      integer, intent(in) :: integerInput
      real, intent(in) :: realInput
   
      if (associated(node%firstChild)) return
      nodeFinal => node
      do while (nodeFinal%isPrimaryProgenitor())
        nodeFinal => nodeFinal%parent
      end do
 
      timeFormation=Dark_Matter_Halo_Formation_Time(nodeFinal,formationMassFraction,self%darkMatterHaloMassAccretionHistory_)

      tau=0.0d0
      timePrevious=0.0d0
      nodeWork => node
      do while (associated(nodeWork))
         ! Check if this node exists after the formation time. (If it is before, we do
         ! not increment τ.)
         basic => nodeWork%basic()
         time=basic%time()
         if (time > timeFormation) then
            ! Compute tc and increment τ.
            tc=get_tc(nodeWork, darkMatterProfileDMO_%circularVelocityMaximum(nodeWork), darkMatterProfileDMO_%radiusCircularVelocityMaximum(nodeWork))
            tau=tau+(time-timePrevious)/tc
         end if
         ! Store the value of τ in the dark matter profile component of this node.
         darkMatterProfile => node%darkMatterProfile()
         call darkMatterProfile%floatRank0MetaPropertySet(self%tauID,tau)
         ! Update the timePrevious variable so that we can compute the time interval
         ! on the next pass through this loop.
         timePrevious=time
         ! Move to the next node.
         if (associated(nodeWork,nodeFinal)) then
            ! We've reached the final node in the branch - set our worker node to null
            ! so we exit the loop.
            nodeWork => null()
         else
            ! Move up the branch.
            nodeWork => nodeWork%parent
         end if
      end do


    end subroutine MyNewClass

  contains
    
    double precision function get_tc(node, Vmax, Rvmax)
      !!{
      Evalluating tc based on Eq. 2.2 from Yang et al. 2024:https://arxiv.org/pdf/2305.16176
      !!}
  
      double precision, intent(in   ) :: Vmax, Rvmax

      double precision                :: reff, rhoeff, GG = 4.30073*1e-6, C = 0.75, pi = 3.1415926535897932384626433832795d0


      reff = Rvmax/2.1626    
      rhoeff = (Vmax/1.648/reff)**2 / GG
      sigmaeff = sigma_eff()

      get_tc = (150.0d0/C) * (1.0d0/(sigmaeff * rhoeff * reff)) * (4.0d0*pi*GG*rhoeff)**(-0.5)


    end function get_tc

    double precision function sigma_eff()

    end function sigma_eff

    double precision function VelocityDependentCrossSectionDifferential(self,theta,velocityRelative)
    !!{
    Return the differential self-interaction cross section, $\mathrm{d}\sigma/\mathrm{d}\theta$, in units of cm$^2$ g$^{-1}$ ster$^{-1}$, of a self-interacting dark matter particle.
    !!}
    implicit none
    class           (darkMatterParticleSIDMVelocityDependent), intent(inout) :: self
    double precision                                         , intent(in   ) :: velocityRelative
    double precision                                         , intent(in   ) :: theta

    ! Currently isotropic scattering is assumed.
    VelocityDependentCrossSectionDifferential=+(self%velocityCharacteristic)**4
* (self%sigma0) *0.5d0 * sin(theta) / ((self%velocityCharacteristic)**2 + 0.5d0
* velocityRelative**2 * (1 - cos(theta)))**2
    return
  end function sidmVelocityDependentCrossSectionSelfInteractionDifferential




  end module MyNewClassModule

