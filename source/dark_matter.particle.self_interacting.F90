!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023
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
Contains a module which implements a selfInteracting dark matter particle class.
!!}

  !![
  <darkMatterParticle name="darkMatterParticleSelfInteractingDarkMatter" abstract="yes">
   <description>Provides a selfInteracting dark matter particle.</description>
  </darkMatterParticle>
  !!]
  type, abstract, extends(darkMatterParticleClass) :: darkMatterParticleSelfInteractingDarkMatter
     !!{
     A selfInteracting dark matter particle class.
     !!}
     private
   contains
     !![
     <methods>
       <method description="Return the self-interaction cross section, $\sigma$, of the dark matter particle in units of cm$^2$ g$^{-1}$." method="crossSectionSelfInteraction" />
       <method description="Return the differential self-interaction cross section, $\mathrm{d}\sigma/\mathrm{d}\Omega$, of the dark matter particle in units of cm$^2$ g$^{-1}$ ster$^{-1}$." method="crossSectionSelfInteractionDifferential"/>
       <method description="Return the momentum transfer self-interaction cross section, $\sigma$, of the dark matter particle in units of cm$^2$ g$^{-1}$." method="crossSectionSelfInteractionMomentumTransfer" />
       <method description="Return the viscosity self-interaction cross section, $\sigma$, of the dark matter particle in units of cm$^2$ g$^{-1}$." method="crossSectionSelfInteractionViscosity" />
     </methods>
     !!]
     procedure(crossSectionSelfInteractionTemplate                ), deferred :: crossSectionSelfInteraction
     procedure(crossSectionSelfInteractionDifferentialTemplate    ), deferred :: crossSectionSelfInteractionDifferential
     procedure(crossSectionSelfInteractionMomentumTransferTemplate), deferred :: crossSectionSelfInteractionMomentumTransfer
     procedure(crossSectionSelfInteractionViscosityTemplate       ), deferred :: crossSectionSelfInteractionViscosity
  end type darkMatterParticleSelfInteractingDarkMatter

  abstract interface
    double precision function crossSectionSelfInteractionTemplate(self,velocityRelative)
      !!{
      Interface for self-interaction cross section, in units of cm$^2$ g$^{-1}$, of a self-interacting dark matter particle.
      !!}
      import darkMatterParticleSelfInteractingDarkMatter
      class(darkMatterParticleSelfInteractingDarkMatter), intent(inout) :: self
      double precision                                  , intent(in   ) :: velocityRelative
      !return
    end function crossSectionSelfInteractionTemplate

    double precision function crossSectionSelfInteractionDifferentialTemplate(self,theta,velocityRelative)
      !!{
      Interface for differential self-interaction cross section, $\mathrm{d}\sigma/\mathrm{d}\theta$, in units of cm$^2$ g$^{-1}$ ster$^{-1}$, of a self-interacting dark matter particle.
      !!}
      import darkMatterParticleSelfInteractingDarkMatter
      class (darkMatterParticleSelfInteractingDarkMatter), intent(inout) :: self
      double precision                                   , intent(in)    :: velocityRelative
      double precision                                   , intent(in   ) :: theta
    end function crossSectionSelfInteractionDifferentialTemplate

    double precision function crossSectionSelfInteractionMomentumTransferTemplate(self,velocityRelative)
      !!{
      Interface for momentum transfer self-interaction cross section, in units of cm$^2$ g$^{-1}$, of a self-interacting dark matter particle.
      !!}
      import darkMatterParticleSelfInteractingDarkMatter
      class(darkMatterParticleSelfInteractingDarkMatter), intent(inout) :: self
      double precision                                  , intent(in   ) :: velocityRelative
    end function crossSectionSelfInteractionMomentumTransferTemplate

    double precision function crossSectionSelfInteractionViscosityTemplate(self,velocityRelative)
      !!{
      Interface for viscosity self-interaction cross section, in units of cm$^2$ g$^{-1}$, of a self-interacting dark matter particle.
      !!}
      import darkMatterParticleSelfInteractingDarkMatter
      class(darkMatterParticleSelfInteractingDarkMatter), intent(inout) :: self
      double precision                                  , intent(in   ) :: velocityRelative
    end function crossSectionSelfInteractionViscosityTemplate

  end interface

