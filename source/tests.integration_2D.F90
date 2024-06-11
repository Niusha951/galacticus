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
Contains a program to test integration routines.
!!}

program Test_Integration2D
  !!{
  Tests that numerical integration routines work.
  !!}
  use :: Display                   , only : displayVerbositySet, verbosityLevelStandard
  use :: Numerical_Constants_Math  , only : Pi
!!  use :: Numerical_Integration     , only : integrator
  use :: SIDM_Effective_Cross_Section_Integrator, only : SIDMEffectiveCrossSectionIntegrator2D
  use :: Test_Integration2D_Functions, only : Integrand1         , Integrand2            , Integrand3          , Integrand4
  use :: Unit_Tests                , only : Assert             , Unit_Tests_Begin_Group, Unit_Tests_End_Group, Unit_Tests_Finish
  implicit none
  double precision                                                        :: integral
  double precision                                       , dimension(2,2) :: boundaries
  type            (SIDMEffectiveCrossSectionIntegrator2D), allocatable    :: integrator_

  ! Set verbosity level.
  call displayVerbositySet(verbosityLevelStandard)

  ! Begin unit tests.
  call Unit_Tests_Begin_Group("Numerical integration")

  ! Test simple integrations.
  allocate(integrator_)
  boundaries (1,:) = (/-1.0d0,1.0d0/)
  boundaries (2,:) = (/-1.0d0,1.0d0/)
  integrator_=SIDMEffectiveCrossSectionIntegrator2D(Integrand1,toleranceRelative=1.0d-6)
  integral=integrator_%integrate(boundaries)
  deallocate(integrator_)
  call Assert("integrate f(x,y)=sin(x^2) cos(y)          from x=-1……1, y=-1?~@??~@?1"          ,integral,1.04433d0             ,relTol=1.0d-6)

  allocate(integrator_)
  boundaries (1,:) = (/-1.0d0,1.0d0/)
  boundaries (2,:) = (/-1.0d0,1.0d0/)
  integrator_=SIDMEffectiveCrossSectionIntegrator2D(Integrand2,toleranceAbsolute=1.0d-6)
  integral=integrator_%integrate(boundaries)
  deallocate(integrator_)
  call Assert("integrate f(x)=x^2 cos(y)     from x=-1?~@??~@?1, y=-1?~@??~@?1"          ,integral,1.12196d0             ,absTol=1.0d-6)

  allocate(integrator_)
  boundaries (1,:) = (/0.0d0,1.0d0/)
  boundaries (2,:) = (/-1.0d0,1.0d0/)
  integrator_=SIDMEffectiveCrossSectionIntegrator2D(Integrand3,toleranceRelative=1.0d-6)
  integral=integrator_%integrate(boundaries)
  deallocate(integrator_)
  call Assert("integrate f(x)=√x y^2      from x=0?~@??~@?1, y=-1?~@??~@?1"          ,integral,0.444444d0,relTol=1.0d-6)

  allocate(integrator_)
  boundaries (1,:) = (/0.0d0,1.0d0/)
  boundaries (2,:) = (/0.0d0,2.0d0/)
  integrator_=SIDMEffectiveCrossSectionIntegrator2D(Integrand4,toleranceRelative=1.0d-6)
  integral=integrator_%integrate(boundaries)
  deallocate(integrator_)
  call Assert("integrate f(x,y)=y/?~H~Zx  from x=0…1π and y=0…2",integral,4.0d0          ,relTol=1.0d-6)

  ! End unit tests.
  call Unit_Tests_End_Group()
  call Unit_Tests_Finish   ()

end program Test_Integration2D
