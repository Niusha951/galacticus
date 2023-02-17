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

program Tests_SIDM_velocity_dependent
  !!{
  Tests velocity dependent cross section calculations.
  !!}
  use :: Dark_Matter_Particles, only : darkMatterParticleCDM         , darkMatterParticleSIDMVelocityDependent
  use :: Display              , only : displayVerbositySet           , verbosityLevelStandard
  use :: Unit_Tests           , only : Assert                        , Unit_Tests_Begin_Group         , Unit_Tests_End_Group, Unit_Tests_Finish
  implicit none
  double precision                                 , dimension(2), parameter :: velocity            =[14.0d0, 100.0d0]
  double precision                                 , dimension(2), parameter :: crossSectionCheck   =[100.0d0, 2.0d0]
  type            (darkMatterParticleCDM          )                          :: darkMatterParticle_
  type            (darkMatterParticleSIDMVelocityDependent)                  :: velocityDependentCrossSection_
  character       (len=1024                       )                          :: message
  integer                                                                    :: i
  double precision                                                           :: crossSection

  ! Set verbosity level.
  call displayVerbositySet(verbosityLevelStandard)
  ! Begin unit tests.
  call Unit_Tests_Begin_Group("SIDM: velocity dependent")
  ! Test growth factor in an Einstein-de Sitter universe. Growth factor should equal the expansion factor.
  !![
  <referenceConstruct object="darkMatterParticle_">
   <constructor>
    darkMatterParticleCDM ()
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="velocityDependentCrossSection_">
   <constructor>
    darkMatterParticleSIDMVelocityDependent (                                      &amp;
     &amp;                          velocityCharacteristic = 24.33d0             , &amp;
     &amp;                          sigma0                 = 147.1d0             , &amp; 
     &amp;                          darkMatterParticle_    = darkMatterParticle_   &amp;
     &amp;                                  )
   </constructor>
  </referenceConstruct>
  !!]
  do i=1,size(velocity)
     crossSection = velocityDependentCrossSection_%crossSectionSelfInteractionMomentumTransfer(velocityRelative = velocity(i))
     write (message,'(a,f6.1,a)') "cross section at velocity [v=",velocity(i),"]"
     call Assert(trim(message),crossSection,crossSectionCheck(i),relTol=1.0d-3)
  end do
  ! End unit tests.
  call Unit_Tests_End_Group()
  call Unit_Tests_Finish   ()
end program Tests_SIDM_velocity_dependent
