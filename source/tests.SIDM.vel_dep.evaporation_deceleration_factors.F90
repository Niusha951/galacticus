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

program Tests_SIDM_evaporation_deceleration_factors
  !!{
  Tests calculation of evaporation and deceleration factors based on Kummer2018.
  !!}
  use :: Dark_Matter_Profiles_DMO    , only : darkMatterProfileDMOIsothermal
  use :: Galactic_Structure          , only : galacticStructureStandard
  use :: Dark_Matter_Halo_Scales     , only : darkMatterHaloScaleVirialDensityContrastDefinition
  use :: Cosmology_Parameters        , only : cosmologyParametersSimple
  use :: Cosmology_Functions         , only : cosmologyFunctionsMatterLambda
  use :: Virial_Density_Contrast     , only : virialDensityContrastSphericalCollapseClsnlssMttrCsmlgclCnstnt
  use :: Dark_Matter_Profiles        , only : darkMatterProfileAdiabaticGnedin2004
  use :: Satellite_Evaporation_SIDM  , only : satelliteEvaporationSIDMKummer2018
  use :: Satellite_Deceleration_SIDM , only : satelliteDecelerationSIDMKummer2018
  use :: Dark_Matter_Profiles_Generic, only : nonAnalyticSolversNumerical
  use :: Dark_Matter_Particles       , only : darkMatterParticleCDM         , darkMatterParticleSIDMVelocityDependent
  use :: Display                     , only : displayVerbositySet           , verbosityLevelStandard
  use :: Unit_Tests                  , only : Assert                        , Unit_Tests_Begin_Group         , Unit_Tests_End_Group, Unit_Tests_Finish
  implicit none
  double precision                                 , dimension(6), parameter    :: velocity               =[14.0d0, 100.0d0, 24.33d0, 24.33d0, 24.33d0, 100.0d0]
  double precision                                 , dimension(6), parameter    :: x                      =[0.5d0, 0.5d0, 0.5d0, 1.5d0, 100.0d0, 100.0d0]
  double precision                                 , dimension(6), parameter    :: evaporationFactorCheck =[0.592163d0, 0.168609d0, 0.555556d0, -0.345794d0, -0.99975d0, -0.998805d0]
  double precision                                 , dimension(6), parameter    :: decelerationFactorCheck =[0.131778d0, 0.0766487d0, 0.129266d0, 0.580017d0, 0.999833d0, 0.998805d0]
  type            (darkMatterProfileDMOIsothermal ), pointer                    :: darkMatterProfileIsothermal_
  type            (galacticStructureStandard      ), pointer                    :: galacticStructureStandard_
  type            (darkMatterParticleSIDMVelocityDependent)                     :: darkMatterParticleSIDM_
  type            (darkMatterParticleCDM)                                       :: darkMatterParticle_
  type            (darkMatterHaloScaleVirialDensityContrastDefinition), pointer :: darkMatterHaloScale_
  type            (cosmologyParametersSimple      ), pointer                    :: cosmologyParameters_
  type            (cosmologyFunctionsMatterLambda ), pointer                    :: cosmologyFunctions_
  type            (virialDensityContrastSphericalCollapseClsnlssMttrCsmlgclCnstnt), pointer :: virialDensityContrast_
  type            (darkMatterProfileAdiabaticGnedin2004), pointer               :: darkMatterProfileAdiabaticGnedin2004_
  type            (satelliteEvaporationSIDMKummer2018), pointer                 :: satelliteEvaporationSIDM_
  type            (satelliteDecelerationSIDMKummer2018), pointer                :: satelliteDecelerationSIDM_
  character       (len=1024                       )                             :: message
  integer                                                                       :: i
  double precision                                                              :: evaporationFactor, decelerationFactor,  fractionBaryons = 0.15d0
  ! Set verbosity level.
  call displayVerbositySet(verbosityLevelStandard)
  ! Begin unit tests.
  call Unit_Tests_Begin_Group("SIDM: evaporation deceleration factors")
  ! Test evaporation and deceleration factors.
  ! Build required objects. Use a cosmology with a baryon fraction of 0.15 to match what was assumed in the reference run of the "contra" code.
  allocate(cosmologyParameters_                 )
  allocate(cosmologyFunctions_                  )
  allocate(virialDensityContrast_               )
  allocate(galacticStructureStandard_           )
  allocate(darkMatterHaloScale_                 )
  allocate(darkMatterProfileIsothermal_         )
  allocate(darkMatterProfileAdiabaticGnedin2004_)
  allocate(satelliteEvaporationSIDM_            )
  allocate(satelliteDecelerationSIDM_            )
  !![
  <referenceConstruct object="darkMatterParticle_">
   <constructor>
    darkMatterParticleCDM ()
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="darkMatterParticleSIDM_">
   <constructor>
    darkMatterParticleSIDMVelocityDependent                       (                                                      &amp;
     &amp;                                   velocityCharacteristic              = 24.33d0                             , &amp;
     &amp;                                   sigma0                              = 147.1d0                             , &amp;
     &amp;                                   darkMatterParticle_                 = darkMatterParticle_                   &amp;
     &amp;                                                        )
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="cosmologyParameters_"                 >
   <constructor>
    cosmologyParametersSimple                                     (                                                       &amp;
     &amp;                                   OmegaMatter                         = 0.3000d0                             , &amp;
     &amp;                                   OmegaBaryon                         = 0.3000d0*fractionBaryons             , &amp;
     &amp;                                   OmegaDarkEnergy                     = 0.7000d0                             , &amp;
     &amp;                                   temperatureCMB                      = 2.7800d0                             , &amp;
     &amp;                                   HubbleConstant                      =69.3000d0                               &amp;
     &amp;                                                        )
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="cosmologyFunctions_"                  >
   <constructor>
    cosmologyFunctionsMatterLambda                                (                                                       &amp;
     &amp;                                  cosmologyParameters_                =cosmologyParameters_                     &amp;
     &amp;                                                        )
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="virialDensityContrast_"               >
   <constructor>
    virialDensityContrastSphericalCollapseClsnlssMttrCsmlgclCnstnt(                                                       &amp;
     &amp;                                   tableStore                          =.true.                                , &amp;
     &amp;                                   cosmologyFunctions_                 =cosmologyFunctions_                     &amp;
     &amp;                                                        )
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="darkMatterHaloScale_"                 >
   <constructor>
    darkMatterHaloScaleVirialDensityContrastDefinition            (                                                       &amp;
     &amp;                                   cosmologyParameters_                =cosmologyParameters_                   ,&amp;
     &amp;                                   cosmologyFunctions_                 =cosmologyFunctions_                    ,&amp;
     &amp;                                   virialDensityContrast_              =virialDensityContrast_                  &amp;
     &amp;                                                        )
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="darkMatterProfileIsothermal_">
   <constructor>
    darkMatterProfileDMOIsothermal           (                                                                            &amp;
     &amp;                                   darkMatterHaloScale_                =darkMatterHaloScale_                    &amp;
     &amp;                                   )
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="galacticStructureStandard_"           >
   <constructor>
    galacticStructureStandard                (                                                                            &amp;
     &amp;                                   cosmologyFunctions_                 =cosmologyFunctions_                    ,&amp;
     &amp;                                   darkMatterHaloScale_                =darkMatterHaloScale_                   ,&amp;
     &amp;                                   darkMatterProfile_                  =darkMatterProfileAdiabaticGnedin2004_   &amp;
     &amp;                                                        )
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="darkMatterProfileAdiabaticGnedin2004_">
   <constructor>
     darkMatterProfileAdiabaticGnedin2004                         (                                                       &amp;
     &amp;                                   A                                   =0.85d0                                , &amp;
     &amp;                                   omega                               =0.80d0                                , &amp;
     &amp;                                   radiusFractionalPivot               =1.00d0                                , &amp;
     &amp;                                   toleranceRelative                   =1.0d-2                                , &amp;
     &amp;                                   nonAnalyticSolver                   =nonAnalyticSolversNumerical           , &amp;
     &amp;                                   cosmologyParameters_                =cosmologyParameters_                  , &amp;
     &amp;                                   darkMatterHaloScale_                =darkMatterHaloScale_                  , &amp;
     &amp;                                   darkMatterProfileDMO_               =darkMatterProfileIsothermal_          , &amp;
     &amp;                                   galacticStructure_                  =galacticStructureStandard_              &amp;
     &amp;                                                        )
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="satelliteEvaporationSIDM_">
   <constructor>
    satelliteEvaporationSIDMKummer2018      (                                                                             &amp;
     &amp;                                  darkMatterProfileDMO_                = darkMatterProfileIsothermal_         , &amp;
     &amp;                                  galacticStructure_                   = galacticStructureStandard_           , &amp; 
     &amp;                                  darkMatterParticle_                  = darkMatterParticleSIDM_                &amp;
     &amp;                                  )
   </constructor>
  </referenceConstruct>
  <referenceConstruct object="satelliteDecelerationSIDM_">
   <constructor>
    satelliteDecelerationSIDMKummer2018      (                                                                            &amp;
     &amp;                                  darkMatterProfileDMO_                = darkMatterProfileIsothermal_         , &amp;
     &amp;                                  galacticStructure_                   = galacticStructureStandard_           , &amp;
     &amp;                                  darkMatterParticle_                  = darkMatterParticleSIDM_                &amp;
     &amp;                                  )
   </constructor>
  </referenceConstruct>
  !!]
  do i=1,size(velocity)
     evaporationFactor = satelliteEvaporationSIDM_%evaporationFactor(x = x(i),speedOrbital = velocity(i))
     decelerationFactor = satelliteDecelerationSIDM_%decelerationFactor(x = x(i),speedOrbital = velocity(i))
     write (message,'(a,f6.1,a,f6.1,a)') "evaporation factor at [x=",x(i),"] and orbital speed [v=",velocity(i),"]"
     call Assert(trim(message),evaporationFactor,evaporationFactorCheck(i),relTol=1.0d-2)
     write (message,'(a,f6.1,a,f6.1,a)') "deceleration factor at [x=",x(i),"] and orbital speed [v=",velocity(i),"]"
     call Assert(trim(message),decelerationFactor,decelerationFactorCheck(i),relTol=1.0d-2)
  end do
  ! End unit tests.
  call Unit_Tests_End_Group()
  call Unit_Tests_Finish   ()
  ! Clean up objects.
  !![
  <objectDestructor name="cosmologyParameters_"        />
  <objectDestructor name="cosmologyFunctions_"         />
  <objectDestructor name="virialDensityContrast_"      />
  <objectDestructor name="darkMatterHaloScale_"        />
  <objectDestructor name="darkMatterProfileIsothermal_"/>
  <objectDestructor name="galacticStructureStandard_"  />
  <objectDestructor name="satelliteEvaporationSIDM_"   />
  <objectDestructor name="satelliteDecelerationSIDM_"   />
  !!]
end program Tests_SIDM_evaporation_deceleration_factors
