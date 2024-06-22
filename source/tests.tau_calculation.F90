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

program Tests_Tau_Calculation
  use :: Display            , only : displayVerbositySet     , verbosityLevelStandard
  use :: Functions_Global_Utilities, only : Functions_Global_Set
  use :: Cosmology_Functions, only : cosmologyFunctionsMatterLambda
  use :: Cosmology_Parameters,only : cosmologyParametersSimple
  use :: Dark_Matter_Halo_Scales, only : darkMatterHaloScaleVirialDensityContrastDefinition
  use :: Dark_Matter_Particles, only : darkMatterParticleSelfInteractingDarkMatter, darkMatterParticleSIDMVelocityDependent
  !use :: Virial_Density_Contrast   , only : virialDensityContrastSphericalCollapseClsnlssMttrCsmlgclCnstnt
  use :: Dark_Matter_Profiles_DMO, only : darkMatterProfileDMOPenarrubia2010, darkMatterProfileDMONFW
  use :: Events_Hooks                , only : eventsHooksInitialize
  use :: Galacticus_Nodes   , only : mergerTree              , treeNode              , treeNodeList, nodeComponentBasic, nodeComponentDarkmatterProfile, nodeComponentBasicStandard, nodeClassHierarchyInitialize
  use :: Node_Components           , only : Node_Components_Initialize, Node_Components_Uninitialize, Node_Components_Thread_Uninitialize
  use :: ISO_Varying_String , only : assignment(=)           , varying_string
  use :: Input_Parameters   , only : inputParameters
  use :: Kind_Numbers       , only : kind_int8
  !use :: Merger_Tree_Walkers, only : mergerTreeWalkerAllNodes
  use :: Unit_Tests         , only : Assert                  , Unit_Tests_Begin_Group, Unit_Tests_End_Group, Unit_Tests_Finish
  use :: Virial_Density_Contrast, only : virialDensityContrastSphericalCollapseClsnlssMttrCsmlgclCnstnt, virialDensityContrastFixed, fixedDensityTypeCritical, virialDensityContrastBryanNorman1998
  use :: tauCalculationClassModule, only : tauCalculation, tauCalculationClass


  implicit none
  type   (varying_string          )            :: parameterFile
  integer(kind=kind_int8          ), parameter :: N = 214
  type   (treeNodeList            )            :: nodes        (N)
  !logical                                      :: nodeFound    (N)
  type   (mergerTree              )            :: tree
  type   (treeNode                ), pointer   :: node
  type   (cosmologyParametersSimple)     , pointer :: cosmologyParameters_
  type   (cosmologyFunctionsMatterLambda), pointer :: cosmologyFunctions_
  type   (darkMatterProfileDMONFW), pointer    :: darkMatterProfileDMONFW_ 
  type   (darkMatterHaloScaleVirialDensityContrastDefinition), pointer :: darkMatterHaloScale_ 
  !type   (virialDensityContrastSphericalCollapseClsnlssMttrCsmlgclCnstnt), pointer :: virialDensityContrast_
  !type   (virialDensityContrastFixed), pointer :: virialDensityContrast_
  type   (virialDensityContrastBryanNorman1998), pointer :: virialDensityContrast_
  type   (darkMatterParticleSelfInteractingDarkMatter) :: darkMatterParticle_
  type   (darkMatterParticleSIDMVelocityDependent) :: darkMatterParticleSIDMVelDep_

  class(nodeComponentBasic),             pointer :: basic
  class(nodeComponentDarkmatterProfile), pointer :: darkMatterProfile
  !class(darkMatterProfileDMONFW),      pointer :: darkMatterProfile

  integer(kind=kind_int8          )            :: i, io_status
  type   (inputParameters         )          :: parameters
  !type   (mergerTreeWalkerAllNodes)            :: treeWalker
  !real(kind=8)                                 :: time_test(N), scale_a(N), mvir_test(N), rvir_test(N), Vmax_test(N), Rvmax_test(N)
  double precision                             :: time_test(N), scale_a(N), mvir_test(N), rvir_test(N), Vmax_test(N), Rvmax_test(N), rs_test(N), velocityMaximum(N), radiusVelocityMaximum(N), radiusScale(N), massVirial(N), radiusVirial(N)
 
  character(len=256)                           :: line
  double precision :: temp1, temp2, temp3, temp4, temp5, temp6, temp7, temp8
  double precision :: h = 0.7

  integer :: output_unit
  type(tauCalculation) :: tauCalc


  ! Set verbosity level.
  call displayVerbositySet(verbosityLevelStandard)
  ! Begin unit tests.
  call Unit_Tests_Begin_Group("SIDM parametric model tests: tau calculation")


  ! Open the parameter file.
  parameterFile='testSuite/parameters/nodes/nodes_modified.xml'
  parameters=inputParameters(parameterFile)

  ! Initialize the Galacticus nodes objects module.
  !call eventsHooksInitialize       (          )
  !call Functions_Global_Set        (          )
  call nodeClassHierarchyInitialize(parameters)
  !call Node_Components_Initialize  (parameters)
  !call Node_Components_Thread_Initialize(parameters)

  allocate(cosmologyParameters_                             )
  allocate(cosmologyFunctions_                              )
  allocate(virialDensityContrast_                           )
  allocate(darkMatterHaloScale_                             )
  allocate(darkMatterProfileDMONFW_                         )
  allocate(darkMatterParticle_                              )
  allocate(darkMatterParticleSIDMVelDep_                    )

  cosmologyParameters_ = cosmologyParametersSimple( &
       OmegaMatter=0.2815d0, &
       OmegaBaryon=0.0465d0, &
       OmegaDarkEnergy=0.7185d0, &
       temperatureCMB=2.7800d0, &
       HubbleConstant=70.0d0)
       !HubbleConstant=69.3000d0)
  cosmologyFunctions_ = cosmologyFunctionsMatterLambda(cosmologyParameters_)
  !virialDensityContrast_ = virialDensityContrastSphericalCollapseClsnlssMttrCsmlgclCnstnt(tableStore=.false., cosmologyFunctions_ = cosmologyFunctions_)
  !virialDensityContrast_ = virialDensityContrastFixed(densityContrastValue = 200.0d0, densityType = fixedDensityTypeCritical, turnAroundOverVirialRadius = 2.0d0, cosmologyParameters_ = cosmologyParameters_, cosmologyFunctions_ = cosmologyFunctions_)
  virialDensityContrast_ = virialDensityContrastBryanNorman1998(allowUnsupportedCosmology = .false., cosmologyParameters_ = cosmologyParameters_, cosmologyFunctions_ = cosmologyFunctions_)
  darkMatterHaloScale_ = darkMatterHaloScaleVirialDensityContrastDefinition(cosmologyParameters_ = cosmologyParameters_, cosmologyFunctions_ = cosmologyFunctions_, virialDensityContrast_ = virialDensityContrast_)
  darkMatterProfileDMONFW_ = darkMatterProfileDMONFW(velocityDispersionUseSeriesExpansion=.false., darkMatterHaloScale_ = darkMatterHaloScale_)
  darkMatterParticle_ = darkMatterParticleSelfInteractingDarkMatter()
  darkMatterParticleSIDMVelDep_ = darkMatterParticleSIDMVelocityDependent(velocityCharacteristic = 24.3289794155754d0, sigma0 = 147.10088d0, darkMatterParticle_ = darkMatterParticle_)

  ! Read the data from file
  open(unit=10, file='data_799_cdm_NFW.txt', status='old', action='read')
  read(10, '(A)', iostat=io_status) ! Skip the header line
  if (io_status /= 0) then
    print *, 'Error reading file'
    stop
  end if
  i = N
  do 
    read(10, '(A)', iostat=io_status) line
    if (io_status /= 0) exit
    read(line, *) temp1, temp2, temp3, temp4
    scale_a(i) = temp1
    mvir_test(i) = temp2/h
    rvir_test(i) = temp3*scale_a(i)/h
    Rvmax_test(i) = temp4*scale_a(i)/h

    time_test(i) = cosmologyFunctions_%cosmicTime(scale_a(i))
    rs_test(i) = Rvmax_test(i)/2.1626d0

    i = i - 1
  end do
  close(10)

  ! Read the data from file
  open(unit=10, file='data_799_cdm.txt', status='old', action='read')
  read(10, '(A)', iostat=io_status) ! Skip the header line
  if (io_status /= 0) then
    print *, 'Error reading file'
    stop
  end if
  i = N
  do
    read(10, '(A)', iostat=io_status) line
    if (io_status /= 0) exit
    read(line, *) temp1, temp2, temp3, temp4, temp5, temp6, temp7, temp8
    !scale_a(i) = temp1
    !mvir_test(i) = temp2/h
    !rvir_test(i) = temp5*scale_a(i)/h
    !rs_test(i) = temp6*scale_a(i)/h
    Vmax_test(i) = temp7
    !Rvmax_test(i) = temp8*scale_a(i)/h

    !time_test(i) = cosmologyFunctions_%cosmicTime(scale_a(i))
    !rs_test(i) = Rvmax_test(i)/2.1626d0

    i = i - 1
  end do
  close(10)

  !print *, 'a_scale from file'
  do i=1, 214
  !   print *, scale_a(i)
  !   print *, time_test(i)
  !   print *, Rvmax_test(i)
  end do  

  ! Create nodes.
  do i=1,N
     nodes(i)%node => treeNode()
  end do

  ! Set child nodes.
  do i=1,N-1
     nodes(i)%node%firstChild => nodes(i+1)%node
     nodes(i+1)%node%parent => nodes(i)%node
  end do

  do i=1,N
     
     basic => nodes(i)%node%basic(autoCreate=.true.)
     call basic%massSet(mvir_test(i))
     call basic%timeSet(time_test(i)) 

     darkMatterProfile => nodes(i)%node%darkMatterProfile(autoCreate=.true.)
     call darkMatterProfile%scaleSet(rs_test(i)/1e3)

     !print *, scale_a(i)    
 
     velocityMaximum(i)=darkMatterProfileDMONFW_%circularVelocityMaximum(nodes(i)%node)
     !call Assert("maximum Velocity comparison:", velocityMaximum(i),Vmax_test(i), relTol=1.0d-3)

     radiusVelocityMaximum(i)=darkMatterProfileDMONFW_%radiusCircularVelocityMaximum(nodes(i)%node)
     radiusScale(i)=darkMatterProfile%scale()
     massVirial(i)= basic%mass()

     radiusVirial(i)=darkMatterHaloScale_%radiusVirial(nodes(i)%node)
     !call Assert("maximum Velocity radius comparison:", radiusVelocityMaximum(i),Rvmax_test(i), relTol=1.0d-3)

     ! Call the tauCalculationClass subroutine
     call tauCalculationClass(tauCalc, nodes(i)%node)

     ! Retrieve the computed values of VmaxSIDM and RmaxSIDM
     call darkMatterProfile%floatRank0MetaPropertyGet(tauCalc%VmaxSIDMID, VmaxSIDM)
     call darkMatterProfile%floatRank0MetaPropertyGet(tauCalc%RmaxSIDMID, RmaxSIDM)


     call Assert("virial radius: ", radiusVirial(i)*1e3, rvir_test(i), relTol=1.0d-2)
     call Assert("mass radius: ", massVirial(i), mvir_test(i), relTol=1.0d-2)
     call Assert("max Velocity: ", velocityMaximum(i), Vmax_test(i), relTol=1.0d-2)
     !call Assert("scale radius: ", radiusScale(i), rs_test(i), relTol=1.0d-3)

  end do

  ! Open the output file
!  open(unit=20, file='output_data.txt', status='replace', action='write')
!  write(20, '(A)') 'time_test rs_test Rs Rvmax_test Rvmax Vmax_test Vmax Rvir_test Rvir Mvir_test Mvir' ! Write header line

!  do i=1,N
     ! Write the data to the output file
!     write(20, '(F20.6, 2X, F20.6, 2X, F20.6, 2X, F20.6, 2X, F20.6, 2X, F20.6, 2X, F20.6, 2X, F20.6, 2X, F20.6, 2X, F20.6, 2X, F20.6)') time_test(i), rs_test(i), radiusScale(i), Rvmax_test(i), radiusVelocityMaximum(i), Vmax_test(i), velocityMaximum(i), rvir_test(i), radiusVirial(i), mvir_test(i), massVirial(i)  
!  end do

 ! Close the output file
!  close(20)

  ! Set the base node of our tree.
  tree%nodeBase => nodes(1)%node

  !call Assert('All nodes walked to',all(nodeFound),.true.)

  ! Destroy nodes.
  do i=1,N
     call nodes(i)%node%destroy()
  end do
  

  ! End unit tests.
  call Unit_Tests_End_Group()
  call Unit_Tests_Finish()
  !call nodeClassHierarchyFinalize         ()

end program Tests_Tau_Calculation


