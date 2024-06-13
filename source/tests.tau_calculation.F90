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
  !use :: Functions_Global_Utilities, only : Functions_Global_Set
  use :: Galacticus_Nodes   , only : mergerTree              , treeNode              , treeNodeList, nodeComponentBasic, nodeComponentBasicStandard, nodeClassHierarchyInitialize
  !use :: Node_Components           , only : Node_Components_Initialize, Node_Components_Uninitialize
  use :: ISO_Varying_String , only : assignment(=)           , varying_string
  use :: Input_Parameters   , only : inputParameters
  use :: Kind_Numbers       , only : kind_int8
  !use :: Merger_Tree_Walkers, only : mergerTreeWalkerAllNodes
  use :: Unit_Tests         , only : Assert                  , Unit_Tests_Begin_Group, Unit_Tests_End_Group, Unit_Tests_Finish
  use :: Cosmology_Parameters,only : cosmologyParametersSimple
  use :: Cosmology_Functions, only : cosmologyFunctionsMatterLambda
  implicit none
  type   (varying_string          )            :: parameterFile
  integer(kind=kind_int8          ), parameter :: N = 214
  type   (treeNodeList            )            :: nodes        (N)
  !logical                                      :: nodeFound    (N)
  type   (mergerTree              )            :: tree
  type   (treeNode                ), pointer   :: node
  type   (cosmologyParametersSimple)     , pointer :: cosmologyParameters_
  type   (cosmologyFunctionsMatterLambda), pointer :: cosmologyFunctions_

  class(nodeComponentBasic),            pointer :: basic

  integer(kind=kind_int8          )            :: i, io_status
  type   (inputParameters         )          :: parameters
  !type   (mergerTreeWalkerAllNodes)            :: treeWalker
  !real(kind=8)                                 :: time_test(N), scale_a(N), mvir_test(N), rvir_test(N), Vmax_test(N), Rvmax_test(N)
  double precision                             :: time_test(N), scale_a(N), mvir_test(N), rvir_test(N), Vmax_test(N), Rvmax_test(N)
 
  character(len=256)                           :: line
  double precision :: temp1, temp2, temp3, temp4, temp5, temp6, temp7, temp8


  ! Set verbosity level.
  call displayVerbositySet(verbosityLevelStandard)
  ! Begin unit tests.
  call Unit_Tests_Begin_Group("SIDM parametric model tests: tau calculation")

  ! Open the parameter file.
  parameterFile='testSuite/parameters/nodes/nodes_modified.xml'
  parameters=inputParameters(parameterFile)

  ! Initialize the Galacticus nodes objects module.
  !call Functions_Global_Set        (          )
  call nodeClassHierarchyInitialize(parameters)
  !call Node_Components_Initialize  (parameters)


  allocate(cosmologyParameters_                             )
  allocate(cosmologyFunctions_                              )

  cosmologyParameters_ = cosmologyParametersSimple( &
       OmegaMatter=0.2815d0, &
       OmegaBaryon=0.0465d0, &
       OmegaDarkEnergy=0.7185d0, &
       temperatureCMB=2.7800d0, &
       HubbleConstant=69.3000d0)
  cosmologyFunctions_ = cosmologyFunctionsMatterLambda(cosmologyParameters_)



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
    scale_a(i) = temp1
    mvir_test(i) = temp2
    rvir_test(i) = temp5
    Vmax_test(i) = temp7
    Rvmax_test(i) = temp8

    time_test(i) = cosmologyFunctions_%cosmicTime(scale_a(i))

    i = i - 1
  end do
  close(10)

  print *, 'a_scale from file'
  do i=1, 214
     print *, scale_a(i)
     print *, time_test(i)
  end do
  

  ! Create nodes.
  do i=1,N
     nodes(i)%node => treeNode()
  end do

  ! Set child nodes.
  do i=1,N-1
     nodes(i)%node%firstChild => nodes(i+1)%node
     nodes(i+1)%node%parent => nodes(i)%node
     
     basic => nodes(i)%node%basic(autoCreate=.true.)
     call basic%massSet(mvir_test(i))
     call basic%timeSet(time_test(i)) 
  end do

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

end program Tests_Tau_Calculation


