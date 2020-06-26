!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019
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

  use :: Atomic_Ionization_Potentials, only : atomicIonizationPotentialClass
  
  !# <radiativeTransferOutputter name="radiativeTransferOutputterContinuuaRates">
  !#  <description>A radiative transfer outputter class which outputs the photon emission rates in the continuua of all specified elements.</description>
  !# </radiativeTransferOutputter>
  type, extends(radiativeTransferOutputterClass) :: radiativeTransferOutputterContinuuaRates
     !% Implementation of a radiative transfer outputter class which outputs the Lyman continuum photon emission rate.
     private
     class           (atomicIonizationPotentialClass), pointer                     :: atomicIonizationPotential_ => null()
     integer                                         , allocatable, dimension(:  ) :: elementIndices
     double precision                                , allocatable, dimension(:,:) :: continuuaRatesEscaping              , continuumLimitWavelength
     integer                                                                       :: countElements                       , atomicNumberMaximum
   contains
     procedure :: reset               => continuuaRatesReset
     procedure :: sourceProperties    => continuuaRatesSourceProperties
     procedure :: photonPacketEscapes => continuuaRatesPhotonPacketEscapes
     procedure :: finalize            => continuuaRatesFinalize
     procedure :: output              => continuuaRatesOutput
  end type radiativeTransferOutputterContinuuaRates

  interface radiativeTransferOutputterContinuuaRates
     !% Constructors for the {\normalfont \ttfamily continuuaRates} radiative transfer outputter packet class.
     module procedure continuuaRatesConstructorParameters
     module procedure continuuaRatesConstructorInternal
  end interface radiativeTransferOutputterContinuuaRates
  
contains

  function continuuaRatesConstructorParameters(parameters) result(self)
    !% Constructor for the {\normalfont \ttfamily continuuaRates} radiative transfer outputter class which takes a parameter set as
    !% input.
    use :: Atomic_Data     , only : Atom_Lookup
    use :: Input_Parameters, only : inputParameters
    implicit none
    type     (radiativeTransferOutputterContinuuaRates)                              :: self
    type     (inputParameters                         ), intent(inout)               :: parameters
    class    (atomicIonizationPotentialClass          ), pointer                     :: atomicIonizationPotential_
    character(len=3                                   ), allocatable  , dimension(:) :: elements
    integer                                            , allocatable  , dimension(:) :: elementIndices
    integer                                                                          :: i                          , elementsCount

    elementsCount=parameters%count('elements',zeroIfNotPresent=.false.)
    allocate(elements      (elementsCount))
    allocate(elementIndices(elementsCount))
    !# <inputParameter>
    !#   <name>elements</name>
    !#   <cardinality>1..*</cardinality>
    !#   <description>The names of the elements to be tracked.</description>
    !#   <source>parameters</source>
    !#   <type>string</type>
    !# </inputParameter>
    do i=1,elementsCount
       elementIndices(i)=Atom_Lookup(shortLabel=elements(i))
    end do
    !# <objectBuilder class="atomicIonizationPotential" name="atomicIonizationPotential_" source="parameters"/>
    self=radiativeTransferOutputterContinuuaRates(elementIndices,atomicIonizationPotential_)
    !# <objectDestructor name="atomicIonizationPotential_"/>
    !# <inputParametersValidate source="parameters"/>
    return
  end function continuuaRatesConstructorParameters

  function continuuaRatesConstructorInternal(elementIndices,atomicIonizationPotential_) result(self)
    !% Internal constructor for the {\normalfont \ttfamily continuuaRates} radiative transfer outputter class.
    use :: Atomic_Data                 , only : Atomic_Number
    use :: Numerical_Constants_Physical, only : speedLight   , plancksConstant
    use :: Numerical_Constants_Units   , only : electronVolt , angstromsPerMeter
    implicit none
    type   (radiativeTransferOutputterContinuuaRates)                              :: self
    integer                                          , intent(in   ), dimension(:) :: elementIndices
    class  (atomicIonizationPotentialClass          ), intent(in   ), target       :: atomicIonizationPotential_
    integer                                                                        :: i                         , j
    !# <constructorAssign variables="elementIndices, *atomicIonizationPotential_"/>

    ! Determine maximum atomic number.
    self%countElements=size(elementIndices)
    self%atomicNumberMaximum=0
    do i=1,self%countElements
       self%atomicNumberMaximum=max(self%atomicNumberMaximum,Atomic_Number(elementIndices(i)))
    end do
    ! Allocate arrays for escape rates and wavelengths.
    allocate(self%continuuaRatesEscaping  (self%countElements,self%atomicNumberMaximum))
    allocate(self%continuumLimitWavelength(self%countElements,self%atomicNumberMaximum))
    ! Initialize escaping rates to zero.
    self%continuuaRatesEscaping=0.0d0
    ! Compute continuum limit wavelengths for each element.
    do i=1,self%countElements
       do j=1,Atomic_Number(elementIndices(i))
          self%continuumLimitWavelength(i,j)=+plancksConstant                                                                                &
               &                             *speedLight                                                                                     &
               &                             /electronVolt                                                                                   &
               &                             *angstromsPerMeter                                                                              &
               &                             /self%atomicIonizationPotential_%potential(                                                     &
               &                                                                        atomicNumber  =Atomic_Number(elementIndices(i))    , &
               &                                                                        electronNumber=Atomic_Number(elementIndices(i))+1-j  &
               &                             )
       end do
    end do
    return
  end function continuuaRatesConstructorInternal

  subroutine continuuaRatesDestructor(self)
    !% Destructor for the {\normalfont \ttfamily continuuaRates} radiative transfer outputter class.
    implicit none
    type(radiativeTransferOutputterContinuuaRates), intent(inout) :: self

    !# <objectDestructor name="self%atomicIonizationPotential_"/>
    return
  end subroutine continuuaRatesDestructor

  subroutine continuuaRatesReset(self)
    !% Reset the accumulated Lyman continuum photon escape rate.
    implicit none
    class(radiativeTransferOutputterContinuuaRates), intent(inout) :: self

    self%continuuaRatesEscaping=0.0d0
    return
  end subroutine continuuaRatesReset

  subroutine continuuaRatesSourceProperties(self,radiativeTransferSource_,outputGroup)
    !% Compute and output the Lyman continuum photon emission rate.
    use :: Atomic_Data             , only : Atomic_Number, Atomic_Short_Label
    use :: IO_HDF5                 , only : hdf5Access
    use :: Numerical_Integration   , only : integrator
    use :: Numerical_Roman_Numerals, only : Roman_Numerals
    implicit none
    class           (radiativeTransferOutputterContinuuaRates), intent(inout)                 :: self
    class           (radiativeTransferSourceClass            ), intent(inout)                 :: radiativeTransferSource_
    type            (hdf5Object                              ), intent(inout)                 :: outputGroup
    double precision                                          , allocatable  , dimension(:,:) :: rateEmitted
    type            (integrator                              )                                :: integrator_
    integer                                                                                   :: i                       , j
    character       (len=30                                  )                                :: label

    integrator_=integrator(                             &
         &                                   integrand, &
         &                 toleranceRelative=1.0d-2     &
         &                )
    allocate(rateEmitted(self%countElements,self%atomicNumberMaximum))
    do i=1,size(self%elementIndices)
       do j=1,Atomic_Number(self%elementIndices(i))
          rateEmitted(i,j)=+integrator_%integrate(                                           &
               &                                  1.0d-6*self%continuumLimitWavelength(i,j), &
               &                                         self%continuumLimitWavelength(i,j)  &
               &                                 )
       end do
    end do
    !$ call hdf5Access%set  ()
    do i=1,size(self%elementIndices)
       do j=1,Atomic_Number(self%elementIndices(i))
          label='rateContinuumEmitted'//trim(Atomic_Short_Label(self%elementIndices(i)))//char(Roman_Numerals(j))
          call outputGroup%writeAttribute(rateEmitted(i,j),trim(label))
       end do
    end do
    !$ call hdf5Access%unset()
    return

  contains

    double precision function integrand(wavelength)
      !% Integrand over the source spectrum.
      use :: Numerical_Constants_Physical    , only : plancksConstant  , speedLight
      use :: Numerical_Constants_Units       , only : angstromsPerMeter
      use :: Numerical_Constants_Astronomical, only : luminositySolar
      implicit none
      double precision, intent(in   ) :: wavelength
      double precision                :: energyPhoton
      
      energyPhoton=+plancksConstant                               &
           &       *speedLight                                    &
           &       *angstromsPerMeter                             &
           &       /wavelength
      integrand   =+radiativeTransferSource_%spectrum(wavelength) &
           &       *luminositySolar                               &
           &       /energyPhoton
      return
    end function integrand

  end subroutine continuuaRatesSourceProperties

  subroutine continuuaRatesPhotonPacketEscapes(self,photonPacket)
    !% Process an escaping photon packet.
    use :: Atomic_Data                     , only : Atomic_Number, Atomic_Short_Label
    use :: Numerical_Constants_Physical    , only : plancksConstant                   , speedLight
    use :: Numerical_Constants_Units       , only : angstromsPerMeter
    use :: Numerical_Constants_Astronomical, only : luminositySolar
    implicit none
    class           (radiativeTransferOutputterContinuuaRates), intent(inout) :: self
    class           (radiativeTransferPhotonPacketClass      ), intent(inout) :: photonPacket
    double precision                                                          :: energyPhoton
    integer                                                                   :: i           , j
    
    do i=1,size(self%elementIndices)
       do j=1,Atomic_Number(self%elementIndices(i))
          if (photonPacket%wavelength() < self%continuumLimitWavelength(i,j)) then
             energyPhoton                    =+plancksConstant                               &
                  &                           *speedLight                                    &
                  &                           *angstromsPerMeter                             &
                  &                           /photonPacket     %wavelength            (   )
             self%continuuaRatesEscaping(i,j)=+self             %continuuaRatesEscaping(i,j) &
                  &                           +photonPacket     %luminosity            (   ) &
                  &                           *luminositySolar                               &
                  &                           /energyPhoton
          end if
       end do
    end do
    return
  end subroutine continuuaRatesPhotonPacketEscapes

  subroutine continuuaRatesFinalize(self)
    !% Finalize the Lyman continuum photon escape rate.
    use :: MPI_Utilities, only : mpiSelf
    implicit none
    class(radiativeTransferOutputterContinuuaRates), intent(inout) :: self

    ! Sum the Lyc rate across all MPI processes.
    self%continuuaRatesEscaping=mpiSelf%sum(self%continuuaRatesEscaping)
    return
  end subroutine continuuaRatesFinalize

  subroutine continuuaRatesOutput(self,outputGroup)
    !% Output the Lyman continuum photon escape rate.
    use :: Atomic_Data             , only : Atomic_Number , Atomic_Short_Label
    use :: IO_HDF5                 , only : hdf5Access
    use :: Numerical_Roman_Numerals, only : Roman_Numerals
    implicit none
    class    (radiativeTransferOutputterContinuuaRates), intent(inout) :: self
    type     (hdf5Object                              ), intent(inout) :: outputGroup
    integer                                                            :: i          , j
    character(len=30                                  )                :: label

    !$ call hdf5Access%set  ()
    do i=1,size(self%elementIndices)
       do j=1,Atomic_Number(self%elementIndices(i))
          label='rateContinuumEscaping'//trim(Atomic_Short_Label(self%elementIndices(i)))//char(Roman_Numerals(j))
          call outputGroup%writeAttribute(self%continuuaRatesEscaping(i,j),trim(label))
       end do
    end do
    !$ call hdf5Access%unset()
    return
  end subroutine continuuaRatesOutput