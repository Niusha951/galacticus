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

module SIDM_Effective_Cross_Section_Integrator
  !!{
  Implements numerical 2D integration.
  !!}
  implicit none
  private
  public :: SIDMEffectiveCrossSectionIntegrator2D

  type :: SIDMEffectiveCrossSectionIntegrator2D
     !!{
     Implementation of a computational domain for 2D Cartesian cells.
     !!}
     private
     procedure        (integrandTemplate), nopass        , pointer:: integrand
     double precision                    , dimension(2,2)         :: boundaries
   contains
     procedure :: integrate => effectiveCrossSection2DIntegrate
  end type SIDMEffectiveCrossSectionIntegrator2D

  interface SIDMEffectiveCrossSectionIntegrator2D
     !!{
     Constructors for the {\normalfont \ttfamily effetiveCrossSection2D} computational domain.
     !!}
     module procedure :: effectiveCrossSection2DConstructor
  end interface SIDMEffectiveCrossSectionIntegrator2D

  ! Integrand interface.
  abstract interface
     double precision function integrandTemplate(x,y)
       double precision, intent(in   ) :: x,y
     end function integrandTemplate
  end interface

contains

  function effetiveCrossSection2DConstructor(integrand) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily effetiveCrossSection2D} computational domain volume integrator class.
    !!}
    implicit none
    type            (SIDMEffectiveCrossSectionIntegrator2D)                                                :: self
    procedure       (integrandTemplate)                                                             :: integrand

    self%integrand         => integrand

    return
  end function effectiveCrossSection2DConstructor

  double precision function effectiveCrossSection2DIntegrate(self, boundaries)
    !!{
    Integrate over the computational domain cell.
    !!}
    use :: Numerical_Integration, only : integrator
    use :: Coordinates          , only : coordinateCartesian
    implicit none
    class    (SIDMEffectiveCrossSectionIntegrator2D)                , intent(inout), target :: self
    double precision                                                , dimension(2,2), intent(in   ) :: boundaries
!    procedure(integrandDenominator                          )                        :: integrand
    type     (integrator                                    )                        :: integrator_
    double precision                                                                 :: x_, y_

    integrator_         = integrator(                                         &
         &                           effectiveCrossSection2DIntegrandX,       &
         &                           toleranceRelative=1.0d-2                 &
         &                          )
    effectiveCrossSection2DIntegrate=+integrator_%integrate(                  &
         &                                                  boundaries(1,1),  &
         &                                                  boundaries(1,2)   &
         &                                                 )  
    return

  contains

    double precision function effectiveCrossSection2DIntegrandX(x)
      !!{
      $x$-integrand over Cartesian 2D computational domain cells.
      !!}
      implicit none
      double precision            , intent(in   ) :: x
      type            (integrator)                :: integrator_

      x_=x
      integrator_         = integrator            (                                         &
           &                                       effectiveCrossSection2DIntegrandY,       &
           &                                       toleranceRelative=1.0d-2                 &
           &                                      )
      effectiveCrossSection2DIntegrandX=+integrator_%integrate(                             &
           &                                                   boundaries(2,1) ,            &
           &                                                   boundaries(2,2)              &
           &                                                  )  
      return
    end function effectiveCrossSection2DIntegrandX

    double precision function effectiveCrossSection2DIntegrandY(y)
      !!{
      $y$-integrand over Cartesian 2D computational domain cells.
      !!}
      implicit none
      double precision            , intent(in   ) :: y

      y_=y
      effectiveCrossSection2DIntegrandY=self%integrand(x_,y_)
      return
    end function effectiveCrossSection2DIntegrandY

  end function effectiveCrossSection2DIntegrate

end module SIDM_Effective_Cross_Section_Integrator
