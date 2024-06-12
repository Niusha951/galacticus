module SIDM_Effective_Cross_Section_Integrator
  implicit none
  private
  public :: SIDMEffectiveCrossSectionIntegrator2D

  type :: SIDMEffectiveCrossSectionIntegrator2D
    private
    procedure(integrandTemplate), nopass, pointer :: integrand
    double precision, dimension(2,2) :: boundaries
  contains
    procedure :: integrate => effectiveCrossSection2DIntegrate
    procedure :: setIntegrand
  end type SIDMEffectiveCrossSectionIntegrator2D

  interface SIDMEffectiveCrossSectionIntegrator2D
    module procedure :: effectiveCrossSection2DConstructor
  end interface SIDMEffectiveCrossSectionIntegrator2D

  abstract interface
    double precision function integrandTemplate(x, y)
      double precision, intent(in) :: x, y
    end function integrandTemplate
  end interface

contains

  function effectiveCrossSection2DConstructor() result(self)
    type(SIDMEffectiveCrossSectionIntegrator2D) :: self
    self%integrand => null()
  end function effectiveCrossSection2DConstructor

  subroutine setIntegrand(self, integrand)
    class(SIDMEffectiveCrossSectionIntegrator2D), intent(inout) :: self
    procedure(integrandTemplate) :: integrand

    self%integrand => integrand
  end subroutine setIntegrand

  double precision function effectiveCrossSection2DIntegrate(self, boundaries)
    use :: Numerical_Integration, only: integrator
    use :: Coordinates, only: coordinateCartesian
    implicit none
    class(SIDMEffectiveCrossSectionIntegrator2D), intent(inout), target :: self
    double precision, dimension(2,2), intent(in) :: boundaries
    type(integrator) :: integrator_
    double precision :: x_, y_

    self%boundaries = boundaries

    integrator_ = integrator(effectiveCrossSection2DIntegrandX, toleranceRelative=1.0d-2)
    effectiveCrossSection2DIntegrate = integrator_%integrate(boundaries(1,1), boundaries(1,2))
    return

  contains

    double precision function effectiveCrossSection2DIntegrandX(x)
      implicit none
      double precision, intent(in) :: x
      type(integrator) :: integrator_

      x_ = x
      integrator_ = integrator(effectiveCrossSection2DIntegrandY,toleranceRelative=1.0d-2)
      effectiveCrossSection2DIntegrandX = integrator_%integrate(self%boundaries(2,1), self%boundaries(2,2))
      return
    end function effectiveCrossSection2DIntegrandX

    double precision function effectiveCrossSection2DIntegrandY(y)
      implicit none
      double precision, intent(in) :: y

      y_ = y
      effectiveCrossSection2DIntegrandY = self%integrand(x_, y_)
      return
    end function effectiveCrossSection2DIntegrandY

  end function effectiveCrossSection2DIntegrate

end module SIDM_Effective_Cross_Section_Integrator

