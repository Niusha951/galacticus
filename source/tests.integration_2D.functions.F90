module Test_Integration2D_Functions
  implicit none
  contains

  double precision function Integrand1(x, y)
    double precision, intent(in) :: x, y
    Integrand1 = sin(x**2) * cos(y)
  end function Integrand1

  double precision function Integrand2(x, y)
    double precision, intent(in) :: x, y
    Integrand2 = x**2 * cos(y)
  end function Integrand2

  double precision function Integrand3(x, y)
    double precision, intent(in) :: x, y
    Integrand3 = sqrt(x) * y**2
  end function Integrand3

  double precision function Integrand4(x, y)
    double precision, intent(in) :: x, y
    Integrand4 = y / sqrt(x)
  end function Integrand4

end module Test_Integration2D_Functions

