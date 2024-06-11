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
Contains a module of integrands for unit tests.
!!}

module Test_Integration2D_Functions
  !!{
  Contains integrands for unit tests.
  !!}
  implicit none
  private
  public :: Integrand1, Integrand2, Integrand3, Integrand4

contains

  double precision function Integrand1(x,y)
    !!{
    Integral for unit testing.
    !!}
    implicit none
    double precision, intent(in   ) :: x, y

    Integrand1=sin(x**2)*cos(y)
    return
  end function Integrand1

  double precision function Integrand2(x,y)
    !!{
    Integral for unit testing.
    !!}
    implicit none
    double precision, intent(in   ) :: x, y

    Integrand2=(x**2) * cos(y) 
    return
  end function Integrand2

  double precision function Integrand3(x,y)
    !!{
    Integral for unit testing.
    !!}
    implicit none
    double precision, intent(in   ) :: x, y

    Integrand3=sqrt(x) * y**2
    return
  end function Integrand3

  double precision function Integrand4(x,y)
    !!{
    Integral for unit testing.
    !!}
    implicit none
    double precision, intent(in   ) :: x, y

    Integrand4 = 1.0d0/sqrt(x) * y
    return
  end function Integrand4

end module Test_Integration2D_Functions
