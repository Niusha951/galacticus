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
Contains a module which displays a banner for \glc.
!!}

module Display_Banner
  !!{
  Displays a banner for \glc.
  !!}
  implicit none
  private
  public :: Display_Banner_Show

contains

  subroutine Display_Banner_Show
    !!{
    Displays the \glc\ banner.
    !!}
    use :: Display      , only : displayVerbosity, verbosityLevelSilent
#ifdef USEMPI
    use :: MPI_Utilities, only : mpiSelf
#endif
    implicit none
#ifdef USEMPI
    if (mpiSelf%rank() == 0) then
#endif
       if (displayVerbosity() > verbosityLevelSilent) then
          write (0,*) '             ##                                     '
          write (0,*) '  ####        #                  #                  '
          write (0,*) ' #   #        #             #                       '
          write (0,*) '#       ###   #  ###   ### ###  ##   ### ## ##   ## '
          write (0,*) '#       #  #  #  #  # #  #  #    #  #  #  #  #  #   '
          write (0,*) '#   ###  ###  #   ### #     #    #  #     #  #   #  '
          write (0,*) ' #   #  #  #  #  #  # #     #    #  #     #  #    # '
          write (0,*) '  ####  #### ### ####  ###   ## ###  ###   #### ##  '
          write (0,*)
          write (0,*) '© 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016,'
          write (0,*) '  2017, 2018, 2019, 2020, 2021, 2022, 2023, 2024'
          write (0,*) '  - Andrew Benson'
          write (0,*)
       end if
#ifdef USEMPI
    end if
#endif
    return
  end subroutine Display_Banner_Show

end module Display_Banner
