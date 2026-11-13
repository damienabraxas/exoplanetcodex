! SPDX-License-Identifier: GPL-3.0-only
! Typed calls to exact pinned native wavelength helpers; no copied formula.
program airvac_probe
  use handyfuncs_mod, only: air2vac,vac2air
  implicit none
  real(8) :: vac,table_coordinate
  integer :: ios
  do
    read(*,*,iostat=ios) vac,table_coordinate
    if(ios<0)exit
    if(ios/=0.or.vac<=0.or.table_coordinate<=0)stop 2
    write(*,'(5(es25.16e3,1x))') vac,table_coordinate,air2vac(table_coordinate), &
      vac2air(vac),air2vac(vac2air(vac))
  end do
end program
