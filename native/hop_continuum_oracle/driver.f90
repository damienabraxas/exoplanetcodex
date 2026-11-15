! SPDX-License-Identifier: GPL-3.0-only
! Read-only caller of pinned OpenMULTI3D HOP continuum interpolation.
program continuum_probe
  use trans_mod, only: cont_type
  use const_mod, only: HZ_TO_AA
  implicit none
  type(cont_type) :: cont
  real(8), allocatable :: wave(:), nu(:)
  real(4), allocatable :: alpha(:)
  integer :: n,m,k
  read(*,*) n,m
  if(n<2.or.n>20000.or.m<2.or.m>20000) stop 1
  allocate(cont%lam_input(n),cont%alpha_input(n),wave(m),nu(m),alpha(m))
  do k=1,n
    read(*,*) cont%lam_input(k),cont%alpha_input(k)
  enddo
  do k=1,m
    read(*,*) wave(k)
  enddo
  if(any(cont%lam_input<=0d0).or.any(wave<=0d0))stop 2
  if(any(cont%lam_input(2:)<=cont%lam_input(:n-1)))stop 3
  nu=HZ_TO_AA/wave
  call cont%hermite_interp(alpha,nu)
  do k=1,m
    write(*,'(es26.17e3)')real(alpha(k),8)
  enddo
end program continuum_probe
