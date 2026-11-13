! SPDX-License-Identifier: GPL-3.0-only
! Diagnostic caller of pinned OpenMULTI3D HOP trans_mod::get_lam.
! Upstream source/license remains unchanged; no radiation or population update.
program line_grid_probe
  use trans_mod, only: line_type
  implicit none
  type(line_type) :: line
  real(8), pointer, contiguous :: lam(:)
  integer :: nr, k, j, nq
  real(8) :: qnorm, qmax, q0
  logical :: legacy
  read(*,*) nr, qnorm, legacy
  if(nr/=256.or.qnorm/=10d0.or.legacy) stop 1
  do k=1,nr
    read(*,*) nq,qmax,q0
    if(nq<5.or.nq>1000.or.mod(nq,2)/=1.or.qmax<=0d0.or.q0<0d0)stop 2
    line%nnu=nq;line%qmax=qmax;line%q0=q0;line%lam0=1d0
    call line%get_lam(qnorm,lam,legacy)
    write(*,'(2(i0,1x))') k,size(lam)
    do j=1,size(lam)
      write(*,'(es26.17e3)')lam(j)-1d0
    enddo
    deallocate(lam)
  enddo
end program line_grid_probe
