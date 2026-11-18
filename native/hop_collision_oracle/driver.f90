! SPDX-License-Identifier: GPL-3.0-only
! Read-only driver of pinned OpenMULTI3D HOP interpolation; no RT/EOS/SE.
program collision_interpolation_probe
  use interpol_mod, only: hermite_derivative, hermite_interpolation
  use const_mod, only: e_k
  implicit none
  integer :: nr,ng,nq,i,j
  real(8),allocatable :: grid(:),query(:),values(:),slopes(:),answer(:)
  read(*,*)nr,ng,nq
  if(nr<1.or.nr>2000.or.ng<4.or.ng>100.or.nq<1.or.nq>256) stop 1
  allocate(grid(ng),query(nq),values(ng),slopes(ng),answer(nq))
  read(*,*)grid
  read(*,*)query
  if(any(grid(2:ng)<=grid(1:ng-1)).or.minval(query)<grid(1).or.maxval(query)>grid(ng)) stop 2
  write(*,'(a,es26.17e3)')'E_K ',e_k
  do i=1,nr
    read(*,*)values
    if(any(values<=0))stop 3
    values=log(values)
    call hermite_derivative(slopes,grid,values)
    call hermite_interpolation(answer,query,grid,values,slopes)
    answer=exp(answer)
    do j=1,nq
      write(*,'(2i6,es26.17e3)')i,j,answer(j)
    enddo
  enddo
end program
