! SPDX-License-Identifier: GPL-3.0-only
! Diagnostic driver of pinned upstream QLG routines; no spectral RT or SE solve.
program profile_probe
  use qlg, only: qlg_init, qlg_bot_sweep, qlg_ref_f64
  implicit none
  real(8) :: a, v, reference(1)
  real(4) :: bot(1)
  integer :: ios
  call qlg_init
  do
    read(*,*,iostat=ios) a,v
    if (ios<0) exit
    if (ios/=0 .or. a<=0 .or. abs(v)>100) stop 2
    call qlg_bot_sweep(real(a,4),real(v,4),0.0,1,bot)
    call qlg_ref_f64(a,v,0.0d0,1,reference)
    write(*,'(4(es25.16e3,1x))') a,v,real(bot(1),8),reference(1)
  end do
end program
