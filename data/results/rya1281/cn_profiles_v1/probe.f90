program probe
  implicit none
  integer count,i
  real a,v,newvoigt
  read(*,*) count
  if (count<1.or.count>20000) stop 2
  do i=1,count
    read(*,*) a,v
    write(*,'(I8,1P,3E25.16)') i,a,v,newvoigt(a,v)
  enddo
end program probe
