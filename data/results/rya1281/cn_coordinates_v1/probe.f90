program probe
 implicit none
 integer count,i
 double precision wave,centre
 real c,xl,width,a,vt,v,newvoigt
 read(*,*) count
 if(count<1.or.count>20000)stop 2
 do i=1,count
  read(*,*)wave,centre,c,xl,width,a
  vt=(wave-centre)*1.d-8
  vt=c*vt/xl**2
  v=vt/width
  write(*,'(I8,1P,2E25.16)')i,v,newvoigt(a,v)
 enddo
end program probe
