PROGRAM subdiffusion
use percroutines
IMPLICIT NONE

integer(kind=long)                    :: i,t,loop
!------------------------------------------------------------------------------------
call initialize
open(1, FILE=name//"_x.dat", STATUS = "REPLACE", ACTION="WRITE")
open(2, FILE=name//"_y.dat", STATUS = "REPLACE", ACTION="WRITE")
do loop=1,Mtra/10
   call cluster
   print*,"run",loop
   do t=2,ze
      call diffuse(t)
      if (mod(t,10000).eq.0) call dynobst
   enddo
   do i=1,10
      write(1,*) dx*pos(i,1:ze:100,1) 
      write(2,*) dx*pos(i,1:ze:100,2) 
   enddo
enddo
close(1)
close(2)
deallocate(orb)
END PROGRAM subdiffusion
