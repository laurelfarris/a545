program test

real*8 :: x=2

100 print*,x

if (x==2) then
 goto 100
else
 goto 200
end if
 
!100 print*, x
200 print*, x+100



end program test
