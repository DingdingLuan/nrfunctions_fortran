program main
implicit none
integer NN
parameter (NN=20)
integer i
real*8 :: a(NN),b(NN)
real*8 :: chsq,prob,df


a=0
b=0
!write(*,*)a(1)
do i=1,NN
	a(i)=1
enddo
!write(*,*)a(1)
do i=1,NN
	b(i)=2
ENDDO

do i=1,NN
	write(*,*)a(i),b(i)
enddo

call chstwo(a,b,NN,0,df,chsq,prob)


write(*,*) 'Chi-squared:',chsq
write(*,*)'Probability:',prob
write(*,*)df

end program main

include 'gser.f90'
include 'gcf.f90'
include 'gammq.f90'
include 'gammln.f90'
include 'chstwo.f90'
include 'expdev.f90'
include 'ran1.f90'
