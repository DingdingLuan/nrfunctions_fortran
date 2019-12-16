PROGRAM xchstwo
!C     driver for routine chstwo
INTEGER n
PARAMETER(n=10)
INTEGER i
REAL chsq,df,prob
real a1(n),a2(n)
! idum=-17
! do j=1,NBINS
!   bins1(j)=0.0
!   bins2(j)=0.0
! enddo
! do i=1,NPTS
!   x=expdev(idum)
!   ibin=x*NBINS/3.0+1
!   if (ibin.le.NBINS) bins1(ibin)=bins1(ibin)+1.0
!   x=expdev(idum)
!   ibin=x*NBINS/3.0+1
!   if (ibin.le.NBINS) bins2(ibin)=bins2(ibin)+1.0
! enddo
a2(i)=0.
do i=1,n
	a1(i)=i+1
	a2(i+2)=i+1
enddo

call chstwo(a1,a2,n,0,df,chsq,prob)
write(*,'(1x,t10,a,t25,a)') 'Dataset 1','Dataset 2'
do i=1,n
  write(*,'(1x,2f15.2)') a1(i),a2(i)
enddo
write(*,'(/1x,t10,a,e12.4)') 'Chi-squared:',chsq
write(*,'(1x,t10,a,e12.4)') 'Probability:',prob
write(*,*)df
END


include 'gser.f90'
include 'gcf.f90'
include 'gammq.f90'
include 'gammln.f90'
include 'chstwo.f90'
include 'ran1.f90'
include 'expdev.f90'