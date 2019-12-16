PROGRAM xpearsn
!C     driver for routine pearsn
INTEGER i
REAL*8 prob,r,z
REAL*8 dose(20),spore(20),aindex(20)
CHARACTER(len=10) head

! DATA dose/56.1,64.1,70.0,66.6,82.,91.3,90.,99.7,115.3,110./
! DATA spore/0.11,0.4,0.37,0.48,0.75,0.66,0.71,1.2,1.01,0.95/

open(17,file='/Users/dingding 1/Desktop/3.15/obsresult/fluence.txt')
open(18,file='/Users/dingding 1/Desktop/3.15/obsresult/output_fluence.txt')
read(18,*)head

do i=1,20
  read(17,*)aindex(i),dose(i)
  read(18,*)spore(i)
enddo

! do i=1,20
!   dose(i)=i**2
!   spore(i)=i**3
! enddo


write(*,'(1x,a)') &
&     'Effect of Gamma Rays on Man-in-the-Moon Marigolds'
write(*,'(1x,a,t29,a)') 'Count Rate (cpm)','Pollen Index'
do i=1,10
  write(*,'(1x,f10.2,f25.2)') dose(i),spore(i)
enddo
call pearsn(dose,spore,10,r,prob,z)
write(*,'(/1x,t24,a,t38,a)') 'PEARSN','Expected'
write(*,'(1x,a,t18,2e15.6)') 'Corr. Coeff.',r,0.906959
write(*,'(1x,a,t18,2e15.6)') 'Probability',prob,0.292650e-3
write(*,'(1x,a,t18,2e15.6/)') 'Fisher''s Z',z,1.51011
END

include 'betacf.f90'
include 'betai.f90'
include 'gammln.f90'
include 'pearsn.f90'