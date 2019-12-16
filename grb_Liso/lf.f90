program main
implicit none
integer*4,parameter :: n=373
integer*4 :: i,j,k,cnt(10)
real*8,dimension(n) :: z,T90,S,L,logL,log10L,Eiso,Ep
real*8 ::dl_m,dl_cm,ss,pi,h,bin,lmin,x,zmax,a1,a2,kz,lmin2
real*8,dimension(n) :: alpha,norm
character(len=80) :: head
character(len=10),dimension(n) :: grbname
external ez,spec

h=0.71 
pi=3.14159265

cnt=0
k=0

open(15,file='Final429.txt')
!read(15,'(A80)') head
write(*,*) 'GRB list:'
do i=1,n
	read(15,*) grbname(i),T90(i),z(i),S(i)
	write(*,"(a10,f8.3,f9.4,e12.3)") grbname(i),T90(i),z(i),S(i)
enddo

open(13,file='Liso.txt')
write(13,*) "#GRBname     z      T90      Eiso       Liso       S"
do i=1,n
      !  calculate luminosity distance:
	call qromb(ez,0.0d0,z(i),ss)
	dl_m=9.26d25/h*(1.0+z(i))*ss  ! in units M
	dl_cm=dl_m*1.0d2 ! in cm
      !  calculate k-correction:
      call qromb(spec,1.5d1,1.5d2,a1)
	call qromb(spec,45.0/(1.0+z(i)),4.5d2/(1.0+z(i)),a2)
	kz=a2/a1
      kz=((1.0+z(i))/3.0)**-0.5
      !  calculate isotropic luminosity:
	Eiso(i)=4.0*pi*dl_cm**2*S(i)*kz/(1.0+z(i))
      L(i)=Eiso(i)/T90(i)*(1.0+z(i))

	write(13,"(a10,f7.3,f9.3,3e11.3)") grbname(i),z(i),T90(i),Eiso(i),L(i),S(i)
	write(*,"(a10,f7.3,2e11.3)") grbname(i),z(i),Eiso(i),L(i)
end do

!open(14,file='limit.txt')
!do x=0.1,8.0,0.02
!	call qromb(ez,0.0d0,x,ss)
!	dl_m=9.26d25/h*(1.0+x)*ss  ! in units M
!	dl_cm=dl_m*1.0d2 ! in cm
!	call qromb(spec,1.0/(1.0+x),1.0d4/(1.0+x),a2)
!	kz=a2/a1
!      kz=((1.0+x)/3.0)**-0.5
!	lmin=4.0*pi*dl_cm**2*1.0d-8*kz
!	write(14,*) x,log10(lmin)
!enddo


end

function ez(zprime)
implicit none
real*8 :: zprime,ez,omegal,omegam
omegal=0.734
omegam=0.266

ez=1.0/sqrt(omegal+omegam*(1.0+zprime)**3)
end function

function spec(E)
use global
implicit none
real*8 spec,E,a,b,Ep

a=-0.5
b=-2.3
Ep=490.0
if (E<=(a-b)*Ep/(2.+a)) then
	spec=(E/100.0)**a*exp(-E*(2+a)/Ep)*E
else
	spec=((a-b)*Ep/(100.0*(2.+a)))**(a-b)*exp(b-a)*(E/100.0)**b*E
endif
end function spec


include 'trapzd.f90'
include 'qromb.f90'
include 'polint.f90'
