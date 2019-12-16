! Monte Carlo simulation: jet-corrected luminosity
include 'num_pars.F90'
module global
implicit none
real*8,save :: Lb,v1,v2,Ep,norm
real*8,save :: AL,sigmaL,logLc
real*8,save :: A_ang,sigma_ang,log_angc
real*8,parameter :: omegam=0.3,omegal=0.7,h=0.71
end module global

program main
use global
use Numerical_Parameters, only : PI,RAD2DEGREES,M2CM
implicit none
integer*4,parameter :: nz=500,nL=500,num=4000000
integer*4 :: i,j,k,LM,idum,cnt(50),cnt2(50),idum2,cnt3(50)
real*8 :: z1,z2,ss,dl_m,dl_mpc,dtdz,dvdz,gammln,gammp,tot,tot2,L1,L2,LF2,sfr,ran1
real*8 :: x,y,k1,k2,Liso,ang1,ang2,ang,L_gam,psi,func,eta_z,eta_a
real*8,dimension(nz) :: z,tz,dndz,bin,f_metal
real*8,dimension(nL) :: logL,L,log_ang
real*8,dimension(num) :: z0,logL0,Fp,ang0
external ez,gammln,gammp,LF2,psi,band,bandE

!Calculate the pdf of redshift:
z1=0.001;z2=10.0
call tzt(z1,z2,nz,z,tz)
dndz=0.0
do i=1,nz
  call qromb(ez,0.0d0,z(i),ss)
  dl_m=9.26d25/h*(1.0+z(i))*ss  ! in units M
  dl_mpc=dl_m*3.24d-23 ! in Mpc
  dtdz=9.78d9/h/((1.0+z(i))*sqrt(omegal+omegam*(1.0+z(i))**3)) ! in yr
  dvdz=4.0*pi*3.065948d-7*dl_mpc**2/(1.0+z(i))*dtdz ! in Mpc^3
!  sfr=0.01*(1.0+z(i))**2.6/(1.0+((1.0+z(i))/3.2)**6.2) !SFR of Madau 2017
  if (z(i)<=1.0) then
    sfr=(1.0+z(i))**3.44
  else
    sfr=2.0**3.44
  end if
  f_metal(i)=gammp(0.84d0,0.4**2*10.0**(0.3*z(i)))
  dndz(i)=f_metal(i)*sfr*dvdz/(1+z(i))
enddo
call trapz(z,dndz,nz,tot)
open(11,file="pdf_redshift.txt")
do i=1,nz
  write(11,"(f7.3,2e11.3)") z(i),dndz(i)/tot,dndz(i)
end do

! Generate redshift z0:
call random_seed()
!idum=-1257
i=1;cnt=0
do while (i<=num)
  call random_number(x)
  call random_number(y)
  x=10.0*x
  y=0.2*y
!    x=10.0*ran1(idum)
!    y=0.2*ran1(idum)
    call locate(z,nz,x,j)
    if (y<=dndz(j)/tot) then
      z0(i)=x
      do k=1,50
        if (x>=(k-1)*0.2 .and. x<k*0.2) then
          cnt(k)=cnt(k)+1
        end if
      end do
      i=i+1
    end if
end do
open(21,file='result/input_z_dis.txt')
do k=1,50
  write(21,*) (k-0.5)*0.2,cnt(k)/0.2/real(num)
end do


! Compute the pdf of the jet-corrected luminosity:
L1=46.0;L2=52.0  !log10(L_gam)
logL=(/((L1+(L2-L1)/(nL-1)*LM),LM=0,nL-1)/)

!1.log-normal function:
! Set parameters:
AL=1.0
sigmaL=0.4
logLc=49.69

call qromb(LF2,L1,L2,tot)
AL=1.0/tot
open(12,file="pdf_lf2.txt")
do i=1,nL
  write(12,*) logL(i),LF2(logL(i))
end do

!2.broken power law:
! Set parameters:
! norm=1.0
! Lb=log10(2.0d52)
! v1=0.2
! v2=1.4
! call qromb(LF,L1,L2,tot)
! norm=1.0/tot

! Generate the jet-corrected luminosity:
!idum=-1358
!call random_seed()
i=1;cnt=0
do while (i<=num)
  call random_number(x)
  call random_number(y)
    x=(L2-L1)*x+L1
!    x=(L2-L1)*ran1(idum)+L1
!    y=ran1(idum)
    call locate(logL,nL,x,j)
    if (y<=LF2(logL(j))) then
      logL0(i)=x   !log10(L_gam)
      do k=1,10
        if (x>=L1+(k-1)*0.5 .and. x<L1+k*0.5) then
          cnt(k)=cnt(k)+1
        end if
      end do
      i=i+1
    end if
end do
open(22,file='result/input_Lgam_dis.txt')
do k=1,10
  write(22,*) L1+(k-0.5)*0.5,cnt(k)/0.5/real(num)
end do

! Compute the pdf of the jet opening angle:
ang1=-2.5;ang2=-0.1  ! log10(theta/rad)
log_ang=(/((ang1+(ang2-ang1)/(nL-1)*LM),LM=0,nL-1)/)

! Set parameters:
A_ang=1.0
sigma_ang=0.6
log_angc=-1.27

call qromb(psi,ang1,ang2,tot)
A_ang=1.0/tot
open(13,file="pdf_angle.txt")
do i=1,nL
  write(13,*) log_ang(i),psi(log_ang(i))
end do

! Generate the jet opening angle::
!idum=-1388
!call random_seed()
i=1;cnt=0
do while (i<=num)
  call random_number(x)
  call random_number(y)
    x=(ang2-ang1)*x+ang1
    y=0.8*y
!    x=(ang2-ang1)*ran1(idum)+ang1
!    y=0.7*ran1(idum)
    call locate(log_ang,nL,x,j)
    if (y<=psi(log_ang(j))) then
      ang0(i)=10.0**x  !angle in rad
      do k=1,12
        if (x>=(ang1+(k-1)*0.2) .and. x<(ang1+k*0.2)) then
          cnt(k)=cnt(k)+1
        end if
      end do
      i=i+1
    end if
end do
open(14,file='result/input_angle.txt')
do k=1,12
  write(14,*) ang1+(k-0.5)*0.2,cnt(k)/0.2/real(num)
end do

write(*,*) "End step input"

! Compute the peak flux using z and Liso:
!idum=-1322
call random_seed()
open(15,file='result/total_jc_sample.txt')
open(16,file='result/exp_jc_sample.txt')
open(17,file='result/exp_ang_dis.txt')
k=0;cnt=0;cnt2=0
iloop: do i=1,num
  ! Compute the luminosity distance at redshift z0:
  call qromb(ez,0.0d0,z0(i),ss)
  dl_m=9.26d25/h*(1.0+z0(i))*ss  !in units M
  ! Calculate Epeak using L-Ep relation:
  L_gam=10.0**logL0(i)
  Liso=L_gam/(1.0-cos(ang0(i)))  !Liso in erg/s
  call random_number(y)
  y=0.9*y+0.1
!  y=0.9*ran1(idum)+0.1
  Ep=200.0*(Liso/1.0d52)**0.5/y/(1.0+z0(i))
  ! Compute k-correction:
  call qromb(band,1.5d1,1.5d2,k1)
  call qromb(bandE,1.0d0/(1.0+z0(i)),1.0d4/(1.0+z0(i)),k2)
  k2=k2*1.60219d-9  ! keV to erg
  Fp(i)=Liso/(4.0*PI*(dl_m*1.0d2)**2)*k1/k2
!  ang=ang0(i)*RAD2DEGREES ! rad to degree   
  write(15,"(f7.3,e11.3,2f9.3)") z0(i),Liso,log10(ang0(i)),Fp(i)

! Calculate the detection probability:
  if (Fp(i)>0.2) then
    call random_number(x)
    if (x<func(Fp(i))) then
      eta_z=0.26+0.032*exp(1.61*log10(Fp(i)))
      call random_number(x)
      if (x<eta_z) then
        eta_a=1.0-cos(ang0(i))
        call random_number(x)
        if (x<0.2*eta_a) then
          k=k+1
          write(*,"(i4,f7.3,e11.3,2f10.3)") k,z0(i),Liso,log10(ang0(i)),Fp(i)
          write(16,"(f7.3,e11.3,2f10.3)") log10(1.0+z0(i)),Liso,log10(ang0(i)),Fp(i)
          do j=1,14
            if(log10(ang0(i))>=(ang1+(j-1)*0.2).and.log10(ang0(i))<(ang1+j*0.2)) then
              cnt(j)=cnt(j)+1
            end if
          end do
          do j=1,12
            if (log10(Liso)>=49.0+(j-1)*0.5 .and. log10(Liso)<49.0+j*0.5) then
              cnt2(j)=cnt2(j)+1
            end if
          end do
          do j=1,25
            if (log10(1.0+z0(i))>=(j-1)*0.04 .and. log10(1.0+z0(i))<j*0.04) then
              cnt3(j)=cnt3(j)+1
            end if
          end do
          if(k==150) exit
        end if
      end if
    end if
  end if
end do iloop

do j=1,14
  write(17,*) ang1+(j-0.5)*0.2,cnt(j)/real(k)
end do
open(18,file='result/exp_L_dis.txt')
do j=1,12
  write(18,*) 49.0+(j-0.5)*0.5,cnt2(j)/real(k)
end do
open(19,file='result/exp_z_dis.txt')
do j=1,25
  write(19,*) (j-0.5)*0.04,cnt3(j)/real(k)
end do

end



!------------------------------Functions----------------------------------
function ez(zprime)
use global, only : omegal, omegam
implicit none
real*8 zprime,ez
ez=1.0/sqrt(omegal+omegam*(1+zprime)**3)
end function


function band(E)
!band funciton for lgrb:
use global, only : Ep
implicit none
real*8 band,E,a,b

a=-1.0
b=-2.25

if (E<=(a-b)*Ep/(2.0+a)) then
	band=(E/100.0)**a*exp(-E*(2.0+a)/Ep)
else
	band=((a-b)*Ep/(100.0*(2.0+a)))**(a-b)*exp(b-a)*(E/100.0)**b
endif
end function

function bandE(E)
!band funciton for lgrb:
use global, only : Ep
implicit none
real*8 bandE,E,a,b

a=-1.0
b=-2.25

if (E<=(a-b)*Ep/(2.0+a)) then
	bandE=(E/100.0)**a*exp(-E*(2.0+a)/Ep)*E
else
	bandE=((a-b)*Ep/(100.0*(2.0+a)))**(a-b)*exp(b-a)*(E/100.0)**b*E
endif
end function

function LF(logL)
! borken power-law luminosity function
use global, only : Lb,v1,v2,norm
implicit none
real*8 :: LF,logL
if (logL<4.9d1) then
  LF=0.0
else if (logL<=Lb) then
  LF=norm*10.0**(-v1*(logL-Lb))
else 
  LF=norm*10.0**(-v2*(logL-Lb))
end if
end function

function LF2(logL)
! log-normal luminosity function
use global, only : AL,sigmaL,logLc
use Numerical_Parameters, only : PI
implicit none
real*8 LF2,logL
LF2=AL/sqrt(2.0*pi)/sigmaL*exp(-(logL-logLc)**2/(2.0*sigmaL**2))
end function

function psi(log_ang)
! the intrinsic distribution of angle
use global, only : A_ang,sigma_ang,log_angc
use Numerical_Parameters, only : PI
implicit none
real*8 psi,log_ang
psi=A_ang/sqrt(2.0*pi)/sigma_ang*exp(-(log_ang-log_angc)**2/(2.0*sigma_ang**2))
end function

function func(P)
! detection probability
implicit none
real*8 func,P
if (P<0.45) then
!  func=5.0*P**3.85
  func=P**2.0
else
  func=0.67*(1.0-0.4/P)**0.52
end if
end function

include 'mathlib/trapzd.f90'
include 'mathlib/trapz.f90'
include 'mathlib/qromb.f90'
include 'mathlib/gammp.f90'
include 'mathlib/gammln.f90'
include 'mathlib/erf.f90'
include 'lib/tzt.f90'
include 'mathlib/polint.f90'
include 'mathlib/gcf.f90'
include 'mathlib/gser.f90'
include 'lib/locate.f90'
include 'mathlib/ran1.f90'
