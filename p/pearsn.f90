SUBROUTINE pearsn(x,y,n,r,prob,z)
INTEGER n
REAL*8 prob,r,z,x(n),y(n),TINY
PARAMETER (TINY=1.e-20)
!CU    USES betai
INTEGER j
REAL*8 ax,ay,df,sxx,sxy,syy,t,xt,yt,betai
ax=0.
ay=0.
do j=1,n
  ax=ax+x(j)
  ay=ay+y(j)
enddo
ax=ax/n
ay=ay/n
sxx=0.
syy=0.
sxy=0.
do j=1,n
  xt=x(j)-ax
  yt=y(j)-ay
  sxx=sxx+xt**2
  syy=syy+yt**2
  sxy=sxy+xt*yt
enddo
r=sxy/(sqrt(sxx*syy)+TINY)
z=0.5*log(((1.+r)+TINY)/((1.-r)+TINY))
df=n-2
t=r*sqrt(df/(((1.-r)+TINY)*((1.+r)+TINY)))
prob=betai(0.5*df,0.5,df/(df+t**2))
!C     prob=erfcc(abs(z*sqrt(n-1.))/1.4142136)
return
END
