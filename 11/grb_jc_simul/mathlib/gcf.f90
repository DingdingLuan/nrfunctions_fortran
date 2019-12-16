SUBROUTINE gcf(gammcf,a,x,gln)
INTEGER*4 ITMAX
REAL*8 a,gammcf,gln,x,EPS,FPMIN
PARAMETER (ITMAX=100,EPS=3.d-7,FPMIN=1.d-30)
!CU    USES gammln
INTEGER*4 i
REAL*8 an,b,c,d,del,h,gammln
EXTERNAL gammln
gln=gammln(a)
b=x+1.-a
c=1./FPMIN
d=1./b
h=d
do i=1,ITMAX
	an=-i*(i-a)
	b=b+2.
	d=an*d+b
	if(abs(d).lt.FPMIN)d=FPMIN
	c=b+an/c
	if(abs(c).lt.FPMIN)c=FPMIN
	d=1./d
	del=d*c
	h=h*del
	if(abs(del-1.).lt.EPS)goto 1
end do
pause 'a too large, ITMAX too small in gcf'
1     gammcf=exp(-x+a*log(x)-gln)*h
return
END
