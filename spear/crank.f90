SUBROUTINE crank(n,w,s)
INTEGER n
REAL*8 s,w(n)
INTEGER j,ji,jt
REAL*8 rank,t
s=0.
j=1
1     if(j.lt.n)then
  if(w(j+1).ne.w(j))then
    w(j)=j
    j=j+1
  else
    do jt=j+1,n
      if(w(jt).ne.w(j))goto 2
    enddo
    jt=n+1
2         rank=0.5*(j+jt-1)
    do ji=j,jt-1
      w(ji)=rank
    enddo
    t=jt-j
    s=s+t**3-t
    j=jt
  endif
goto 1
endif
if(j.eq.n)w(n)=n
return
END
