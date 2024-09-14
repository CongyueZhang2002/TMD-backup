c--------------------------------------------------------------------------
c               / xf
c               |
c       value = | dx func(x)         n > 1, number of divisions for xf-xi
c               |
c              / xi
c--------------------------------------------------------------------------

      function qgauss(func,xi,xf,n)
      implicit none
      real*8 qgauss,func,xi,xf,xn,value,x1,x2,val
      integer i,n
      external func

      if(n.le.1) then           ! same as n=1
         x1=xi
         x2=xf
         call gauss(func,x1,x2,val)
         qgauss=val
         return
      endif
	
      xn=(xf-xi)/float(n)
      value=0d0
      x2=xi
      Do 100 i=1,n
         x1=x2
         x2=x1+xn
         call gauss(func,x1,x2,val)
         value=value+val
 100  continue

      qgauss=value
      return
      end

c-------------------------------------------------
      subroutine gauss(func,xi,xf,value)
      implicit none
      real*8 func,xi,xf,value,xm,xr,dx,x(8),w(8)
      real*8 eps
      integer j
      data eps /1.0d-25/
      data w
     1   / 0.02715 24594 11754 09485 17805 725D0,
     2     0.06225 35239 38647 89286 28438 370D0,
     3     0.09515 85116 82492 78480 99251 076D0,
     4     0.12462 89712 55533 87205 24762 822D0,
     5     0.14959 59888 16576 73208 15017 305D0,
     6     0.16915 65193 95002 53818 93120 790D0,
     7     0.18260 34150 44923 58886 67636 680D0,
     8     0.18945 06104 55068 49628 53967 232D0 /
      DATA X
     1   / 0.98940 09349 91649 93259 61541 735D0,
     2     0.94457 50230 73232 57607 79884 155D0,
     3     0.86563 12023 87831 74388 04678 977D0,
     4     0.75540 44083 55003 03389 51011 948D0,
     5     0.61787 62444 02643 74844 66717 640D0,
     6     0.45801 67776 57227 38634 24194 430D0,
     7     0.28160 35507 79258 91323 04605 015D0,
     8     0.09501 25098 37637 44018 53193 354D0 /
      
      xm=0.5d0*(xf+xi)
      xr=0.5d0*(xf-xi)
      if (abs(xr).lt.eps) print *,
     >     'WARNING: Too high accuracy required for QGAUSS!'
 
      value=0d0

      Do 100 j=1,8
         dx=xr*x(j)
         value=value+w(j)*(func(xm+dx)+func(xm-dx))
 100  continue
	
      value=xr*value
      return
      end
