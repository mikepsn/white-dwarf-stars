c-------------------------------------------------------------
c PROGRAM         : WHITE DWARF STARS
c-------------------------------------------------------------
      program main
      implicit none
      real*8 h, rhoc, rmax

      rmax = 2.55d0
      write(*,*)'Enter rhoc'
      read(*,*) rhoc
c	call euler(h, rhoc, rmax)
      call runge_kutta(h, rhoc, rmax)

      end
c-------------------------------------------------------------
c SUBROUTINE      : EULER
c-------------------------------------------------------------

      subroutine euler(h, rhoc, rmax)
      implicit none
      real*8 h, rhoc, rmax
      real*8 f, g, rho0, m0
      real*8 r, rho, m
      real*8 ri, rhoi, mi
      integer*4 n, i 
      
      open(unit=1, file='e_rho.dat', status='unknown')
      open(unit=2, file='e_m.dat', status='unknown')

      n = nint(rmax/h)

      r = h
      m = m0(rhoc, h)
      rho = rho0(rhoc, h)

      do 100 i = 0, n, +1
            ri = r
            rhoi = rho
            mi = m

            r = ri + h
            rho = rhoi + h*f(ri,rhoi,mi)
            m = mi + h*g(ri, rhoi, mi)

            write(1,*)r, rho
            write(2,*)r, m
100   continue

      close(unit=1)
      close(unit=2)

      end

c-------------------------------------------------------------
c SUBROUTINE      : RUNGE_KUTTA
c-------------------------------------------------------------

      subroutine runge_kutta(h, rhoc, rmax)
      implicit none
      real*8 h, rhoc, rmax
      real*8 f, g, rho0, m0
      real*8 r, rho, m
      real*8 ri, rhoi, mi, hon2
      real*8 f1, f2, f3, f4
      real*8 g1, g2, g3, g4
      integer*4 n, i 
      
      open(unit=1, file='rho.dat', status='unknown')
      open(unit=2, file='m.dat', status='unknown')

      n = nint(rmax/h)
      hon2 = h/2.0d0

      r = h
      m = m0(rhoc, h)
      rho = rho0(rhoc, h)

      do 100 i = 0, n, +1
            ri = r
            rhoi = rho
            mi = m

            f1 = f(ri, rhoi, mi)
            g1 = g(ri, rhoi, mi)

            f2 = f(ri + hon2, rhoi + hon2*f1, mi + hon2*g1)
            g2 = g(ri + hon2, rhoi + hon2*f1, mi + hon2*g1)

            f3 = f(ri + hon2, rhoi + hon2*f2, mi + hon2*g2)
            g3 = g(ri + hon2, rhoi + hon2*f2, mi + hon2*g2)

            f4 = f(ri + h, rhoi + h*f3, mi + h*g3)
            g4 = g(ri + h, rhoi + h*f3, mi + h*g3)

            r = ri + h
            rho = rhoi + (h/6.0d0)*(f1 + 2*f2 + 2*f3 + f4)
            m = mi + (h/6.0d0)*(g1 + 2*g2 + 2*g3 + g4)

            write(1,*)r, rho
            write(2,*)r, m, (r*0.005553956), (r*0.005154071d0)
100   continue

      close(unit=1)
      close(unit=2)

      end

c-------------------------------------------------------------
c FUNCTION        : f
c-------------------------------------------------------------

      double precision function f(r, rho, m)
      implicit none
      real*8 r, rho, m, x
      real*8 gama

	if (rho .lt. 0.0d0) then
		write(*,*)'Rho < 0 : Stopping Program'
		stop
	endif 

      f = -(m*rho)/(gama(rho**(1.0d0/3.0d0))*r**2.0d0)
      return

      end

c-------------------------------------------------------------
c FUNCTION        : g
c-------------------------------------------------------------

      double precision function g(r, rho, m)
      implicit none
      real*8 r, rho, m

      g = (r**2.0d0)*rho
      return

      end

c-------------------------------------------------------------
c FUNCTION        : gama
c-------------------------------------------------------------

      double precision function gama(x)
      implicit none
      real*8 x

      gama = (x**2.0d0)/(3.0d0*dsqrt(1.0d0 + x**2.0d0))
      return

      end

c-------------------------------------------------------------
c FUNCTION        : m0
c                 : initial value of m
c                 : rhoc = central density
c-------------------------------------------------------------

      double precision function m0(rhoc, h)
      implicit none
      real*8 rhoc, h

      m0 = rhoc*(h**3.0d0)/3.0d0
      return

      end

c-------------------------------------------------------------
c FUNCTION        : rho0
c                 : initial value of rho
c                 : rhoc = central density
c-------------------------------------------------------------

      double precision function rho0(rhoc, h)
      implicit none
      real*8 rhoc, h, x
      real*8 gama

      x = rhoc**(1.0d0/3.0d0)
      rho0 = rhoc*(1 - ((h**2.0d0)*rhoc)/(6*gama(x)) )
      return

      end
