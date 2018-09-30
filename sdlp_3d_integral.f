      subroutine sdlp_3d_integral 
     +
     +  (x0,y0,z0
     +  ,u0,v0,w0
     +  ,i0
     +  ,k
     +  ,mint
     +  ,u
     +  ,uxel,uyel,uzel
     +  )

c==========================================
c FDLIB, BEMLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved.
c
c This program is to be used only under the
c stipulations of the licensing agreement
c==========================================

c---------------------------------------
c integrate the dlp over the kth element
c---------------------------------------

      Implicit Double Precision (a-h,o-z)
      Dimension   p(1026,3)
      Dimension  ne(1026,7)
      Dimension vna(1026,3)
      Dimension   u(1026,3)
      Dimension   a(3078,3078)
 
      Dimension     n(512,6),nbe (512,3)
      Dimension alpha(512),  beta(512), gamma(512)

      Dimension xiq(20),etq(20),wq(20)

c--------------
c common blocks
c--------------

      common/points/p,ne
      common/elmnts/n,nbe
      common/albega/alpha,beta,gamma

      common/veloc3/a
      
      common/geo2/vna

      common/var/shrt,wall

      common/trq/xiq,etq,wq

      common/ph/ph1,ph2,ph3,ph4,ph5,ph6

c------
c flags
c------

      Iopt_sgf = 2    ! for the Green's function

c----------------
c  prepare to run
c----------------

      i1 = n(k,1)
      i2 = n(k,2)
      i3 = n(k,3)
      i4 = n(k,4)
      i5 = n(k,5)
      i6 = n(k,6)

      al = alpha(k)
      be =  beta(k)
      ga = gamma(k)

      uxel = 0.0D0
      uyel = 0.0D0
      uzel = 0.0D0

c---
c triangle quadrature
c---


      Do i=1,mint

        xi  = xiq(i)
        eta = etq(i)

        
        call sdlp_3d_interp 
     +
     +   (p(i1,1),p(i1,2),p(i1,3)
     +   ,p(i2,1),p(i2,2),p(i2,3)
     +   ,p(i3,1),p(i3,2),p(i3,3)
     +   ,p(i4,1),p(i4,2),p(i4,3)
     +   ,p(i5,1),p(i5,2),p(i5,3)
     +   ,p(i6,1),p(i6,2),p(i6,3)
     +
     +   ,vna(i1,1),vna(i1,2),vna(i1,3)
     +   ,vna(i2,1),vna(i2,2),vna(i2,3)
     +   ,vna(i3,1),vna(i3,2),vna(i3,3)
     +   ,vna(i4,1),vna(i4,2),vna(i4,3)
     +   ,vna(i5,1),vna(i5,2),vna(i5,3)
     +   ,vna(i6,1),vna(i6,2),vna(i6,3)
     +
     +   ,u(i1,1),u(i1,2),u(i1,3)
     +   ,u(i2,1),u(i2,2),u(i2,3)
     +   ,u(i3,1),u(i3,2),u(i3,3)
     +   ,u(i4,1),u(i4,2),u(i4,3)
     +   ,u(i5,1),u(i5,2),u(i5,3)
     +   ,u(i6,1),u(i6,2),u(i6,3)
     +
     +   ,al,be,ga
     +   ,xi,eta
     +
     +   ,x,y,z
     +   ,vnx,vny,vnz
     +   ,ux,uy,uz
     +   ,hs
     +       )

c        write (6,*) hs,ph1,ph2,ph3

c---
c compute the Green's function
c---

      call sgf_3d_w 
     +
     +  (Iopt_sgf
     +  ,y,z,x
     +  ,y0,z0,x0
     +  ,wall
     +  ,gyy,gyz,gyx
     +  ,gzy,gzz,gzx
     +  ,gxy,gxz,gxx
     +  ,py,pz,px
     +  ,Tyyy,Tyyz,Tyyx,Tzyz,Tzyx,Txyx
     +  ,Tyzy,Tyzz,Tyzx,Tzzz,Tzzx,Txzx
     +  ,Tyxy,Tyxz,Tyxx,Tzxz,Tzxx,Txxx
     +  )

      Txyy = Tyyx
      Txyz = Tzyx
      Tzyy = Tyyz

      Txzy = Tyzx
      Txzz = Tzzx
      Tzzy = Tyzz

      Txxy = Tyxx
      Txxz = Tzxx
      Tzxy = Tyxz

      ux = ux-u0 
      uy = uy-v0
      uz = uz-w0

      cf = 0.5D0*hs*wq(i)

c     ux = 0.0    ! for debugging
c     uy = 1.0
c     uz = 0.0

      uxel = uxel +((ux*Txxx + uy*Tyxx + uz*Tzxx)*vnx
     +             +(ux*Txxy + uy*Tyxy + uz*Tzxy)*vny
     +             +(ux*Txxz + uy*Tyxz + uz*Tzxz)*vnz)*cf

      uyel = uyel +((ux*Txyx + uy*Tyyx + uz*Tzyx)*vnx
     +             +(ux*Txyy + uy*Tyyy + uz*Tzyy)*vny
     +             +(ux*Txyz + uy*Tyyz + uz*Tzyz)*vnz)*cf

      uzel = uzel +((ux*Txzx + uy*Tyzx + uz*Tzzx)*vnx
     +             +(ux*Txzy + uy*Tyzy + uz*Tzzy)*vny
     +             +(ux*Txzz + uy*Tyzz + uz*Tzzz)*vnz)*cf

      ! u(i,1) <- u(i1~i6,1)
      a(3*i0-2,3*i1-2) =  a(3*i0-2,3*i1-2)
     +     + cf*ph1*(Txxx*vnx + Txxy*vny + Txxz*vnz)
      
      a(3*i0-2,3*i2-2) =  a(3*i0-2,3*i2-2)
     +     + cf*ph2*(Txxx*vnx + Txxy*vny + Txxz*vnz)

      a(3*i0-2,3*i3-2) =  a(3*i0-2,3*i3-2)
     +     + cf*ph3*(Txxx*vnx + Txxy*vny + Txxz*vnz)

      a(3*i0-2,3*i4-2) =  a(3*i0-2,3*i4-2)
     +     + cf*ph4*(Txxx*vnx + Txxy*vny + Txxz*vnz)

      a(3*i0-2,3*i5-2) =  a(3*i0-2,3*i5-2)
     +     + cf*ph5*(Txxx*vnx + Txxy*vny + Txxz*vnz)

      a(3*i0-2,3*i6-2) =  a(3*i0-2,3*i6-2)
     +     + cf*ph6*(Txxx*vnx + Txxy*vny + Txxz*vnz)

      ! u(i,1) <- u(i1~6,2)
      a(3*i0-2,3*i1-1) =  a(3*i0-2,3*i1-1)
     +     + cf*ph1*(Tyxx*vnx + Tyxy*vny + Tyxz*vnz)

      a(3*i0-2,3*i2-1) =  a(3*i0-2,3*i2-1)
     +     + cf*ph2*(Tyxx*vnx + Tyxy*vny + Tyxz*vnz)

      a(3*i0-2,3*i3-1) =  a(3*i0-2,3*i3-1)
     +     + cf*ph3*(Tyxx*vnx + Tyxy*vny + Tyxz*vnz)

      a(3*i0-2,3*i4-1) =  a(3*i0-2,3*i4-1)
     +     + cf*ph4*(Tyxx*vnx + Tyxy*vny + Tyxz*vnz)

      a(3*i0-2,3*i5-1) =  a(3*i0-2,3*i5-1)
     +     + cf*ph5*(Tyxx*vnx + Tyxy*vny + Tyxz*vnz)

      a(3*i0-2,3*i6-1) =  a(3*i0-2,3*i6-1)
     +     + cf*ph6*(Tyxx*vnx + Tyxy*vny + Tyxz*vnz)

      ! u(i,1) <- u(i1~6,3)
      a(3*i0-2,3*i1) =  a(3*i0-2,3*i1)
     +     + cf*ph1*(Tzxx*vnx + Tzxy*vny + Tzxz*vnz)

      a(3*i0-2,3*i2) =  a(3*i0-2,3*i2)
     +     + cf*ph2*(Tzxx*vnx + Tzxy*vny + Tzxz*vnz)

      a(3*i0-2,3*i3) =  a(3*i0-2,3*i3)
     +     + cf*ph3*(Tzxx*vnx + Tzxy*vny + Tzxz*vnz)

      a(3*i0-2,3*i4) =  a(3*i0-2,3*i4)
     +     + cf*ph4*(Tzxx*vnx + Tzxy*vny + Tzxz*vnz)

      a(3*i0-2,3*i5) =  a(3*i0-2,3*i5)
     +     + cf*ph5*(Tzxx*vnx + Tzxy*vny + Tzxz*vnz)

      a(3*i0-2,3*i6) =  a(3*i0-2,3*i6)
     +     + cf*ph6*(Tzxx*vnx + Tzxy*vny + Tzxz*vnz)

      
      ! u(i,2) <- u(i1~i6,1)
      a(3*i0-1,3*i1-2) =  a(3*i0-1,3*i1-2)
     +     + cf*ph1*(Txyx*vnx + Txyy*vny + Txyz*vnz)
      
      a(3*i0-1,3*i2-2) =  a(3*i0-1,3*i2-2)
     +     + cf*ph2*(Txyx*vnx + Txyy*vny + Txyz*vnz)

      a(3*i0-1,3*i3-2) =  a(3*i0-1,3*i3-2)
     +     + cf*ph3*(Txyx*vnx + Txyy*vny + Txyz*vnz)

      a(3*i0-1,3*i4-2) =  a(3*i0-1,3*i4-2)
     +     + cf*ph4*(Txyx*vnx + Txyy*vny + Txyz*vnz)

      a(3*i0-1,3*i5-2) =  a(3*i0-1,3*i5-2)
     +     + cf*ph5*(Txyx*vnx + Txyy*vny + Txyz*vnz)

      a(3*i0-1,3*i6-2) =  a(3*i0-1,3*i6-2)
     +     + cf*ph6*(Txyx*vnx + Txyy*vny + Txyz*vnz)

      ! u(i,2) <- u(i1~6,2)
      a(3*i0-1,3*i1-1) =  a(3*i0-1,3*i1-1)
     +     + cf*ph1*(Tyyx*vnx + Tyyy*vny + Tyyz*vnz)

      a(3*i0-1,3*i2-1) =  a(3*i0-1,3*i2-1)
     +     + cf*ph2*(Tyyx*vnx + Tyyy*vny + Tyyz*vnz)

      a(3*i0-1,3*i3-1) =  a(3*i0-1,3*i3-1)
     +     + cf*ph3*(Tyyx*vnx + Tyyy*vny + Tyyz*vnz)

      a(3*i0-1,3*i4-1) =  a(3*i0-1,3*i4-1)
     +     + cf*ph4*(Tyyx*vnx + Tyyy*vny + Tyyz*vnz)

      a(3*i0-1,3*i5-1) =  a(3*i0-1,3*i5-1)
     +     + cf*ph5*(Tyyx*vnx + Tyyy*vny + Tyyz*vnz)

      a(3*i0-1,3*i6-1) =  a(3*i0-1,3*i6-1)
     +     + cf*ph6*(Tyyx*vnx + Tyyy*vny + Tyyz*vnz)

      ! u(i,2) <- u(i1~6,3)
      a(3*i0-1,3*i1) =  a(3*i0-1,3*i1)
     +     + cf*ph1*(Tzyx*vnx + Tzyy*vny + Tzyz*vnz)

      a(3*i0-1,3*i2) =  a(3*i0-1,3*i2)
     +     + cf*ph2*(Tzyx*vnx + Tzyy*vny + Tzyz*vnz)

      a(3*i0-1,3*i3) =  a(3*i0-1,3*i3)
     +     + cf*ph3*(Tzyx*vnx + Tzyy*vny + Tzyz*vnz)

      a(3*i0-1,3*i4) =  a(3*i0-1,3*i4)
     +     + cf*ph4*(Tzyx*vnx + Tzyy*vny + Tzyz*vnz)

      a(3*i0-1,3*i5) =  a(3*i0-1,3*i5)
     +     + cf*ph5*(Tzyx*vnx + Tzyy*vny + Tzyz*vnz)

      a(3*i0-1,3*i6) =  a(3*i0-1,3*i6)
     +     + cf*ph6*(Tzyx*vnx + Tzyy*vny + Tzyz*vnz)

      
      ! u(i,3) <- u(i1~i6,1)
      a(3*i0,3*i1-2) =  a(3*i0,3*i1-2)
     +     + cf*ph1*(Txzx*vnx + Txzy*vny + Txzz*vnz)
      
      a(3*i0,3*i2-2) =  a(3*i0,3*i2-2)
     +     + cf*ph2*(Txzx*vnx + Txzy*vny + Txzz*vnz)

      a(3*i0,3*i3-2) =  a(3*i0,3*i3-2)
     +     + cf*ph3*(Txzx*vnx + Txzy*vny + Txzz*vnz)

      a(3*i0,3*i4-2) =  a(3*i0,3*i4-2)
     +     + cf*ph4*(Txzx*vnx + Txzy*vny + Txzz*vnz)

      a(3*i0,3*i5-2) =  a(3*i0,3*i5-2)
     +     + cf*ph5*(Txzx*vnx + Txzy*vny + Txzz*vnz)

      a(3*i0,3*i6-2) =  a(3*i0,3*i6-2)
     +     + cf*ph6*(Txzx*vnx + Txzy*vny + Txzz*vnz)

      ! u(i,3) <- u(i1~6,2)
      a(3*i0,3*i1-1) =  a(3*i0,3*i1-1)
     +     + cf*ph1*(Tyzx*vnx + Tyzy*vny + Tyzz*vnz)

      a(3*i0,3*i2-1) =  a(3*i0,3*i2-1)
     +     + cf*ph2*(Tyzx*vnx + Tyzy*vny + Tyzz*vnz)

      a(3*i0,3*i3-1) =  a(3*i0,3*i3-1)
     +     + cf*ph3*(Tyzx*vnx + Tyzy*vny + Tyzz*vnz)

      a(3*i0,3*i4-1) =  a(3*i0,3*i4-1)
     +     + cf*ph4*(Tyzx*vnx + Tyzy*vny + Tyzz*vnz)

      a(3*i0,3*i5-1) =  a(3*i0,3*i5-1)
     +     + cf*ph5*(Tyzx*vnx + Tyzy*vny + Tyzz*vnz)

      a(3*i0,3*i6-1) =  a(3*i0,3*i6-1)
     +     + cf*ph6*(Tyzx*vnx + Tyzy*vny + Tyzz*vnz)

      ! u(i,3) <- u(i1~6,3)
      a(3*i0,3*i1) =  a(3*i0,3*i1)
     +     + cf*ph1*(Tzzx*vnx + Tzzy*vny + Tzzz*vnz)

      a(3*i0,3*i2) =  a(3*i0,3*i2)
     +     + cf*ph2*(Tzzx*vnx + Tzzy*vny + Tzzz*vnz)

      a(3*i0,3*i3) =  a(3*i0,3*i3)
     +     + cf*ph3*(Tzzx*vnx + Tzzy*vny + Tzzz*vnz)

      a(3*i0,3*i4) =  a(3*i0,3*i4)
     +     + cf*ph4*(Tzzx*vnx + Tzzy*vny + Tzzz*vnz)

      a(3*i0,3*i5) =  a(3*i0,3*i5)
     +     + cf*ph5*(Tzzx*vnx + Tzzy*vny + Tzzz*vnz)

      a(3*i0,3*i6) =  a(3*i0,3*i6)
     +     + cf*ph6*(Tzzx*vnx + Tzzy*vny + Tzzz*vnz)
      
      


      

      End Do
c      write (6,*) i0,k

c-----
c done
c-----

      return
      end
