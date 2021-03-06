      subroutine sgf_3d_sw 
     +
     +   (Iopt
     +   ,x,y,z
     +   ,x0,y0,z0
     +   ,wall
     +
     +   ,Gxx,Gxy,Gxz
     +   ,Gyx,Gyy,Gyz
     +   ,Gzx,Gzy,Gzz
     +   ,px,py,pz
     +   ,Txxx,Txxy,Txxz,Tyxy,Tyxz,Tzxz
     +   ,Txyx,Txyy,Txyz,Tyyy,Tyyz,Tzyz
     +   ,Txzx,Txzy,Txzz,Tyzy,Tyzz,Tzzz
     +   )

c=========================================
c Green's function for semi-infinite flow
c bounded by a plane wall located at x = wall
c The wall is perfectly slip consition
c      
c Iopt =  1 Compute only the Green's function
c      ne 1 Compute the Green's function,
c           pressure, and stress
c=========================================

      Implicit Double Precision (a-h,o-z)

      pi = 3.14159 265358 979 32384 D0

      piq = 0.25D0*pi
      pih = 0.50D0*pi
      pi2 = 2.00D0*pi
      pi4 = 4.00D0*pi
      pi6 = 6.00D0*pi
      pi8 = 8.00D0*pi

c-------------
c primary part (Compute original source point)
c-------------

      dx = x-x0
      dy = y-y0
      dz = z-z0

      dxx = dx*dx
      dxy = dx*dy
      dxz = dx*dz
      dyy = dy*dy
      dyz = dy*dz
      dzz = dz*dz

      r   = Dsqrt(dxx+dyy+dzz)
      r3  = r*r*r
      r5  = r3*r*r
      ri  = 1.0D0/r
      ri3 = 1.0D0/r3
      ri5 = 1.0D0/r5
! GA
      Sxx = ri + dxx*ri3
      Sxy =      dxy*ri3
      Sxz =      dxz*ri3
      Syy = ri + dyy*ri3
      Syz =      dyz*ri3
      Szz = ri + dzz*ri3

c------------------------------
      if(Iopt.ne.1) then    ! compute the stress and pressure
c$$$
c$$$! compute dG[ij]/dx[k]
c$$$         TA0xxx = dx *ri3 - 3.0D0 *dxx *dx *ri5
c$$$         TA0xxy = - dy *ri3 - 3.0D0 *dxy *dx *ri5
c$$$         TA0xxz = - dz *ri3 - 3.0D0 *dxz *dx *ri5
c$$$         TA0yxy = dx *ri3 - 3.0D0 *dxy *dy *ri5
c$$$         TA0yxz = - 3.0D0 *dxy *dz *ri5
c$$$         TA0zxz = dx *ri3 - 3.0D0 *dxz *dz *ri5
c$$$
c$$$         TA0xyx = dy *ri3 - 3.0D0 *dxy *dx *ri5
c$$$         TA0xyy = dx *ri3 - 3.0D0 *dxy *dy *ri5
c$$$         TA0xyz = - 3.0D0 *dxy *dz *ri5
c$$$         TA0yyy = dy *ri3 - 3.0D0 *dyy *dy *ri5
c$$$         TA0yyz = - dz *ri3 - 3.0D0 *dyz *dy *ri5
c$$$         TA0zyz = dy *ri3 - 3.0D0 *dyz *dz *ri5
c$$$
c$$$         TA0xzx = dz *ri3 - 3.0D0 *dxz *dx *ri5
c$$$         TA0xzy = - 3.0D0 *dxy *dz *ri5
c$$$         TA0xzz = dx *ri3 - 3.0D0 *dxz *dz *ri5
c$$$         TA0yzy = dz *ri3 - 3.0D0 *dyz *dy *ri5
c$$$         TA0yzz = dy *ri3 - 3.0D0 *dyz *dz *ri5
c$$$         TA0zzz = dz *ri3 - 3.0D0 *dzz *dz *ri5
c$$$!     compute dG[kj]/dx[i]
c$$$         TAxxx = TA0xxx + TA0xxx
c$$$         TAxxy = TA0xxy + dy *ri3 - 3.0D0 *dxy *dy *ri5
c$$$         TAxxz = TA0xxz + dz *ri3 - 3.0D0 *dxz *dx *ri5
c$$$         TAyxy = TA0yxy + TA0yxy
c$$$         TAyxz = TA0yxz - 3.0D *dxy *dz *ri5
c$$$         TAzxz = TA0zxz + TA0zxz
c$$$
c$$$         TAxyx = TA0xyx + TA0xyx
c$$$         TAxyy = TA0xyy - dx *ri3 - 3.0D0 *dxy *dy *ri5
c$$$         TAxyz = TA0xyz - 3.0D0 *dxy *dz *ri5
c$$$         TAyyy = TA0yyy + TA0xyz
c$$$         TAyyz = TA0yyz + dz *ri3 - 3.0D0 *dyz *dy *ri5
c$$$         TAzyz = TA0zyz + TA0zyz
c$$$
c$$$         TAxzx = TA0xzx +TA0xzx
c$$$         TAxzy = TA0xzy - 3.0D0 *dxy *dz *ri5
c$$$         TAxzz = TA0xzz - dx *ri3 - 3.0D0 *dxz *dz *ri5
c$$$         TAyzy = TA0yzy + TA0yzy
c$$$         TAyzz = TA0yzz - dy *ri3 - 3.0D0 *dyz *dz *ri5
c$$$         TAzzz = TA0zzz + TA0zzz

         
      cf = -6.0D0*ri5

      TSxxx = dxx*dx * cf
      TSxxy = dxy*dx * cf
      TSxxz = dxz*dx * cf
      TSyxy = dyy*dx * cf
      TSyxz = dyz*dx * cf
      TSzxz = dzz*dx * cf

      TSyyy = dyy*dy * cf
      TSyyz = dyz*dy * cf
      TSzyz = dzz*dy * cf

      TSyzy = dyy*dz * cf
      TSyzz = dyz*dz * cf
      TSzzz = dzz*dz * cf

      cf = 2.0D0*ri3

      psx = dx * cf
      psy = dy * cf
      psz = dz * cf

         
         
      

      End If
c------------------------------

c-------------
c image system
c-------------

      x0im = 2.0D0*wall-x0
      dx   = x-x0im

      dxx = dx*dx
      dxy = dx*dy
      dxz = dx*dz
      
      r   = Dsqrt(dxx+dyy+dzz)
      r3  = r*r*r
      r5  = r3*r*r
      ri  = 1.0D0/r
      ri3 = 1.0D0/r3
      ri5 = 1.0D0/r5

      r1i = 1.0D0/(r + dx)

c----------------
c image stokeslet : GB
c----------------

      Sxx = Sxx - ri - dxx*ri3
      Sxz = Sxz      - dxz*ri3
      Syy = Syy - ri - dyy*ri3
      Syz = Syz      - dyz*ri3
      Szz = Szz - ri - dzz*ri3

      Syx = Sxy
      Szx = Sxz
      Szy = Syz

c-----------
      if(Iopt.ne.1) then   ! compute the stress and pressure

      cf = -6.0D0*ri5

      TSxxx = TSxxx - dxx*dx * cf
      TSxxy = TSxxy - dxy*dx * cf
      TSxxz = TSxxz - dxz*dx * cf
      TSyxy = TSyxy - dyy*dx * cf
      TSyxz = TSyxz - dyz*dx * cf
      TSzxz = TSzxz - dzz*dx * cf

      TSxyx = TSxxy
      TSxyy = TSyxy
      TSxyz = TSyxz
      TSyyy = TSyyy - dyy*dy * cf
      TSyyz = TSyyz - dyz*dy * cf
      TSzyz = TSzyz - dzz*dy * cf

      TSxzx = TSxxz
      TSxzy = TSyxz
      TSxzz = TSzxz
      TSyzy = TSyzy - dyy*dz * cf
      TSyzz = TSyzz - dyz*dz * cf
      TSzzz = TSzzz - dzz*dz * cf

      cf = 2.0D0*ri3

      psx = psx - dx * cf
      psy = psy - dy * cf
      psz = psz - dz * cf

      end if
c-----------

c-----------------------
c image potential dipole : GD
c-----------------------

      PDxx = - ri3 + 3.0D0*dxx*ri5 
      PDyx =         3.0D0*dxy*ri5
      PDzx =         3.0D0*dxz*ri5

      PDxy = - PDyx
      PDyy =   ri3 - 3.0D0*dyy*ri5
      PDzy =       - 3.0D0*dyz*ri5

      PDxz = - PDzx
      PDyz =   PDzy
      PDzz =   ri3 - 3.0D0*dzz*ri5

c-----------------------
      if(Iopt.ne.1) then

      r7  = r5*r*r
      ri7 = 1.0D0/r7

      cf  =  6.0D0*ri5
      cf1 = 30.0D0*ri7

      TPDxxx =  (dx+dx+dx) * cf - dxx*dx * cf1
      TPDxxy =   dy        * cf - dxy*dx * cf1
      TPDxxz =   dz        * cf - dxz*dx * cf1
      TPDyxy =   dx        * cf - dyy*dx * cf1
      TPDyxz =                  - dyz*dx * cf1
      TPDzxz =   dx        * cf - dzz*dx * cf1

      TPDxyx = - dy        * cf + dxx*dy * cf1
      TPDxyy = - dx        * cf + dxy*dy * cf1
      TPDxyz =                    dxz*dy * cf1
      TPDyyy = -(dy+dy+dy) * cf + dyy*dy * cf1
      TPDyyz = - dz        * cf + dyz*dy * cf1
      TPDzyz = - dy        * cf + dzz*dy * cf1

      TPDxzx = - dz        * cf + dxx*dz * cf1
      TPDxzy =                    dxy*dz * cf1
      TPDxzz = - dx        * cf + dxz*dz * cf1
      TPDyzy = - dz        * cf + dyy*dz * cf1
      TPDyzz = - dy        * cf + dyz*dz * cf1
      TPDzzz = -(dz+dz+dz) * cf + dzz*dz * cf1

c     no pressure contribution

      end if
c-----------------------

c-----------------------
c image stokeslet dipole : GC
c-----------------------

      SDxx = dx * PDxx
      SDyx = dx * PDyx - dy*ri3
      SDzx = dx * PDzx - dz*ri3

      SDxy = dx * PDxy - dy*ri3
      SDyy = dx * PDyy
      SDzy = dx * PDzy

      SDxz = dx * PDxz - dz*ri3
      SDyz = dx * PDyz
      SDzz = dx * PDzz

c-----------------------
      if(Iopt.ne.1) then

      cf = 6.0D0*ri5

      TSDxxx = dx *TPDxxx 
      TSDxxy = dx *TPDxxy - (   -dxy) * cf
      TSDxxz = dx *TPDxxz - (   -dxz) * cf
      TSDyxy = dx *TPDyxy - (dxx-dyy) * cf
      TSDyxz = dx *TPDyxz - (   -dyz) * cf
      TSDzxz = dx *TPDzxz - (dxx-dzz) * cf

      TSDxyx = dx *TPDxyx + dxy * cf
      TSDxyy = dx *TPDxyy 
      TSDxyz = dx *TPDxyz 
      TSDyyy = dx *TPDyyy + dxy * cf
      TSDyyz = dx *TPDyyz 
      TSDzyz = dx *TPDzyz + dxy * cf

      TSDxzx = dx *TPDxzx + dxz * cf
      TSDxzy = dx *TPDxzy 
      TSDxzz = dx *TPDxzz 
      TSDyzy = dx *TPDyzy + dxz * cf
      TSDyzz = dx *TPDyzz
      TSDzzz = dx *TPDzzz + dxz * cf

      PSDx =  - 2.0D0*ri3 + 6.0D0*dxx*ri5
      PSDy =              - 6.0D0*dxy*ri5
      PSDz =              - 6.0D0*dxz*ri5

      end if
c-----------------------


c---------
c GE
c---------
      GExx = 0.0D0
      GEyx = 0.0D0
      GEzx = 0.0D0

      GExy = - pi8 *(- dxy *ri3 + dy *ri *r1i)
      GEyy = - pi8 *(ri - dyy *ri3)
      GEzy = - pi8 *(- dyz *ri3)

      GExz = - pi8 *(- dxz *ri3 + dz *ri *r1i)
      GEyz = GEzy
      GEzz = - pi8 *(ri - dzz *ri3)

! compute dG[ij]/dx[k]
      TE0xxx = 0.0D0
      TE0xxy = 0.0D0
      TE0xxz = 0.0D0
      TE0yxy = 0.0D0
      TE0yxz = 0.0D0
      TE0zxz = 0.0D0

      TE0xyx = - 2.0D0 *dy *ri3 + 3.0D0 *dx *dxy *ri5
      TE0xyy = - dx*ri3 + 3.0D0 *dy *dxy *ri5 + ri *r1i
     +        - dy *dy *ri *ri *r1i *r1i - dy *dy *ri3 *r1i
      TE0xyz = - 3.0D0 *dxy *dz *ri5 - dyz *ri3 *r1i
     +        - dyz *ri *ri *r1i *r1i
      TE0yyy = - 3.0D0 *dy *ri3 + 3.0D0 *dyy *dy *ri5
      TE0yyz = - dz *ri3 + 3.0D0 *dyz *dy *ri5
      TE0zyz = - dy *ri3 + 3.0D0 *dyz *dz *ri5

      TE0xzx = - 2.0D0 *dz *ri3 + 3.0D0 *dx *dxz *ri5
      TE0xzy = 3.0D0 *dxy *dz *ri5 - dyz *ri3 *r1i - dyz *ri *ri *r1i
      TE0xzz = - dx *ri3 + 3.0D0 *dxz *dz *ri5 + ri *r1i
     +     - dzz *dz *ri3 *r1i - dzz *dz *ri *ri *r1i *r1i
      TE0yzy = -dz *ri3 + 3.0D0 *dyz *dy *ri5
      TE0yzz = -dy *ri3 + 3.0D0 *dyz *dz *ri5
      TE0zzz = - 3.0D0 *dz *ri3 - 3.0D0 *dzz *dz *ri5

! compute dG[kj]/dx[i]
      TExxx = TE0xxx + TE0xxx
      TExxy = TE0xxy
      TExxz = TE0xxz
      TEyxy = TE0yxy
      TEyxz = TE0yxz
      TEzxz = TE0zxz

      TExyx = TE0xyx + TE0xyx
      TExyy = TE0xyy - dx *ri3 + 3.0D0 *dxy *dy *ri5
      TExyz = TE0xyz + 3.0D0 *dxy *dz *ri5
      TEyyy = TE0yyy + TE0yyy
      TEyyz = TE0yyz - dz *ri3 + 3.0D0 *dyz *dy *ri5
      TEzyz = TE0zyz + TE0zyz

      TExzx = TE0xzx + TE0xzx
      TExzy = TE0xzy + 3.0D0 *dxy *dz *ri5
      TExzz = TE0xzz - dx *ri3 + 3.0D0 *dxz *dz *ri5
      TEyzy = TE0yzy + TE0yzy
      TEyzz = TE0yzz - dy *ri3 + 3.0D0 *dyz *dz *ri5
      TEzzz = TE0zzz + TE0zzz

!     *-pi8
      TExxx = - pi8 *TExxx
      TExxy = - pi8 *TExxy
      TExxz = - pi8 *TExxz
      TEyxy = - pi8 *TEyxy
      TEyxz = - pi8 *TEyxz
      TEzxz = - pi8 *TEzxz

      TExyx = - pi8 *TExyx
      TExyy = - pi8 *TExyy
      TExyz = - pi8 *TExyz
      TEyyy = - pi8 *TEyyy
      TEyyz = - pi8 *TEyyz
      TEzyz = - pi8 *TEzyz

      TExzx = - pi8 *TExzx
      TExzy = - pi8 *TExzy
      TExzz = - pi8 *TExzz
      TEyzy = - pi8 *TEyzy
      TEyzz = - pi8 *TEyzz
      TEzzz = - pi8 *TEzzz

c---------
c GF
c---------
      GFxx = 0.0D0
      GFyx = 0.0D0
      GFzx = 0.0D0

      GFxy = - pi8 *(dy *ri /(r + dx))
      GFyy = 0.0D0
      GFzy = 0.0D0

      GFxz = - pi8 *(dz *ri /(r + dx))
      GFyz = 0.0D0
      GFzz = 0.0D0

! compute dG[ij]/dx[k]
      TF0xxx = 0.0D0
      TF0xxy = 0.0D0
      TF0xxz = 0.0D0
      TF0yxy = 0.0D0
      TF0yxz = 0.0D0
      TF0zxz = 0.0D0

      TF0xyx = - dy *ri3
      TF0xyy = ri *r1i - dyy *ri3 *r1i - dyy *ri *ri *r1i *r1i
      TF0xyz = - dyz *ri3 *r1i - dyz *ri *ri *r1i *r1i
      TF0yyy = 0.0D0
      TF0yyz = 0.0D0
      TF0zyz = 0.0D0

      TF0xzx = - dz *ri3
      TF0xzy = - dyz *ri3 *r1i - dyz *ri *ri *r1i *r1i
      TF0xzz = ri *r1i - dzz *ri3 *r1i - dzz *ri *ri *r1i *r1i
      TF0yzy = 0.0D0
      TF0yzz = 0.0D0
      TF0zzz = 0.0D0
      
!     compute +dG[kj]/dx[i]
      TFxxx = TF0xxx + TF0xxx
      TFxxy = TF0xxy
      TFxxz = TF0xxz
      TFyxy = TF0yxy + TF0yxy
      TFyxz = TF0yxz
      TFzxz = TF0zxz + TF0zxz

      TFxyx = TF0xyx + TF0xyx
      TFxyy = TF0xyy
      TFxyz = TF0xyz
      TFyyy = TF0yyy + TF0yyy
      TFyyz = TF0yyz
      TFzyz = TF0zyz + TF0zyz

      TFxzx = TF0xzx + TF0xzx
      TFxzy = TF0xzy
      TFxzz = TF0xzz
      TFyzy = TF0yzy + TF0yzy
      TFyzz = TF0yzz
      TFzzz = TF0zzz + TF0zzz

!     *-pi8
      TFxxx = - pi8 *TFxxx
      TFxxy = - pi8 *TFxxy
      TFxxz = - pi8 *TFxxz
      TFyxy = - pi8 *TFyxy
      TFyxz = - pi8 *TFyxz
      TFzxz = - pi8 *TFzxz

      TFxyx = - pi8 *TFxyx
      TFxyy = - pi8 *TFxyy
      TFxyz = - pi8 *TFxyz
      TFyyy = - pi8 *TFyyy
      TFyyz = - pi8 *TFyyz
      TFzyz = - pi8 *TFzyz

      TFxzx = - pi8 *TFxzx
      TFxzy = - pi8 *TFxzy
      TFxzz = - pi8 *TFxzz
      TFyzy = - pi8 *TFyzy
      TFyzz = - pi8 *TFyzz
      TFzzz = - pi8 *TFzzz      
c---------
c GG
c---------
      GGxx = 0.0D0
      GGyx = 0.0D0
      GGzx = 0.0D0

      GGxy = 0.0D0
      GGyy = - pi8 *(ri /(2 * pi))
      GGzy = 0.0D0

      GGxz = 0.0D0
      GGyz = 0.0D0
      GGzz = - pi8 *(ri /(2 * pi))

!     compute dG[ij]/dx[k]
      TG0xxx = 0.0D0
      TG0xxy = 0.0D0
      TG0xxz = 0.0D0
      TG0yxy = 0.0D0
      TG0yxz = 0.0D0
      TG0zxz = 0.0D0

      TG0xyx = 0.0D0
      TG0xyy = 0.0D0
      TG0xyz = 0.0D0
      TG0yyy = 4.0D0 *dy *ri3
      TG0yyz = 4.0D0 *dz *ri3
      TG0zyz = 0.0D0

      TG0xzx = 0.0D0
      TG0xzy = 0.0D0
      TG0xzz = 0.0D0
      TG0yzy = 0.0D0
      TG0yzz = 0.0D0
      TG0zzz = 4.0D0 *dz *ri3

!     compute dG[kj]/dx[i]
      TGxxx = TG0xxx + TG0xxx
      TGxxy = TG0xxy
      TGxxz = TG0xxz
      TGyxy = TG0yxy
      TGyxz = TG0yxz
      TGzxz = TG0zxz

      TGxyx = TG0xyx + TG0xyx
      TGxyy = TG0xyy + 4.0D0 *dx *ri3
      TGxyz = TG0xyz
      TGyyy = TG0yyy + 4.0D0 *dy *ri3
      TGyyz = TG0yyz
      TGzyz = TG0zyz

      TGxzx = TG0xzx + TG0xzz
      TGxzy = TG0xzy
      TGxzz = TG0xzz + 4.0D0 *dx *ri3
      TGyzy = TG0yzy + TG0yzy
      TGyzz = TG0yzz + 4.0D0 *dy *ri3
      TGzzz = TG0zzz + 4.0D0 *dz *ri3
      
c---------
c assemble
c---------

      h0   = x0-wall
      h02  = 2.0D0*h0
      h0s2 = 2.0D0*h0*h0

      h0   = x0-wall
      h04  = 4.0D0*h0
      h0s4 = 4.0D0*h0*h0


      Gxx = Sxx + h0s4 * PDxx - h04 * SDxx + GExx + GFxx + GGxx 
      Gxy = Sxy + h0s4 * PDxy - h04 * SDxy + GExy + GFxy + GGxy
      Gxz = Sxz + h0s4 * PDxz - h04 * SDxz + GExz + GFxz + GGxz 

      Gyx = Syx + h0s4 * PDyx - h04 * SDyx + GEyx + GFyx + GGyx
      Gyy = Syy + h0s4 * PDyy - h04 * SDyy + GEyy + GFyy + GGyy
      Gyz = Syz + h0s4 * PDyz - h04 * SDyz + GEyz + GFyz + GGyz
     

      Gzx = Szx + h0s4 * PDzx - h04 * SDzx + GEzx + GFzx + GGzx
      Gzy = Szy + h0s4 * PDzy - h04 * SDzy + GEzy + GFzy + GFzy
      Gzz = Szz + h0s4 * PDzz - h04 * SDzz + GEzz + GFzz + GFzz

      

c-----------------------
      if(Iopt.ne.1) then


      Hx = h02 *(- ri3 + 3.0D0 +dxx *ri5)
      Hy = - h02 *3.0D0 *dxy *ri5
      Hz = - h02 *3.0D0 *dxz *ri5

      Txxx = TSxxx + h0s4 *TPDxxx - h04 *TSDxxx + TExxx + TFxxx + TGxxx
     +       - Hx 
      Txxy = TSxxy + h0s4 *TPDxxy - h04 *TSDxxy + TExxy + TFxxy + TGxxy
      Txxz = TSxxz + h0s4 *TPDxxz - h04 *TSDxxz + TExxz + TFxxz + TGxxz
      Tyxy = TSyxy + h0s4 *TPDyxy - h04 *TSDyxy + TEyxy + TFyxy + TGyxy
     +       - Hx
      Tyxz = TSyxz + h0s4 *TPDyxz - h04 *TSDyxz + TEyxz + TFyxz + TGyxz
      Tzxz = TSzxz + h0s4 *TPDzxz - h04 *TSDzxz + TEzxz + TFzxz + TGzxz
     +       - Hx


      
      Txyx = TSxyx + h0s4 *TPDxyx - h04 *TSDxyx + TExyx + TFxyx + TGxyx
     +       - Hy
      Txyy = TSxyy + h0s4 *TPDxyy - h04 *TSDxyy + TExyy + TFxyy + TGxyy
      Txyz = TSxyz + h0s4 *TPDxyz - h04 *TSDxyz + TExyz + TFxyz + TGxyz
      Tyyy = TSyyy + h0s4 *TPDyyy - h04 *TSDyyy + TEyyy + TFyyy + TGyyy
     +       - Hy
      Tyyz = TSyyz + h0s4 *TPDyyz - h04 *TSDyyz + TEyyz + TFyyz + TGyyz
      Tzyz = TSzyz + h0s4 *TPDzyz - h04 *TSDzyz + TEzyz + TFzyz + TGzyz
     +       - Hy


      
      Txzx = TSxzx + h0s4 *TPDxzx - h04 *TSDxzx + TExzx + TFxzx + TGxzx
     +       - Hz
      Txzy = TSxzy + h0s4 *TPDxzy - h04 *TSDxzy + TExzy + TFxzy + TGxzy
      Txzz = TSxzz + h0s4 *TPDxzz - h04 *TSDxzz + TExzz + TFxzz + TGxzz
      Tyzy = TSyzy + h0s4 *TPDyzy - h04 *TSDyzy + TEyzy + TFyzy + TGyzy
     +       - Hz
      Tyzz = TSyzz + h0s4 *TPDyzz - h04 *TSDyzz + TEyzz + TFyzz + TGyzz
      Tzzz = TSzzz + h0s4 *TPDzzz - h04 *TSDzzz + TEzzz + TFzzz + TGzzz
     +       - Hz


      
      px = psx - h02 * PSDx
      py = psy - h02 * PSDy
      pz = psz - h02 * PSDz

      end if
c----------------------
c$$$      write (6,*)
c$$$      write (6,*) Gxx,Gxy,Gxz
c$$$      write (6,*) Gyx,Gyy,Gyz
c$$$      write (6,*) Gzx,Gzy,Gzz
c$$$      write (6,*) GEyy,GFyy,GGyy,ri,ri3,dyy

c$$$      write (6,*)
c$$$      write (6,*) Txxx,Txxy,Txxz
c$$$      write (6,*) Tyxy,Tyxz,Tzxz
c-----
c done
c-----

      return
      end
