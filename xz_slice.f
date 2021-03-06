      subroutine xz_slice 
     +
     +  (npts
     +  ,cx,cy,cz
     +  ,nxz,jxz
     +  )

c========================================
c FDLIB, BEMLIB
c
c Copyright by C. Pozrikidis, 1999
c All rights reserved
c
c This program is to be used only under the
c stipulations of the licensing agreement
c========================================

c---------------------------
c identify nodes with z = cz
c---------------------------

      Implicit Double Precision (a-h,o-z)

      Dimension  p(1026,3)
      Dimension ne(1026,7)

      Dimension n(512,6),nbe(512,3)
      Dimension jxz(100)

      Logical switch

      common/points/p,ne
      common/elmnts/n,nbe

c---
c introduce an integer array jxz that contains
c the labels of all points in the z = cz plane
c---

      nxz = 0   ! counter

      Do i=1,npts
        if(abs(p(i,2)-cy).lt.0.000001) then
          nxz = nxz+1
          jxz(nxz) = i
        end if
      End Do

c---
c sort marker points
c according to meridional angle
c using bubble-sort
c---

      switch = .true.

      Do while (switch)

        switch = .false.

        Do i=1,nxz-1
          i1 = i+1
          phii  = atan2(p(jxz(i ),1)-cx,
     +                  p(jxz(i ),3)-cz)
          phip1 = atan2(p(jxz(i1),1)-cx,
     +                  p(jxz(i1),3)-cz)

          if(phii.gt.phip1) then
             jhold = jxz(i)
             jxz(i)  = jxz(i1)
             jxz(i1) = jhold
             switch = .true.
            end if
         End Do

      End Do

c---
c close loop
c---

      nxz = nxz+1
      jxz(nxz) = jxz(1)

c-----
c done
c-----

      return
      end
