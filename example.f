
      program example

      implicit double precision(a-h,o-z)
      dimension a(3,3)
      dimension ipiv(3)
      dimension b(3)

      a(1,1) = 1.0d0
      a(1,2) = 2.0d0
      a(1,3) = 0.0d0
      a(2,1) = 3.0d0
      a(2,2) = 2.89372d0
      a(2,3) = 0.0d0
      a(3,1) = 0.0d0
      a(3,2) = 0.0d0
      a(3,3) = 0.0d0

      b(1) = 1.0d0
      b(2) = 1.3d0
      b(3) = 0.0d0
      
      call dgesv (3,1,a,3,ipiv,3,3,info)

      if(info.eq.0) then
         write (6,*) b(1),b(2),b(3)
      else
         write (6,*) 'failed',info
      end if
      
      end program
      
