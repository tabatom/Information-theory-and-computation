      PROGRAM Ex2_b
      
      real*4 pi4,sq2_4,res_4
      real*8 pi8,sq2_8,res_8
      
      pi4 = ACOS(-1.) * 10e32
      pi8 = ACOS(-1.) * 10e32
      
      sq2_4 = SQRT(2.) * 10e21
      sq2_8 = SQRT(2.) * 10e21
      
      res_4 = pi4+sq2_4
      res_8 = pi8+sq2_8
      
      print*,res_4
      print*,res_8
      
      STOP
      END PROGRAM Ex2_b
