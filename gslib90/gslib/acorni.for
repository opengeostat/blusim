      double precision function acorni(idum)
c-----------------------------------------------------------------------
c
c Fortran implementation of ACORN random number generator of order less
c than or equal to 12 (higher orders can be obtained by increasing the
c parameter value MAXORD).
c
c
c NOTES: 1. The variable idum is a dummy variable. The common block
c           IACO is used to transfer data into the function.
c
c        2. Before the first call to ACORN the common block IACO must
c           be initialised by the user, as follows. The values of
c           variables in the common block must not subsequently be
c           changed by the user.
c
c             KORDEI - order of generator required ( must be =< MAXORD)
c
c             MAXINT - modulus for generator, must be chosen small
c                      enough that 2*MAXINT does not overflow
c
c             ixv(1) - seed for random number generator
c                      require 0 < ixv(1) < MAXINT
c
c             (ixv(I+1),I=1,KORDEI)
c                    - KORDEI initial values for generator
c                      require 0 =< ixv(I+1) < MAXINT
c
c        3. After initialisation, each call to ACORN generates a single
c           random number between 0 and 1.
c
c        4. An example of suitable values for parameters is
c
c             KORDEI   = 10
c             MAXINT   = 2**30
c             ixv(1)   = an odd integer in the (approximate) range 
c                        (0.001 * MAXINT) to (0.999 * MAXINT)
c             ixv(I+1) = 0, I=1,KORDEI
c
c
c
c Author: R.S.Wikramaratna,                           Date: October 1990
c-----------------------------------------------------------------------
c      implicit double precision (a-h,o-z)
      integer KORDEI, MAXOP1, MAXINT,idum,i
      parameter (KORDEI=12,MAXOP1=KORDEI+1,MAXINT=2**30)
      integer ixv(MAXOP1)
      common/iaco/ ixv
c      common/iaco/ ixv(MAXOP1)

      do i=1,KORDEI
            ixv(i+1)=(ixv(i+1)+ixv(i))
            if(ixv(i+1).ge.MAXINT) ixv(i+1)=ixv(i+1)-MAXINT
      end do
      acorni=dble(ixv(KORDEI+1))/MAXINT
      return
      end




            double precision function acorniopt(idum)

c      implicit double precision (a-h,o-z)
      integer KORDEI, MAXOP1, MAXINT,idum,i,tmpixv
      parameter (KORDEI=12,MAXOP1=KORDEI+1,MAXINT=2**30)
      integer ixv(MAXOP1)
      common/iaco/ ixv

c      common/iaco/ ixv(MAXOP1)

      tmpixv=(ixv(2 )+ixv(1 ))
      ixv(2 )=tmpixv-shiftl(shiftr(MAXINT-tmpixv,31),30)
      tmpixv=(ixv(3 )+ixv(2 ))
      ixv(3 )=tmpixv-shiftl(shiftr(MAXINT-tmpixv,31),30)
      tmpixv=(ixv(4 )+ixv(3 ))
      ixv(4 )=tmpixv-shiftl(shiftr(MAXINT-tmpixv,31),30)
      tmpixv=(ixv(5 )+ixv(4 ))
      ixv(5 )=tmpixv-shiftl(shiftr(MAXINT-tmpixv,31),30)
      tmpixv=(ixv(6 )+ixv(5 ))
      ixv(6 )=tmpixv-shiftl(shiftr(MAXINT-tmpixv,31),30)
      tmpixv=(ixv(7 )+ixv(6 ))
      ixv(7 )=tmpixv-shiftl(shiftr(MAXINT-tmpixv,31),30)
      tmpixv=(ixv(8 )+ixv(7 ))
      ixv(8 )=tmpixv-shiftl(shiftr(MAXINT-tmpixv,31),30)
      tmpixv=(ixv(9 )+ixv(8 ))
      ixv(9 )=tmpixv-shiftl(shiftr(MAXINT-tmpixv,31),30)
      tmpixv=(ixv(10)+ixv(9 ))
      ixv(10)=tmpixv-shiftl(shiftr(MAXINT-tmpixv,31),30)
      tmpixv=(ixv(11)+ixv(10))
      ixv(11)=tmpixv-shiftl(shiftr(MAXINT-tmpixv,31),30)
      tmpixv=(ixv(12)+ixv(11))
      ixv(12)=tmpixv-shiftl(shiftr(MAXINT-tmpixv,31),30)
      tmpixv=(ixv(13)+ixv(12))
      ixv(13)=tmpixv-shiftl(shiftr(MAXINT-tmpixv,31),30)

      acorniopt=dble(ixv(KORDEI+1))/MAXINT
      return
      end
