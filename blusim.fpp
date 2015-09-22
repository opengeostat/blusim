      program main
c-----------------------------------------------------------------------
c
c             Block LUSIM of a 3-D Rectangular Grid
c             *************************************
c
c
c AUTHOR: Julian M. Ortiz                      DATE: Dec. 2006
c OpenMP add-ons: Oscar Peredo                 DATE: Nov. 2013
c-----------------------------------------------------------------------
c 
c - msflib was commented and the dynamic arrays were defined inside of 
c   the main routine (also commented, on top of the code)
c
c-----------------------------------------------------------------------
c
c      use msflib
c      use geostat

#ifdef _OPENMP
      use omp_lib
#endif

#ifdef TRACE
      use extrae_module
#endif

      include  'blusim.inc'
      parameter (MV=20,KORDEI=12,MAXOP1=KORDEI+1,MAXINT=2**30,MAXTHR=10)
#ifdef _OPENMP
      parameter (MAXTHREADS=16)
#endif
      real      var(MV),ltpar,utpar,dtest,thres(MAXTHR),
     +mean(MAXTHR),prop(MAXTHR),pavg,pvar,pgvar
      character datafl*512,transfl*512,smthfl*512,outfl*512,dbgfl*512,
     +          str*512,tmpfl*512
#ifdef _OPENMP
      character outflThreads(43,MAXTHREADS) , outfltmp*43
#endif
      integer ltail,utail,tipx,tipy,tipz,ibx,iby,ibz,tibx,tiby,tibz,
     +nthr
      logical   testfl,coloc
      real*8  p,cp,w, rotmat(MAXROT,3,3)
      real*4 timeIni, timeFin
      integer indpn,iz,iy,ix,sumsim,ic,nclose,infoct,i,j,k,na,nu
      real xloc,yloc,zloc

#ifdef _OPENMP
      integer threadId,numThreads,ithread
      integer loutThreads(MAXTHREADS)
#endif

      integer,allocatable :: nisb(:),ixsbtosr(:),iysbtosr(:),izsbtosr(:)
      real,allocatable    :: vrtr(:),vrgtr(:),x(:),y(:),z(:),vr(:),
     +wt(:),tmp(:),closeArray(:),xa(:),ya(:),za(:),
     +vra(:),xdb(:),ydb(:),zdb(:),c11(:,:),
     +c12(:,:),c22(:,:),u12(:,:),c22i(:,:),
     +l11(:,:),l11inv(:,:),l21v1(:),l22(:,:),
     +le(:),slsk(:),v1(:),w2(:),y2(:),simsmu(:),
     +xg(:),yg(:),zg(:),simbl(:),sdbg(:),pdbg(:),
     +sumbl(:),sqdbl(:)

#ifdef OPTI
      integer,allocatable :: xyzdb(:,:)
#endif

c
c For random number generator:
c
      real*8  acorni
      common /iaco/ ixv(MAXOP1)
      data   ixv/MAXOP1*0.0/


#ifdef TRACE
      call extrae_init()
#endif

#ifdef _OPENMP
c$omp parallel
      numThreads = OMP_get_num_threads()
      print *,'numThreads=',numThreads
c$omp end parallel
#endif

      do i=1,50000
          p=acorni(idum)
      end do
c
c FORTRAN Units:
c
      lin   = 1
      ldbg  = 3
      lout  = 4
#ifdef _OPENMP
      do i=1,MAXTHREADS
            loutThreads(i)=lout*100+i
      end do
#endif
      lext  = 7
      ljack = 8
      ldbg2 = 5
      ldbg3 = 6

c
c READ PARAMETERS
c
c
c Note VERSION number:
c
      write(*,9999) VERSION
 9999 format(/' BLUSIM Version: ',f5.3/)
c
c Get the name of the parameter file - try the default name if no input:
c
      do i=1,512
            str(i:i) = ' '
      end do
      call getarg(1,str)
      if(str(1:1).eq.' ')then
            write(*,*) 'Which parameter file do you want to use?'
            read (*,'(a)') str
      end if
      if(str(1:1).eq.' ') str(1:20) = 'blusim.par            '
      inquire(file=str,exist=testfl)
      if(.not.testfl) then
            write(*,*) 'ERROR - the parameter file does not exist,'
            write(*,*) '        check for the file and try again  '
            write(*,*)
            if(str(1:20).eq.'blusim.par            ') then
                  write(*,*) '        creating a blank parameter file'
                  call makepar
                  write(*,*)
            end if
            stop
      endif
      open(lin,file=str,status='OLD')
c
c Find Start of Parameters:
c
 1    read(lin,'(a4)',end=98) str(1:4)
      if(str(1:4).ne.'STAR') go to 1
      write(*,*) ' *** 1 *** Reading parameters'
c
c Read Input Parameters:
c
      read(lin,'(a512)',err=98) datafl
      call chknam(datafl,512)
      write(*,*) ' data file = ',datafl(1:40)

      read(lin,*,err=98) ixl,iyl,izl,ivrl,iwtl
      write(*,*) ' columns = ',ixl,iyl,izl,ivrl,iwtl

      read(lin,*,err=98) tmin,tmax
      write(*,*) ' trimming limits = ',tmin,tmax

      read(lin,*,err=98) itrans
      write(*,*) ' transformation flag = ',itrans

      read(lin,'(a512)',err=98) transfl
      call chknam(transfl,512)
      write(*,*) ' transformation file = ',transfl(1:40)

      read(lin,*,err=98) ismooth
      write(*,*) ' consider smoothed distribution (1=yes) = ',ismooth

      read(lin,'(a512)',err=98) smthfl
      call chknam(smthfl,512)
      write(*,*) ' file with smoothed distribution = ',smthfl(1:40)

      read(lin,*,err=98) isvr,iswt
      write(*,*) ' columns = ',isvr,iswt

      read(lin,*,err=98) zmin,zmax
      write(*,*) ' data limits (tails) = ',zmin,zmax

      read(lin,*,err=98) ltail,ltpar
      write(*,*) ' lower tail = ',ltail,ltpar

      read(lin,*,err=98) utail,utpar
      write(*,*) ' upper tail = ',utail,utpar

      read(lin,*,err=98) idbg
      write(*,*) ' debugging level = ',idbg

      read(lin,'(a512)',err=98) dbgfl
      call chknam(dbgfl,512)
      write(*,*) ' debugging file = ',dbgfl(1:40)
      open(ldbg,file=dbgfl,status='UNKNOWN')

      read(lin,'(a512)',err=98) outfl
      call chknam(outfl,512)
      write(*,*) ' output file = ',outfl(1:40)

#ifdef _OPENMP
      do i=1,numThreads
            write(outfltmp,"(A40,I3)") trim(outfl(1:40)),
     + loutThreads(i)
            outfltmp = adjustl(trim(outfltmp))
            write(*,*) ' output file (threads) = ',
     + outfltmp
      end do
#endif

      read(lin,*,err=98) nsim
      write(*,*) ' number of realizations = ',nsim

      read(lin,*,err=98) nx,xmn,xsiz
      write(*,*) ' Pannel: nx, xmn, xsiz = ',nx,xmn,xsiz

      read(lin,*,err=98) ny,ymn,ysiz
      write(*,*) '         ny, ymn, ysiz = ',ny,ymn,ysiz

      read(lin,*,err=98) nz,zmn,zsiz
      write(*,*) '         nz, zmn, zsiz = ',nz,zmn,zsiz

      read(lin,*,err=98) nbx,nby,nbz
      write(*,*) ' SMUs per Pannel:',nbx,nby,nbz

      read(lin,*,err=98) nxdis,nydis,nzdis
      write(*,*) ' block discretization:',nxdis,nydis,nzdis

      read(lin,*,err=98) nthr
      write(*,*) ' number of cutoffs for reporting:',nthr
      if(nthr.gt.MAXTHR) stop ' nthr too large!'

      read(lin,*,err=98) (thres(j),j=1,nthr)
      write(*,*) ' cutoffs:',(thres(j),j=1,nthr)

      read(lin,*,err=98) ixv(1)
      write(*,*) ' random number seed = ',ixv(1)

      read(lin,*,err=98) ndmin,ndmax
      write(*,*) ' ndmin,ndmax = ',ndmin,ndmax

      read(lin,*,err=98) noct
      write(*,*) ' max per octant = ',noct

      read(lin,*,err=98) radius,radius1,radius2
      write(*,*) ' search radii = ',radius,radius1,radius2
      if(radius.lt.EPSLON) stop 'radius must be greater than zero'
      radsqd = radius  * radius
      sanis1 = radius1 / radius
      sanis2 = radius2 / radius

      read(lin,*,err=98) sang1,sang2,sang3
      write(*,*) ' search anisotropy angles = ',sang1,sang2,sang3

      read(lin,*,err=98) nst(1),c0(1)
      write(*,*) ' nst, c0 = ',nst(1),c0(1)

      if(nst(1).le.0) then
            write(*,9997) nst(1)
 9997       format(' nst must be at least 1, it has been set to ',i4,/,
     +             ' The c or a values can be set to zero')
            stop
      endif

      do i=1,nst(1)
            read(lin,*,err=98) it(i),cc(i),ang1(i),ang2(i),ang3(i)
            read(lin,*,err=98) aa(i),aa1,aa2
            anis1(i) = aa1 / max(aa(i),EPSLON)
            anis2(i) = aa2 / max(aa(i),EPSLON)
            write(*,*) ' it,cc,ang[1,2,3]; ',it(i),cc(i),
     +                   ang1(i),ang2(i),ang3(i)
            write(*,*) ' a1 a2 a3: ',aa(i),aa1,aa2
            if(it(i).eq.4) then
                  if(aa(i).lt.0.0) stop ' INVALID power variogram'
                  if(aa(i).gt.2.0) stop ' INVALID power variogram'
            end if
      end do

      close(lin)
c
c Done reading parameters 
c
c
c Find the needed parameters:
c
      MAXSBX = 1
      if(nx.gt.1)then
            MAXSBX = int(nx/2.00)
            if(MAXSBX.gt.50)MAXSBX=50
      end if
c
      MAXSBY = 1
      if(ny.gt.1)then
            MAXSBY = int(ny/2.00)
            if(MAXSBY.gt.50)MAXSBY=50
      end if
c
      MAXSBZ = 1
      if(nz.gt.1)then
            MAXSBZ = int(nz/2.00)
            if(MAXSBZ.gt.50)MAXSBZ=50
      end if
c
      MAXSB = MAXSBX*MAXSBY*MAXSBZ
      MXSXY = 4 * MAXSBX * MAXSBY
      MXSX  = 2 * MAXSBX
      MAXSAM = ndmax
      MAXDIS = nxdis*nydis*nzdis
c
c TRANFORMATION TABLE
c
c
c Find MAXDAT:
c
      inquire(file=datafl,exist=testfl)
      if(testfl)then
            open(lin,file=datafl,status='UNKNOWN')
            read(lin,*,err=99)
            read(lin,*,err=99) nvari
            do i=1,nvari
                  read(lin,*)
            end do
            MAXDAT = 0
 33         read(lin,*,end=66,err=99)(var(j),j=1,nvari)
            MAXDAT = MAXDAT + 1
            go to 33
 66         continue
            rewind(lin)
            close(lin)
      else
            stop ' MISSING DATA FILE'
      end if
      write(*,*) 'MAXDAT',MAXDAT
c
c Decide which file to use for establishing the transformation table:
c
      if(ismooth.eq.1) then
            tmpfl  = smthfl
            icolvr = isvr
            icolwt = iswt
      else
            tmpfl  = datafl
            icolvr = ivrl
            icolwt = iwtl
      end if
      inquire(file=tmpfl,exist=testfl)
      if(.not.testfl) stop 'ERROR: Transformation'
c
c Find MAXTMP:
c
      inquire(file=tmpfl,exist=testfl)
      if(testfl)then
            open(lin,file=tmpfl,status='UNKNOWN')
            read(lin,*,err=99)
            read(lin,*,err=99) nvari
            do i=1,nvari
                  read(lin,*)
            end do
            MAXTMP = 0
 34         read(lin,*,end=67,err=99)(var(j),j=1,nvari)
            MAXTMP = MAXTMP + 1
            go to 34
 67         continue
            rewind(lin)
            close(lin)
      else
            stop ' MISSING TRANSF FILE'
      end if
      write(*,*) 'MAXTMP',MAXTMP
c
c Some parameters
c
      nxy   = nx*ny
      nxyz  = nx*ny*nz
      nsmup = nbx*nby*nbz
      irepo = max(1,min((nxyz/10),10000))
c
c Allocate memory
c
      allocate(vrtr(MAXTMP),stat = test)
            if(test.ne.0)stop ' ALLOCATION ERROR'
      allocate(vrgtr(MAXTMP),stat = test)
            if(test.ne.0)stop ' ALLOCATION ERROR'
      allocate(x(MAXDAT),stat = test)
            if(test.ne.0)stop ' ALLOCATION ERROR'
      allocate(y(MAXDAT),stat = test)
            if(test.ne.0)stop ' ALLOCATION ERROR'
      allocate(z(MAXDAT),stat = test)
            if(test.ne.0)stop ' ALLOCATION ERROR'
      allocate(vr(MAXDAT),stat = test)
            if(test.ne.0)stop ' ALLOCATION ERROR'
      allocate(wt(MAXDAT),stat = test)
            if(test.ne.0)stop ' ALLOCATION ERROR'
      allocate(tmp(MAXDAT),stat = test)
            if(test.ne.0)stop ' ALLOCATION ERROR'
      allocate(nisb(MAXSB),stat = test)
            if(test.ne.0)stop ' ALLOCATION ERROR'
      allocate(ixsbtosr(8 * MAXSB),stat = test)
            if(test.ne.0)stop ' ALLOCATION ERROR'
      allocate(iysbtosr(8 * MAXSB),stat = test)
            if(test.ne.0)stop ' ALLOCATION ERROR'
      allocate(izsbtosr(8 * MAXSB),stat = test)
            if(test.ne.0)stop ' ALLOCATION ERROR'
      allocate(closeArray(MAXDAT),stat = test)
            if(test.ne.0)stop ' ALLOCATION ERROR'
      allocate(xa(MAXSAM),stat = test)
            if(test.ne.0)stop ' ALLOCATION ERROR'
      allocate(ya(MAXSAM),stat = test)
            if(test.ne.0)stop ' ALLOCATION ERROR'
      allocate(za(MAXSAM),stat = test)
            if(test.ne.0)stop ' ALLOCATION ERROR'
      allocate(vra(MAXSAM),stat = test)
            if(test.ne.0)stop ' ALLOCATION ERROR'
#ifdef OPTI
      allocate(xyzdb(3,MAXDIS),stat = test)
            if(test.ne.0)stop ' ALLOCATION ERROR'
#else
      allocate(xdb(MAXDIS),stat = test)
            if(test.ne.0)stop ' ALLOCATION ERROR'
      allocate(ydb(MAXDIS),stat = test)
            if(test.ne.0)stop ' ALLOCATION ERROR'
      allocate(zdb(MAXDIS),stat = test)
            if(test.ne.0)stop ' ALLOCATION ERROR'
#endif
      allocate(xg(MAXDIS),stat = test)
            if(test.ne.0)stop ' ALLOCATION ERROR'
      allocate(yg(MAXDIS),stat = test)
            if(test.ne.0)stop ' ALLOCATION ERROR'
      allocate(zg(MAXDIS),stat = test)
            if(test.ne.0)stop ' ALLOCATION ERROR'
      allocate(simbl(nsim),stat = test)
            if(test.ne.0)stop ' ALLOCATION ERROR'
      allocate(sdbg(nxyz*nsmup),stat = test)
            if(test.ne.0)stop ' ALLOCATION ERROR'
      allocate(pdbg(MAXDIS*nxyz),stat = test)
            if(test.ne.0)stop ' ALLOCATION ERROR'
      allocate(simsmu(nsim),stat = test)
            if(test.ne.0)stop ' ALLOCATION ERROR'
      allocate(sumbl(nsim),stat = test)
            if(test.ne.0)stop ' ALLOCATION ERROR'
      allocate(sqdbl(nsim),stat = test)
            if(test.ne.0)stop ' ALLOCATION ERROR'
      allocate(c11(MAXSAM,MAXSAM),stat = test)
            if(test.ne.0)stop ' ALLOCATION ERROR'
      allocate(c12(MAXSAM,MAXDIS),stat = test)
            if(test.ne.0)stop ' ALLOCATION ERROR'
      allocate(c22(MAXDIS,MAXDIS),stat = test)
            if(test.ne.0)stop ' ALLOCATION ERROR'
      allocate(u12(MAXSAM,MAXDIS),stat = test)
            if(test.ne.0)stop ' ALLOCATION ERROR'
      allocate(c22i(MAXDIS,MAXDIS),stat = test)
            if(test.ne.0)stop ' ALLOCATION ERROR'
      allocate(l11(MAXSAM,MAXSAM),stat = test)
            if(test.ne.0)stop ' ALLOCATION ERROR'
      allocate(l11inv(MAXSAM,MAXSAM),stat = test)
            if(test.ne.0)stop ' ALLOCATION ERROR'
      allocate(l21v1(MAXDIS),stat = test)
            if(test.ne.0)stop ' ALLOCATION ERROR'
      allocate(l22(MAXDIS,MAXDIS),stat = test)
            if(test.ne.0)stop ' ALLOCATION ERROR'
      allocate(le(MAXSAM),stat = test)
            if(test.ne.0)stop ' ALLOCATION ERROR'
      allocate(slsk(MAXDIS),stat = test)
            if(test.ne.0)stop ' ALLOCATION ERROR'
      allocate(v1(MAXSAM),stat = test)
            if(test.ne.0)stop ' ALLOCATION ERROR'
      allocate(w2(MAXDIS),stat = test)
            if(test.ne.0)stop ' ALLOCATION ERROR'
      allocate(y2(MAXDIS),stat = test)
            if(test.ne.0)stop ' ALLOCATION ERROR'
c
c Check to make sure the data file exists:
c
      nd = 0
      av = 0.0
      ss = 0.0
      gav = 0.0
      gss = 0.0
c
c Establish the reference histogram for the simulation (provided that
c we have data, and we are transforming the data):
c
      write(*,*) ' *** 2 *** Transformation Table'
      if(itrans.eq.1) then
            write(*,*) 'Setting up transformation table'
c
c Open up the file with reference distribution:
c
            open(lin,file=tmpfl,status='UNKNOWN')
            read(lin,'(a40)',err=99) str(1:40)
            read(lin,*,err=99) nvari
            do i=1,nvari
                  read(lin,*,err=98)
            end do
c
c Now, read in the actual data:
c
            nt     = 0
            ntr    = 0
            twt    = 0.0
 3          read(lin,*,end=4,err=99) (var(j),j=1,nvari)
c
c Trim this data?
c
            if(var(icolvr).lt.tmin.or.var(icolvr).ge.tmax) then
                  nt = nt + 1
                  go to 3
            endif
            ntr = ntr + 1
c
c Error in column number?
c
            if(icolvr.gt.nvari.or.icolwt.gt.nvari) then
                  write(*,*) ' ERROR: too few columns in ref data '
                  stop
            endif
c
c Keep this data: Assign the data value:
c
            vrtr(ntr) = var(icolvr)
            if(icolwt.le.0) then
                  vrgtr(ntr) = 1.0
            else
                  vrgtr(ntr) = var(icolwt)
            endif
            if(vrgtr(ntr).le.0.0) then
                  ntr = ntr - 1
                  nt  = nt  + 1
                  go to 3
            end if
            twt = twt + vrgtr(ntr)
c
c Go back for another datum:
c
            go to 3
 4          close(lin)
            if(ntr.le.1) then
                  write(*,*) 'ERROR: too few data for transformation'
                  stop
            endif
c
c Write transformation table:
c
            open(lout,file=transfl,status='UNKNOWN')

c
c Sort data by value:
c
            istart = 1
            iend   = ntr
            call sortem(istart,iend,vrtr,1,vrgtr,c,d,e,f,g,h)
c
c Compute the cumulative probabilities and write transformation table
c
            twt   = max(twt,EPSLON)
            oldcp = 0.0
            cp    = 0.0
            do j=istart,iend
                  cp =  cp + dble(vrgtr(j)/twt)
                  w  = (cp + oldcp)*0.5
                  call gauinv(w,vrg,ierr)
                  if(ierr.eq.1) vrg = UNEST
                  write(lout,201) vrtr(j),vrg
 201              format(f12.5,1x,f12.5)
                  oldcp =  cp
c
c Now, reset the weight to the normal scores value:
c
                  vrgtr(j) = vrg
            end do
            close(lout)
      end if
c
c READ AND TRANSFORM DATA:
c
      write(*,*) ' *** 3 *** Transforming Data'
      inquire(file=datafl,exist=testfl)
      if(testfl) then
            write(*,*) 'Reading input data'
            open(lin,file=datafl,status='OLD')
            read(lin,*,err=99)
            read(lin,*,err=99) nvari
            do i=1,nvari
                  read(lin,*,err=99)
            end do
            if(ixl.gt.nvari.or.iyl.gt.nvari.or.izl.gt.nvari.or.
     +         ivrl.gt.nvari.or.iwtl.gt.nvari) then
                  write(*,*) 'ERROR: you have asked for a column number'
                  write(*,*) '       greater than available in file'
                  stop
            end if
c
c Read all the data until the end of the file:
c
            twt = 0.0
            nd  = 0
            nt  = 0
 5          read(lin,*,end=6,err=99) (var(j),j=1,nvari)
            if(var(ivrl).lt.tmin.or.var(ivrl).ge.tmax) then
                  nt = nt + 1
                  go to 5
            end if
            nd = nd + 1
c
c Acceptable data, assign the value, X, Y, Z coordinates, and weight:
c
            vr(nd) = var(ivrl)
            if(ixl.le.0) then
                  x(nd) = xmn
            else
                  x(nd) = var(ixl)
            endif
            if(iyl.le.0) then
                  y(nd) = ymn
            else
                  y(nd) = var(iyl)
            endif
            if(izl.le.0) then
                  z(nd) = zmn
            else
                  z(nd) = var(izl)
            endif
            if(iwtl.le.0) then
                  wt(nd) = 1.0
            else
                  wt(nd) = var(iwtl)
            endif
c
c Normal scores transform?
c
            if(itrans.eq.1) then
                  vrr = vr(nd)
                  call locate(vrtr,ntr,1,ntr,vrr,j)
                  j   = min(max(1,j),(ntr-1))
                  vrg = powint(vrtr(j),vrtr(j+1),vrgtr(j),vrgtr(j+1),
     +                         vrr,1.0)
                  if(vrg.lt.vrgtr(1)  ) vrg = vrgtr(1)
                  if(vrg.gt.vrgtr(ntr)) vrg = vrgtr(nd)
                  vr(nd) = vrg
            end if
            twt = twt + wt(nd)
            av  = av  + var(ivrl)*wt(nd)
              gav = gav + vrg*wt(nd)
            ss  = ss  + var(ivrl)*var(ivrl)*wt(nd)
              gss = gss + vrg*vrg*wt(nd)
            go to 5
 6          close(lin)
c
c Compute the averages and variances as an error check for the user:
c
            av = av / max(twt,EPSLON)
            ss =(ss / max(twt,EPSLON)) - av * av
              gav = gav / max(twt,EPSLON)
              gss = (gss /max(twt,EPSLON))- gav * gav
            write(ldbg,111) nd,nt,av,ss,gav,gss
            write(*,   111) nd,nt,av,ss,gav,gss
 111  format(/,' Data for SGSIM: Number of acceptable data  = ',i8,/,
     +         '                 Number trimmed             = ',i8,/,
     +         '                 Weighted Average           = ',f12.4,/,
     +         '                 Weighted Variance          = ',f12.4,/,
     +         '                 Weighted Gaussian Average  = ',f12.4,/,
     +         '                 Weighted Gaussian Variance = ',f12.4,/)
      endif
      open(lout,file=outfl,status='UNKNOWN')
#ifdef _OPENMP
c each thread will write in the files file.out401, file.out402, so on.
c up to MAXTHREADS
      do i=1,numThreads
            write(outfltmp,"(A40,I3)") trim(outfl(1:40)),
     + loutThreads(i)
            outfltmp = adjustl(trim(outfltmp))
            open(100*lout+i,file=outfltmp,status='UNKNOWN')
       end do
#endif

      write(lout,801)3+2*nthr
 801  format('BLUSIM Output',/,i5,/,'block average',/,'block variance'
     + ,/,'block global variance')
      do i=1,nthr
            write(lout,*) ' Mean above cutoff',thres(i)
      end do
      do i=1,nthr
            write(lout,*) ' Proportion above cutoff',thres(i)
      end do
c      do k=1,nsim
c      write(lout,*) 'Realization',k
c      end do
      open(ldbg2,file='dbgsim.out')
      write(ldbg2,*) ' Debugging file- Full Realization (last one)'
      write(ldbg2,*) ' 1'
      write(ldbg2,*) ' value - point support'
c      open(ldbg3,file='dbgsmu.out')
c      write(ldbg3,*) ' Debugging file- Full Realization 1'
c      write(ldbg3,*) ' 1'
c      write(ldbg3,*) ' value - smu support'

c
c Set up the rotation/anisotropy matrices that are needed for the
c variogram and search.  Also compute the maximum covariance for
c the rescaling factor:
c
      write(*,*) ' *** 4 *** Setting up the kriging/LUSIM'
      write(*,*) 'Setting up rotation matrices for variogram and search'
      radsqd = radius * radius
      PMX    = 999.0
      covmax = c0(1)
      do is=1,nst(1)
            call setrot(ang1(is),ang2(is),ang3(is),anis1(is),anis2(is),
     +                  is,MAXROT,rotmat)
            if(it(is).eq.4) then
                  covmax = covmax + PMX 
            else
                  covmax = covmax + cc(is)
            endif
      end do
      isrot = MAXNST + 1
      call setrot(sang1,sang2,sang3,sanis1,sanis2,isrot,MAXROT,rotmat)
      resc=1.0
c
c Set up for super block searching:
c
      write(*,*) 'Setting up super block search strategy'
      nsec = 0
      call setsupr(nx,xmn,xsiz,ny,ymn,ysiz,nz,zmn,zsiz,nd,x,y,z,
     +             vr,tmp,nsec,sec1,sec2,sec3,MAXSBX,MAXSBY,MAXSBZ,nisb,
     +             nxsup,xmnsup,xsizsup,nysup,ymnsup,ysizsup,nzsup,
     +             zmnsup,zsizsup)
      write(*,*) ' DONE WITH SETSUPR'
      call picksup(nxsup,xsizsup,nysup,ysizsup,nzsup,zsizsup,
     +             isrot,MAXROT,rotmat,radsqd,nsbtosr,ixsbtosr,
     +             iysbtosr,izsbtosr)
      write(*,*) ' DONE WITH PICKSUP'

c
c Set up the discretization points per block.  Figure out how many
c are needed, the spacing, and fill the xdb,ydb, and zdb arrays with
c the offsets relative to the block center (this only gets done once):
c
      write(*,*) ' *** 5 *** Setting up discretization points'
      if(nxdis.lt.1) nxdis = 1
      if(nydis.lt.1) nydis = 1
      if(nzdis.lt.1) nzdis = 1
      nds = MAXDIS/nsmup
      xdis = xsiz  / max(real(nxdis),1.0)
      ydis = ysiz  / max(real(nydis),1.0)
      zdis = zsiz  / max(real(nzdis),1.0)
      ii    = 0
      xloc = -0.5*(xsiz+xdis)
      do ix =1,nxdis
            xloc = xloc + xdis
            yloc = -0.5*(ysiz+ydis)
            do iy=1,nydis
                  yloc = yloc + ydis
                  zloc = -0.5*(zsiz+zdis)
                  do iz=1,nzdis
                        zloc = zloc + zdis
                        ii = ii+1
#ifdef OPTI
                        xyzdb(1,ii) = xloc
                        xyzdb(2,ii) = yloc
                        xyzdb(3,ii) = zloc
      write(*,*) '       *** ix,iy,iz,ii',ix,iy,iz,ii
      write(*,*) '         * xdb,ydb,zdb',xyzdb(1,ii),xyzdb(2,ii),
     +xyzdb(3,ii)
#else
                        xdb(ii) = xloc
                        ydb(ii) = yloc
                        zdb(ii) = zloc
      write(*,*) '       *** ix,iy,iz,ii',ix,iy,iz,ii
      write(*,*) '         * xdb,ydb,zdb',xdb(ii),ydb(ii),zdb(ii)
#endif
                  end do
            end do
      end do
c
c Ready to krige
c
      sdbg=0
      write(*,*) ' *** 6 *** Kriging'


c
c MAIN LOOP OVER ALL THE PANNELS IN THE GRID:
c

      call cpu_time(timeIni)
      mean=0.0
      prop=0.0
      sumsim=0

c$omp parallel default(shared) 
c$omp&     firstprivate(threadId,indpn,iz,iy,ix,xloc,yloc,zloc,
c$omp& 			nclose,closeArray,infoct,i,j,k,na,xa,ya,za,vra,
c$omp& 			nu,ii,coloc,dtest,
#ifdef OPTI
c$omp&                  xyzdb,
#else
c$omp&                  xdb,ydb,zdb,
#endif
c$omp&                  xg,yg,zg,
c$omp& 			cmax,sill,cov, 
c$omp& 			c22i,c11,l11,l11inv,ierr,v1,
c$omp& 			sle,le,c12,u12,c22,l22,l21v1,slsk,
c$omp& 			avgbl,simbl,isim,w2,p,xp,y2,
c$omp& 			ic,sumsim,zz,yy,xx,simval,
c$omp& 			tipz,tipy,tipx,ipx,ipy,ipz,indpt,pdbg,
c$omp& 			tibz,tiby,tibx,ibx,iby,ibz,indb,indsmu,
c$omp&                  sumbl,sqdbl,prop,mean,pavg,pvar,pgvar,
c$omp&                  MAXDIS,sum,
c$omp&                  ntr,vrtr,vrgtr,zmin,
c$omp&                  zmax,ltail,ltpar,utail,utpar,
cc cova3 variables
c$omp&                  istart,is,hsqd,h,hr,
cc sqrdist variables
c$omp&                  dx,dy,dz,sqdist,cont,
cc srchsupr variables
c$omp&                  index,inflag,isup,ixsup,iysup,izsup,
c$omp&                  nums,tmp,inoct,nt,iq,
cc sortem variables
c$omp&                  lt,ut,m,q,iring,ta,tb,tc,td,te,tf,tg,th,
c$omp&                  xb,xc,xd,xe,xf,xh)
#ifdef _OPENMP
      threadId = int(OMP_get_thread_num())
#else
      threadId = int(0)
#endif
      print *,'threadId=',threadId
      print *,'size(y2)=',size(y2)
c$omp do schedule(static) 
      do indpn=1,nxyz
c            print *,'indpn=',indpn
            if((int(indpn/irepo)*irepo).eq.indpn) 
     +                    write(*,100) int(threadId),indpn
 100            format(i2,':   currently on node ',i9)
c
c Where are we making an estimate?
c
            iz   = int((indpn-1)/nxy) + 1
            iy   = int((indpn-(iz-1)*nxy-1)/nx) + 1
            ix   = indpn - (iz-1)*nxy - (iy-1)*nx
            xloc = xmn + real(ix-1)*xsiz
            yloc = ymn + real(iy-1)*ysiz
            zloc = zmn + real(iz-1)*zsiz
c      write(*,*) '         * ix   iy   iz  ',ix,iy,iz
c      write(*,*) '         * xloc yloc zloc',xloc,yloc,zloc

c
c Find the nearest samples:
c
            call srchsupr(xloc,yloc,zloc,radsqd,isrot,MAXROT,rotmat,
     +              nsbtosr,
     +              ixsbtosr,iysbtosr,izsbtosr,noct,nd,x,y,z,tmp,
     +              nisb,nxsup,xmnsup,xsizsup,nysup,ymnsup,ysizsup,
     +              nzsup,zmnsup,zsizsup,nclose,closeArray,infoct)
c      write(*,*) '         * DONE WITH SRCHSUPR '
c
c Load the nearest data in xa,ya,za,vra:
c

            na = 0
            do i=1,nclose
                  ind    = int(closeArray(i)+0.5)
                  if(na.lt.ndmax) then
                        na = na + 1
                        xa(na)  = x(ind) - xloc
                        ya(na)  = y(ind) - yloc
                        za(na)  = z(ind) - zloc
                        vra(na) = vr(ind)
                  end if
            end do

c      write(*,*) '         * Nearest data',na
c      do i=1,na
c            write(*,*) '         * ',i,xa(i),ya(i),za(i),vra(i)
c      end do
c
c Establish positions of grid points (may be fewer than nx*ny*nz):
c
c                 nxyzu = number actually used
c
            nu = 0
            do ii=1, MAXDIS
                  coloc = .false.
                  do i=1,nd
#ifdef OPTI
                        dtest = abs(xyzdb(1,ii)-x(i)+xloc)+ 
     +                          abs(xyzdb(2,ii)-y(i)+yloc)+ 
     +                          abs(xyzdb(3,ii)-z(i)+zloc)
#else
                        dtest = abs(xdb(ii)-x(i)+xloc)+ 
     +                          abs(ydb(ii)-y(i)+yloc)+ 
     +                          abs(zdb(ii)-z(i)+zloc)
#endif
                        if(dtest.le.EPSLON) coloc = .true.
                  end do
c
c Only simulate this grid node if not colocated:
c
                  if(.not.coloc) then
                        nu     = nu + 1
#ifdef OPTI
                        xg(nu) = xyzdb(1,ii)
                        yg(nu) = xyzdb(2,ii)
                        zg(nu) = xyzdb(3,ii)
#else
                        xg(nu) = xdb(ii)
                        yg(nu) = ydb(ii)
                        zg(nu) = zdb(ii)
#endif
                        if(idbg.ge.3) write(ldbg,101) 
     +                                           nu,xg(nu),yg(nu),zg(nu)
 101              format(' Node ',i4,': x = ',f9.3,' y = ',f9.3,
     +                                             ' z = ',f9.3)
                  endif
            end do

c      write(*,*) '         * Nodes to simulate',nu
c      do i=1,nu
c            write(*,*) '         * ',i,xg(i),yg(i),zg(i)
c      end do
c
c Compute C22i: first get all the covariances:
c


            call cova3(xg(1),yg(1),zg(1),xg(1),yg(1),zg(1),
     +               1,nst,MAXNST,
     +           c0,it,cc,aa,1,MAXROT,rotmat,cmax,sill)
c      write(*,*) '         * Node to node covariances'
            do i = 1,nu
                  do j = i,nu
c            write(*,*) '         * xg,yg,zg i',xg(i),yg(i),zg(i)
c            write(*,*) '         * xg,yg,zg j',xg(j),yg(j),zg(j)
                        call cova3(xg(i),yg(i),zg(i),xg(j),yg(j),zg(j),
     +                        1,nst,MAXNST,
     +                        c0,it,cc,aa,1,MAXROT,rotmat,cmax,cov)
                        c22i(i,j) = cov
                        c22i(j,i) = cov
                        if(idbg.ge.3) write(ldbg,102) i,j,cov
c              if(idbg.ge.3) write(*,*) '         * i,j,cov',i,j,cov
 102  format(' Covariance between Node ',i4,' and ',i4,' = ',f9.6)
                  end do
            end do
c
c Compute C11, L11, inv(L11) AND inv(L11).Z1: First compute c11:
c
c      write(*,*) '         * Data to data covariances'
            do i=1,na
                  c11(i,i) = sill
                  do j=i+1,na
c            write(*,*) '         * xa,ya,za i',xa(i),ya(i),za(i)
c            write(*,*) '         * xa,ya,za j',xa(j),ya(j),za(j)
                        call cova3(xa(i),ya(i),za(i),xa(j),ya(j),za(j),
     +                             1,nst,MAXNST,
     +                             c0,it,cc,aa,1,MAXROT,rotmat,cmax,cov)
                        c11(i,j) = cov
                        c11(j,i) = cov
                        if(idbg.ge.3) write(ldbg,103) i,j,cov
c                 if(idbg.ge.3) write(*,*) '         * i,j,cov',i,j,cov
                  end do
            end do
 103  format(' Data-Data Covariance: ',i4,' and ',i4,' = ',f9.6)

c
c Compute l11:
c
            call chol(c11,l11,na,MAXSAM,ierr)
c
c Compute inv(l11):
c
            call linv(l11,l11inv,na,MAXSAM)
c
c Compute l11inv.z1:
c
            do i=1,na
                  v1(i) = 0.0
                  do k=1,i
                        v1(i) = v1(i) + l11inv(i,k)*vra(k)
                  end do
            end do
c
c Computation of LAMBDAe:
c
            sle = 0.0
            do i = 1,na
                  le(i) = 0.0
                  do k = 1,na
                        do j = 1,na
                              le(i) = le(i) + l11inv(j,i)*l11inv(j,k)
                        end do
                  end do
                  sle = sle + le(i)
            end do
c
c Compute C12 :
c
c      write(*,*) '         * Data to node covariances'
            do i = 1,na
                  do j = 1,nu
c            write(*,*) '         * xa,ya,za i',xa(i),ya(i),za(i)
c            write(*,*) '         * xg,yg,zg j',xg(j),yg(j),zg(j)
                        call cova3(xa(i),ya(i),za(i),xg(j),yg(j),zg(j),
     +                             1,nst,
     +                             MAXNST,c0,it,cc,aa,1,MAXROT,rotmat,
     +                             cmax,cov)
                        c12(i,j) = cov
                        if(idbg.ge.3) write(ldbg,104) i,j,cov
c                  if(idbg.ge.3) write(*,*) '         * i,j,cov',i,j,cov
            end do
            end do
 104  format(' Covariance between data ',i4,' and node ',i4,' = ',f9.6)
c
c Compute U12 = L11INV.C12:
c
c      write(*,*) '         * Compute U12'
            do i = 1,na
            do j = 1,nu
                  u12(i,j) = 0.0
                  do k = 1,i
                        u12(i,j) = u12(i,j) + l11inv(i,k)*c12(k,j)
                  end do
            end do
            end do
c
c Compute C22 - L21.U12 = C22 - U12^T.U12 ---> C22:
c
c      write(*,*) '         * Compute C22'
            do i = 1,nu
            do j = 1,nu
                  c22(i,j) = c22i(i,j)
                  do k = 1,na
                        c22(i,j) = c22(i,j) - u12(k,i)*u12(k,j)
                  end do
            end do
            end do
c
c Compute L22:
c
c      write(*,*) '         * Compute L22'
            call chol(c22,l22,nu,MAXDIS,ierr)


c
c Compute L21.L11INV.Z1 = U12^T.V1 ---> L21V1:
c
c      write(*,*) '         * Compute L21v1'
            do i = 1,nu
            l21v1(i) = 0.0
            do k = 1,na
                l21v1(i) = l21v1(i) + u12(k,i)*v1(k)
            end do
            end do
c
c Compute LAMBDAsk:
c
c      write(*,*) '         * Compute LAMBDAsk'
            do i = 1,nu
            slsk(i) = 0.0
            do k = 1,na
                  do j = k,na
                        slsk(i) = slsk(i) + u12(j,i)*l11inv(j,k)
                  end do
            end do
            end do
c
c Local bias correction:
c
c      write(*,*) '         * Compute Local bias correction',sle
            do i = 1,nu
            do j = 1,na
c                  print *,'in: l21v1(',i,')=',l21v1(i)
c                  print *,'in: slsk(',i,')=',slsk(i)
c                  print *,'in: le(',j,')=',le(j)
c                  print *,'in: vra(',j,')=',vra(j)
                  l21v1(i) = l21v1(i) + (1-slsk(i))*le(j)*vra(j)/sle
            end do
c                  print *,'out: l21v1(',i,')=',l21v1(i)
            end do
c
c LOOP OVER THE NUMBER OF SIMULATIONS      FOR THE PANNEL
c
            avgbl=0.0
            sumbl=0.0
            sqdbl=0.0
            simbl=0.0

c      write(*,*) '       *** Number realizations',nsim
            do isim=1,nsim
c      write(*,*) '       *** Realization',isim
c
c Generate (nxyz) Gaussian Random Numbers:
c
                  do i=1,nu
                        w2(i) = real(acorni(idum))
c                        w2(i) = 0.5
c      write(*,*) '         * rand',i,w2(i)
                  end do
                  do i = 1,nu
 112                          p = dble(w2(i))
                              call gauinv(p,xp,ierr)
                              w2(i) = xp
c      write(*,*) '         * gaus',i,w2(i)
                              if(w2(i).lt.(-6.).or.w2(i).gt.(6.)) then
                                         w2(i) = real(acorni(idum))
c                                         w2(i) = 0.5
                                    go to 112
                              endif
                  end do
c
c Compute the simulation at each grid point:
c

                  do i=1,nu
c                  write(*,*) '         * simu',i,l21v1(i)
                              y2(i) = l21v1(i)
                              do k=1,i
                                    y2(i) = y2(i) + l22(i,k)*w2(k)
                              end do
                  end do

c
c Write the simulation out to file:
c
                  ic = 0
                  sumsim=0
                  simbl=0

                  do ii=1,MAXDIS

#ifdef OPTI
                        zz = zloc + xyzdb(1,ii)
                        yy = yloc + xyzdb(2,ii)
                        xx = xloc + xyzdb(3,ii)
                        ic = ic + 1
                        dtest = abs(xyzdb(1,ii)-xg(ic)) + 
     +                          abs(xyzdb(2,ii)-yg(ic)) +
     +                          abs(xyzdb(3,ii)-zg(ic))

#else
                        zz = zloc + xdb(ii)
                        yy = yloc + ydb(ii)
                        xx = xloc + zdb(ii)
                        ic = ic + 1
                        dtest = abs(xdb(ii)-xg(ic)) + 
     +                          abs(ydb(ii)-yg(ic)) +
     +                          abs(zdb(ii)-zg(ic))
#endif
c
c Was this grid node simulated or was it co-located with a datum:
c
                        if(dtest.le.EPSLON) then

                              simval = y2(ic)
c      write(*,*) '         * ii,dtest,ic',ii,dtest,ic 
                        else
                              ic = ic - 1
                              do i=1,nd
                                    dtest = abs(xx-x(i)) + 
     +                                      abs(yy-y(i)) +
     +                                      abs(zz-z(i))
                                    if(dtest.le.EPSLON) then
                                          simval = vr(i)
c      write(*,*) '    coll * ii,dtest,ic',ii,dtest,ic 
                                          go to 2
                                    end if
                              end do
 2                            continue
                        end if
                        simval = backtr(simval,ntr,vrtr,vrgtr,zmin
     +                                  ,zmax,ltail,ltpar,utail,utpar)

                        if(simval.lt.zmin) simval = zmin
                        if(simval.gt.zmax) simval = zmax

c Write out the point simulation:
                        tipz=1+int((ii-1)/(nxdis*nydis))
                        tipy=1+int((ii-(tipz-1)*nxdis*nydis-1)/nxdis)
                        tipx=real(ii)-(tipz-1)*nxdis*nydis
     +                       -(tipy-1)*nxdis
                        ipx=tipx+(ix-1)*nxdis
                        ipy=tipy+(iy-1)*nydis
                        ipz=tipz+(iz-1)*nzdis
                        indpt=ipx+(ipy-1)*nxdis*nx 
     +                        +(ipz-1)*nx*nxdis*ny*nydis
c                        pdbg(indpt)=simval
c Compute the SMU values

                        tibz=int((nbz)*(tipz-1)/(nzdis))+1
                        tiby=int((nby)*(tipy-1)/(nydis))+1
                        tibx=int((nbx)*(tipx-1)/(nxdis))+1
                        ibx=tibx+(ix-1)*nbx
                        iby=tiby+(iy-1)*nby
                        ibz=tibz+(iz-1)*nbz
                        indb=ibx+(iby-1)*nbx*nx +(ibz-1)*nbx*nby*nx*ny
c                        if(isim.eq.nsim) sdbg(indb)=sdbg(indb)+simval
                        indsmu=tibx+(tiby-1)*nbx+(tibz-1)*nbx*nby
c Accumulate the simulated value within the pannel
                        simbl(indsmu)=simbl(indsmu)+simval
                  end do
c Accumulate to get the mean of SMU within the pannel for each realization
c and the variance of SMU grades wrt the pannel average for each realization

                  sumbl(isim)=0.0
                  sqdbl(isim)=0.0
                  do i=1,nsmup
                        simbl(i)=simbl(i)/MAXDIS*nsmup
                        sumbl(isim)=sumbl(isim)+simbl(i)
                        sqdbl(isim)=sqdbl(isim)+simbl(i)*simbl(i)
                        do j=1,nthr
                              if(simbl(i).gt.thres(j)) then
                                    prop(j)=prop(j)+1
                                    mean(j)=mean(j)+simbl(i)
                              end if
                        end do
                  end do
                  simbl=0.0

c
c END LOOP OVER SIMULATIONS:
c
            end do


            pavg = 0.0
            pvar = 0.0
            pgvar = 0.0

            do isim=1,nsim
                  sumbl(isim)=sumbl(isim)/nsmup
                  sqdbl(isim)=sqdbl(isim)/nsmup
                  pavg=pavg+sumbl(isim)
                  pvar=pvar+sqdbl(isim)-sumbl(isim)*sumbl(isim)
                  pgvar=pgvar+sqdbl(isim)
            end do

            sumbl=0.0
            sqdbl=0.0
            do j=1,nthr
                  if (prop(j).ne.0) then
                        mean(j)=mean(j)/prop(j)
                  else
                        mean(j)=0.0
                  end if
                  prop(j)=prop(j)/(nsmup*nsim)
            end do
            pavg=pavg/nsim
            pvar=pvar/nsim
            pgvar=pgvar/nsim-pavg*pavg

#ifdef _OPENMP
cc            write(100*lout+threadId+1,'(I6,30(f12.4))') 
            write(100*lout+threadId+1,'(30(f12.4))') 
     +                          pavg,pvar,pgvar,
     +                             (mean(j),j=1,nthr),
     +                             (prop(j),j=1,nthr)
#else
c            write(lout,'(I6,30(f12.4))') indpn,pavg,pvar,pgvar,
            write(lout,'(30(f12.4))') pavg,pvar,pgvar,
     +                             (mean(j),j=1,nthr),
     +                             (prop(j),j=1,nthr)
#endif

            pavg=0.0
            pvar=0.0
            pgvar=0.0
            mean=0.0
            prop=0.0
c
c END LOOP
c
      end do
c$omp end do
c$omp end parallel
     
#ifdef TRACE
      call extrae_fini()
#endif

      call cpu_time(timeEnd)

      write(*,*)'Time loop total=',(timeEnd-timeIni)
      write(*,*)'Time loop iter=',(timeEnd-timeIni)/(real(nxyz))

c      do i=1,nxyz*MAXDIS
c            write(ldbg2,'(f12.4)') pdbg(i)
c      end do
c Write out the SMU simulation
      do i=1,nxyz*nsmup
            sdbg(i)=sdbg(i)*nsmup/MAXDIS
c            write(ldbg3,'(f12.4)') sdbg(i)
      end do
      close(ldbg2)
c      close(ldbg3)

c
c Finished:
c
      close(ldbg)
      close(lout)
#ifdef _OPENMP
      do i=1,numThreads
            close(lout*100+i)
      end do
#endif

      write(*,9998) VERSION

 9998 format(/' BLUSIM Version: ',f5.3, ' Finished'/)
      stop
 98   stop 'ERROR in parameter file!'
 99   stop 'ERROR in data file!'
      end


      subroutine chol(a,t,n,ndim,ierr)
c-----------------------------------------------------------------------
c
c                      Cholesky Decomposition
c                      **********************
c
c This subroutine calculates the lower triangular matrix T which, when
c multiplied by its own transpose, gives the symmetric matrix A. (from
c "Numerical Analysis of Symmetric Matrices,"  H.R. Schwarz et al.,
c p. 254)
c
c
c
c INPUT VARIABLES:
c
c   a(n,n)           Symmetric positive definite matrix to be
c                      decomposed (destroyed in the calculation of t)
c   t(n,n)           Lower triangular matrix solution
c   n                Dimension of the system you're decomposing
c   ndim             Dimension of matrix a (Note: In the main program,
c                      matrix a may have been dimensioned larger than
c                      necessary, i.e. n, the size of the system you're
c                      decomposing, may be smaller than ndim.)
c   ierr             Error code:  ierr=0 - no errors; ierr=1 - matrix a
c                      is not positive definite
c
c
c
c NO EXTERNAL REFERENCES:
c-----------------------------------------------------------------------
      dimension a(ndim,ndim),t(ndim,ndim)
      ierr = 0
c
c Check for positive definiteness:
c
      do ip=1,n
            if(a(ip,ip).le.0.0) then
                  write(*,'(a)') 'WARNING: chol - not positive definite'
                  ierr = 1
                  go to 1
            endif
            t(ip,ip) = sqrt (a(ip,ip))
            if(ip.ge.n) return
            do k = ip+1,n
                  t(k,ip) = a(ip,k)/t(ip,ip)
            end do
            do i = ip+1,n
                  do k = i,n
                        a(i,k) = a(i,k) - t(i,ip) * t(k,ip)
                  end do
            end do
 1          continue
      end do
c
c Finished:
c
      return
      end
 
 
 
      subroutine linv(a,b,n,ndim)
c-----------------------------------------------------------------------
c
c                Inverse of a Lower Triangular Matrix
c                ************************************
c
c This subroutine finds the inverse of a lower triangular matrix A and
c stores the answer in B. (from "Numerical Analysis of Symmetric
c Matrices,"  H.R. Schwarz et al.,)
c
c
c
c INPUT VARIABLES:
c
c   a(n,n)           Lower triangular matrix to be inverted
c   b(n,n)           the inverse
c   n                Dimension of the matrix you're inverting
c   ndim             Dimension of matrix a (Note: In the main program,
c                      matrix a may have been dimensioned larger than
c                      necessary, i.e. n, the size of the system you're
c                      decomposing, may be smaller than ndim.)
c
c-----------------------------------------------------------------------
      dimension a(ndim,ndim),b(ndim,ndim)
      do i = 1,n
            if(i.gt.1) then
                  do k = 1,i-1
                        sum=0.0
                        do j = k,i-1
                              sum = sum + a(i,j)*b(j,k)
                        end do
                        b(i,k) = -sum/a(i,i)
                  end do
            end if
            b(i,i) = 1./a(i,i)
      end do
c
c Finished:
c
      return
      end


      subroutine setrot(ang1,ang2,ang3,anis1,anis2,ind,MAXROT,rotmat)
c-----------------------------------------------------------------------
c
c              Sets up an Anisotropic Rotation Matrix
c              **************************************
c
c Sets up the matrix to transform cartesian coordinates to coordinates
c accounting for angles and anisotropy (see manual for a detailed
c definition):
c
c
c INPUT PARAMETERS:
c
c   ang1             Azimuth angle for principal direction
c   ang2             Dip angle for principal direction
c   ang3             Third rotation angle
c   anis1            First anisotropy ratio
c   anis2            Second anisotropy ratio
c   ind              matrix indicator to initialize
c   MAXROT           maximum number of rotation matrices dimensioned
c   rotmat           rotation matrices
c
c
c NO EXTERNAL REFERENCES
c
c
c-----------------------------------------------------------------------
      parameter(DEG2RAD=3.141592654/180.0,EPSLON=1.e-20)
      real*8    rotmat(MAXROT,3,3),afac1,afac2,sina,sinb,sint,
     +          cosa,cosb,cost
c
c Converts the input angles to three angles which make more
c  mathematical sense:
c
c         alpha   angle between the major axis of anisotropy and the
c                 E-W axis. Note: Counter clockwise is positive.
c         beta    angle between major axis and the horizontal plane.
c                 (The dip of the ellipsoid measured positive down)
c         theta   Angle of rotation of minor axis about the major axis
c                 of the ellipsoid.
c
      if(ang1.ge.0.0.and.ang1.lt.270.0) then
            alpha = (90.0   - ang1) * DEG2RAD
      else
            alpha = (450.0  - ang1) * DEG2RAD
      endif
      beta  = -1.0 * ang2 * DEG2RAD
      theta =        ang3 * DEG2RAD
      write(*,*) ' a,b,t',alpha,beta,theta
c
c Get the required sines and cosines:
c
      sina  = dble(sin(alpha))
      sinb  = dble(sin(beta))
      sint  = dble(sin(theta))
      cosa  = dble(cos(alpha))
      cosb  = dble(cos(beta))
      cost  = dble(cos(theta))
      write(*,*) sina,cosa
      write(*,*) sinb,cosb
      write(*,*) sint,cost
c
c Construct the rotation matrix in the required memory:
c
      afac1 = 1.0 / dble(max(anis1,EPSLON))
      afac2 = 1.0 / dble(max(anis2,EPSLON))
      rotmat(ind,1,1) =       (cosb * cosa)
      rotmat(ind,1,2) =       (cosb * sina)
      rotmat(ind,1,3) =       (-sinb)
      rotmat(ind,2,1) = afac1*(-cost*sina + sint*sinb*cosa)
      rotmat(ind,2,2) = afac1*(cost*cosa + sint*sinb*sina)
      rotmat(ind,2,3) = afac1*( sint * cosb)
      rotmat(ind,3,1) = afac2*(sint*sina + cost*sinb*cosa)
      rotmat(ind,3,2) = afac2*(-sint*cosa + cost*sinb*sina)
      rotmat(ind,3,3) = afac2*(cost * cosb)
c
c Return to calling program:
c
      return
      end



      subroutine makepar
c-----------------------------------------------------------------------
c
c                      Write a Parameter File
c                      **********************
c
c
c
c-----------------------------------------------------------------------
      lun = 99
      open(lun,file='blusim.par',status='UNKNOWN')
      write(lun,10)
 10   format('                  Parameters for BLUSIM',/,
     +       '                  *********************',/,/,
     +       'START OF PARAMETERS:')

      write(lun,11)
 11   format('cluster.dat                      ',
     +       '-file with data')
      write(lun,12)
 12   format('1  2  3  4  5                    ',
     +       '-   columns for X,Y,Z,var,wt')
      write(lun,13)
 13   format('-1.0e21   1.0e21                 ',
     +       '-   trimming limits')
      write(lun,14)
 14   format('1                             ',
     +       '-transform the data (0=no, 1=yes)')
      write(lun,15)
 15   format('blusim.trn                    ',
     +       '-  file for output trans table')
      write(lun,16)
 16   format('0                             ',
     +       '-  consider ref. dist (0=no, 1=yes)')
      write(lun,17)
 17   format('histsmth.out                  ',
     +       '-  file with ref. dist distribution')
      write(lun,18)
 18   format('1  2                             ',
     +       '-  columns for vr and wt')
      write(lun,19)
 19   format('0.0    15.0                      ',
     +       '-  zmin,zmax(tail extrapolation)')
      write(lun,20)
 20   format('1       0.0                      ',
     +       '-  lower tail option, parameter')
      write(lun,21)
 21   format('1      15.0                      ',
     +       '-  upper tail option, parameter')
      write(lun,22)
 22   format('3                                ',
     +       '-debugging level: 0,1,2,3')
      write(lun,23)
 23   format('blusim.dbg                       ',
     +       '-file for debugging output')
      write(lun,24)
 24   format('blusim.out                       ',
     +       '-file for kriged output')
      write(lun,25)
 25   format('100                              ',
     +       '-number of realizations to generate')
      write(lun,26)
 26   format('50   0.5    100.0                ',
     +       '-nx,xmn,xsiz - Pannel Size')
      write(lun,27)
 27   format('50   0.5    100.0                ',
     +       '-ny,ymn,ysiz')
      write(lun,28)
 28   format('1    0.5    1.0                  ',
     +       '-nz,zmn,zsiz')
      write(lun,29)
 29   format('5   5   5                        ',
     +       '-nbx,nby,nbz- SMUs per Pannel')
      write(lun,32)
 32   format('3    3   3                       ',
     +       '-x,y and z block discretization')
      write(lun,33)
 33   format('5                                ',
     +       '-number of cutoffs for reporting')
      write(lun,34)
 34   format('0.2 0.5 0.8 1.0 2.0              ',
     +       '-cutoffs')
      write(lun,35)
 35   format('9784585                          ',
     +       '-random number seed')
      write(lun,36)
 36   format('4    8                           ',
     +       '-min, max data for kriging')
      write(lun,37)
 37   format('0                                ',
     +       '-max per octant (0-> not used)')
      write(lun,38)
 38   format('20.0  20.0  20.0                 ',
     +       '-maximum search radii')
      write(lun,39)
 39   format(' 0.0   0.0   0.0                 ',
     +       '-angles for search ellipsoid')
      write(lun,40)
 40   format('1    0.2                         ',
     +       '-nst, nugget effect')
      write(lun,41)
 41   format('1    0.8  0.0   0.0   0.0        ',
     +       '-it,cc,ang1,ang2,ang3')
      write(lun,42)
 42   format('         10.0  10.0  10.0        ',
     +       '-a_hmax, a_hmin, a_vert')

      close(lun)
      return
      end
