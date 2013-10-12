


cccccccccccccccccccccccccc
c
c As stand-alone program:
c
c     1) Set PRODUCTION 0, 1, or 2
c     2) To compile only:
c            gfortran -Wall -cpp -g -O2 test.f -o test.e
c        To compile/run with debug info, do:
c            cp ../test.f . ; gfortran -Wall -cpp -g -O2 test.f -o test.e ; ./test.e > test.e.out
c        To compile/run with debug info and gdb easy reading, do:
c            cp ../test.f . ; gfortran -Wall -cpp -g -O0 test.f -o test.e ; ./test.e > test.e.out
c        If crashes, then do:
c            gdb ./test.e core

c     To test for bad warnings:
c     cp ../test.f . ; gfortran -Wall -cpp -g -O2 test.f -o test.e &> test.e.makelog ; grep -v "Unused dummy" test.e.makelog | grep -v "Unused variable" | grep -i warning 


cccccccccccccccccccccccccc
c
c As harm subroutine
c
c 1) Set PRODUCTION 3
c 2) fpp -P test.f > testfpp.temp.f ; f2c -Wall -f -P testfpp.temp.f ; sed 's/static//g' testfpp.temp.c > testfpp.c ;sed 's/static//g' testfpp.temp.P > testfpp.P
c    i.e. testfpp.c from f2c: MUST remove static in front of variables! 


c about pre-processor directives:
c http://gcc.gnu.org/onlinedocs/gfortran/Preprocessing-Options.html
c
c 0 means normal full output
c 1 means very simple Jon's version of err and sol output
c 2 means only Jon's err output
c 3 means no output or input (harm mode)
c#define PRODUCTION 0
c#define PRODUCTION 1
c#define PRODUCTION 2
#define PRODUCTION 3

#define VEL4 0
#define VEL3 1
#define VELREL4 2
c 0 : iterate 4-velocities u^i=prim(2,3,4) (radiation is inverted separately and values for u^i_{rad} computed separately)
c 1 : iterate 3-velocities (not used)
c 2 : iterate relative 4-velocities \tilde{u}^i = prim(2,3,4) ("")

c     WHICHVEL should be set same as in harm
#define WHICHVEL VELREL4
c#define WHICHVEL VEL4

c error below which will consider BAD solution and not even compute final solution
c#define FAILLIMIT (1D-6)
c Choose actual harm choice even if pretty strict.
c#define FAILLIMIT (1D-8)
#define FAILLIMIT (tolallow)

#define NUMARGS (211+11)
c 11 vars, failcode, error, iterations
#define NUMRESULTS 15

#define SMALL (1D-300)

#if(PRODUCTION<3)
      program testradiation
      double precision args(NUMARGS)
c results for energy
      double precision resultseng(NUMRESULTS)
c results for entropy
      double precision resultsent(NUMRESULTS)
c      args just dummy here, not used except by external call to rameshsolver subroutine
      call rameshsolver(args,resultseng,resultsent)
      stop
      end
#endif







      subroutine rameshsolver(args
     &     ,resultseng,resultsent)
c     Reads in Jon's error file fails.dat and tests our ideas for
c     applying the radiation source term

c     variables with 'p' at the end refer to the previous time step,
c     those with 'i' at the end refer to the current time step
c     pre-radiation, and those with 'f' at the end are for the current
c     time step post-radiation.

c     variables with 'c' at the end are conserved quantities
c     corresponding to the problem currently being solved

c     variables with 'final' at the end correspond to the final solution
      implicit double precision (a-h,o-z)
c     inputs
      dimension args(NUMARGS)
c     results
      double precision resultseng(NUMRESULTS)
      double precision resultsent(NUMRESULTS)
c     internals
      real*8 itertot
      dimension isc(4),gn(4,4),gv(4,4),hp(4,4),hf(4,4)
      dimension vgasp(4),vgasf(4),B1(5),B2(5),
     &     BBp(4),BBf(4),BBc(4),vradp(4),vradf(4),s(5),
     &     ugasconf(4),ugascovf(4),ugasconp(4),ugascovp(4),
     &     uradconf(4),uradcovf(4),uradconp(4),uradcovp(4),
     &     ugasconi(4),ugascovi(4),uradconi(4),uradcovi(4),
     &     bconp(4),bcovp(4),bconf(4),bcovf(4),bconi(4),bcovi(4),
     &     Tmunup(4,4),Tmunuf(4,4),Rmunup(4,4),Rmunuf(4,4),
     &     Tmunui(4,4),Rmunui(4,4),Ttcovi(4),Rtcovi(4),
     &     Ttcovc(4),Rtcovc(4),Ttotc(4),prim(4),error(4),
     &     primeng(4),priment(4),
     &     ugasconfinal(4),ugascovfinal(4),uradconfinal(4),
     &     uradcovfinal(4),bconfinal(4),bcovfinal(4),
     &     Ttcovfinal(4),Rtcovfinal(4),
     &     Tmunufinal(4,4),Rmunufinal(4,4)
      common/metric/gn,gv
      common/accuracy/eps,epsbis,dvmin,tol,uminfac,dlogmax,
     &     gammaradceiling,itermax
      common/conserved/Gam,Gam1,en,en1,rhou0c,sc,Ttcovc,Rtcovc,
     &     BBc,dt,
     &     FNUMARGSHARM,FNUMRESULTSHARM,WHICHVELRAMESHHARM,
     &     GAMMAMAXRAD,ERADLIMIT,toltry,tolallow,
     &     ARAD_CODE,OKAPPA_ES_CODE11,OKAPPA_FF_CODE11,OKAPPA_BF_CODE11
      external funcMHD1,funcMHD2,funcrad1,funcrad2
      character filehead*30,fileext*4,infile*60
      character conhead*3,confile*60
      character solhead*3,solfile*60
      character errhead*3,errfile*60
      double precision mymax,myabs,mydiv
      integer myisnan
      integer which,showboth
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc



ccccccccccccccccccccccccccccccccccccccccccccccccccc
c     origintype=0 means read from file
c     origintype=1 means read from function args as if called by C code using latest number of columns, in which case don't write anything, just pass back result.
      if(PRODUCTION==3) then
         origintype=1
      else
         origintype=0
      endif
ccccccccccccccccccccccccccccccccccccccccccccccccccc
c     idatatype=1 means Jon's old data format with 134 numbers
c     idatatype=2 means Jon's new data format with 181 numbers
c     idatatype=3 means Jon's new data format with 211+11 numbers

c      write (*,*) ' which type data file? old(1) new(2) '
c      read (*,*) idatatype
      idatatype=3

ccccccccccccccccccccccccccccccccccccccccccccccccccc
c     iMHD=0 means don't do initial MHD inversion before radiation
c     iMHD=1 means do initial MHD inversion before radiation

      iMHD=0
      itertot=0
      iproblem=0



c     Read in data from Jon's datafile


#if(PRODUCTION<=2)
c     CHOOSE FILENAME HEAD TO READ
      filehead='fails'
c      filehead='fails1'
c      filehead='failsearly'
c      
c     CAN CHOOSE EXTENSION
      fileext='.txt'
      infile=trim(filehead)//fileext
c
c     CAN CHOOSE HEADER FOR OUTPUTS
      conhead='con'
      confile=conhead//'_'//trim(filehead)//fileext
      solhead='sol'
      solfile=solhead//'_'//trim(filehead)//fileext
      errhead='err'
      errfile=errhead//'_'//trim(filehead)//fileext
c
      open (11,file=infile)
#endif
#if(PRODUCTION<=0)
      open (12,file=confile)
#endif
#if(PRODUCTION<=1)
      open (13,file=solfile)
#endif
#if(PRODUCTION<=2)
      open (14,file=errfile)
#endif
#if(PRODUCTION<=2)
      write (14,"(9X,1A)",advance="no") 'NUM HITER HTOTITER HARMERR'
      write (14,"(15X,1A)",advance="no") 'PREVBESTHARMERR'
      write (14,"(9X,1A)",advance="no")
     &     'HARMEOMTYPE HARMITERMODE'
      write (14,"(2X,1A)",advance="no") 'HARMSTAT'
      write (14,"(1X,1A)",advance="no") 'RAMENGSTAT'
      write (14,"(4X,1A)",advance="no") 'RAMENGITER   RAMENGERR'
      write (14,"(11X,1A)",advance="no") 'RAMENTSTAT'
      write (14,"(5X,1A)") 'RAMENTITER RAMENTERR'
#endif

#if(PRODUCTION==3)
      do i=1,1
#else
      do i=1,1000000
#endif

c     If ientropy=0, we will try the energy equation. If it is 1, we
c     proceed directly to the entropy equation

#if(PRODUCTION<=2)
         write (14,"(1I10)",advance="no") i
#endif
#if(PRODUCTION<=1)
         if(WHICHVEL.eq.VEL4) then
            write (13,*)
            write (13,"(54X,22X,A,22X,A)",advance="no")
     &           '|rho','|u_g'
            write (13,"(22X,A,22X,A,22X,A,22X,A)",advance="no")
     &           '|u^t','u^1','u^2','u^3|'
            write (13,"(22X,A)",advance="no")
     &           '|Erf'
            write (13,"(22X,A,22X,A,22X,A,22X,A)")
     &           '|u^t','u^1','u^2','u^3|'
            write (13,*) ' JON PROBLEM NUMBER: ',i
         else if(WHICHVEL.eq.VELREL4) then
            write (13,*)
            write (13,"(54X,22X,A,22X,A)",advance="no")
     &           '|rho','|u_g'
            write (13,"(22X,A,22X,A,22X,A,22X,A)",advance="no")
     &           '|tu^t','tu^1','tu^2','tu^3|'
            write (13,"(22X,A)",advance="no")
     &           '|Erf'
            write (13,"(22X,A,22X,A,22X,A,22X,A)")
     &           '|tu^t','tu^1','tu^2','tu^3|'
            write (13,*) ' JON PROBLEM NUMBER: ',i            
         else
            stop
         endif
#endif
         ientropy=0
         
         if (idatatype.eq.1) then

         call readtype1(gn,gv,rhof,rhop,rhou0i,rhou0f,rhou0p,
     &        uf,up,T00i,T00f,T00p,vgasf,vgasp,T01i,T01f,T01p,
     &        T02i,T02f,T02p,T03i,T03f,T03p,BBf,BBp,
     &        Ef,Ep,R00i,R00f,R00p,vradf,vradp,R01i,R01f,R01p,
     &        R02i,R02f,R02p,R03i,R03f,R03p,s,si,sf,sp,
     &        uradconf,uradcovf,ugasconf,ugascovf,
     &        uradconp,uradcovp,ugasconp,ugascovp,ifinish,errorabs)
         FNUMARGSHARM=NUMARGS
         FNUMRESULTSHARM=NUMRESULTS
         WHICHVELRAMESHHARM=WHICHVEL
         Gam=4.d0/3.d0
         GAMMAMAXRAD=50.0d0
         ERADLIMIT=10.0D0*1d-300
         toltry=1d-12
         tolallow=1d-8
         ARAD_CODE=1.18316d17
         OKAPPA_ES_CODE11=5.90799d5
         OKAPPA_FF_CODE11=3.46764d-17
         OKAPPA_BF_CODE11=6.93528d-15

         else if (idatatype.eq.2) then

         call readtype2(gn,gv,rhof,rhop,rhou0i,rhou0f,rhou0p,
     &        uf,up,T00i,T00f,T00p,vgasf,vgasp,T01i,T01f,T01p,
     &        T02i,T02f,T02p,T03i,T03f,T03p,BBf,BBp,
     &        Ef,Ep,R00i,R00f,R00p,vradf,vradp,R01i,R01f,R01p,
     &        R02i,R02f,R02p,R03i,R03f,R03p,s,si,sf,sp,
     &        uradconf,uradcovf,ugasconf,ugascovf,
     &        uradconp,uradcovp,ugasconp,ugascovp,ifinish,errorabs)

         FNUMARGSHARM=NUMARGS
         FNUMRESULTSHARM=NUMRESULTS
         WHICHVELRAMESHHARM=WHICHVEL
c         Gam=4.d0/3.d0
         GAMMAMAXRAD=50.0d0
         ERADLIMIT=10.0D0*1D-300
         toltry=1d-12
         tolallow=1d-8
         ARAD_CODE=1.18316d17
         OKAPPA_ES_CODE11=5.90799d5
         OKAPPA_FF_CODE11=3.46764d-17
         OKAPPA_BF_CODE11=6.93528d-15

         else if (idatatype.eq.3) then

         call readtype3(origintype,args,
     &        gn,gv,rhof,rhop,rhou0i,rhou0f,rhou0p,
     &        uf,up,T00i,T00f,T00p,vgasf,vgasp,T01i,T01f,T01p,
     &        T02i,T02f,T02p,T03i,T03f,T03p,BBf,BBp,
     &        Ef,Ep,R00i,R00f,R00p,vradf,vradp,R01i,R01f,R01p,
     &        R02i,R02f,R02p,R03i,R03f,R03p,s,si,sf,sp,
     &        uradconf,uradcovf,ugasconf,ugascovf,
     &        uradconp,uradcovp,ugasconp,ugascovp,ifinish,errorabs)

         endif

c     Initialize constants and read in data from Jon's datafile

c     eps is fractional shift in primitives for numerical derivatives
c     epsbib is fractional accuracy desired in bisection
c     dvmin: if abs(uconmu)<1d-4 then shift used is dvmin*eps
c     tol: all fractional errors must be below tol for convergence
c     uminfac: minimum u_g allowed is uminfrac*rho
c     dlogmax: minimum log() in entropy is log() conserved - dlogmin
c     itermax is maximum no. of iterations with u_g, entropy below minimum

         eps=1.d-6
         epsbis=1.d-3
         dvmin=1.d-4
c     tol=1.d-10
c     harm-like error
         tol=toltry
         uminfac=1.d-10
         dlogmax=log(2.d0)
         itermax=3
         gammaradceiling=GAMMAMAXRAD


c     Gam=adiabatic index Gamma, Gam1=Gamma-1
c     en=polytropic index = 1/Gamma, en1 = n+1

         Gam1=Gam-1.d0
         en=1.d0/Gam1
         en1=en+1.d0

         

#if(PRODUCTION==0)
c         write (*,*)
c         write (*,*) (isc(j),j=1,4),dt,
c     &        ((gn(j,k),k=1,4),j=1,4),
c     &        ((gv(j,k),k=1,4),j=1,4),
c     &        rhof,rhop,rhou0i,rhou0f,rhou0p,
c     &        uf,up,T00i,T00f,T00p,
c     &        vgasf(2),vgasp(2),T01i,T01f,T01p,
c     &        vgasf(3),vgasp(3),T02i,T02f,T02p,
c     &        vgasf(4),vgasp(4),T03i,T03f,T03p,
c     &        BBf(2),BBp(2),scr,scr,scr,
c     &        BBf(3),BBp(3),scr,scr,scr,
c     &        BBf(4),BBp(4),scr,scr,scr,
c     &        Ef,Ep,R00i,R00f,R00p,
c     &        vradf(2),vradp(2),R01i,R01f,R01p,
c     &        vradf(3),vradp(3),R02i,R02f,R02p,
c     &        vradf(4),vradp(4),R03i,R03f,R03p,
c     &        (s(j),j=1,2),si,sf,sp,
c     &        (uradconf(j),uradcovf(j),j=1,4),
c     &        (ugasconf(j),ugascovf(j),j=1,4),
c     &        (uradconp(j),uradcovp(j),j=1,4),
c     &        (ugasconp(j),ugascovp(j),j=1,4)
#endif         
         if (ifinish.eq.1) go to 10
         iproblem=iproblem+1

c     Calculate b^mu and b_mu for 'p' and 'f'

         call calcbconbcov(BBp,ugasconp,ugascovp,bconp,bcovp,bsqp)
         call calcbconbcov(BBf,ugasconf,ugascovf,bconf,bcovf,bsqf)

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

#if(PRODUCTION==0)
         write (*,*)
         write (*,*)
         write (*,*) ' JON PROBLEM NUMBER: ',i
         write (12,*)
         write (12,*) ' JON PROBLEM NUMBER: ',i
#endif

c     Set up initial guess primitives ('p' quantities from prevous time
c     step). Make sure up is reasonable. If not, reset up to
c     umin=uminfac*rhop. 

c      CHOOSE to use previous primitve or harm solution as starting point.
c     0 = use original prims (normal mode)
c     1 = use harm prims and still do normal Ramesh stages
c     2 = use harm prims and go straight to 4D iterations to avoid mis-steps (works to get some entropy solutions not otherwise found by Ramesh code, but misses *many* energy solutions even though I'm providing very close guess.)
         guesstype=0
c         guesstype=2
c         guesstype=1

         if(guesstype.eq.0.or.errorabs.gt.1E-9) then
#if(PRODUCTION==0)
            write(*,*) 'Using original primitives as guess: guesstype='
     &           ,guesstype
#endif
            prim(1)=mymax(up,uminfac*rhop)
            if(WHICHVEL.eq.VEL4) then
               prim(2)=ugasconp(2)
               prim(3)=ugasconp(3)
               prim(4)=ugasconp(4)
            else if(WHICHVEL.eq.VELREL4) then
               prim(2)=vgasp(2)
               prim(3)=vgasp(3)
               prim(4)=vgasp(4)               
            else
               stop
            endif
            guessstype=0
         else
#if(PRODUCTION==0)
            write(*,*) 'Using harm primitives as guess: guesstype='
     &           ,guesstype
#endif
            prim(1)=mymax(uf,uminfac*rhof)
            if(WHICHVEL.eq.VEL4) then
               prim(2)=ugasconf(2)
               prim(3)=ugasconf(3)
               prim(4)=ugasconf(4)
            else if(WHICHVEL.eq.VELREL4) then
               prim(2)=vgasf(2)
               prim(3)=vgasf(3)
               prim(4)=vgasf(4)
            else
               stop
            endif
         endif


c     put any overrides for guesses here
c         prim(1)=1.1081585946596330690836004556281327032809D-8
c         prim(2)=-0.075584957687290353262066821240292555500
c         prim(3)=-0.003678584407270331563545276169512108269
c         prim(4)=0.0923803797665372578214676962826934511966
         



#if(PRODUCTION==0)
         write (*,*) ' prim: ',(prim(j),j=1,4)
#endif

c     Set up conserved quantities 'c' of the current time step ('i'
c     state). The values go into common block /conserved/


         rhou0c=rhou0i
         sc=si

         Ttcovc(1)=T00i
         Ttcovc(2)=T01i
         Ttcovc(3)=T02i
         Ttcovc(4)=T03i

         Rtcovc(1)=R00i
         Rtcovc(2)=R01i
         Rtcovc(3)=R02i
         Rtcovc(4)=R03i

         BBc(1)=0.d0
         BBc(2)=BBf(2)
         BBc(3)=BBf(3)
         BBc(4)=BBf(4)

      if (iMHD.eq.1) then

c     Carry out initial MHD inversion of the 'i' state

#if(PRODUCTION==0)
         write (*,*) 
         write (*,*) ' MHD INVERSION OF PRE-RADIATION SOLUTION '
#endif       
         call MHDinvert(prim,iflag,jflag,ientropy,itertot,guesstype)

c     MHD inversion done. Update all the relevant quantities
c     corresponding to the 'i' solution

c     Calculate u, u^0 and rho

         ui=prim(1)
         call solveucon(prim,ugasconi)
         call contocov(ugasconi,ugascovi)
         rhoi=rhou0i/ugasconi(1)

#if(PRODUCTION==0)
c         write (*,*)
c         write (*,*) ' final solution: rho, u, u^mu: ',i,
c     &        rhoi,ui,(ugasconi(j),j=1,4)
c         write (*,*)
c         write (*,*) ' original T^t_mu: ',(Ttcovc(j),j=1,4)
c         write (*,*) ' original R^t_mu: ',(Rtcovc(j),j=1,4)
#endif

c     Calculate b^mu and b_mu, and update T^t_mu and R^t_mu so as to be
c     consistent with primitives

         do j=1,4
            Ttotc(j)=Ttcovc(j)+Rtcovc(j)
         enddo

         call calcbconbcov(BBc,ugasconi,ugascovi,bconi,bcovi,bsqi)

         call calcTmunu(rhoi,ui,bsqi,Gam,ugasconi,ugascovi,
     &        bconi,bcovi,Tmunui)

         do j=1,4
            Ttcovi(j)=Tmunui(1,j)
            Ttcovc(j)=Ttcovi(j)
            Rtcovi(j)=Ttotc(j)-Ttcovi(j)
            Rtcovc(j)=Rtcovi(j)
         enddo

#if(PRODUCTION==0)
c         write (*,*)
c         write (*,*) ' updated T^t_mu: ',(Ttcovc(j),j=1,4)
c         write (*,*) ' updated R^t_mu: ',(Rtcovc(j),j=1,4)
#endif

c     Initial MHD inversion done

      else

c     No MHD inversion done. Proceed directly to solve the radiation
c     problem.

#if(PRODUCTION==0)
         write (*,*) 
         write (*,*) ' SKIP MHD INVERSION OF PRE-RADIATION SOLUTION '
#endif

      endif

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c     Now it is time to apply the implicit radiation source term and
c     solve for the post-radiation primitives

c     common block /conserved/ already contains the relevant conserved
c     quantities for the pre-radiation state. For the primitives, use
c     the 'i' state as the initial guess. Then solve the lab frame
c     energy equation using Newton-Raphson

#if(PRODUCTION==0)
         write (*,*)
         write (*,*) ' SOLVE WITH IMPLICIT RADIATION SOURCE TERM '
#endif
         call radsource(prim,iter,iflag,jflag,ientropy
     &        ,itertot,guesstype,primeng,priment
     &        ,resultseng,resultsent)

c     Radiation inversion done. Calculate relevant quantities
c     corresponding to the final solutions
         showboth=nint(resultseng(12))+nint(resultsent(12))
         if(nint(resultseng(12)).eq.0) then
            which=3
            call getfinal(which,showboth,resultseng,
     &     isc,gn,gv,hp,hf,
     &     vgasp,vgasf,B1,B2,
     &     BBp,BBf,BBc,vradp,vradf,s,
     &     ugasconf,ugascovf,ugasconp,ugascovp,
     &     uradconf,uradcovf,uradconp,uradcovp,
     &     ugasconi,ugascovi,uradconi,uradcovi,
     &     bconp,bcovp,bconf,bcovf,bconi,bcovi,
     &     Tmunup,Tmunuf,Rmunup,Rmunuf,
     &     Tmunui,Rmunui,Ttcovi,Rtcovi,
c     prim->primeng
     &     Ttcovc,Rtcovc,Ttotc,primeng,error,
     &     primeng,priment,
     &     ugasconfinal,ugascovfinal,uradconfinal,
     &     uradcovfinal,bconfinal,bcovfinal,
     &     Ttcovfinal,Rtcovfinal,
     &     Tmunufinal,Rmunufinal,rhou0c,Gam)
         endif
         if(nint(resultsent(12)).eq.0) then
            which=2
            call getfinal(which,showboth,resultsent,
     &     isc,gn,gv,hp,hf,
     &     vgasp,vgasf,B1,B2,
     &     BBp,BBf,BBc,vradp,vradf,s,
     &     ugasconf,ugascovf,ugasconp,ugascovp,
     &     uradconf,uradcovf,uradconp,uradcovp,
     &     ugasconi,ugascovi,uradconi,uradcovi,
     &     bconp,bcovp,bconf,bcovf,bconi,bcovi,
     &     Tmunup,Tmunuf,Rmunup,Rmunuf,
     &     Tmunui,Rmunui,Ttcovi,Rtcovi,
c     prim->priment
     &     Ttcovc,Rtcovc,Ttotc,priment,error,
     &     primeng,priment,
     &     ugasconfinal,ugascovfinal,uradconfinal,
     &     uradcovfinal,bconfinal,bcovfinal,
     &     Ttcovfinal,Rtcovfinal,
     &     Tmunufinal,Rmunufinal,rhou0c,Gam)
        endif



c     enddo below is loop over cases
      enddo

 10   continue
      problem=iproblem

#if(PRODUCTION==0)
      write (*,*)
      write (*,*)
      write (*,*) ' iproblem: ',iproblem
      write (*,*) ' itertot: ',itertot
      write (*,*) ' iterations per problem: ',itertot/problem
      write (12,*)
      write (12,*) ' iproblem: ',iproblem
      write (12,*) ' itertot: ',itertot
      write (12,*) ' iterations per problem: ',itertot/problem
      write (*,*)
#endif
#if(PRODUCTION<=2)
      close (11)
#endif
#if(PRODUCTION<=0)
      close (12)
#endif
#if(PRODUCTION<=1)
      close (13)
#endif
#if(PRODUCTION<=2)
      close (14)
#endif
      end



c     Once radiation inversion done, calculate relevant quantities
c     corresponding to the final solutions
      subroutine getfinal(which,showboth,results,
     &     isc,gn,gv,hp,hf,
     &     vgasp,vgasf,B1,B2,
     &     BBp,BBf,BBc,vradp,vradf,s,
     &     ugasconf,ugascovf,ugasconp,ugascovp,
     &     uradconf,uradcovf,uradconp,uradcovp,
     &     ugasconi,ugascovi,uradconi,uradcovi,
     &     bconp,bcovp,bconf,bcovf,bconi,bcovi,
     &     Tmunup,Tmunuf,Rmunup,Rmunuf,
     &     Tmunui,Rmunui,Ttcovi,Rtcovi,
c     prim that's used
     &     Ttcovc,Rtcovc,Ttotc,prim,error,
     &     primeng,priment,
     &     ugasconfinal,ugascovfinal,uradconfinal,
     &     uradcovfinal,bconfinal,bcovfinal,
     &     Ttcovfinal,Rtcovfinal,
     &     Tmunufinal,Rmunufinal,rhou0c,Gam)
      implicit double precision (a-h,o-z)
      integer which,showboth
      double precision results(NUMRESULTS)
      integer radinvmod
      dimension isc(4),gn(4,4),gv(4,4),hp(4,4),hf(4,4)
      double precision vgasp(4),vgasf(4),B1(5),B2(5),
     &     BBp(4),BBf(4),BBc(4),vradp(4),vradf(4),s(5),
     &     ugasconf(4),ugascovf(4),ugasconp(4),ugascovp(4),
     &     uradconf(4),uradcovf(4),uradconp(4),uradcovp(4),
     &     ugasconi(4),ugascovi(4),uradconi(4),uradcovi(4),
     &     bconp(4),bcovp(4),bconf(4),bcovf(4),bconi(4),bcovi(4),
     &     Tmunup(4,4),Tmunuf(4,4),Rmunup(4,4),Rmunuf(4,4),
     &     Tmunui(4,4),Rmunui(4,4),Ttcovi(4),Rtcovi(4),
     &     Ttcovc(4),Rtcovc(4),Ttotc(4),prim(4),error(4),
     &     primeng(4),priment(4),
     &     ugasconfinal(4),ugascovfinal(4),uradconfinal(4),
     &     uradcovfinal(4),bconfinal(4),bcovfinal(4),
     &     Ttcovfinal(4),Rtcovfinal(4),
     &     Tmunufinal(4,4),Rmunufinal(4,4),rhou0c,Gam
      double precision error0(4),errornorm(4),err4,err3
      integer iflag,jflag
      integer fixcons
      double precision turadconfinal(4)
      double precision tugasconfinal(4)

#if(PRODUCTION==0)
      write (*,*)
      write (*,*)
      write (*,*) ' GETTING FINAL SOLUTION: ',which,showboth
#endif

      ufinal=prim(1)
      call solveucon(prim,ugasconfinal)
      call contocov(ugasconfinal,ugascovfinal)
      rhou0final=rhou0c
      rhofinal=rhou0c/ugasconfinal(1)

#if(PRODUCTION==0)
c     write (*,*)
c     write (*,*) ' final solution: rho, u, u^mu: ',
c     &        rhofinal,ufinal,(ugasconfinal(j),j=1,4)
#endif
c     Update T^mu_nu so as to be consistent with the final
c     primitives. Adjust R^t_mu so that the total energy and momentum
c     density are unchaned, then calculate the full R^mu_nu

      do j=1,4
         Ttotc(j)=Ttcovc(j)+Rtcovc(j)
      enddo

      call calcbconbcov(BBc,ugasconfinal,ugascovfinal,
     &     bconfinal,bcovfinal,bsqfinal)

      call calcTmunu(rhofinal,ufinal,bsqfinal,Gam,ugasconfinal,
     &     ugascovfinal,bconfinal,bcovfinal,Tmunufinal)

c     whether to fix cons or not
      do j=1,4
         Ttcovfinal(j)=Tmunufinal(1,j)
         Rtcovfinal(j)=Ttotc(j)-Ttcovfinal(j)
      enddo

      call Rmunuinvert(Rtcovfinal,Efinal,uradconfinal,
     &     uradcovfinal,ugasconfinal,Rmunufinal,radinvmod)

      alpha=1.d0/sqrt(-gn(1,1))
      gammagas=ugasconfinal(1)*alpha
      gammarad=uradconfinal(1)*alpha
#if(PRODUCTION==0)
      write (*,*)
      write (*,*)
      write (*,*) ' FINAL SOLUTION: ',which,showboth
      write (*,*)
      write (*,*) ' conserved rho u^0: ',rhou0final
      write (*,*) ' orig conserved T^t_mu: ',(Ttcovc(j),j=1,4)
      write (*,*) ' new  conserved T^t_mu: ',(Ttcovfinal(j),j=1,4)
      write (*,*) ' orig conserved R^t_mu: ',(Rtcovc(j),j=1,4)
      write (*,*) ' new  conserved R^t_mu: ',(Rtcovfinal(j),j=1,4)
      write (*,*) ' rho, u_g, E: ',rhofinal,ufinal,Efinal
      write (*,*) ' ugascon: ',(ugasconfinal(j),j=1,4)
      write (*,*) ' uradcon: ',(uradconfinal(j),j=1,4)
      write (*,*) ' gammagas, gammarad: ',gammagas,gammarad
      write (*,*) ' BBc, bsqfinal: ',BBc,bsqfinal
      write (*,*) ' Gam: ',Gam
c      Check on final error
      if(which.eq.3) then
         call funcrad1(prim,error0,errornorm,err4,err3,iflag,jflag)
         write (*,*) ' Final error: ',(errornorm(j),j=1,4),err4,err3
      endif
c      Check on final error
      if(which.eq.2) then
         call funcrad2(prim,error0,errornorm,err4,err3,iflag,jflag)
         write (*,*) ' Final error: ',(errornorm(j),j=1,4),err4,err3
      endif
#endif

      if(WHICHVEL.eq.VELREL4) then
c     get \tilde{u}^\mu_{\rm rad}
         call uconrel(uradconfinal,turadconfinal)
c     get \tilde{u}^\mu_{\rm gas}
         call uconrel(ugasconfinal,tugasconfinal)
      endif


#if(PRODUCTION<=1)
      if(showboth.eq.0.and.which.eq.2) then
         write (13,"(18X)",advance="no")
      endif
      if(WHICHVEL.eq.VEL4) then
      write (13,*) 'RAMESH(type) rho, u_g, u^mu Erf urad^mu: '
     &     ,which,rhofinal,
     &     ufinal,
     &     (ugasconfinal(j),j=1,4),
     &     Efinal,
     &     (uradconfinal(j),j=1,4)
      else if(WHICHVEL.eq.VELREL4) then
      write (13,*) 'RAMESH(type) rho, u_g, tu^mu Erf turad^mu: '
     &     ,which,rhofinal,
     &     ufinal,
     &     (tugasconfinal(j),j=1,4),
     &     Efinal,
     &     (turadconfinal(j),j=1,4)
      else
         stop
      endif


#endif

c     HARM order
      results(1) = rhofinal
      results(2) = ufinal
      if(WHICHVEL.eq.VEL4) then
         results(3) = ugasconfinal(1)
         results(4) = ugasconfinal(2)
         results(5) = ugasconfinal(3)
         results(6) = ugasconfinal(4)
      else if(WHICHVEL.eq.VELREL4) then
         results(3) = tugasconfinal(1)
         results(4) = tugasconfinal(2)
         results(5) = tugasconfinal(3)
         results(6) = tugasconfinal(4)
      else
         stop
      endif
      results(7) = Efinal
      if(WHICHVEL.eq.VEL4) then
         results(8) = uradconfinal(1)
         results(9) = uradconfinal(2)
         results(10) = uradconfinal(3)
         results(11) = uradconfinal(4)
      else if(WHICHVEL.eq.VELREL4) then
         results(8) = turadconfinal(1)
         results(9) = turadconfinal(2)
         results(10) = turadconfinal(3)
         results(11) = turadconfinal(4)
      else
         stop
      endif         
      results(15) = DBLE(radinvmod)

c     DEBUG TEST
c      write (13,*) 'RAMESH2: '
c     &     ,results(1),results(2),results(3),results(4),results(5)
c     &     ,results(6),results(7),results(8),results(9),results(10)
c     &     ,results(11)
      

      return
      end





cccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine readtype1(gn,gv,rhof,rhop,rhou0i,rhou0f,rhou0p,
     &        uf,up,T00i,T00f,T00p,vgasf,vgasp,T01i,T01f,T01p,
     &        T02i,T02f,T02p,T03i,T03f,T03p,BBf,BBp,
     &        Ef,Ep,R00i,R00f,R00p,vradf,vradp,R01i,R01f,R01p,
     &        R02i,R02f,R02p,R03i,R03f,R03p,s,si,sf,sp,
     &        uradconf,uradcovf,ugasconf,ugascovf,
     &        uradconp,uradcovp,ugasconp,ugascovp,ifinish,errorabs)

c     Read in data in Jon's old format (134 numbers)

      implicit double precision (a-h,o-z)
      dimension isc(4),gn(4,4),gv(4,4),vgasp(4),vgasf(4),
     &     BBp(4),BBf(4),vradp(4),vradf(4),s(5),
     &     ugasconf(4),ugascovf(4),ugasconp(4),ugascovp(4),
     &     uradconf(4),uradcovf(4),uradconp(4),uradcovp(4)

      ifinish=0

#if(PRODUCTION<=2)
         read (11,*,end=10) (isc(j),j=1,4),dt,
     &        ((gn(j,k),k=1,4),j=1,4),
     &        ((gv(j,k),k=1,4),j=1,4),
     &        rhof,rhop,rhou0i,rhou0f,rhou0p,
     &        uf,up,T00i,T00f,T00p,
     &        vgasf(2),vgasp(2),T01i,T01f,T01p,
     &        vgasf(3),vgasp(3),T02i,T02f,T02p,
     &        vgasf(4),vgasp(4),T03i,T03f,T03p,
     &        BBf(2),BBp(2),scr,scr,scr,
     &        BBf(3),BBp(3),scr,scr,scr,
     &        BBf(4),BBp(4),scr,scr,scr,
     &        Ef,Ep,R00i,R00f,R00p,
     &        vradf(2),vradp(2),R01i,R01f,R01p,
     &        vradf(3),vradp(3),R02i,R02f,R02p,
     &        vradf(4),vradp(4),R03i,R03f,R03p,
     &        (s(j),j=1,2),si,sf,sp,
     &        (uradconf(j),uradcovf(j),j=1,4),
     &        (ugasconf(j),ugascovf(j),j=1,4),
     &        (uradconp(j),uradcovp(j),j=1,4),
     &        (ugasconp(j),ugascovp(j),j=1,4)
#endif

#if(PRODUCTION==0)
c         write (*,*) (isc(j),j=1,4),dt,
c     &        ((gn(j,k),k=1,4),j=1,4),
c     &        ((gv(j,k),k=1,4),j=1,4),
c     &        rhof,rhop,rhou0i,rhou0f,rhou0p,
c     &        uf,up,T00i,T00f,T00p,
c     &        vgasf(2),vgasp(2),T01i,T01f,T01p,
c     &        vgasf(3),vgasp(3),T02i,T02f,T02p,
c     &        vgasf(4),vgasp(4),T03i,T03f,T03p,
c     &        BBf(2),BBp(2),scr,scr,scr,
c     &        BBf(3),BBp(3),scr,scr,scr,
c     &        BBf(4),BBp(4),scr,scr,scr,
c     &        Ef,Ep,R00i,R00f,R00p,
c     &        vradf(2),vradp(2),R01i,R01f,R01p,
c     &        vradf(3),vradp(3),R02i,R02f,R02p,
c     &        vradf(4),vradp(4),R03i,R03f,R03p,
c     &        (s(j),j=1,2),si,sf,sp,
c     &        (uradconf(j),uradcovf(j),j=1,4),
c     &        (ugasconf(j),ugascovf(j),j=1,4),
c     &        (uradconp(j),uradcovp(j),j=1,4),
c     &        (ugasconp(j),ugascovp(j),j=1,4)
c         write (*,*) ' rho: ',rhof,rhop,rhou0i,rhou0f,rhou0p
#endif
         errorabs=1.0
         errorabsbestexternal=1.0
         tolallow=1E-8
c      set \tilde{u}^t=0
         vgasf(1)=0.0d0
         vgasp(1)=0.0d0
         vradf(1)=0.0d0
         vradp(1)=0.0d0
         return

 10   ifinish=1

      return
      end



      subroutine readtype2(gn,gv,rhof,rhop,rhou0i,rhou0f,rhou0p,
     &        uf,up,T00i,T00f,T00p,vgasf,vgasp,T01i,T01f,T01p,
     &        T02i,T02f,T02p,T03i,T03f,T03p,BBf,BBp,
     &        Ef,Ep,R00i,R00f,R00p,vradf,vradp,R01i,R01f,R01p,
     &        R02i,R02f,R02p,R03i,R03f,R03p,s,si,sf,sp,
     &        uradconf,uradcovf,ugasconf,ugascovf,
     &        uradconp,uradcovp,ugasconp,ugascovp,ifinish,errorabs)

c     Read in data in Jon's new format (181 numbers)

      implicit double precision (a-h,o-z)
      dimension isc(4),gn(4,4),gv(4,4),vgasp(4),vgasf(4),
     &     BBp(4),BBf(4),vradp(4),vradf(4),s(5),
     &     ugasconf(4),ugascovf(4),ugasconp(4),ugascovp(4),
     &     uradconf(4),uradcovf(4),uradconp(4),uradcovp(4),
     &     ugasconb(4),ugascovb(4),uradconb(4),uradcovb(4)
      dimension Ttcovc(4),Rtcovc(4),BBc(4)
      common/conserved/Gam,Gam1,en,en1,rhou0c,sc,Ttcovc,Rtcovc,
     &     BBc,dt,
     &     FNUMARGSHARM,FNUMRESULTSHARM,WHICHVELRAMESHHARM,
     &     GAMMAMAXRAD,ERADLIMIT,toltry,tolallow,
     &     ARAD_CODE,OKAPPA_ES_CODE11,OKAPPA_FF_CODE11,OKAPPA_BF_CODE11
      integer eomtype

      ifinish=0

#if(PRODUCTION<=2)
         read (11,*,end=10) (isc(j),j=1,4),
     &        errorabs,iters,dt,nstep,steppart,Gam,
     &        ((gn(j,k),k=1,4),j=1,4),
     &        ((gv(j,k),k=1,4),j=1,4),
     &        rhof,scr,rhob,rhop,rhou0i,rhou0f,rhou0p,
     &        uf,scr,ub,up,T00i,T00f,T00p,
     &        vgasf(2),scr,scr,vgasp(2),T01i,T01f,T01p,
     &        vgasf(3),scr,scr,vgasp(3),T02i,T02f,T02p,
     &        vgasf(4),scr,scr,vgasp(4),T03i,T03f,T03p,
     &        BBf(2),scr,scr,BBp(2),scr,scr,scr,
     &        BBf(3),scr,scr,BBp(3),scr,scr,scr,
     &        BBf(4),scr,scr,BBp(4),scr,scr,scr,
     &        Ef,scr,Eb,Ep,R00i,R00f,R00p,
     &        vradf(2),scr,scr,vradp(2),R01i,R01f,R01p,
     &        vradf(3),scr,scr,vradp(3),R02i,R02f,R02p,
     &        vradf(4),scr,scr,vradp(4),R03i,R03f,R03p,
     &        s(1),scr,scr,s(2),si,sf,sp,
     &        (uradconf(j),uradcovf(j),j=1,4),
     &        (ugasconf(j),ugascovf(j),j=1,4),
     &        (uradconb(j),uradcovb(j),j=1,4),
     &        (ugasconb(j),ugascovb(j),j=1,4),
     &        (uradconp(j),uradcovp(j),j=1,4),
     &        (ugasconp(j),ugascovp(j),j=1,4)
#endif
         eomtype=0
         tolallow=1E-8
c      set \tilde{u}^t=0
         vgasf(1)=0.0d0
         vgasp(1)=0.0d0
         vradf(1)=0.0d0
         vradp(1)=0.0d0

#if(PRODUCTION<=2)
         write (14,"(2X,1I5,2X,1E21.15,2X,1I8)",advance="no")
     & iters,errorabs,0
#endif
         if(errorabs.lt.FAILLIMIT) then
#if(PRODUCTION<=1)
            write (13,"(1X,A,9X)",advance="no") ' GOOD   '
#endif
#if(PRODUCTION<=2)
            write (14,"(1X,A)",advance="no") ' GOOD   '
#endif
         else
#if(PRODUCTION<=1)
            write (13,"(1X,A,9X)",advance="no") '  BAD   '
#endif
#if(PRODUCTION<=2)
            write (14,"(1X,A)",advance="no") '  BAD   '
#endif
         endif
        
#if(PRODUCTION<=1)
      if(WHICHVEL.eq.VEL4) then
         write (13,*) 'HARMJM(type) rho, u_g, u^mu Erf urad^mu: '
     &        ,eomtype
     &        ,rhof,uf,
     &        (ugasconf(j),j=1,4),
     &        Ef,
     &        (uradconf(j),j=1,4)
      else if(WHICHVEL.eq.VELREL4) then
         write (13,*) 'HARMJM(type) rho, u_g, tu^mu Erf turad^mu: '
     &        ,eomtype
     &        ,rhof,uf,
     &        (vgasf(j),j=1,4),
     &        Ef,
     &        (vradf(j),j=1,4)
      else
         stop
      endif
#endif
#if(PRODUCTION==0)
c         write (*,*) ' gn: ',((gn(i,j),j=1,4),i=1,4)
c         write (*,*) ' gv: ',((gv(i,j),j=1,4),i=1,4)
c         write (*,*) ' rhou0i, rhou0f, rhou0p: ',rhou0i,rhou0f,rhou0p
c         write (*,*) ' uf, up, T001, T00f, T00p: ',uf,up,T00i,T00f,T00p
c         write (*,*) ' ugasconp: ',(ugasconp(j),j=1,4)
c         write (*,*) ' ugasconf: ',(ugasconf(j),j=1,4)
c         write (*,*) ' BBp: ',(BBp(j),j=2,4)
c         write (*,*) ' BBf: ',(BBf(j),j=2,4)
         write (*,*) ' Ef, Ep, R00i, R00f, R00p: ',Ef,Ep,R00i,R00f,R00p
c         write (*,*) ' uradconp: ',(uradconp(j),j=1,4)
c         write (*,*) ' uradconf: ',(uradconf(j),j=1,4)
         write (*,*) ' si, sf, sp: ',si,sf,sp

c         write (*,*) (isc(j),j=1,4),dt,
c     &        ((gn(j,k),k=1,4),j=1,4),
c     &        ((gv(j,k),k=1,4),j=1,4),
c     &        rhof,rhop,rhou0i,rhou0f,rhou0p,
c     &        uf,up,T00i,T00f,T00p,
c     &        vgasf(2),vgasp(2),T01i,T01f,T01p,
c     &        vgasf(3),vgasp(3),T02i,T02f,T02p,
c     &        vgasf(4),vgasp(4),T03i,T03f,T03p,
c     &        BBf(2),BBp(2),scr,scr,scr,
c     &        BBf(3),BBp(3),scr,scr,scr,
c     &        BBf(4),BBp(4),scr,scr,scr,
c     &        Ef,Ep,R00i,R00f,R00p,
c     &        vradf(2),vradp(2),R01i,R01f,R01p,
c     &        vradf(3),vradp(3),R02i,R02f,R02p,
c     &        vradf(4),vradp(4),R03i,R03f,R03p,
c     &        (s(j),j=1,2),si,sf,sp,
c     &        (uradconf(j),uradcovf(j),j=1,4),
c     &        (ugasconf(j),ugascovf(j),j=1,4),
c     &        (uradconp(j),uradcovp(j),j=1,4),
c     &        (ugasconp(j),ugascovp(j),j=1,4)
#endif
         return

 10   ifinish=1

      return
      end


      subroutine readtype3(origintype,args,
     &        gn,gv,rhof,rhop,rhou0i,rhou0f,rhou0p,
     &        uf,up,T00i,T00f,T00p,vgasf,vgasp,T01i,T01f,T01p,
     &        T02i,T02f,T02p,T03i,T03f,T03p,BBf,BBp,
     &        Ef,Ep,R00i,R00f,R00p,vradf,vradp,R01i,R01f,R01p,
     &        R02i,R02f,R02p,R03i,R03f,R03p,s,si,sf,sp,
     &        uradconf,uradcovf,ugasconf,ugascovf,
     &        uradconp,uradcovp,ugasconp,ugascovp,ifinish,errorabs)

c     Read in data in Jon's new format (211+11 numbers)

      implicit double precision (a-h,o-z)
      integer orgintype
      dimension args(NUMARGS)
      integer eomtype,itermode,iters,totaliters
      dimension isc(4),gn(4,4),gv(4,4),vgasp(4),vgasf(4),
     &     BBp(4),BBf(4),vradp(4),vradf(4),s(5),
     &     ugasconf(4),ugascovf(4),ugasconp(4),ugascovp(4),
     &     uradconf(4),uradcovf(4),uradconp(4),uradcovp(4),
     &     ugasconb(4),ugascovb(4),uradconb(4),uradcovb(4)
      dimension Ttcovc(4),Rtcovc(4),BBc(4)
      common/conserved/Gam,Gam1,en,en1,rhou0c,sc,Ttcovc,Rtcovc,
     &     BBc,dt,
     &     FNUMARGSHARM,FNUMRESULTSHARM,WHICHVELRAMESHHARM,
     &     GAMMAMAXRAD,ERADLIMIT,toltry,tolallow,
     &     ARAD_CODE,OKAPPA_ES_CODE11,OKAPPA_FF_CODE11,OKAPPA_BF_CODE11
      integer na

      ifinish=0

#if(PRODUCTION<=2)
c     read (11,*,end=10) (isc(j),j=1,4),
         read (11,*,end=10) ,
     &     FNUMARGSHARM,FNUMRESULTSHARM,WHICHVELRAMESHHARM,
     &        failtype,myid,failnum,gotfirstnofail,
     &        eomtype,itermode,
     &        errorabs,errorabsbestexternal,iters,totaliters,
     &        dt,nstep,steppart,Gam,
     &        GAMMAMAXRAD,ERADLIMIT,toltry,tolallow,
     &        ARAD_CODE,
     &        OKAPPA_ES_CODE11,OKAPPA_FF_CODE11,OKAPPA_BF_CODE11,
     &        ((gn(j,k),k=1,4),j=1,4),
     &        ((gv(j,k),k=1,4),j=1,4),
     &        rhof,scr,rhob,rhop,src,src,rhou0i,rhou0f,rhou0p,
     &        uf,scr,ub,up,src,src,T00i,T00f,T00p,
     &        vgasf(2),scr,scr,vgasp(2),src,src,T01i,T01f,T01p,
     &        vgasf(3),scr,scr,vgasp(3),src,src,T02i,T02f,T02p,
     &        vgasf(4),scr,scr,vgasp(4),src,src,T03i,T03f,T03p,
     &        BBf(2),scr,scr,BBp(2),src,src,scr,scr,scr,
     &        BBf(3),scr,scr,BBp(3),src,src,scr,scr,scr,
     &        BBf(4),scr,scr,BBp(4),src,src,scr,scr,scr,
     &        Ef,scr,Eb,Ep,src,src,R00i,R00f,R00p,
     &        vradf(2),scr,scr,vradp(2),src,src,R01i,R01f,R01p,
     &        vradf(3),scr,scr,vradp(3),src,src,R02i,R02f,R02p,
     &        vradf(4),scr,scr,vradp(4),src,src,R03i,R03f,R03p,
     &        s(1),scr,scr,s(2),src,src,si,sf,sp,
     &        (uradconf(j),uradcovf(j),j=1,4),
     &        (ugasconf(j),ugascovf(j),j=1,4),
     &        (uradconb(j),uradcovb(j),j=1,4),
     &        (ugasconb(j),ugascovb(j),j=1,4),
     &        (uradconp(j),uradcovp(j),j=1,4),
     &        (ugasconp(j),ugascovp(j),j=1,4)
#endif
#if(PRODUCTION==3)
c     then pull from array of doubles from args
         na=0
c
         na=na+1
         FNUMARGSHARM=args(na)
         na=na+1
         FNUMRESULTSHARM=args(na)
         na=na+1
         WHICHVELRAMESHHARM=args(na)
         na=na+1
c
         failtype=args(na)
         na=na+1
         myid=args(na)
         na=na+1
         failnum=args(na)
         na=na+1
         gotfirstnofail=args(na)
         na=na+1
         eomtype=args(na)
         na=na+1
         itermode=args(na)
         na=na+1
         errorabs=args(na)
         na=na+1
         errorabsbestexternal=args(na)
         na=na+1
         iters=args(na)
         na=na+1
         totaliters=args(na)
         na=na+1
         dt=args(na)
         na=na+1
         nstep=args(na)
         na=na+1
         steppart=args(na)
         na=na+1
         Gam=args(na)
         na=na+1
c
         GAMMAMAXRAD=args(na)
         na=na+1
         ERADLIMIT=args(na)
         na=na+1
         toltry=args(na)
         na=na+1
         tolallow=args(na)
         na=na+1
         ARAD_CODE=args(na)
         na=na+1
         OKAPPA_ES_CODE11=args(na)
         na=na+1
         OKAPPA_FF_CODE11=args(na)
         na=na+1
         OKAPPA_BF_CODE11=args(na)
c
         na=na+1
         gn(1,1)=args(na)
         na=na+1
         gn(1,2)=args(na)
         na=na+1
         gn(1,3)=args(na)
         na=na+1
         gn(1,4)=args(na)
         na=na+1
         gn(2,1)=args(na)
         na=na+1
         gn(2,2)=args(na)
         na=na+1
         gn(2,3)=args(na)
         na=na+1
         gn(2,4)=args(na)
         na=na+1
         gn(3,1)=args(na)
         na=na+1
         gn(3,2)=args(na)
         na=na+1
         gn(3,3)=args(na)
         na=na+1
         gn(3,4)=args(na)
         na=na+1
         gn(4,1)=args(na)
         na=na+1
         gn(4,2)=args(na)
         na=na+1
         gn(4,3)=args(na)
         na=na+1
         gn(4,4)=args(na)
         na=na+1
         gv(1,1)=args(na)
         na=na+1
         gv(1,2)=args(na)
         na=na+1
         gv(1,3)=args(na)
         na=na+1
         gv(1,4)=args(na)
         na=na+1
         gv(2,1)=args(na)
         na=na+1
         gv(2,2)=args(na)
         na=na+1
         gv(2,3)=args(na)
         na=na+1
         gv(2,4)=args(na)
         na=na+1
         gv(3,1)=args(na)
         na=na+1
         gv(3,2)=args(na)
         na=na+1
         gv(3,3)=args(na)
         na=na+1
         gv(3,4)=args(na)
         na=na+1
         gv(4,1)=args(na)
         na=na+1
         gv(4,2)=args(na)
         na=na+1
         gv(4,3)=args(na)
         na=na+1
         gv(4,4)=args(na)
         na=na+1
         rhof=args(na)
         na=na+1
         scr=args(na)
         na=na+1
         rhob=args(na)
         na=na+1
         rhop=args(na)
         na=na+1
         src=args(na)
         na=na+1
         src=args(na)
         na=na+1
         rhou0i=args(na)
         na=na+1
         rhou0f=args(na)
         na=na+1
         rhou0p=args(na)
         na=na+1
         uf=args(na)
         na=na+1
         scr=args(na)
         na=na+1
         ub=args(na)
         na=na+1
         up=args(na)
         na=na+1
         src=args(na)
         na=na+1
         src=args(na)
         na=na+1
         T00i=args(na)
         na=na+1
         T00f=args(na)
         na=na+1
         T00p=args(na)
         na=na+1
         vgasf(2)=args(na)
         na=na+1
         scr=args(na)
         na=na+1
         scr=args(na)
         na=na+1
         vgasp(2)=args(na)
         na=na+1
         src=args(na)
         na=na+1
         src=args(na)
         na=na+1
         T01i=args(na)
         na=na+1
         T01f=args(na)
         na=na+1
         T01p=args(na)
         na=na+1
         vgasf(3)=args(na)
         na=na+1
         scr=args(na)
         na=na+1
         scr=args(na)
         na=na+1
         vgasp(3)=args(na)
         na=na+1
         src=args(na)
         na=na+1
         src=args(na)
         na=na+1
         T02i=args(na)
         na=na+1
         T02f=args(na)
         na=na+1
         T02p=args(na)
         na=na+1
         vgasf(4)=args(na)
         na=na+1
         scr=args(na)
         na=na+1
         scr=args(na)
         na=na+1
         vgasp(4)=args(na)
         na=na+1
         src=args(na)
         na=na+1
         src=args(na)
         na=na+1
         T03i=args(na)
         na=na+1
         T03f=args(na)
         na=na+1
         T03p=args(na)
         na=na+1
         BBf(2)=args(na)
         na=na+1
         scr=args(na)
         na=na+1
         scr=args(na)
         na=na+1
         BBp(2)=args(na)
         na=na+1
         src=args(na)
         na=na+1
         src=args(na)
         na=na+1
         scr=args(na)
         na=na+1
         scr=args(na)
         na=na+1
         scr=args(na)
         na=na+1
         BBf(3)=args(na)
         na=na+1
         scr=args(na)
         na=na+1
         scr=args(na)
         na=na+1
         BBp(3)=args(na)
         na=na+1
         src=args(na)
         na=na+1
         src=args(na)
         na=na+1
         scr=args(na)
         na=na+1
         scr=args(na)
         na=na+1
         scr=args(na)
         na=na+1
         BBf(4)=args(na)
         na=na+1
         scr=args(na)
         na=na+1
         scr=args(na)
         na=na+1
         BBp(4)=args(na)
         na=na+1
         src=args(na)
         na=na+1
         src=args(na)
         na=na+1
         scr=args(na)
         na=na+1
         scr=args(na)
         na=na+1
         scr=args(na)
         na=na+1
         Ef=args(na)
         na=na+1
         scr=args(na)
         na=na+1
         Eb=args(na)
         na=na+1
         Ep=args(na)
         na=na+1
         src=args(na)
         na=na+1
         src=args(na)
         na=na+1
         R00i=args(na)
         na=na+1
         R00f=args(na)
         na=na+1
         R00p=args(na)
         na=na+1
         vradf(2)=args(na)
         na=na+1
         scr=args(na)
         na=na+1
         scr=args(na)
         na=na+1
         vradp(2)=args(na)
         na=na+1
         src=args(na)
         na=na+1
         src=args(na)
         na=na+1
         R01i=args(na)
         na=na+1
         R01f=args(na)
         na=na+1
         R01p=args(na)
         na=na+1
         vradf(3)=args(na)
         na=na+1
         scr=args(na)
         na=na+1
         scr=args(na)
         na=na+1
         vradp(3)=args(na)
         na=na+1
         src=args(na)
         na=na+1
         src=args(na)
         na=na+1
         R02i=args(na)
         na=na+1
         R02f=args(na)
         na=na+1
         R02p=args(na)
         na=na+1
         vradf(4)=args(na)
         na=na+1
         scr=args(na)
         na=na+1
         scr=args(na)
         na=na+1
         vradp(4)=args(na)
         na=na+1
         src=args(na)
         na=na+1
         src=args(na)
         na=na+1
         R03i=args(na)
         na=na+1
         R03f=args(na)
         na=na+1
         R03p=args(na)
         na=na+1
         s(1)=args(na)
         na=na+1
         scr=args(na)
         na=na+1
         scr=args(na)
         na=na+1
         s(2)=args(na)
         na=na+1
         src=args(na)
         na=na+1
         src=args(na)
         na=na+1
         si=args(na)
         na=na+1
         sf=args(na)
         na=na+1
         sp=args(na)
c con/cov in 1,2,3,4
         na=na+1
         uradconf(1)=args(na)
         na=na+1
         uradcovf(1)=args(na)
         na=na+1
         uradconf(2)=args(na)
         na=na+1
         uradcovf(2)=args(na)
         na=na+1
         uradconf(3)=args(na)
         na=na+1
         uradcovf(3)=args(na)
         na=na+1
         uradconf(4)=args(na)
         na=na+1
         uradcovf(4)=args(na)
c
         na=na+1
         ugasconf(1)=args(na)
         na=na+1
         ugascovf(1)=args(na)
         na=na+1
         ugasconf(2)=args(na)
         na=na+1
         ugascovf(2)=args(na)
         na=na+1
         ugasconf(3)=args(na)
         na=na+1
         ugascovf(3)=args(na)
         na=na+1
         ugasconf(4)=args(na)
         na=na+1
         ugascovf(4)=args(na)
c
         na=na+1
         uradconb(1)=args(na)
         na=na+1
         uradcovb(1)=args(na)
         na=na+1
         uradconb(2)=args(na)
         na=na+1
         uradcovb(2)=args(na)
         na=na+1
         uradconb(3)=args(na)
         na=na+1
         uradcovb(3)=args(na)
         na=na+1
         uradconb(4)=args(na)
         na=na+1
         uradcovb(4)=args(na)
c
         na=na+1
         ugasconb(1)=args(na)
         na=na+1
         ugascovb(1)=args(na)
         na=na+1
         ugasconb(2)=args(na)
         na=na+1
         ugascovb(2)=args(na)
         na=na+1
         ugasconb(3)=args(na)
         na=na+1
         ugascovb(3)=args(na)
         na=na+1
         ugasconb(4)=args(na)
         na=na+1
         ugascovb(4)=args(na)
c
         na=na+1
         uradconp(1)=args(na)
         na=na+1
         uradcovp(1)=args(na)
         na=na+1
         uradconp(2)=args(na)
         na=na+1
         uradcovp(2)=args(na)
         na=na+1
         uradconp(3)=args(na)
         na=na+1
         uradcovp(3)=args(na)
         na=na+1
         uradconp(4)=args(na)
         na=na+1
         uradcovp(4)=args(na)
c
         na=na+1
         ugasconp(1)=args(na)
         na=na+1
         ugascovp(1)=args(na)
         na=na+1
         ugasconp(2)=args(na)
         na=na+1
         ugascovp(2)=args(na)
         na=na+1
         ugasconp(3)=args(na)
         na=na+1
         ugascovp(3)=args(na)
         na=na+1
         ugasconp(4)=args(na)
         na=na+1
         ugascovp(4)=args(na)
c
         if(na.ne.NUMARGS) then
#if(PRODUCTION==0)
            write(*,*) 'Wrong number of args'
#endif
            stop
         endif
#endif


c      set \tilde{u}^t=0
         vgasf(1)=0.0d0
         vgasp(1)=0.0d0
         vradf(1)=0.0d0
         vradp(1)=0.0d0


c     HARMEOMTYPE=2 is entropy, 3 is energy.  If 3 fails, could have reverted to entropy.
c     If itermode=0, then in default harm mode, this was not the last attempt for the solver, so just looking at why this strategy failed -- not failure of harm ultimately.
c     If itermode=1 shows up, then check whether PREVBESTHARMERR was ok/good enough even if not <tol.  If itermode=1 and both current and best error is bad, actually BAD case for harm.
#if(PRODUCTION==0)
         write (*,"(A,1X,1F21.15)") 'TEST',ugascovp(4)
#endif
#if(PRODUCTION<=2)
         write (14,"(1I5,2X,1I5,5X,1E21.15,
     &         2X,1E21.15,2X,1I8,8X,1I1,8X)"
     &        ,advance="no")
     &        iters,totaliters,errorabs,errorabsbestexternal
     &        ,eomtype,itermode
#endif
         if(errorabs.lt.FAILLIMIT) then
#if(PRODUCTION<=1)
            write (13,"(1X,A,9X)",advance="no") ' GOOD   '
#endif
#if(PRODUCTION<=2)
            write (14,"(1X,A)",advance="no") ' GOOD   '
#endif
         else
#if(PRODUCTION<=1)
            write (13,"(1X,A,9X)",advance="no") '  BAD   '
#endif
#if(PRODUCTION<=2)
            write (14,"(1X,A)",advance="no") '  BAD   '
#endif
         endif
         
#if(PRODUCTION<=1)
      if(WHICHVEL.eq.VEL4) then
         write (13,*) 'HARMJM(type) rho, u_g, u^mu Erf urad^mu: '
     &        ,eomtype
     &        ,rhof,uf,
     &        (ugasconf(j),j=1,4),
     &        Ef,
     &        (uradconf(j),j=1,4)
      else if(WHICHVEL.eq.VELREL4) then
         write (13,*) 'HARMJM(type) rho, u_g, tu^mu Erf turad^mu: '
     &        ,eomtype
     &        ,rhof,uf,
     &        (vgasf(j),j=1,4),
     &        Ef,
     &        (vradf(j),j=1,4)
      else
         stop
      endif
#endif

#if(PRODUCTION==0)
c         write (*,*) ' gn: ',((gn(i,j),j=1,4),i=1,4)
c         write (*,*) ' gv: ',((gv(i,j),j=1,4),i=1,4)
c         write (*,*) ' rhou0i, rhou0f, rhou0p: ',rhou0i,rhou0f,rhou0p
c         write (*,*) ' uf, up, T001, T00f, T00p: ',uf,up,T00i,T00f,T00p
c         write (*,*) ' ugasconp: ',(ugasconp(j),j=1,4)
c         write (*,*) ' ugasconf: ',(ugasconf(j),j=1,4)
c         write (*,*) ' BBp: ',(BBp(j),j=2,4)
c         write (*,*) ' BBf: ',(BBf(j),j=2,4)
         write (*,*) ' Ef, Ep, R00i, R00f, R00p: ',Ef,Ep,R00i,R00f,R00p
c         write (*,*) ' uradconp: ',(uradconp(j),j=1,4)
c         write (*,*) ' uradconf: ',(uradconf(j),j=1,4)
         write (*,*) ' si, sf, sp: ',si,sf,sp

c         write (*,*) (isc(j),j=1,4),dt,
c     &        ((gn(j,k),k=1,4),j=1,4),
c     &        ((gv(j,k),k=1,4),j=1,4),
c     &        rhof,rhop,rhou0i,rhou0f,rhou0p,
c     &        uf,up,T00i,T00f,T00p,
c     &        vgasf(2),vgasp(2),T01i,T01f,T01p,
c     &        vgasf(3),vgasp(3),T02i,T02f,T02p,
c     &        vgasf(4),vgasp(4),T03i,T03f,T03p,
c     &        BBf(2),BBp(2),scr,scr,scr,
c     &        BBf(3),BBp(3),scr,scr,scr,
c     &        BBf(4),BBp(4),scr,scr,scr,
c     &        Ef,Ep,R00i,R00f,R00p,
c     &        vradf(2),vradp(2),R01i,R01f,R01p,
c     &        vradf(3),vradp(3),R02i,R02f,R02p,
c     &        vradf(4),vradp(4),R03i,R03f,R03p,
c     &        (s(j),j=1,2),si,sf,sp,
c     &        (uradconf(j),uradcovf(j),j=1,4),
c     &        (ugasconf(j),ugascovf(j),j=1,4),
c     &        (uradconp(j),uradcovp(j),j=1,4),
c     &        (ugasconp(j),ugascovp(j),j=1,4)
#endif
         return

 10   ifinish=1

      return
      end


      

      double precision function mymax(x,y)
      double precision x,y

      mymax=max(x,y)
      if(x.ne.x) then
         mymax=x
      endif
      if(y.ne.y) then
         mymax=y
      endif

      if(myisnan(x).eq.1) then
         mymax=x
      endif
      if(myisnan(y).eq.1) then
         mymax=y
      endif

      return
      end


      double precision function myabs(x)
      double precision x

      myabs=abs(x)
      if(myisnan(x).eq.1) then
         myabs=x
      endif

      return
      end


      integer function myisnan(x)
      double precision x

c for C / harm
#if(PRODUCTION==3)
      if(x.ne.x) then
#else
c for fortran
      if(isnan(x)) then
#endif
         myisnan=1
      else
         myisnan=0
      endif

      return
      end


      double precision function mydiv(x,y)
      double precision x,y

      double precision myabs
c     Avoid division by zero
      mydiv=x*dsign(1.0d0,y)/(SMALL+myabs(y))
      if(x.ne.x) then
         mydiv=x
      endif
      if(y.ne.y) then
         mydiv=y
      endif

      if(myisnan(x).eq.1) then
         mydiv=x
      endif
      if(myisnan(y).eq.1) then
         mydiv=y
      endif

      return
      end



      subroutine solveucon (prim,con)

c     Calculates full u^\mu given the primitives

      implicit double precision (a-h,o-z)
      dimension con(4), tcon(4), gn(4,4),gv(4,4)
      dimension prim(4)
      common/metric/gn,gv

      if(WHICHVEL.eq.VEL4) then
         do j=2,4
            con(j)=prim(j)
         enddo
         call solveu0old(con)
      else if(WHICHVEL.eq.VELREL4) then
         tcon(1)=0.0d0
         do j=2,4
            tcon(j)=prim(j)
         enddo         
         call solveu0new(tcon,con)
      else
         stop
      endif

      return
      end


      subroutine uconrel (con,tcon)

c     Calculates \tilde{u}^i from full u^\mu

      implicit double precision (a-h,o-z)
      dimension con(4), tcon(4), gn(4,4),gv(4,4)
      common/metric/gn,gv

c     \alpha = 1/sqrt(-g^{tt})
      alphalapse = 1.0/sqrt(abs(-gn(1,1)))

      tcon(1)=0.0d0
      do j=2,4
         tcon(j) = con(j) + con(1)* (gn(1,j)*alphalapse*alphalapse)
      enddo

      return
      end


      subroutine solveu0old (con)

c     Calculates u^0 given the other three components of the 4-velocity,
c     u^1, u^2, u^3

      implicit double precision (a-h,o-z)
      dimension con(4), gn(4,4),gv(4,4)
      common/metric/gn,gv

      au0=gv(1,1)
      bu0=0.d0
      cu0=1.d0
      do j=2,4
         bu0=bu0+2.d0*gv(1,j)*con(j)
      do k=2,4
         cu0=cu0+gv(j,k)*con(j)*con(k)
      enddo
      enddo
         
      con(1)=(-bu0-sqrt(bu0*bu0-4.d0*au0*cu0))/(2.d0*au0)

      return
      end


      subroutine solveu0new (tcon,con)

c     Calculates u^\mu given \tilde{u}^i

      implicit double precision (a-h,o-z)
      double precision qsq,gamma,alphalapse
      dimension tcon(4),con(4), gn(4,4),gv(4,4)
      common/metric/gn,gv

c      q^2 = \tilde{u}^\mu \tilde{u}^\nu g_{\mu\nu} and \tilde{u}^t=0
      tcon(1)=0.0
      qsq = 0.0
      do j=2,4
         do k=2,4
            qsq = qsq + gv(j,k)*tcon(j)*tcon(k)
         enddo
      enddo

c     \gamma = \sqrt{1+q^2}
      gamma = sqrt(1.0+qsq)
c     \alpha = 1/sqrt(-g^{tt})
      alphalapse = 1.0/sqrt(abs(-gn(1,1)))

c     u^t = \gamma/\alpha
      con(1) = gamma/alphalapse

c     u^j = \tilde{u}^j - u^t \beta^j
      do j=2,4
         con(j) = tcon(j) - con(1)* (gn(1,j)*alphalapse*alphalapse)
      enddo


      return
      end



      subroutine contocov(vecin,vecout)

c     Converts a contravariant 4-vector to a covariant 4-vector

      implicit double precision (a-h,o-z)
      dimension vecin(4),vecout(4),gn(4,4),gv(4,4)
      common/metric/gn,gv

      do i=1,4
         vecout(i)=0.d0
      do j=1,4
         vecout(i)=vecout(i)+gv(i,j)*vecin(j)
      enddo
      enddo

      return
      end



      subroutine covtocon(vecin,vecout)

c     Converts a covariant 4-vector to a contravariant 4-vector

      implicit double precision (a-h,o-z)
      dimension vecin(4),vecout(4),gn(4,4),gv(4,4)
      common/metric/gn,gv

      do i=1,4
         vecout(i)=0.d0
      do j=1,4
         vecout(i)=vecout(i)+gn(i,j)*vecin(j)
      enddo
      enddo

      return
      end



      function condotcov(con,cov)

c     Calculates the dot product con^mu cov_mu

      implicit double precision (a-h,o-z)
      dimension con(4),cov(4),gn(4,4),gv(4,4)
      common/metric/gn,gv

      condotcov=0.d0
      do i=1,4
         condotcov=condotcov+con(i)*cov(i)
      enddo

      return
      end



      subroutine calcbconbcov(BB,ucon,ucov,bcon,bcov,bsq)

c     Given B(1), B(2), B(3), and gas 4-velocity, calculates the
c     four-vectors b^mu and b_mu, and bsq = b^mub_mu

      implicit double precision (a-h,o-z)
      dimension BB(4),ucon(4),ucov(4),bcon(4),bcov(4),gn(4,4),gv(4,4)
      common/metric/gn,gv

      BB(1)=0.d0
      udotB=condotcov(BB,ucov)
      do i=1,4
         bcon(i)=(BB(i)+udotB*ucon(i))/ucon(1)
      enddo

      call contocov(bcon,bcov)

      bsq=condotcov(bcon,bcov)

      return
      end



      subroutine calcTmunu(rho,u,bsq,Gam,ucon,ucov,
     &     bcon,bcov,Tmunu)

c     Given rho, u, 4-velocity, magnetic 4-vector, calculates the gas
c     stress energy tensor T^mu_nu

      implicit double precision (a-h,o-z)
      dimension ucon(4),ucov(4),bcon(4),bcov(4),Tmunu(4,4),
     &     gn(4,4),gv(4,4),delta(4,4)
      common/metric/gn,gv

      call makedelta(delta)

      p=(Gam-1.d0)*u
      bsq2=0.5d0*bsq

      
      do i=1,4
      do j=1,4

      if (i.eq.1.and.j.eq.1) then

c     Replace T^0_0 -> T^0_0 + rho*u^0                                             

         Tmunu(i,j)=rho*ucon(i)*(1.d0+ucov(j))
     &        +(u+p+bsq)*ucon(i)*ucov(j)
     &        +(p+bsq2)*delta(i,j)
     &        -bcon(i)*bcov(j)

      else

         Tmunu(i,j)=(rho+u+p+bsq)*ucon(i)*ucov(j)
     &        +(p+bsq2)*delta(i,j)
     &        -bcon(i)*bcov(j)

      endif

      enddo
      enddo


      return
      end



      subroutine calcRmunu(E,ucon,ucov,Rmunu)

c     Given the radiation frame energy density E and 4-velocity,
c     calculates the radiation stress-energy tensor R^mu_nu

      implicit double precision (a-h,o-z)
      dimension ucon(4),ucov(4),Rmunu(4,4),
     &     gn(4,4),gv(4,4),delta(4,4)
      common/metric/gn,gv

      call makedelta(delta)

      four3=4.d0/3.d0
      one3=1.d0/3.d0

      do i=1,4
      do j=1,4

         Rmunu(i,j)=four3*E*ucon(i)*ucov(j)
     &        +one3*E*delta(i,j)

      enddo
      enddo

      return
      end



      subroutine makedelta(delta)

c     Calculate Kronecker delta

      implicit double precision (a-h,o-z)
      dimension delta (4,4)

      do i=1,4
      do j=1,4

         if (i.eq.j) then
            delta(i,j)=1.d0
         else
            delta(i,j)=0.d0
         endif

      enddo
      enddo

      return
      end


      subroutine Rmunuinvert(Rtcov,E,ucon,ucov,ugascon,Rmunu,radinvmod)
      implicit double precision (a-h,o-z)
      integer radinvmod
      dimension Rtcov(4),Rtcon(4),ucon(4),ucov(4),Rmunu(4,4),
     &     ugascon(4),gn(4,4),gv(4,4)
      common/accuracy/eps,epsbis,dvmin,tol,uminfac,dlogmax,
     &     gammaradceiling,itermax
      common/metric/gn,gv

c      call Rmunuinvert0(Rtcov,E,ucon,ucov,ugascon,Rmunu,radinvmod)
      call Rmunuinvert1(Rtcov,E,ucon,ucov,ugascon,Rmunu,radinvmod)


      return
      end



c     Old routine
      subroutine Rmunuinvert0(Rtcov,E,ucon,ucov,ugascon,Rmunu,radinvmod)

c     Given the row R^t_mu of the radiation tensor, solves for the
c     radiation frame energy density E and 4-velocity and calculates the
c     full tensor R^mu_nu

      implicit double precision (a-h,o-z)
      integer radinvmod
      dimension Rtcov(4),Rtcon(4),ucon(4),ucov(4),Rmunu(4,4),
     &     ugascon(4),gn(4,4),gv(4,4)
      common/accuracy/eps,epsbis,dvmin,tol,uminfac,dlogmax,
     &     gammaradceiling,itermax
      common/metric/gn,gv

c     Convert R^0_mu to R^0^mu

      call covtocon(Rtcov,Rtcon)
#if(PRODUCTION==0)
c      write (*,*) ' Rtcov: ',(Rtcov(i),i=1,4)
c      write (*,*) ' Rtcon: ',(Rtcon(i),i=1,4)
#endif
c     Set up and solve quadratic equation for E

c     default is no ceiling hit on gamma
      radinvmod=0
c      
      aquad=gn(1,1)
      bquad=-2.d0*Rtcon(1)
      cquad=0.d0
      do i=1,4
      do j=1,4
         cquad=cquad-gv(i,j)*Rtcon(i)*Rtcon(j)
      enddo
      enddo
      cquad=3.d0*cquad
      disc=bquad*bquad-4.d0*aquad*cquad

c     Make sure signs are okay for a physical solution. If not, go to 10
c     for alternative calculation

      if (aquad*cquad.ge.0.d0.or.disc.lt.0.d0) then
#if(PRODUCTION==0)
         write(*,*) 'not physical: ',aquad,cquad,disc
#endif
         go to 10
      endif

c     Use negative sign of discriminant and solve for E and ucon

      E=2.d0*cquad/(-bquad+sqrt(disc))
#if(PRODUCTION==0)
c      write (*,*) ' negative sign of discriminant '
      write (*,*) ' a, b, c, disc, E: ',aquad,bquad,cquad,disc,E
#endif

c     Make sure E is positive. If not, go to 10

      if (E.lt.0.d0) then
#if(PRODUCTION==0)
         write(*,*) 'E is not positive: ',E
#endif
         go to 10
      endif

      ucon(1)=0.5d0*sqrt(3.d0*Rtcon(1)/E-gn(1,1))
      do i=2,4
         ucon(i)=(3.d0*Rtcon(i)-E*gn(1,i))/(4.d0*E*ucon(1))
      enddo

c     Make sure the Lorentz factor is below the ceiling. If not go to 10

      if (ucon(1).le.gammaradceiling*sqrt(-gn(1,1))) then
         go to 20
      else
         go to 10
      endif

c     This segment is for problem cases. We set gamma_radiation equal to
c     its ceiling value and solve for E and u^i without using R^00.

 10   ucon(1)=gammaradceiling*sqrt(-gn(1,1))
      radinvmod=1
#if(PRODUCTION==0)
      write (*,*) ' gamma_rad hit ceiling: g^tt, urad^t = ',
     &     gn(1,1),ucon(1)
      write (*,*) ' Rtcon: ',(Rtcon(j),j=1,4)
#endif
      aquad=0.d0
      bquad=0.d0
      cquad=1.d0+gv(1,1)*ucon(1)*ucon(1)
      do i=2,4
         bquad=bquad+1.5d0*gv(1,i)*Rtcon(i)
         cquad=cquad-0.5d0*gv(1,i)*gn(1,i)
      do j=2,4
         aquad=aquad+9.d0*gv(i,j)*Rtcon(i)*Rtcon(j)/
     &        (16.d0*ucon(1)*ucon(1))
         bquad=bquad-3.d0*gv(i,j)*(Rtcon(i)*gn(1,j)+
     &        Rtcon(j)*gn(1,i))/(16.d0*ucon(1)*ucon(1))
         cquad=cquad+gv(i,j)*gn(1,i)*gn(1,j)/
     &        (16.d0*ucon(1)*ucon(1))
      enddo
      enddo

      disc=bquad*bquad-4.d0*aquad*cquad

c     Check if disc > 0. If not, we have trouble again, and we need to
c     use yet another scheme!

      if (disc.lt.0.d0) go to 30

      E=(-bquad-sqrt(disc))/(2.d0*cquad)
#if(PRODUCTION==0)
      write (*,*) ' a, b, c, disc, E: ',aquad,bquad,cquad,disc,E
#endif
      do i=2,4
         ucon(i)=(3.d0*Rtcon(i)-E*gn(1,i))/(4.d0*E*ucon(1))
      enddo
#if(PRODUCTION==0)
c      write (*,*) ' ucon: ',(ucon(i),i=1,4)
#endif
      sum=0.d0
      do i=1,4
      do j=1,4
         sum=sum+gv(i,j)*ucon(i)*ucon(j)
      enddo
      enddo
#if(PRODUCTION==0)
c      write (*,*) ' check norm: ',sum
#endif
      go to 20

c     Third try! What should we do here?

 30   continue

c     Last-ditch effort. Set radiation velocity equal to gas velocity

#if(PRODUCTION==0)
      write (*,*) ' Error: No solution for radiation '
#endif
      do i=1,4
         ucon(i)=ugascon(i)
      enddo

      E=3.d0*Rtcon(1)/(4.d0*ucon(1)**2-gn(1,1))

c     Calculate urad_mu

 20   call contocov(ucon,ucov)

c     Calculate the full tensor R^mu_nu

      call calcRmunu(E,ucon,ucov,Rmunu)

      return
      end



c     New routine
      subroutine Rmunuinvert1(Rtcov,E,ucon,ucov,ugascon,Rmunu,radinvmod)

c     Given the row R^t_mu of the radiation tensor, solves for the              
c     radiation frame 4-velocity u^mu and energy density E and then             
c     calculates the full tensor R^mu_nu                                        

      implicit double precision (a-h,o-z)
      integer radinvmod
      dimension Rtcov(4),Rtcon(4),ucon(4),ucov(4),Rmunu(4,4),
     &     ugascon(4),gn(4,4),gv(4,4),etacov(4),etacon(4)
      common/accuracy/eps,epsbis,dvmin,tol,uminfac,dlogmax,
     &     gammaradceiling,itermax
      common/metric/gn,gv
      dimension Ttcovc(4),Rtcovc(4),BBc(4)
      common/conserved/Gam,Gam1,en,en1,rhou0c,sc,Ttcovc,Rtcovc,
     &     BBc,dt,
     &     FNUMARGSHARM,FNUMRESULTSHARM,WHICHVELRAMESHHARM,
     &     GAMMAMAXRAD,ERADLIMIT,toltry,tolallow,
     &     ARAD_CODE,OKAPPA_ES_CODE11,OKAPPA_FF_CODE11,OKAPPA_BF_CODE11

c     default is no ceiling hit on gamma
      radinvmod=0

c      arad=1.18316d17

c     Convert R^0_mu to R^0^mu                                                  

      call covtocon(Rtcov,Rtcon)
c      write (*,*) ' Rtcov: ',(Rtcov(i),i=1,4)                                  
c      write (*,*) ' Rtcon: ',(Rtcon(i),i=1,4)                                  

c     Calculate lapse alpha and ZAMO four velocity eta_mu and eta^mu

      alphasq=1.d0/(-gn(1,1))
      alpha=sqrt(alphasq)

      etacov(1)=-alpha
      etacov(2)=0.d0
      etacov(3)=0.d0
      etacov(4)=0.d0

      call covtocon(etacov,etacon)

      if (Rtcon(1).lt.0.d0) then

c     Check if R^tt is positive. If not, set energy density to a very
c     small value and set 4-velocity equal to ZAMO velocity.

         ucon(1)=1.d0/alpha
         ucon(2)=etacon(2)
         ucon(3)=etacon(3)
         ucon(4)=etacon(4)

         E = ERADLIMIT
c         E=arad*1.d-36

#if(PRODUCTION<=0)
         write (*,*) ' R^tt<0: ucon, E ',(ucon(i),i=1,4),E
#endif

      else

c     Compute |R|^2, xisqr=|R|^2*g^tt/(R^tt)^2, and xisqrceiling.  Check
c     which regime we are in and choose appropriate solution.

      Rsqr=0.d0
      do i=1,4
      do j=1,4
         Rsqr=Rsqr+gv(i,j)*Rtcon(i)*Rtcon(j)
      enddo
      enddo

      xisqr=Rsqr*gn(1,1)/Rtcon(1)**2
      gc=gammaradceiling
      xisqrc=(8.d0*gc**2+1.d0)/(16.d0*gc**4-8.d0*gc**2+1.d0)

#if(PRODUCTION<=0)
      write (*,*) ' R00, Rsqr, xisqr, xisqrceiling: ',
     &     Rtcon(1),Rsqr,xisqr,xisqrc
#endif

      if (xisqr.ge.1.d0) then
         radinvmod=1

c     If xisqr>=1, set gammarad=1, solve for E, and obtain remaining
c     urad^i.

         ucon(1)=1.d0/alpha
         ucon(2)=etacon(2)
         ucon(3)=etacon(3)
         ucon(4)=etacon(4)

         E=3.d0*Rtcon(1)/(4.d0*ucon(1)**2+gn(1,1))

#if(PRODUCTION<=0)
         write (*,*) ' xisqr>=1: ucon, E ',(ucon(i),i=1,4),E
#endif

      elseif (xisqr.le.xisqrc) then

c     If xisqr<=xisqrceiling, set gammarad=gammaradceiling and compute
c     the rest of the solution using our old method
         radinvmod=1

         ucon(1)=gammaradceiling*sqrt(-gn(1,1))

         aquad=0.d0
         bquad=0.d0
         cquad=1.d0+gv(1,1)*ucon(1)*ucon(1)
         do i=2,4
            bquad=bquad+1.5d0*gv(1,i)*Rtcon(i)
            cquad=cquad-0.5d0*gv(1,i)*gn(1,i)
         do j=2,4
            aquad=aquad+9.d0*gv(i,j)*Rtcon(i)*Rtcon(j)/
     &        (16.d0*ucon(1)*ucon(1))
            bquad=bquad-3.d0*gv(i,j)*(Rtcon(i)*gn(1,j)+
     &        Rtcon(j)*gn(1,i))/(16.d0*ucon(1)*ucon(1))
            cquad=cquad+gv(i,j)*gn(1,i)*gn(1,j)/
     &        (16.d0*ucon(1)*ucon(1))
         enddo
         enddo

         disc=bquad*bquad-4.d0*aquad*cquad
         E=(-bquad-sqrt(disc))/(2.d0*cquad)
#if(PRODUCTION<=0)
         write (*,*) ' a, b, c, disc, E: ',aquad,bquad,cquad,disc,E
#endif

         do i=2,4
            ucon(i)=(3.d0*Rtcon(i)-E*gn(1,i))/(4.d0*E*ucon(1))
         enddo

#if(PRODUCTION<=0)
         write (*,*) ' xisqr<=xisqrceil: ucon, E ',(ucon(i),i=1,4),E
#endif

      else

c     We have a physically valid situation with xisqr within the
c     acceptable range. Calculate the solution.

         gammaradsqr=((1.d0+xisqr)+sqrt(1.d0+3.d0*xisqr))/
     &        (4.d0*xisqr)

         ucon(1)=sqrt(gammaradsqr)/alpha
         E=3.d0*Rtcon(1)/(4.d0*ucon(1)**2+gn(1,1))

         do i=2,4
            ucon(i)=(3.d0*Rtcon(i)-E*gn(1,i))/(4.d0*E*ucon(1))
         enddo

#if(PRODUCTION<=0)
         write (*,*) ' physical solution: ucon, E ',(ucon(i),i=1,4),E
#endif

      endif

      endif

#if(PRODUCTION<=0)
c     Check if the 4-velocity is properly normalized                            
      sum=0.d0
      do i=1,4
      do j=1,4
         sum=sum+gv(i,j)*ucon(i)*ucon(j)
      enddo
      enddo

      write (*,*) ' check norm: ',sum
#endif

c     Calculate the full tensor R^mu_nu                                         
      call contocov(ucon,ucov)
      call calcRmunu(E,ucon,ucov,Rmunu)

      return
      end




      subroutine MHDinvert(prim,iflag,jflag,ientropy,itertot,guesstype)

c     Given initial primitives prim(4), solves for primitives that
c     satisfy the initial pre-radiation conserved quantities: rhou0c, s,
c     Ttcovc(4), Rtcovc(4)

      implicit double precision (a-h,o-z)
      real*8 itertot
      dimension prim(4),primsave(4),Ttcovc(4),Rtcovc(4),BBc(4),
     &     ucon(4)
      common/accuracy/eps,epsbis,dvmin,tol,uminfac,dlogmax,
     &     gammaradceiling,itermax
      common/conserved/Gam,Gam1,en,en1,rhou0c,sc,Ttcovc,Rtcovc,
     &     BBc,dt,
     &     FNUMARGSHARM,FNUMRESULTSHARM,WHICHVELRAMESHHARM,
     &     GAMMAMAXRAD,ERADLIMIT,toltry,tolallow,
     &     ARAD_CODE,OKAPPA_ES_CODE11,OKAPPA_FF_CODE11,OKAPPA_BF_CODE11
      external funcMHD1,funcMHD2,uMHD1,uMHD2

#if(PRODUCTION==0)
      write (*,*)
      write (*,*) ' initial primitives: ',(prim(j),j=1,4)
#endif

      do i=1,4
         primsave(i)=prim(i)
      enddo


      if(guesstype.le.1) then
c     First do one round of Newton-Raphson just on the velocities

      call Newton3(prim,iter,iflag,jflag,funcMHD1)
      itertot=itertot+iter
#if(PRODUCTION==0)
      write (12,*) ' MHD inversion, velocities only ',
     &        iter      
      write (*,*)
#endif
c     Next work only on u_g using the energy equation

      call usolveMHD(prim,iter,iflag,jflag,funcMHD1,uMHD1)
#if(PRODUCTION==0)
      write (*,*) ' usolveMHD done, u: ',prim(1)
      write (*,*)
      write (12,*) ' MHD inversion, u_g only (energy) '
#endif
      do i=1,4
         primsave(i)=prim(i)
      enddo

      endif

c     Carry out the full Newton-Raphson MHD inversion via the
c     energy equation

      err4=1.0d5
      call Newton4(prim,iter,iflag,jflag,funcMHD1,err4)
#if(PRODUCTION==0)
      write (12,*) ' MHD inversion, energy equation ',
     &     iter
#endif
      itertot=itertot+iter

c     Check if the iterations converged. If not, re-solve using the
c     entropy equation.

#if(PRODUCTION==0)
      if (iflag.eq.0) then
         write (*,*)
         write (*,*) ' Energy equation converged: iter ',
     &        iter
         write (*,*) ' primitives: ',(prim(j),j=1,4)
      else
         write (*,*)
         write (*,*) ' ERROR: no convergence with energy equation: ',
     &        iter
         write (*,*) ' switching to entropy equation '
      endif
#endif

      if (iflag.eq.0.and.prim(1).lt.0.d0) then
c     Check if the internal energy is okay. If not, re-solve using the
c     entropy equation.
#if(PRODUCTION==0)
         write (*,*)
         write (*,*) ' ERROR: u is negative!! '
         write (*,*) ' switching to entropy equation '
#endif
         iflag=2

      else

c     Check if entropy looks okay. If not, re-solve using the entropy
c     equation.

         pressure=Gam1*prim(1)
         call solveucon(prim,ucon)
         rhoi=rhou0c/ucon(1)

         entr=log(pressure**en/
     &        rhoi**en1)
         if ((entr-sc/rhou0c).lt.(-dlogmax)) then
#if(PRODUCTION==0)
            write (*,*)
            write (*,*) ' ERROR: s is too low!! '
            write (*,*) ' (entropy,sc)/(rho u^0): ',entr,sc/rhou0c
            write (*,*) ' switching to entropy equation '
#endif
            iflag=3
         endif

      endif

      if (iflag.ge.1) then

c     Restore saved primitives to post-3+1 solution before proceeding to
c     solving the entropy equation
         
c     Reset jflag
         jflag=0

            

         do i=1,4
            prim(i)=primsave(i)
         enddo
         ientropy=1

c     Next work only on u_g using the entropy equation

c      call usolveMHD(prim,iter,iflag,jflag,funcMHD2,uMHD2)
#if(PRODUCTION==0)
c      write (*,*)
c      write (*,*) ' usolveMHD done, u: ',prim(1)
c      write (*,*)
c      write (12,*) ' Radiation inversion, u_g only (entropy) '
#endif

c     Carry out full Newton-Raphson inversion with the entropy equation

      err4=1.0d6
      call Newton4(prim,iter,iflag,jflag,funcMHD2,err4)
      itertot=itertot+iter

      if (iflag.eq.9) then
#if(PRODUCTION==0)
         write (*,*) ' ERROR: u is negative! '
         write (*,*) ' reverting to previous solution '
#endif
         do j=1,4
            prim(j)=primsave(j)
         enddo
      else
#if(PRODUCTION==0)
         write (12,*) ' MHD inversion, entropy equation ',
     &        iter
         write (*,*)
#endif
      endif

      if (iflag.eq.0) then
#if(PRODUCTION==0)
         write (*,*) ' entropy equation converged: iter ',
     &        iter
         write (*,*) ' primitives: ',(prim(j),j=1,4)
#endif
      else
#if(PRODUCTION==0)
         write (*,*) ' ERROR: no convergence with entropy equation ',
     &        iter
         write (*,*) ' reverting to saved solution '
#endif
         do j=1,4
            prim(j)=primsave(j)
         enddo
         iflag=9
      endif

      endif

      return
      end



      subroutine radsource(prim,iter,iflag,jflag,ientropy,
     &     itertot,guesstype,primeng,priment,
     &     resultseng,resultsent)

c     Given initial primitives prim(4), solves for primitives that
c     satisfy the post-radiation equations

      implicit double precision (a-h,o-z)
      double precision resultseng(NUMRESULTS)
      double precision resultsent(NUMRESULTS)
      double precision primeng(4),priment(4)
      real*8 itertot
      dimension prim(4),primsave(4)
      common/accuracy/eps,epsbis,dvmin,tol,uminfac,dlogmax,
     &     gammaradceiling,itermax
      dimension Ttcovc(4),Rtcovc(4),BBc(4)
      common/conserved/Gam,Gam1,en,en1,rhou0c,sc,Ttcovc,Rtcovc,
     &     BBc,dt,
     &     FNUMARGSHARM,FNUMRESULTSHARM,WHICHVELRAMESHHARM,
     &     GAMMAMAXRAD,ERADLIMIT,toltry,tolallow,
     &     ARAD_CODE,OKAPPA_ES_CODE11,OKAPPA_FF_CODE11,OKAPPA_BF_CODE11
      external funcrad1,funcrad2,uerr1,uerr2
      dimension ugascon(4)

      itereng=0
      erreng=1.d3
      iterent=0
      errent=1.d3

#if(PRODUCTION==0)
      write (*,*)
      write (*,*) ' initial primitives: ',(prim(j),j=1,4)
      write (*,*)
#endif

      do i=1,4
         primsave(i)=prim(i)
      enddo

      if(guesstype.le.1) then

c     First do one round of Newton-Raphson just on the velocities

      call Newton3(prim,iter,iflag,jflag,funcrad1)
      itertot=itertot+iter
#if(PRODUCTION==0)
      write (12,*) ' Radiation inversion, velocities only ',
     &        iter      
      write (*,*)
#endif

      do i=1,4
         primsave(i)=prim(i)
      enddo

      if (ientropy.eq.0) then

c     ientropy=0, so we will first try the energy equation. Initially
c     work only on u_g using the energy equation. If ientropy=1, we skip
c     all this and go directly to working with the entropy equation.

      call usolverad(prim,iter,iflag,jflag,funcrad1,uerr1)
#if(PRODUCTION==0)
      write (*,*) ' usolverad done, u: ',prim(1)
      write (*,*)
      write (12,*) ' Radiation inversion, u_g only (energy) '
#endif
      if (prim(1).gt.0.d0) primsave(1)=prim(1)

      endif

c   endif guesstype.le.1
      endif

c     Carry out full Newton-Raphson on all four primitives using
c     Newton-Raphson

      envconv=0
      err4=1.0d8
      call Newton4(prim,iter,iflag,jflag,funcrad1,err4)
      itertot=itertot+iter
      itereng=iter
      erreng=err4
#if(PRODUCTION==0)
      write (*,*)
      write (12,*) ' Radiation inversion, energy equation ',
     &        iter,err4,isnan(erreng)
#endif

c     default is failed
      resultseng(12)=1.0d0
      if (iflag.eq.0.and.erreng.lt.FAILLIMIT) then
#if(PRODUCTION==0)
         write (*,*) ' Energy equation converged: iter ',
     &        iter
         write (*,*) ' primitives: ',(prim(j),j=1,4)
#endif
         call solveucon(prim,ugascon)

         if (erreng.lt.FAILLIMIT
     &        .and.prim(1).gt.0.0
     &        .and.jflag.eq.0
     &        .and.prim(1).eq.prim(1)
     &        .and.prim(2).eq.prim(2)
     &        .and.prim(3).eq.prim(3)
     &        .and.prim(4).eq.prim(4)
     &        .and.ugascon(1).eq.ugascon(1)
     &        ) then
#if(PRODUCTION<=2)
            write (14,"(1A)",advance="no") '  GOOD   '
#endif
#if(PRODUCTION<=1)
            write (13,"(1A)",advance="no") '  GOOD   '
#endif
         resultseng(12)=0.0d0
         engconv=1
         else if (erreng.lt.FAILLIMIT
     &        .and.prim(1).eq.prim(1)
     &        .and.prim(2).eq.prim(2)
     &        .and.prim(3).eq.prim(3)
     &        .and.prim(4).eq.prim(4)
     &        .and.ugascon(1).eq.ugascon(1)
     &        ) then
#if(PRODUCTION<=2)
            write (14,"(1A)",advance="no") '   BADNEG'
#endif
#if(PRODUCTION<=1)
            write (13,"(1A)",advance="no") '   BADNEG'
#endif
c     only bad because negative something
         resultseng(12)=2.0d0
         engconv=1
         else
#if(PRODUCTION<=2)
            write (14,"(1A)",advance="no") '   BAD   '
#endif
#if(PRODUCTION<=1)
            write (13,"(1A)",advance="no") '   BAD   '
#endif
c            treat as bad because error not small or nan'ed out
         resultseng(12)=1.0d0
         engconv=0
         endif
c     commenting below to see entropy even if energy converged.
c         return
      else
         resultseng(12)=1.0d0
         engconv=0

c     If energy equation does not converge, calculate using the entropy
c     equation.

#if(PRODUCTION==0)
       write (*,*) ' ERROR: no convergence with energy equation: ',
     &     iter
      write (*,*) ' Proceed to entropy equation '
      write (*,*)
      write (12,*) ' ERROR: no convergence with energy equation: '
#endif
#if(PRODUCTION<=2)
      write (14,"(1A)",advance="no") '   BAD   '
#endif
#if(PRODUCTION<=1)
      write (13,"(1A)",advance="no") '   BAD   '
#endif

      endif


      resultseng(13)=erreng
      resultseng(14)=DBLE(itereng)
#if(PRODUCTION<=2)
      if(1.eq.1) then
          write (14,"(9X,1I3,3X,1E21.15)",advance="no")
     &        itereng,erreng
      else
         write (14,*) itereng,erreng,iterent,errent
      endif
#endif

      do i=1,4
         primeng(i)=prim(i)
      enddo

c     Restore saved primitives and do Newton-Raphson with the entropy
c     equation

      do i=1,4
         prim(i)=primsave(i)
      enddo


c     Next work only on u_g using the entropy equation. Since we are
c     caurrently using the entropy equation even for the energy equation
c     step, we do not need to repeat the 1D search. If and when we
c     modify uerr1 to do the proper energy equation, we need this second
c     round of 1D search using uerr2.

c      call usolverad(prim,iter,iflag,jflag,funcrad2,uerr2)
#if(PRODUCTION==0)
c      write (*,*)
c      write (*,*) ' usolverad done, u: ',prim(1)
c      write (*,*)
c      write (12,*) ' Radiation inversion, u_g only (entropy) '
#endif

c      default is ent failed
      resultsent(12)=1.0d0
      err4=1.0d8
      call Newton4(prim,iter,iflag,jflag,funcrad2,err4)
      itertot=itertot+iter
      iterent=iter
      errent=err4
      if (iflag.eq.0.and.errent.lt.FAILLIMIT) then
#if(PRODUCTION==0)
      write (*,*)
      write (12,*) ' Radiation inversion, entropy equation ',
     &        iter
#endif
      call solveucon(prim,ugascon)
      if ((erreng.lt.FAILLIMIT.or.errent.lt.FAILLIMIT)
     &     .and.prim(1).gt.0.0
     &     .and.jflag.eq.0
     &     .and.prim(1).eq.prim(1)
     &     .and.prim(2).eq.prim(2)
     &     .and.prim(3).eq.prim(3)
     &     .and.prim(4).eq.prim(4)
     &     .and.ugascon(1).eq.ugascon(1)
     &     ) then
#if(PRODUCTION<=2)
         write (14,"(1A)",advance="no") '  GOOD   '
#endif
#if(PRODUCTION<=1)
         write (13,"(1A)",advance="no") '  GOOD   '
#endif
      resultsent(12)=0.0d0
      else if ((erreng.lt.FAILLIMIT.or.errent.lt.FAILLIMIT)
     &     .and.prim(1).eq.prim(1)
     &     .and.prim(2).eq.prim(2)
     &     .and.prim(3).eq.prim(3)
     &     .and.prim(4).eq.prim(4)
     &     .and.ugascon(1).eq.ugascon(1)
     &     ) then
#if(PRODUCTION<=2)
         write (14,"(1A)",advance="no") '   BADNEG'
#endif
#if(PRODUCTION<=1)
         write (13,"(1A)",advance="no") '   BADNEG'
#endif
      resultsent(12)=2.0d0
      else
#if(PRODUCTION<=2)
         write (14,"(1A)",advance="no") '   BAD   '
#endif
#if(PRODUCTION<=1)
         write (13,"(1A)",advance="no") '   BAD   '
#endif
      resultsent(12)=1.0d0
      endif
#if(PRODUCTION==0)
      if (iflag.eq.0.and.errent.lt.FAILLIMIT) then
         write (*,*) ' Entropy equation converged: iter ',
     &        iter
         write (*,*) ' primitives: ',(prim(j),j=1,4)
c         return
      endif
#endif
      else
      resultsent(12)=1.0d0

#if(PRODUCTION==0)
      write (*,*) ' ERROR: no convergence with entropy equation: ',
     &     iter
      write (12,*) ' ERROR: no convergence with entropy equation: '
#endif
#if(PRODUCTION<=2)
      write (14,"(1A)",advance="no") '   BAD   '
#endif
#if(PRODUCTION<=1)
      write (13,"(1A)",advance="no") '   BAD   '
#endif

      endif
      
      resultsent(13)=errent
      resultsent(14)=DBLE(iterent)
#if(PRODUCTION<=2)
      if(1.eq.1) then
         write (14,*) iterent,errent
      else
         write (14,*) itereng,erreng,iterent,errent
      endif
#endif


      do i=1,4
         priment(i)=prim(i)
      enddo

c     We should never reach this point. Unclear what to do in this case!
c     We could keep the solution or restore the saved primitives.


c      do i=1,4
c         prim(i)=primsave(i)
c      enddo

      return
      end



      subroutine Newton4 (prim0,iter,iflag,jflag,func,err4)

c     Calculates iteratively by the Newton-Raphson technique the
c     solution to a set of four non-linear equations. The initial
c     primitives are in the array prim0, and the final solution is
c     returned in the same array. func is the function that computes the
c     four error terms for a given set of primitives.

      implicit double precision (a-h,o-z)
      dimension prim0(4),error0(4),errornorm(4),prim(4),error(4),
     &     AJac(4,4),indx(4),Ttcov(4),Rtcov(4),BB(4),
     &     primold(4)
      common/accuracy/eps,epsbis,dvmin,tol,uminfac,dlogmax,
     &     gammaradceiling,itermax
      common/conserved/Gam,Gam1,en,en1,rhou0,s,Ttcov,Rtcov,BB,dt,
     &     FNUMARGSHARM,FNUMRESULTSHARM,WHICHVELRAMESHHARM,
     &     GAMMAMAXRAD,ERADLIMIT,toltry,tolallow,
     &     ARAD_CODE,OKAPPA_ES_CODE11,OKAPPA_FF_CODE11,OKAPPA_BF_CODE11
      external func
      double precision mymax,myabs,mydiv
      integer myisnan

c     niter is the maximum number of Newton-Raphson iterations
c     iflag=0 means that a good solution was found

      niter=100
c      niter=20
      iflag=0
      jflag=0

c     Do Newton-Raphson until err is smaller than tolerance tol

      iter=0
      do i=1,niter

c     Make sure u == prim0(1) is in safe territory

      prim0(1)=mymax(prim0(1),uminfac*rhou0)

      call func(prim0,error0,errornorm,err4,err3,iflag,jflag)
#if(PRODUCTION==0)
      write (*,*) ' Newton-Raphson-4 '
      write (*,*) ' i, primitives, error, errornorm, err4: ',i,
     &     (prim0(j),j=1,4),(error0(j),j=1,4),
     &     (errornorm(j),j=1,4),err4
#endif
c     iflag=9 means a serious error in the value of u. This can be
c     tolerated for a few steps (up to iter=itermax). After that, return
c     with iflag=9.

      if (iflag.eq.9.and.iter.ge.itermax) then
#if(PRODUCTION==0)
         write (*,*) 'Hit iflag=',iflag,'iter=',iter
#endif
         return
      endif

c     If all the four equations give fractional errors less than tol,
c     return

      if (err4.lt.tol) then
#if(PRODUCTION==0)
         write (*,*) 'Hit err4=',err4
#endif
         return
      endif

      if (i.gt.1) then
         sum=0.d0
         do j=1,4
            sum=sum+myabs(prim0(j)-primold(j))
         enddo
         if (sum.eq.0.d0) then
#if(PRODUCTION==0)
            write (12,*) ' sum = 0! '
#endif
            return
         endif
      endif

      do j=1,4
         primold(j)=prim0(j)
      enddo

c     Calculate the Jacobian numerically by shifting the primitives one
c     by one by a fraction eps of their current values and calculating
c     the errors.

      do j=1,4
         
      do k=1,4
         prim(k)=prim0(k)
      enddo

         if (j.eq.1.or.myabs(prim0(j)).gt.dvmin) then
            dprim=prim0(j)*eps
         else
            dprim=dvmin*eps
         endif
         prim(j)=prim0(j)+dprim
c         prim(j)=prim0(j)-myabs(dprim)

         call func(prim,error,errornorm,err4,err3,iflag,jflag)
#if(PRODUCTION==0)
c         write (*,*) ' error: ',j,(error(k),k=1,4)
#endif
      do k=1,4
         AJac(k,j)=(error(k)-error0(k))/dprim
      enddo

      enddo
#if(PRODUCTION==0)
c      write (*,*) ' Jacobian: ',((AJac(k,j),j=1,4),k=1,4)
#endif

c     Invert the Jacobian using subroutines from Numerical Recipes and
c     compute the shifts to the primitives

      call ludcmp(AJac,4,4,indx,d,retval)
      if(retval.ne.0) then
         err4=128.0
         return
      endif
      call lubksb(AJac,4,4,indx,error0)

#if(PRODUCTION==0)
c      do j=1,4
c         write (*,*) ' j, prim0(j), dprim(j): ',j,prim0(j),
c     &        error0(j)
c      enddo
#endif

c     Apply the Newton-Raphson shifts

      do j=1,4
         prim0(j)=prim0(j)-error0(j)
      enddo

      iter=iter+1

      enddo

      iflag=1

      return
      end



      subroutine Newton3 (prim0,iter,iflag,jflag,func)

c     Calculates iteratively by the Newton-Raphson technique the
c     solution to a set of three non-linear equations. The initial
c     primitives are in the 4-array prim0, of which the first element is
c     not varied and only the other three are solved for. The final
c     solution is returned in the same 4-array with the first element
c     unchanged. func is the function that computes the four error terms
c     for a given set of primitives. Only the last three errors are
c     used.

      implicit double precision (a-h,o-z)
      dimension prim0(4),prim30(3),error30(3),prim3(3),error3(3),
     &     AJac(3,3),indx(3),Ttcov(4),Rtcov(4),BB(4),errornorm(3)
      common/accuracy/eps,epsbis,dvmin,tol,uminfac,dlogmax,
     &     gammaradceiling,itermax
      common/conserved/Gam,Gam1,en,en1,rhou0,s,Ttcov,Rtcov,BB,dt,
     &     FNUMARGSHARM,FNUMRESULTSHARM,WHICHVELRAMESHHARM,
     &     GAMMAMAXRAD,ERADLIMIT,toltry,tolallow,
     &     ARAD_CODE,OKAPPA_ES_CODE11,OKAPPA_FF_CODE11,OKAPPA_BF_CODE11
      external func
      double precision mymax,myabs,mydiv
      integer myisnan

c     niter is the maximum number of Newton-Raphson
c     iterations. Currently we do only one iteration of Newton3.
c     iflag=0 measn that a good solution was found

      niter=1
c      niter=5
      iflag=0
      jflag=0

c     Make sure u == primsave(1) is in safe territory

      primsave=mymax(prim0(1),uminfac*rhou0)
      do i=1,3
         prim30(i)=prim0(i+1)
      enddo

c     Do Newton-Raphson until err is less than tolerance tol

      iter=0
      do i=1,100

      call func3(primsave,prim30,error30,errornorm
     &        ,err4,err3,iflag,jflag,func)
#if(PRODUCTION==0)
      write (*,*) ' Newton-Raphson-3 '
      write (*,*) ' iter, primitives, error, errornorm, err3: ',iter,
     &     (prim30(j),j=1,3),(error30(j),j=1,3),
     &     (errornorm(j),j=1,3),err3
#endif

      if (iter.ge.niter) then
         prim0(1)=primsave
         do j=1,3
            prim0(j+1)=prim30(j)
         enddo
         return
      endif

c     If all the three equations give fractional errors less than tol,
c     return

      if (err3.lt.tol) then
         prim0(1)=primsave
         do j=1,3
            prim0(j+1)=prim30(j)
         enddo
         return
      endif

c     Calculate the Jacobian numerically by shifting the primitives one
c     by one by a fraction eps of their current values and calculating
c     the errors.

      do j=1,3
         
      do k=1,3
         prim3(k)=prim30(k)
      enddo

         if (myabs(prim30(j)).gt.dvmin) then
            dprim=prim30(j)*eps
         else
            dprim=dvmin*eps
         endif
         prim3(j)=prim30(j)+dprim
c         prim3(j)=prim30(j)-myabs(dprim)

         call func3(primsave,prim3,error3,errornorm,
     &        err4,err3,iflag,jflag,func)
#if(PRODUCTION==0)
c         write (*,*) ' error: ',j,(error3(k),k=1,3),err3
#endif
      do k=1,3
         AJac(k,j)=(error3(k)-error30(k))/dprim
      enddo

      enddo
#if(PRODUCTION==0)
c      write (*,*) ' Jacobian: ',((AJac(k,j),j=1,3),k=1,3)
#endif

c     Invert the Jacobian using subroutines from Numerical Recipes and
c     compute the shifts to the primitives

      call ludcmp(AJac,3,3,indx,d,retval)
      if(retval.ne.0) then
         err3=1.0d4
         return
      endif
      call lubksb(AJac,3,3,indx,error30)

#if(PRODUCTION==0)
c      do j=1,3
c         write (*,*) ' j, prim30(j), dprim(j): ',j,prim0(j),
c     &        error30(j)
c      enddo
#endif
c     Apply the Newton-Raphson shifts

      do j=1,3
         prim30(j)=prim30(j)-error30(j)
      enddo

      iter=iter+1

      enddo

      iflag=1
#if(PRODUCTION==0)
c      write (*,*) ' too many iterations! '
#endif

      prim0(1)=primsave
      do i=1,3
         prim0(i+1)=prim30(i)
      enddo

      return
      end



      subroutine usolveMHD(prim,iter,iflag,jflag,funcMHD,uMHD)

c     Solves a 1D equation for u, using either the energy or entropy
c     equation without radiation source term

      implicit double precision (a-h,o-z)
      dimension prim(4),error(4),errornorm(4),Ttcov(4),Rtcov(4),BB(4),
     &     ucon(4),ucov(4),bcon(4),bcov(4),Tmunu(4,4),
     &     urcon(4),urcov(4),Rmunu(4,4),Gcon(4),Gcov(4)
      common/accuracy/eps,epsbis,dvmin,tol,uminfac,dlogmax,
     &     gammaradceiling,itermax
      common/conserved/Gam,Gam1,en,en1,rhou0,s,Ttcov,Rtcov,BB,dt,
     &     FNUMARGSHARM,FNUMRESULTSHARM,WHICHVELRAMESHHARM,
     &     GAMMAMAXRAD,ERADLIMIT,toltry,tolallow,
     &     ARAD_CODE,OKAPPA_ES_CODE11,OKAPPA_FF_CODE11,OKAPPA_BF_CODE11
      common/funcMHDd/ffkap,eskap,arad,ucon,ucov,rho,bcon,bcov,bsq,
     &     Tmunu,E,Ehat,urcon,urcov,Rmunu,Tgas,Trad,B4pi,gamma,dtau,
     &     Gcon,Gcov
      external funcMHD,uMHD
      double precision mymax,myabs,mydiv
      integer myisnan

c     Call funcMHD and obtain basic parameters needed for solving the 1D
c     energy equation: rho, Ehat

      u0=prim(1)
      call funcMHD(prim,error,errornorm,err4,err3,iflag,jflag)
      rho0=rho
      Ehat0=Ehat
#if(PRODUCTION==0)
c      write (*,*) ' usolveMHD: ',(prim(i),i=1,4),rho0,Ehat0
#endif

c     Calculate initial error

      call uMHD(u0,u0,rho0,Ehat0,ffkap,arad,dtau,err0)
#if(PRODUCTION==0)
      write (*,*) ' Bracketing the solution '
      write (*,*) ' u, err: ',u0,err0
#endif

c     Decide whether u is too low or too high and search accordingly to
c     bracket the solution for u

      if (err0.lt.0.d0) then

         ur=u0
         er=err0

         do i=1,50
           ul=ur
           el=er
           ur=2.d0*ur
           call uMHD(ur,u0,rho0,Ehat0,ffkap,arad,dtau,er)
#if(PRODUCTION==0)
           write (*,*) ' u, err: ',ur,er
#endif
           if (er.ge.0.d0) go to 10
         enddo
#if(PRODUCTION==0)
         write (*,*) ' No bracket! '
#endif
         return

      else

         ul=u0
         el=err0

         do i=1,50
           ur=ul
           er=el
           ul=0.5d0*ul
           call uMHD(ul,u0,rho0,Ehat0,ffkap,arad,dtau,el)
#if(PRODUCTION==0)
           write (*,*) ' u, err: ',ul,el
#endif
           if (er.lt.0.d0) go to 10
         enddo
#if(PRODUCTION==0)
         write (*,*) ' No bracket! '
#endif
         return

      endif

c     Bracketing is done. Now solve for u by bisection

 10   continue
#if(PRODUCTION==0)
      write (*,*)
      write (*,*) ' Solving by bisection '
#endif

      do i=1,50
         umid=0.5d0*(ul+ur)
         call uMHD(umid,u0,rho0,Ehat0,ffkap,arad,dtau,emid)
#if(PRODUCTION==0)
         write (*,*) ' u, err: ',umid,emid
#endif
         if (emid.ge.0.d0) then
            ur=umid
            er=emid
         else
            ul=umid
            el=emid
         endif
         if (myabs(ur-ul).lt.epsbis*ul) then
            prim(1)=0.5d0*(ul+ur)
            return
         endif

      enddo

      prim(1)=0.5d0*(ul+ur)
#if(PRODUCTION==0)
      write (*,*) ' bisection did not converge! '
#endif

      return
      end



      subroutine uMHD1(u,u0,rho,Ehat,ffkap,arad,dtau,err)

c     Error function for 1D search in u for the energy equation without
c     radiation source term. Currently, this has the entropy equation.

      implicit double precision (a-h,o-z)
      dimension Ttcov(4),Rtcov(4),BB(4)
      common/conserved/Gam,Gam1,en,en1,rhou0,s,Ttcov,Rtcov,BB,dt,
     &     FNUMARGSHARM,FNUMRESULTSHARM,WHICHVELRAMESHHARM,
     &     GAMMAMAXRAD,ERADLIMIT,toltry,tolallow,
     &     ARAD_CODE,OKAPPA_ES_CODE11,OKAPPA_FF_CODE11,OKAPPA_BF_CODE11

c     Compute entropy and compute deviation from conserved entropy
      
      Tgas=Gam1*u/rho
c      B4pi=arad*Tgas**4
c      ff=ffkap*rho*rho/Tgas**(3.5d0)
c      Gdt=ff*(Ehat-B4pi)*dt
      entropy=rhou0*log((Gam1*u)**en/rho**en1)
      err=entropy-s
c      err=entropy-s-Gdtau
#if(PRODUCTION==0)
c      write (*,*) ' rhou0, rho, u, entropy, s, err: ',
c     &     rhou0,rho,u,entropy,s,err
#endif

      return
      end



      subroutine uMHD2(u,u0,rho,Ehat,ffkap,arad,dtau,err)

c     Error function for 1D search in u for the entropy equation without
c     radiation source term.

      implicit double precision (a-h,o-z)
      dimension Ttcov(4),Rtcov(4),BB(4)
      common/conserved/Gam,Gam1,en,en1,rhou0,s,Ttcov,Rtcov,BB,dt,
     &     FNUMARGSHARM,FNUMRESULTSHARM,WHICHVELRAMESHHARM,
     &     GAMMAMAXRAD,ERADLIMIT,toltry,tolallow,
     &     ARAD_CODE,OKAPPA_ES_CODE11,OKAPPA_FF_CODE11,OKAPPA_BF_CODE11
      
c     Compute entropy and compute deviation from conserved entropy
      
      Tgas=Gam1*u/rho
c      B4pi=arad*Tgas**4
c      ff=ffkap*rho*rho/Tgas**(3.5d0)
c      Gdt=ff*(Ehat-B4pi)*dt
      entropy=rhou0*log((Gam1*u)**en/rho**en1)
      err=entropy-s
c      err=entropy-s-Gdtau

      return
      end



      subroutine usolverad(prim,iter,iflag,jflag,funcrad,uerr)

c     Solves a 1D equation for u, using either the energy or entropy
c     equation including the radiation source term

      implicit double precision (a-h,o-z)
      dimension prim(4),error(4),Ttcov(4),Rtcov(4),BB(4),
     &     ucon(4),ucov(4),bcon(4),bcov(4),Tmunu(4,4),
     &     urcon(4),urcov(4),Rmunu(4,4),Gcon(4),Gcov(4),
     &     errornorm(4)
      common/accuracy/eps,epsbis,dvmin,tol,uminfac,dlogmax,
     &     gammaradceiling,itermax
      common/conserved/Gam,Gam1,en,en1,rhou0,s,Ttcov,Rtcov,BB,dt,
     &     FNUMARGSHARM,FNUMRESULTSHARM,WHICHVELRAMESHHARM,
     &     GAMMAMAXRAD,ERADLIMIT,toltry,tolallow,
     &     ARAD_CODE,OKAPPA_ES_CODE11,OKAPPA_FF_CODE11,OKAPPA_BF_CODE11
      common/funcradd/ffkap,eskap,arad,ucon,ucov,rho,bcon,bcov,bsq,
     &     Tmunu,E,Ehat,urcon,urcov,Rmunu,Tgas,Trad,B4pi,gamma,dtau,
     &     Gcon,Gcov,u0,iuerr
      external funcrad,uerr,uuerr

c     Call funcrad and obtain basic parameters needed for solving the 1D
c     energy equation: rho, Ehat

      u0=prim(1)
      call funcrad(prim,error,errornorm,err4,err3,iflag,jflag)
      rho0=rho
      Ehat0=Ehat

c     Calculate initial error

#if(PRODUCTION==0)
      write (*,*) ' Bracketing solution for u_g '
#endif
      iuerr=0
      call uuerr(u0,err0,derr)
c      call uerr(u0,u0,rho0,Ehat0,ffkap,arad,dtau,err0)
#if(PRODUCTION==0)
c      write (*,*) ' u, err: ',u0,err0
#endif

c     Decide whether u is too low or too high and search accordingly to
c     bracket the solution for u

      if (err0.lt.0.d0) then

         ur=u0
         er=err0

         do i=1,50
           ul=ur
           el=er
           ur=3.d0*ur
           call uuerr(ur,er,derr)
#if(PRODUCTION==0)
c           write (*,*) ' u, err: ',ur,er
#endif
           if (er.ge.0.d0) go to 10
         enddo
#if(PRODUCTION==0)
         write (*,*) ' No bracket! '
#endif
         return

      else

         ul=u0
         el=err0

         do i=1,50
           ur=ul
           er=el
           ul=0.3d0*ul
           call uuerr(ul,el,derr)
#if(PRODUCTION==0)
c           write (*,*) ' u, err: ',ul,el
#endif
           if (el.lt.0.d0) go to 10
         enddo
#if(PRODUCTION==0)
         write (*,*) ' No bracket! '
#endif
         return

      endif

c     Bracketing is done. Now solve for u using rtsafe

 10   continue

      iuerr=1
#if(PRODUCTION==0)
      write (*,*)
      write (*,*) ' Solving for u_g using rtsafe '
#endif

      uacc=eps*min(ul,ur)
      ugas=rtsafe(uuerr,ul,ur,el,er,uacc)
      prim(1)=ugas

#if(PRODUCTION==0)
c      write (*,*)
c      write (*,*) ' Solving for u_g by bisection '
#endif

c      do i=1,50
c         umid=0.5d0*(ul+ur)
c         call uerr(umid,u0,rho0,Ehat0,ffkap,arad,dtau,emid)
#if(PRODUCTION==0)
c         write (*,*) ' u, err: ',umid,emid
#endif
c         if (emid.ge.0.d0) then
c            ur=umid
c            er=emid
c         else
c            ul=umid
c            el=emid
c         endif
c         if (myabs(ur-ul).lt.epsbis*ul) then
c            prim(1)=0.5d0*(ul+ur)
c            return
c         endif

c      enddo

c      prim(1)=0.5d0*(ul+ur)
#if(PRODUCTION==0)
c      write (*,*) ' bisection did not converge! '
#endif

      return
      end



      subroutine uerr1(u,u0,rho,Ehat0,ffkap,arad,dtau,err)

c     Error function for 1D search in u for the energy equation
c     including the radiation source term. Currently, this has the
c     entropy equation.

      implicit double precision (a-h,o-z)
      dimension Ttcov(4),Rtcov(4),BB(4)
      common/conserved/Gam,Gam1,en,en1,rhou0,s,Ttcov,Rtcov,BB,dt,
     &     FNUMARGSHARM,FNUMRESULTSHARM,WHICHVELRAMESHHARM,
     &     GAMMAMAXRAD,ERADLIMIT,toltry,tolallow,
     &     ARAD_CODE,OKAPPA_ES_CODE11,OKAPPA_FF_CODE11,OKAPPA_BF_CODE11

c     Compute entropy and compute deviation from conserved entropy
      
      Tgas=Gam1*u/rho
      B4pi=arad*Tgas**4
      ff=ffkap*rho*rho/Tgas**(3.5d0)

      Ehat=u0+Ehat0-u
      Gdtau=ff*(Ehat-B4pi)*dtau

      entropy=rhou0*log((Gam1*u)**en/rho**en1)
      err=entropy-s-Gdtau

      return
      end



      subroutine uerr2(u,u0,rho,Ehat0,ffkap,arad,dtau,err)

c     Error function for 1D search in u for the entropy equation
c     including the radiation source term

      implicit double precision (a-h,o-z)
      dimension Ttcov(4),Rtcov(4),BB(4)
      common/conserved/Gam,Gam1,en,en1,rhou0,s,Ttcov,Rtcov,BB,dt,
     &     FNUMARGSHARM,FNUMRESULTSHARM,WHICHVELRAMESHHARM,
     &     GAMMAMAXRAD,ERADLIMIT,toltry,tolallow,
     &     ARAD_CODE,OKAPPA_ES_CODE11,OKAPPA_FF_CODE11,OKAPPA_BF_CODE11
      
c     Compute entropy and compute deviation from conserved entropy

      Tgas=Gam1*u/rho
      B4pi=arad*Tgas**4
      ff=ffkap*rho*rho/Tgas**(3.5d0)
      
      Ehat=u0+Ehat0-u
      Gdtau=ff*(Ehat-B4pi)*dtau

      entropy=rhou0*log((Gam1*u)**en/rho**en1)
      err=entropy-s-Gdtau

      return
      end



      subroutine uuerr(u,err,derr)

      implicit double precision (a-h,o-z)
      dimension Ttcov(4),Rtcov(4),BB(4),ucon(4),ucov(4),bcon(4),
     &     bcov(4),Tmunu(4,4),urcon(4),urcov(4),Rmunu(4,4),Gcon(4),
     &     Gcov(4)
      common/conserved/Gam,Gam1,en,en1,rhou0,s,Ttcov,Rtcov,BB,dt,
     &     FNUMARGSHARM,FNUMRESULTSHARM,WHICHVELRAMESHHARM,
     &     GAMMAMAXRAD,ERADLIMIT,toltry,tolallow,
     &     ARAD_CODE,OKAPPA_ES_CODE11,OKAPPA_FF_CODE11,OKAPPA_BF_CODE11
      common/funcradd/ffkap,eskap,arad,ucon,ucov,rho,bcon,bcov,bsq,
     &     Tmunu,E,Ehat0,urcon,urcov,Rmunu,Tgas,Trad,B4pi,gamma,dtau,
     &     Gcon,Gcov,u0,iuerr

c     The fluid frame error is:
c        rho*ln(p^n/rho^(n+1)) - s - (ffkap*rho^2/T^3.5)*dtau*(Ehat-B4pi)
c     We compute coefficients corresponding to these terms and then
c     compute the error and its derivative wrt u

      c1=Gam1/rho
      Tgas=c1*u

      c2=arad*c1**4
      B4pi=c2*u**4

      c3=ffkap*rho*rho/c1**(3.5d0)
      ff=c3/u**(3.5d0)
#if(PRODUCTION==0)
c      write (*,*) ' Tgas, B4pi, ff: ',Tgas,B4pi,ff
#endif

      c4=u0+Ehat0
      Ehat=c4-u

      Gdtau=(c3*dtau)*u**(-3.5d0)*(c4-u-c2*u**4)

      entropy=rho*log((Gam1*u)**en/rho**en1)

#if(PRODUCTION==0)
c      write (*,*) ' Ehat, Gdtau, entropy: ',Ehat,Gdtau,entropy
#endif

      err=entropy-(s/ucon(1))-Gdtau
#if(PRODUCTION==0)
c      write (*,*) ' s, ucon(1), err: ',s,ucon(1),err
#endif

      if (iuerr.eq.1) then

         derr=(rho*en/u)+3.5d0*c3*dtau*c4*u**(-4.5d0)-
     &        2.5d0*c3*dtau*u**(-3.5d0)+
     &        0.5d0*c3*dtau*c2*u**(-0.5d0)
#if(PRODUCTION==0)
         write (*,10) u,err,derr
#endif
 10      format (' u, err, derr: ',3e20.12)

      else

#if(PRODUCTION==0)
         write (*,*) ' u, err: ',u,err
#endif

      endif

      return
      end



      subroutine funcMHD1(prim,error,errornorm,err4,err3,iflag,jflag)

c     This subroutine calculates errors for the MHD inversion problem
c     without radiation source term using the energy equation

      implicit double precision (a-h,o-z)
      dimension prim(4),error(4),errornorm(4),Ttcov(4),Rtcov(4),BB(4),
     &     ucon(4),ucov(4),bcon(4),bcov(4),Tmunu(4,4),
     &     urcon(4),urcov(4),Rmunu(4,4),Gcon(4),Gcov(4)
      dimension gn(4,4),gv(4,4)
      common/metric/gn,gv
      common/accuracy/eps,epsbis,dvmin,tol,uminfac,dlogmax,
     &     gammaradceiling,itermax
      common/conserved/Gam,Gam1,en,en1,rhou0,s,Ttcov,Rtcov,BB,dt,
     &     FNUMARGSHARM,FNUMRESULTSHARM,WHICHVELRAMESHHARM,
     &     GAMMAMAXRAD,ERADLIMIT,toltry,tolallow,
     &     ARAD_CODE,OKAPPA_ES_CODE11,OKAPPA_FF_CODE11,OKAPPA_BF_CODE11
      common/funcMHDd/ffkap,eskap,arad,ucon,ucov,rho,bcon,bcov,bsq,
     &     Tmunu,E,Ehat,urcon,urcov,Rmunu,Tgas,Trad,B4pi,gamma,dtau,
     &     Gcon,Gcov
      double precision mymax,myabs,mydiv
      integer myisnan

c     Save u, compute rho, ucon, ucov, bcon, bcov, Tmunu

      u=prim(1)
      call solveucon(prim,ucon)
      call contocov(ucon,ucov)

      rho=rhou0/ucon(1) 

      if(rho.lt.0.0.or.myisnan(rho).eq.1) then
         jflag=1
      else
         jflag=0
      endif

      call calcbconbcov(BB,ucon,ucov,bcon,bcov,bsq)

      call calcTmunu(rho,u,bsq,Gam,ucon,ucov,bcon,bcov,Tmunu)

c     Compute normalized err3 from the three momentum equations and the
c     normalized err4 from all four equations

      err3=0.d0
      do i=2,4
         error(i)=Tmunu(1,i)-Ttcov(i)
         errornorm(i)=SMALL+myabs(mydiv(error(i),(myabs(error(i))+
     &        myabs(Tmunu(1,i))+myabs(Ttcov(i)))))
         err3=mymax(err3,errornorm(i))
      enddo

c     Use lab frame energy equation

      error(1)=Tmunu(1,1)-Ttcov(1)
      errornorm(1)=SMALL+myabs(mydiv(error(1),(myabs(error(1))+
     &     myabs(Tmunu(1,1))+myabs(Ttcov(1)))))
      err4=mymax(err3,errornorm(1))

#if(PRODUCTION==0)
c      write (*,*)
c      write (*,*) ' target Ttcov: ',(Ttcov(i),i=1,4)
c      write (*,*) ' current Ttcov: ',(Tmunu(1,i),i=1,4)
      write (*,*) 'B error: ',(error(i),i=1,4)
      write (*,*) 'B err3, err4: ',err3,err4
#endif

c     If u has an unreasonable value, set iflag=9 and reset u

      if (u.lt.0.d0) then
         iflag=9
#if(PRODUCTION==0)
         write (*,*) 'BSet iflag=',iflag
#endif
      else
         entropy=log((Gam1*u)**en/rho**en1)
         if (entropy.lt.((s/rhou0)-dlogmax)) iflag=9
      endif
      if (iflag.eq.9) prim(1)=mymax(prim(1),uminfac*rho)

#if(PRODUCTION==0)
c      write (*,*) ' error: ',(error(i),i=1,4),err3
#endif

      return 
      end



      subroutine funcMHD2(prim,error,errornorm,err4,err3,iflag,jflag)

c     This subroutine calculates errors for the MHD inversion problem
c     without radiation source term using the entropy equation

      implicit double precision (a-h,o-z)
      dimension prim(4),error(4),errornorm(4),Ttcov(4),Rtcov(4),BB(4),
     &     ucon(4),ucov(4),bcon(4),bcov(4),Tmunu(4,4)
      dimension gn(4,4),gv(4,4)
      common/metric/gn,gv
      common/accuracy/eps,epsbis,dvmin,tol,uminfac,dlogmax,
     &     gammaradceiling,itermax
      common/conserved/Gam,Gam1,en,en1,rhou0,s,Ttcov,Rtcov,BB,dt,
     &     FNUMARGSHARM,FNUMRESULTSHARM,WHICHVELRAMESHHARM,
     &     GAMMAMAXRAD,ERADLIMIT,toltry,tolallow,
     &     ARAD_CODE,OKAPPA_ES_CODE11,OKAPPA_FF_CODE11,OKAPPA_BF_CODE11
      double precision mymax,myabs,mydiv
      integer myisnan

c     Save u, compute rho, ucon, ucov, bcon, bcov, Tmunu

      u=prim(1)
      call solveucon(prim,ucon)
      call contocov(ucon,ucov)

      rho=rhou0/ucon(1) 

      if(rho.lt.0.0.or.myisnan(rho).eq.1) then
         jflag=1
      else
         jflag=0
      endif

      call calcbconbcov(BB,ucon,ucov,bcon,bcov,bsq)

      call calcTmunu(rho,u,bsq,Gam,ucon,ucov,bcon,bcov,Tmunu)

c     Compute normalized err3 from the three momentum equations and the
c     normalized err4 from all four equations

      err3=0.d0
      do i=2,4
         error(i)=Tmunu(1,i)-Ttcov(i)
         errornorm(i)=SMALL+myabs(mydiv(error(i),(myabs(error(i))+
     &        myabs(Tmunu(1,i))+myabs(Ttcov(i)))))
         err3=mymax(err3,errornorm(i))
      enddo

c     Use entropy equation

      entropy=ucon(1)*rho*log((Gam1*u)**en/rho**en1)
      error(1)=entropy-s
      errornorm(1)=SMALL+myabs(mydiv(error(1),(myabs(error(1))+
     &     myabs(entropy)+myabs(s))))
      err4=mymax(err3,errornorm(1))

#if(PRODUCTION==0)
c      write (*,*)
c      write (*,*) ' target entropy: ',s
c      write (*,*) ' current entropy: ',entropy
#endif

c     If u has an unreasonable value, set iflag=9 and reset u

      if (prim(1).lt.0.d0) then
         iflag=9
#if(PRODUCTION==0)
         write (*,*) 'CSet iflag=',iflag
#endif
         prim(1)=mymax(prim(1),uminfac*rho)
      endif

#if(PRODUCTION==0)
c      write (*,*) ' primitives: ',(prim(i),i=1,4)
c      write (*,*) ' target Ttcov: ',(Ttcov(i),i=2,4)
c      write (*,*) ' current Ttcov: ',(Tmunu(1,i),i=2,4)
      write (*,*) 'A error: ',(error(i),i=1,4)
      write (*,*) 'A err3, err4: ',err3,err4
#endif

      return
      end



      subroutine func3(primsave,prim3,error3,errornorm3,
     &     err4,err3,iflag,jflag,func)

c     This function is called by Newton3. It takes a 3-array with
c     primitives prim3(3), transfers the value primsave to prim4(1) and
c     the three velocity components to prim4(2-4) and computes the
c     4-array of errors error4(4). Then transfers appropriate elements
c     to error3(3) and returns this along with overall error err3.

      implicit double precision (a-h,o-z)
      dimension prim3(3),error3(3),prim4(4),error4(4),
     &     errornorm4(4),errornorm3(3)
      external func

c     Transfer primsave and prim3(3) to prim4(4)

      prim4(1)=primsave
      do i=1,3
         prim4(i+1)=prim3(i)
      enddo

c     Call appropriate function to calculate error4(4)

      call func(prim4,error4,errornorm4,err4,err3,iflag,jflag)

c     Transfer the momentum equation errors to error3(3) and return

      do i=1,3
         error3(i)=error4(i+1)
         errornorm3(i)=errornorm4(i+1)
      enddo

      return
      end



      subroutine funcrad1(prim,error,errornorm,err4,err3,iflag,jflag)

c     This subroutine calculates errors for the radiation inversion
c     problem including the radiation source term using the energy
c     equation

      implicit double precision (a-h,o-z)
      dimension prim(4),error(4),errornorm(4),Ttcov(4),Rtcov(4),BB(4),
     &     ucon(4),ucov(4),bcon(4),bcov(4),Tmunu(4,4),
     &     Tt(4),Rt(4),urcon(4),urcov(4),Rmunu(4,4),
     &     gn(4,4),gv(4,4),Gcon(4),Gcov(4)
      common/metric/gn,gv
      common/accuracy/eps,epsbis,dvmin,tol,uminfac,dlogmax,
     &     gammaradceiling,itermax
      common/conserved/Gam,Gam1,en,en1,rhou0,s,Ttcov,Rtcov,BB,dt,
     &     FNUMARGSHARM,FNUMRESULTSHARM,WHICHVELRAMESHHARM,
     &     GAMMAMAXRAD,ERADLIMIT,toltry,tolallow,
     &     ARAD_CODE,OKAPPA_ES_CODE11,OKAPPA_FF_CODE11,OKAPPA_BF_CODE11
      common/funcradd/ffkap,eskap,arad,ucon,ucov,rho,bcon,bcov,bsq,
     &     Tmunu,E,Ehat,urcon,urcov,Rmunu,Tgas,Trad,B4pi,gamma,dtau,
     &     Gcon,Gcov,u0,iuerr
      double precision mymax,myabs,mydiv
      integer myisnan
      integer radinvmod

c      ffkap=3.46764d-17
c      eskap=5.90799d5
c      arad=1.18316d17

      ffkap=OKAPPA_FF_CODE11
      eskap=OKAPPA_ES_CODE11
      arad=ARAD_CODE

c     Save u and solve for lab frame conserved quantities corresponding
c     to the given primitives: rho, ucon, ucov

      u=prim(1)
c     Calculate u^0 and rho
      call solveucon(prim,ucon)
      call contocov(ucon,ucov)

      rho=rhou0/ucon(1) 
#if(PRODUCTION==0)
c      write (*,*)
c      write (*,*) ' rhou0, rho, u: ',rhou0,rho,u
#endif

c     Calculate b^mu, b_mu, bsq

      call calcbconbcov(BB,ucon,ucov,bcon,bcov,bsq)

#if(PRODUCTION==0)
      write(*,*) 'BB bsq=',BB,bsq
#endif

c     Calculate the full gas stress-energy tensor T^mu_nu

      call calcTmunu(rho,u,bsq,Gam,ucon,ucov,bcon,bcov,Tmunu)

c     Evaluate the first row of the radiation stress-energy tensor:
c     R^0_mu

      do i=1,4
         Tt(i)=Tmunu(1,i)
         Rt(i)=Ttcov(i)+Rtcov(i)-Tt(i)
      enddo
#if(PRODUCTION==0)
c      write (*,*) ' Tt: ',(Tt(i),i=1,4)
c      write (*,*) ' Rt: ',(Rt(i),i=1,4)
#endif

c     Calculate the full R^mu_nu tensor

      call Rmunuinvert(Rt,E,urcon,urcov,ucon,Rmunu,radinvmod)
#if(PRODUCTION==0)
c      write (*,*) ' E, urcon, urcov: ',E,(urcon(i),urcov(i),i=1,4)
      write (*,*) ' upda conserved T^t_mu: ',(Tt(j),j=1,4)
      write (*,*) ' upda conserved R^t_mu: ',(Rt(j),j=1,4)
#endif

c     Calculate radiation energy density in the gas frame \hat{E}, and
c     the gas and radiation temperatures

      Ehat=0.d0
      do i=1,4
      do j=1,4
         Ehat=Ehat+Rmunu(i,j)*ucov(i)*ucon(j)
      enddo
      enddo
#if(PRODUCTION==0)
      write (*,*) ' \hat{E}: ',Ehat
#endif

      Trad=(Ehat/arad)**(0.25d0)
      Tgas=Gam1*u/rho
#if(PRODUCTION==0)
      write (*,*) ' Tgas, Trad: ',Tgas,Trad
#endif

      if(Trad.lt.0.0.or.myisnan(Trad).eq.1) then
         jflag=1
      else
         jflag=0
      endif

c     Calculate quantities needed to compute the radiation source term

      ff=ffkap*rho*rho/Tgas**(3.5d0)
      es=eskap*rho
      B4pi=arad*Tgas**4
#if(PRODUCTION==0)
c      write (*,*) ' ff, es, B4pi: ',ff,es,B4pi
#endif

c     The following side calculation is to estimate quantities that are
c     relevant for deciding whether the radiation source term can be
c     handled explicitly. Currently, all calculations are done
c     implicitly.

      alpha=1.d0/sqrt(-gn(1,1))
      gamma=ucon(1)*alpha

c     Decide how to set dtau: either dt/gamma or dt/ucon(1)

c      dtau=dt/gamma
      dtau=dt/ucon(1)

      chi1=ff*dtau*(1.d0+4.d0*B4pi/u)
      chi2=(ff+es)*dtau*(1.d0+Ehat/(rho+Gam*u))
#if(PRODUCTION==0)
c      write (*,*) ' chi1, chi2: ',chi1,chi2
#endif

c     Compute the radiation source term G^mu and G_mu directly in the
c     lab frame (formula taken from Jon's draft of the paper, with some
c     corrections)

      do i=1,4
         Gcon(i)=-(ff*B4pi+es*Ehat)*ucon(i)
      do j=1,4
         Gcon(i)=Gcon(i)-(ff+es)*Rmunu(i,j)*ucon(j)
      enddo
      enddo

      call contocov(Gcon,Gcov)
#if(PRODUCTION==0)
c      write (*,*)
c      write (*,*) ' lab Gcon: ',(Gcon(i),i=1,4)
c      write (*,*) ' lab Gcov: ',(Gcov(i),i=1,4)
#endif

c     Compute normalized err3 from the three momentum equations and the
c     normalized err4 from all four equations

      err3=0.d0
      do i=2,4
         error(i)=Tt(i)-Ttcov(i)-Gcov(i)*dt
         errornorm(i)=SMALL+myabs(mydiv(error(i),(myabs(Gcov(i)*dt)+
     &        myabs(Tt(i))+myabs(Ttcov(i)))))
         err3=mymax(err3,errornorm(i))
      enddo

c     Calculate error(1) from lab frame energy equation

         error(1)=Tt(1)-Ttcov(1)-Gcov(1)*dt
         errornorm(1)=myabs(mydiv(error(1),
     &        (myabs(Gcov(1)*dt)
     &        +myabs(Tt(1))+myabs(Ttcov(1))
     &        +myabs((Gam*u+bsq)*ucon(1)*ucov(1))
     &        )))
         err4=mymax(err3,errornorm(1))

#if(PRODUCTION==0)
c      write (*,*)
c      write (*,*) ' target Ttcov: ',(Ttcov(i),i=1,4)
c      write (*,*) ' current Ttcov: ',(Tt(i),i=1,4)
      write (*,*) 'C error: ',(error(i),i=1,4)
      write (*,*) 'C error2: ',(myabs(error(i)),i=1,4)
      write (*,*) 'C errornorm: ',(errornorm(i),i=1,4)
      write (*,*) 'C errornorm2: ',
     & (myabs(Gcov(1)*dt)+myabs(Tt(1))+myabs(Ttcov(1)))
      write (*,*) 'C err3, err4: ',err3,err4
#endif
c     If u has an unreasonable value, set iflag=9 and reset u

      if (u.lt.0.d0) then
         iflag=9
#if(PRODUCTION==0)
         write (*,*) 'DSet iflag=',iflag
#endif
      else
         entropy=log((Gam1*u)**en/rho**en1)
         if (entropy.lt.((s/rhou0)-dlogmax)) iflag=9
      endif
      if (iflag.eq.9) prim(1)=mymax(prim(1),uminfac*rho)

      u0=prim(1)

      return
      end



      subroutine funcrad2(prim,error,errornorm,err4,err3,iflag,jflag)

c     This subroutine calculates errors for the radiation inversion
c     problem including the radiation source term using the entropy
c     equation

      implicit double precision (a-h,o-z)
      dimension prim(4),error(4),errornorm(4),Ttcov(4),Rtcov(4),BB(4),
     &     ucon(4),ucov(4),bcon(4),bcov(4),Tmunu(4,4),
     &     Tt(4),Rt(4),urcon(4),urcov(4),Rmunu(4,4),
     &     gn(4,4),gv(4,4),Gcon(4),Gcov(4)
      common/metric/gn,gv
      common/conserved/Gam,Gam1,en,en1,rhou0,s,Ttcov,Rtcov,BB,dt,
     &     FNUMARGSHARM,FNUMRESULTSHARM,WHICHVELRAMESHHARM,
     &     GAMMAMAXRAD,ERADLIMIT,toltry,tolallow,
     &     ARAD_CODE,OKAPPA_ES_CODE11,OKAPPA_FF_CODE11,OKAPPA_BF_CODE11
      common/funcradd/ffkap,eskap,arad,ucon,ucov,rho,bcon,bcov,bsq,
     &     Tmunu,E,Ehat,urcon,urcov,Rmunu,Tgas,Trad,B4pi,gamma,dtau,
     &     Gcon,Gcov,u0,iuerr
      double precision mymax,myabs,mydiv
      integer myisnan
      integer radinvmod

c      ffkap=3.46764d-17
c      eskap=5.90799d5
c      arad=1.18316d17

      ffkap=OKAPPA_FF_CODE11
      eskap=OKAPPA_ES_CODE11
      arad=ARAD_CODE

c     Save u and solve for lab frame conserved quantities corresponding
c     to the given primitives: rho, ucon, ucov

      u=prim(1)
      if (prim(1).lt.0.d0) then
         iflag=9
#if(PRODUCTION==0)
         write (*,*) 'FSet iflag=',iflag
#endif
      endif

c     Calculate u^0 and rho
      call solveucon(prim,ucon)
      call contocov(ucon,ucov)

      rho=rhou0/ucon(1) 
#if(PRODUCTION==0)
c      write (*,*)
c      write (*,*) ' rho, u: ',rho,u
#endif

c     Calculate b^mu, b_mu, bsq

      call calcbconbcov(BB,ucon,ucov,bcon,bcov,bsq)

c     Calculate the full gas stress-energy tensor T^mu_nu

      call calcTmunu(rho,u,bsq,Gam,ucon,ucov,bcon,bcov,Tmunu)

c     Evaluate the first row of the radiation stress-energy tensor:
c     R^0_mu

      do i=1,4
         Tt(i)=Tmunu(1,i)
         Rt(i)=Ttcov(i)+Rtcov(i)-Tt(i)
      enddo

c     Calculate the full R^mu_nu tensor

      call Rmunuinvert(Rt,E,urcon,urcov,ucon,Rmunu,radinvmod)

c     Calculate radiation energy density in the gas frame \hat{E}, and
c     the gas and radiation temperatures

      Ehat=0.d0
      do i=1,4
      do j=1,4
         Ehat=Ehat+Rmunu(i,j)*ucov(i)*ucon(j)
      enddo
      enddo
#if(PRODUCTION==0)
c      write (*,*) ' \hat{E}: ',Ehat
#endif
      Trad=(Ehat/arad)**(0.25d0)
      Tgas=Gam1*u/rho
#if(PRODUCTION==0)
c      write (*,*) ' Tgas, Trad: ',Tgas,Trad
#endif

c     Calculate quantities needed to compute the radiation source term

      ff=ffkap*rho*rho/Tgas**(3.5d0)
      es=eskap*rho
      B4pi=arad*Tgas**4
#if(PRODUCTION==0)
c      write (*,*) ' ff, es, B4pi: ',ff,es,B4pi
#endif

      if(Trad.lt.0.0.or.myisnan(Trad).eq.1) then
         jflag=1
      else
         jflag=0
      endif

c     The following side calculation is to estimate quantities that are
c     relevant for deciding whether the radiation source term can be
c     explicitly. Currently, all calculations are done implicitly.

      alpha=1.d0/sqrt(-gn(1,1))
      gamma=ucon(1)*alpha

c     Decide how to set dtau: either dt/gamma or dt/ucon(1)

c      dtau=dt/gamma
      dtau=dt/ucon(1)

      Ghatdtau=ff*(Ehat-B4pi)*dtau

      chi1=ff*dtau*(1.d0+4.d0*B4pi/u)
      chi2=(ff+es)*dtau*(1.d0+Ehat/(rho+Gam*u))
#if(PRODUCTION==0)
c      write (*,*) ' chi1, chi2: ',chi1,chi2
#endif

c     Compute the radiation source term G^mu and G_mu directly in the
c     lab frame (formula taken from Jon's draft of the paper, with some
c     corrections)

      do i=1,4
         Gcon(i)=-(ff*B4pi+es*Ehat)*ucon(i)
      do j=1,4
         Gcon(i)=Gcon(i)-(ff+es)*Rmunu(i,j)*ucon(j)
      enddo
      enddo

      call contocov(Gcon,Gcov)
#if(PRODUCTION==0)
c      write (*,*)
c      write (*,*) ' lab Gcon: ',(Gcon(i),i=1,4)
c      write (*,*) ' lab Gcov: ',(Gcov(i),i=1,4)
#endif
c     Compute normalized err3 from the three momentum equations and the
c     normalized err4 from all four equations

      err3=0.d0
      do i=2,4
         error(i)=Tt(i)-Ttcov(i)-Gcov(i)*dt
         errornorm(i)=SMALL+myabs(mydiv(error(i),(myabs(Gcov(i)*dt)+
     &        myabs(Tt(i))+myabs(Ttcov(i)))))
         err3=mymax(err3,errornorm(i))
      enddo

c     Calculate error(1) corresponding to the lab frame entropy equation

      entropy=ucon(1)*rho*log((Gam1*prim(1))**en/rho**en1)
      Gt=Gcov(1)*dt
      error(1)=entropy-s-Gt
      errornorm(1)=SMALL+myabs(mydiv(error(1),(myabs(entropy)+
     &     myabs(s)+myabs(Gt))))

c      entropy=rho*log((Gam1*prim(1))**en/rho**en1)
c      Gt=Ghatdtau
c      error(1)=entropy-(s/ucon(1))-Gt
c      errornorm(1)=SMALL+myabs(mydiv(error(1),(myabs(entropy)+
c     &     myabs(s/ucon(1))+myabs(Gt))))

c     Alternatively, use the fluid frame entropy equation

c      entropy=ucon(1)*rho*log((Gam1*prim(1))**en/rho**en1)
c      Gdtau=ff*(Ehat-B4pi)*dtau
c      Gt=Gdtau
c      error(1)=entropy-s-Gdtau
      err4=mymax(err3,errornorm(1))
      u0=prim(1)
#if(PRODUCTION==0)
c      write (*,*) ' entropy, s, Gdtau: ',entropy,s,Gdtau
c      write (*,*) ' error(1), err3, err4 ',error(1),err3,err4
c      write (*,*)
c      write (*,*) ' target Ttcov: ',(Ttcov(i),i=1,4)
c      write (*,*) ' current Ttcov: ',(Tt(i),i=1,4)
      write (*,*) 'E error: ',(error(i),i=1,4)
      write (*,*) 'E errornorm: ',(errornorm(i),i=1,4)
      write (*,*) 'E err3, err4: ',err3,err4
#endif


      return
      end



      SUBROUTINE lubksb(a,n,np,indx,b)

c     Matrix inversion subroutine 1 from Numerical Recipes

      implicit double precision (a-h,o-z)
      INTEGER n,np,indx(n)
      REAL*8 a(np,np),b(n)
      INTEGER i,ii,j,ll
      REAL*8 sum
      ii=0
      do 12 i=1,n
        ll=indx(i)
        sum=b(ll)
        b(ll)=b(i)
        if (ii.ne.0)then
          do 11 j=ii,i-1
            sum=sum-a(i,j)*b(j)
11        continue
        else if (sum.ne.0.d0) then
          ii=i
        endif
        b(i)=sum
12    continue
      do 14 i=n,1,-1
        sum=b(i)
        do 13 j=i+1,n
          sum=sum-a(i,j)*b(j)
13      continue
        b(i)=sum/a(i,i)
14    continue
      return
      END



      SUBROUTINE ludcmp(a,n,np,indx,d,retval)

c     Matrix inversion subroutine 2 from Numerical Recipes

      implicit double precision (a-h,o-z)
      INTEGER n,np,indx(n),NMAX
      REAL*8 d,a(np,np),TINY
      PARAMETER (NMAX=500,TINY=1.0d-80)
      INTEGER i,imax,j,k
      REAL*8 aamax,dum,sum,vv(NMAX)
      double precision mymax,myabs,mydiv
      integer myisnan
      d=1.d0
      imax=-100
      do 12 i=1,n
        aamax=0.d0
        do 11 j=1,n
          if (myabs(a(i,j)).gt.aamax) aamax=myabs(a(i,j))
11      continue
        if (aamax.eq.0.d0) then
#if(PRODUCTION==0)
           write (*,*) 'singular matrix in ludcmp'
#endif
           retval=1
           return
        endif
        vv(i)=1.d0/aamax
12    continue
      do 19 j=1,n
        do 14 i=1,j-1
          sum=a(i,j)
          do 13 k=1,i-1
            sum=sum-a(i,k)*a(k,j)
13        continue
          a(i,j)=sum
14      continue
        aamax=0.d0
        do 16 i=j,n
          sum=a(i,j)
          do 15 k=1,j-1
            sum=sum-a(i,k)*a(k,j)
15        continue
          a(i,j)=sum
          dum=vv(i)*myabs(sum)
          if (dum.ge.aamax) then
            imax=i
            aamax=dum
          endif
16      continue
        if (imax.gt.n.or.imax.lt.1.or.imax.eq.-100) then
#if(PRODUCTION==0)
           write (*,*) 'Bad imax in ludcmp: ',imax
#endif
           retval=2
           return
        endif
        if (j.ne.imax)then
          do 17 k=1,n
            dum=a(imax,k)
            a(imax,k)=a(j,k)
            a(j,k)=dum
17        continue
          d=-d
          vv(imax)=vv(j)
        endif
        indx(j)=imax
        if(a(j,j).eq.0.d0)a(j,j)=TINY
        if(j.ne.n)then
          dum=1.d0/a(j,j)
          do 18 i=j+1,n
            a(i,j)=a(i,j)*dum
18        continue
        endif
19    continue
      retval=0
      return
      END



      FUNCTION rtsafe(funcd,x1,x2,fl,fr,xacc)
      implicit double precision (a-h,o-z)
      INTEGER MAXIT
      REAL*8 rtsafe,x1,x2,xacc
      EXTERNAL funcd
      PARAMETER (MAXIT=100)
      INTEGER j
      REAL*8 df,dx,dxold,f,fh,fl,temp,xh,xl
      double precision mymax,myabs,mydiv
      integer myisnan
c      call funcd(x1,fl,df)
c      call funcd(x2,fh,df)
#if(PRODUCTION==0)
c      if((fl.gt.0.d0.and.fh.gt.0.d0).or.(fl.lt.0.d0.and.fh.lt.0.d0))
c     &     write (*,*) 'root is not bracketed in rtsafe! '
#endif
c      if(fl.eq.0.d0)then
c        rtsafe=x1
c        return
c      else if(fh.eq.0.d0)then
c        rtsafe=x2
c        return
      if(fl.lt.0.d0)then
        xl=x1
        xh=x2
      else
        xh=x1
        xl=x2
      endif
      rtsafe=.5d0*(x1+x2)
      call funcd(rtsafe,f,df)
      if (f.lt.0.d0) then
         xl=rtsafe
      else
         xh=rtsafe
      endif
      dxold=myabs(xh-xl)
      dx=dxold
      do 11 j=1,MAXIT
        if(((rtsafe-xh)*df-f)*((rtsafe-xl)*df-f).ge.0.d0.or. 
     &        myabs(2.d0*f).gt.myabs(dxold*df) ) then
          dxold=dx
          dx=0.5d0*(xh-xl)
          rtsafe=xl+dx
          if(xl.eq.rtsafe)return
        else
          dxold=dx
          dx=f/df
          temp=rtsafe
          rtsafe=rtsafe-dx
          if(temp.eq.rtsafe)return
        endif
        if(myabs(dx).lt.xacc) return
        call funcd(rtsafe,f,df)
        if(f.lt.0.d0) then
          xl=rtsafe
        else
          xh=rtsafe
        endif
11    continue
#if(PRODUCTION==0)
      write (*,*) ' rtsafe has exceeded maximum iterations '
#endif
      return
      END
