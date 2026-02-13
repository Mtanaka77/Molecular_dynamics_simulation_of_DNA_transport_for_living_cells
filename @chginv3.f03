!*-------------------------------------------- Update: 2025/01/17 --*
!*                                                                  *
!*   ## Charge Inversion of Macroions by Counter/Co-Ions ##         *
!*                                                                  *
!*   Author: Motohiko Tanaka, Ph.D., Chikusa, Nagoya 464, Japan     *
!*                                                                  *
!*  References                                                      *
!*  1. M.Tanaka and A.Yu Grosberg, J.Chem.Phys., vol.115,567 (2001).*
!*  2. M.Tanaka and A.Yu Grosberg, Euro.Phys.J., E7, 371 (2002).    *
!*  3. M.Tanaka, Phys.Reviews., E68, 061501 (2003).                 *
!*                                                                  *
!*  Equations of motion:                                            *
!*    Original r,v,t,m are normalized as R_i,V_i,hat(t),M_i         *
!*                                                                  * 
!*        dV_i   (t*e)^2 qq'R   t^2  48*eps'R   r0       1 r0       *
!*    M_i ---- = -------*---- + ----*-------- [(--)^12 - -(--)^6]   * 
!*         dt     ma^3   R^3    ma^2   R^2      R        2 R        *
!*                                                                  *
!*                t^2*e  Edc(V/cm)                                  *
!*              + ------ ---------  -(mue0*a^2/t) mue*ag(i)*V_i     *
!*                  ma     300         Langevin thermostat          *
!*    dR_i                                                          * 
!*    ---- = V_i                                                    *
!*     dt                                                           * 
!*                                                                  * 
!*    ccel = 48.d0*pref_eps*epsav*snt*(snt-0.5d0)*sqrt(rsi/r2)      *
!*    forceV/r = prefactor*ch(i)*ch(j)* &                           *
!*                    (erfc/r +2*alpha/sqrtpi)*exp(-ar**2)/r2       *
!*                                                                  *  
!*------------------------------------------------------------------*
!*  Main and subroutines (ca. 5000 lines):                          *
!*    program charge_inv                                            *
!*                                                                  *  
!*    RUN_MD                                                        *
!*    Equations of motion - moldyn                                  *
!*    Short-range - realteil_s                                      *
!*    Long-range - p3m_perform, p3m_init, perform_aliasing_sums,... *
!*      Note: The p3m routines in C were written by Dr.M.Deserno    *
!*       and Dr.C.Holm, University of Mainz, Germany in Jan.1999,   *
!*       and rewritten by Fortran 90 by M.Tanaka in Dec.1999.       *
!*                                                                  *
!*    init, ggauss, fun                                             *
!*    LPLOT1/HPLOT1, lplots, ppl3da,...                             *
!*    gopen graphic package                                         *
!*                                                                  *
!********************************************************************
!*  To get a free format (f90, f03), convert f77 into:              *
!*    :%s/^c/!/ and "tr 'A-Z' 'a-z' <@chginv3.f >@chginv3.f03"      *
!*  Also, "use, intrinsic :: iso_c_binding"                         *
!*  The @ character is not permitted in the Intel (LX) system.      *
!*                                                                  *
!*$ mpif90 -mcmodel=medium -fPIC -O2 -o a.out @chginv3.f03 -I/opt/fftw3/include -L/opt/fftw3/lib -lfftw3 &> log *
!********************************************************************
!*-- 12/23/1999 --------------------------------------- 7/07/2001 --*
!
      program  charge_inv
!
      use, intrinsic :: iso_c_binding
      implicit none
!
      include 'param_inv3.h'
      include 'mpif.h'
!
      real(C_DOUBLE) ctime1,ctime2,ctime
      integer(C_int) rank,size,ierror,ipar,igrp, &
                     nframe,cl_first
!
      integer(C_int) io_pe,num_proc
      common/sub_proc/ io_pe,num_proc
!
      call mpi_init (ierror)
      call mpi_comm_rank (mpi_COMM_WORLD,rank,ierror)
      call mpi_comm_size (mpi_COMM_WORLD,size,ierror)
!
      num_proc= size
!
      ipar = 1 + rank             ! PE #
      igrp = 1 + rank/size        ! Group #
!
      io_pe= 0
      if(ipar.eq.1) io_pe= 1      ! plot if io_pe= 1
!
      nframe= 4
      call FLOPEN (nframe,igrp)
!
      if(io_pe.eq.1) then
        write(11,*) 'My process is #',ipar,'of group #',igrp
        write(11,*) '  rank, size=',rank,size
        write(11,*)
      end if
!********************************************
!
      cl_first= 1
      call clocks (ctime1,cl_first)
!
! -------------------------------------------------
      call Run_MD (ipar,igrp,size)
! -------------------------------------------------
!
      cl_first= 2
      call clocks (ctime2,cl_first)
!
      ctime= ctime2 -ctime1
      if(io_pe.eq.1) then
        write(11,*)
        write(11,'(/,"*ipar, ctime(sec)=",i3,f10.3)') ipar,ctime
      end if
!*
      call mpi_finalize  (ierror)
!
      stop
      end program charge_inv
!
!
!----------------------------------------------------------------------
      subroutine clocks (tt,cl_first)
!----------------------------------------------------------------------
!  For this job, t1: beginning, t2: current
!
      use, intrinsic :: iso_c_binding
      implicit none
!
      include 'param_inv3.h'
      include 'mpif.h'
!
      real(C_DOUBLE)  tt,unif1,unif2
      integer(C_int)  cl_first,ierror
      integer(Selected_Int_Kind(18)),save :: t1, t2, t_rate, t_max, diff
!
      integer(C_int) io_pe,num_proc,it,is
      common/sub_proc/ io_pe,num_proc
      common/parm1/ it,is
!
      if(cl_first.eq.1) then
        call system_clock (t1)   ! 開始時を記録
      end if
!
      call system_clock (t2, t_rate, t_max)   ! 終了時を記録
!
      if ( t2 < t1 ) then
        diff = (t_max - t1) + t2 + 1
      else
        diff = t2 - t1
      endif
!
      unif2= diff/dble(t_rate)
      call mpi_allreduce (unif2,unif1,1,mpi_real8,mpi_sum, &
                          mpi_comm_world,ierror)
      tt= unif1/num_proc
!
      return
      end subroutine clocks
!
!
!----------------------------------------------------------------------
      subroutine Run_MD (ipar,igrp,size)
!----------------------------------------------------------------------
      use, intrinsic :: iso_c_binding
      implicit none
!
      include    'param_inv3.h'
      include    'mpif.h'
!
      integer(C_int) ipar,igrp,size,io_pe,num_proc
      common/sub_proc/ io_pe,num_proc
!
      real(C_DOUBLE),dimension(npqr0) :: &  !! include all atoms
                    xg,yg,zg,vx,vy,vz,ch,am,ag,ep
!
      real(C_DOUBLE) t_unit,a_unit,w_unit,e_unit,  &
                     kjoule,kcal,mol,kbT
      common/units/ t_unit,a_unit,w_unit,e_unit
      common/unit2/ kjoule,kcal,mol,kbT
!
      integer(C_int) np,nq,nr,npqr
!*
!  Single precision for the diagnostic routines.
      real(C_DOUBLE) pi,dt,rbmax,vth,vth1,vth2,vth3,tmax,    &
                     xmax,ymax,zmax,t8,                      &
                     tmax0,qfrac,Qcore,Rmac,Wmac,Zcp,Zcn,Edc
      common/time3/  t8
      integer(C_int) it,is,nsg,nseg,istop,iwa,iwb,iwc
!
      character(len=8) label,cdate*10,ctime*8
      real(C_float)  t,xleng,phi,tht,dtwr,dtwr2,cptot
      common/headr1/  label,cdate,ctime
      common/headr2/  t,xleng
!
      common/parm1/ it,is
      common/parm2/ pi,dt,rbmax,vth,vth1,vth2,vth3,tmax
      common/parm3/ xmax,ymax,zmax,phi,tht,dtwr,dtwr2
      common/parm9/ cptot
      common/psegm/ nsg(30),nseg
      common/abterm/ istop
!
      common/imemo/ iwa,iwb,iwc
      common/cntion/ qfrac,Qcore,Rmac,Wmac,Zcp,Zcn
      common/EBfild/ Edc
!
      real(C_DOUBLE)    rcutpme,rcutLJ,r_fene,k_fene,r_fene2,     &
                        alpha,length,Bjerrum,rcutlj2,temperat,    &
                        dmesh,prefactor,pref_eps,econv,           &
                        scale_c,surreps,gamma
!
      common/CONFdataR/ rcutpme,rcutLJ,r_fene
      common/fenepara/  k_fene,r_fene2
      common/elsta/     Bjerrum
      common/cutoffLJ/  rcutlj2,temperat
      common/ewald1/    alpha,length,dmesh,prefactor,pref_eps,econv
!
      real(C_float),dimension(nhist) :: &
                    ekin,ekn2,ekn3,ecr,elj,ep3m,pot,time,          &
                    vxio,vxco,vxct,vxmc,vxtot,xmcc,vxwat,etot,     &
                    xmac,ymac,zmac,xcon,ycon,zcon,xcoi,ycoi,zcoi,  &
                    xwat,ywat,zwat,dmcc
      common/ehist/ ekin,ekn2,ekn3,ecr,elj,ep3m, pot,time,         &
                    vxio,vxco,vxct,vxmc,vxtot,xmcc,vxwat,etot,     &
                    xmac,ymac,zmac,xcon,ycon,zcon,xcoi,ycoi,zcoi,  &
                    xwat,ywat,zwat,dmcc
!
!--------------------------------------------------------------
!* Read basic data in /READ_CONF/, constants for realteil and p3m.
!    >> Bjerrum
!    >> dt, length, alpha, rcutpme, rcutlj
!    >> r_fene, k_fene 
!    >> cptot, dtwr, dtwr2, rmax.
!
!   cptot (in min.)
!     cptot= 570.0  <--- 10 H
!     dtwr=    50.
!     dtwr2=  500.
!     rmax =   50.  <--- Initial size of Polymer.
!
!    Non-neutral PA................ ifqq= 1
!                                        qfrac = dq/q
!    Valence of counterions........ Zcp, Zcn
!      Salt (Zcp, -Zcn) after sufficient counterions.
!*********************
!     ifqq= 1
!     qfrac= 0.30
!     Qcore= -28.d0
!     Rmac = 3.d0
!     Wmac = 300.d0
!     Zcp=   3.
!     Zcn=  -1.
!*********************
!*  rbmax: Maximum bond length for /sprmul/.
!
      label='P3M.CINV'
      call date_and_time_7 (cdate,ctime)
!
      if(io_pe.eq.1) then
        write(11,600) label,cdate,ctime
  600   format('<< Charge inversion of moving macroions -- P3M >> ', &
               a8,/,' Date=',a10,'  time=',a8,//)
      end if
!
!  Four basic quantities of one's in CGS
      t_unit = 1.0000d-15  ! 1 fs
      a_unit = 1.0000d-08  ! 1 Angstrom
      w_unit = 1.6726d-24  ! proton, g
      e_unit = 4.8033d-10  ! charge, esu
!
      if(io_pe.eq.1) then
        write(11,30) 
   30   format('*Four basic quantities must be used in equations...',/,&
               '  t_unit = 1.0000d-15',/,     &
               '  a_unit = 1.0000d-08',/,     &
               '  w_unit = 1.6605d-24',/,     &
               '  e_unit = 4.8033d-10',/)
      end if
!
      rbmax= 1.5d0 
!
!*  nseg :  defined in /READ_CONF/.
!       -> np, alpha
      call READ_CONF (xg,yg,zg,vx,vy,vz,ch,am,ag,np)
!
      kjoule = 1.00d+10             ! in erg
      kcal = 4.1868d0 *kjoule       ! 1 cal= 4.18 J/cal, in erg
      mol  = 6.0220d+23
      kbT  = 1.3807d-16 *temperat   ! READ_CONF, erg in Kelvin
!
!  Four basic quantities of one's in CGS
!     t_unit = 1.0000d-15  ! 0.01 ps
!     a_unit = 1.0000d-08  ! 1 Angstrom
!     w_unit = 1.6605d-24  ! m_proton
!     e_unit = 4.8033d-10  ! charge esu
!
!!    eps= k_B*T for LJ eq. and epsilon=78 for C.coupling are different !
!!    scale_c  = (t_unit*e_unit)**2 /(w_unit*a_unit**3)
!!    pref_eps = t_unit**2/(w_unit*a_unit**2)
!
!*     ^  dv    (t*e)^2 qq'R   t^2  48*eps'R   r0       1 r0        *
!*     m ---- = -------*---- + ----*-------- [(--)^12 - -(--)^6]    * 
!*        dt     ma^3   r^3    ma^2   r*r      r        2 r         *
!*                                                                  *
!*                t^2*e  Edc0(V/cm)                                 *
!*              + ------ ----------  -(mue0* a^2/t) mue*ag(i)*v     *
!*                  ma     300         Langevin thermostat          *
!
!       ccel = 48.d0*pref_eps*epsav*snt*(snt-0.5d0)*sqrt(rsi/r2)
!       forceV = prefactor*ch(i)*ch(j)*                            &
!                          (erfc/r +2*alpha/sqrtpi)*exp(-ar**2)/r2
!
!       econv= t_unit**2*e_unit /(w_unit*a_unit *299.98d0)
!       Edc  = econv * Edc0 
!         maximum drift at 1.2d4 statV/cm= 3.8d6 V/cm 
!           Edc0= 1.d6 V/cm -> Edc= econv*Edc0= 1.27d-04 
!
!  it is surrounded by water, which suppresses by eps= 78
!       epsil_j = 3.1879d-14  ! erg
!
      scale_c  = (t_unit*e_unit)**2 /(w_unit*a_unit**3)
      pref_eps = t_unit**2/(w_unit*a_unit**2)
      surreps  = 78.d0
!
      gamma= scale_c /surreps 
      prefactor = gamma  !<- include surreps
!
      dmesh= mesh
      xleng= length     
!
      if(io_pe.eq.1) then
        write(11,230) length,scale_c,gamma,pref_eps, &
                      prefactor*abs(3.d0*1.d0),      &
                      pref_eps*3.19d-14/3.2,         &
                      Edc
  230   format('*length =',1pd12.5,/,          &
               ' scale_c (in vacuum)     =',d12.5,/, &
               ' prefactor (surreps=78)  =',d12.5,/, &
               ' pref_eps (LJ force)     =',d12.5,/, &
               '   prefactor*ch(i)*ch(j) =',d12.5,/, &
               '   pref_eps*3.19d-14/3.2 =',d12.5,/, &
               ' Edc = econv*Edc0 =',d12.5,/)
      end if
!
!--------------------------------------------------------------
      istop = 0    ! Signal for termination: istop= 1
!--------------------------
!****************************************
!*  Prepare for graphic output.         *
!****************************************
!*  For 3-D plot of particles /ppl3d/.
!
      phi= -60.
      tht=  15.
!
      pi = 4.d0*atan(1.d0)  !<- generic name atan()
      call ggauss 
!
!-----------------------------------------------------
!*  System size (0., xmax), and Ewald sum parameter.
!
      xmax=  length
      ymax=  length
      zmax=  length
!
      if(io_pe.eq.1) then
        write(11,*) '# 3-D cube: xmax, alpha (~/length) = ',xmax,alpha
      end if
!
!************************************
!*   Step 0 : Initial loading.      *
!************************************
!* Generate all the charge sequences (cases) at first.
!  imix: orikaesi.... different random numbers.
!
!************************************
!*   Step 1 : Molecular dynamics.   *
!************************************
!       -> xg - ep, nq - npqr
!
      call init (xg,yg,zg,vx,vy,vz,ch,am,ag,ep,np,nq,nr,npqr)
!     *************************
      call Initialisierung
!
        if(io_pe.eq.1) then
          write(11,*) '>> np,nq,nr=',np,nq,nr
        end if
!
      call p3m_init
!     *************
!
      if(kstart.eq.0) then
        if(io_pe.eq.1) then
          write(11,*) '*Now kstart= 0'
        end if
!
!* Restart data.
      else
        if(io_pe.eq.1) then
          write(11,667) praefixc//'.12'//suffix1
  667     format('Restart data are loaded from FT12.....', &
                 ' file=',a34,/)
        end if
!                                     now start with suffix1 
        OPEN (unit=12,file=praefixc//'.12'//suffix1,           &
                               status='old',form='unformatted')
!
        read(12) it,is,size,np,nq,nr,npqr
        read(12) xg,yg,zg,vx,vy,vz,ch,am,ag,ep
        read(12) pi,dt,rbmax,vth,vth1,vth2,vth3,tmax0,t8 
        read(12) t,phi,tht,dtwr,dtwr2 
        read(12) ekin,ekn2,ekn3,ecr,elj,ep3m,pot,time,           &
                 vxio,vxco,vxct,vxmc,vxtot,xmcc,vxwat,etot,      &
                 xmac,ymac,zmac,xcon,ycon,zcon,xcoi,ycoi,zcoi,   &
                 xwat,ywat,zwat,dmcc                   !*added
        read(12) iwa,iwb,iwc
        read(12) nsg,nseg
        read(12) alpha,xmax,ymax,zmax
        read(12) qfrac,Qcore,Rmac,Wmac,Zcp,Zcn
        read(12) Edc
        close(12)
!
        if(io_pe.eq.1) then
          write(11,*) 'Present time t, is=',t,is
        end if
      end if
!
!************************************
!*   Step 1 : Molecular dynamics.   *
!************************************
!
      if(io_pe.eq.1) then
        write(11,603) alpha
  603   format('*alpha=',1pd15.6)
      end if
!
!*-------------------------------------------------------------
      call moldyn (xg,yg,zg,vx,vy,vz,ch,am,ag,ep,ipar,igrp, &
                   size,np,nq,nr,npqr)
!*-------------------------------------------------------------
!************************************************************
!* Step 3: Restart data.
!
      if(io_pe.eq.1) then                 !  finished
        OPEN (unit=12,file=praefixc//'.12'//suffix2,        &
                                         form='unformatted')
        write(12) it,is,size,np,nq,nr,npqr
        write(12) xg,yg,zg,vx,vy,vz,ch,am,ag,ep
        write(12) pi,dt,rbmax,vth,vth1,vth2,vth3,tmax,t8 
        write(12) t,phi,tht,dtwr,dtwr2
        write(12) ekin,ekn2,ekn3,ecr,elj,ep3m,pot,time,           &
                  vxio,vxco,vxct,vxmc,vxtot,xmcc,vxwat,etot,      &
                  xmac,ymac,zmac,xcon,ycon,zcon,xcoi,ycoi,zcoi,   &
                  xwat,ywat,zwat,dmcc
        write(12) iwa,iwb,iwc
        write(12) nsg,nseg
        write(12) alpha,xmax,ymax,zmax
        write(12) qfrac,Qcore,Rmac,Wmac,Zcp,Zcn
        write(12) Edc
        close(12)
!------------------------------
        open (unit=77,file=praefixc//'.77'//suffix2//'.ps',      & 
            status='unknown',position='append',form='formatted')
!
        call lplots
!------------------------------
        call plote
        close(77)
      end if
!**********************************************************
!
      return
      end subroutine Run_MD
!
!
!------------------------------------------------------------------------
      subroutine date_and_time_7 (date_now,time_now)
!------------------------------------------------------------------------
      use, intrinsic :: iso_c_binding 
      implicit none
!
      integer, dimension(8) :: ipresent_time
      character(len=10) :: date_now
      character(len=8)  :: time_now

      call date_and_time(values=ipresent_time)

      write(time_now,'(i2,":",i2,":",i2)') ipresent_time(5:7)
      write(date_now,'(i4,"/",i2,"/",i2)') &
               ipresent_time(1),ipresent_time(2),ipresent_time(3)
!
      return
      end subroutine date_and_time_7
!
!
!----------------------------------------------------------------------
      subroutine FLOPEN (nframe,igrp)
!----------------------------------------------------------------------
      use, intrinsic :: iso_c_binding
      implicit none
!
      include  'param_inv3.h'

      integer(C_int) nframe,igrp,io_pe,num_proc
      common/sub_proc/ io_pe,num_proc
!
      if(io_pe.eq.1) then
        OPEN (unit=11,file=praefixc//'.06'//suffix2)
!
        OPEN (unit=77,file=praefixc//'.77'//suffix2//'.ps')
          call gopen (nframe)
        close(77)
!
        write(11,*) 'igrp=',igrp,' uses ',            &
                     praefixs//'_config.START'//suffix0
      end if
!
      return
      end subroutine FLOPEN
!
!
!------------------------------------------------------------------
      subroutine moldyn (xg,yg,zg,vx,vy,vz,ch,am,ag,ep,ipar,igrp, &
                         size,np,nq,nr,npqr)
!------------------------------------------------------------------
!*  Double precision.
      use, intrinsic :: iso_c_binding
      implicit none
!
      include    'param_inv3.h'
      include    'mpif.h'
!
      real(C_DOUBLE),dimension(npqr0) :: &
                     xg,yg,zg,vx,vy,vz,ch,am,ag,ep,  &
                     frx,fry,frz,fsx,fsy,fsz
      integer(C_int) np,nq,nCLp,nr,npqr,ipar,igrp,size,   &
                     io_pe,num_proc,npq
!
      common/sub_proc/ io_pe,num_proc
! 
      real(C_DOUBLE) pi,dt,rbmax,vth,vth1,vth2,vth3,tmax, &
                     xmax,ymax,zmax,t8,E_Coulomb_P3M,mue, &
                     e_c_s,e_c_pme,e_c_r,e_LJ,e_fene,     &
                     red,vcn,vss0
      common/time3/  t8
      common/Energy/ e_c_s,e_c_pme,e_c_r,e_LJ,e_fene
!
      real(C_float),dimension(npio) :: &  !! Plot: charged only
                    x4,y4,z4,vx4,vy4,vz4,ch4,am4,ag4,ep4
!
      character(len=8) label,cdate*10,ctime*8
      common/headr1/  label,cdate,ctime
      real(C_float)  t,xleng
      common/headr2/ t,xleng
!
      real(C_float),dimension(nhist) :: &
                    ekin,ekn2,ekn3,ecr,elj,ep3m,pot,time,          &
                    vxio,vxco,vxct,vxmc,vxtot,xmcc,vxwat,etot,     &
                    xmac,ymac,zmac,xcon,ycon,zcon,xcoi,ycoi,zcoi,  &
                    xwat,ywat,zwat,dmcc
      common/ehist/ ekin,ekn2,ekn3,ecr,elj,ep3m, pot,time,         &
                    vxio,vxco,vxct,vxmc,vxtot,xmcc,vxwat,etot,     &
                    xmac,ymac,zmac,xcon,ycon,zcon,xcoi,ycoi,zcoi,  &
                    xwat,ywat,zwat,dmcc
!
      integer(C_int) it,is,iwa,iwb,iwc,iwrt1,iwrt2,iwrt3, &
                     nsg,nseg,istop,i,k
      real(C_float)  phi,tht,dtwr,dtwr2,cptot
!
      common/parm1/ it,is
      common/parm2/ pi,dt,rbmax,vth,vth1,vth2,vth3,tmax
      common/parm3/ xmax,ymax,zmax,phi,tht,dtwr,dtwr2
      common/parm9/ cptot
      common/imemo/ iwa,iwb,iwc
      common/iotim/ iwrt1,iwrt2,iwrt3
      common/psegm/ nsg(30),nseg
      common/abterm/ istop
!
      real(C_DOUBLE) rcutlj2,temperat,Bjerrum,alpha,length, &
                     dmesh,prefactor,pref_eps,econv,    &
                     qfrac,Qcore,Rmac,Wmac,Zcp,Zcn,Edc,awat
!
      common/cutoffLJ/ rcutlj2,temperat
      common/elsta/  Bjerrum
      common/ewald1/ alpha,length,dmesh,prefactor,pref_eps,econv
      common/cntion/ qfrac,Qcore,Rmac,Wmac,Zcp,Zcn
      common/EBfild/ Edc
!
      real(C_DOUBLE) dx,dy,dz,r2,r
      real(C_float) qfrac4,qcore4,rmac4,wmac4,zcp4,zcn4,bjerrum4, &
                    alpha4,Edc4,xmax4,ymax4,zmax4,exc,            &  
                    xdd,vmax2,xa,ya,za,xd,yd,zd,dgaus2,           &
                    dtm,vm,s1,s2,s3,vsq
!
      integer(C_int) i1,j
      character(len=2) spc(npq0)  !! Plot: charged only
      logical      first_23,first_p3m
      data         first_23/.true./,first_p3m/.true./
!
      real(C_DOUBLE) cpu0,cpu1,cpu2,cpu3
      integer(C_int) iwrta,iwrtb,iwrtc,nmac,ncti,ncoi,nwat,ncc,cl_first
      real(C_float)  sx1,sy1,sz1,sx2,sy2,sz2,sx3,sy3,sz3,sx4,sy4,sz4, &
                     svx1,svy1,svz1,svx2,svy2,svz2,svx3,svy3,svz3,    &
                     svx4,svy4,svz4,ss
!
      real(C_DOUBLE) t_unit,a_unit,w_unit,e_unit,  &
                     kjoule,kcal,mol,kbT
      common/units/ t_unit,a_unit,w_unit,e_unit
      common/unit2/ kjoule,kcal,mol,kbT
      logical       if_write06/.true./,if_mue/.true./
!
!--------------------------
!*  Initial condition.
!--------------------------
      nCLp= np + nq
!
      if(kstart.le.0) then
        t8= 0.d0
        t = t8
!
        is= 0
        it= 0
        iwa= -1
        iwb= -1
        iwc= -1
      end if
!
!-------------------------------------------------------
      if(io_pe.eq.1) then
        qfrac4 = qfrac
        Qcore4 = Qcore
        Rmac4  = Rmac
        Wmac4  = Wmac
        Zcp4   = Zcp
        Zcn4   = Zcn
!
        Bjerrum4 = Bjerrum
        alpha4 = alpha
        Edc4   = Edc
        xmax4  = xmax
        ymax4  = ymax
        zmax4  = zmax
!
        OPEN (unit=13,file=praefixc//'.13'//suffix2,           &
                           status='unknown',form='unformatted')
!  
        write(13) np,nq,npio,nseg,qfrac4,Qcore4,Rmac4,Wmac4,   & !*added
                  Zcp4,Zcn4,Bjerrum4,alpha4,Edc4,              &
                  xmax4,ymax4,zmax4
!
        do i= 1,npio  !! Plot: charged only
        ch4(i)= ch(i)
        am4(i)= am(i)
        ag4(i)= ag(i)
        ep4(i)= ep(i)
        end do
!
        write(13) ch4,am4,ag4,ep4
        close(13)
      end if
!-------------------------------------------------------
!
 1000 t8= t8 +dt
      t=  t8
      it= it +1
!
!-------------------------------
!  Define write-out intervals.
!-------------------------------
!
      iwrt1= iwrta(t,dtwr )
      iwrt2= iwrtb(t,dtwr2)
      iwrt3= iwrtc(t,5000.) ! 50 ps
!
      cl_first= 2
      call clocks (cpu1,cl_first)
!
      if(t8.gt.tmax) go to 2000
      if((cpu1/60.).gt.cptot) go to 2000
!
      if(istop.ge.1) then
        write(6,*) 'Abnormal termination (istop=1)... '
        write(6,*) '  ipar,t8=',ipar,t8
        go to 2000
      end if
!
!***********************************************
!*  Thermal bath - neutrals  (NOT for ions).   *
!********************************* 9/07/2001 ***
!-------------------------
!*   Velocity update.
!-------------------------
!
      do i = 1,npqr
      frx(i) = 0.d0
      fry(i) = 0.d0
      frz(i) = 0.d0
!
      fsx(i) = 0.d0
      fsy(i) = 0.d0
      fsz(i) = 0.d0
      end do
!
! Short-range forces
      call realteil_s (xg,yg,zg,ch,am,ag,ep,frx,fry,frz, &
                       ipar,size,np,nq,npqr)
!
        cl_first= 2
        call clocks (cpu2,cl_first)
!  
! Long-range forces
      call p3m_perform (xg,yg,zg,ch,fsx,fsy,fsz,npqr,    &
                        E_Coulomb_P3M,first_p3m)
!
      do i= 1,npqr
      frx(i)= frx(i) +fsx(i)
      fry(i)= fry(i) +fsy(i)
      frz(i)= frz(i) +fsz(i)
      end do
!
!************************************************
!*  Free relaxation in the initial stage.       *
!************************************************
!  Temp below is relative to the initial one.
!     Temp= 1.d0/2.d0**(t/2500.d0)
!
!*   <(1/2)*m v**2>= (3/2) kT, n a**3=1. in /init/
!--------------------------------------------
!     vth  = sqrt(2.d0*kbT/(Wmac*w_unit)) /(a_unit/t_unit)
!     vth1 = sqrt(2.d0*kbT/(am(np+nq1)*w_unit)) /(a_unit/t_unit)
!     vth2 = sqrt(2.d0*kbT/(am(np+nq2)*w_unit)) /(a_unit/t_unit)
!     vth3 = sqrt(2.d0*kbT/(am(np+nq+1)*w_unit))/(a_unit/t_unit)
!--------------------------------------------
!
!     epsilj= epsilj0/Temp
!     gamma= Bjerrum0 /(aLJ*Temp)
!
!     gamma= scale_c /surreps   include in L.300
!     prefactor = gamma 
!
      if(t8.lt.10.d0) then
        Exc= 0.d0
      else
        Exc= Edc  ! Converted Edc= econv*Edc0, by READ_CONF 
      end if      ! econv= t_unit**2*e_unit/(w_unit*a_unit *299.98d0)
!
!-------------------------
!*   Velocity update.
!-------------------------
!  the Langevin thermostat 
!    -(mue0* a^2/t) mue*ag(i)*v -> 10^-1/10^2/100= 1 ps
      mue= (a_unit**2/t_unit/100.d0) /50.d0
!                                     50 ps
      if(io_pe.eq.1 .and. if_mue) then
        if_mue= .false.
        write(11,*) 'mue=',mue
      end if  
!
      do i= 1,npqr
      dtm= dt/am(i)
!
      vx(i)= vx(i) +(frx(i) +ch(i)*exc)*dtm -mue*ag(i)*vx(i)*dtm
      vy(i)= vy(i) +(fry(i)           )*dtm -mue*ag(i)*vy(i)*dtm
      vz(i)= vz(i) +(frz(i)           )*dtm -mue*ag(i)*vz(i)*dtm
!
      xg(i)= xg(i) +dt*vx(i)
      yg(i)= yg(i) +dt*vy(i)
      zg(i)= zg(i) +dt*vz(i)
      end do
!
!*  Fold back positions: (xmin,xmax) = (-L/2, L/2). 
!
      if(it.eq.1) then
        do i= 1,npqr
        xg(i) = xg(i) - anint(xg(i)/xmax)*xmax ! real of dnint()
        yg(i) = yg(i) - anint(yg(i)/ymax)*ymax
        zg(i) = zg(i) - anint(zg(i)/zmax)*zmax
        end do
      end if
!
        cl_first= 2
        call clocks (cpu3,cl_first)
!
      if(.false.) then
!     if(io_pe.eq.1 .and. mod(it,100).eq.1) then
        write(11,800) it,cpu3-cpu1,cpu2-cpu1,cpu3-cpu2,cpu1
  800   format('it=',i6,' cpu=',3f10.5,2x,' cpu1=',f10.5,' sec')
      end if
!
!------------------------------
!*  Diagnosis section.
!------------------------------
! 1. Energy.
! --------------------------------------- on major nodes --------------
      if(iwrt1.eq.0 .and. io_pe.eq.1) then
!
        close(11)
        open (unit=11,file=praefixc//'.06'//suffix2,              & 
              status='unknown',position='append',form='formatted')
!      
        is= is +1
        if(is.ge.nhist) then
          write(6,*) ' Stop: hist: is > 10000...'
          istop= 1
        end if
!
        if(if_write06) then
          write(11,600)
  600     format(/,' time:      E_kin    E_kin(cou+co)  E_kin(wat)  ', &
                   'E_pot       e_c_r       e_LJ        e_p3m      ', &
                   ' e_tot       cpu_123')
          if_write06= .false.
        end if
!
        vm= 0
        s1= 0
        s2= 0
        s3= 0
!
        do i= 1,np
        vsq= vx(i)**2 +vy(i)**2 +vz(i)**2
        s1= s1 +0.5d0*am(i)*vsq
!       vm= max(sqrt(vsq), vm)
        end do
!
        do i= np+1,np+nq
        s2= s2 +0.5d0*am(i)*(vx(i)**2 +vy(i)**2 +vz(i)**2)
        end do
!
        do i= nCLp+1,npqr
        s3= s3 +0.5d0*am(i)*(vx(i)**2 +vy(i)**2 +vz(i)**2)
        end do
!
!
        time(is)= t
!
        ekin(is) = s1 
        ekn2(is) = s2 
        ekn3(is) = s3
!
        ecr (is) = e_c_r
        elj (is) = e_LJ 
        ep3m(is) = E_Coulomb_P3M 
        pot (is) = e_c_r +e_LJ +E_Coulomb_P3M
        etot(is) = s1 +s2 +s3 +e_c_r +e_LJ +E_Coulomb_P3M
!
!
        svx1= 0
        svx2= 0
        svx3= 0
        svx4= 0
        sx1 = 0
        sy1 = 0
        sz1 = 0
        sx2 = 0
        sy2 = 0
        sz2 = 0
        sx3 = 0
        sy3 = 0
        sz3 = 0
        sx4 = 0
        sy4 = 0
        sz4 = 0
        nmac= 0
        ncti= 0
        ncoi= 0
        nwat= 0
!
        do i= 1,np
        svx1= svx1 +am(i)*vx(i)
        sx1 = sx1  +xg(i)
        sy1 = sy1  +yg(i)
        sz1 = sz1  +zg(i)
        nmac= nmac +1
        end do
!
        do i= np+1,np+nq
        if(ch(i).gt.0.) then
          svx2= svx2 +am(i)*vx(i)
          sx2 = sx2  +xg(i)
          sy2 = sy2  +yg(i)
          sz2 = sz2  +zg(i)
          ncti= ncti +1
        else
          svx3= svx3 +am(i)*vx(i)
          sx3 = sx3  +xg(i)
          sy3 = sy3  +yg(i)
          sz3 = sz3  +zg(i)
          ncoi= ncoi +1
        end if
        end do
!
        do i= nCLp+1,npqr
        svx4= svx4 +am(i)*vx(i)
        svy4= svy4 +am(i)*vy(i)
        svz4= svz4 +am(i)*vz(i)
        sx4 = sx4  +xg(i)
        sy4 = sy4  +yg(i)
        sz4 = sz4  +zg(i)
        nwat= nwat +1
        end do
!
        vxio(is)= svx1  !/np
        vxct(is)= svx2  !/ncti
        vxco(is)= svx3  !/ncoi
        vxwat(is)= svx4 !/nwat
!
        xmac(is)=  sx1/np
        ymac(is)=  sy1/np
        zmac(is)=  sz1/np
!
        xcon(is)=  sx2/ncti
        ycon(is)=  sy2/ncti
        zcon(is)=  sz2/ncti
!
        xcoi(is)=  sx3/ncoi
        ycoi(is)=  sy3/ncoi
        zcoi(is)=  sz3/ncoi
!
        xwat(is)=  sx4/nwat
        ywat(is)=  sy4/nwat
        zwat(is)=  sz4/nwat
!
!
        svx1= 0
        svx2= 0
        sx1 = 0
        ncc = 0
!
        ss = 0
        do i= 1,np
        ss = ss +ch(i)
        end do
!
        do i= 1,np
        do j= np+1,np+nq
        dx = xg(j) -xg(i)
        dy = yg(j) -yg(i)
        dz = zg(j) -zg(i)
!
        dx = dx - anint(dx/length)*length
        dy = dy - anint(dy/length)*length
        dz = dz - anint(dz/length)*length
        r2 = dx**2 + dy**2 + dz**2
        r  = sqrt(r2)
!
        if(r.lt.Rmac +8.d0) then
          ss = ss +ch(j)  ! charge disribution
!
          sx1 = sx1 +1
          svx1= svx1 +am(j)*vx(j) !-am(i)*vx(i)
          ncc = ncc +1
          if(ch(j).gt.0.d0) then
            svx2= svx2 +am(j)*vx(j) !-am(i)*vx(i)
          end if
        end if
        end do
        end do
!
        dmcc(is) = ss    ! charge distribution f(r)
        xmcc(is) = ncc   !! particles within Rmac+5 A
        vxtot(is)= svx1  !  counter + co
        vxmc(is) = svx2  !  counter only
!
        write(11,610) time(is),s1,s2,s3,pot(is),          &
                      ecr(is),elj(is),ep3m(is),etot(is),  &
                      cpu3-cpu1,cpu2-cpu1,cpu3-cpu2,cpu1/60.
  610   format('t=',f8.1,1p8d12.4,2x,3d12.4,1x,0pf11.3,' M') 
      end if
!
!---------------------------------------------------------------------
!* The routine /redefP/ redefines the monomer positions by shifting
!  each chain so as to fall in the primary domain.
      if(iwrt1.eq.0 .and. io_pe.eq.1) then
!
        OPEN (unit=13,file=praefixc//'.13'//suffix2,            &
              status='old',position='append',form='unformatted')
!
        do i= 1,npio  !! charged only
        x4(i)= xg(i)
        y4(i)= yg(i)
        z4(i)= zg(i)
        end do
!
        write(13) t,x4,y4,z4
        close(13)
      end if
!
!---------------------------------------------------------------------
!*  write(23)
!   make a plot cimv23.xyz
      if(iwrt1.eq.0 .and. io_pe.eq.1) then
!
        if(first_23) then
          first_23= .false.
          open (unit=23,file=praefixc//suffix2//'.xyz',    &
                          status='replace',form='formatted')
        else  
          open (unit=23,file=praefixc//suffix2//'.xyz',    &
                status='unknown',position='append',form='formatted')
        end if
!
        do i= 1,np
        x4(i) = xg(i) 
        y4(i) = yg(i)
        z4(i) = zg(i)
        spc(i)= 'B ' ! 'Me'
        end do
!
        do i= np+1,nCLp  !np+nq
        if(ch4(i).gt.0.) then
          x4(i) = xg(i) 
          y4(i) = yg(i)
          z4(i) = zg(i)
          spc(i) = 'K '
        else if(ch4(i).lt.0.) then
          x4(i) = xg(i) 
          y4(i) = yg(i)
          z4(i) = zg(i)
          spc(i) = 'Cl'
        end if
        end do
!
        do i= 1,np+nq
        x4(i)= x4(i) -anint(x4(i)/xmax)*xmax
        y4(i)= y4(i) -anint(y4(i)/ymax)*ymax
        z4(i)= z4(i) -anint(z4(i)/zmax)*zmax
        end do
!
        write(23,371) np+nq
  371   format(i5,/)

        do i= 1,nCLp !np+nq
        write(23,372) spc(i),x4(i),y4(i),z4(i)
  372   format(a2,3f12.6)
        end do
!
        close(23)
      end if
!
!************************************************************
!  ppl3d/vdistr plots - Charged only !
!    file.77 cimv23.771
      if(iwrt2.eq.0 .and. io_pe.eq.1) then 
!
        open (unit=77,file=praefixc//'.77'//suffix2//'.ps',      & 
              status='unknown',position='append',form='formatted')
!
        i1= 0
!
        do i= 1,nCLp  !!np+nq, charged only
        i1= i1 +1
!
        x4(i1)= xg(i1) -anint(xg(i1)/xmax)*xmax
        y4(i1)= yg(i1) -anint(yg(i1)/ymax)*ymax
        z4(i1)= zg(i1) -anint(zg(i1)/zmax)*zmax
!
        ch4(i1)= ch(i)
        ag4(i1)= ag(i)
        ep4(i1)= ep(i)
        end do
!
        call ppl3da (x4,y4,z4,ch4,ag4,np,nq) 
!
!
!       do j= 1,nCLp  !np+nq
!       vx4(j)= vx(j)
!       vy4(j)= vy(j)
!       vz4(j)= vz(j)
!       end do
!
!       call vdistr (vx4,vy4,vz4,np,nq) 
!
        close(77)
      end if
!
!*  Restart
!
      if(iwrt2.eq.0 .and. io_pe.eq.1) then 
!
        OPEN (unit=12,file=praefixc//'.12'//suffix2,        &
                                         form='unformatted')
        write(12) it,is,size,np,nq,nr,npqr
        write(12) xg,yg,zg,vx,vy,vz,ch,am,ag,ep
        write(12) pi,dt,rbmax,vth,vth1,vth2,vth3,tmax,t8  !*added
        write(12) t,phi,tht,dtwr,dtwr2
        write(12) ekin,ekn2,ekn3,ecr,elj,ep3m,pot,time,           &
                  vxio,vxco,vxct,vxmc,vxtot,xmcc,vxwat,etot,      &
                  xmac,ymac,zmac,xcon,ycon,zcon,xcoi,ycoi,zcoi,   &
                  xwat,ywat,zwat,dmcc
        write(12) iwa,iwb,iwc
        write(12) nsg,nseg
        write(12) alpha,xmax,ymax,zmax
        write(12) qfrac,Qcore,Rmac,Wmac,Zcp,Zcn
        write(12) Edc
        close(12)
      end if     
! --------------------------------------- on major nodes --------------
      go to 1000
!
 2000 continue
      if(io_pe.eq.1) then
        write(11,*) ' final: t8, tmax=',t8,tmax
        write(11,*) ' final: cpu1/60., cptot=',cpu1/60.,cptot
      end if
!
      return
      end subroutine moldyn
!
!
!-------------------------------------------------------------
      subroutine ppl3da (x0,y0,z0,ch,ag,np,nq) !! charged only
!-------------------------------------------------------------
      use, intrinsic :: iso_c_binding
      implicit none
      include   'param_inv3.h'
!
      real(C_float)  x0(npio),y0(npio),z0(npio),ch(npio),ag(npio), &
                     xg(npio),yg(npio),zg(npio),phi,tht,dtwr,dtwr2
      integer(C_INT) np,nq
!
      real(C_DOUBLE) alpha,length,dmesh,prefactor,pref_eps,econv,  &
                     xmax,ymax,zmax,dx,dy,dz,rr,                   &
                     Bjerrum,qfrac,Qcore,Rmac,Wmac,Zcp,Zcn,Edc
      common/ewald1/ alpha,length,dmesh,prefactor,pref_eps,econv
      common/parm3/  xmax,ymax,zmax,phi,tht,dtwr,dtwr2
      common /elsta/ Bjerrum
      common/cntion/ qfrac,Qcore,Rmac,Wmac,Zcp,Zcn
      common/EBfild/ Edc
!
      character(len=8) label,cdate*10,ctime,cax*1
      common/headr1/   label,cdate,ctime
      real(C_float)  time,xleng
      common/headr2/ time,xleng
!
      integer(C_INT) i,j,k,nph1,nph2,nhi
      real(C_float) fsize,rmax1,xx,yy,ps,dd1,dd2,                     &
                    HH,Qcore4,Zcp4,Zcn4,Edc4,Bjerrum4,Rmac4,Wmac4,    &
                    HL,VD,pi,pha,tha,cph,sph,cth,sth,xp,yp,zp,        &
                    xpp,ypp,zpp,x1,y1,z1,dph,dth,th,ph
!
      call newcolor (0, 1.,0.,0.)
!
      HH= 0.8
      CALL SYMBOL (0.5,17.0,HH,label,0.,8)
      CALL SYMBOL ( 0.5,1.0,HH,cdate, 0.,10)
      CALL SYMBOL (15.9,1.0,HH,'t=  ', 0.,4)
      CALL VALUES (999.,1.0,HH,time,0.,101)
!
      Qcore4= Qcore
      Zcp4  = Zcp
      Zcn4  = Zcn
      Edc4  = Edc
!
      CALL SYMBOL (0.5,16.0,HH,'Qcore=', 0.,6)
      CALL VALUES (3.0,16.0,HH,Qcore4,0.,101)
      CALL SYMBOL (7.5,16.0,HH,'Zcp=', 0.,4)
      CALL VALUES (9.3,16.0,HH,Zcp4,0.,101)
      CALL SYMBOL (13.8,16.0,HH,'Zcn=', 0.,4)
      CALL VALUES (15.8,16.0,HH,Zcn4,0.,101)
      CALL SYMBOL (19.8,16.0,HH,'Edc=', 0.,4)
      CALL VALUES (21.8,16.0,HH,Edc4,0.,101)
!
!     Bjerrum4 = Bjerrum
      Rmac4    = Rmac
      Wmac4    = Wmac
      CALL SYMBOL ( 0.5,15.0,HH,'Rmac=', 0.,5)
      CALL VALUES ( 2.8,15.0,HH,Rmac4,0.,101)
      CALL SYMBOL ( 7.0,15.0,HH,'npq=', 0.,4)
      CALL VALUES ( 9.1,15.0,HH,float(npio),0.,101)
      CALL SYMBOL (13.8,15.0,HH,'Wmac=', 0.,5)
      CALL VALUES (15.8,15.0,HH,Wmac4,0.,101)
!
      fsize= 8.
      HL= 11.
      VD=  6.
!
      pi=  4.0*atan(1.0)
      pha= pi*phi/180.
      tha= pi*tht/180.
!
      cph= cos(pha)
      sph= sin(pha)
      cth= cos(tha)
      sth= sin(tha)
!
      xp= xmax*cph -ymax*sph
      yp= xmax*sph +ymax*cph
      zp= zmax
!
      ypp=  yp
      zpp= -xp*sth +zp*cth
!
      rmax1= sqrt(ypp**2 +zpp**2)
      ps= fsize/rmax1
!
!     fsize= 13.
!     rmax1= 0.
!     do i= 1,np+nq  ! for all atoms
!     xx = - x(i)     ! arrange the x-direction with ppl3da plot
!     yy =   y(i)
!     rmax1= max(rmax1, sqrt(xx**2 +yy**2))
!!    rmax1= 2.11d-2     ! fixed in time
!     end do
!     ps= 0.5*fsize/rmax1
!
!**********************
!*  Draw Axes.        *
!**********************
!
      do i= 1,3
      if(i.eq.1) then
         x1= rmax1
         y1= 0.
         z1= 0.
         cax='X'
      else if(i.eq.2) then
         x1= 0.
         y1= rmax1
         z1= 0.
         cax='Y'
      else if(i.eq.3) then
         x1= 0.
         y1= 0.
         z1= rmax1
         cax='Z'
      end if
!
      xp= x1*cph -y1*sph
      yp= x1*sph +y1*cph
      zp= z1
!
      ypp= yp
      zpp= -xp*sth +zp*cth
!
      xx= ps*ypp  +HL
      yy= ps*zpp  +VD
!      
      CALL PLOT (HL,VD,3)
      CALL PLOT (xx,yy,2)
!
      CALL SYMBOL (xx-0.7,yy-0.5,HH,cax,0.,1)
      end do
!
!********************************
!*  Draw a sphere (Macroion).   *
!********************************
!* Draw the Big Globe.
!
      do i= 1,np
      xp= x0(i)*cph -y0(i)*sph
      yp= x0(i)*sph +y0(i)*cph
      zp= z0(i)
!
      xg(i)=  xp*cth +zp*sth
      yg(i)=  yp
      zg(i)= -xp*sth +zp*cth
      end do
!
!
      dph= 2.*pi/360.
      dth=    pi/180.
!
      nph1= -30
      nph2= nph1 +361
!
      do 300 i= 1,np
!* 1. Latitude circles.
!*
      do k= 30,181,30
      nhi= 0
      th= dth * k
!
      do j= nph1,nph2,5
      ph= dph*(j-1)
!
      x1= x0(i) +ag(i)*sin(th)*cos(ph)
      y1= y0(i) +ag(i)*sin(th)*sin(ph)
      z1= z0(i) +ag(i)*cos(th)
!
      xp= x1*cph -y1*sph
      yp= x1*sph +y1*cph
      zp= z1
!
      xpp=  xp*cth +zp*sth
      ypp=  yp
      zpp= -xp*sth +zp*cth
!
!     if(xpp.lt.xg(i)) then  !! if not visible
!       CALL PLOT (xx,yy,3)
!       nhi= nhi +1
!       go to 330
!     end if
!
      xx= ps*ypp +HL
      yy= ps*zpp +VD
!
      if(nhi.eq.0) CALL PLOT (xx,yy,3)
      nhi= nhi +1
      CALL PLOT (xx,yy,2)
      end do
      end do
  330 continue
!
!* 2. Longitude circles.
!
      do j= nph1,nph2,30
      nhi= 0
      ph= dph*(j-1)
!
      do k= 1,181,3
      th= dth * k
!
      x1= x0(i) +ag(i)*sin(th)*cos(ph)
      y1= y0(i) +ag(i)*sin(th)*sin(ph)
      z1= z0(i) +ag(i)*cos(th)
!
      xp= x1*cph -y1*sph
      yp= x1*sph +y1*cph
      zp= z1
!
      xpp=  xp*cth +zp*sth
      ypp= yp
      zpp= -xp*sth +zp*cth
!
!     if(xpp.lt.xg(i)) then
!       CALL PLOT (xx,yy,3)
!       nhi= 0
!       go to 360
!     end if
!
      xx= ps*ypp +HL
      yy= ps*zpp +VD
!
      if(nhi.eq.0) CALL PLOT (xx,yy,3)
      nhi= nhi +1
      CALL PLOT (xx,yy,2)
      end do
      end do
  360 continue
!*
  300 continue
!
!     CALL SYMBOL (15.9,0.1,HH,'RMAX=', 0.,5)
!     CALL VALUES (999.,0.1,HH,rmax1,0.,101)
!
!
!* Draw counterions and coions
      do 400 i= np+1,np+nq
      dx= x0(i) -x0(1)
      dy= y0(i) -y0(1)
      dz= z0(i) -z0(1)
!
      dx = dx - anint(dx/length)*length
      dy = dy - anint(dy/length)*length
      dz = dz - anint(dz/length)*length
      rr = sqrt(dx**2 + dy**2 + dz**2)
      if(rr.gt.21.) go to 400
!
      xp= x0(i)*cph -y0(i)*sph
      yp= x0(i)*sph +y0(i)*cph
      zp= z0(i)
!
      xpp=  xp*cth +zp*sth
      ypp= yp
      zpp= -xp*sth +zp*cth
!
!     do j= 1,np
!     rr= sqrt((ypp -yg(j))**2 +(zpp -zg(j))**2) 
!
!     if(.false.) then  !! now transparent 
!!    if(rr.lt.ag(j)) then
!       if(xpp.lt.xg(j)) go to 400    ! Hide if behind the globe.
!     end if
!     end do
!
      xx= ps*ypp +HL
      yy= ps*zpp +VD
!
      if(ch(i).gt.0.) then
        dd1= 3*ps*ag(i) ! 7*
        call symbol (xx-0.15,yy-0.30,0.7,'+',0.,1)
!
        call newcolor (3, 0.,0.,1.)  ! blue
        call circle (xx-0.12,yy-0.12,dd1,2)
!
      else if(ch(i).lt.0.) then
        dd2= 3*ps*ag(i) ! 7*
        call symbol (xx-0.15,yy-0.30,0.7,'*',0.,1)
!
        call newcolor (3, 0.,1.,0.)  ! green
        call circle (xx-0.12,yy-0.30,dd2,1)
      end if
  400 continue
!
!---------------------
      CALL CHART
!---------------------
!                  grey black.
      call newcolor (0, 1.,0.,0.)
!
      return
      end subroutine ppl3da
!
!
!-------------------------------------------------
       subroutine circle (x,y,d,ic)
!-------------------------------------------------
!*  Open circle centered at (x,y) /or outer edge.
!
      use, intrinsic :: iso_c_binding
!
      write(77,*) " 3.0 setlinewidth"
!
      pi= 3.1415927
      nc= 13
      dth= 2.*pi/nc
      a= d/2.
!
      x0= x +a
      y0= y
      call plot (x0,y0,3)
!
      do j= 1,nc
      th= dth*j
!
      x1= x +a*cos(th)
      y1= y +a*sin(th)
!
      call plot (x1,y1,2)
      end do
!
      call plot (x1,y1,3)
      write(77,*) " 1.0 setlinewidth"
!
      if(ic.eq.1) return
!------------------------------------
!*  Filled circle centered at (x,y).
!------------------------------------
!
      write(77,*) " 3.0 setlinewidth"
!
      nc= 5
      dth= pi/(2*nc +1)
!
      do j= -nc,nc
      th= 0.5*pi +dth*j
!
      x1= x +a*cos(th)
      y1= y +a*sin(th)
!
      x2= x1
      y2= 2.*y -y1
!
      call plot (x1,y1,3)
      call plot (x2,y2,2)
      end do
!
      call plot (x2,y2,3)
      write(77,*) " 1.0 setlinewidth"
!
      return
      end subroutine circle
!
!
!---------------------------------------------------------------------
      subroutine realteil_s (xg,yg,zg,ch,am,ag,ep,frx,fry,frz,  &
                             ipar,size,np,nq,npqr)
!---------------------------------------------------------------------
      use, intrinsic :: iso_c_binding
      implicit none
!*
      include    'param_inv3.h'
      include    'mpif.h'
!
      real(C_DOUBLE),dimension(npqr0) :: &
                    xg,yg,zg,ch,am,ag,ep,frx,fry,frz,ffx,ffy,ffz
      integer(C_int) ipar,size,np,nq,nCLp,npqr
!
      real(C_DOUBLE) t_unit,a_unit,w_unit,e_unit,  &
                     kjoule,kcal,mol,kbT
      common/units/ t_unit,a_unit,w_unit,e_unit
      common/unit2/ kjoule,kcal,mol,kbT
!
      integer(C_int) i,j,ierror
      real(C_DOUBLE) e_cl1,e_lj1,unif1(2),unif2(2),          &
                    dx,dy,dz,r2,r,rlj,asc2,rsi,snt,ccel,     &
                    addpot,forceV,erfc,ar,t,                 &
                    sqrtpi,driwu2,driwu,                     &
                    r_cut,rlj_m,rlj0,pref3,alpha,length,     &
                    dmesh,prefactor,pref_eps,econv,          &
                    pi,dt,rbmax,vth,vth1,vth2,vth3,tmax,     &
                    xmax,ymax,zmax, &
                    rcutlj2,rcutlj,temperat,rcutpme2,rcps2,  &
                    rcutpme,e_c_s,e_c_pme,e_c_r,e_LJ,e_fene, &
                    qfrac,Qcore,Rmac,Wmac,Zcp,Zcn,rr,epsav
      real(C_float) phi,tht,dtwr,dtwr2
!-----------
      common/ewald1/ alpha,length,dmesh,prefactor,pref_eps,econv
      common/parm2/  pi,dt,rbmax,vth,vth1,vth2,vth3,tmax
      common/parm3/  xmax,ymax,zmax,phi,tht,dtwr,dtwr2
      common/cntion/ qfrac,Qcore,Rmac,Wmac,Zcp,Zcn
!
      common /cutoffLJ/ rcutlj2,temperat
      common /cutoffel/ rcutpme2,rcps2
      common /Energy/ e_c_s,e_c_pme,e_c_r,e_LJ,e_fene
!
      real(C_DOUBLE) acount,acoion,epsil_p,epsil_n,epsil_j,epsil_w
      common/ionsiz/ acount,acoion,epsil_p,epsil_n,epsil_j,epsil_w   
      character*2    spc(npq0)  !! Plot: charged only
!
      integer(C_int) it,is,io_pe,num_proc
      common/parm1/ it,is
      common/sub_proc/ io_pe,num_proc
      real(C_float) time,xleng
      common/headr2/ time,xleng
!
!*---------------------------------------------------------------
      pi = 4.d0*atan(1.d0)  !<- generic name
      nCLp= np +nq
!
      driwu2 = 1.25992104989487316476721060728d0  ! 2**(1/3)
      driwu  = sqrt(driwu2)                      ! 2**(1/6)
      sqrtpi = sqrt(pi)
!
      rcutlj =  sqrt(rcutlj2)
      rcutpme = sqrt(rcutpme2)
!
      e_cl1 = 0
      e_lj1 = 0
!
      rlj_m = 0.6d0    ! r^12= 1/(0.6)^12= 459 >> r^2, x0.7d0
      addpot =  1.d0/rcutlj2**6 - 1.d0/rcutlj2**3
!
      if(it.eq.1 .and. io_pe.eq.1) then
        write(11,*) '# rlj_m= ',rlj_m,' for LJ stability rlj !!'
        write(11,*)
      end if
!
!$OMP PARALLEL DEFAULT(NONE)                        &
!$OMP SHARED(ipar,np,nq,nCLp,npqr,size,xg,yg,zg,    &
!$OMP        ch,ep,length,ag,alpha,prefactor,       &
!$OMP        pref_eps,rcutlj,addpot,sqrtpi,         &
!$OMP        rcutpme,rlj_m)                         &
!$OMP PRIVATE(i,j,dx,dy,dz,r2,r,rr,ccel,e_lj,       &
!$OMP        epsav,rlj,rsi,snt,ar,t,erfc,forceV,    &
!$OMP        r_cut,rlj0,pref3)                      &
!$OMP REDUCTION(+:frx,fry,frz,e_cl1,e_lj1)
!$OMP DO SCHEDULE(STATIC,1)
!
      do i= ipar,npqr,size
      do j= 1,npqr
      if(j.eq.i) go to 300
!
      dx = xg(i) -xg(j)
      dy = yg(i) -yg(j)
      dz = zg(i) -zg(j)
!
      dx = dx - anint(dx/length)*length  !<- real of dnint
      dy = dy - anint(dy/length)*length
      dz = dz - anint(dz/length)*length
      r2 = dx**2 + dy**2 + dz**2
      r  = sqrt(r2)
!
!!    scale_c  = (t_unit*e_unit)**2 /(w_unit*a_unit**3)
!!    pref_eps = t_unit**2/(w_unit*a_unit**2)
!          
!*     ^  dv    (t*e)^2 qq' R   t^2  48*eps R   r0       1 r0        *
!*     m ---- = -------*--- - + ----*------ - [(--)^12 - -(--)^6]    * 
!*        dt     ma^3   r^2 r   ma^2   r    r   r        2 r         *
!*                                           snt=1/rlj^6, rlj=r/(ag_i+ag_j)
!*                t^2*e  Edc0(V/cm)                                  *
!*              + ------ ----------  -(mue0*a^2/t) mue*ag(i)*v       *
!*                  ma     300         Langevin thermostat           *
!
!*  Coulomb force: 
!     prefactor - used in realteil, p3m_perform; defined in L.230
!
      r_cut= ag(i)+ag(j) + rcutpme 
      if(r.lt.r_cut) then
!
!  if (i or j) > np+nq
        if(i.gt.nCLp .or. j.gt.nCLp) then
          forceV= 0
!
        else
          ar   = alpha*r  ! <- 1/length
          t    = 1.d0 / ( 1.d0 + PP * ar )
          erfc = t*(A1+t*(A2+t*(A3+t*(A4+t*A5))))
!
          forceV = prefactor*ch(i)*ch(j)*                          &
                       (erfc/r +2*alpha/sqrtpi)*exp(-ar**2)/(r*r)  ! r*r2
          e_cl1 = e_cl1 +prefactor*ch(i)*ch(j)*erfc*exp(-ar**2)/r
        end if
!
!  To prevent close particles from touching
        rlj= r/(ag(i)+ag(j))        ! Non-dimension: hard contact
!
        if(rlj.gt.rcutlj) then      ! rcutlj= driwu
          ccel= 0
!
        else
          rlj0= max(rlj,rlj_m)    ! rlj > rlj_m= 0.6 <- 459^1/12
!
          rsi = 1.d0/rlj0**2        ! 
          snt = rsi*rsi*rsi
!
          epsav = sqrt(ep(i)*ep(j))
          ccel  = 48*(epsav*pref_eps/r) *snt*(snt -0.5d0) /r   ! r*r
          e_lj1 = e_lj1 + 4*epsav*pref_eps *(snt*(snt -1.d0) -addpot)
        end if
!
        frx(i) = frx(i) +(forceV +ccel)*dx
        fry(i) = fry(i) +(forceV +ccel)*dy
        frz(i) = frz(i) +(forceV +ccel)*dz
      end if
!*
  300 end do
      end do
!$OMP END DO
!$OMP END PARALLEL
!
! --------------------
!*   Final results stored in fsx().
! ------------------------------------
      call mpi_allreduce (frx,ffx,npqr,mpi_real8,mpi_sum, & 
                          mpi_comm_world,ierror)
      call mpi_allreduce (fry,ffy,npqr,mpi_real8,mpi_sum, &
                          mpi_comm_world,ierror)
      call mpi_allreduce (frz,ffz,npqr,mpi_real8,mpi_sum, &
                          mpi_comm_world,ierror)
!
      do i= 1,npqr
      frx(i)= ffx(i)
      fry(i)= ffy(i)
      frz(i)= ffz(i)
      end do
!
      unif2(1) = e_cl1
      unif2(2) = e_lj1
      call mpi_allreduce (unif2,unif1,2,mpi_real8,mpi_sum, &
                          mpi_comm_world,ierror)
      e_c_r  = unif1(1)
      e_lj   = unif1(2)
!***
      return
      end subroutine realteil_s
!
!
!CCCCCCCCCCCCCC   P3M (FORTRAN 77)  CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!                                 26.10.1999 
!                                     Motohiko Tanaka, Christian Holm
!  Calling sequences:
!
!     call  p3m_init    (length,npartcls,alpha0Bjerrum)
!     call  p3m_perform (coox,cooy,cooz,Q,FX,FY,FZ,E_Coulomb_P3M)
!
!/*---------------------------------------------------------------------
! SUBUNIT:  p3m_v2.f  (FORTRAN 90)
! 
!       Version   20.10.1999 (CH) 
!       corrected 26.10.1999 (MT)
!                 23.06.2000 (MT)
!
! VERSION:  20 January 1999, written in C
! AUTHOR:   Markus Deserno
!
!    Brillouin is now a parameter
!    MaxInterpol --> MINTPOL
!    floor --> Defined: must be stepwise at x= 0.
!    dround --> = nint (Round off to nearest integer). 
!/*---------------------------------------------------------------------
!
! void p3m_perform(double coo[][3], double *Q, double For[][3], 
! double *E_Coulomb_P3M)
!/*---------------------------------------------------------------------
      subroutine p3m_perform (xg,yg,zg,ch,fsx,fsy,fsz,npqr,  &
                              E_Coulomb_P3M,first_p3m)
!/*---------------------------------------------------------------------
      use, intrinsic :: iso_c_binding
      use omp_lib
      implicit none
!
      include    'param_inv3.h'
!     include    'aslfftw3.f03'                 ! SX
!     include    'fftw3.f03'                    ! LX
      include    '/opt/fftw3/include/fftw3.f03' ! physique
!
      real(C_DOUBLE),dimension(0:npqr0-1) :: xg,yg,zg,ch,fsx,fsy,fsz
      real(C_DOUBLE),dimension(0:npqr0-1,0:2) :: fek
!
!     integer n_thread,ddd
      type(C_PTR),save :: plan, pinv 
      integer(C_INT)   :: npqr,ierror,io_pe,num_proc
      common/sub_proc/ io_pe,num_proc
!
      real(C_DOUBLE),dimension(0:mesh-1,0:mesh-1,0:mesh-1)         &
                                  :: qq,phi_x,phi_y,phi_z
      complex(C_DOUBLE_COMPLEX),                                   &
                     dimension(0:mesh/2,0:mesh-1,0:mesh-1)         &
                                  :: qq_c,phi_x_c,phi_y_c,phi_z_c 
!-----------
      real(C_DOUBLE) alpha,length,dmesh,prefactor,pref_eps,econv,  &
                     E_Coulomb_P3M,ecl,wupi,pi
      common/ewald1/ alpha,length,dmesh,prefactor,pref_eps,econv
!
      real(C_DOUBLE) fft_scale2
      complex(kind(0d0)) ei,eigq  ! complex*16
!
      integer(C_int) xarg(0:npq0-1),yarg(0:npq0-1),   &
                  zarg(0:npq0-1),i,j,k,m,MESHMASK, &
                  Gi0,Gi1,Gi2, xpos,ypos,zpos,     &
                  assignshift, QZahl,m0
      real(C_DOUBLE) d1,Hi,MI2,modadd1,modadd2, T1,T2,T3, &
                  sum_q_2, sum_q2
!-----------
      real(C_DOUBLE) coop,qp,ql,meshift,Dn,Ghat,intCAF
      integer(C_int) global,G
!-----------
      common/coordi/ coop(0:npq0-1,0:2),qp(0:npq0-1)
      common/pindex/ global(0:npq0-1),G(0:npq0-1,0:2)
      common/prmesh/ ql(0:P_max**3-1,0:npq0-1)
      common/mesh01/ meshift(0:mesh-1),Dn(0:mesh-1)
!
      common/influf/ Ghat(0:mesh/2,0:mesh-1,0:mesh-1)
      common/intcaf/ intCAF(0:P_max,0:2*MINTPOL+1)
!
!-----------
      integer(C_int) it,is
      common/parm1/  it,is
!
      integer(C_int) For_flag, E_Coulomb_P3M_flag
      common/flags1/ For_flag, E_Coulomb_P3M_flag
!
      integer(C_int) iwrt1,iwrt2,iwrt3
      common/iotim/  iwrt1,iwrt2,iwrt3
!
      logical :: first_p3m
!
! ---------------------------
!*  Prepare for FFTW calls
! ---------------------------
      if(first_p3m) then
        first_p3m= .false.
!
        ierror= 0
!       n_thread= 1 ! 8 
!       ddd= fftw_init_threads (ierror)
!       call fftw_plan_with_nthreads (n_thread)
!       call fftw_plan_with_nthreads (omp_get_max_threads())
!       plan= fftw_plan_dft_r2c_3d & 
!                 (mesh,mesh,mesh,qq,qq_c,FFTW_ESTIMATE)
!                  N, M, L
!       pinv= fftw_plan_dft_c2r_3d &
!                 (mesh,mesh,mesh,qq_c,qq,FFTW_ESTIMATE)
!                  N/2+1, M, L
!
!  SX case only
!       call fftw_plan_with_nthreads (omp_get_max_threads()) 
!
!       call dfftw_plan_dft_r2c_3d  &
        plan= fftw_plan_dft_r2c_3d  &
                (mesh,mesh,mesh,qq,qq_c,FFTW_ESTIMATE)
!                   n0,m0,l0 ---------> (n0/2+1) complex
!
!       call dfftw_plan_dft_c2r_3d  &
        pinv= fftw_plan_dft_c2r_3d  &
                (mesh,mesh,mesh,phi_x_c,phi_x,FFTW_ESTIMATE)
      end if
!
      ei = cmplx(0.d0,1.d0,kind(0d0))
! -----------------------------------
      pi = 4.d0*atan(1.d0)
      MESHMASK = mesh-1     
!
      dmesh = dfloat(mesh)
      Hi = dmesh / length
      MI2 = 2.d0*dfloat(MINTPOL)
      assignshift = mesh -(IP-1)/2
!
      qzahl = 0
      sum_q_2 = 0
      sum_q2  = 0
!
      do i= 0,npq0-1
      if (abs(ch(i)) .gt. 1.d-5) then 
        coop(qzahl, 0) = xg(i) - anint(xg(i)/length -0.5d0)*length ! real of dnint
        coop(qzahl, 1) = yg(i) - anint(yg(i)/length -0.5d0)*length
        coop(qzahl, 2) = zg(i) - anint(zg(i)/length -0.5d0)*length
!
        qzahl= qzahl + 1
        qp(qzahl) = ch(i)
        global(qzahl) = i
!
        sum_q_2 = sum_q_2 + qp(qzahl)
        sum_q2  = sum_q2  + qp(qzahl)**2
      end if
      end do
      sum_q_2 = sum_q_2 **2
!
!
      do i= 0,mesh-1 
      do j= 0,mesh-1 
      do k= 0,mesh-1 
      qq(i,j,k) = 0.d0
      end do
      end do
      end do
!
      if(IP.eq.2 .or. IP.eq.4 .or. IP.eq.6) then
         modadd1 =  0.5d0
         modadd2 = -0.5d0
      else 
         if(IP.eq.1 .or. IP.eq.3 .or. IP.eq.5 .or. IP.eq.7) then
            modadd1 = 0.0d0
            modadd2 = 0.5d0
         else
            write(06,*) "Error in function 'P3M_perform':"
            write(06,*) "Charge assignment order IP=",IP," unknown."
            write(06,*) "Program terminated."
            call exit(1)
         end if
      end if
!
!***
!
      do i= 0, qzahl-1
      d1  = coop(i,0)*Hi + modadd1
      Gi0 = int(d1 + modadd2) + assignshift
      G(i,0) = Gi0 
      xarg(i) = int( (d1 - anint(d1) + 0.5d0)*MI2 ) !! real of dnint()
!      
      d1  = coop(i,1)*Hi + modadd1 
      Gi1 = int(d1 + modadd2) + assignshift
      G(i,1) = Gi1
      yarg(i) = int( (d1 - anint(d1) + 0.5d0)*MI2 )
!      
      d1  = coop(i,2)*Hi + modadd1 
      Gi2 = int(d1 + modadd2) + assignshift
      G(i,2) = Gi2 
      zarg(i) = int( (d1 - anint(d1) + 0.5d0)*MI2 )
!
      m0= -1
!*
      do j = 0, IP-1
      do k = 0, IP-1
      do m = 0, IP-1
      t1 = qp(i) * intCAF(j,xarg(i))
      t2 = t1 *    intCAF(k,yarg(i))
      t3 = t2 *    intCAF(m,zarg(i))
!    
      m0= m0 + 1
      ql(m0,i) = t3          ! assignment factor.
!      
      xpos = IAND( (G(i,0) + j), MESHMASK)
      ypos = IAND( (G(i,1) + k), MESHMASK)
      zpos = IAND( (G(i,2) + m), MESHMASK)
!    
      qq(xpos,ypos,zpos) = qq(xpos,ypos,zpos) + ql(m0,i)
      end do
      end do
      end do
      end do
!
!* Scale factor for paired (-1,+1) FFTW execution.
!
      fft_scale2= mesh**3  
!     call dfftw_execute_dft_r2c (plan,qq,qq_c)
      call fftw_execute_dft_r2c (plan,qq,qq_c)
!
!
      do k= 0,mesh-1
      do j= 0,mesh-1
      do i= 0,mesh/2
      if(i.eq.0) then
        phi_x_c(0,j,k) = 0.d0  !! <-- zeros at i= 0
        phi_y_c(0,j,k) = 0.d0
        phi_z_c(0,j,k) = 0.d0
      else
      eigq = ei*Ghat(i,j,k)*qq_c(i,j,k)/fft_scale2
!
      phi_x_c(i,j,k) = eigq*Dn(i)
      phi_y_c(i,j,k) = eigq*Dn(j)
      phi_z_c(i,j,k) = eigq*Dn(k)
      end if
!*
      end do 
      end do 
      end do 
!
!*  For diagnosis.
!    
      if( For_flag.ne.0 .and. (iwrt1.eq.0)) then
         do i= 0,mesh/2
         do j= 0,mesh-1
         do k= 0,mesh-1
         qq_c(i,j,k) = Ghat(i,j,k)*qq_c(i,j,k)/fft_scale2
         end do
         end do
         end do
!
         call fftw_execute_dft_c2r (pinv,qq_c,qq)
      end if
!
!     call dfftw_execute_dft_c2r (pinv,phi_x_c,phi_x)
      call fftw_execute_dft_c2r (pinv,phi_x_c,phi_x)
      call fftw_execute_dft_c2r (pinv,phi_y_c,phi_y)
      call fftw_execute_dft_c2r (pinv,phi_z_c,phi_z)
! 
!
      do i = 0, qzahl-1
      m0= -1
!
      do j = 0, IP-1
      do k = 0, IP-1
      do m = 0, IP-1
      xpos = IAND( (G(i,0) + j), MESHMASK)
      ypos = IAND( (G(i,1) + k), MESHMASK)
      zpos = IAND( (G(i,2) + m), MESHMASK)
!    
      m0 = m0 +1
      d1 = prefactor * QL(m0,i)
!          +++++++++ 
      fsx(global(i)) = fsx(global(i)) - d1*phi_x(xpos,ypos,zpos)
      fsy(global(i)) = fsy(global(i)) - d1*phi_y(xpos,ypos,zpos)
      fsz(global(i)) = fsz(global(i)) - d1*phi_z(xpos,ypos,zpos)
      end do
      end do
      end do
      end do
!
!*
      if (E_Coulomb_P3M_flag.ne.0 .and. (iwrt1.eq.0)) then
        ecl = 0.d0
!
        do i= 0, mesh/2
        do j= 0, mesh-1
        do k= 0, mesh-1
        ecl = ecl + Ghat(i,j,k)*abs(qq_c(i,j,k))**2/fft_scale2 
        end do
        end do
        end do
!
        wupi = sqrt(pi)
        E_Coulomb_P3M = prefactor * ecl* (1/dmesh) * length/(4.d0*pi) &
                       - length *prefactor *( sum_q2 *alpha /wupi     &
                             + sum_q_2 *pi /(2.d0*length**2*alpha**2) )
!       write(6,*) 'E_Coulomb (all)  =',E_Coulomb_P3M
      end if
!
      return
      end subroutine p3m_perform
!
!
!/*---------------------------------------------------------------------
      subroutine p3m_init
!/*---------------------------------------------------------------------
      use, intrinsic :: iso_c_binding
      implicit none
!
      include  'param_inv3.h'
!
      real(C_DOUBLE) t_unit,a_unit,w_unit,e_unit,  &
                     kjoule,kcal,mol,kbT,Bjerrum
      common/units/ t_unit,a_unit,w_unit,e_unit
      common/unit2/ kjoule,kcal,mol,kbT
      common /elsta/ Bjerrum
!
!     integer*4   P_max,MINTPOL --> param_inv3.h
      integer(C_int) io_pe,num_proc
      common/sub_proc/ io_pe,num_proc
!
      real(C_DOUBLE) alpha,length,dmesh,prefactor,pref_eps,econv
      common/ewald1/ alpha,length,dmesh,prefactor,pref_eps,econv
!
      real(C_DOUBLE) pi,dt,rbmax,vth,vth1,vth2,vth3,tmax,  &
                     rcutlj2,temperat
      common/parm2/  pi,dt,rbmax,vth,vth1,vth2,vth3,tmax
      common /cutoffLJ/ rcutlj2,temperat
!
      integer(C_int) For_flag, E_Coulomb_P3M_flag
      common/flags1/ For_flag, E_Coulomb_P3M_flag
!
!
      pi = 4.d0*atan(1.d0)
!
!  Use flags instead of 0-addresses.
!  *********************************
      For_flag      = 1
      E_Coulomb_P3M_flag = 1
!  *********************************
      dmesh = dfloat(mesh)
!
!  P3M is initialized
      call  interpol_charge_assign_function
      call  calculate_meshift
      call  calculate_differential_operator
      call  calculate_influence_function
!
      if(io_pe.eq.1) then
        write(11,*) 'P3M is successfully initialized !!'
        write(11,*)
      end if
!
      return
      end subroutine P3M_init
!
!
!/*---------------------------------------------------------------------
      subroutine perform_aliasing_sums (nx,ny,nz,nominatorX,nominatorY, &
                                        nominatorZ,denominator)
!/*---------------------------------------------------------------------
      use, intrinsic :: iso_c_binding
      implicit none
!
      include  'param_inv3.h'
!
!  Brillouin --> param_inv3.h
      integer(C_int) io_pe,num_proc
      common/sub_proc/ io_pe,num_proc
!
      real(C_DOUBLE) alpha,length,dmesh,prefactor,pref_eps,econv, & 
                     meshift,Dn,pi
      common/ewald1/ alpha,length,dmesh,prefactor,pref_eps,econv
      common/mesh01/ meshift(0:mesh-1),Dn(0:mesh-1)
!
      integer(C_int) nx,ny,nz,nnx,nny,nnz
      real(C_DOUBLE) nominatorX,nominatorY,nominatorZ,denominator, &
                  S1,S2,S3,fak1,fak2,fak3,NMX,NMY,NMZ,NM2,expo,    &
                  exponent_limit,sinc
!
      exponent_limit = 30.d0
!
      pi = 4.d0*atan(1.d0)
      fak1 = 1.d0/dmesh
      fak2 = (pi/(alpha*length))**2  ! <- READ_-CONF
!
      nominatorX = 0.d0
      nominatorY = 0.d0
      nominatorZ = 0.d0
      denominator = 0.d0
!
      do nnx = -Brillouin, Brillouin
      NMX = meshift(nx) + dmesh*nnx
      S1  = sinc(fak1*NMX)**(2*IP)
!
        do nny = -Brillouin, Brillouin
        NMY = meshift(ny) + dmesh*nny
        S2  = S1* sinc(fak1*NMY)**(2*IP)
!
          do nnz = -Brillouin, Brillouin
          NMZ = meshift(nz) + dmesh*nnz
          S3  = S2* sinc(fak1*NMZ)**(2*IP)
!
          denominator = denominator + S3
          NM2 = NMX**2 + NMY**2 + NMZ**2
!
          expo= fak2*NM2
          if (expo .lt. exponent_limit) then
              fak3 =  S3* dexp(-expo)/NM2 
          else
              fak3 = 0.d0
          end if
!
          nominatorX = nominatorX + fak3 * NMX
          nominatorY = nominatorY + fak3 * NMY
          nominatorZ = nominatorZ + fak3 * NMZ
          end do
        end do
      end do

      return
      end subroutine perform_aliasing_sums
!
!
!/*---------------------------------------------------------------------
      subroutine  calculate_differential_operator
!/*---------------------------------------------------------------------
      use, intrinsic :: iso_c_binding
      implicit none
!
      include  'param_inv3.h'
!
!     integer*4   P_max,MINTPOL --> param_inv3.h
      integer(C_int) io_pe,num_proc,i
      common/sub_proc/ io_pe,num_proc
!  
      real(C_DOUBLE) alpha,length,dmesh,prefactor,pref_eps,econv, & 
                     meshift,Dn
      common/ewald1/ alpha,length,dmesh,prefactor,pref_eps,econv
      common/mesh01/ meshift(0:mesh-1),Dn(0:mesh-1)
!                                      ******  defined here.
!
      if(io_pe.eq.1) then
        write(11,*) " - calculating differential operator"
      end if
!
      do i= 0,mesh-1
      Dn(i) =  real(i) - anint(dfloat(i)/dmesh)*dmesh
!     Dn(i) = float(i) - anint(dfloat(i)/dmesh)*dmesh
      end do 
!
      Dn(mesh/2) = 0.d0
!
      return
      end subroutine  calculate_differential_operator
!
!
!/*---------------------------------------------------------------------
      subroutine  calculate_influence_function
!/*---------------------------------------------------------------------
      use, intrinsic :: iso_c_binding
      implicit none
!
      include  'param_inv3.h'
!
!     integer*4   P_max,MINTPOL --> param_inv3.h
      real(C_DOUBLE) alpha,length,dmesh,prefactor,pref_eps,econv, &
                     meshift,Dn,Ghat,pi
      common/ewald1/ alpha,length,dmesh,prefactor,pref_eps,econv
      common/mesh01/ meshift(0:mesh-1),Dn(0:mesh-1)
      common/influf/ Ghat(0:mesh/2,0:mesh-1,0:mesh-1)
!                    **************************  defined here.
!
      integer(C_int) nx,ny,nz
      real(C_DOUBLE) Dnx,Dny,Dnz,Dn2,fak1,fak2,fak3,            &
                     nominatorX,nominatorY,nominatorZ,denominator
      integer(C_int) io_pe,num_proc
      common/sub_proc/ io_pe,num_proc
!
      if(io_pe.eq.1) then
        write(11,*) " - calculating influence function with parameters."
        write(11,601) mesh,IP,alpha,length
  601   format('   mesh=',i5,'  ,IP=',i5,/,     &
               '   alpha=',d18.10,',  length=',d18.10)
      end if
! 
       pi = 4.d0*atan(1.d0)
       fak1  = dmesh*dmesh*dmesh * 2.d0 / length**2
       fak2  = (pi/(alpha*length))**2  ! <--READ_CONF
!
!     do NX = 0, mesh-1
      do nx = 0, mesh/2
       do ny = 0, mesh-1
        do nz = 0, mesh-1
        if ( (nx.eq.0).and.(ny.eq.0).and.(nz.eq.0)) then
           Ghat(nx,ny,nz)= 0.d0 
        else
           call perform_aliasing_sums (nx,ny,nz,nominatorX,nominatorY, &
                                       nominatorZ,denominator)
           Dnx = Dn(nx)
           Dny = Dn(ny)
           Dnz = Dn(nz)
!
           Dn2 = Dnx**2 + Dny**2 + Dnz**2
!  
           if (Dn2 .gt. 1.d-7) then
             fak3 = Dnx*nominatorX + Dny*nominatorY + Dnz*nominatorZ
             fak3 = fak3 /( Dn2 * denominator**2 )
             Ghat(nx,ny,nz) = fak3*fak1
           else 
             Ghat(nx,ny,nz) = 0.d0
           end if
        end if
!
        end do
       end do
      end do
!
      return
      end subroutine  calculate_influence_function
!
!
!/*---------------------------------------------------------------------
      subroutine interpol_charge_assign_function
!/*---------------------------------------------------------------------
      use, intrinsic :: iso_c_binding
      implicit none
!
      include  'param_inv3.h'
!
!     integer*4   P_max,MINTPOL --> param_inv3.h
      integer(C_int) io_pe,num_proc,i
      common/sub_proc/ io_pe,num_proc
!
      real(C_DOUBLE) intCAF, dInterpol, x
      common/intcaf/ intCAF(0:P_max,0:2*MINTPOL+1)
!                    ***************************** defined here.
!
      dInterpol= dfloat(MINTPOL)
!
      if(io_pe.eq.1) then
        write(11,601) IP
  601   format(/,'- interpolating the order-',i2, &
               ' charge assignment function')
      end if
!
      if (IP.eq.1) then
      do 100 i= -MINTPOL, MINTPOL
      x= i/(2.d0*dInterpol)
      intCAF(0, i+MINTPOL) = 1.d0
  100 continue
!
      else if (IP.eq.2) then
      do 200 i= -MINTPOL, MINTPOL
      x= i/(2.d0*dInterpol)
      intCAF(0, i+MINTPOL) = 0.5d0 -x
      intCAF(1, i+MINTPOL) = 0.5d0 +x
  200 continue
!
      else if (IP.eq.3) then
      do 300 i= -MINTPOL, MINTPOL
      x= i/(2.d0*dInterpol)
      intCAF(0, i+MINTPOL) = 0.50d0*(0.5d0 - x)**2
      intCAF(1, i+MINTPOL) = 0.75d0 - x*x
      intCAF(2, i+MINTPOL) = 0.50d0*(0.5d0 + x)**2
  300 continue
!
!     else if (IP.eq.4) then
!     do 400 i= -MINTPOL, MINTPOL
!     x= i/(2.d0*dInterpol)
!     intCAF(0, i+MINTPOL) = 
!    *                       ( 1.d0+x*( -6.d0+x*( 12.d0-x* 8.d0)))/48.d0
!     intCAF(1, i+MINTPOL) = 
!    *                       (23.d0+x*(-30.d0+x*(-12.d0+x*24.d0)))/48.d0
!     intCAF(2, i+MINTPOL) = 
!    *                       (23.d0+x*( 30.d0+x*(-12.d0-x*24.d0)))/48.d0
!     intCAF(3, i+MINTPOL) =  
!    *                       ( 1.d0+x*(  6.d0+x*( 12.d0+x* 8.d0)))/48.d0
! 400 continue
!
!     else if (IP.eq.5) then
!     do 500 i= -MINTPOL, MINTPOL
!     x= i/(2.d0*dInterpol)
!     intCAF(0, i+MINTPOL) = 
!    *         (  1.d0+x*( -8.d0+x*(  24.d0+x*(-32.d0+x*16.d0))))/384.d0
!     intCAF(1, i+MINTPOL) = 
!    *         ( 19.d0+x*(-44.d0+x*(  24.d0+x*( 16.d0-x*16.d0))))/ 96.d0
!     intCAF(2, i+MINTPOL) = 
!    *         (115.d0+x*        x*(-120.d0+x*        x*48.d0))  /192.d0
!     intCAF(3, i+MINTPOL) = 
!    *         ( 19.d0+x*( 44.d0+x*(  24.d0+x*(-16.d0-x*16.d0))))/ 96.d0
!     intCAF(4, i+MINTPOL) = 
!    *         (  1.d0+x*(  8.d0+x*(  24.d0+x*( 32.d0+x*16.d0))))/384.d0
! 500 continue
!
!     else if (IP.eq.6) then
!     do 600 i= -MINTPOL, MINTPOL
!     x= i/(2.d0*dInterpol)
!     intCAF(0, i+MINTPOL) = 
!    *     (  1.d0+x*( -10.d0+x*(  40.d0+x*( -80.d0+x*(  80.d0-x* 32.d0)
!    *                                                      ))))/3840.d0
!     intCAF(1, i+MINTPOL) = 
!    *     (237.d0+x*(-750.d0+x*( 840.d0+x*(-240.d0+x*(-240.d0+x*160.d0)
!    *                                                      ))))/3840.d0
!     intCAF(2, i+MINTPOL) = 
!    *     (841.d0+x*(-770.d0+x*(-440.d0+x*( 560.d0+x*(  80.d0-x*160.d0)
!    *                                                      ))))/1920.d0
!     intCAF(3, i+MINTPOL) = 
!    *     (841.d0+x*(+770.d0+x*(-440.d0+x*(-560.d0+x*(  80.d0+x*160.d0)
!    *                                                      ))))/1920.d0
!     intCAF(4, i+MINTPOL) = 
!    *     (237.d0+x*( 750.d0+x*( 840.d0+x*( 240.d0+x*(-240.d0-x*160.d0)
!    *                                                      ))))/3840.d0
!     intCAF(5, i+MINTPOL) = 
!    *     (  1.d0+x*(  10.d0+x*(  40.d0+x*(  80.d0+x*(  80.d0+x* 32.d0)
!    *                                                      ))))/3840.d0
! 600 continue
!
!     else if (IP.eq.7) then
!     do 700 i= -MINTPOL, MINTPOL
!     x= i/(2.d0*dInterpol)
!     intCAF(0, i+MINTPOL) = 
!    *            (    1.d0+x*(   -12.d0+x*(   60.d0+x*( -160.d0+x*( 
!    *                         240.d0+x*(-192.d0+x* 64.d0))))))/46080.d0
!     intCAF(1, i+MINTPOL) = 
!    *            (  361.d0+x*( -1416.d0+x*( 2220.d0+x*(-1600.d0+x*(
!    *                         240.d0+x*( 384.d0-x*192.d0))))))/23040.d0
!     intCAF(2, i+MINTPOL) = 
!    *            (10543.d0+x*(-17340.d0+x*( 4740.d0+x*( 6880.d0+x*(
!    *                       -4080.d0+x*(-960.d0+x*960.d0))))))/46080.d0
!     intCAF(3, i+MINTPOL) = 
!    *            ( 5887.d0+x*          x*(-4620.d0+x*         x*( 
!    *                         1680.d0-x*        x*320.d0)))   /11520.d0
!     intCAF(4, i+MINTPOL) = 
!    *            (10543.d0+x*( 17340.d0+x*( 4740.d0+x*(-6880.d0+x*(
!    *                       -4080.d0+x*( 960.d0+x*960.d0))))))/46080.d0
!     intCAF(5, i+MINTPOL) = 
!    *            (  361.d0+x*(  1416.d0+x*( 2220.d0+x*( 1600.d0+x*(
!    *                         240.d0+x*(-384.d0-x*192.d0))))))/23040.d0
!     intCAF(6, i+MINTPOL) = 
!    *            (    1.d0+x*(    12.d0+x*(   60.d0+x*(  160.d0+x*(  
!    *                         240.d0+x*( 192.d0+x* 64.d0))))))/46080.d0
! 700 continue
!
!     else if (IP.gt.7 .or. IP.lt.1) then
      else if (IP.gt.3 .or. IP.lt.1) then
        if(io_pe.eq.1) then
          write(11,*) "Error in function ", &
                      "'interpolate_charge_assignment_function':"
          write(11,*) IP
  611     format("Charge assignment order",i2," unknown.",/, &
                 "Program terminated.")
        end if
        call exit(1)
      end if
!
      return
      end subroutine interpol_charge_assign_function
!
!
!/*---------------------------------------------------------------------
      subroutine  calculate_meshift
!/*---------------------------------------------------------------------
      use, intrinsic :: iso_c_binding
      implicit none
!
      include  'param_inv3.h'
!
!     integer*4   P_max,MINTPOL --> param_inv3.h
      integer(C_int) io_pe,num_proc,i
      common/sub_proc/ io_pe,num_proc
!
      real(C_DOUBLE) alpha,length,dmesh,prefactor,pref_eps,econv,  &
                     meshift,Dn
      common/ewald1/ alpha,length,dmesh,prefactor,pref_eps,econv
      common/mesh01/ meshift(0:mesh-1),Dn(0:mesh-1)
!-----------
! 
      if(io_pe.eq.1) then
        write(11,*) " - calculating mesh-shift"
      end if
!  
      do 100 i= 0, mesh-1
      meshift(i) = i - anint(i/dmesh)*dmesh  !! real of dnint() 
  100 continue
!
      return
      end subroutine  calculate_meshift
!
!
!/*---------------------------------------------------------------------
      function sinc (d)
!/*---------------------------------------------------------------------
      use, intrinsic :: iso_c_binding
      implicit none
!
      real(C_DOUBLE) sinc,d,epsi,c2,c4,c6,c8,pi,pid,pid2
!
      epsi =  0.1d0
      c2 = -0.1666666666667d-0
      c4 =  0.8333333333333d-2
      c6 = -0.1984126984127d-3
      c8 =  0.2755731922399d-5
!
      pi = 4.d0*atan(1.d0)
      pid = pi*d
!
      if (abs(d).gt.epsi) then
         sinc = dsin(pid) / pid
      else 
         pid2 = pid*pid
         sinc = 1.d0 + pid2*( c2 + pid2*( c4 + pid2*(c6 + pid2*c8) ) )
      end if
!
      if(abs(sinc).lt.1.d-100) sinc= 0.d0
!
      return
      end function sinc
!
!
!---------------------------------------------------------------------
      function floor (x)
!---------------------------------------------------------------------
      use, intrinsic :: iso_c_binding
      implicit none
!
      real(C_DOUBLE) floor, x, xlim
      integer(C_int) io_pe,num_proc
      common/sub_proc/ io_pe,num_proc
!
       xlim= 100000.d0
       if(abs(x).lt.xlim) then
         floor = int(x + xlim) - xlim
       else
         if(io_pe.eq.1) then
           write(11,*) " FLOOR: Argument too large -- Run terminated."
         end if
         call exit (1)
       end if
!
       return
       end function floor
!
!
!  READ /WRITE configuration data.
!------------------------------------------------------------------
      subroutine READ_CONF (xg,yg,zg,vx,vy,vz,ch,am,ag,np)
!------------------------------------------------------------------
!  Four quantities - a(length), m(mass), t(time) and e(charge)
!    1 Angstrom, proton's mass, 10^-14 s, and electron charge
!
      use, intrinsic :: iso_c_binding
      implicit none
!
      include    'param_inv3.h'
!
      integer(C_int) io_pe,num_proc
      common/sub_proc/ io_pe,num_proc
!
      real(C_DOUBLE) t_unit,a_unit,w_unit,e_unit,  &
                     kjoule,kcal,mol,kbT
      common/units/ t_unit,a_unit,w_unit,e_unit
      common/unit2/ kjoule,kcal,mol,kbT
!
      character  praefix*6,text1*40
      common/CONFdataS/ praefix
      integer(C_int) maxpol
      parameter  (maxpol=128)  !,maxmesh=64)
!
      real(C_DOUBLE),dimension(npqr0) :: &
                     xg,yg,zg,vx,vy,vz,ch,am,ag
!
      real(C_DOUBLE) alpha,length,dmesh,prefactor,pref_eps,econv, &
                     rcutpme,rcutLJ,r_fene,rcutlj2,temperat,rcutpme2, &
                     rcps2,k_fene,r_fene2,Bjerrum
!
      common/ewald1/ alpha,length,dmesh,prefactor,pref_eps,econv
      common/CONFdataR/ rcutpme,rcutLJ,r_fene
      common/cutoffLJ/  rcutlj2,temperat
      common/cutoffel/  rcutpme2,rcps2
      common/fenepara/  k_fene,r_fene2
      common/elsta/ Bjerrum
!
      integer(C_int) N_P,MPC,LADABST,MAXION,PLUSION,      &
                     i_Steps,MeasStep,ConfStep,           &
                     N_LP,v_G,n_SMol,v_SP,v_SM,seed
      real(C_DOUBLE) skin,skin2
      common/POLY/  N_P,MPC,LADABST
      common/ION/ MAXION,PLUSION
      common/confdatai/ i_Steps,MeasStep,ConfStep,   &
                        N_LP,v_G,n_SMol,v_SP,v_SM,seed
      common/einles/ skin,skin2

      integer(C_int) np,nsg,nseg,npartcls,isize,I,Verletfix
!----------------------------------------------------------------
      real(C_DOUBLE) pi,dt,rbmax,vth,vth1,vth2,vth3,tmax,tmin, &
                     xmax,ymax,zmax,qfrac,Qcore,Rmac,Wmac,     &
                     Zcp,Zcn,Edc0,Edc
      real(C_float)  phi,tht,dtwr,dtwr2,cptot
!
      common/parm2/ pi,dt,rbmax,vth,vth1,vth2,vth3,tmax
      common/psegm/ nsg(30),nseg
      common/parm3/ xmax,ymax,zmax,phi,tht,dtwr,dtwr2
      common/parm9/ cptot
      common/cntion/ qfrac,Qcore,Rmac,Wmac,Zcp,Zcn
      common/EBfild/ Edc
!
      real(C_DOUBLE) acount,acoion,epsil_p,epsil_n,epsil_j,epsil_w
      common/ionsiz/ acount,acoion,epsil_p,epsil_n,epsil_j,epsil_w   
!----------------------------------------------------------------
!
      if(io_pe.eq.1) then
        write(11,*) 'READ_CONF: Parameter read... start'
      end if
!*                                        this is suffix0 !!
      OPEN (unit=08,file=praefixs//'_config.START'//suffix0, &
                                            form='formatted')
!
      read (8,'(a40,a6)')   text1, praefix  ! String der simulationserkennung
!     read (8,'(a40,i12)')   text1,ifqq     ! Non-neutral if ifqq= 1
      read (8,'(a40,f20.0)') text1,cptot    ! Maximum cpu time for each run
      read (8,'(a40,f20.0)') text1,dt       ! Zeitschritt
      read (8,'(a40,f20.0)') text1,tmax     ! Zeit zum Abbruch
      read (8,'(a40,f20.0)') text1,dtwr     ! Write out interval for IWRT1 
      read (8,'(a40,f20.0)') text1,dtwr2    ! Write out for IWRT2 
!
      if(io_pe.eq.1) then
        write(11,*) 'cptot (min)=',cptot
        write(11,*) 'tmax =',tmax
        write(11,*) 'dt   =',dt
      end if
!
      read (8,'(a40,f20.0)') text1,Qcore    ! charge of a macroion 
      read (8,'(a40,f20.0)') text1,Rmac     ! radius of a macroion 
      read (8,'(a40,f20.0)') text1,Wmac     ! mass of a macroion  *added
!
      if(io_pe.eq.1) then
        write(11,*) 'Qcore =',Qcore
        write(11,*) 'Rmac  =',Rmac
        write(11,*) 'Wmac  =',Wmac
      end if
!
      read (8,'(a40,i12)') text1,np         ! Zahl der Macroionen
      read (8,'(a40,i12)') text1,N_P        ! Zahl der Polymere
      read (8,'(a40,i12)') text1,N_LP       ! Zahl der Ladungen pro Polymer
      Ladabst= 1
!
      if(io_pe.eq.1) then
        write(11,*) 'np (Macroionen) =',np
        write(11,*) 'N_P (Polymere)  =',N_P
        write(11,*) 'N_LP(per Polymere)=',N_LP
      end if
!
      read (8,'(a40,f20.0)') text1,qfrac    ! Fraction of Polymer's non-neutrality 
      read (8,'(a40,i12)') text1,v_G        ! Valenz der Gegen(salz)ionen
      read (8,'(a40,i12)') text1,v_SM       ! Valenz der negativen Salzionen
      Zcp = v_G  
      Zcn = v_SM 
!
!  /init/ handles them 
!     read (8,'(a40,i12)') text1,npartcls   ! charged atoms: MAXPAR
!     npq= npartcls 
!     nqr= npartcls - (N_P*N_LP + np)
!     nseg = N_P
!
      if(io_pe.eq.1) then
        write(11,*) 'Zcp =',Zcp
        write(11,*) 'Zcn =',Zcn
        write(11,*)
      end if
!
      read (8,'(a40,f20.0)') text1,length
      read (8,'(a40,i12)')   text1,isize
      read (8,'(a40,f20.0)') text1,skin
!
      if(io_pe.eq.1) then
        write(11,*) 'length  =',length
      end if
!
      read (8,'(a40,f20.0)') text1,rcutpme
      read (8,'(a40,f20.0)') text1,acount
      read (8,'(a40,f20.0)') text1,acoion
      read (8,'(a40,f20.0)') text1,epsil_p
      read (8,'(a40,f20.0)') text1,epsil_n
      read (8,'(a40,f20.0)') text1,temperat
      read (8,'(a40,f20.0)') text1,rcutLJ
!
      epsil_p = 1.661d-14*epsil_p
      epsil_n = 1.661d-14*epsil_n
!
      if(io_pe.eq.1) then
        write(11,*) 'rcutpme =',rcutpme
        write(11,*) 'acount  =',acount
        write(11,*) 'acoion  =',acoion
        write(11,*) 'epsil_p (erg) =',epsil_p
        write(11,*) 'epsil_n (erg) =',epsil_n
        write(11,*) 'rcutLJ  =',rcutLJ
      end if
!
      read (8,'(a40,f20.0)') text1,alpha
!     read (8,'(a40,i12)')   text1,MESH  ! parm_inv3.h
!     read (8,'(a40,i12)')   text1,IP    ! parm_inv3.h
      read (8,'(a40,f20.0)') text1,Bjerrum
      read (8,'(a40,f20.0)') text1,Edc0
!
      econv= t_unit**2*e_unit /(w_unit*a_unit *299.98d0)
      Edc  = econv * Edc0
!
      if(io_pe.eq.1) then
        write(11,*) 'alpha (~/length) =',alpha  ! ~ 1/length
        write(11,*) 'Bjerrum      =',Bjerrum
        write(11,*) 'econv * Edc0 =',Edc
      end if
!     
      close(08)
!
      if(io_pe.eq.1) then
        write(11,*) 'READ_CONF: Parameter read... end'
        write(11,*)
      end if
!
!     ORTE und Geschwindigkeiten einlesen
!  ***********************
!
!     HIER EVENTUELL KONSISTENZCHECK
!
      if (mesh.ne.maxmesh) then
        write(06,*) ' mesh, maxmesh= ',mesh,maxmesh
        print*, 'Die programminterne Fouriergitter-Groesse'
        print*, 'und die eingelesene stimmen nicht ueberein'
        print*, 'Programm abgebrochen!'
        stop
      end if
      goto 1000
!
 999  print*, 'Configuration has less particles than allowed'
      print*,' program terminated'
      stop
!
 1000 continue
      return
      end subroutine READ_CONF
!
!
!-------------------------------------------------------------
      subroutine Initialisierung
!-------------------------------------------------------------
      use, intrinsic :: iso_c_binding
      implicit none
!
      include  'param_inv3.h'
!
      real(C_DOUBLE) alpha,length,dmesh,prefactor,pref_eps,econv
      integer(C_INT) npartcls
      common/ewald1/ alpha,length,dmesh,prefactor,pref_eps,econv
      common/ewald2/ npartcls
!
      real(C_DOUBLE) rcutlj2,temperat,rcutpme2,rcps2,         &
                     rcutpme,rcutLJ,r_fene,k_fene,r_fene2,  &
                     Bjerrum
      integer(C_INT) i_Steps,MeasStep,ConfStep,             &
                     N_LP,v_G,n_SMol,v_SP,v_SM,seed
      common /cutoffLJ/ rcutlj2,temperat
      common /cutoffel/ rcutpme2,rcps2
      common /CONFdataR/ rcutpme,rcutLJ,r_fene
      common /fenepara/ k_fene,r_fene2
      common /elsta/ Bjerrum
      common /confdatai/ i_Steps,MeasStep,ConfStep,   &
                         N_LP,v_G,n_SMol,v_SP,v_SM,seed

      real(C_DOUBLE) e_c_s,e_c_pme,e_c_r,e_LJ,e_fene
      common /Energy/ e_c_s,e_c_pme,e_c_r,e_LJ,e_fene
!
      integer    N_G,N_SP,N_SM,N_SI,N_P,MPC,LADABST,  &
                 qlj,MAXION,PLUSION
      common /INIDATI/ N_G,N_SP,N_SM,N_SI
      common /POLY/ N_P,MPC,LADABST
      common /ION/ MAXION,PLUSION
!
      real(C_DOUBLE) qwert,fmesh,skin,skin2
      common /INIDATR/ qwert,fmesh
      common /chargeLJ/ qlj(npq0)
      common /einles/ skin,skin2
!
      real(C_DOUBLE) rcps,wupi
      integer Verletfix,i,isize,IVERG
!
      real(C_DOUBLE) pi,dt,rbmax,vth,vth1,vth2,vth3,tmax
      common/parm2/ pi,dt,rbmax,vth,vth1,vth2,vth3,tmax
!----------------------------------------------------------------

      wupi   = 1.77245385090551602729816748334D0

      N_G   = N_P*N_LP/v_G              ! Zahl der Polymer-Gegenionen 
      N_SP  = N_SMol*v_SM               ! Zahl der positiven Salzionen 
      N_SM  = N_SMol*v_SP               ! Zahl der negativen Salzionen 
      N_SI  = N_SP+N_SM                 ! Gesamtzahl der Salzionen 
      MPC   = (N_LP - 1)*Ladabst+1      ! Zahl der Monomere pro 
!                                         Polymerkette 
      qwert = 1.0D0*(N_P*N_LP)        & ! Summe der Quadrate aller Ladungen 
               +dfloat(v_G)**2  * ((N_P*N_LP)/v_G) & 
               +dfloat(v_SP)**2 * N_SP             &
               +dfloat(v_SM)**2 * N_SM

      MAXION  = N_G + N_SP + N_SM
      PLUSION = N_G + N_SP
  
      rcps      = rcutpme + skin
      fmesh     = dfloat(mesh)
  
      skin2     = skin   * skin
      rcutpme2  = rcutpme* rcutpme
      r_fene2   = r_fene * r_fene
      rcutLJ2   = rcutLJ * rcutLJ
      rcps2     = rcps   * rcps
  
      do i=1,N_P*MPC
          qlj(i) = 1 
      end do
      do i=N_P*MPC +1,npartcls
          qlj(i) = 0
      end do

      IF (rcutpme.LT.rcutLJ) THEN
        print*,'PME-CUT-OFF ZU KLEIN:',rcutpme,'->rcutLJ',rcutLJ
        rcutpme = rcutLJ
      END IF
!
      if (rcps.gt.(0.5*length)) then
        print*,'rcps > 0.5*length: Minimum convention not satisfied'
      end if
!
!   Coulomb Self Energy Term
      e_c_s = alpha*Bjerrum*qwert/wupi
!
      return
      end subroutine Initialisierung
!
!
!--------------------------------------------------------------
      subroutine init (x,y,z,vx,vy,vz,ch,am,ag,ep,    &
                       np,nq,nr,npqr)     ! neutrality assumed
!--------------------------------------------------------------
!************************
!*  Periodic version.   *
!************************
      use, intrinsic :: iso_c_binding
      implicit none
      include 'param_inv3.h'
!
      real(C_DOUBLE),dimension(npqr0) :: &
                       x,y,z,vx,vy,vz,ch,am,ag,ep
      integer(C_int) np,nq,nCLp,nr,npqr,i,j,k,nq1,nq2
      integer(C_int) io_pe,num_proc
      common/sub_proc/ io_pe,num_proc
!
      real(C_DOUBLE) t8
      common/time3/  t8
!
      real(C_DOUBLE) t_unit,a_unit,w_unit,e_unit, &
                     kjoule,kcal,mol,kbT
      common/units/ t_unit,a_unit,w_unit,e_unit
      common/unit2/ kjoule,kcal,mol,kbT
!
      real(C_DOUBLE) pi,dt,rbmax,vth,vth1,vth2,vth3,tmax, &
                     xmax,ymax,zmax,pi2
      real(C_float) ranff,vmax1,vmax2,dgaus2,s0,s1,s2, &
                    phi,tht,dtwr,dtwr2,t
!
      integer(C_INT) it,is,nsg,nseg
      common/parm1/ it,is
      common/parm2/ pi,dt,rbmax,vth,vth1,vth2,vth3,tmax
      common/parm3/ xmax,ymax,zmax,phi,tht,dtwr,dtwr2
      common/psegm/ nsg(30),nseg
!
      real(C_DOUBLE) qfrac,Qcore,Rmac,Wmac,Zcp,Zcn,awat
      common/cntion/ qfrac,Qcore,Rmac,Wmac,Zcp,Zcn
!
      real(C_DOUBLE) acount,acoion,epsil_p,epsil_n,epsil_j,epsil_w
      common/ionsiz/ acount,acoion,epsil_p,epsil_n,epsil_j,epsil_w   
      character*2    spc(npq0)  !! charged only
!
      real(C_DOUBLE) rcutpme,rcutLJ,r_fene
      common /CONFdataR/ rcutpme,rcutLJ,r_fene
!
      real(C_DOUBLE) Qmacro,Qnegative,Qpositive,QQQ, &
                     gx,gy,gz,hx,hy,hz,hxi,hyi,hzi,svx1,svy1,svz1
      integer(C_INT) ncti
      common/ktable/ gx(mx+1),gy(my+1),gz(mz+1),  &
                     hx,hy,hz,hxi,hyi,hzi
!
      integer(C_INT) pxr,pxc,pxl,pyr,pyc,pyl,pzr,pzc,pzl
      common/ptable/ pxr(mx1),pxc(mx1),pxl(mx1),pyr(my1),pyc(my1),&
                     pyl(my1),pzr(mz1),pzc(mz1),pzl(mz1)
!
      pi2= 2.*pi
!
      hx= xmax/mx
      hy= ymax/my
      hz= zmax/mz
!
      hxi= 0.99999999d0/hx
      hyi= 0.99999999d0/hy
      hzi= 0.99999999d0/hz
!
!-------------------------
!*  Periodic tables.
!-------------------------
!
      do 110 i= 1,mx+1
      gx(i)= hx*(i-1)
!
      pxr(i)= i+1
      pxc(i)= i
      pxl(i)= i-1
  110 continue
!
      pxr(mx)=  1
      pxl( 1)= mx
      pxr(mx+1)= 2
      pxc(mx+1)= 1
!
!
      do 120 j= 1,my+1
      gy(j)= hy*(j-1)
!
      pyr(j)= j+1
      pyc(j)= j
      pyl(j)= j-1
  120 continue
!
      pyr(my)=  1
      pyl( 1)= my
      pyr(my+1)= 2
      pyc(my+1)= 1
!
!
      do 130 k= 1,mz+1
      gz(k)= hz*(k-1)
!
      pzr(k)= k+1
      pzc(k)= k
      pzl(k)= k-1
  130 continue
!
      pzr(mz)=  1
      pzl( 1)= mz
      pzr(mz+1)= 2
      pzc(mz+1)= 1
!
!----------------------------------------------------------
!* Physical constants: Run_MD  realteil READ_CONF init 
!                      L.210-, L.1000-  L.1840-,  L.2200-
!----------------------------------------------------------
!  Domain (-L/2, L/2)
!
!* DNA (Adenosine) .... M= 312/e 
!       (5+5)*C + 5*O +P +(7+5)*H +5*N; C=12, O=14, N=14, P=30
!  Al  .... 26/3e,  Ca .... 38/2e7
!  New data: in erg (Matsumoto group, 2020)
!     1 mol=  6.0221 10^23 particles
!
      epsil_j = 1.661d-14   ! 1.6606 10^-14 erg
      epsil_w = 3.188d-14   ! TIP4P Ewald sum
!
      if(io_pe.eq.1) then
        write(11,*) 'epsil_j = 1.66d-14    ! 1.66 10^-14 erg'
        write(11,*) 'epsil_w = 3.188d-14   ! TIP4P-2005 Ewald sum'
        write(11,*)
      end if
!
      if(kstart.gt.0) return
!     **********************
!
      np = np0     ! macroions
!
      nq1= nq10    ! 60 (counterions against macro) +100 counter 
      nq2= nq20    ! 300 coions
      nq = nq1 +nq2
!
      nr= nr0      ! 8000 water molecules
!
      npqr= np + nq + nr
      nCLp= np + nq 
!     ++++++++++++++++++
!
!* Macroions
      do i= 1,np
      x(i) = (-0.5d0 +ranff(0.))*xmax
      y(i) = (-0.5d0 +ranff(0.))*ymax 
      z(i) = (-0.5d0 +ranff(0.))*zmax
!
      ch(i)= Qcore
      am(i)= Wmac  ! 34* 200.d0 
      ag(i)= Rmac  !  5a*2.1 A
      ep(i)= 3 *epsil_j 
!
      if(io_pe.eq.1) then
        write(11,330) i,ch(i),am(i),ag(i),ep(i)
  330   format(' i,ch,am,ag,ep(i)=',i5,1p4d12.4)
      end if
      end do
!
!* Counter and co-ions
      do i= np+1,nCLp
      x(i) = (-0.5d0 +ranff(0.))*xmax
      y(i) = (-0.5d0 +ranff(0.))*ymax 
      z(i) = (-0.5d0 +ranff(0.))*zmax
!
      if(i.le.np+nq1) then
        ch(i)= Zcp 
        am(i)= 22.d0   ! counterion small
        ag(i)= acount
        ep(i)= epsil_j
        spc(i)= 'K '
!
      else if(i.gt.np+nq1) then
        ch(i)= Zcn
        am(i)= 34.d0   ! coion larger
        ag(i)= acoion
        ep(i)= epsil_j
        spc(i)= 'Cl'
      end if
!
      if(io_pe.eq.1) then
        write(11,330) i,ch(i),am(i),ag(i),ep(i)
      end if
      end do
!
!* Water molecules
      do i= nCLp+1,npqr
      x(i) = (-0.5d0 +ranff(0.))*xmax
      y(i) = (-0.5d0 +ranff(0.))*ymax 
      z(i) = (-0.5d0 +ranff(0.))*zmax
!
      ch(i)= 0
      am(i)= 18.d0 
      ag(i)= 3.16d0/2.d0  ! radius 
      ep(i)= epsil_w  !! for water, epsil_w
!
      if(io_pe.eq.1) then
        if(i.gt.nCLp .and. i.le.nCLp+5) then
        write(11,330) i,ch(i),am(i),ag(i),ep(i)
        end if
      end if
      end do
!
!---------------------------------
!   (First) Give charge states.
!---------------------------------
!* For counterions (Qcore <0, Zcp > 0).
!    different valence for negative particles.
!         ****
      Qmacro   =  Qcore *np 
      Qpositive=  Zcp *nq1   !! Zcp > 0
      Qnegative=  Zcn *nq2   !! Zcn < 0
!
      QQQ= Qmacro + Qnegative + Qpositive
      if(abs(QQQ).lt.1.d-8) then
        if(io_pe.eq.1) then
          write(11,*) 'Macroions Qmacro, np=',Qmacro,np
          write(11,*) 'Qpositive nq1=',Qpositive,nq1
          write(11,*) 'Qnegative nq2=',Qnegative,nq2
          write(11,*) 'QQQ is ok as QQQ=',QQQ
          write(11,*)
        end if
        go to 300
      else
        if(io_pe.eq.1) then
          write(11,*) 'Macroions Qmacro, np=',Qmacro,np
          write(11,*) 'Qnegative Qnegative, nq1=',Qnegative,nq1
          write(11,*) 'Qpositive Qpositive, nq2=',Qpositive,nq2
          write(11,*) 'QQQ is not good as QQQ=',QQQ
        end if
        return
      end if
  300 continue
!
! Within the (-L/2,+L/2) system
      do i= 1,npqr
      x(i) = x(i) -anint(x(i)/xmax)*xmax ! real of dnint()
      y(i) = y(i) -anint(y(i)/ymax)*ymax
      z(i) = z(i) -anint(z(i)/zmax)*zmax
      end do
!
      if(io_pe.eq.1) then
        write(11,*) 'Number of Macroion, counter, coion, water:'
        write(11,*) '  np, nq, nr=',np,nq,nr
        write(11,*) '  npqr=',npqr
!
        write(11,550) Qcore,Rmac,Wmac,Zcp,Zcn,QQQ
  550   format(' Qcore, Rmac, Wmac =',3f7.1,/,    & 
               ' Zcp, Zcn, QQQ= ',2f7.1,d15.6)
        write(11,551) acount,acoion
  551   format(' Size of counterion=',f7.1,  &
             '           coion     =',f7.1)
      end if
!
!***************
!*  Velocity.  *
!***************
!*   <(1/2)*m v**2>= (3/2) kT, n a**3=1.
!--------------------------------------------
      t8   = 0.d0
      t    = t8
!
      vth  = sqrt(2.d0*kbT/(Wmac*w_unit)) /(a_unit/t_unit)
      vth1 = sqrt(2.d0*kbT/(am(np+nq1)*w_unit)) /(a_unit/t_unit)
      vth2 = sqrt(2.d0*kbT/(am(np+nq2)*w_unit)) /(a_unit/t_unit)
      vth3 = sqrt(2.d0*kbT/(am(np+nq+1)*w_unit))/(a_unit/t_unit)
!--------------------------------------------
!
      do i= 1,np
      vmax1= vth
      vx(i)= 0 !dgaus2(vmax1)
      vy(i)= 0 !dgaus2(vmax1)
      vz(i)= 0 !dgaus2(vmax1)
      end do
!
      do i= np+1,nCLp
      if(i.le.np+nq1) then
        vmax2= vth1
      else
        vmax2= vth2
      end if
!
      vx(i)= dgaus2(vmax2)
      vy(i)= dgaus2(vmax2)
      vz(i)= dgaus2(vmax2)
      end do
!
      do i= nCLp+1,npqr
      vmax2= vth3
      vx(i)= dgaus2(vmax2)
      vy(i)= dgaus2(vmax2)
      vz(i)= dgaus2(vmax2)
      end do
!
!* Nullify counterions and coions
      svx1= 0
      svy1= 0
      svz1= 0
      ncti= 0
!
      do i= np+1,np+nq1  ! Counterions
      svx1= svx1 +vx(i)
      svy1= svy1 +vy(i)
      svz1= svz1 +vz(i)
      ncti= ncti +1
      end do
!
      svx1= svx1/ncti
      svy1= svy1/ncti
      svz1= svz1/ncti
!
      do i= np+1,np+nq1
      vx(i)= vx(i) -svx1
      vy(i)= vy(i) -svy1
      vz(i)= vz(i) -svz1
      end do
!
!
      svx1= 0
      svy1= 0
      svz1= 0
      ncti= 0
!
      do i= np+nq1,npq0  !<-- Coions
      svx1= svx1 +vx(i)
      svy1= svy1 +vy(i)
      svz1= svz1 +vz(i)
      ncti= ncti +1
      end do
!
      svx1= svx1/ncti
      svy1= svy1/ncti
      svz1= svz1/ncti
!
      do i= np+nq1,npq0  !<--
      vx(i)= vx(i) -svx1
      vy(i)= vy(i) -svy1
      vz(i)= vz(i) -svz1
      end do
!
!* Nulify water molecules
      svx1= 0
      svy1= 0
      svz1= 0
      ncti= 0
!
      do i= nCLp+1,npqr
      svx1= svx1 +vx(i)
      svy1= svy1 +vy(i)
      svz1= svz1 +vz(i)
      ncti= ncti +1
      end do
!
      svx1= svx1/ncti
      svy1= svy1/ncti
      svz1= svz1/ncti
!
      do i= nCLp+1,npqr
      vx(i)= vx(i) -svx1
      vy(i)= vy(i) -svy1
      vz(i)= vz(i) -svz1
      end do
!
!
      if(io_pe.eq.1) then
        s0= 0
        do i= 1,np
        s0= s0 +am(i)*(vx(i)**2 +vy(i)**2 +vz(i)**2)
        end do
!
        s1= 0
        do i= np+1,nCLp 
        s1= s1 +am(i)*(vx(i)**2 +vy(i)**2 +vz(i)**2)
        end do
!
        s2= 0
        do i= nCLp+1,npqr
        s2= s2 +am(i)*(vx(i)**2 +vy(i)**2 +vz(i)**2)
        end do
!
        write(11,*) '*Energy per atom and neutral=',  &
                     s0/np,s1/nq,s2/nr
        write(11,*)
      end if
!
      return
      end subroutine init
!
!
!---------------------------------------------------------------
      subroutine ggauss
!---------------------------------------------------------------
      use, intrinsic :: iso_c_binding
      implicit none
!
      real(C_float) fv,vv0,dv,vv,s,sdv,fun
      integer(C_INT) i,j,k,ns,k2
      common/gaus1/ fv(51),vv0,dv
!
      FV(1)=0.0
!
      VV0= -3.
      DV= 2.*ABS(VV0)/50.0
!
      VV= VV0
      DO J=1,50
      S=0.0
      NS=1000
      K2=NS/2
      SDV=DV/FLOAT(NS)
!
      DO K=1,K2
      VV=VV +2.0*SDV
      S=S +4.0*FUN(VV-SDV) +2.0*FUN(VV)
      end do
!
      S= (S +4.0*FUN(VV+SDV) +FUN(VV+2.0*SDV))*SDV/3.0
      FV(J+1)= FV(J)+S
      end do
!
      DO I=1,51
      FV(I)=FV(I)/FV(51)
      end do
!
      return
      end subroutine ggauss
!
!
!---------------------------------------------------------------
      function dgaus2 (vmax)
!---------------------------------------------------------------
      use, intrinsic :: iso_c_binding
!
      integer(C_INT) k,k2
      real(C_float) dgaus2,vmax,fv,vv0,dv,eps,ranff
      common/gaus1/ fv(51),vv0,dv
!
      eps= ranff(0.)
      DO K=1,51
      K2=K
      IF(FV(K).GT.eps) GO TO 200
      end do
!
  200 Y1= FV(K2-1)
      Y2= FV(K2)
      X2= (EPS-Y2)/(Y2-Y1)+K2
      dgaus2= vmax*(VV0 +DV*(X2-1.0))
!
      return
      end function dgaus2 
!
!
!---------------------------------------------------------------
      function fun (v)
!---------------------------------------------------------------
      use, intrinsic :: iso_c_binding
      implicit none
      real(C_float) fun,v
!
      fun= exp(-v**2/2.)
!
      return
      end function fun
!
!
!-----------------------------------------------------------
      subroutine vdistr (vx,vy,vz,am,np,nq)
!-----------------------------------------------------------
      use, intrinsic :: iso_c_binding
      implicit none
      include  'param_inv3.h'
!
      real(C_DOUBLE) vx(npq0),vy(npq0),vz(npq0),am(npq0)
      real(C_float)  fvx(51),fvy(51),fvz(51),xsc(51)
      integer(C_INT) np,nq
!
!     character*1  char(8**3)
      real(C_DOUBLE) pi,dt,rbmax,vth,vth1,vth2,vth3,tmax,awat
      common/parm2/ pi,dt,rbmax,vth,vth1,vth2,vth3,tmax
!
      real(C_float) vmax1,vmax2,aiv,fmax1,fmax2,fmax3, &
                    fmin1,fmin2,fmin3
      integer(C_INT) ix,iy,iz,i,k,ILN
!
!*******************************
!*  For the chain parts only.  *
!*******************************
!
      vmax1= 6.*vth /sqrt(am(1))
      vmax2= 6.*vth1/sqrt(am(np+1))
!
!* (1) Coulomb particles.
      aiv= 25./vmax1
!
      do k= 1,51
      xsc(k)= (k -26)/aiv
      fvx(k)= 0.
      fvy(k)= 0.
      fvz(k)= 0.
      end do
!
      do i= 1,np
      ix= aiv*vx(i) +25.5
      iy= aiv*vy(i) +25.5
      iz= aiv*vz(i) +25.5
!
      if(iabs(ix-26).gt.25)  go to 230
        fvx(ix)= fvx(ix) +1.
!
  230 if(iabs(iy-26).gt.25)  go to 250
        fvy(iy)= fvy(iy) +1.
!
  250 if(iabs(iz-26).gt.25)  go to 200
        fvz(iz)= fvz(iz) +1.
  200 continue
      end do
!
      call lplmax (fvx,fmax1,fmin1,51)
      call lplmax (fvy,fmax2,fmin2,51)
      call lplmax (fvz,fmax3,fmin3,51)
!
      ILN= 1
      call hplot1 (2,4,51,xsc,fvx,fmax1,fmin1,ILN,'FVX(mac)',8, &
                   '   VX   ',8,'        ',8)
      call hplot1 (2,5,51,xsc,fvy,fmax2,fmin2,ILN,'FVY(mac)',8, &
                   '   VY   ',8,'        ',8)
      call hplot1 (2,6,51,xsc,fvz,fmax3,fmin3,ILN,'FVZ(mac)',8, &
                   '   VZ   ',8,'        ',8)
!
!* (2) Water particles.
      aiv= 25./vmax2
!
      do k= 1,51
      xsc(k)= (k -26)/aiv
      fvx(k)= 0.
      fvy(k)= 0.
      fvz(k)= 0.
      end do
!
      do i= np+1,np+nq
      ix= aiv*vx(i) +25.5
      iy= aiv*vy(i) +25.5
      iz= aiv*vz(i) +25.5
!
      if(iabs(ix-26).gt.25)  go to 430
        fvx(ix)= fvx(ix) +1.
!
  430 if(iabs(iy-26).gt.25)  go to 450
        fvy(iy)= fvy(iy) +1.
!
  450 if(iabs(iz-26).gt.25)  go to 400
        fvz(iz)= fvz(iz) +1.
  400 continue
      end do
!
      call lplmax (fvx,fmax1,fmin1,51)
      call lplmax (fvy,fmax2,fmin2,51)
      call lplmax (fvz,fmax3,fmin3,51)
!
      ILN= 1
      call hplot1 (3,4,51,xsc,fvx,fmax1,fmin1,ILN,'FVX(nqr)',8, &
                   '   VX   ',8,'        ',8)
      call hplot1 (3,5,51,xsc,fvy,fmax2,fmin2,ILN,'FVY(nqr)',8, &
                   '   VY   ',8,'        ',8)
      call hplot1 (3,6,51,xsc,fvz,fmax3,fmin3,ILN,'FVZ(nqr)',8, &
                   '   VZ   ',8,'        ',8)
!---------------------
      call chart
!---------------------
      return
      end subroutine vdistr
!
!
!------------------------------------------------------
      subroutine averg1 (q,qav,is)
!------------------------------------------------------
      use, intrinsic :: iso_c_binding
      dimension  q(4096)
!
      qav= 0.
!
      do i= is-9,is
      qav= qav +q(i)
      end do
!
      qav= qav/10.
!
      return
      end subroutine averg1 
!
!
!------------------------------------------------------
      subroutine lplots
!------------------------------------------------------
      use, intrinsic :: iso_c_binding
      implicit none
!
      include 'param_inv3.h'
      real(C_float),dimension(nhist) :: &
                   ekin,ekn2,ekn3,ecr,elj,ep3m,pot,time,           &
                   vxio,vxco,vxct,vxmc,vxtot,xmcc,vxwat,etot,      &
                   xmac,ymac,zmac,xcon,ycon,zcon,xcoi,ycoi,zcoi,   &
                   xwat,ywat,zwat,dmcc
      common/ehist/ ekin,ekn2,ekn3,ecr,elj,ep3m, pot,time,         &
                    vxio,vxco,vxct,vxmc,vxtot,xmcc,vxwat,etot,     &
                    xmac,ymac,zmac,xcon,ycon,zcon,xcoi,ycoi,zcoi,  &
                    xwat,ywat,zwat,dmcc
!
!--------------------------------------------------------------
!* Read basic data in /READ_CONF/, constants for realteil and P3M.
!    >> Bjerrum
!    >> dt, length, alpha, rcutpme, rcutlj
!    >> r_fene, k_fene 
!    >> cptot, dtwr, dtwr2, rmax.
!--------------------------------------------------------------
!  pgf95 @lplot240.f03     
!
      character(len=8) label,cdate*10,ctime*8
      common/headr1/  label,cdate,ctime
!
      real(C_float) t,xleng
      integer(C_INT) it,is,nsg,nseg,ILN,ILG,i
      common/headr2/ t,xleng
      common/parm1/ it,is
      common/psegm/ nsg(30),nseg
!
      real(C_float) emax1,emin1,emax2,emin2,pmax,pmin,etmax,etmin, &
                    ecmax,ecmin,elmax,elmin,e3max,e3min,           &
                    vxmc1,vxmc2,xmcc1,xmcc2,vxt1,vxt2,vxio1,vxio2, &
                    vxct1,vxct2,vxco1,vxco2,emax
!
      ekin(1)= ekin(3)
      ekin(2)= ekin(3)
      pot (1)=  pot(3)
      pot (2)=  pot(3)
      etot(1)= etot(3)
      etot(2)= etot(3)
      ep3m(1)= ep3m(3)
      ep3m(2)= ep3m(3)
!
        do i= 1,5
        elj(i)= elj(5)
        ep3m(i)= ep3m(5)
        end do
!
      call lplmax (ekin,emax1,emin1,is)
      call lplmax (ekn2,emax2,emin2,is)
      call lplmax ( pot,pmax,pmin,is)
      call lplmax ( ecr,ecmax,ecmin,is)
!
      call lplmax ( elj,elmax,elmin,is)
      call lplmax (ep3m,e3max,e3min,is)
      call lplmax (etot,etmax,etmin,is)
!
      emax= max(emax1,emax2)
!
      ILN= 1
      ILG= 2
!
!  pgf95 @lplot240.f03     
!     call lplot1 (2,4,is,time,ekin,emax1,emin1,ILN,'Kin Macr',8, &
      call lplot1 (2,4,is,time,ekin,emax1,0.,ILN,'Kin Macr',8, &
                 '        ',8,'        ',8)
!     call lplot1 (2,5,is,time,ekn2,emax2,emin2,ILN,'K.cou+co',8, &
      call lplot1 (2,5,is,time,ekn2,emax2,0.,ILN,'K.cou+co',8, &
                 '        ',8,'        ',8)
      call lplot1 (2,6,is,time,etot,etmax,etmin,ILN,'E.total ',8, &
                 '  time  ',8,'        ',8)
!
      call lplot1 (3,4,is,time, ecr,ecmax,ecmin,ILN,'E_realtl',8, &
                 '        ',8,'        ',8)
      call lplot1 (3,5,is,time, elj,elmax,elmin,ILN,'E_LJ    ',8, &
                 '        ',8,'        ',8)
      call lplot1 (3,6,is,time,ep3m,e3max,e3min,ILN,'E_p3m   ',8, &
                 '  time  ',8,'        ',8)
!------------------------
      call chart
!------------------------
!       dmcc(is) = ss    ! charge distribution with Rmac+5 A
!       vxtot(is)= svx1  ! counter + co
!       vxmc(is) = svx2  ! counter only
!
!     call lplmax (xmcc,xmcc1,xmcc2,is)
      call lplmax (dmcc,xmcc1,xmcc2,is)
      call lplmax (vxmc,vxmc1,vxmc2,is)
      call lplmax (vxtot,vxt1,vxt2, is)
!
      call lplmax (vxio,vxio1,vxio2,is)
      call lplmax (vxct,vxct1,vxct2,is)
      call lplmax (vxco,vxco1,vxco2,is)
!
      call lplot1 (2,4,is,time,dmcc,xmcc1,xmcc2,ILN,'S.Rch<8A',8, & ! distr
                 '        ',8,'        ',8)
      call lplot1 (2,5,is,time,vxmc,vxmc1,vxmc2,ILN,'S.Vx Mc ',8, & ! Vx Mc
                 '        ',8,'        ',8)
      call lplot1 (2,6,is,time,vxtot, vxt1,vxt2,ILN,'S.Vx Mcc',8, &
                 '  time  ',8,'        ',8)
!
      call lplot1 (3,4,is,time,vxio,vxio1,vxio2,ILN,'VX(MACI)',8, & ! vx M 
                 '        ',8,'        ',8)
      call lplot1 (3,5,is,time,vxct,vxct1,vxct2,ILN,'VX(CNTI)',8, & ! vx count
                 '        ',8,'        ',8)
      call lplot1 (3,6,is,time,vxco,vxco1,vxco2,ILN,'VX(COIO)',8, & ! vx co
                 '  time  ',8,'        ',8)
!------------------------
      call chart
!------------------------
!
!  pgf95 @lplot240.f03     
      call lplmax (xmac,xmcc1,xmcc2,is)
      call lplmax (ymac,vxmc1,vxmc2,is)
      call lplmax (zmac,vxt1,vxt2,  is)
      call lplmax (xcon,vxio1,vxio2,is)
      call lplmax (ycon,vxct1,vxct2,is)
      call lplmax (zcon,vxco1,vxco2,is)
!
      call lplot1 (2,4,is,time,xmac,xmcc1,xmcc2,ILN,'X.MAC   ',8, &
                 '        ',8,'        ',8)
      call lplot1 (2,5,is,time,ymac,vxmc1,vxmc2,ILN,'Y.MAC   ',8, &
                 '        ',8,'        ',8)
      call lplot1 (2,6,is,time,zmac, vxt1,vxt2,ILN, 'Z.MAC    ',8, &
                 '  time  ',8,'        ',8)
!
      call lplot1 (3,4,is,time,xcon,vxio1,vxio2,ILN,'X.COUNT ',8, &
                 '        ',8,'        ',8)
      call lplot1 (3,5,is,time,ycon,vxct1,vxct2,ILN,'Y.COUNT ',8, &
                 '        ',8,'        ',8)
      call lplot1 (3,6,is,time,zcon,vxco1,vxco2,ILN,'Z.COUNT ',8, &
                 '  time  ',8,'        ',8)
!------------------------
      call chart
!------------------------
!
      call lplmax (xcoi,xmcc1,xmcc2,is)
      call lplmax (ycoi,vxmc1,vxmc2,is)
      call lplmax (zcoi,vxt1,vxt2,  is)
      call lplmax (xwat,vxio1,vxio2,is)
      call lplmax (ywat,vxct1,vxct2,is)
      call lplmax (zwat,vxco1,vxco2,is)
!
      call lplot1 (2,4,is,time,xcoi,xmcc1,xmcc2,ILN,'X.COION ',8, &
                 '        ',8,'        ',8)
      call lplot1 (2,5,is,time,ycoi,vxmc1,vxmc2,ILN,'Y.COION ',8, &
                 '        ',8,'        ',8)
      call lplot1 (2,6,is,time,zcoi, vxt1,vxt2,ILN, 'Z.COION ',8, &
                 '  time  ',8,'        ',8)
!
      call lplot1 (3,4,is,time,xwat,vxio1,vxio2,ILN,'X.WATER ',8, &
                 '        ',8,'        ',8)
      call lplot1 (3,5,is,time,ywat,vxct1,vxct2,ILN,'Y.WATER ',8, &
                 '        ',8,'        ',8)
      call lplot1 (3,6,is,time,zwat,vxco1,vxco2,ILN, 'Z.WATER ',8, &
                 '  time  ',8,'        ',8)
!------------------------
      call chart
!------------------------
!
      write(6,*) 'lplots has been plotted'
!
      return
      end subroutine lplots
!
!
!------------------------------------------------------
      subroutine lplmax (f,fmax,fmin,is)
!------------------------------------------------------
      use, intrinsic :: iso_c_binding
      integer(C_INT) is
      dimension  f(is)
!
      fmax= -1.e10
      fmin=  1.e10
!
      do i= 1,is
      fmax= max(fmax,f(i))
      fmin= min(fmin,f(i))
      end do
!
      return
      end subroutine lplmax 
!
!
!------------------------------------------------
      function iwrta (t,twr)
!------------------------------------------------
      use, intrinsic :: iso_c_binding
      integer(C_INT) iwrta
      
      common/imemo/ iwa,iwb,iwc
!
      iw= t/twr 
      if(iw.gt.iwa) then
        iwa= iw
        iwrta= 0
      else
        iwrta= 1
      end if
!
      return
      end function iwrta 
!
!
!------------------------------------------------
      function iwrtb (t,twr)
!------------------------------------------------
      use, intrinsic :: iso_c_binding
      integer(C_INT) iwrtb
      common/imemo/ iwa,iwb,iwc
!
      iw= t/twr 
      if(iw.gt.iwb) then
        iwb= iw
        iwrtb= 0
      else
        iwrtb= 1
      end if
!
      return
      end function iwrtb 
!
!
!------------------------------------------------
      function iwrtc (t,twr)
!------------------------------------------------
      use, intrinsic :: iso_c_binding
      integer(C_INT) iwrtc
      common/imemo/ iwa,iwb,iwc
!
      iw= t/twr 
      if(iw.gt.iwc) then
        iwc= iw
        iwrtc= 0
      else
        iwrtc= 1
      end if
!
      return
      end function iwrtc
!
!
!------------------------------------------------
      block data
!------------------------------------------------
      common/ranfff/ ir,iq
      data  ir/3021/,iq/7331/    ! original
!     data  ir/17331/,iq/37711/
      end
!
!
!------------------------------------------------
      function ranff (x)
!------------------------------------------------
!*  ranf= (0,1)
      use, intrinsic :: iso_c_binding
!
      real(C_float)  x
      common/ranfff/ ir,iq
!
      real(C_DOUBLE) INVM
      parameter  (MASK=2**30+(2**30-1),INVM= 0.5D0**31)
      parameter  (LAMBDA=48828125)
!
      IR= IAND( LAMBDA*IR, MASK)
      ranff= IR*INVM
!
!     ask= 371597.
!     ambda= sqrt(ask)
!     qq= 0.3713*ask
!
!     ir= amod( ambda*ir +qq, ask)
!     ranff= ir/ask
!
      return
      end function ranff
!
!
!-----------------------------------------------------------------------
      subroutine LPLOT1 (IX,IY,npt1,x,y,YMAX,YMIN,IL,LAB1,N1,LAB2,N2,&
                         LAB3,N3)
!-----------------------------------------------------------------------
!  <<Warning>>  Order and number of arguments /LPLOT/ have been changed.
!               Also, X (time) is defined for all range.
!               Date: 5/18/96 at MIT.
!***********************************************************************
!   IL=1................ LINEAR PLOT OF (X,Y)
!   IL=2................ LOG10 PLOT OF (X,LOG Y)
!***********************************************************************
      use, intrinsic :: iso_c_binding
!
      integer(C_INT) npt1
      DIMENSION  x(npt1),y(npt1),U(4096),V(4096)
      DIMENSION  XCM(6),YCM(6),PL(6),PR(6),QL(6),QR(6)
!
      character(len=8) LAB1,LAB2,LAB3
      character(len=8) label,cdate*10,ctime*8,cax*1
      common/headr1/  label,cdate,ctime
!
      common/headr2/ time,xleng
      common/PPLCOM/ NFINE,PL1(10),PR1(10),QL1(10),QR1(10), &
                     XMIN1(10),XMAX1(10),YMIN1(10),YMAX1(10)
!
!   FOR FUJITSU.
!     data  XCM/18.46,2*9.867,3*6.18/,
!    *      YCM/16.85,2*7.435,3*4.381/,
!    *      PL/2*2.00,15.132,2.00,8.00,18.20/,
!    *      QL/1.95,10.885,1.95,13.832,7.891,1.95/
!
!   FOR NEC.
      data  XCM/21.0, 2*9.00, 3*6.00/,        &
            YCM/15.0, 2*6.80, 3*3.90/,        &
            PL/2.0,  2.0,14.0, 2.0,9.0,16.0/, &
            QL/2.3, 10.5,2.3, 12.9,7.6,2.3/
!
      IPLOT=1
      GO TO 1
!
!-----------------------------------------------------------------------
      ENTRY HPLOT1 (IX,IY,NPT1,X,Y,ymax,ymin,IL,LAB1,N1,LAB2,N2,LAB3,N3)
!-----------------------------------------------------------------------
      IPLOT=2
!
    1 NPT= NPT1
      ISC= 1
!
      DO I=1,6
      PR(I)= PL(I) +XCM(I)
      end do
!
      DO J=1,6
      QR(J)= QL(J) +YCM(J)
      end do
!
!                 ******************************************************
!*                **  MAKE A COPY BEFORE THE TOP-LEFT FRAME IS DRAWN. **
!                 ******************************************************
      HH= 0.60
      I1= IABS(IX)
      J1= IABS(IY)
      IF(I1.GE.3) GO TO 10
      IF(J1.EQ.3.OR.J1.GE.5) GO TO 10
!                                              ************************
!                                              ** LABEL OF THE PAGE. **
!                                              ************************
      CALL SYMBOL (0.1,18.7,HH,label,0.,8)
      CALL SYMBOL (3.1,18.7,HH,cdate, 0.,10)
      CALL SYMBOL (15.9,0.1,HH,'t =',0.,3)
      CALL VALUES (999.0,999.0,HH,time,0.,101)
!
   10 CONTINUE
!
      DO I=1,NPT
      U(I)= X(I)
      end do
!
      XMAX= U(NPT)
      XMIN= U(1)
!                             ************************************
!                             ** THREE-POINT AVERAGE IF IL > 0  **
!                             ************************************
      IF(IL.GT.0) THEN
        V(1)=   Y(1)
        V(NPT)= Y(NPT)
        DO I=2,NPT-1
        V(I)= 0.33333*(Y(I-1)+Y(I)+Y(I+1))
        end do
      ELSE
        DO I=1,NPT
        V(I)= Y(I)
        end do
      END IF
!                                                *****************
!                                                **  LOG. SCALE **
!                                                *****************
      IF(IABS(IL).EQ.2) THEN
         DO I=1,NPT
         IF(V(I).GT.0.) THEN
            V(I)= ALOG10(V(I))
         ELSE
            V(I)= -10.
         END IF
         end do
      END IF
!                                **************************************
!                                ** SET A NEW SCALE AND DRAW A FRAME.**
!                                **************************************
      IF(IPLOT.EQ.2) THEN
         ymax= -1.e10
         ymin=  1.e10
!
         do i= 1,npt
         ymax= max(ymax,v(i))
         ymin= min(ymin,v(i))
         end do
!
         if(ymin.ge.0.) then
           ymax= 1.1*ymax
           ymin= 0.
         else
           ymax= max(0.,ymax)
           ymin= 1.1*ymin
         end if
      END IF
!
      IF(YMAX.LE.YMIN) YMAX= YMIN+1.0
      IF(IABS(IL).EQ.2) THEN
         IF(YMAX.GT.0.0) YMAX= YMAX+1.0
      END IF
!
      DX= (XMAX-XMIN)/XCM(I1)
      DY= (YMAX-YMIN)/YCM(J1)
      X0= XMIN
      Y0= YMIN
!
      CALL SCALEX (PL(I1),QL(J1),X0,Y0,DX,DY,ISC)
!
      PL1(ISC)= PL(I1)
      PR1(ISC)= PR(I1)
      QL1(ISC)= QL(J1)
      QR1(ISC)= QR(J1)
      XMIN1(ISC)= XMIN
      XMAX1(ISC)= XMAX
      YMAX1(ISC)= YMAX
      YMIN1(ISC)= YMIN
!                                                      *************
!                                                      **  FRAME. **
!                                                      *************
      CALL PLOT (PL(I1),QL(J1),3)
      CALL PLOT (PL(I1),QR(J1),2)
      CALL PLOT (PR(I1),QR(J1),2)
      CALL PLOT (PR(I1),QL(J1),2)
      CALL PLOT (PL(I1),QL(J1),2)
!                                                    ******************
!                                                    **  TICK MARKS. **
!                                                    ******************
      SCX= XCM(I1)/4.0
      SCY= YCM(J1)/4.0
!
      X0= PL(I1)
      Y1= QL(J1)
      Y4= QR(J1)
      Y2= Y1 +0.25
      Y3= Y4 -0.25
!
      DO K=1,3
      X0= X0 +SCX
      CALL PLOT (X0,Y1,3)
      CALL PLOT (X0,Y2,2)
      CALL PLOT (X0,Y3,3)
      CALL PLOT (X0,Y4,2)
      end do
!
      Y0= QL(J1)
      X1= PL(I1)
      X4= PR(I1)
      X2= X1 +0.25
      X3= X4 -0.25
!
      DO K=1,3
      Y0= Y0 +SCY
      CALL PLOT (X1,Y0,3)
      CALL PLOT (X2,Y0,2)
      CALL PLOT (X3,Y0,3)
      CALL PLOT (X4,Y0,2)
      end do
!                                                     **************
!                                                     ** NUMBERS. **
!                                                     **************
!
      CALL NUMBER (PL(I1)-0.45,QL(J1)-0.45,0.30,XMIN,0.,101)
      CALL NUMBER (PR(I1)-0.45,QL(J1)-0.45,0.30,XMAX,0.,101)
!
      CALL NUMBER (PL(I1)-1.30,QL(J1)     ,0.30,YMIN,0.,101)
      CALL NUMBER (PL(I1)-1.30,QR(J1)-0.30,0.30,YMAX,0.,101)
!
!                                                     **************
!                                                     **  LABELS. **
!                                                     **************
      XC= 0.5*(PL(I1)+PR(I1))
      XU= XC -1.60
      XD= XC -0.20*N2/2
!
      YR= QR(J1)+0.15
      YL= QL(J1)-0.70
!
      CALL SYMBOL (XU,YR,HH,LAB1,0.,N1)
      CALL SYMBOL (XD,YL,HH,LAB2,0.,N2)
!
      XL= PL(I1)-1.50
      YC= 0.5*(QL(J1)+QR(J1))
      CALL SYMBOL (XL,YC,HH,LAB3,0.,N3)
!                                     **********************************
!                                     **  NO PLOT IS MADE IF NPT1 < 0 **
!                                     **********************************
   70 IF(NPT1.LT.0) RETURN
!
      CALL PLOTL (U(1),V(1),ISC,3)
!**
      IF(IPLOT.EQ.1) THEN
         DO I=1,NPT
         CALL PLOTL (U(I),V(I),ISC,2)
         end do
      ELSE
         DO I=1,NPT-1
         CALL PLOTL (U(I+1),V(I)  ,ISC,2)
         CALL PLOTL (U(I+1),V(I+1),ISC,2)
         end do
      END IF
!**
      CALL PLOTL (U(NPT),V(NPT),ISC,3)
!
      RETURN
      END subroutine LPLOT1 
!
!
!-----------------------------------------------------------------------
      subroutine CPLOT3 (Q,xmax8,ymax8,zmax8,CHAR,NC)
!-----------------------------------------------------------------------
!***********************************************************************
!*   CONTOUR PLOTS OF SCALAR QUANTITIES.                               *
!***********************************************************************
      use, intrinsic :: iso_c_binding
      include 'param_inv3.h'
!
      character(len=8) label,cdate*10,ctime*8,char
      common/headr1/  label,cdate,ctime
      common/headr2/ time,xleng
      integer(C_INT) pxr,pxc,pxl,pyr,pyc,pyl,pzr,pzc,pzl
      common/ptable/ pxr(mx1),pxc(mx1),pxl(mx1),pyr(my1),pyc(my1),&
                     pyl(my1),pzr(mz1),pzc(mz1),pzl(mz1)
      real(C_DOUBLE) Q(mx,my,mz),xmax8,ymax8,zmax8
      real(C_float)  a(2048),b(2048),WW(2048),cut(200,4)
!
      J0= my/2 +1
      K0= mz/2 +1
      xmax= xmax8
      ymax= ymax8
      zmax= zmax8
!
!* 1. PLOT AT K= K0: SUBSCRIPT J FIRST.
!
      npx= 0
      IJ= 0
      qc = 1./16.
!***
      do i= 1,mx
      npx= npx +1
      IR= PXR(I)
      IL= PXL(I)
!
      npy= 0
      DO j= 1,my
      npy= npy +1
      JR= PYR(J)
      JL= PYL(J)
!
      IJ= IJ+1
      a(IJ)= qc*(   Q(IR,JR,K0) +2.*Q(IR,J,K0)    +Q(IR,JL,K0)  &
                +2.*Q(I ,JR,K0) +4.*Q(I ,J,K0) +2.*Q(I ,JL,K0)  &
                +   Q(IL,JR,K0) +2.*Q(IL,J,K0)    +Q(IL,JL,K0) )
      end do
      end do
!
!* 2. PLOT AT J= J0: SUBSCRIPT K FIRST.
!
      npx= 0
      IJ= 0
      qc = 1./16.
!***
      do i= 1,mx
      npx= npx +1
      IR= PXR(I)
      IL= PXL(I)
!
      npz= 0
      DO k= 1,mz
      npz= npz +1
      KR= PZR(K)
      KL= PZL(K)
!
      IJ= IJ+1
      b(IJ)= qc*(   Q(IR,J0,KR) +2.*Q(IR,J0,K)    +Q(IR,J0,KL) &
                +2.*Q(I ,J0,KR) +4.*Q(I ,J0,K) +2.*Q(I ,J0,KL) &
                +   Q(IL,J0,KR) +2.*Q(IL,J0,K)    +Q(IL,J0,KL) )
      end do
      end do
!
!
      HH = 0.70
      CALL SYMBOL (0.1,18.2,HH,label,0.,8)
      CALL SYMBOL (13.0,0.7,HH,cdate,0.,10)
      CALL SYMBOL (13.0,0.1,HH,'t =',0.,3)
      CALL VALUES (999.0,999.0,HH,time,0.,101)
!
      XL1=  1.8
      XR1=  9.3
      XL2= 10.0
      XR2= 17.5
!
      ZL=  1.0
      ZR=  ZL +(XR1 -XL1)*ZMAX/XMAX
      if(zr.gt.10.) zr= 10.
!                          <--- Limit elongated y-length.
!
      YL=  ZR +1.
      YR=  YL +(XR1 -XL1)*YMAX/XMAX
      if(yr.gt.25.) yr= 25.
!                          <--- Limit elongated y-length.
!
      XC1= 0.5*(XR1+XL1)
      XC2= 0.5*(XR2+XL2)
      YC=  0.5*(YR+YL)
      ZC=  0.5*(ZR+ZL)
!
!---------------------------------------------
!*  **MAXIMUM OF THE VECTORS**
!---------------------------------------------
!
      AM2= 0.
      AM4= 0.
!
      DO 100 IJ= 1,npx*npy
      AM2= AMAX1(AM2,ABS(A(IJ)))
  100 CONTINUE
!
      DO 200 IJ= 1,npx*npz
      AM4= AMAX1(AM4,ABS(B(IJ)))
  200 CONTINUE
!
      AMS= AMAX1(AM2,AM4)
      IF(AMS.LT.1.E-10) AMS=999.0
!
      CALL SYMBOL (ZL,0.10,HH,'SCALAR.MAX= ',0.,12)
      CALL VALUES (999.0,999.0,HH,AMS,0.,101)
!
!---------------------------------------------
!*  (1): Contours in (x,z) plane.
!---------------------------------------------
!
      CALL SETSCL (0.,0.,ZMAX,XMAX,ZL,XL2,ZR,XR2,GDZ,GDX, &
                   NC,CHAR,6,' (X-Z)',0.4,                &
                   1,'Z',0.4,1,'X',0.4,1)
!
      CALL VALUES (ZL-0.45,XL2-0.5,HH,0.0,0.,101)
      CALL VALUES (ZL-1.3,XR2-0.3, HH,XMAX,0.,101)
      CALL VALUES (ZR-1.3,XL2-0.5, HH,ZMAX,0.,101)
!
      NXZ= npx*npz
      CALL DAISHO (B,NXZ,WAMIN,WAMAX)
!
      NCONTR= 11
      CALL EQCNTR (B,WW,npz,npx,ZL,XL2,ZR,XR2,WAMIN,0.0,WAMAX, &
                   NCONTR,1)
!
!---------------------------------------------
!*  (2): Contours in (x,y) plane.
!---------------------------------------------
!
      CALL SETSCL (0.,0.,YMAX,XMAX,YL,XL2,YR,XR2,GDY,GDX, &
                   1,' ',1,' ',0.4,                       &
                   1,'Y',0.4,1,'X',0.4,1)
!
      CALL VALUES (YL-0.45,XL2-0.5,HH,0.0,0.,101)
      CALL VALUES (YL-1.3,XR2-0.3, HH,XMAX,0.,101)
      CALL VALUES (YR-0.3,XL2-0.5, HH,YMAX,0.,101)
!
      NXY= npx*npy
      CALL DAISHO (A,NXY,WAMIN,WAMAX)
!
      NCONTR= 11
      CALL EQCNTR (A,WW,npy,npx,YL,XL2,YR,XR2,WAMIN,0.0,WAMAX, &
                   NCONTR,1)
!
!---------------------------------------------
!*  (3): Cut plots.
!---------------------------------------------
!
      do jj= 1,4
      j= (my/4)*(jj-1) +1
!
      do i= 1,mx
      cut(i,jj)= q(i,j,k0)
      end do
      end do
!
!
      amax7= -1.e+10
      amin7=  1.e+10
!
      do jj= 1,4
      do i= 1,mx
      amax7= max(cut(i,jj),amax7)
      amin7= min(cut(i,jj),amin7)
      end do
      end do
!
      if(amax7.lt.0.) amax7= 0.
      if(amin7.gt.0.) amin7= 0.
!
!
      dd= amax7 -amin7
      dx= (YR -YL)/6.
      dy= XR1 -XL1
!
      do jj= 1,4
      xo= YL +1.5*(jj-1)*dx
      xu= xo +dx
      call plot (xo,XL1,3)
      call plot (xu,XL1,2)
      call plot (xu,XR1,2)
      call plot (xo,XR1,2)
      call plot (xo,XL1,2)
!
!* zero line.
      x1= xo +dx*(0. -amin7)/dd
      call plot (x1,XL1,3)
      call plot (x1,XR1,2)
!
      x1= xo +dx*(cut(1,jj) -amin7)/dd
      y1= XL1 
      call plot (x1,y1,3)
!
      do i= 1,mx
      x1= xo  +dx*(cut(i,jj) -amin7)/dd
      y1= XL1 +dy*(i-1)/float(mx-1)
!
      call plot (x1,y1,2)
      end do
      end do
!
!---------------------
      call chart
!---------------------
!
      RETURN
      END subroutine CPLOT3 
!
!
!-----------------------------------------------------------------------
      subroutine EQCNTR(U,W,NX,NY,XL,YL,XR,YR,UMIN,UBUND,UMAX, &
                        LANK,IWAKU)
!-----------------------------------------------------------------------
!  << EQCNTR >>
!         PRESENTED BY KUNIHIKO.WATANABE 14.NOV.1989
!         REVICED   BY HISANORI.TAKAMARU 16.MAR.1990
!-----------------------------------------------------------------------
!     1. FUNCTION
!        (1) TO DRAW TOKOSEN
!     2. ARGUMENTS   (SIZE)   (I/O)     (MEANING)
!        (1) U       NX,NY     (I)       WORLD VALUE
!        (2) W       NX,NY     (I)       WORK ARRAY (real*4)
!        (3) XL,XR,YL,YR       (I)       ABSOLUTE COORDINATE VALUE
!        (4) UMIN,UMAX         (I)       HEIGHT OF MAX & MIN
!                                        UMIN>UMAX : AUTOMATIC CONTROL
!        (5) UBUND             (I)       DRAW DASH LINE (U < UBUND)
!        (6) LANK              (I)       NUMBER OF DRAW LINES
!        (7) IWAKU             (I)       =1 : DRAW FRAME
!     3. CALLED BY
!             (** NOTHING **)
!     4. CALLS
!             (** PLOT   **)
!-----------------------------------------------------------------------
      DIMENSION U(1),W(1)
!
      IF (NX.LT.2) RETURN
      IF (NY.LT.2) RETURN
      IF (XR.LT.XL) RETURN
      IF (YR.LT.YL) RETURN
!
      NXY = NX*NY
      NXM1 = NX - 1
      NYM1 = NY - 1
!
      DX = (XR-XL)/ FLOAT(NXM1)
      DY = (YR-YL)/ FLOAT(NYM1)
!
      UMAX1 = UMAX
      UMIN1 = UMIN
!
      IF(UMAX1.GT.(1.000001*UMIN1)) THEN
!
        DO 10 I = 1 , NXY
          W(I) = U(I) - UMIN1
          IF(U(I).GT.UMAX1) W(I) = UMAX1 - UMIN1
          IF(U(I).LT.UMIN1) W(I) = 0.
   10   CONTINUE
!
      ELSE
!
        UMAX1=-1.E+30
        UMIN1= 1.E+30
        DO 20 I = 1 , NXY
          UMAX1=AMAX1(UMAX1,U(I))
          UMIN1=AMIN1(UMIN1,U(I))
   20   CONTINUE
        DO 25 I = 1 , NXY
          W(I) = U(I) - UMIN1
   25   CONTINUE
!
      ENDIF
!
!------------------------------------------------
      IF(UMAX1.LE.(1.000001*UMIN1))  RETURN
!------------------------------------------------
!
      IF(IWAKU.EQ.1) THEN
        CALL PLOT(XL,YL,3)
        CALL PLOT(XR,YL,2)
        CALL PLOT(XR,YR,2)
        CALL PLOT(XL,YR,2)
        CALL PLOT(XL,YL,2)
        CALL PLOT(XL,YL,3)
      ENDIF
!
      ULD = FLOAT(LANK+1) / (UMAX1-UMIN1)
      EPS = 1.0E-8
!
      NXYM1 = NXM1*NYM1
      DO 9000  IJNXY1 = 1,NXYM1
        J = (IJNXY1-1)/NXM1 + 1
        I = IJNXY1 - (J-1)*NXM1
!
          I1 = I + NX * (J - 1)
          I2 = I1 + 1
          I3 = I1 + 1 + NX
          I4 = I1 + NX
!
          U1 =  W(I1) * ULD
          U2 =  W(I2) * ULD
          U3 =  W(I3) * ULD
          U4 =  W(I4) * ULD
!
          K1 = IFIX(U1)
          K2 = IFIX(U2)
          K3 = IFIX(U3)
          K4 = IFIX(U4)
!
          J1 = IABS(K2-K1)
          J2 = IABS(K3-K2)
          J3 = IABS(K4-K3)
!
          IF(J1.NE.0) THEN
            DO 1000 LL = 1 , J1
              U0 = FLOAT(LL) + FLOAT(MIN0(K1,K2))
                UJOUGE = U0/ULD + UMIN1
                IF (UJOUGE.LT.UBUND) THEN
                  JOUGE = 4
                ELSE
                  JOUGE = 1
                END IF
!
              IF(ABS(U2-U1).LT.EPS)                 GO TO 1000
!
              X1 = XL + DX * ( (U0-U1)/(U2-U1) + FLOAT(I-1) )
              Y1 = YL + DY * FLOAT(J-1)
!
              IF( ((U3-U0)*(U2-U0)).GT.0. )         GO TO 1100
              IF( ( (U0-U2).GT.0. ).AND.( (U0-U4).GT.0. ) ) GO TO 1100
              IF( ABS(U3-U2).LT.EPS )               GO TO 1100
!
                X2 = XL + DX * FLOAT(I)
                Y2 = YL + DY * ( (U0-U2)/(U3-U2) + FLOAT(J-1) )
!
                CALL WDASH(X1,Y1,X2,Y2,JOUGE)
!
 1100         CONTINUE
              IF( ((U4-U0)*(U3-U0)).GT.0. )         GO TO 1200
              IF( ((U1-U0)*(U3-U0)).GT.0. )         GO TO 1200
              IF( ((U2-U0)*(U4-U0)).GT.0. )         GO TO 1200
              IF( ABS(U4-U3).LT.EPS )               GO TO 1200
!
                X2 = XL + DX * ( (U0-U4)/(U3-U4) + FLOAT(I-1) )
                Y2 = YL + DY * FLOAT(J)
                CALL WDASH(X1,Y1,X2,Y2,JOUGE)
!
 1200         CONTINUE
              IF( ((U1-U0)*(U4-U0)).GT.0. )         GO TO 1300
              IF( ( (U0-U1).GT.0. ).AND.( (U0-U3).GT.0. ) ) GO TO 1300
              IF( ABS(U1-U4).LT.EPS )               GO TO 1300
!
                X2 = XL + DX * FLOAT(I-1)
                Y2 = YL + DY*((U0-U1)/(U4-U1)+FLOAT(J-1))
                CALL WDASH(X1,Y1,X2,Y2,JOUGE)
 1300         CONTINUE
 1000       CONTINUE
!
          ENDIF
!
          IF(J2.NE.0) THEN
!
            DO 2000 LL = 1 , J2
              U0 = FLOAT(LL) + FLOAT(MIN0(K2,K3))
                UJOUGE = U0/ULD + UMIN1
                IF (UJOUGE.LT.UBUND) THEN
                  JOUGE = 4
                ELSE
                  JOUGE = 1
                END IF
              IF( ABS(U3-U2).LT.EPS )               GO TO 2000
!
              X1 = XL + DX * FLOAT(I)
              Y1 = YL + DY * ( (U0-U2)/(U3-U2) + FLOAT(J-1) )
!
              IF( ((U4-U0)*(U3-U0)).GT.0. )         GO TO 2100
              IF( ( (U0-U1).GT.0. ).AND.( (U0-U3).GT.0. ) ) GO TO 2100
              IF( ABS(U4-U3).LT.EPS )               GO TO 2100
!
                X2 = XL + DX * ( (U0-U4)/(U3-U4) + FLOAT(I-1) )
                Y2 = YL + DY * FLOAT(J)
!
                CALL WDASH(X1,Y1,X2,Y2,JOUGE)
!
 2100         CONTINUE
              IF( ((U1-U0)*(U4-U0)).GT.0. )         GO TO 2200
              IF( ((U1-U0)*(U3-U0)).GT.0. )         GO TO 2200
              IF( ((U2-U0)*(U4-U0)).GT.0. )         GO TO 2200
              IF( ABS(U1-U4).LT.EPS )               GO TO 2200
!
                X2 = XL + DX * FLOAT(I-1)
                Y2 = YL + DY * ( (U0-U1)/(U4-U1)+FLOAT(J-1) )
                CALL WDASH(X1,Y1,X2,Y2,JOUGE)
 2200         CONTINUE
 2000       CONTINUE
!
          ENDIF
!
          IF(J3.NE.0) THEN
!
            DO 3000 LL = 1 , J3
              U0 = FLOAT(LL) + FLOAT(MIN0(K3,K4))
                UJOUGE = U0/ULD + UMIN1
                IF (UJOUGE.LT.UBUND) THEN
                  JOUGE = 4
                ELSE
                  JOUGE = 1
                END IF
              IF( ABS(U4-U3).LT.EPS )               GO TO 3000
!
              X1 = XL + DX * ( (U0-U4)/(U3-U4) + FLOAT(I-1) )
              Y1 = YL + DY * FLOAT(J)
!
              IF( ((U1-U0)*(U4-U0)).GT.0. )         GO TO 3100
              IF( ( (U0-U2).GT.0. ).AND.( (U0-U4).GT.0. ) ) GO TO 3100
              IF( ABS(U1-U4).LT.EPS )               GO TO 3100
!
                X2 = XL + DX * FLOAT(I-1)
                Y2 = YL + DY * ( (U0-U1)/(U4-U1) + FLOAT(J-1) )
                CALL WDASH(X1,Y1,X2,Y2,JOUGE)
 3100         CONTINUE
 3000       CONTINUE
          ENDIF
 9000 CONTINUE
!
      RETURN
      END
!
!
!-----------------------------------------------------------------------
      subroutine SETSCL (WMINX,WMINY,WMAXX,WMAXY, XL,YL,XR,YR, GDX,GDY,&
                 N1,CHAR1, N2,CHAR2, HIGHT1,                           &
                 NNX,CHARX, HIGHTX, NNY,CHARY, HIGHTY, IWAKU)
!-----------------------------------------------------------------------
!  << SETSCL >>                   /CHAR1/
!                          WMAXY  +--------------------+  (XL,YL)
!                               Y |             (XR,YR)|  (XR,YR) ON 0
!                               R |                    |
!                               A |                    |
!    (WMINX,WMINY)              H |                    |
!    (WMAXX,WMAXY) ON IS        C |                    |
!                                 |(XL,YL)             |
!                          WMINY  +--------+--+--------+
!                                 WMINX  /CHARX/       WMAXX
!-----------------------------------------------------------------------
!
!     SETSCL
!
!     1. FUNCTION
!        (1) TO SCALE THE GRAPHICS BY CALCOMP SPECIFICATIONS
!     2. ARGUMENTS            (I/O)     (MEANING)
!        (1) WMINX,WMAXX,
!            WMINY,WMAXY       (I)       WORLD COORDINATE VALUE
!        (2) XL,XR,YL,YR       (I)       ABSOLUTE COORDINATE VALUE
!        (3) GDX,GDY           (O)       SCALING FACTOR OF COORDINATE
!                                        FROM WORLD TO ABSOLUTE
!        (4) CHAR1,CHARX,CAHRY (I)       TITLE ON GRAPH,X-AXIS,Y-AXIS
!        (5) IWAKU             (I)       DRAW FRAME (0:OFF ; 1:ON)
!                                         999 : WRITE OUT ONLY TITLE,
!                                                    NOT DRAW OTHERWISE
!     3. CALLED BY
!             (** NOTHING **)
!     4. CALLS
!             (** PLOT   **)
!             (** SYMBOL **)
!             (** NUMBER **)
!-----------------------------------------------------------------------
      character*1  CHAR1(1),CHAR2(1),CHARX(1),CHARY(1)
!
      IF (WMAXX.LE.WMINY) GOTO 9999
      IF (WMAXX.LE.WMINY) GOTO 9999
      IF (XR.LE.XL)       GOTO 9999
      IF (YR.LE.YL)       GOTO 9999
!
      GDX= (XR-XL)/(WMAXX-WMINX)
      GDY= (YR-YL)/(WMAXY-WMINY)
!
      XC = 0.5*( XR + XL )
      YC = 0.5*( YR + YL )
!
      IF (N1 .GT.0) THEN
        IF (HIGHT1.GT.0) THEN
          XS1= XC -0.5*N1*HIGHT1
          XS2= XS1 +(N1+1)*HIGHT1
          CALL SYMBOL(XS1,YR+0.1,HIGHT1,CHAR1(1),0.,N1)
          CALL SYMBOL(XS2,YR+0.1,HIGHT1,CHAR2(1),0.,N2)
        END IF
      END IF
!-----------------------------------------------------------------------
      IF (IWAKU.EQ.999) RETURN
!-----------------------------------------------------------------------
!
      IF (IWAKU.EQ.1) THEN
        CALL PLOT (XL,YL,3)
        CALL PLOT (XL,YR,2)
        CALL PLOT (XR,YR,2)
        CALL PLOT (XR,YL,2)
        CALL PLOT (XL,YL,2)
        CALL PLOT (999.,999.0,3)
      END IF
!
      IF (NNX.GT.0) THEN
        IF (HIGHTX.GT.0) THEN
          CALL SYMBOL(XC-0.5*HIGHTX*NNX,YL-0.5,HIGHTX,CHARX(1),0.,1)
          DO NNX1=2,NNX
          CALL SYMBOL(999.0,999.0,HIGHTX,CHARX(NNX1),0.,1)
          end do
        END IF
      END IF
      IF (NNY.GT.0) THEN
        IF (HIGHTY.GT.0) THEN
          CALL SYMBOL(XL-0.5,YC-0.5*HIGHTY*NNY,HIGHTY,CHARY(1),0.,1)
          DO NNY1=2,NNY
          CALL SYMBOL(999.0,999.0,HIGHTY,CHARY(NNY1),0.,1)
          end do
        END IF
      ELSE IF(NNY.LT.0) THEN
        IF (HIGHTY.GT.0) THEN
          CALL SYMBOL(XC-0.5*HIGHTY*NNY,YC,HIGHTY,CHARY(1),0.,1)
          DO NNY1=2,NNY
          CALL SYMBOL(999.0,999.0,HIGHTY,CHARY(NNY1),0.,1)
          end do
        END IF
      END IF
!
      RETURN
!
!-----------------------------------------------------------------------
!
 9999 CONTINUE
      WRITE(6,*) '**********  ABNORMAL WORLD COORDINATE ********'
      WRITE(6,*) '      '
      WRITE(6,*) '    WMAXX =',WMAXX,' WMINX = ',WMINX
      WRITE(6,*) '    WMAXY =',WMAXY,' WMINY = ',WMINY
      WRITE(6,*) '    XL,YL,XR,YR =',XL,YL,XR,YR
      WRITE(6,*) '    FCTR  =',FCTR
      WRITE(6,*) '      '
      CALL CHART
      CALL SYMBOL(1.0,10.0,0.2,' ABNORMAL WORLD COORDINATE CALL',0.,31)
      CALL SYMBOL(1.0,09.0,0.2,' WMAXX =',0.,8)
      CALL NUMBER(999.0,999.0,0.2,WMAXX,0.,2)
      CALL SYMBOL(1.0,08.5,0.2,' WMINX =',0.,8)
      CALL NUMBER(999.0,999.0,0.2,WMINY,0.,2)
      CALL SYMBOL(1.0,08.0,0.2,' WMAXY =',0.,8)
      CALL NUMBER(999.0,999.0,0.2,WMAXY,0.,2)
      CALL SYMBOL(1.0,07.5,0.2,' WMINY =',0.,8)
      CALL NUMBER(999.0,999.0,0.2,WMINY,0.,2)
      CALL SYMBOL(1.0,07.0,0.2,' FCTR  =',0.,8)
      CALL NUMBER(999.0,999.0,0.2,FCTR,0.,2)
      CALL SYMBOL(1.0,06.5,0.2,' XLEFT =',0.,8)
      CALL NUMBER(999.0,999.0,0.2,XL,0.,2)
      CALL SYMBOL(1.0,06.0,0.2,' YLEFT =',0.,8)
      CALL NUMBER(999.0,999.0,0.2,YL,0.,2)
      CALL SYMBOL(1.0,05.5,0.2,' XRIGHT=',0.,8)
      CALL NUMBER(999.0,999.0,0.2,XR,0.,2)
      CALL SYMBOL(1.0,05.0,0.2,' YRIGHT=',0.,8)
      CALL NUMBER(999.0,999.0,0.2,YR,0.,2)
      CALL EXIT
      RETURN
      END
!
!
!-----------------------------------------------------------------------
      subroutine SCALEX (XCM,YCM,X00,Y00,DX,DY,ISC)
!-----------------------------------------------------------------------
      common/GSCALE/ X0(10),Y0(10),XL(10),YL(10),DXI(10),DYI(10)
!
      X0(ISC)= X00
      Y0(ISC)= Y00
      DXI(ISC)= 1./DX
      DYI(ISC)= 1./DY
!
      XL(ISC)= XCM
      YL(ISC)= YCM
!
      RETURN
      END
!
!
!-----------------------------------------------------------------------
      subroutine PLOTL (X,Y,ISC,IPL)
!-----------------------------------------------------------------------
      common/GSCALE/ X0(10),Y0(10),XL(10),YL(10),DXI(10),DYI(10)
!
      XCM= XL(ISC) +DXI(ISC)*(X -X0(ISC))
      YCM= YL(ISC) +DYI(ISC)*(Y -Y0(ISC))
!
      CALL PLOT (XCM,YCM,IPL)
!
      RETURN
      END
!
!
!-----------------------------------------------------------------------
      subroutine VALUES (X,Y,HEIGHT,VAL,THETA,IFMAT)
!-----------------------------------------------------------------------
!  << VALUES >>
!     1. FUNCTION
!        (1) TO DRAW VARIABLE
!     2. ARGUMENTS   (SIZE)   (I/O)     (MEANING)
!        (1) X,Y               (I)       ABSOLUTE COORDINATE VALUE
!        (2) HEIGHT            (I)       DRAW OUT SIZE ON PAPER
!        (3) VAL               (I)       VARIABLE
!        (4) THETA             (I)       ANGLE
!        (5) IFMAT             (I)       FORMAT TYPE
!     3. CALLED BY
!             (** NOTHING **)
!     4. CALLS
!             (** NUMBER **)
!             (** SYMBOL **)
!-----------------------------------------------------------------------
!        IFMAT = (N100)*100 + KETA
!        N100 = 0 : integer FORMAT
!        N100 = 1 : F FORMAT ::  NUMBER(X,Y,HEIGHT,VAL,THETA,KETA)
!        N100 = 2 : E FORMAT ::
!        N100 = 3 : POWER OF TEN FORMAT
!        N100 = OTHEWISE : NOT WRITE OUT
!-----------------------------------------------------------------------
      use, intrinsic :: iso_c_binding
!
      real(C_float) VAL
      character CHR13*13,CHR12*12,CHR3*3
      character(len=1) MINUS,ZERO,BLANK
      parameter(RATIO = 6./7. )
      data MINUS/'-'/,ZERO/'0'/,BLANK/' '/
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      IF (IFMAT.LT.0) RETURN
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      N100 = IFMAT/100
      KETA = IFMAT - N100*100
!
      IF (N100.EQ.0) THEN
        CALL NUMBER(X,Y,HEIGHT,VAL,THETA,-1)
      ELSE IF (N100.EQ.1) THEN
        CALL NUMBER(X,Y,HEIGHT,VAL,THETA,KETA)
      ELSE IF (N100.EQ.2) THEN
        CHR13 = '             '
        CHR12 = '            '
        IF (KETA.EQ.0) THEN
          WRITE(CHR13,'(1PE13.6)') VAL
          CHR12(1:4) = CHR13(1:3)//'E'
          NUMSYM = 4
        ELSE
          KETA = KETA + 1
          IF (VAL.LT.0.) THEN
            CHRVAL = VAL - 5.*10**FLOAT(-KETA)
            WRITE(CHR13,'(1PE13.6)') CHRVAL
            CHR12(1:KETA+3) = CHR13(1:KETA+2)//'E'
            NUMSYM = KETA + 3
          ELSE IF (VAL.EQ.0) THEN
            CHRVAL = VAL
            WRITE(CHR13,'(1PE13.6)') CHRVAL
            CHR12(1:KETA+3) = CHR13(1:KETA+2)//'E'
            NUMSYM = KETA + 3
          ELSE
            CHRVAL = VAL + 5.*10**FLOAT(-KETA)
            WRITE(CHR13,'(1PE13.6)') CHRVAL
            CHR12(1:KETA+2) = CHR13(2:KETA+2)//'E'
            NUMSYM = KETA + 2
          END IF
        END IF
        CHR3 = '   '
!
        IF (CHR13(11:11) .EQ. MINUS) THEN
          IF (CHR13(12:12) .EQ. ZERO  .OR. &
              CHR13(12:12) .EQ. BLANK) THEN
            CHR3(1:2) = '-'//CHR13(13:13)
          ELSE
            CHR3(1:3) = '-'//CHR13(12:13)
          END IF
          NUMSY1 = 3
        ELSE
          IF (CHR13(12:12) .EQ. ZERO  .OR. &
              CHR13(12:12) .EQ. BLANK) THEN
            CHR3(1:1) = CHR13(13:13)
            NUMSY1 = 1
          ELSE
            CHR3(1:2) = CHR13(12:13)
            NUMSY1 = 2
          END IF
        END IF
        AKAKU = 2. * 3.1415927 / 360.
        COST = COS(THETA*AKAKU)
        CALL SYMBOL(X,Y,HEIGHT,CHR12,THETA,NUMSYM)
        CALL SYMBOL(999.,999.,HEIGHT,CHR3,THETA,NUMSY1)
      ELSE IF (N100.EQ.3) THEN
        CHR13 = '             '
        CHR12 = '            '
        IF (KETA.EQ.0) THEN
          WRITE(CHR13,'(1PE13.6)') VAL
          CHR12(1:6) = CHR13(1:3)//'X10'
          NUMSYM = 6
        ELSE
          KETA = KETA + 1
          IF (VAL.LT.0.) THEN
            CHRVAL = VAL - 5.*10**FLOAT(-KETA)
            WRITE(CHR13,'(1PE13.6)') CHRVAL
            CHR12(1:KETA+5) = CHR13(1:KETA+2)//'X10'
            NUMSYM = KETA + 5
          ELSE
            CHRVAL = VAL + 5.*10**FLOAT(-KETA)
            WRITE(CHR13,'(1PE13.6)') CHRVAL
            CHR12(1:KETA+4) = CHR13(2:KETA+2)//'X10'
            NUMSYM = KETA + 4
          END IF
        END IF
        CHR3 = '   '
!
        IF (CHR13(11:11) .EQ. MINUS) THEN
          IF (CHR13(12:12) .EQ. ZERO  .OR. &
              CHR13(12:12) .EQ. BLANK) THEN
            CHR3(1:2) = '-'//CHR13(13:13)
          ELSE
            CHR3(1:3) = '-'//CHR13(12:13)
          END IF
          NUMSY1 = 3
        ELSE
          IF (CHR13(12:12) .EQ. ZERO  .OR. &
              CHR13(12:12) .EQ. BLANK) THEN
            CHR3(1:1) = CHR13(13:13)
            NUMSY1 = 1
          ELSE
            CHR3(1:2) = CHR13(12:13)
            NUMSY1 = 2
          END IF
        END IF
        AKAKU = 2. * 3.1415927 / 360.
        COST = COS(THETA*AKAKU)
        SINT = SIN(THETA*AKAKU)
        CALL SYMBOL(X,Y,HEIGHT,CHR12,THETA,NUMSYM)
!
!                                             *******************
!                                             ** EXPONENT PART **
!                                             *******************
!
        H2 = HEIGHT * 5./7.
        X1 = (NUMSYM+1)* HEIGHT * RATIO
        Y1 = HEIGHT * 4./7.
        IF (ABS(THETA).LT.1E-04) THEN
          X1 = X + X1
          Y1 = Y + Y1
        ELSE
          X2 =     X1 * COST - Y1 * SINT
          Y1 = Y + X1 * SINT + Y1 * COST + H2*COST
          X1 = X + X2                    - H2*SINT
        END IF
        CALL SYMBOL(X1,Y1,H2,CHR3,THETA,NUMSY1)
      END IF
      RETURN
      END
!
!
!-----------------------------------------------------------------------
      subroutine WDASH (X1,Y1,X2,Y2,IPEN )
!-----------------------------------------------------------------------
!  << WDASH  >>                      VER 2.00   16.MAR.1990
!
!     1. FUNCTION
!        (1) TO DRAW LINE FROM (X1,Y1) TO (X2,Y2) BY WDASH
!                            IN ABSOLUTE COORDINATE
!     2. ARGUMENTS            (I/O)     (MEANING)
!        (1) X1,X2,Y1,Y2       (I)       ABSOLUTE COORDINATE VALUE
!        (2) IPEN              (I)       PEN TYPE OF 'WDASH'
!     3. CALLED BY
!             (** EQCNTR  **)
!             (** WDASHL  **)
!     4. CALLS
!             (** PLOT   **)
!-----------------------------------------------------------------------
!       IPEN : MEANING           - : 0.05 (CM)
!        1   :       LINE     -------------------
!        2   :  DASH LINE     --- --- --- --- ---
!        3   :  DASH LINE     -- -- -- -- -- -- --
!        4   :  DASH LINE     - - - - - - - - - -
!        5   :  1 POINT DASH  ---- - ---- - ---- -
!        6   :  2 POINT DASH  --2.0-- - - --2.0--
!   OTHERWISE:  LINE          ---------------------
!-----------------------------------------------------------------------
!
      H1  =  0.05
      H2  =  2.0 * H1
      H3  =  3.0 * H1
      H4  =  4.0 * H1
      H20 = 20.0 * H1
      CALL PLOT ( X1 , Y1 , 3 )
      K = - 1
      IF(IPEN.LT.2) THEN
        GO TO 999
      ELSE IF(IPEN.EQ.2) THEN
        HH1 = H3
        HH2 = H1
      ELSE IF (IPEN.EQ.3) THEN
        HH1 = H2
        HH2 = H1
      ELSE IF (IPEN.EQ.4) THEN
        HH1 = H1
        HH2 = H1
      ELSE IF (IPEN.EQ.5) THEN
        HH1 = H4
        HH2 = H1
        HH3 = H1
        HH4 = H1
      ELSE IF (IPEN.EQ.6) THEN
        HH1 = H20
        HH2 = H1
        HH3 = H1
        HH4 = H1
        HH5 = H1
        HH6 = H1
      END IF
      IF(IPEN.LT.5) THEN
        RLENG = SQRT ( ( X2 - X1 ) **2 + ( Y2 - Y1 ) **2 )
        IF(RLENG.LT.1.0E-5) GOTO 999
        IF(RLENG.LT.HH1) GOTO 999
        COSTH = ( X2 - X1 ) / RLENG
        SINTH = ( Y2 - Y1 ) / RLENG
        D = HH1
        X = X1 + D * COSTH
        Y = Y1 + D * SINTH
        CALL PLOT ( X , Y , ( 5 + K ) / 2 )
        K = - K
        D = D + HH2
        HHH = HH1
        HH1 = HH2
        HH2 = HHH
  200   IF(D.LE.RLENG) THEN
          X = X1 + D * COSTH
          Y = Y1 + D * SINTH
          CALL PLOT ( X , Y , ( 5 + K ) / 2 )
          K = - K
          HHH = HH1
          HH1 = HH2
          HH2 = HHH
          D=D+HH1
          GOTO 200
        END IF
      ELSE IF (IPEN.EQ.5) THEN
        RLENG = SQRT ( ( X2 - X1 ) **2 + ( Y2 - Y1 ) **2 )
        IF(RLENG.LT.1.0E-5) GOTO 999
        IF(RLENG.LT.HH1) GOTO 999
        COSTH = ( X2 - X1 ) / RLENG
        SINTH = ( Y2 - Y1 ) / RLENG
        D = HH1
        X = X1 + D * COSTH
        Y = Y1 + D * SINTH
        CALL PLOT ( X , Y , ( 5 + K ) / 2 )
        K = - K
        D = D + HH2
        HHH = HH1
        HH1 = HH2
        HH2 = HH3
        HH3 = HH4
        HH4 = HHH
  500   IF(D.LE.RLENG) THEN
          X = X1 + D * COSTH
          Y = Y1 + D * SINTH
          CALL PLOT ( X , Y , ( 5 + K ) / 2 )
          K = - K
          HHH = HH1
          HH1 = HH2
          HH2 = HH3
          HH3 = HH4
          HH4 = HHH
          D=D+HH1
          GOTO 500
        END IF
      ELSE IF (IPEN.EQ.6) THEN
        RLENG = SQRT ( ( X2 - X1 ) **2 + ( Y2 - Y1 ) **2 )
        IF(RLENG.LT.1.0E-5) GOTO 999
        IF(RLENG.LT.HH1) GOTO 999
        COSTH = ( X2 - X1 ) / RLENG
        SINTH = ( Y2 - Y1 ) / RLENG
        D = HH1
        X = X1 + D * COSTH
        Y = Y1 + D * SINTH
        CALL PLOT ( X , Y , ( 5 + K ) / 2 )
        K = - K
        D = D + HH2
        HHH = HH1
        HH1 = HH2
        HH2 = HH3
        HH3 = HH4
        HH4 = HH5
        HH5 = HH6
        HH6 = HHH
  600   IF(D.LE.RLENG) THEN
          X = X1 + D * COSTH
          Y = Y1 + D * SINTH
          CALL PLOT ( X , Y , ( 5 + K ) / 2 )
          K = - K
          HHH = HH1
          HH1 = HH2
          HH2 = HH3
          HH3 = HH4
          HH4 = HH5
          HH5 = HH6
          HH6 = HHH
          D=D+HH1
          GOTO 600
        END IF
      END IF
  999 CALL PLOT ( X2 , Y2 , ( 5 + K ) / 2 )
      CALL PLOT ( X2 , Y2 , 3)
      RETURN
      END
!
!
!-----------------------------------------------------------------------
      subroutine DAISHO(X  ,NX,XMIN1,XMAX1)
!-----------------------------------------------------------------------
      DIMENSION X(1)
!
      XMAX1= X(1)
      XMIN1= X(1)
      DO 100 I=2,NX
      XMAX1= AMAX1(XMAX1,X(I) )
      XMIN1= AMIN1(XMIN1,X(I) )
  100 CONTINUE
      RETURN
      END
!
!
!***************************************************************
!*     This program package generates a UNIX postscript        *
!*     graphic file when called by calcomp-compatible          *
!*     /plot23.f/.                                             *
!***************************************************************
!----------------------------------------------------------
!      PostScript header by fortran
!        T. Ogino (Nagoya University) February 27, 1992
!      Modified to conform GSIPP commands
!        Motohiko Tanaka (NIFS)       November 23, 1993
!
!----------------------------------------------- 5/27/1996 -----
!   This PS-Adobe-2.0 header allows us full paging features in
!   the Ghostview. To scroll up the page (backward), click the 
!   page number and press two buttons of mouse simultaneously.
!
!   Consult: A.Saitou (Kyoto U.)  The definition of "@eop"  
!   needs stroke for line drawings (not in the TeX header).
!---------------------------------------------------------------
       subroutine gopen (nframe)
!---------------------------------------------------------------
       common/convsn/ fmag,x0,y0,h0,n0
       common/pages/  ipage,nfrm
!
!*  This is an Adobe-2.0 postscript file.
!
       write(77,10)
   10  format('%!PS-Adobe-2.0',/      &
              '%%Pages: (atend)',/    &
              '%%PageOrder: Ascend',/ &
              '%%EndComments',/       &
              '%%BeginDocument')
!
!%%%%%%%%%%%%%%%%%%% Procedure Defintions %%%%%%%%%%%%%%%%%%%%%%%%%%
!
!     write(77,11) 
!  11 format('%%BoundingBox: 150. 400. 550. 600.')
!
      write(77,21) 
   21 format('/l {lineto} bind def  % x y l -- line to position',/ &
             '/m {moveto} bind def  % x y m -- move to position')
!
      write(77,23) 
   23 format('/tr {/Times-Roman findfont} bind def',/ &
             '/sf {scalefont} bind def',/             &
             '/se {setfont} bind def',/               &
             '/ro {rotate}  bind def',/               &
             '/tl {translate} bind def',/             &
             '/sc {scale} bind def')
!
      write(77,24) 
   24 format('/@bop          % @bop -- begin the a new page',/ &
             '{erasepage newpath initgraphics',/ &
             '/SaveImage save def',/             &
             '} bind def')
!
      write(77,25) 
   25 format('/@eop          % @eop -- end a page',/ &
             '{stroke showpage',/   &
             ' SaveImage restore',/ &
             '} bind def')
!
      write(77,26) 
   26 format('/@end          % @end -- done the whole shebang',/ &
             ' /end load def')
!
      write(77,27) 
   27 format('/dir 0 def')
!
      write(77,29) 
   29 format('/s             % string s -- show the string',/ &
             '{dir 1 eq',/                                    &
             ' {gsave currentpoint translate 90 rotate 0 0 moveto',/&
             ' show grestore}',/ &
             ' {show} ifelse',/  &
             '} bind def')
!
      write(77,31)
   31 format('%%EndDocument',/ &
             '%%EndProlog',/   &
             '%%BeginSetup',/  &
             '/Resolution 300 def',/ &
             '/#copies 1 def',/ &
             '%%EndSetup')
!
!%%%%%%%%%%%%%%%%%%% End of the header %%%%%%%%%%%%%%%%%%%%%%%%%%
!
!*  initiate the page one.
!
       nfrm = nframe
!
       ipage = 1
       write(77,12) ipage,ipage
   12  format('%%Page:',1x,i2,1x,i2)
!
       write(77,30) 
   30  format('%%BeginPageSetup',/ &
              '%%EndPageSetup',/   &
              '@bop')
!
!
!*  Set magnifying factor (GSIPP to Sun coordinate).
!   Rotate and translate to output on A4-L paper.
!      Left corner ...... (  0.,  0.)
!      Right corner ..... (600.,780.)
!
       xcm=  25.
       xwc= 700.
       fmag= xwc/xcm
!
       write(77,*) '90.0 ro'
       write(77,*) '50.0 -550.0 tl'
!
!*  If nfrm=4, four frames in a page (top-left frame).
!
       if(nfrm.eq.1) then
          write(77,*) '1.00 1.00 sc'
       else
          write(77,*) '0.50 0.50 sc'
          write(77,*) '0.0 550.0 tl'
       end if
!
       return
       end
!
!
!-----------------------------
       subroutine gclose
!-----------------------------
       call plote
       return
       end
!
!
!-----------------------------
       subroutine plote
!-----------------------------
       write(77,10) 
   10  format('@eop')
       return
       end
!
!
!-----------------------------------------
       subroutine chart
!-----------------------------------------
!*     Four frames in a page (if nfrm=4).
       common/pages/ ipage,nfrm
!
!
       ipage = ipage +1
       loc= mod(ipage-1,nfrm)
!
!*  Frame 1: open a new page.
!
       if(loc.eq.0) then
          call plote
!
          if(nfrm.eq.1) lpage= ipage
          if(nfrm.ne.1) lpage= (ipage+3)/4
!
          write(77,10) 
   10     format('%%PageTrailer    % Need for the page count')
!
          write(77,20) lpage,lpage
   20     format('%%Page:',1x,i2,1x,i2)
!
          write(77,30) 
   30     format('%%BeginPageSetup',/ &
                 '%%EndPageSetup',/   &
                 '@bop')
!
          write(77,*) '90.0 ro'
          write(77,*) '50.0 -550.0 tl'
!
          if(nfrm.eq.1) then
             write(77,*) '1.00 1.00 sc'
          else
             write(77,*) '0.50 0.50 sc'
             write(77,*) '0.0  550.0 tl'
          end if
!
          return
       end if
!
!-----------------------------------------------------
!      First cancel the previous translation, then
!      make a new translation (scale factor alive).
!-----------------------------------------------------
!*   Frames 2-4:
!
       if(loc.eq.1) then
          write(77,*) '  0.0 -550.0 tl'
          write(77,*) '700.0  550.0 tl'
       end if
!
       if(loc.eq.2) then
          write(77,*) '-700.0 -550.0 tl'
          write(77,*) '   0.0    0.0 tl'
       end if
!
       if(loc.eq.3) then
          write(77,*) '  0.0 0.0 tl'
          write(77,*) '700.0 0.0 tl'
       end if
!
       return
       end
!
!
!------------------------------------
       subroutine factor(fct)
!------------------------------------
       write(77,10) fct,fct
   10  format(f6.2,1x,f6.2,' sc')
       return
       end
!
!
!------------------------------------
       subroutine newpen (ip)
!------------------------------------
       i1=(ip-1)/2
       i2=ip-2*i1
       write(77,*) 'sn'
       pi1=0.40*float(i1-1)
       write(77,30) pi1
   30  format(f3.1,' sl')
       if(i2.ne.1) then
       write(77,*) '[2 2] 0 sd'
       endif
       return
       end
!
!
!---------------------------------------
       subroutine newcolor (ic,r,g,b)
!---------------------------------------
!  ic= 3 tri-color
!  ic= 0 gray scale, r= 0. for black
!
       write(77,*) 'stroke'
!
       if(ic.eq.0) then
         write(77,10) 1.-r  !! (0,1.,0.,0.) for black
   10    format(f4.1,' setgray')
       end if
!
       if(ic.eq.3) then
         write(77,30) r,g,b
   30    format(3f4.1,' setrgbcolor')
       end if
!
       return
       end subroutine newcolor
!
!
!-----------------------------
       subroutine linee
!-----------------------------
       write(77,*) 'st'
       return
       end
!
!
!------------------------------------
       subroutine plot (x0,y0,ip)
!------------------------------------
!
       x= x0
       y= y0
       h= 0.
       n= 777
       call sunscl (x,y,h,n)
!
       if(ip.eq.3)  write(77,10) x,y
       if(ip.eq.2)  write(77,20) x,y
       if(ip.eq.-3) write(77,30) x,y
       if(ip.eq.-2) write(77,40) x,y,x,y
   10  format(f5.1,1x,f5.1,' m')
   20  format(f5.1,1x,f5.1,' l')
   30  format(f5.1,1x,f5.1,' tl')
   40  format(f5.1,1x,f5.1,' l sn',1x,f5.1,1x,f5.1,' tl')
!       write(77,*) 'st'
       return
       end
!
!
!-------------------------------------------------
       subroutine symbol (x0,y0,h0,isymb,ang,n0)
!-------------------------------------------------
       character ica*80,ich(80)*1
       character(*) isymb
       equivalence (ica,ich(1))
!
       x= x0
       y= y0
       h= h0
       n= n0
       call sunscl (x,y,h,n)
!
       write(77,*) 'tr'
       write(77,10) h
   10  format(f5.1,' sf')
       write(77,*) 'se'
       write(77,20) x,y
   20  format(f5.1,1x,f5.1,' m')
       write(77,30) ang
   30  format(f5.1,' ro')
!*
       ica= isymb
       write(77,*) '(',(ich(i),i=1,n),') s'
!
       return
       end
!
!
!-----------------------------------------------
       subroutine number (x0,y0,h0,anu,ang,n0)
!-----------------------------------------------
       character  isymb*9
!
       x= x0
       y= y0
       h= h0
       n= 777
       call sunscl (x,y,h,n)
!
       write(77,*) 'tr'
       write(77,10) h
   10  format(f5.1,' sf')
       write(77,*) 'se'
!
       write(77,20) x,y
   20  format(f5.1,1x,f5.1,' m')
       write(77,30) ang
   30  format(f5.1,' ro')
!
       write(isymb,40) anu
   40  format(1pe9.2)
       write(77,*) '(',isymb,') s'
!
       return
       end
!
!
!---------------------------------------------------
       subroutine sunscl (x,y,h,n)
!---------------------------------------------------
       common/convsn/ fmag,x0,y0,h0,n0
!
       if(x.eq.999.) then
         x= x0 +iabs(n0)*h0
       else
         x= fmag*x
         x0= x
       end if
!
       if(y.eq.999.) then
         y= y0
       else
         y= fmag*y
         y0= y
       end if
!
       h= fmag*h
       h0= h
       if(n.ne.777) n0= n
!
       return
       end
