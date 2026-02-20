!************************************************** 2025/7/03 ****
!*                                                               *
!*   ## Molecular Dynamics of Electrostatic Living Cells ##      *
!*     @nanoporAPG.f03 - Short-range Coulomb and LJ forces and   *
!*     long-range Poisson forces for electrostatic living cells  *
!*                                                               *
!*   Author: Motohiko Tanaka, Ph.D.                              *
!*         Nature and Science Applications, Nagoya 464, Japan.   *
!*   This code is permitted by GNU General Public License v3.0.  *
!*                                                               *
!*     The short-range Coulomb forces of /moldyn/, L.1220, and   *
!*   the long-range electrostatic effects by Poisson equation,   *
!*   L.1275, are treated in this simulation code.                *
!*                                                               *
!*     This molecular dynamics simulation of nanoscale pores     *
!*   in 3D electrostatic effects has been updated in April 2025. *
!*                                                               *
!*   References:                                                 *
!*   1) Y.Rabin and M.Tanaka, DNA in nanopores: Counterion cond- * 
!*    ensation and..., Phys.Rev.Letters, vol.94, 148103 (2005).  *
!*   2) M.Tanaka, https://github.com/Mtanaka77/Electrostatic_    *
!*    molecular_dynamics_of_living_human_cells, 2025, USA.       *
!*                                                               *
!*---------------------------------------------------------------*
!*  >> Note: isize must be chosen such that the sub-box size     *
!*            Lx/isize > ag(counter)+ag(water)                   *
!*                                                               *
!*     ch(1),  am(1),  ag(1)                                     *
!*     charge  mass    radius                                    *
!*                                                               *
!*   Restriction:  p_max = 3                                     *
!*                                                               *
!*     x----x-----x-----x    sub-boxes can be used if            *
!*     |    |   q |     |     1) L_s > rcut_pme (avoid mishit)   *
!*     x----x-----x-----x     2) L >= 4 * L_s   (efficiency)     *
!*     |    | p   |     |                                        *
!*     x----x-----x-----x                                        *
!*     |    |     |     |                                        *
!*     x----x-----x-----x                                        *
!*            L_s                                                *
!*                                                               *
!*  System units:                                                *
!*     length...... a= 1.0 Ang= 1.0 10^-8 cm                     *
!*     mass........ m= 1.67 10^-24 g                             *
!*     time........ t= 0.01 ps= 10^-14 s                         *
!*     charge...... 4.8 10^-10 esu= 1.60 10^-19 C                *
!*          with (1/2)M*(a/tau)**2= kT                           *
!*                                                               *
!*  Equation of motion:                                          *
!*                                                               *
!*      dv    Gamma*q'q grad r                                   *
!*   m ---- = --------- ------ - fgm*(2*r(i)-r(i+1)-r(i-1))      *
!*      dt       r^2      r                                      *
!*                                                               *
!*                      epsLJ       sigma       sigma            *
!*               + 48*----- grad[(-----)^12- (-----)^6] -qE      *
!*                       kT         r_ij        r_ij             *
!*                                                               *
!*  Poisson equation:                                            *
!*   div(eps(i,j,k) [grad pot(i,j,k)]) = - 4*pi*Gamma*rho(i,j,k) *
!*                                                               *
!*   Gamma = Bjerrum/(a*kT) = e**2/(epsLJ*aLJ*kT)                *
!*     The electrostatic coupling constant, Bjerrum=7 at T=300 K *
!*                                                               *
!*****************************************************************
!*  Main program and subroutines:                                *
!*                                                               *
!*   Program nanopore  MPI setup -> setups /Run_MD/ -> /moldyn/  *
!*    param_APG.h (parameter), PORW21_config.start3 (config)     *     
!*                                                               *
!*   /moldyn/     Time cycles, Coulomb and EM fields, L.685-     *
!*   /sht_forces/ Coulomb forces and LJ potential, L.1225, 1815- *
!*   /sprmul/     Spring forces, L.1230, 2290-                   * 
!*   /reflect_endpl/ Particles boundary, L.1455, 2680-           *
!*                                                               *
!*   /init/       Setups from /RUN_MD/, L.385,3320-              *
!*   /poissn/     Poisson equation, L.1300, 5080-                *
!*   /escof3/     ES forces, closed boundary, L.5150,5210-       *
!*     /bound_s/    for it > 1, L.5550                           * 
!*   /cresmd/-/avmult/  Conjugate residual method, L.5170,6320-  *
!*    Graphics    /gopen/ (Adobe-2.0 postscript)                 *
!*                                                               *
!*   First version : 2004/10/25                                  *
!*   Second version: 2006/12/18 (Fortran 95)                     *
!*   Third version : 2025/06/23 (Fortran 2003)                   *
!*                                                               *
!*****************************************************************
!*                                                               *
!*  To get a free format of Fortan f90 or f03, convert f77 into  *
!*    :%s/^c/!/  Change '~c' of ALL columns to '!'               *
!*    :%s/^C/!/                                                  *
!*    tr 'A-Z' 'a-z' <@nanopor.f >@nanopor.f03                   *
!*                                                               *
!*  For subroutines, write "use, intrinsic :: iso_c_binding",    *
!*  "implicit none" is recommended for minimizing typoerrors.    *
!*                                                               *
!*  Compilation by Linux (choose -O0, -O1,...):                  *
!*  % mpif90 -mcmodel=medium -fpic -O2 -o a.out @nanoAPGa.f03 -I/opt/fftw3/include -L/opt/fftw3/lib -lfftw3 &> log       *
!*  % mpiexec -n <proc> a.out &  <proc>=6 or more                                       *
!*                                                               *
!*****************************************************************
!
      program nanopore
      use, intrinsic :: iso_c_binding 
!
      include    'mpif.h'
      include    'paramAPGa.h'
!
      integer(C_INT) size,rank,ierror,ipar,cl_first
      real(C_DOUBLE) wtime,walltime1,walltime2
      logical        first_recyc
!
      integer(C_INT) io_pe
      common/sub_proc/ io_pe
!
!  -----------------
!     num_proc = 6  !  <- size of job: paramAPGa.h
!  -----------------
!
      call mpi_init (ierror)
      call mpi_comm_rank (mpi_comm_world,rank,ierror)
      call mpi_comm_size (mpi_comm_world,size,ierror)
!
      ipar = 1 + rank             ! this is my rank number 
!
      io_pe= 0
      if(ipar.eq.1) io_pe= 1
!*
      if(kstart.eq.0) then
        first_recyc = .true.      ! at fresh start
      else
        first_recyc = .false.
      end if
!
      if(io_pe.eq.1) then
        OPEN (unit=11,file=praefixc//'.06'//suffix2,form='formatted')
!
        write(11,'("My process is #",i3," #")') ipar
        close(11)
      end if
!
      cl_first= 1
      call clocks (walltime1,size,cl_first)
! --------------------------------------------------------------
      call Run_MD (rank,ipar,first_recyc,size,cl_first)
! --------------------------------------------------------------
      cl_first= 2
      call clocks (walltime2,size,cl_first)
!
      wtime= walltime2 -walltime1
!
      if(io_pe.eq.1) then
        OPEN (unit=11,file=praefixc//'.06'//suffix2,             &
              status='unknown',position='append',form='formatted')
!
        write(11,'("*ipar, wtime(sec)=",i3,f10.1)') ipar,wtime
        close(11)
      end if
!*
      call mpi_finalize (ierror)
!*
      stop
      end program nanopore
!
!
!--------------------------------------------------------------------
      subroutine Run_MD (rank,ipar,first_recyc,size,cl_first)
!--------------------------------------------------------------------
      use, intrinsic :: iso_c_binding 
      implicit none
!
      include    'mpif.h'
      include    'paramAPGa.h'
!
      integer(C_INT) rank,ipar,ifrgrod,ifrodco,kstart0,nframe, &
                     np,nq,nr,nCLp,npqr,size,cl_first,i
!
      logical        first_recyc,ende
!
      integer(C_INT) io_pe
      common/sub_proc/ io_pe
!
!  Three particle categories of npqr0, npq0 and np0
      real(C_DOUBLE),dimension(npqr0) :: &
                                xg,yg,zg,vx,vy,vz,ch,am,ag,ep, &
                                vxs,vys,vzs
      real(C_DOUBLE),dimension(npq0) :: vxc,vyc,vzc
!
      real(C_DOUBLE) xmax,ymax,zmax,xmin,ymin,zmin,Temp,      &
                     xleng,yleng,zleng,aLJ,wwat,awat,tmax0
!
      integer(C_INT) ifqq
      real(C_DOUBLE) qfrac,Rpore,Hpore,Zci,Zcp,Zcn
      common/cntion/ qfrac,Rpore,Hpore,Zci,Zcp,Zcn,ifqq
!
      real(C_DOUBLE) xgc,ygc,zgc,vxg,vyg,vzg,rod_leng,Rmac,   &
                     Rhelix,q1,q2,q3,q4,angx,angy,angz,       &
                     ipx,ipy,ipz,tau_rot,gdump
      integer(C_INT) n_rodp,ifbase
!
      common/macroion/ xgc,ygc,zgc,vxg,vyg,vzg,rod_leng,Rmac, &
                       Rhelix,n_rodp
      common/macroiov/ q1,q2,q3,q4,angx,angy,angz,ipx,ipy,ipz
      common/macrotat/ tau_rot,gdump
!
      real(C_DOUBLE),dimension(np0) :: &
                       dox,doy,doz,doxc,doyc,dozc
      common/macroin1/ dox,doy,doz,doxc,doyc,dozc
!
      real(C_DOUBLE) t8,cptot
      common/headr3/ t8,cptot
!*
!  Single precision for the diagnostic routines.
!
      real(C_float)  ecr1,ekin1,elj1,pot1
      real(C_float)  t,t00,xp_leng
      character*8    label,cdate*10,ctime
      common/headr1/ label,cdate
      common/headr2/ t,t00,xp_leng
!
      integer(C_INT) it,is
      real(C_DOUBLE) pi,dt,axi,Gamma,rbmax,vth,tmax
      real(C_float) phi,tht,dtwr1,dtwr2,dtwr3
      common/parm1/ it,is
      common/parm2/ pi,dt,axi,Gamma,rbmax,vth,tmax
      common/parm3/ xmax,ymax,zmax,xmin,ymin,zmin
      common/parm4/ phi,tht,dtwr1,dtwr2,dtwr3
      common/parm8/ xleng,yleng,zleng
!
      integer(C_INT) nsg,nseg,istop
      common/psegm/  nsg(30),nseg
      common/abterm/ istop
!
      integer(C_INT) iwa,iwb,iwc,iwd
      common/imemo/  iwa,iwb,iwc,iwd
      common/shrgf/ aLJ
      common/waterm/ wwat,awat
!                            double real*8
      real(C_DOUBLE)  t_pe,t_poisn
      common/timings/ t_pe,t_poisn
!
      real(C_DOUBLE) Vtop0,Vbot0,Vtop,Vbot,t_recyc,rpr
      integer(C_INT) fixedbc_x,fixedbc_y,itermax,filtx,filty,filtz,ipr
      common/espot/  Vtop0,Vbot0,Vtop,Vbot,fixedbc_x,fixedbc_y, &
                     itermax,filtx,filty,filtz
      common/cresm3/ ipr(10),rpr(10)
!
      real(C_DOUBLE)  rcut_clf,rcutlj,r_fene,k_fene, &
                      r_fene2,Bjerrum
      common /confdatar/ rcut_clf,rcutlj,r_fene
      common /fenepara/  k_fene,r_fene2
      common /elsta/     Bjerrum
!
      integer(C_INT) np00,nnb,ist1,ist2
      common/pbase/  np00,nnb,ist1(100),ist2(100)
!
      real(C_DOUBLE) diel2,dielpr,aa
      common/dielec/ diel2,dielpr,aa
      common/ewald2/ nCLp     ! ip0, paramAPGa.h L.4920
!
      integer(C_INT) ifLJ
      real(C_DOUBLE) epsLJ,eps_K,eps_Cl
      common/parmlj/ epsLJ,eps_K,eps_Cl
      common/parmlj2/ ifLJ
!
      real(C_float),dimension(ntmax) :: &
                    ekin,ppot,ekn2,etot,z_pe,vzpe,vzco,vzct, &
                    vdtm,vpor,time,ecr,elj,espr
      common/ehist/ ekin,ppot,ekn2,etot,z_pe,vzpe,vzco,vzct, &
                    vdtm,vpor,time,ecr,elj,espr
!
      namelist/inp1/ phi,tht
      namelist/inp2/ dt,rbmax,aLJ,tmax,dtwr1,dtwr2,dtwr3
!
      label='+Recycle'
      call date_and_time_7 (cdate,ctime)
!
!**************************************************************
      xp_leng= 7.
      nframe= 4
!
      call FLOPEN (nframe)
!
!     nCLp= npq0   ! np + nq; defined in /init/, L.4920
!     npqr= npqr0  ! total number of particles, READ_CONF in L.3785
!**************************************************************
!--------------------------------------------------------------
!* Read basic data and constants in /READ_CONF/
!    >> Bjerrum, axi, 
!    >> dt, xmax, rcut_clf, rcutlj
!    >> r_fene 
!    >> kstart, cptot, dtwr1, dtwr2, dtwr3
!    >> ifqq, qfrac, Rpore, Hpore, Zcp, Zcn.
!
!   cptot (in min.)
!     cptot= 60.0  <--- 1 h
!     dtwr1=   10.
!     dtwr2=  500.
!     dtwr3=   50.
!     kstart= -1    if particle data read in /READ_CONF/.
!
!    non-neutral PE................ ifqq= 1
!                                        qfrac = dq/q
!    Valence of counterions........ Zcp, Zcn
!      salt (Zcp, -Zcn) after sufficient counterions.
!*********************
!     ifqq= 1
!     qfrac= 0.30
!     Zcp=   1.
!     Zcn=  -1.
!*********************
!*  rbmax: maximum bond length for /sprmul/.
!
      rbmax= 1.5
      aLJ  = 1.4d0
!
!  Water particles.
! ++++++++++++++++++++++++++++++++++++++++++++
      wwat = 18.d0         ! Water mass
      awat = 1.4d0         !  radius
! ++++++++++++++++++++++++++++++++++++++++++++
!
!*  np,nq,nseg :  defined in /READ_CONF/.
!   Vtop0,Vbot0,rpr(10) in READ_CONF
!
      call READ_CONF (np,nq,npqr,ifbase)
!     +++++++++++++++++++++++++++++++++++++++++
!
!* Bjerrum = Gamma aLJ*kT = e**2/epsLJ  
      Temp= 1.d0   ! T= 300 K
      Gamma= Bjerrum/(awat*Temp) 
!
      istop = 0    ! signal for termination: istop= 1
!-----------------------------------------------------
!****************************************
!*  Prepare for graphic output.         *
!****************************************
!*  for 3-D plot of particles /pplt3d/.
!
      phi= -60.
      tht=  15.
!
      pi = 4.d0*atan(1.d0)
!
      if(io_pe.eq.1) then
        OPEN (unit=11,file=praefixc//'.06'//suffix2,             &
              status='unknown',position='append',form='formatted')
!
        write(11,600) label,cdate,ctime
!
        write(11,*) '# np, nq=',np,nq
        write(11,*) '&inp1  phi,tht.....'
        write(11,inp1)
!
        write(11,*) '&inp2  rbmax .....'
        write(11,inp2)
!*
        write(11,610)   t_pe,t_poisn,Vtop0,Vbot0
!
!*  System size (-l/2, l/2)
        write(11,*) '# Lx/2, Ly/2, Lz/2 = ',xmax,ymax,zmax
!
  600   format(/,'<< Nanopore problem (isolated system) >> ', &
               a8,/,' date=',a10,'  time=',a8,/)
  610   format(' %Timings: t_pe    =',f7.1,/, &
               '           t_poisn =',f7.1,/  &
               ' %Potential bias Vtop0, Vbot0=',2f8.3,/)
!
        close(11)
      end if
!-----------------------------------------------------
!************************************
!*   Step 1 : Molecular dynamics.   *
!************************************
!*   <(1/2)*mv**2>= (3/2)t, na**3=1.
!--------------------------------------------
      vth= sqrt(epsLJ)
!--------------------------------------------
! 
! Create the unity Gauss distribution 
!   Result is used -> /gaus1/ fv(101) -> dgaus2(vmax)
!
      call ggauss
!
      call init (xg,yg,zg,vx,vy,vz,ch,am,ag,ep,ifrgrod,ifrodco, &
                 ifbase,ipar,np,nq,nCLp,nr,npqr)
!                             !  np,nq are given at /READ_CONF/
!
!  ------------------------------------------
      call Initialisierung                 !<- in nCLp used
      call Interpol_charge_assign_function
!  ------------------------------------------
!
      if(kstart.le.0) then
        if(io_pe.eq.1) then
          OPEN (unit=11,file=praefixc//'.06'//suffix2,             &
                status='unknown',position='append',form='formatted')
!
          write(11,*) '# Polyelectrolyte #'
          write(11,620) (ch(i),i=1,np)
!
          if(nq.ne.0) then
             write(11,*) '# cations and anions #'
             write(11,620) (ch(i),i=np+1,np+nq)
  620        format(10f5.1)
          end if
!
          close(11)
        end if
      else
!
!* Restart data: kstart= 1.
!    these data overwrite definitions made above.
!
        if(io_pe.eq.1) then
          OPEN (unit=11,file=praefixc//'.06'//suffix2,             &
                status='unknown',position='append',form='formatted')
!
          write(11,630) 
  630     format('Restart data are loaded from Ft12.....',/)
!
          close(11)
        end if
!
!------------------------------------------
!* FT12 must be mounted on all nodes
!------------------------------------------
        open (unit=12,file=praefixi//'.12'//suffix1,status='old', &
                                                form='unformatted')
!                *******                  ****  do not alter.
        read(12) kstart0,it,np,nq,nCLp,nr,npqr
        read(12) ifrgrod,ifrodco,ifbase
! 
        read(12) t8,xg,yg,zg,vxs,vys,vzs,ch,am,ag,ep
        read(12) vxc,vyc,vzc
!
        read(12) pi,dt,Gamma,rbmax,vth,tmax0
        read(12) t,phi,tht,dtwr1,dtwr2,dtwr3,is
        read(12) ekin,ppot,ekn2,etot,z_pe,vzpe,vzco,vzct, &
                 vdtm,vpor,time,ecr,elj,espr
!       read(12) qzd,qzu,cjz,qtop,qpor,qbot
!       iead(12) z0,cj1,cj2
!
        read(12) iwa,iwb,iwc,iwd
        read(12) nsg,nseg
!       read(12) nbox,list
!
        read(12) np00,nnb,ist1,ist2
        read(12) diel2,aa
!                **** **** **** **** 
!       read(12) xmax,ymax,zmax,xmin,ymin,zmin
        read(12) xleng,yleng,zleng
        read(12) aLJ,epsLJ,eps_K,eps_Cl,Vtop0,Vbot0,ifLJ
        read(12) qfrac,Rpore,Hpore,Zci,Zcp,Zcn,ifqq
!
!   If there are rods
        if(ifrgrod.eq.1) then
          read(12) xgc,ygc,zgc,vxg,vyg,vzg,rod_leng,Rmac, &
                   dox,doy,doz,Rhelix,doxc,doyc,dozc,n_rodp
          read(12) q1,q2,q3,q4,angx,angy,angz,ipx,ipy,ipz
          read(12) tau_rot,gdump
        end if
!
        close(12)
!
  101   continue
! ------------------------
!       dtwr1=  10.  by /READ_CONF/
!       dtwr2= 500.
!       dtwr3=  50.
!
        if(first_recyc) then
          t8= 0
          t= 0
!
          is= 0
          it= 0
          iwa= -1
          iwb= -1
          iwc= -1
          iwd= -1
        end if
!
!    nCLp= np +nq  !<-- defined in /init/ L.4920
! ------------------------
        if(io_pe.eq.1) then
          OPEN (unit=11,file=praefixc//'.06'//suffix2,             &
                status='unknown',position='append',form='formatted')
!
          write(11,771) t8
          write(11,772) diel2,dielpr,aa,t00
  771     format(' t8=',f7.1)
  772     format(' diel2(pore),dielpr(water)=',2f7.4,/,'  aa,t00=',2f7.1)
!
          write(11,773) np,nq,nCLp,npqr
  773     format(' np,nq,nCLp,npqr=',4i5)
!
          close(11)
        end if
!
!       call mpi_broadcast ()
!
        if(io_pe.eq.1) then
          OPEN (unit=11,file=praefixc//'.06'//suffix2,             &
                status='unknown',position='append',form='formatted')
!
          write(11,*) 'Present time t, is=',t,is
          close(11)
        end if
      end if
!
!************************************
!*  Step 2 : Molecular dynamics.    *
!************************************
!
      if(io_pe.eq.1) then
        OPEN (unit=11,file=praefixc//'.06'//suffix2,             &
              status='unknown',position='append',form='formatted')
!
        write(11,640) Gamma,awat
  640   format('* Gamma, awat=',2f7.3)
        close(11)
      end if
!
      if(io_pe.eq.1) then
        OPEN (unit=11,file=praefixe//'.06'//suffix2,             &
              status='unknown',position='append',form='formatted')
!
        write(11,*)
        write(11,*) '--- Subroutine WRITE_CONF ------------------------'
        close(11)
!
        OPEN (unit=07,file=praefixe//'.06'//suffix2,             &
              status='unknown',position='append',form='formatted')
!
        ende = .true.
        call WRITE_CONF (ende,ifbase,np,nq,npqr)
!
        close (07)
!
        OPEN (unit=11,file=praefixe//'.06'//suffix2,             &
              status='unknown',position='append',form='formatted')
!
        write(11,*) '--- Subroutine WRITE_CONF ------------------------'
        write(11,*)
        close(11)
!
      end if
!
!**************************************************************
!*  Step 3: Subroutine /moldyn/                               *
!**************************************************************
!*-------------------------------------------------------------
      call moldyn (xg,yg,zg,ch,am,ag,ep,vxs,vys,vzs,vxc,vyc,vzc, &
                   ipar,ifrgrod,ifrodco,ifbase,first_recyc,      &
                   np,nq,nCLp,nr,npqr,size,cl_first)
!*-------------------------------------------------------------
!
!**************************************************************
!*  Step 4: Restart data.                                     *
!**************************************************************
!
      if(io_pe.eq.1) then
!
        OPEN (unit=12,file=praefixc//'.12'//suffix2,status='unknown', &
                                                    form='unformatted')
!
        write(12) kstart,it,np,nq,nCLp,nr,npqr
        write(12) ifrgrod,ifrodco,ifbase
!
        write(12) t8,xg,yg,zg,vxs,vys,vzs,ch,am,ag,ep
        write(12) vxc,vyc,vzc
!
        tmax0= tmax
        write(12) pi,dt,Gamma,rbmax,vth,tmax0
        write(12) t,phi,tht,dtwr1,dtwr2,dtwr3,is
        write(12) ekin,ppot,ekn2,etot,z_pe,vzpe,vzco,vzct, &
                  vdtm,vpor,time,ecr,elj,espr
!       write(12) qzd,qzu,cjz,qtop,qpor,qbot
!       write(12) z0,cj1,cj2
!
        write(12) iwa,iwb,iwc,iwd
        write(12) nsg,nseg
!       write(12) nbox,list
!
        write(12) np00,nnb,ist1,ist2
        write(12) diel2,aa
!                 **** *** **** **** 
!       write(12) xmax,ymax,zmax,xmin,ymin,zmin
        write(12) xleng,yleng,zleng
        write(12) aLJ,epsLJ,eps_K,eps_Cl,Vtop0,Vbot0,ifLJ
        write(12) qfrac,Rpore,Hpore,Zci,Zcp,Zcn,ifqq
!
        if(ifrgrod.eq.1) then
          write(12) xgc,ygc,zgc,vxg,vyg,vzg,rod_leng,Rmac, &
                    dox,doy,doz,Rhelix,doxc,doyc,dozc,n_rodp
          write(12) q1,q2,q3,q4,angx,angy,angz,ipx,ipy,ipz
          write(12) tau_rot,gdump
        end if
!
        close(12)
!**********************************************************
!
        call averg1 (ekin,ekin1,is)
        call averg1 (ppot, pot1,is)
        call averg1 ( ecr, ecr1,is)
        call averg1 ( elj, elj1,is)
!
        OPEN (unit=11,file=praefixc//'.06'//suffix2,             &
              status='unknown',position='append',form='formatted')
!
        write(11,650) ekin1,pot1
  650   format(/,'@@ average @@',/,       &
               ' e_kin, e_grid= ',1p2e11.3) 
        close(11)
!
        OPEN (unit=77,file=praefixc//'.77'//suffix2//'.ps',      &
              status='unknown',position='append',form='formatted')
!
        call lplots
!
        call plote
        close(77)
      end if
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
      subroutine FLOPEN (nframe)
!----------------------------------------------------------------------
      use, intrinsic :: iso_c_binding 
!
      include    'paramAPGa.h'
!
      integer(C_INT) io_pe
      common/sub_proc/ io_pe
!
!   --------------------------
      if(io_pe.ne.1) return
!   --------------------------
!
!     if(io_pe.eq.1) then
        OPEN (unit=11,file=praefixc//'.06'//suffix2,             &
              status='unknown',position='append',form='formatted')
!
        write(11,*) ' Uses ',praefixs//'_config.start'//suffix0
        close(11)
!     end if
!
      OPEN (unit=77,file=praefixc//'.77'//suffix2//'.ps', &
                                          form='formatted')
!
      call gopen (nframe)
      close (77)
!
!     open (unit=12,file=praefixc//'.12'//suffix2,
!     open (unit=13,file=praefixc//'.13'//suffix2,
!
      return
      end subroutine FLOPEN
!
!
!-----------------------------------------------------------------------
      subroutine moldyn (xg,yg,zg,ch,am,ag,ep,vxs,vys,vzs,vxc,vyc,vzc, &
                         ipar,ifrgrod,ifrodco,ifbase,first_recyc,      &
                         np,nq,nCLp,nr,npqr,size,cl_first)
!-----------------------------------------------------------------------
!*  Double precision.
!
      use, intrinsic :: iso_c_binding 
      implicit none
!
      include    'mpif.h'
      include    'paramAPGa.h'
!
      real(C_DOUBLE),dimension(npqr0) :: &
                                      xg,yg,zg,vx,vy,vz,ch,am,ag,ep, &
                                      vxs,vys,vzs
      real(C_DOUBLE),dimension(npqr0) :: fsx,fsy,fsz
      real(C_DOUBLE),dimension(npq0) ::  fcx,fcy,fcz,fgx,fgy,fgz, &
                                         vxc,vyc,vzc
      real(C_DOUBLE),dimension(np0)  ::  ftx,fty,ftz,fpx,fpy,fpz
!
      integer(C_INT) ipar,np,nq,nCLp,nr,npqr
      integer(C_INT) io_pe,cnt_recv,disp_recv,i0,i2,i3,i4, &
                     cnt_recv2,disp_recv2,size,cl_first
      common/sub_proc/ io_pe
      common/dat_typ0/ i0(30),i2(30),cnt_recv(30),disp_recv(30)
      common/dat_typ2/ i3(30),i4(30),cnt_recv2(30),disp_recv2(30)
!
      real(C_DOUBLE) t8,cptot
      common/headr3/ t8,cptot
!
!  The space cell index uses 0,mx-1
      real*8   t1,t2,t3,d1,xx,yy,zz
      real(C_DOUBLE),dimension(0:ip0**3-1,0:npq0-1) :: ql
      real(C_DOUBLE),dimension(0:mx-1,0:my-1,0:mz-1) :: &
                                 rho,ex,ey,ez,dec2,pot,ez1
      real(C_DOUBLE) diel
      common/potsav/ pot,dec2
      logical     first_recyc 
!
      integer(C_INT) ierr
      real(C_DOUBLE) sres,avex,rsdl,conv
      common/crconv/ sres,avex,rsdl,conv,ierr
!
      integer(C_INT) i,j,k,m,m0,xpos,ypos,zpos,          &
                     g(0:npq0-1,0:2),mx1,mx2,my1,my2,mz1,mz2
!
      real(C_float),dimension(npio) :: x4,y4,z4,vx4,vy4,vz4,  &
                                       ch4,am4,ag4 
      real(C_float) t,t00,xp_leng
      character*8   label,cdate*10
      common/headr1/ label,cdate
      common/headr2/ t,t00,xp_leng
!
      real(C_float),dimension(ntmax) :: &
                    ekin,ppot,ekn2,etot,z_pe,vzpe,vzco,vzct, &
                    vdtm,vpor,time,ecr,elj,espr
      common/ehist/ ekin,ppot,ekn2,etot,z_pe,vzpe,vzco,vzct, &
                    vdtm,vpor,time,ecr,elj,espr
!
      real(C_float) vpr1,vpr2,rr
      integer(C_INT) npor,kp1,kp2
!
      integer(C_INT) it,is,iwa,iwb,iwc,iwd,iwrt1,iwrt2,iwrt3,iwrt4, &
                     nsg,nseg,istop,ifqq,iwrta,iwrtb,iwrtc,iwrtd,   &
                     ndim,iterp,ncti,ncoi,ntimes,istep
!
      integer(C_INT) np00,nnb,ist1,ist2,nc_DNA
      common/pbase/ np00,nnb,ist1(100),ist2(100)
      real(C_DOUBLE) diel2,dielpr,aa
      common/dielec/ diel2,dielpr,aa
!
      real(C_float) phi,tht,dtwr1,dtwr2,dtwr3,                 &
                    qfrac4,Rpore4,Hpore4,Zcp4,Zcn4,Bjerrum4,  &
                    xmax4,ymax4,zmax4,walltime,               &
                    dtm,vm,s0,s1,s2,vsq,ekin0,ekin1,          &
                    ekin2,sz1,svz1,svz2,svz3,ef2,ef20
      real(C_DOUBLE) pi,dt,axi0,Gamma,rbmax,vth,tmax,         &
                    xmax,ymax,zmax,xmin,ymin,zmin,symp,symp2, &
                    xleng,yleng,zleng,dgaus2,vmax0,vmax1,vmax2, &
                    e_c_s,e_grid,e_c_r,e_lj,e_elas,q0,ss,dt0, &
                    fv,vv,dv,Temp
      common/gaus1/ fv(101),vv,dv
!
      common/parm1/ it,is
      common/parm2/ pi,dt,axi0,Gamma,rbmax,vth,tmax
      common/parm3/ xmax,ymax,zmax,xmin,ymin,zmin
      common/parm4/ phi,tht,dtwr1,dtwr2,dtwr3
      common/parm8/ xleng,yleng,zleng
!
      common/imemo/ iwa,iwb,iwc,iwd
      common/iotim/ iwrt1,iwrt2,iwrt3,iwrt4
      common/psegm/ nsg(30),nseg
      common/abterm/ istop
      common/energy/ e_c_s,e_grid,e_c_r,e_lj,e_elas
!
      real(C_DOUBLE)  t_pe,t_poisn
      common/timings/ t_pe,t_poisn
!
      real(C_DOUBLE)  rcut_clf,rcutlj,r_fene,Bjerrum,aLJ,qfrac, &
                      Rpore,Hpore,Zci,Zcp,Zcn,wwat,awat
      common /confdatar/ rcut_clf,rcutlj,r_fene
      common /elsta/ Bjerrum
      common/shrgf/  aLJ
      common/cntion/ qfrac,Rpore,Hpore,Zci,Zcp,Zcn,ifqq
      common/waterm/ wwat,awat
!
      real(C_DOUBLE) gx,gy,gz,ghx,ghy,ghz,ghx2,ghy2,ghz2, &
                     ghxsq,ghysq,ghzsq,hxi,hyi,hzi
      integer(C_INT) ptx,pty,ptz
      common/xtable/ gx(0:mx),ghx(0:mx),ghx2(0:mx),ghxsq(0:mx), &
                     gy(0:my),ghy(0:my),ghy2(0:my),ghysq(0:my), &
                     gz(0:mz),ghz(0:mz),ghz2(0:mz),ghzsq(0:mz), &
                     hxi,hyi,hzi
      common/xtabl2/ ptx(-10:3000),pty(-10:3000),ptz(-10:4000)
!
      integer(C_INT) pxr,pxc,pxl,pyr,pyc,pyl,pzr,pzc,pzl
      common/ptable/ pxr(-1:mx+1),pxc(-1:mx+1),pxl(-1:mx+1), &
                     pyr(-1:my+1),pyc(-1:my+1),pyl(-1:my+1), &
                     pzr(-1:mz+1),pzc(-1:mz+1),pzl(-1:mz+1)
!
      real(C_DOUBLE) Vtop0,Vbot0,Vtop,Vbot,t_recyc,rpr
      integer(C_INT) fixedbc_x,fixedbc_y,itermax,filtx,filty,filtz,ipr
      common/espot/  Vtop0,Vbot0,Vtop,Vbot,fixedbc_x,fixedbc_y, &
                     itermax,filtx,filty,filtz
      common/trecyc/ t_recyc
      common/cresm3/ ipr(10),rpr(10)
!
      real(C_DOUBLE) xgc,ygc,zgc,vxg,vyg,vzg,rod_leng,Rmac,  &
                     Rhelix,q1,q2,q3,q4,angx,angy,angz,      &
                     ipx,ipy,ipz,tau_rot,gdump,ranff
      integer(C_INT) n_rodp,ifrgrod,ifrodco,ifrecyc,ifbase
      common/macroion/ xgc,ygc,zgc,vxg,vyg,vzg,rod_leng,Rmac,  &
                       Rhelix,n_rodp
      common/macroiov/ q1,q2,q3,q4,angx,angy,angz,ipx,ipy,ipz
      common/macrotat/ tau_rot,gdump
!
      real(C_DOUBLE),dimension(np0) :: &
                       dox,doy,doz,doxc,doyc,dozc
      common/macroin1/ dox,doy,doz,doxc,doyc,dozc
!
      real(C_DOUBLE) wtime0,wtime1,wtime2,wtime3,wtime4,wtime5
!
      integer(C_INT) ifLJ
      real(C_float)  epsLJ4
      real(C_DOUBLE) epsLJ,eps_K,eps_Cl
      common/parmlj/ epsLJ,eps_K,eps_Cl
      common/parmlj2/ ifLJ
!
      logical    first_11/.true./,first_15/.true./, &
                 first_ppl/.true./
!
! -------------------------
!*  Paralell task by MPI 
! -------------------------
!
      if(num_proc.gt.8) then
        if(io_pe.eq.1) then
        OPEN (unit=11,file=praefixc//'.06'//suffix2,             &
              status='unknown',position='append',form='formatted')
!
        write(11,*) ' num_proc= 1,6,8 is assumed ! '
        close(11)
        end if
!
        return
      end if
!
!++++++++++++++++++++++++++++
!* 1) Particles  i0(1)= 1   +
!++++++++++++++++++++++++++++
!
      do k= 1,8        
      if(num_proc.eq.1) then
          i0(k)= 1    !* np+1
          i2(k)= npqr
      end if 
!
      if(num_proc.eq.6) then
        if(k.eq.1) then
          i0(k)= 1    !* np+1
          i2(k)= npqr/6 -1
        else if(k.eq.2) then
          i0(k)= npqr/6
          i2(k)= 2*npqr/6 -1
        else if(k.eq.3) then
          i0(k)= 2*npqr/6
          i2(k)= 3*npqr/6 -1
        else if(k.eq.4) then
          i0(k)= 3*npqr/6
          i2(k)= 4*npqr/6 -1
        else if(k.eq.5) then
          i0(k)= 4*npqr/6
          i2(k)= 5*npqr/6 -1
        else if(k.eq.6) then
          i0(k)= 5*npqr/6
          i2(k)= npqr
        end if
      end if
!
      if(num_proc.eq.8) then
        if(k.eq.1) then
          i0(k)= 1    !* np+1
          i2(k)= npqr/8 -1
        else if(k.eq.2) then
          i0(k)= npqr/8
          i2(k)= 2*npqr/8 -1
        else if(k.eq.3) then
          i0(k)= 2*npqr/8
          i2(k)= 3*npqr/8 -1
        else if(k.eq.4) then
          i0(k)= 3*npqr/8
          i2(k)= 4*npqr/8 -1
        else if(k.eq.5) then
          i0(k)= 4*npqr/8
          i2(k)= 5*npqr/8 -1
        else if(k.eq.6) then
          i0(k)= 5*npqr/8
          i2(k)= 6*npqr/8 -1
        else if(k.eq.7) then
          i0(k)= 6*npqr/8
          i2(k)= 7*npqr/8 -1
        else if(k.eq.6) then
          i0(k)= 7*npqr/8
          i2(k)= npqr
        end if
      end if
!
      disp_recv(k)= i0(k) -1         !! = 0 if no displacement 
      cnt_recv(k) = i2(k) -i0(k) +1
      end do
!
!++++++++++++++++++++++++
!* 2) Fields. i3(1)= 0  +
!++++++++++++++++++++++++
!
      do k= 1,8
      if(num_proc.eq.1) then
        i3(k)= 0
        i4(k)= mxyz-1
      end if
!
      if(num_proc.eq.6) then
        if(k.eq.1) then
          i3(k)= 0
          i4(k)= mxyz/6-2
        else if(k.eq.2) then
          i3(k)= mxyz/6-1
          i4(k)= 2*mxyz/6-2
        else if(k.eq.3) then
          i3(k)= 2*mxyz/6-1
          i4(k)= 3*mxyz/6-2
        else if(k.eq.4) then
          i3(k)= 3*mxyz/6-1
          i4(k)= 4*mxyz/6-2
        else if(k.eq.5) then
          i3(k)= 4*mxyz/6-1
          i4(k)= 5*mxyz/6-2
        else if(k.eq.6) then
          i3(k)= 5*mxyz/6-1
          i4(k)= mxyz -1
        end if
      end if
!
      if(num_proc.eq.8) then
        if(k.eq.1) then
          i3(k)= 0
          i4(k)= mxyz/8-2
        else if(k.eq.2) then
          i3(k)= mxyz/8-1
          i4(k)= 2*mxyz/8-2
        else if(k.eq.3) then
          i3(k)= 2*mxyz/8-1
          i4(k)= 3*mxyz/8-2
        else if(k.eq.4) then
          i3(k)= 3*mxyz/8-1
          i4(k)= 4*mxyz/8-2
        else if(k.eq.5) then
          i3(k)= 4*mxyz/8-1
          i4(k)= 5*mxyz/8-2
        else if(k.eq.6) then
          i3(k)= 5*mxyz/8-1
          i4(k)= 6*mxyz/8-2
        else if(k.eq.7) then
          i3(k)= 6*mxyz/8-1
          i4(k)= 7*mxyz/8-2
        else if(k.eq.8) then
          i3(k)= 7*mxyz/8-1
          i4(k)= mxyz -1
        end if
      end if
!
      disp_recv2(k)= i3(k)
      cnt_recv2(k) = i4(k) -i3(k) +1
!
      if(io_pe.eq.1) then
        OPEN (unit=11,file=praefixc//'.06'//suffix2,             &
              status='unknown',position='append',form='formatted')
!
        write(11,117) k,i3(k),i4(k)
  117   format('# k=',i6,'  i3,i4=',2i8)
!
        close(11)
      end if
      end do
!
!--------------------------
!*  Settings
!--------------------------
!
      dt0= dt
      if(kstart.eq.0) then
        t8= -dt 
        t = t8
!
        is= 0
        it= 0
        iwa= -1
        iwb= -1
        iwc= -1
        iwd= -1
!
        istep= 0  !<-- write for the Poisson solver
!
        do i= 1,nCLp
        vxc(i)= 0
        vyc(i)= 0
        vzc(i)= 0
        end do
!
        do i= 1,npqr0
        vxs(i)= 0
        vys(i)= 0
        vzs(i)= 0
        end do
      end if
!
!
      do k= 0,mz-1
      do j= 0,my-1
      do i= 0,mx-1
      xx= gx(i)
      yy= gy(j)
      zz= gz(k)
!
      dec2(i,j,k)= diel(xx,yy,zz)
      end do
      end do
      end do
! 700 continue
!
!-------------------------------------------------------
      if(io_pe.eq.1) then
        qfrac4 = qfrac
        Rpore4 = Rpore
        Hpore4 = Hpore
        Zcp4   = Zcp
        Zcn4   = Zcn
!
        Bjerrum4 = Bjerrum
        epsLJ4   = epsLJ
        xmax4  = xmax
        ymax4  = ymax
        zmax4  = zmax
!
        OPEN (unit=13,file=praefixc//'.13'//suffix2, &
                      status='unknown',form='unformatted')
!
!                            ** ****** ******
        write(13) np,nq,nCLp,nr,ifbase,epsLJ4,Rpore4,Hpore4, &
                  Zcp4,Zcn4,Bjerrum4,xmax4,ymax4,zmax4
!
        do i= 1,npio
        ch4(i)= ch(i)
        am4(i)= am(i)
        ag4(i)= ag(i)
        end do
!
        write(13) ch4,am4,ag4
        close(13)
      end if
!
!***********************************************
!* >> Reset ion velocities at t= t_vreset      *
!********************************* 9/07/2001 ***
!*  Read data in /READ_CONF/
!   Values are suggested:
!     t_pe    =    3. ! PE motion starts at t_pe
!     t_poisn =    1. ! solve Poisson eq after t_poisn
!
!* For ions.
!
!   ww1  = 38.d0  ! K
!   ww2  = 34.d0  ! Cl
!
      vth = sqrt(epsLJ)
!     vmax0= vth/sqrt(218.d0/wwat)  !<-- mass of sugar, 218./wwat
      vmax0= vth/sqrt( 96.d0/wwat)   !<-- mass of PO_4
      vmax1= vth/sqrt( 38.d0/wwat)   !<-- mass of K(+)
      vmax2= vth/sqrt( 18.d0/wwat)   !<-- water
!
      if(kstart.eq.0) then
!
      do i= 1,np
      vxs(i)= dgaus2(vmax0)
      vys(i)= dgaus2(vmax0)
      vzs(i)= dgaus2(vmax0)
      end do
!
      do i= np+1,nCLp
      vxs(i)= dgaus2(vmax1)
      vys(i)= dgaus2(vmax1)
      vzs(i)= dgaus2(vmax1)
      end do
!
      do i= nCLp+1,npqr
      vxs(i)= dgaus2(vmax2)
      vys(i)= dgaus2(vmax2)
      vzs(i)= dgaus2(vmax2)
      end do
!
      do i= 1,nCLp
      fgx(i)= 0.d0
      fgy(i)= 0.d0
      fgz(i)= 0.d0
      end do
!
      end if
!
!***************************
!* Step 1: Initialization  *
!***************************
!
 1000 continue
!    ++++++++++++
      t8= t8 +dt
      t = t8
      it= it +1
!    ++++++++++++
!
      cl_first= 2
      call clocks (wtime0,size,cl_first)
!
!     if(t8.gt.tmax) go to 2000
!     if((wtime0/60.d0).gt.cptot) go to 2000
!
      if(istop.ge.1) then
        if(io_pe.eq.1) then
        OPEN (unit=11,file=praefixc//'.06'//suffix2,             &
              status='unknown',position='append',form='formatted')
!
        write(11,*) ' Uses Abnormal termination (istop=1)... '
        write(11,*) '  ipar,t8=',ipar,t8
        close(11)
        end if
       
        go to 2000
      end if
!
!-------------------------------
!  Define write-out intervals.
!-------------------------------
      iwrt1= iwrta(t,dtwr1)   ! energy data, dtwr1=10.
      iwrt2= iwrtb(t,dtwr2)   ! restart data, dtwr2=200
      iwrt3= iwrtc(t,dtwr3)   ! field data, dtwr3=50
      iwrt4= iwrtd(t,  10.)   ! write(13), 10.
!*
      if(t8.le.t00) then
        Vtop = 0.d0
        Vbot = 0.d0
!
        do k= 0,mz-1
        do j= 0,my-1
        do i= 0,mx-1
        ez1(i,j,k)= 0
        end do
        end do
        end do
!
      else if(t8.gt.t00) then
        Vtop = Vtop0 *(1.d0 -exp(-(t8-t00)/100.d0))
        Vbot = Vbot0 *(1.d0 -exp(-(t8-t00)/100.d0))
!
        do k= 0,mz-1
        do j= 0,my-1
        do i= 0,mx-1
        ez1(i,j,k)=  - (Vtop-Vbot)/zleng  !<- ez1< 0, increse -10
        end do                            !   (-e)ez1> 0
        end do
        end do
      end if
!
      if(t8.le.t_recyc) then
        ifrecyc = 0         ! Equilibration phase
        first_recyc= .false.
!
      else if(t8.gt.t_recyc) then
        ifrecyc = 1         ! Recycle phase
      end if
!
      if(iwrt1.eq.0) then
        is= is +1
        if(is.ge.ntmax) then
          call rehist
        end if
!
      end if
!
      Temp= 1.d0
      Gamma= Bjerrum/(awat*Temp) 
!
!*****************************
!*  Step 2: Velocity update  *
!*****************************
!
      e_lj   = 0.d0
      e_c_r  = 0.d0
!
      do i = 1,np
      fpx(i)= 0.d0
      fpy(i)= 0.d0
      fpz(i)= 0.d0
!
      ftx(i)= 0.d0
      fty(i)= 0.d0
      ftz(i)= 0.d0
      end do
!
      do i = 1,nCLp
      fcx(i)= 0.d0
      fcy(i)= 0.d0
      fcz(i)= 0.d0
      end do
!
      do i = 1,npqr
      fsx(i)= 0.d0
      fsy(i)= 0.d0
      fsz(i)= 0.d0
      end do
!
!* Complete sht_forces, and sprmul
!
        cl_first= 2
        call clocks (wtime1,size,cl_first)
!
      call sht_forces (xg,yg,zg,ch,ag,ep,fcx,fcy,fcz,  &
                       fsx,fsy,fsz,ipar,np,nq,nCLp,npqr)
!
        cl_first= 2
        call clocks (wtime2,size,cl_first)
!
      call sprmul (xg,yg,zg,ch,ag,fpx,fpy,fpz,np)
!
        cl_first= 2
        call clocks (wtime3,size,cl_first)
!
!
!**************************************************
!*  Step 3: Coulomb forces on grids               *
!**************************************************
!   The space cell index follows i= 0,mx-1 as /charges/.
!
!*  -4*pi*rho = div*(eps*grad pot)
!    pot = pot_q + pot_g
!        pot_q: direct sum of coulomb forces
!        pot_g: geometry correction (boundary)
!                     
! ------------------------------------
      if(t8.lt.t_poisn) go to 370 
! ------------------------------------
      ntimes= 10 
      if(mod(it,ntimes).eq.0) then
!**                       +++
        istep= istep +1  
!         increment of time: dt*ntimes= 0.01x10 = 0.1
!
        do i= 1,nCLp
        fgx(i)= 0.d0
        fgy(i)= 0.d0
        fgz(i)= 0.d0
        end do
!
        call charges (rho,xg,yg,zg,ch,ql,g,ipar,nCLp)
!       call water_dens (xg,yg,zg,ipar,nCLp,npqr)
!
        do i= 0,mx-1
        do j= 0,my-1
        do k= 0,mz-1
        rho(i,j,k)= 0
        pot(i,j,k)= 0
        end do
        end do
        end do
!
        do k= 1,mz-2   ! inner points for i= 1,mx-2
        do j= 1,my-2
        do i= 1,mx-2
        rho(i,j,k)= -4.d0*pi*Gamma*rho(i,j,k)/dec2(i,j,k)
        pot(i,j,k)= 1.d-3*(ranff(0.d0) -0.5)  ! random seed
        end do 
        end do
        end do
!
        symp = -1.d0
        symp2=  1.d0
        call filt1p (rho,filtx,filty,filtz,symp,symp2)
!            +++++++++++     
!
!* Abandon the solution if not converged,
!  which occurs when rho() is large locally
!
        ndim= 3
        call poissn (rho,pot,ndim,itermax,iterp,ipar)
!       +++++++++++++++++++++++++++++++++++++++++++++
!
!       symp = -1.d0
!       symp2=  1.d0
!       call filt1p (pot,filtx,filty,filtz,symp,symp2)
!
        e_grid = 0.d0
!
        do k= 1,mz-2    ! inner points i= 1,mx-2
        do j= 1,my-2
        do i= 1,mx-2
        ex(i,j,k)= -(pot(i+1,j,k) -pot(i-1,j,k))/ghx2(i)
        ey(i,j,k)= -(pot(i,j+1,k) -pot(i,j-1,k))/ghy2(j)
        ez(i,j,k)= -(pot(i,j,k+1) -pot(i,j,k-1))/ghz2(k) +ez1(i,j,k)
!
        e_grid = e_grid + &
                    sqrt(ex(i,j,k)**2 +ey(i,j,k)**2 +ez(i,j,k)**2)
        end do
        end do
        end do
!
        e_grid = e_grid/float(mx*my*mz)
!
! ----------------------------------------------------
!         increment of time: dt*ntimes= 0.01x10 = 0.1
        if(mod(istep,50).eq.0 .and. io_pe.eq.1) then   !  istep=50 for Dt=5 
!
          OPEN (unit=11,file=praefixc//'.06'//suffix2,             &
                status='unknown',position='append',form='formatted')
!
          write(11,341) t8,ierr,e_grid,rsdl
  341     format('#(bcg-ps) t8, ierr=',f8.2,i6,'  e_grid,rsdl=',1p2d12.3)
          close(11)
        end if
! ----------------------------------------------------
!
        do i= 1,nCLp  !<-- only ions 
        m0= -1
!
        do j = -1,1   ! symmetric
        do k = -1,1
        do m = -1,1
        xpos = pxc(g(i,0) + j)   ! use fresh g()
        ypos = pyc(g(i,1) + k)
        zpos = pzc(g(i,2) + m)
!    
        m0 = m0 +1
        d1 = ql(m0,i)
!    
        fgx(i) = fgx(i) + d1*ex(xpos,ypos,zpos)
        fgy(i) = fgy(i) + d1*ey(xpos,ypos,zpos)
        fgz(i) = fgz(i) + d1*ez(xpos,ypos,zpos)
        end do
        end do
        end do
        end do
!
!* Coulomb forces are calculated under a large stride.
!
        do i= 1,nCLp
        dtm= ntimes*dt/am(i)
!
        vxc(i)= vxc(i) +fgx(i)*dtm  !<-- Coulomb forces
        vyc(i)= vyc(i) +fgy(i)*dtm  !    once in ntimes steps
        vzc(i)= vzc(i) +fgz(i)*dtm
        end do
!
        if(t8.le.t_pe) then
          do i= 1,np
          vxc(i)= 0.d0
          vyc(i)= 0.d0
          vzc(i)= 0.d0
          end do
        end if
!
      end if
  370 continue  !<-- From L.1480-1620
!
        cl_first= 2
        call clocks (wtime4,size,cl_first)
!
!
!* Forces for Coulomb forces in a short stride.
!
      if(t8.ge.t_pe) then
        do i= 1,np
        dtm= dt/am(i)
!                                  Coulomb and PO_4
        vxs(i)= vxs(i) +(fcx(i) +fpx(i) +fsx(i))*dtm  !<- Macroions 
        vys(i)= vys(i) +(fcy(i) +fpy(i) +fsy(i))*dtm
        vzs(i)= vzs(i) +(fcz(i) +fpz(i) +fsz(i))*dtm
        end do
      end if
!
      do i= np+1,nCLp
      dtm= dt/am(i)
!
      vxs(i)= vxs(i) +(fcx(i) +fsx(i))*dtm  !<- Counter/co-ions
      vys(i)= vys(i) +(fcy(i) +fsy(i))*dtm
      vzs(i)= vzs(i) +(fcz(i) +fsz(i))*dtm
      end do
!
      if(t8.le.t_pe+10.d0) then   ! retain salt until DNA contracts
        do i= np+1,nCLp           !             5/28/2006
        if(abs(zg(i)).lt.0.5*Hpore) then
          vzs(i)= 0.d0
        end if
        end do
      end if
!
      do i= nCLp+1,npqr
      dtm= dt/am(i)
!
      vxs(i)= vxs(i) +fsx(i)*dtm  !<- water (ball) 
      vys(i)= vys(i) +fsy(i)*dtm
      vzs(i)= vzs(i) +fsz(i)*dtm
      end do
!
!  vxc() +vxs()
      do i= 1,nCLp
      vx(i)= vxc(i) +vxs(i)  !<-- Coulomb + LJ
      vy(i)= vyc(i) +vys(i)
      vz(i)= vzc(i) +vzs(i)
      end do
!
      do i= nCLp+1,npqr
      vx(i)= vxs(i)          !<-- LJ
      vy(i)= vys(i)
      vz(i)= vzs(i)
      end do
!
! ------------------------
!* Pin down DNA center
! ------------------------
!     if(np.ne.0) then
!        nc_DNA= np00/2
!
!        vx(nc_DNA)= 0.d0
!        vy(nc_DNA)= 0.d0
!        vz(nc_DNA)= 0.d0
!     end if
!
      do i= 1,npqr
      xg(i)= xg(i) +dt*vx(i)  !<- Update the positions
      yg(i)= yg(i) +dt*vy(i)
      zg(i)= zg(i) +dt*vz(i)
      end do
!
!---------------------
!* Reflection walls
!---------------------
! 
      call reflect_endpl (xg,yg,zg,vx,vy,vz,ch,am,ag,np,nq,npqr)
!
!* After the reflection, separate vx() -vxc() -> vxs()
!
      do i= 1,nCLp
      vxs(i)= vx(i) -vxc(i)  !<- Charged species
      vys(i)= vy(i) -vyc(i)
      vzs(i)= vz(i) -vzc(i)
      end do
!
      do i= nCLp+1,npqr
      vxs(i)= vx(i)          !<- Water (ball)
      vys(i)= vy(i)
      vzs(i)= vz(i)
      end do
!
        cl_first= 2
        call clocks (wtime5,size,cl_first)
!
      if(iwrt1.eq.0 .and. io_pe.eq.1) then
        OPEN (unit=11,file=praefixc//'.06'//suffix2,             &
              status='unknown',position='append',form='formatted')
!
        write(11,612) t8,wtime5-wtime0,wtime1-wtime0,wtime2-wtime1, &
                      wtime3-wtime2,wtime4-wtime3,wtime5-wtime4
  612   format(' t8,wtime=',f6.1,f8.2,' 1-5=',5f7.4)
!*
!       write(11,*)
!       write(11,*) '# t8=',t8
!
!       write(11,*) 'xg - zg, vx - vz...'
!       do i= 1,np+nq
!       write(11,991) i,xg(i),yg(i),zg(i),vx(i),vy(i),vz(i)
! 991   format('i=',i5,3f7.1,2x,1p3d12.3)
!       end do
!
        write(11,*) '##### t8,Vtop,Vbot=',t8,Vtop,Vbot
        close(11)
      end if
!
!*******************************
!*  Step 4: Diagnosis section  *
!*******************************
!------------------------------
! 1. Energy.
!
      if(iwrt1.eq.0 .and. io_pe.eq.1) then
!
        OPEN (unit=11,file=praefixc//'.06'//suffix2,             &
              status='unknown',position='append',form='formatted')
!
! --------------------------------------- On major node --------------
!
        vm= 0
        s0= 0
        s1= 0
        s2= 0
!
        if(ifrgrod.eq.1) then
          ekin0= 0.5*am(1)*(vxg**2 +vyg**2 +vzg**2)
        else
!
          do i= 1,np
          s0= s0 +0.5*am(i)*(vx(i)**2 +vy(i)**2 +vz(i)**2)
          end do
        end if
!
        do i= np+1,nCLp
        vsq= vx(i)**2 +vy(i)**2 +vz(i)**2
        s1= s1 +0.5*am(i)*vsq
        vm= max(sqrt(vsq), vm)
        end do
!
        do i= nCLp+1,npqr
        s2= s2 +0.5*am(i)*(vx(i)**2 +vy(i)**2 +vz(i)**2)
        end do
!
        Ekin0= s0/(np +1.e-5)
        Ekin1= s1/(nq +1.e-5)
        Ekin2= s2/(nr +1.e-5)
!
        time(is)= t
        vdtm(is)= vm*dt
!
        ekin(is) = ekin1
        ekn2(is) = ekin2
        ecr (is) = e_c_r  !/nCLp
        elj (is) = e_lj   !/npqr
        espr(is) = e_elas !/np
        ppot(is) = e_grid !/nCLp
        etot(is) = 0.5d0*( s0 +s1 +s2 ) & ! /npqr 
                   +ecr(is) +elj(is) +espr(is) +ppot(is)
!
        sz1 = 0.
        svz1= 0.
        svz2= 0.
        svz3= 0.
        ncti= 0
        ncoi= 0
!
        do i= 1,np
        svz1= svz1 +vz(i)
        sz1 = sz1  +zg(i)
        end do
!
        do i= np+1,nCLp
        if(ch(i).gt.0.) then
           svz2= svz2 +vz(i)
           ncti= ncti +1
        else
           svz3= svz3 +vz(i)
           ncoi= ncoi +1
        end if
        end do
!
        vzpe(is)= svz1/(np +0.001)
        z_pe(is)=  sz1/(np +0.001)
!
        vzct(is)= svz2/ncti
        vzco(is)= svz3/ncoi
!
        npor= 0
        vpr1= 0.
        vpr2= 0.
        do k= 0,mz-1
        if(gz(k).lt. Hpore/2) kp2= k
        if(gz(k).lt.-Hpore/2) kp1= k+1
        end do
!
        do j= my/2-5,my/2+5
        do i= mx/2-5,mx/2+5
        rr= sqrt(gx(i)**2 +gy(j)**2)
        if(rr.lt.Rpore) then
          npor= npor +1
          vpr1= vpr1 +pot(i,j,kp1)  
          vpr2= vpr2 +pot(i,j,kp2)  
        end if
        end do
        end do
!
        vpor(is) = (vpr2 -vpr1)/npor
        if(it.eq.1.and.io_pe.eq.1) then
          write(11,*) 'Vpore: zp1,zp2,npor=',gz(kp1),gz(kp2),npor
        end if
!
        if(first_11) then
          first_11= .false.
!
          if(io_pe.eq.1) then
          write(11,600)
  600     format(/,'  time:    E_kin0/np  E_kin1/nq  E_kin2/nr  ', &
                   'E_coulb    E_grid     E_lj       E_elas     ', &
                   'E_tot    iterp rsdl      wall(sec)')
          end if
        end if
!
        if(io_pe.eq.1) then
          write(11,610) t8,ekin0,ekin1,ekin2,ecr(is),ppot(is),       &
                        elj(is),espr(is),etot(is),ipr(2),rpr(2),wtime5 
  610     format('t8=',f7.1,1p8e11.3,i5,e9.2,0pf9.1) 
        end if
!
        close(11)
      end if
!
!---------------------------------------------------------------------
      if(iwrt4.eq.0 .and. io_pe.eq.1) then
!
        OPEN (unit=13,file=praefixc//'.13'//suffix2,               &
              status='unknown',position='append',form='unformatted')
!
        do i= 1,nCLp
        x4(i) = xg(i)
        y4(i) = yg(i)
        z4(i) = zg(i)
!
        vx4(i)= vx(i)
        vy4(i)= vy(i)
        vz4(i)= vz(i)
        end do
!
        write(13)  t,x4,y4,z4,vx4,vy4,vz4
        close(13)
      end if
! --------------------------------------- On major node --------------
!
! 2. Particle plots.
! --------------------------- 
!*  Total potential: cplot3
! --------------------------- 
      if(iwrt3.eq.0) then
!
        if(t8.ge.t_poisn) then
!
          call charges (rho,xg,yg,zg,ch,ql,g,ipar,nCLp)
!
          do k= 0,mz-1
          do j= 0,my-1
          do i= 0,mx-1
          rho(i,j,k)= -4.d0*pi*Gamma *rho(i,j,k)/dec2(i,j,k)
          pot(i,j,k)= 1.d-3*( ranff(0.d0) -0.5)  ! random seed
          end do
          end do
          end do
!
!* All PE must call /poissn/ (MPI wait!)
!
          ndim= 3
          call poissn (rho,pot,ndim,itermax,iterp,ipar)
!
          if(io_pe.eq.1) then
!
            OPEN (unit=77,file=praefixc//'.77'//suffix2//'.ps',      &
                  status='unknown',position='append',form='formatted')
!
!           call cplot3 (rho,xmax,ymax,zmax,'rho(all)',8)
            call cplot3 (pot,xmax,ymax,zmax,'pot(all)',8)
!
!!          call cplot3 (dec2,xmax,ymax,zmax,'dielect.',8)
            close(77)
!
!
            OPEN (unit=16,file=praefixc//'.16'//suffix2,             &
                status='unknown',position='append',form='unformatted')
!
            mx1= (mx-1)/2-10
            mx2= (mx-1)/2+10
            my1= (my-1)/2-10
            my2= (my-1)/2+10
            mz1= 0
            mz2= mz-1
!
            write(16) t,mx1,mx2,my1,my2,mz1,mz2
            write(16) (((rho(i,j,k),i= mx1,mx2),j=my1,my2),k=mz1,mz2)
            write(16) (((pot(i,j,k),i= mx1,mx2),j=my1,my2),k=mz1,mz2)
!
            close(16)
          end if
        end if
!
!**
        if(iwrt3.eq.0 .and. io_pe.eq.1) then
!
          OPEN (unit=77,file=praefixc//'.77'//suffix2//'.ps',      &
                status='unknown',position='append',form='formatted')
!
          if(t8.gt.1.d0) then
            call ppl3da (xg,yg,zg,ch,ag,rod_leng,Rmac,np,nq,nCLp,npqr, &
                         first_ppl)
!
            call vdistr (vx,vy,vz,np,nq,nCLp,nr,npqr)
          end if

          close(77)
        end if
      end if
!
! ---------------------------
!*  Dump all particle data.
! ---------------------------
      if(iwrt3.eq.0 .and. io_pe.eq.1) then
!     real(C_float),dimension(npio) :: x4,y4,z4,vx4,vy4,vz4,  &
!                                      ch4,am4,ag4,xp,yp,zp
        do i= 1,npio
        x4(i) =  xg(i)
        y4(i) =  yg(i)
        z4(i) =  zg(i)
        vx4(i)=  vx(i)
        vy4(i)=  vy(i)
        vz4(i)=  vz(i)
!
        ch4(i)=  ch(i)
        am4(i)=  am(i)
        ag4(i)=  ag(i)
        end do
!
        if(io_pe.eq.1) then
          OPEN (unit=15,file=praefixc//'.15'//suffix2, &
                status='unknown',form='unformatted')
!
          if(first_15) then
            first_15= .false.
            write(15)  np,nq,nCLp,nr,qfrac4,Rpore4,Hpore4,  &
                       Zcp4,Zcn4,Bjerrum4,xmax4,ymax4,zmax4
            write(15)  ch4,am4,ag4
          end if
!
          write(15)  t,x4,y4,z4,vx4,vy4,vz4
          close(15)
        end if
      end if
!
! -------------------
!*  Restart data.
! -------------------
!
      if(iwrt2.eq.0 .and. io_pe.eq.1) then
        OPEN (unit=12,file=praefixc//'.12'//suffix2,status='replace', &
                                                    form='unformatted')
!
        write(12) kstart,it,np,nq,nCLp,nr,npqr
        write(12) ifrgrod,ifrodco,ifbase
!
        write(12) t8,xg,yg,zg,vxs,vys,vzs,ch,am,ag,ep
        write(12) vxc,vyc,vzc
!
        write(12) pi,dt,Gamma,rbmax,vth,tmax
        write(12) t,phi,tht,dtwr1,dtwr2,dtwr3,is
        write(12) ekin,ppot,ekn2,etot,z_pe,vzpe,vzco,vzct, &
                  vdtm,vpor,time,ecr,elj,espr
!       write(12) qzd,qzu,cjz,qtop,qpor,qbot
!       write(12) z0,cj1,cj2
!
        write(12) iwa,iwb,iwc,iwd
        write(12) nsg,nseg
!       write(12) nbox,list
!
        write(12) np00,nnb,ist1,ist2
        write(12) diel2,aa
!                 **** *** **** **** 
!       write(12) xmax,ymax,zmax,xmin,ymin,zmin
        write(12) xleng,yleng,zleng
        write(12) aLJ,epsLJ,eps_K,eps_Cl,Vtop0,Vbot0,ifLJ
        write(12) qfrac,Rpore,Hpore,Zci,Zcp,Zcn,ifqq
!
!   If rods are present
        if(ifrgrod.eq.1) then
          write(12) xgc,ygc,zgc,vxg,vyg,vzg,rod_leng,Rmac, &
                    dox,doy,doz,Rhelix,doxc,doyc,dozc,n_rodp
          write(12) q1,q2,q3,q4,angx,angy,angz,ipx,ipy,ipz
          write(12) tau_rot,gdump
        end if
!
        close(12)
      end if
! --------------------------------------- On major node --------------
      if(mod(it,10).eq.0) then
        if(t8.gt.tmax) go to 2000
        if((wtime0/60.d0).gt.cptot) go to 2000
      end if
      go to 1000
!
 2000 continue
      if(io_pe.eq.1) then
        OPEN (unit=11,file=praefixc//'.06'//suffix2,             &
              status='unknown',position='append',form='formatted')
!
        write(11,*) ' Final: t, tmax=',t8,tmax
        write(11,*) ' Final: wtime5/60., cptot=',wtime5/60.,cptot
!
        close(11)
      end if
!
      return
      end subroutine moldyn
!
!
!----------------------------------------------------------------
      subroutine sht_forces (xg,yg,zg,ch,ag,ep,fcx,fcy,fcz,  &
                             fsx,fsy,fsz,ipar,np,nq,nCLp,npqr)
!----------------------------------------------------------------
      use, intrinsic :: iso_c_binding 
      implicit  none
!*
      include   'mpif.h'
      include   'paramAPGa.h'
!
      real(C_DOUBLE),dimension(npqr0) :: xg,yg,zg,ch,ag,ep,   &    
                                         fsx,fsy,fsz,frx,fry,frz
      real(C_DOUBLE),dimension(npq0) ::  fcx,fcy,fcz,fdx,fdy,fdz 
      integer(C_INT) ipar,np,nq,nCLp,npqr,istop
      common/abterm/ istop
!
      real(C_DOUBLE) rcutcl2,rsq 
      real(C_DOUBLE) dx,dy,dz,r2,r,forceV,ff,f_cut,diel
      common/rsqcut/ rcutcl2(npq0),rsq(npq0,npq0)
!
      real(C_DOUBLE) t8,cptot
      common/headr3/ t8,cptot
!
      real(C_float)  time,t00,xp_leng
      common/headr2/ time,t00,xp_leng
!
      real(C_DOUBLE) pi,dt,axi,Gamma,rbmax,vth,tmax,       &
                     rcut_clf2,rcps2,unif1,unif2, &
                     e_c_s,e_c_pme,e_c_r,e_lj,e_elas
!
      integer(C_INT) i,j,k,l,ibox,ibx(npqr0),itabl,neigh,  &
                     ierror,io_pe,ix,iy,iz
      common/sub_proc/ io_pe
!-----------
      real(C_DOUBLE) gx,gy,gz,ghx,ghy,ghz,ghx2,ghy2,ghz2,  &
                     ghxsq,ghysq,ghzsq,hxi,hyi,hzi
      integer(C_INT) ptx,pty,ptz
      common/xtable/ gx(0:mx),ghx(0:mx),ghx2(0:mx),ghxsq(0:mx), &
                     gy(0:my),ghy(0:my),ghy2(0:my),ghysq(0:my), &
                     gz(0:mz),ghz(0:mz),ghz2(0:mz),ghzsq(0:mz), &
                     hxi,hyi,hzi
      common/xtabl2/ ptx(-10:3000),pty(-10:3000),ptz(-10:4000)
!
      integer(C_INT) iwrt1,iwrt2,iwrt3,iwrt4
      common/iotim/ iwrt1,iwrt2,iwrt3,iwrt4
!
      common/parm2/  pi,dt,axi,Gamma,rbmax,vth,tmax
      common /cutoffel/ rcut_clf2,rcps2
      common /energy/   e_c_s,e_c_pme,e_c_r,e_lj,e_elas
!
!
      integer(C_INT),dimension(30) :: i0,i2,cnt_recv,disp_recv
      integer(C_INT) ifrodco,cnt_send,i00,jj
      common/dat_typ0/ i0,i2,cnt_recv,disp_recv
!
      integer(C_INT) isize,isizeZ,isize2,isize4,nc3
      parameter  (isize=13,isizeZ=25,nc3=isize**2*isizeZ)
      parameter  (isize2=isize*isize,isize4=isize2+isize)
!
!                                         nbxs in paramAPGa.h
      integer(C_INT)  ncel,lcel,nipl,lipl,ibind
      common/LJ_list/ ncel(nc3), lcel(nbxs,nc3), &
                      nipl(npqr0),lipl(nbxs,npqr0)
      common /boxind/ ibind(27,nc3)
!
      real(C_DOUBLE) aLJ,asc2
      common/shrgf/  aLJ
!
      real(C_DOUBLE) rcut_clf,rcutlj,r_fene,eps_lj,addpot
      common /confdatar/ rcut_clf,rcutlj,r_fene
!
      real(C_DOUBLE) xmax,ymax,zmax,xmin,ymin,zmin, &
                     xleng,yleng,zleng
      real(C_float) phi,tht,dtwr1,dtwr2,dtwr3
      common/parm3/ xmax,ymax,zmax,xmin,ymin,zmin
      common/parm4/ phi,tht,dtwr1,dtwr2,dtwr3
      common/parm8/ xleng,yleng,zleng
!
      real(C_DOUBLE) driwu,driwu2,rlj,rsi,snt,ccel
!
      integer(C_INT) ifLJ
      real(C_DOUBLE) rljcut,rljcut2
      real(C_DOUBLE) epsLJ,eps_K,eps_Cl
      common/parmlj/ epsLJ,eps_K,eps_Cl
      common/parmlj2/ ifLJ
!
      real(C_DOUBLE) xgc,ygc,zgc,vxg,vyg,vzg,rod_leng,Rmac,Rhelix
      integer(C_INT) n_rodp
      common/macroion/ xgc,ygc,zgc,vxg,vyg,vzg,rod_leng,Rmac,Rhelix, &
                       n_rodp
!
      real(C_DOUBLE),dimension(np0) :: &
                       dox,doy,doz,doxc,doyc,dozc
      common/macroin1/ dox,doy,doz,doxc,doyc,dozc
!*---------------------------------------------------------------
!     rcut_clf= sqrt(rcut_clf2) <- from READ_CONF
!     rcutlj <- common/confdatar/ rcut_clf,rcutlj, READ_CONF
!
      driwu2 = 1.25992104989487316476721060728d0  ! 2**(1/3)
      driwu  = sqrt(driwu2)
!
      asc2 = (0.85d0*aLJ)**2
!
      do i= 1,npqr
      rcutcl2(i)= (ag(i) * rcut_clf)**2 
      end do
!            a_phos: ag(i)= 4.1/(2*1.5)= 1.37
!            rcut_clf= 5 -> rcutcl2(i)= 1.37*5 = 6.83
!
!*---------------------------------------------------------
!*  Update interaction table in every 5 steps. 
!*---------------------------------------------------------
!
      itabl= t8/dt +0.1d0
      if(mod(itabl,5).eq.1) then
!     ++++++++++++++++++++++++++
!
!* Step 1.
! ------------------------------------------------
!* Register all particles (include rod charges).
! ------------------------------------------------
! 
!     -----------------
         call Labels
!     -----------------
!
         do k= 1,nc3
         ncel(k)= 0     ! Clear the (ix,iy,iz) cell registers.
         end do
!
         do i = 1,npqr
         ix= int(isize* (xg(i)-xmin)/xleng +1.5001)
         iy= int(isize* (yg(i)-ymin)/yleng +1.5001)
         iz= int(isizeZ*(zg(i)-zmin)/zleng +1.5001)
         if(ix.le.1 .or. ix.ge.isize ) go to 100
         if(iy.le.1 .or. iy.ge.isize ) go to 100
         if(iz.le.1 .or. iz.ge.isizeZ) go to 100
!
         ibox = ix + isize*(iy-1 + isize*(iz-1))
! 
         ncel(ibox)= ncel(ibox) +1     ! Number of targets is registered
         lcel(ncel(ibox),ibox)= i      ! Large bin is used
!
         if(ncel(ibox).gt.nbxs) then
           if(io_pe.eq.1) then
           OPEN (unit=11,file=praefixc//'.06'//suffix2,             &
                 status='unknown',position='append',form='formatted')
!
           write(11,*) 'full: ncel> nbxs !  i,ibox=',i,ibox
           write(11,930) i,xg(i),yg(i),zg(i)
  930      format('i=',i5,' xg,yg,zg=',1p3d12.3)
!
           close(11)
           end if
!
           stop
         end if
  100    end do
!
!* Step 2.
!
         do i= i0(ipar),i2(ipar)
         nipl(i)= 0
!
         ix= int(isize* (xg(i)-xmin)/xleng +1.5001)
         iy= int(isize* (yg(i)-ymin)/yleng +1.5001)
         iz= int(isizeZ*(zg(i)-zmin)/zleng +1.5001)
         if(ix.le.1 .or. ix.ge.isize ) go to 210
         if(iy.le.1 .or. iy.ge.isize ) go to 210
         if(iz.le.1 .or. iz.ge.isizeZ) go to 210
!
         ibox = ix + isize*(iy-1 + isize*(iz-1))
!
         do k = 1,27
         neigh= ibind(k,ibox)
         if(neigh.le.0) go to 230     ! Far away
!
         do l= 1,ncel(neigh)          ! Find ions in the boxes around i-th
         j= lcel(l,neigh)             ! j-th belongs to this box.
!
         if(i.eq.j) go to 240
! 
         dx= xg(i) -xg(j)  ! No folding
         dy= yg(i) -yg(j)
         dz= zg(i) -zg(j)
!
!***
         r2 = dx**2 +dy**2 +dz**2
         if(r2.lt.rcutcl2(i)) then
!
           nipl(i)= nipl(i) +1 
           lipl(nipl(i),i)= j 
!
         end if
  240    end do
  230    end do
!
!        if(nipl(i).eq.0) go to 210
  210    end do
!
!*  ---------------------------------
         do i = i0(ipar),i2(ipar)
         if(nipl(i).gt.nbxs) then
!
           if(io_pe.eq.1) then
           OPEN (unit=11,file=praefixc//'.06'//suffix2,             &
                 status='unknown',position='append',form='formatted')
!
           write(11,*) 'Stop: f_solvent, @ nipl> nbxs @  ipar=',ipar
           write(11,*) '       i, nipl(i)=',i,nipl(i)
           write(11,*) ' common/LJ_list/ lipl(nbxs,npqr0)'
           close(11)
           end if
!
           istop= 1
           return
         end if
         end do
!*  ---------------------------------
      end if
!
!**********************************
!* Step 4: The Coulomb forces.    *
!**********************************
!*  Round-robin partitioning is best for 'allreduce'
!
      do i= ipar,nCLp,num_proc
      do j= i+1,nCLp
      dx = xg(i) -xg(j)
      dy = yg(i) -yg(j)
      dz = zg(i) -zg(j)
!
      r2 = dx**2 + dy**2 + dz**2
      if(r2.le.rcutcl2(i)) then
!
      r  = sqrt(r2)
!
!             Short-range forces
      forceV = f_cut(r,rcut_clf) *Gamma*ch(i)*ch(j)/r2        &
                *diel((xg(i)+xg(j))/2.d0,(yg(i)+yg(j))/2.d0,  &
                                         (zg(i)+zg(j))/2.d0)
!
      fcx(i) = fcx(i) + forceV*dx/r
      fcy(i) = fcy(i) + forceV*dy/r
      fcz(i) = fcz(i) + forceV*dz/r
!
      fcx(j) = fcx(j) - forceV*dx/r
      fcy(j) = fcy(j) - forceV*dy/r
      fcz(j) = fcz(j) - forceV*dz/r
!
      e_c_r = e_c_r + Gamma*ch(i)*ch(j)/(diel((xg(i)+xg(j))/2.d0, &
                        (yg(i)+yg(j))/2.d0,(zg(i)+zg(j))/2.d0)*r) 
      end if
      end do
      end do
!
! --------------------
!*  Unify the force.
! --------------------
!  fcx() spreads on all domains - must take "allreduce".....
!
      if(num_proc.ne.1) then
!
        call mpi_allreduce (fcx,fdx,nCLp,mpi_real8,mpi_sum, &
                            mpi_comm_world,ierror)
        call mpi_allreduce (fcy,fdy,nCLp,mpi_real8,mpi_sum, &
                            mpi_comm_world,ierror)
        call mpi_allreduce (fcz,fdz,nCLp,mpi_real8,mpi_sum, &
                            mpi_comm_world,ierror)
!
        do i= 1,nCLp
        fcx(i)= fdx(i)
        fcy(i)= fdy(i)
        fcz(i)= fdz(i)
        end do
!
        unif1= e_c_r
        call mpi_allreduce (unif1,unif2,1,mpi_real8,mpi_sum, &
                            mpi_comm_world,ierror)
        e_c_r= unif2
!
        if(iwrt1.eq.0 .and. io_pe.eq.1) then
          OPEN (unit=11,file=praefixc//'.06'//suffix2,             &
                status='unknown',position='append',form='formatted')
          write(11,*)
          write(11,*) 'e_c_r/(np+nq)=',e_c_r/(np+nq)
          close(11)
        end if
      end if
!
!*************************
!* Step 5: LJ potential  *
!*************************
!
      do i = i0(ipar),i2(ipar) 
      do jj= 1,nipl(i) 
      j= lipl(jj,i)
!
!  ----------------------------------------
      rljcut= driwu         ! no dimension, driwu=1.12
!
      if(ifLJ.eq.1) then
        if(i.le.nCLp .or. j.le.nCLp) then
          rljcut= rcutlj    ! no dimension, sigma/Ang
        end if
      end if
!  ----------------------------------------
!
      rljcut2= rljcut**3
      addpot =  1.d0/rljcut2**4 - 1.d0/rljcut2**2

      dx = xg(i) -xg(j)
      dy = yg(i) -yg(j)
      dz = zg(i) -zg(j)
!
!*  Lennard-Jones force.
!   >> Energy unit: epsLJ (LJ energy).
!   >> rlj= 1.0 when i and j are touching.
!
      r2 = dx**2 + dy**2 + dz**2
      rlj = sqrt(r2)/(ag(i)+ag(j))
!
      if(rlj.le.rljcut) then
        rsi = 1.d0/max(rlj**2,asc2)
        snt = rsi*rsi*rsi
!
        eps_lj= sqrt(ep(i)*ep(j))
        ccel  = 48.d0*eps_lj*snt*(snt-0.5d0)/r2
!
        fsx(i) = fsx(i) + ccel*dx
        fsy(i) = fsy(i) + ccel*dy
        fsz(i) = fsz(i) + ccel*dz
!
        e_lj = e_lj  + 4.d0*eps_lj*(snt*(snt - 1.d0) - addpot)
      end if
!
      end do
      end do
!
! --------------------
!*  Unify the force.
! --------------------
!* Define packed arrays for mpi_allgatherv.
!
      if(num_proc.ne.1) then
!
        call mpi_allreduce (fsx,frx,npqr0,mpi_real8,mpi_sum, &
                            mpi_comm_world,ierror)
        call mpi_allreduce (fsy,fry,npqr0,mpi_real8,mpi_sum, &
                            mpi_comm_world,ierror)
        call mpi_allreduce (fsz,frz,npqr0,mpi_real8,mpi_sum, &
                            mpi_comm_world,ierror)
!
        do i= 1,npqr
        fsx(i)= frx(i)
        fsy(i)= fry(i)
        fsz(i)= frz(i)
        end do
!
!
        unif1= e_lj
        call mpi_allreduce (unif1,unif2,1,mpi_real8,mpi_sum, &
                            mpi_comm_world,ierror)
        e_lj= unif2
!
        if(iwrt1.eq.0 .and. io_pe.eq.1) then
          OPEN (unit=11,file=praefixc//'.06'//suffix2,             &
                status='unknown',position='append',form='formatted')
          write(11,*) 'e_lj/npqr=',e_lj/npqr
          close(11)
        end if
      end if
!
      return
      end subroutine sht_forces 
!
!
!----------------------------------------------------------------------
      function f_cut (r,rcut_clf)
!----------------------------------------------------------------------
      use, intrinsic :: iso_c_binding 
!
      implicit none
      real(C_DOUBLE) f_cut,r,rcut_clf
!
      if(r.lt.rcut_clf) then
        f_cut= 1.d0
      else
        f_cut= exp(1.d0 -(r/rcut_clf)**2)
      end if         ! rcut_clf=5 -> r/5=...  
!
      return
      end function f_cut
!
!
!---------------------------------------------------------------
      function diel (xg,yg,zg)
!------------------------------------------------ 8/31/2006 ----
      use, intrinsic :: iso_c_binding 
      implicit none
!
      include 'paramAPGa.h'
!
      integer(C_INT) ifqq
      real(C_DOUBLE) xg,yg,zg,rr,diel,dielv,diel_m,     &
                     qfrac,Rpore,Hpore,Hpore2,Zci,Zcp,Zcn
      real(C_DOUBLE) diel2,dielpr,dielp,tau,ddv,aa,bb,r1,z1
!
      common/cntion/ qfrac,Rpore,Hpore,Zci,Zcp,Zcn,ifqq
      common/dielec/ diel2,dielpr,aa
!
      integer(C_INT) io_pe
      common/sub_proc/ io_pe
      logical        first
      data           first /.true./
!
! -------------------------
      aa =   2.0d0   ! 1.5 Ang ?
      bb =   3.0d0   ! vertical depth
!
      dielv = 1.0d0
      dielp = dielpr/79.d0
      diel_m= diel2 /79.d0   !<-- diel_m= diel2/79= 2/79 at /READ_CONF/
! -------------------------
!
      if(first) then
        if(io_pe.eq.1) then
          OPEN (unit=11,file=praefixc//'.06'//suffix2,             &
                status='unknown',position='append',form='formatted')
!
          write(11,*) ' dielp=',dielp
          write(11,*) ' aa,diel[v,p,m]=',aa,dielv,dielpr/79.,diel_m
          close(11)
        end if
!
        first= .false.
      end if
!
      Hpore2= 0.5d0*Hpore
      r1 =  Rpore  -aa
      z1 =  Hpore2 +bb
!
! 1) Water
      diel= dielv  
      rr= sqrt(xg**2 +yg**2)
!
      if(abs(zg).lt.z1) then
!
! 2) Center axis
        diel= dielv -(dielv -dielp)*(z1-abs(zg))/z1
!
! 3) Side of a membrane
        if(abs(zg).lt.Hpore2) then
          if(rr.gt.r1) then
            ddv = dielv -(dielv -dielp)*(z1-abs(zg))/z1
            diel= ddv -(ddv -diel_m)*min((rr-r1)/aa, 1.d0) !<-- pore
          end if
        end if
!
! 4) Topside of a membrane
        if(abs(zg).ge.Hpore2) then
          if(rr.gt.Rpore) then
            diel= dielv -(dielv -diel_m)*(z1-abs(zg))/bb
!
! 5) in-between
          else if(rr.gt.r1) then
            diel= ( dielv -(dielv -diel_m)*(z1-abs(zg))/bb     &
                  + dielv -(dielv -dielp)*(z1-abs(zg))/z1 )/2.d0
          end if
        end if
!
      end if
!
      return
      end function diel 
!
!
!------------------------------------------------------------
      subroutine sprmul (xg,yg,zg,ch,ag,fpx,fpy,fpz,np)
!------------------------------------------------------------
!  >>> Inextensible chains + bond angle rigidity (ggg).
!
      use, intrinsic :: iso_c_binding 
      implicit none
!
      include   'paramAPGa.h'
!    
      real(C_DOUBLE),dimension(npqr0) :: xg,yg,zg,ch,ag
      real(C_DOUBLE),dimension(np0) ::   fpx,fpy,fpz
      real(C_DOUBLE) xleng,yleng,zleng,  &
                     fgm,fgm0,rr,rrmax
      common/parm8/  xleng,yleng,zleng
!
      integer(C_INT) nsg,nseg,np
      common/psegm/  nsg(30),nseg
!
      real(C_DOUBLE) e_c_s,e_c_pme,e_c_r,e_lj,e_elas
      common/energy/ e_c_s,e_c_pme,e_c_r,e_lj,e_elas
!
      real(C_DOUBLE) t8,cptot
      common/headr3/ t8,cptot
!
      integer(C_INT) ia,ib,ns,i,j,jj,np00,nnb,ist1,ist2
      real(C_DOUBLE) bond_ps,bond_ss
      common/pbase/ np00,nnb,ist1(100),ist2(100)
      common/pbond/ bond_ps,bond_ss
!
      integer(C_INT) io_pe
      common/sub_proc/ io_pe
!
      integer(C_INT) iwrt1,iwrt2,iwrt3,iwrt4
      common/iotim/ iwrt1,iwrt2,iwrt3,iwrt4
!
!   fgm0= 1.5d0/2= 0.75d0 
      fgm0= 0.75d0
!
      do ns= 1,nseg     !<-- nseg=1 single, nseg=2 double stranded
      ia= nsg(ns) +1    
      ib= nsg(ns+1)
!
      do i= ia,ib 
      fpx(i)= 0
      fpy(i)= 0
      fpz(i)= 0
      end do
      end do
!                            nps = np/nseg
      do ns= 1,nseg     !<-- nseg=1 single, nseg=2 double stranded
      ia= nsg(ns) +1    
      ib= nsg(ns+1)
!
      do i= ia,ib       !<- ia--ib
      fgm= fgm0
!     rr= sqrt((xg(i)-xg(i+1))**2 +(yg(i)-yg(i+1))**2 +(zg(i)-zg(i+1))**2)
!     rrmax= bond_ps
!     fgm= fgm0*exp(-(rr/rrmax)**2) 
!
      if(i.eq.ia) then
      fpx(i)= fpx(i) -fgm*(xg(i)-xg(i+1))
      fpy(i)= fpy(i) -fgm*(yg(i)-yg(i+1))
      fpz(i)= fpz(i) -fgm*(zg(i)-zg(i+1))
!
      else if(i.eq.ib) then
      fpx(i)= fpx(i) -fgm*(xg(i)-xg(i-1))
      fpy(i)= fpy(i) -fgm*(yg(i)-yg(i-1))
      fpz(i)= fpz(i) -fgm*(zg(i)-zg(i-1))
!
      else 
      fpx(i)= fpx(i) -fgm*(2.d0*xg(i)-xg(i+1)-xg(i-1))
      fpy(i)= fpy(i) -fgm*(2.d0*yg(i)-yg(i+1)-yg(i-1))
      fpz(i)= fpz(i) -fgm*(2.d0*zg(i)-zg(i+1)-zg(i-1))
      end if
      end do
      end do 
!
!  -a-base (base attached to a sugar ring
!
      do jj= 1,nnb
      i= ist1(jj)
      j= ist2(jj)
!
      fgm= fgm0
!     rr= sqrt((xg(i)-xg(j))**2 +(yg(i)-yg(j))**2 +(zg(i)-zg(j))**2)
!     rrmax= ag(i) +ag(j) 
!     fgm= fgm0*exp(-(rr/rrmax)**2)
!
      fpx(i)= fpx(i) -fgm*(xg(i) -xg(j))
      fpy(i)= fpy(i) -fgm*(yg(i) -yg(j))
      fpz(i)= fpz(i) -fgm*(zg(i) -zg(j))
!
      fpx(j)= fpx(j) -fgm*(xg(j) -xg(i))
      fpy(j)= fpy(j) -fgm*(yg(j) -yg(i))
      fpz(j)= fpz(j) -fgm*(zg(j) -zg(i))
      end do
!
!
      e_elas= 0.d0
!
      do ns= 1,nseg   !<-- do ns
      ia= nsg(ns) +1 
      ib= nsg(ns+1)
!
      do i= ia,ib-1  
      fgm= fgm0
      e_elas= e_elas +0.5d0*fgm* &
                        ((xg(i) -xg(i+1))**2 +(yg(i) -yg(i+1))**2 &
                        +(zg(i) -zg(i+1))**2) 
      end do
      end do
!
      if(iwrt1.eq.0 .and. io_pe.eq.1) then
        OPEN (unit=11,file=praefixc//'.06'//suffix2,             &
              status='unknown',position='append',form='formatted')
        write(11,*) 'e_elas/np=',e_elas/np
        close(11)
      end if
!
      return
      end subroutine sprmul 
!
!
!----------------------------------------------------------------
      subroutine charges (rho,x,y,z,ch,ql,g,ipar,nCLp)
!----------------------------------------------------------------
!* Even meshes are always used.
      use, intrinsic :: iso_c_binding 
      implicit none
!
      include    'mpif.h'
      include    'paramAPGa.h'
!
      integer(C_INT) ipar 
!
      real(C_DOUBLE),dimension(0:mx-1,0:my-1,0:mz-1) :: rho,rho2 
      real(C_DOUBLE),dimension(0:npqr0-1) :: x,y,z,ch
      real(C_DOUBLE),dimension(0:ip0**3-1,0:npq0-1) :: ql
      real(C_DOUBLE),dimension(-1:1) :: ffx,ffy,ffz
      real(C_DOUBLE) t1,t2,t3,d1,d3,xx,yy,zz
!
      integer(C_INT),dimension(0:npq0-1,0:2) :: g
      integer(C_INT) i,j,k,m,m0,ik,xpos,ypos,zpos,nclp,xps,yps,zps
!
      real(C_DOUBLE) xmax,ymax,zmax,xmin,ymin,zmin
      real(C_float)  phi,tht,dtwr1,dtwr2,dtwr3
      common/parm3/  xmax,ymax,zmax,xmin,ymin,zmin
      common/parm4/  phi,tht,dtwr1,dtwr2,dtwr3
!
      real(C_DOUBLE) gx,gy,gz,ghx,ghy,ghz,ghx2,ghy2,ghz2, &
                     ghxsq,ghysq,ghzsq,hxi,hyi,hzi
      integer(C_INT) ptx,pty,ptz
      common/xtable/ gx(0:mx),ghx(0:mx),ghx2(0:mx),ghxsq(0:mx), &
                     gy(0:my),ghy(0:my),ghy2(0:my),ghysq(0:my), &
                     gz(0:mz),ghz(0:mz),ghz2(0:mz),ghzsq(0:mz), &
                     hxi,hyi,hzi
      common/xtabl2/ ptx(-10:3000),pty(-10:3000),ptz(-10:4000)
!
      integer(C_INT) pxr,pxc,pxl,pyr,pyc,pyl,pzr,pzc,pzl
      common/ptable/ pxr(-1:mx+1),pxc(-1:mx+1),pxl(-1:mx+1), &
                     pyr(-1:my+1),pyc(-1:my+1),pyl(-1:my+1), &
                     pzr(-1:mz+1),pzc(-1:mz+1),pzl(-1:mz+1)
!
      integer(C_INT) io_pe,ierror
      common/sub_proc/ io_pe
!
!* Non-periodic case
!
      do k= 0,mz-1
      do j= 0,my-1 
      do i= 0,mx-1
      rho(i,j,k) = 0.d0
      end do
      end do
      end do
!
!  ip0= 3, Non-periodic case in (-l/2,l/2).
!  definition of d1, g() are modified.
!
!* Need all definition for ql() and g() in /moldyn/
!
      do i= 0,nCLp-1  
!
      d1  = (x(i) -xmin)*hxi
      xps = ptx(int(d1 + 0.5d0))
      g(i,0) = pxc(xps)
!
      xx = (x(i) -gx(xps))/ghx(xps)
      ffx(-1)= 0.50d0*(0.5d0-xx)*(0.5d0-xx)
      ffx( 0)= 0.75d0-xx*xx
      ffx( 1)= 0.50d0*(0.5d0+xx)*(0.5d0+xx)
!*      
      d1  = (y(i) -ymin)*hyi
      yps = pty(int(d1 + 0.5d0))
      g(i,1) = pyc(yps)
!
      yy = (y(i) -gy(yps))/ghy(yps) 
      ffy(-1)= 0.50d0*(0.5d0-yy)*(0.5d0-yy)
      ffy( 0)= 0.75d0-yy*yy
      ffy( 1)= 0.50d0*(0.5d0+yy)*(0.5d0+yy)
!*      
      d3  = (z(i) -zmin)*hzi
      zps = ptz(int(d3 + 0.5d0))
      g(i,2) = pzc(zps)
!
      zz = (z(i) -gz(zps))/ghz(zps) 
      ffz(-1)= 0.50d0*(0.5d0-zz)*(0.5d0-zz)
      ffz( 0)= 0.75d0-zz*zz
      ffz( 1)= 0.50d0*(0.5d0+zz)*(0.5d0+zz)
!
      m0= -1 
!
!* Symmetric
!
      do j = -1,1            ! must be centered at 0
      do k = -1,1
      do m = -1,1
      t1 = ch(i) * ffx(j)
      t2 = t1 *    ffy(k)
      t3 = t2 *    ffz(m)
!    
      m0= m0 + 1
      ql(m0,i) = t3          ! assignment factor.
      end do
      end do
      end do
!
      end do
!      
!* The index of this loop must be partial (no double count).
!  round-robin partitioning is used for nclp particles.
!
!           ++++++ ++++++
      do i= ipar-1,nCLp-1,num_proc 
      m0= -1
!
      do j = -1,1
      do k = -1,1
      do m = -1,1
      xpos = pxc(g(i,0) + j)
      ypos = pyc(g(i,1) + k)
      zpos = pzc(g(i,2) + m)
!    
      m0= m0 +1
      rho(xpos,ypos,zpos) = rho(xpos,ypos,zpos) + ql(m0,i)
      end do
      end do
      end do
!
      end do
!
!
      if(num_proc.ne.1) then
!
        do k= 0,mz-1
        do j= 0,my-1 
        do i= 0,mx-1
        rho2(i,j,k)= rho(i,j,k)
        end do
        end do
        end do
!
        call mpi_allreduce (rho2,rho,mxyz,mpi_real8,mpi_sum, &
                            mpi_comm_world,ierror)
      end if
!
      return
      end subroutine charges 
!
!
!----------------------------------------------------------
      subroutine Labels
!----------------------------------------------------------
!*  Indices of the neighboring (27) sub-boxes.
!   Bounded case 
!
      use, intrinsic :: iso_c_binding 
      implicit none
!
      include  'paramAPGa.h'
      integer(C_INT) i,j,k,n,ip,im,jp,jm,kp,km,icmax
!
      integer(C_INT) isize,isizeZ,isize2,isize4,nc3
      integer(C_INT) ibind
      parameter  (isize=13,isizeZ=25,nc3=isize**2*isizeZ)
      parameter  (isize2=isize*isize,isize4=isize2+isize)
      common /boxind/ ibind(27,nc3)
!
      icmax= nc3 + 1
!
      do i = 1,isize
      ip = i + 1
      im = i - 1
      if (i.eq.isize) ip = isize ! 1     ! must be excluded 
      if (i.eq.1) im = 1     ! isize ! - otherwise, double counted
!
      do j = 1,isize
      jp = j + 1
      jm = j - 1
      if (j.eq.isize) jp= isize 
      if (j.eq.1) jm= 1 
!
      do k = 1,isizeZ
      kp = k + 1
      km = k - 1
      if (k.eq.isizeZ) kp= isizeZ
      if (k.eq.1) km = 1
!
      n = i + isize*(j-1) + isize2*(k-1)
!
      ibind( 1,n) = i  +  isize*( j-1) +   k*isize2 - isize4 
      ibind( 2,n) = ip +  isize*( j-1) +   k*isize2 - isize4
      ibind( 3,n) = ip +  isize*(jm-1) +   k*isize2 - isize4
      ibind( 4,n) = ip +  isize*(jp-1) +   k*isize2 - isize4
      ibind( 5,n) = i  +  isize*(jp-1) +   k*isize2 - isize4 
!
      ibind( 6,n) = i  +  isize*( j-1) +   kp*isize2 - isize4 
      ibind( 7,n) = im +  isize*( j-1) +   kp*isize2 - isize4
      ibind( 8,n) = ip +  isize*( j-1) +   kp*isize2 - isize4
      ibind( 9,n) = i  +  isize*(jm-1) +   kp*isize2 - isize4
      ibind(10,n) = im +  isize*(jm-1) +   kp*isize2 - isize4
      ibind(11,n) = ip +  isize*(jm-1) +   kp*isize2 - isize4 
      ibind(12,n) = i  +  isize*(jp-1) +   kp*isize2 - isize4
!
      ibind(13,n) = im +  isize*(jp-1) + kp*isize2 - isize4
      ibind(14,n) = ip +  isize*(jp-1) + kp*isize2 - isize4
      ibind(15,n) = i  +  isize*(jm-1) + km*isize2 - isize4
      ibind(16,n) = im +  isize*(jm-1) + km*isize2 - isize4
      ibind(17,n) = ip +  isize*(jm-1) + km*isize2 - isize4
      ibind(18,n) = i  +  isize*(jp-1) + km*isize2 - isize4
      ibind(19,n) = im +  isize*(jp-1) + km*isize2 - isize4
      ibind(20,n) = ip +  isize*(jp-1) + km*isize2 - isize4
!
      ibind(21,n) = ip +  isize*( j-1) + km*isize2 - isize4
      ibind(22,n) = im +  isize*( j-1) + km*isize2 - isize4
      ibind(23,n) = i  +  isize*( j-1) + km*isize2 - isize4
      ibind(24,n) = im +  isize*(jp-1) +  k*isize2 - isize4
      ibind(25,n) = im +  isize*(jm-1) +  k*isize2 - isize4
      ibind(26,n) = im +  isize*( j-1) +  k*isize2 - isize4
      ibind(27,n) = i  +  isize*(jm-1) +  k*isize2 - isize4
      end do
      end do
      end do
!
! Periodic case
!     ibind( 1,n) = i  +  isize*( j-1 + isize*(k-1))
!     ibind(26,n) = im +  isize*( j-1 + isize*(k-1))
!     ibind( 2,n) = ip +  isize*( j-1 + isize*(k-1))
!     ibind(27,n) = i  +  isize*(jm-1 + isize*(k-1))
!     ibind(25,n) = im +  isize*(jm-1 + isize*(k-1))
!     ibind( 3,n) = ip +  isize*(jm-1 + isize*(k-1))
!     ibind( 5,n) = i  +  isize*(jp-1 + isize*(k-1))
!     ibind(24,n) = im +  isize*(jp-1 + isize*(k-1))
!     ibind( 4,n) = ip +  isize*(jp-1 + isize*(k-1))
!
!     ibind( 6,n) = i  +  isize*( j-1 + isize*(kp-1))
!     ibind( 7,n) = im +  isize*( j-1 + isize*(kp-1))
!     ibind( 8,n) = ip +  isize*( j-1 + isize*(kp-1))
!     ibind( 9,n) = i  +  isize*(jm-1 + isize*(kp-1))
!     ibind(10,n) = im +  isize*(jm-1 + isize*(kp-1))
!     ibind(11,n) = ip +  isize*(jm-1 + isize*(kp-1))
!     ibind(12,n) = i  +  isize*(jp-1 + isize*(kp-1))
!     ibind(13,n) = im +  isize*(jp-1 + isize*(kp-1))
!     ibind(14,n) = ip +  isize*(jp-1 + isize*(kp-1))
!
!     ibind(23,n) = i  +  isize*( j-1 + isize*(km-1))
!     ibind(22,n) = im +  isize*( j-1 + isize*(km-1))
!     ibind(21,n) = ip +  isize*( j-1 + isize*(km-1))
!     ibind(15,n) = i  +  isize*(jm-1 + isize*(km-1))
!     ibind(16,n) = im +  isize*(jm-1 + isize*(km-1))
!     ibind(17,n) = ip +  isize*(jm-1 + isize*(km-1))
!     ibind(18,n) = i  +  isize*(jp-1 + isize*(km-1))
!     ibind(19,n) = im +  isize*(jp-1 + isize*(km-1))
!     ibind(20,n) = ip +  isize*(jp-1 + isize*(km-1))
!
!* Nullify if ibind() is out of bounds
!
      do n= 1,nc3
      do k= 1,27
      if(ibind(k,n).gt.nc3) ibind(k,n)= 0
      end do
      end do

      return
      end subroutine Labels
!
!
!---------------------------------------------------------------------
      subroutine reflect_endpl (xg,yg,zg,vx,vy,vz,ch,am,ag,np,nq,npqr)
!---------------------------------------------------------------------
      use, intrinsic :: iso_c_binding 
      implicit none
!
      include   'paramAPGa.h'
!
      real(C_DOUBLE),dimension(npqr0) :: xg,yg,zg,ch,vx,vy,vz,am,ag 
      real(C_DOUBLE) xmax3,ymax3,zmax3,xmin3,ymin3,zmin3
!
      integer(C_INT) isize,isizeZ,isize2,isize4,nc3
      parameter  (isize=13,isizeZ=25,nc3=isize**2*isizeZ)
      parameter  (isize2=isize*isize,isize4=isize2+isize)
!
      integer(C_INT) io_pe
      common/sub_proc/ io_pe
!
      real(C_DOUBLE) t8,cptot
      common/headr3/ t8,cptot
!
      integer(C_INT) np,nq,npqr,i,ifqq,ifrecyc,is
      real(C_DOUBLE) xmax,ymax,zmax,xmin,ymin,zmin,      &
                     xleng,yleng,zleng,                  & 
                     Hpore2,rmax,ddx,ddy,ddz,dd,         &
                     qfrac,Rpore,Hpore,Zci,Zcp,Zcn,      &
                     rr,rr2,zpor1,zpor2,vpara
!
      real(C_float)  phi,tht,dtwr1,dtwr2,dtwr3
      common/parm3/  xmax,ymax,zmax,xmin,ymin,zmin
      common/parm4/  phi,tht,dtwr1,dtwr2,dtwr3
      common/parm8/  xleng,yleng,zleng
      common/cntion/ qfrac,Rpore,Hpore,Zci,Zcp,Zcn,ifqq
!
      real(C_DOUBLE) pi,dt,axi,Gamma,rbmax,vth,tmax
      common/parm2/  pi,dt,axi,Gamma,rbmax,vth,tmax
!
      real(C_float)  t,t00,xp_leng
      common/headr2/ t,t00,xp_leng
!
!-----------------------------------------
!* Center of the system is (0,0,0)
!-----------------------------------------
!* Particle box is smaller than the field box.
!* grids i= 0, ..., (mx-1)
!        ix= int(isize* (xg(i)-xmin)/xleng +1.5001)
!        if(ix.le.1 .or. ix.ge.isize ) go to 210
!
      ddx = 0.5d0* xleng/isize
      ddy = 0.5d0* yleng/isize
      ddz = 0.5d0* zleng/isizeZ
!
      do i= 1,npqr
      Hpore2 = 0.5d0*Hpore           ! also modify /init/
!
! 1) Outside the pore region, hitting the box sides.
! 
      if(abs(zg(i)).gt.Hpore2) then  ! Outside region
!
        xmax3 =  xmax -ddx !-ag(i)
        ymax3 =  ymax -ddy !-ag(i)
        zmax3 =  zmax -ddz
!
        xmin3 =  xmin +ddx !+ag(i)
        ymin3 =  ymin +ddy !+ag(i)
        zmin3 =  zmin +ddz
!
!* X,Y sides are closed
!
        if(xg(i).lt.xmin3) then
          xg(i) = 2.d0*xmin3 -xg(i)
          vx(i)= -vx(i)
!
        else if(xg(i).gt.xmax3) then
          xg(i) = 2.d0*xmax3 -xg(i)
          vx(i)= -vx(i)
        end if
!
        if(yg(i).lt.ymin3) then
          yg(i) = 2.d0*ymin3 -yg(i)
          vy(i)= -vy(i)
!
        else if(yg(i).gt.ymax3) then
          yg(i) = 2.d0*ymax3 -yg(i)
          vy(i)= -vy(i)
        end if
!*
!* Z boundaries are closed for i=1,np+nq, but 
!  they are open in i=np+nq+1,np+nq+nr
!
!       if(i.le.np+nq) then
!
          if(zg(i).lt.zmin3) then
            zg(i) = 2.d0*zmin3 -zg(i)
            vz(i)= -vz(i)
!
          else if(zg(i).gt.zmax3) then
            zg(i) = 2.d0*zmax3 -zg(i)
            vz(i)= -vz(i)
          end if
!
!       else if(i.gt.np+nq) then
!
!         if(zg(i).lt.zmin3) then
!           zg(i)= zg(i) +zmax3
!
!         else if(zg(i).gt.zmax3) then
!           zg(i)= zg(i) +zmin3
!         end if
!       end if
!*
!
! 2) Within the pore region
!     Specular reflection: v'= v - 2*v_para
!
      else if(abs(zg(i)).le.Hpore2) then
!
        rr= sqrt(xg(i)**2 +yg(i)**2)
!
        Rmax = Rpore -ag(i)
        if(rr.gt.Rmax) then
!                     
          rr2= 2.d0*Rmax -rr
          xg(i)= rr2* xg(i)/rr
          yg(i)= rr2* yg(i)/rr
!
          vpara= (vx(i)*xg(i) +vy(i)*yg(i))/rr**2
          vx(i)= vx(i) -2.d0*vpara*xg(i)
          vy(i)= vy(i) -2.d0*vpara*yg(i)
!
!         else
!           if(zg(i).lt.0.d0) then    ! bottom side
!             zpor1 = -Hpore2 
!             zg(i)= 2.d0*zpor1 -zg(i)
!             vz(i)= -vz(i)
!           else
!             zpor2 = Hpore2
!             zg(i)= 2.d0*zpor2 -zg(i) ! top side
!             vz(i)= -vz(i)
!           end if
!         end if
       end if
!*
      end if
      end do
!
      return
      end subroutine reflect_endpl 
!
!
!  Read /write configuration data.
!------------------------------------------------------------------
      subroutine READ_CONF (np,nq,npqr,ifbase)
!------------------------------------------------------------------
      use, intrinsic :: iso_c_binding 
      implicit none
!
      include  'paramAPGa.h'
!
      integer(C_INT) kstart1,np,nq,npqr,ifbase
      integer(C_INT) io_pe
      common/sub_proc/ io_pe
!
      character praefix*6,text1*40
      integer(C_INT) nsg,nseg,mpc,ladabst,n_p,n_lp,v_g,n_smol, &
                     v_sp,v_sm,seed,Nzions,npartcls,           &
                     SubBox  
!
      real(C_DOUBLE) rcut_clf,rcutlj,r_fene,     &
                     skin,skin2,rcut_clf2,rcps2,         &
                     k_fene,Bjerrum,r_fene2
      common /poly/    n_p,mpc,ladabst
      common /cutoffel/ rcut_clf2,rcps2
      common /fenepara/ k_fene,r_fene2
      common /elsta/ Bjerrum
      common /einles/ skin,skin2
!
      common /confdatai/ n_lp,v_g,n_smol,v_sp,v_sm,seed
      common /confdatar/ rcut_clf,rcutlj,r_fene
      common /confdatas/ praefix
!----------------------------------------------------------------
      real(C_DOUBLE) pi,dt,axi,Gamma,rbmax,vth,tmax,tmin, &
                     xmax,ymax,zmax,xmin,ymin,zmin,       &
                     xleng,yleng,zleng,qfrac,Rpore,Hpore, &
                     Zci,Zcp,Zcn,acount,acoion
      real(C_float)  phi,tht,dtwr1,dtwr2,dtwr3
      integer(C_INT) ifqq
!
      common/parm2/  pi,dt,axi,Gamma,rbmax,vth,tmax
      common/psegm/  nsg(30),nseg
      common/parm3/  xmax,ymax,zmax,xmin,ymin,zmin
      common/parm4/  phi,tht,dtwr1,dtwr2,dtwr3
      common/parm8/  xleng,yleng,zleng
      common/cntion/ qfrac,Rpore,Hpore,Zci,Zcp,Zcn,ifqq
      common/ionsiz/ acount,acoion   
!
      integer(C_INT) ifLJ,ifrgrod,ifrodco
      common/parmlj2/ ifLJ
      common/parmif/ ifrgrod,ifrodco
!
      real(C_float)  t,t00,xp_leng
      real(C_DOUBLE) t8,cptot
      common/headr2/ t,t00,xp_leng
      common/headr3/ t8,cptot
!
      real(C_DOUBLE)  t_pe,t_poisn
      common/timings/ t_pe,t_poisn
!
      real(C_DOUBLE) Vtop0,Vbot0,Vtop,Vbot,t_recyc
      integer(C_INT) fixedbc_x,fixedbc_y,itermax,filtx,filty,filtz
      common/espot/  Vtop0,Vbot0,Vtop,Vbot,fixedbc_x,fixedbc_y, &
                     itermax,filtx,filty,filtz
      common/trecyc/ t_recyc
!
      real(C_DOUBLE) xgc,ygc,zgc,vxg,vyg,vzg,rod_leng,Rmac, &
                     Rhelix
      integer(C_INT) n_rodp
      common/macroion/ xgc,ygc,zgc,vxg,vyg,vzg,rod_leng,Rmac,  &
                       Rhelix,n_rodp
      real(C_DOUBLE),dimension(np0) :: &
                                dox,doy,doz,doxc,doyc,dozc
      common/macroin1/ dox,doy,doz,doxc,doyc,dozc
!
      real(C_DOUBLE) diel2,dielpr,aa
      common/dielec/ diel2,dielpr,aa
!
!----------------------------------------------------------------
!  All ranks must do:
      OPEN (unit=08,file=praefixs//'_config.start'//suffix0, &
                                            form='formatted')
!
      if(io_pe.eq.1) then
        OPEN (unit=11,file=praefixc//'.06'//suffix2,             &
              status='unknown',position='append',form='formatted')
!
        write(11,*) 'READ_CONF: parameter read... start'
        close(11)
      end if
!*
      read(08,'(a40,a6)') text1, praefix    ! String der Simulationserkennung
      kstart1= 0
      ifqq   = 1 !!
!
!     read(08,'(a40,i12)') text1,kstart     !*Restart sequence number 0
      read(08,'(a40,i12)') text1,ifqq       !*Non-neutral Pl (if=1)   1 
!
      read(08,'(a40,f20.0)') text1,cptot    ! Maximum cpu time Min.  1430.
      read(08,'(a40,f20.0)') text1,tmax     ! Physikalische Endzeit
      read(08,'(a40,f20.0)') text1,dt       ! Diskretisierungs-Schritt
      read(08,'(a40,f20.0)') text1,dtwr1    ! Write out for iwrt1, by DDt=1 
      read(08,'(a40,f20.0)') text1,dtwr2    ! Write out for iwrt2 
      read(08,'(a40,f20.0)') text1,dtwr3    ! Write out for iwrt3 
!
      read(08,'(a40,f20.0)') text1,Rpore    ! Radius of a nanopore  5.4 Ang
      read(08,'(a40,f20.0)') text1,Hpore    ! Height of a nanopore 56.0 Ang
      read(08,'(a40,f20.0)') text1,diel2    ! Dielectric constant of membrane diel2=2
!
      read(08,'(a40,i12)') text1,n_p        ! Zahl der Polymere 0,1 or 2
      read(08,'(a40,i12)') text1,n_lp       ! Zahl der Ladungen pro Polymer 24
      read(08,'(a40,i12)') text1,ifbase     ! =1 if neutral base is attached
!
      read(08,'(a40,f20.0)') text1,Rmac     ! radius of a rod  Rmac= 3.9d0
      read(08,'(a40,f20.0)') text1,rod_leng ! length of a rod  rod_leng= 42.0d0
      read(08,'(a40,f20.0)') text1,Rhelix   ! Radius for rod charges Rhelix=4.9Ang 
      read(08,'(a40,f20.0)') text1,dielpr   ! Dielectric constant of water dielpr=79 
!
      Ladabst= 1
      qfrac= 0.d0 
!
      read(08,'(a40,i12)') text1,v_g        ! Valenz der Gegen(salz)ionen     1.
      read(08,'(a40,i12)') text1,v_sm       ! Valenz der Negativen Salzionen -1.
      read(08,'(a40,i12)') text1,NZions     ! Teilchenzahl (including Z:1 salt) 200
!  ***********************
!  (npqr) is redefined in /init/
!
      if(kstart.eq.0) then
        np  = n_p * n_lp
        nseg= n_p
!
        npartcls = np + Nzions
        npqr= npartcls
        nq  = Nzions
      end if
!  ***********************
      Zcp = v_g 
      Zcn = v_sm
!
      read(08,'(a40,f20.0)') text1,xleng    ! Box size: X  81 Ang
      read(08,'(a40,f20.0)') text1,yleng    ! Box size: Y
      read(08,'(a40,f20.0)') text1,zleng    ! Box size: Z 168 Ang
!
!  -----------------------
!   Symmetric of xmax and xmin 
      xmax=  0.5d0*xleng
      ymax=  0.5d0*yleng
      zmax=  0.5d0*zleng
!
      xmin= -0.5d0*xleng
      ymin= -0.5d0*yleng
      zmin= -0.5d0*zleng
!  -----------------------
!
      read(08,'(a40,f20.0)') text1,Vtop0    ! Potential of the top plate 0.
      read(08,'(a40,f20.0)') text1,Vbot0    ! Potential of the bot plate -0.
!
      if(Vtop0.lt.Vbot) then
!
        if(io_pe.eq.1) then
        OPEN (unit=11,file=praefixc//'.06'//suffix2,             &
              status='unknown',position='append',form='formatted')
!
        write(11,*) '# error - Vtop0 < Vbot0, stop...'
        close(11)
        end if
!
        stop
      end if
!
      read(08,'(a40,i12)') text1,fixedbc_x  ! Boundary: 0 (fixed), 1 (smooth)
      fixedbc_y =  fixedbc_x
!
      read(08,'(a40,i12)')   text1,ifLJ     ! ifLJ=1 if attractive LJ 
      read(08,'(a40,i12)')   text1,ifrgrod  ! ifrgrod=1 if rod is true
      read(08,'(a40,i12)')   text1,ifrodco  ! ifrodco=1 if rodco
      read(08,'(a40,f20.0)') text1,t00      ! Vtop starts at this time
      read(08,'(a40,f20.0)') text1,t_recyc  ! t_recyc starts at this time
!
      read(08,'(a40,f20.0)') text1,t_pe     ! Time to start moving PE 3.
      read(08,'(a40,f20.0)') text1,t_poisn  ! Time to start solving Poisson eq  1.
!
      read(08,'(a40,i12)')   text1,itermax  ! Maximum of Poisson iterations 200
      read(08,'(a40,i12)')   text1,filtx    ! Number of spatial fileterings 1
      filty= filtx
      filtz= filtx
!
      read(08,'(a40,i12)')   text1,SubBox   ! Subboxen pro Boxlaenge 1
      read(08,'(a40,f20.0)') text1,skin     ! Skin 1.166
!
!!      rcut_clf= no dimension               "position space"
      read(08,'(a40,f20.0)') text1,rcut_clf ! Ortsraum-Cutoff. = 7 Ang 
      read(08,'(a40,f20.0)') text1,acount   ! Size of counterion 1.33 Ang 0.95
      read(08,'(a40,f20.0)') text1,acoion   ! Size of coion      1.82 Ang  1.3
!
!!      rcutlj  = no dimension
      read(08,'(a40,f20.0)') text1,rcutlj   ! LJ-Cutoff, no dimension=1.122
      read(08,'(a40,f20.0)') text1,Bjerrum  ! Bjerrum-Laenge 7.0 Ang
!
      close (unit=08)
!*
      if(io_pe.eq.1) then
        OPEN (unit=11,file=praefixc//'.06'//suffix2,             &
              status='unknown',position='append',form='formatted')
!
        write(11,*) 'xleng,zleng (Ang)=',xleng,zleng
        write(11,*) 'Rpore,Hpore (Ang)=',Rpore,Hpore
        write(11,*) 'diel2=',diel2
        write(11,*) 'Rmac,rod_leng (Ang)=',Rmac,rod_leng
        write(11,*) 'acount,acoion (Ang)=',acount,acoion
        write(11,*) 'rcut_clf=',rcut_clf
        write(11,*) 'Bjerrum (Ang)=',Bjerrum
!
        write(11,*) 'READ_CONF: parameter read... end'
        close(11)
      end if
!
      return
      end subroutine READ_CONF 
!
!
!----------------------------------------------------------------
      subroutine WRITE_CONF (ende,ifbase,np,nq,npqr)
!----------------------------------------------------------------
      use, intrinsic :: iso_c_binding 
      implicit none
!
      include    'paramAPGa.h'
      integer(C_INT) ifbase,np,nq,npqr
      integer(C_INT) ifqq,mpc,ladabst,n_p,n_lp,v_g,n_smol, &
                     v_sp,v_sm,seed
      real(C_DOUBLE) rcut_clf,rcutlj,r_fene,       &
                     rcut_clf2,rcps2,k_fene,skin,  &
                     Bjerrum,r_fene2,skin2
      character   praefix*6
      logical     ende
!
      common /poly/     n_p,mpc,ladabst
      common /cutoffel/ rcut_clf2,rcps2
      common /fenepara/ k_fene,r_fene2
      common /elsta/ Bjerrum
      common /einles/ skin,skin2
!
      common /confdatai/ n_lp,v_g,n_smol,v_sp,v_sm,seed
      common /confdatar/ rcut_clf,rcutlj,r_fene
      common /confdatas/ praefix
!
      real(C_DOUBLE) qfrac,Rpore,Hpore,Zci,Zcp,Zcn
      common/cntion/ qfrac,Rpore,Hpore,Zci,Zcp,Zcn,ifqq
!
      real(C_DOUBLE) diel2,dielpr,aa
      common/dielec/ diel2,dielpr,aa
!
      real(C_DOUBLE) Vtop0,Vbot0,Vtop,Vbot,t_recyc,rpr
      integer(C_INT) fixedbc_x,fixedbc_y,itermax,filtx,filty,filtz,ipr
      common/espot/  Vtop0,Vbot0,Vtop,Vbot,fixedbc_x,fixedbc_y, &
                     itermax,filtx,filty,filtz
!----------------------------------------------------------------
      real(C_DOUBLE) pi,dt,axi,Gamma,rbmax,vth,tmax, &
                     xmax,ymax,zmax,xmin,ymin,zmin,  &
                     xleng,yleng,zleng,acount,acoion   
      common/parm8/  xleng,yleng,zleng
      real(C_float)  phi,tht,dtwr1,dtwr2,dtwr3
!
      common/parm2/ pi,dt,axi,Gamma,rbmax,vth,tmax
      common/parm3/ xmax,ymax,zmax,xmin,ymin,zmin
      common/parm4/ phi,tht,dtwr1,dtwr2,dtwr3
      common/ionsiz/ acount,acoion   
!
      real(C_DOUBLE) t8,cptot
      common/headr3/ t8,cptot
!
      real(C_DOUBLE) xgc,ygc,zgc,vxg,vyg,vzg,rod_leng,Rmac, &
                     Rhelix
      integer(C_INT) n_rodp
      common/macroion/ xgc,ygc,zgc,vxg,vyg,vzg,rod_leng,Rmac,  &
                       Rhelix,n_rodp
      real(C_DOUBLE),dimension(np0) :: &
                                dox,doy,doz,doxc,doyc,dozc
      common/macroin1/ dox,doy,doz,doxc,doyc,dozc
!----------------------------------------------------------------
!*********************************************************************
      write(07,'("# praefix-string................: ",a6)')praefix
      write(07,'("# Restart sequence number.......: ",i12)')kstart
      write(07,'("# Non-neutral polymer = 1.......: ",i12)')ifqq
      write(07,'("# Maximum cpu time of run.......: ",f20.12)')cptot
      write(07,'("# Physikalische endzeit.........: ",f20.12)')tmax
      write(07,'("# Diskretisierungs-schritt......: ",f20.12)')dt
      write(07,'("# write out interval for iwrt1..: ",f20.12)')dtwr1
      write(07,'("# write out interval for iwrt2..: ",f20.12)')dtwr2
      write(07,'("# write out interval for iwrt3..: ",f20.12)')dtwr3
!
      write(07,'("# Radius of the pore............: ",f20.12)')Rpore
      write(07,'("# Width of the pore.............: ",f20.12)')Hpore
      write(07,'("# Dielectric of membrane diel2..: ",f20.12)')diel2
!
      write(07,'("# Number of polymer segments....: ",i12)')n_p     
      write(07,'("# Number of monomers on backbone: ",i12)')n_lp
      write(07,'("# =1 if neutral base is attached: ",i12)')ifbase
!
      write(07,'("# Radius of a rod Rmac..........: ",f20.12)')Rmac
      write(07,'("# Length of a rod  rod_leng.....: ",f20.12)')rod_leng
      write(07,'("# Radius for rod charges Rhelix.: ",f20.12)')Rhelix
      write(07,'("# Dielectric constant of water..: ",f20.12)')dielpr
!
      v_g  = Zcp
      v_sm = Zcn
      write(07,'("# Valence of counterions........: ",i12)')v_g
      write(07,'("# Valence of negative salt ion..: ",i12)')v_sm
!
      write(07,'("# Boxlaenge (x).................: ",f20.12)')xleng
      write(07,'("# Boxlaenge (y).................: ",f20.12)')yleng
      write(07,'("# Boxlaenge (z).................: ",f20.12)')zleng
      write(07,'("# potential of the top plate....: ",f20.12)')Vtop0
      write(07,'("# boundary at the sides.........: ",i12)')fixedbc_x
      write(07,'("# Number of poisson iterations..: ",i12)')itermax
      write(07,'("# Skin..........................: ",f20.12)')skin
      write(07,'("# Ortsraum-cutoff...............: ",f20.12)')rcut_clf
      write(07,'("# Size of counterions...........: ",f20.12)')acount
      write(07,'("# Size of coions................: ",f20.12)')acoion
      write(07,'("# LJ-cutoff, no dimension.......: ",f20.12)')rcutlj
      write(07,'("# Bjerrum-laenge Ang............: ",f20.12)')Bjerrum
!*********************************************************************
!
      return
      end subroutine WRITE_CONF 
!
!
!-------------------------------------------------------------
      subroutine Initialisierung 
!-------------------------------------------------------------
      use, intrinsic :: iso_c_binding 
      implicit none
!
      include   'paramAPGa.h'
!
      integer(C_INT) n_g,n_sp,n_sm,n_si,mpc,n_p,ladabst, &
                     n_lp,v_g,n_smol,v_sp,v_sm,seed,     &
                     verletfix,qlj,i,iverg,i_steps,      &
                     measstep,confstep,nCLp
      real(C_DOUBLE) qwert,rcutlj,rcut_clf,skin,rcps2,     &
                     rcps,fmesh,skin2,rcut_clf2,r_fene2,rcutlj2, &
                     k_fene,r_fene,Bjerrum,wupi,           &
                     e_c_s,e_c_pme,e_c_r,e_lj,e_elas
!-----------
!                    **** 
      common/ewald2/ nCLp   ! ,ip0 paramAPGa.h
!
      common /inidati/ n_g,n_sp,n_sm,n_si
      common /inidatr/ qwert,fmesh
      common /chargelj/ qlj(npq0)
      common /poly/ n_p,mpc,ladabst
      common /einles/ skin,skin2
      common /cutoffel/ rcut_clf2,rcps2
      common /fenepara/ k_fene,r_fene2
      common /confdatai/ n_lp,v_g,n_smol,v_sp,v_sm,seed

      common /confdatar/ rcut_clf,rcutlj,r_fene
      common /energy/ e_c_s,e_c_pme,e_c_r,e_lj,e_elas
      common /elsta/ Bjerrum
!----------------------------------------------------------------
      real(C_DOUBLE) pi,dt,axi,Gamma,rbmax,vth,tmax, &
                    xmax,ymax,zmax,xmin,ymin,zmin,   &
                    xleng,yleng,zleng,acount,acoion   
      real(C_float) phi,tht,dtwr1,dtwr2,dtwr3
      common/parm2/ pi,dt,axi,Gamma,rbmax,vth,tmax
      common/parm3/ xmax,ymax,zmax,xmin,ymin,zmin
      common/parm4/ phi,tht,dtwr1,dtwr2,dtwr3
      common/parm8/ xleng,yleng,zleng
!----------------------------------------------------------------

      wupi   = 1.77245385090551602729816748334d0

      n_g   = n_p*n_lp/v_g              ! zahl der polymer-gegenionen 
      n_sp  = n_smol*v_sm               ! zahl der positiven salzionen 
      n_sm  = n_smol*v_sp               ! zahl der negativen salzionen 
      n_si  = n_sp+n_sm                 ! gesamtzahl der salzionen 
      mpc   = (n_lp - 1)*ladabst+1      ! zahl der monomere pro 
!                                         polymerkette 
      qwert = 1.0d0*(n_p*n_lp) &        ! summe der quadrate aller ladungen 
               +dfloat(v_g)**2  * ((n_p*n_lp)/v_g) &
               +dfloat(v_sp)**2 * n_sp             &
               +dfloat(v_sm)**2 * n_sm

      rcps      = rcut_clf + skin
  
      skin2     = skin   * skin
      rcut_clf2  = rcut_clf* rcut_clf
      r_fene2   = r_fene * r_fene
      rcutlj2   = rcutlj * rcutlj
      rcps2     = rcps   * rcps
  
      do i= 1,n_p*mpc
      qlj(i) = 1 
      end do
!
      do i= n_p*mpc +1,nCLp
      qlj(i) = 0
      end do

      if (rcut_clf.lt.rcutlj) then
        print*,'pme-cut-off zu klein:',rcut_clf,'->rcutlj',rcutlj
        rcut_clf = rcutlj
      end if
!
      if (rcps.gt.(0.5*xleng)) then
        print*,'rcps > 0.5*Len: minimum convention not satisfied'
      end if
!
      return
      end subroutine initialisierung
!
!
!/*---------------------------------------------------------------------
      subroutine Interpol_charge_assign_function
!/*---------------------------------------------------------------------
      use, intrinsic :: iso_c_binding 
      implicit    none
!
      include    'paramAPGa.h'
!
!     integer(C_INT) ip0,mintpol
!     parameter  (ip0=3,mintpol=4*50048)
!
      integer(C_INT) io_pe
      common/sub_proc/ io_pe
!
      integer(C_INT) i
      real(C_DOUBLE) intcaf, dinterpol, x
      common/intcaf/ intcaf(-1:1,-mintpol:mintpol)
!                    **************************** defined here.
!
      dinterpol= dfloat(mintpol)
!
      if(io_pe.eq.1) then
        OPEN (unit=11,file=praefixc//'.06'//suffix2,             &
              status='unknown',position='append',form='formatted')
!
        write(11,601) ip0 
  601   format(/,'- interpolating the order-',i2, &
               ' charge assignment function')
        close(11)
      end if
!
      if (ip0.eq.3) then
!        ++++++++++
!
        do i= -mintpol, mintpol
        x= i/(2.d0*dinterpol)
        intcaf(-1,i) = 0.50d0*(0.5d0 - x)**2
        intcaf( 0,i) = 0.75d0 - x*x
        intcaf( 1,i) = 0.50d0*(0.5d0 + x)**2
        end do
!
      else
        if(io_pe.eq.1) then
          OPEN (unit=11,file=praefixc//'.06'//suffix2,             &
                status='unknown',position='append',form='formatted')
!
          write(11,*) "Error in function ", &
                      "'interpolate_charge_assignment_function':"
          write(11,*) ip0
  611     format("charge assignment order",i2," unknown.",/, &
                 "program terminated.")
          close(11)
        end if
!
        stop
      end if
!
      return
      end subroutine interpol_charge_assign_function
!
!
!-----------------------------------------------------------
      subroutine foldbk (xg,yg,zg,x4,y4,z4,np,nq,npqr)
!-----------------------------------------------------------
!*******************************
!*  Periodic version.          *
!*******************************
      use, intrinsic :: iso_c_binding 
!
      include  'paramAPGa.h'
!
      real(C_DOUBLE),dimension(npqr0) :: xg,yg,zg
      real(C_float),dimension(npio) ::  x4,y4,z4
!
      do l= 1,npqr
      x4(l) = xg(l)
      y4(l) = yg(l)
      z4(l) = zg(l)
      end do
!
      return
      end subroutine foldbk
!
!
!--------------------------------------------------------------------
      subroutine init (xg,yg,zg,vx,vy,vz,ch,am,ag,ep,ifrgrod,ifrodco, &
                       ifbase,ipar,np,nq,nCLp,nr,npqr)
!--------------------------------------------------------------------
      use, intrinsic :: iso_c_binding 
      implicit none
!
      include   'mpif.h'
      include   'paramAPGa.h'
!
      integer(C_INT) ifrgrod,ifrodco,ifbase,ipar, &
                     np,nq,nCLp,nr,npqr
!
      integer(C_INT) io_pe
      common/sub_proc/ io_pe
!
      real(C_DOUBLE),dimension(npqr0) :: xg,yg,zg,vx,vy,vz, &
                                         ch,am,ag,ep
      real(C_DOUBLE) ch8,abin,pi,dt,axi,Gamma,rbmax, &
                     vth,tmax,                       &
                     xmax,ymax,zmax,xmin,ymin,zmin,  &
                     xleng,yleng,zleng,xcent,ycent,zcent, &
                     ww1,ww2,ddz,awat_diam,          &
                     KJ,KCal,mol,kbT,ranff   
!
      real(C_DOUBLE) dgaus2,vmax1,vmax2
      real(C_float)  phi,tht,dtwr1,dtwr2,dtwr3
      integer(C_INT) it,is,i,j,k,l,nsg,nseg,ifqq,nps,k1,       &
                     Nzi,Nci,Ncc_pore,Ncc_wide,n0,n1,n2,n3,n4, &
                     nn3,nn4,ntried,i0,i1,ll,kl,km,kn,istp
!
      common/parm1/ it,is
      common/parm2/ pi,dt,axi,Gamma,rbmax,vth,tmax
      common/parm3/ xmax,ymax,zmax,xmin,ymin,zmin
      common/parm4/ phi,tht,dtwr1,dtwr2,dtwr3
      common/parm8/ xleng,yleng,zleng
      common/psegm/ nsg(30),nseg
!
      real(C_DOUBLE) qfrac,Rpore,Hpore,Zci,Zcp,Zcn,Hpore2, &
                     wwat,awat,acount,acoion
      common/cntion/ qfrac,Rpore,Hpore,Zci,Zcp,Zcn,ifqq
      common/waterm/ wwat,awat
      common/ionsiz/ acount,acoion   
!
      real(C_DOUBLE) pi2,q_PE,aLJ,cci,dnint,rr,x0,y0,z0, &
                     rcut_clf,rcutlj,r_fene
      common/shrgf/ aLJ
      common/confdatar/ rcut_clf,rcutlj,r_fene
!
      real(C_DOUBLE)fv,vv,dv
      common/gaus1/ fv(101),vv,dv
!
      integer(C_INT) pxr,pxc,pxl,pyr,pyc,pyl,pzr,pzc,pzl
      common/ptable/ pxr(-1:mx+1),pxc(-1:mx+1),pxl(-1:mx+1), &
                     pyr(-1:my+1),pyc(-1:my+1),pyl(-1:my+1), &
                     pzr(-1:mz+1),pzc(-1:mz+1),pzl(-1:mz+1)
!
      real(C_DOUBLE) gx,gy,gz,ghx,ghy,ghz,ghx2,ghy2,ghz2,  &
                     ghxsq,ghysq,ghzsq,hxi,hyi,hzi,        &
                     ghx8(0:mx),ghy8(0:my),ghz8(0:mz),     &
                     ghx0,ghy0,ghz0,hx8,hy8,hz8,hz0,hhz0,  &
                     zz0,zm,hx0,hhx0,xx0,xm,hy0,hhy0,wg
      integer(C_INT) ptx,pty,ptz,icentx,icentz,ismall,k0
      common/xtable/ gx(0:mx),ghx(0:mx),ghx2(0:mx),ghxsq(0:mx), &
                     gy(0:my),ghy(0:my),ghy2(0:my),ghysq(0:my), &
                     gz(0:mz),ghz(0:mz),ghz2(0:mz),ghzsq(0:mz), &
                     hxi,hyi,hzi
      common/xtabl2/ ptx(-10:3000),pty(-10:3000),ptz(-10:4000)
!
      real(C_DOUBLE) xx,yy,zz,diel,xgc,ygc,zgc,vxg,vyg,vzg,    &
                     rod_leng,rod_len(100),Rmac,Rhelix,        &
                     th,ph,ph3,ss,airc,rxg(np0),ryg(np0),rzg(np0), &
                     q1,q2,q3,q4,angx,angy,angz,ipx,ipy,ipz,   &
                     ipz1,xyz1,xyz3,rr0
!
      integer(C_INT) n_rodp,jj,np00,nnb,ist1,ist2
      common/pbase/  np00,nnb,ist1(100),ist2(100)
!
      real(C_DOUBLE) bond_ps,bond_ss,a_phos,a_sugar,a_base1
      common/pbond/  bond_ps,bond_ss
!
      common/macroion/ xgc,ygc,zgc,vxg,vyg,vzg,rod_leng,Rmac,  &
                       Rhelix,n_rodp
      real(C_DOUBLE),dimension(np0) :: &
                                dox,doy,doz,doxc,doyc,dozc
      common/macroin1/ dox,doy,doz,doxc,doyc,dozc
      common/macroiov/ q1,q2,q3,q4,angx,angy,angz,ipx,ipy,ipz
!
      real(C_DOUBLE) Vtop0,Vbot0,Vtop,Vbot,t_recyc
      integer(C_INT) fixedbc_x,fixedbc_y,itermax,filtx,filty,filtz
      common/espot/  Vtop0,Vbot0,Vtop,Vbot,fixedbc_x,fixedbc_y, &
                     itermax,filtx,filty,filtz
!
      integer(C_INT) ifLJ,Ncc_1,Ncc_2,Ncc_3
      real(C_DOUBLE) epsLJ,eps_K,eps_Cl,rhow
      common/parmlj/ epsLJ,eps_K,eps_Cl
      common/parmlj2/ ifLJ
      common/water1/ rhow(0:mx-1,0:my-1,0:mz-1)
!
!-------------------------
!*  Grid quantity.
!-------------------------
!
      pi2= 2.d0*pi
!
! ----------------------------------
!* ptable: used in Poisson solver.
! ----------------------------------
!
      do i= 0,mx-1
      pxr(i)= i +1
      pxc(i)= i 
      pxl(i)= i -1
      end do
!
      pxr(mx-1)= mx-1
      pxl(0)   =  0
!
      pxr(mx)= mx-1
      pxc(mx)= mx-1
      pxl(mx)= mx-1
      pxr(-1) = 0
      pxc(-1) = 0
      pxl(-1) = 0
!
      do j= 0,my-1
      pyr(j)= j +1
      pyc(j)= j 
      pyl(j)= j -1
      end do
!
      pyr(my-1)= my-1
      pyl(0)   =  0
!
      pyr(my)= my-1
      pyc(my)= my-1
      pyl(my)= my-1
      pyr(-1) = 0
      pyc(-1) = 0
      pyl(-1) = 0
!
!
      do k= 0,mz-1
      pzr(k)= k +1
      pzc(k)= k 
      pzl(k)= k -1
      end do
!
      pzr(mz-1)= mz-1
      pzl(0)   =  0
!
      pzr(mz)= mz-1
      pzc(mz)= mz-1
      pzl(mz)= mz-1
      pzr(-1) = 0
      pzc(-1) = 0
      pzl(-1) = 0
!
!***********************************************************************
!  Define uneven mesh quantities.                                      *
!***********************************************************************
!-----------------------------------------------------------------------
!* Generate the grid spacing.
!
!          0                    zmax
!  grid#: .0..1..2.. ..(mz-2)..(mz-1)
!
!     cell 0 and (mz-1) are only halfways in the domain
!-----------------------------------------------------------------------
!  -----------------------
!     xmax=  0.5d0*xleng  !<- READ_CONF
!     ymax=  0.5d0*yleng
!     zmax=  0.5d0*zleng
!
!     xmin= -0.5d0*xleng
!     ymin= -0.5d0*yleng
!     zmin= -0.5d0*zleng
!  -----------------------
!
!  *** Z-direction ***     
!
      do k= 0,mz-1
      ghz8(k)= 1.0d0
      end do
!
      if(io_pe.eq.1) then
        OPEN (unit=11,file=praefixc//'.06'//suffix2,             &
              status='unknown',position='append',form='formatted')
!
        write(11,*)
        write(11,*) '--- Subroutine init (start) ----------------------'
        write(11,*)
        write(11,*) 'Even meshes are used.....'
        write(11,*) ' mx= my=',mx,' mz=',mz
        close(11)
      end if
!
      wg= 0.5d0*ghz8(0)
      do k= 1,mz-2
      wg= wg +ghz8(k)
      end do
!
      wg= wg +0.5d0*ghz8(mz-1)
!
!---------------------------------
      hz8= dble(zleng)/wg  !<-- hz8
      hz0= (1.d0/9.d0)*hz8
!---------------------------------
!*  Define the guard cell (mz).
!
      do k= 0,mz-1
      ghz8(k)= hz8*ghz8(k)
      ghz(k) = ghz8(k)
      end do
!
      gz(0)= zmin
      wg   = zmin +0.5d0*ghz8(0)
!
      do k= 1,mz-2
      wg= wg +0.5d0*ghz8(k)
      gz(k)= wg
      wg= wg +0.5d0*ghz8(k)
      end do
!
      gz(mz-1)= wg +0.5d0*ghz8(mz-1)
!
!
!  *** X,Y-directions ***     
!
      if(my.ne.mx) then
        if(io_pe.eq.1) then
          OPEN (unit=11,file=praefixc//'.06'//suffix2,             &
                status='unknown',position='append',form='formatted')
!
          write(11,*) 'Error in /init/: my.ne.mx, stop.....'
          close(11)
        end if
!
        stop
      end if
!
      do i= 0,mx-1
      ghx8(i)= 1.0d0
      end do
!
      wg= 0.5d0*ghx8(0)
      do i= 1,mx-2
      wg= wg +ghx8(i)
      end do
!
      wg= wg +0.5d0*ghx8(mx-1)
!
!---------------------------------
      hx8= dble(xleng)/wg  !<-- hx8
      hx0= (1.d0/9.d0)*hx8
!
      hy8= hx8
      hy0= hx0
!---------------------------------
!*  Define the guard cell (mx).
!
      do i= 0,mx-1
      ghx8(i)= hx8*ghx8(i)
      ghx(i) = ghx8(i)
!
      ghy8(i)= ghx8(i)   !<-- ghy = ghx
      ghy(i) = ghx(i) 
      end do
!
      gx(0)= xmin
      gy(0)= xmin
!
      wg   = xmin +0.5d0*ghx8(0)
!
      do i= 1,mx-2
      wg= wg +0.5d0*ghx8(i)
      gx(i)= wg
      gy(i)= gx(i)
      wg= wg +0.5d0*ghx8(i)
      end do
!
      gx(mx-1)= wg +0.5d0*ghx8(mx-1)
      gy(mx-1)= gx(mx-1)
!
!---------------------------------
!           * ****
      do i= 1,mx-2
      ghx2(i) =  2.d0*ghx(i)
      ghxsq(i)=  ghx(i)**2
      end do
!
      do j= 1,my-2
      ghy2(j) = ghx2(j)
      ghysq(j)= ghxsq(j)
      end do
!
      do k= 1,mz-2
      ghz2(k) =  2.d0*ghz(k)
      ghzsq(k)=  ghz(k)**2
      end do
!
      if(io_pe.eq.1) then
        OPEN (unit=11,file=praefixc//'.06'//suffix2,             &
              status='unknown',position='append',form='formatted')
!
        write(11,101) (i,gx(i),gy(i),gz(i),ghx(i),ghy(i),ghz(i),i= 0,mx-1)
        write(11,102) (i,gz(i),ghz(i),i=mx,mz-1)
  101   format('i=',i3,' gx,gy,gz,hgx,hgy,hgz=',6f10.5)
  102   format('i=',i3,' gx,gy,gz,hgx.hgy,hgz=',20x,f10.5,20x,f10.5)
        close(11)
      end if
!
!--------------------------------
!* Find the nearest grid.
!--------------------------------
! X
      i0 = 0
      ptx(0)= 0
      hhx0 = 0.3d0 * hx0
!
      do i= 1,3000
      xx0= xmin + hhx0*i
      if(xx0.ge.xmax) go to 120
!
      i1= i
      if(i0.le.mx-2) then
         xm= 0.5d0*(gx(i0) +gx(i0+1))
      else
         xm= xmax
      end if
!
      if(xx0.le.xm) then
         ptx(i)= i0
      else
         i0= i0 +1
         ptx(i)= min(i0, mx-1)
      end if
  120 continue
      end do
!
!
      do i= 1,10
      ptx(i1+i)= mx-1
      ptx(-i  )= 0
!
      if(i1+i.gt.3000) then
        if(io_pe.eq.1) then
          OPEN (unit=11,file=praefixc//'.06'//suffix2,             &
                status='unknown',position='append',form='formatted')
!
          write(11,*) ' # ptx(i) i >3000 ... xx0=',xx0
          close(11)
        end if
!
        stop
      end if
      end do
!
      if(io_pe.eq.1) then
        OPEN (unit=11,file=praefixc//'.06'//suffix2,             &
              status='unknown',position='append',form='formatted')
!
        write(11,*) '## ptx(imax): imax=',i1
        do i= i1-10,i1
        write(11,130) i,ptx(i)
  130   format('i=',i5,' ptx=',i7)
        end do
!
        close(11)
      end if
!
! Y: from j= -10 to j= i1+10
!
      hhy0 = 0.3d0 * hy0
!
      do j= -10,i1+10
      pty(j)= ptx(j)
      end do
!
! -----------------------------------------
! Z:
      k0 = 0
      ptz(0)= 0
      hhz0 = 0.3d0 * hz0
!
      do k= 1,4000    !<<-- 4000
      zz0= zmin + hhz0*k
      if(zz0.ge.zmax) go to 140
!
      k1= k
      if(k0.le.mz-2) then
        zm= 0.5d0*(gz(k0) +gz(k0+1))
      else
        zm= zmax
      end if
!
      if(zz0.le.zm) then
        ptz(k)= k0
      else
        k0= k0 +1
        ptz(k)= min(k0, mz-1)
      end if
!
  140 continue
      end do
!
!  ptz() from k1+1 to k1+10, and -10 to -1
      do k= 1,10
      ptz(k1+k)= mz-1
      ptz(-k  )= 0
!
      if(k1+k.gt.4000) then
        if(io_pe.eq.1) then
          OPEN (unit=11,file=praefixc//'.06'//suffix2,             &
                status='unknown',position='append',form='formatted')
!
          write(11,*) ' # ptz(k) k >4000 ... zz0=',zz0
          close(11)
        end if
!
        stop
      end if
      end do
!
      if(io_pe.eq.1) then
        OPEN (unit=11,file=praefixc//'.06'//suffix2,             &
              status='unknown',position='append',form='formatted')
!
        write(11,*) '## ptz(kmax): kmax=',k1
        do k= k1-10,k1
        write(11,150) k,ptz(k)
  150   format('k=',i5,' ptz=',i7)
        end do
!
        close(11)
      end if
!
! -----------------------------------------
!*  For particles only (must be uniform).
! -----------------------------------------
!* Non-periodic case.
!* Fine mesh
!
      hxi= 0.999999d0/hhx0
      hyi= 0.999999d0/hhy0
      hzi= 0.999999d0/hhz0
!
      if(ipar.eq.1) then
        OPEN (unit=11,file=praefixc//'.06'//suffix2,             &
              status='unknown',position='append',form='formatted')
!
        write(11,170) hxi,hyi,hzi
  170   format(/,' hxi, hyi, hzi=',3f12.7)
        close(11)
      end if
!
!-----------------------------
!* 1. DNA (Polyelectrolyte)
!-----------------------------
!* (Adenosine) .... m= 312/e 
!       (5+5)*C + 5*O +P +(7+5)*H +5*N; C=12, N=14, O=16, P=30
!  Al  .... 26/3e,  Ca .... 40/2e,
!  Na  .... 22/e ,  K  .... 38/e
!  Cl  .... 34/-e
! -----------------------
!
      ww1  = 38.d0  ! K,  2.27 Ang, K+ 2.98 Ang
      ww2  = 34.d0  ! Cl, 1.75 Ang, Cl- 1.81 Ang
!
      a_phos = 4.1d0   ! core 3.2 Ang, poc 2.6 Ang
      a_sugar= 3.6d0   ! core 2.2 Ang, csc 5.0 Ang
      a_base1= 4.3d0
!
      bond_ps= 4.2d0 +0.1d0 ! P-sugar ring (3.5+5.0)/2
      bond_ss= 5.3d0 +0.1d0 ! 4.5d0 +0.3d0 ! two sugar rings
!
!* Lennard-Jones parameters
!
      KJ  = 1.d+10              ! 1 KJ= 1.d10, per erg
      KCal= 4.1868d0 *KJ        ! 4.18 J/Cal
      mol = 6.0220d23
      kbT = 1.3807d-16 * 300.   ! at 300 Kelvin
!
      epsLJ = 1.d+3  *kJ/mol     ! KJ/mol for water= 6.6e-3 E.volt
      eps_K =  418.8 *kJ/mol       ! 0.0148 * 4.1868d0  ! 4.18 x Kcal/mol= KJ/mol
      eps_Cl= 1251.2 *kJ/mol     ! 0.1064 * 4.1868d0  ! 
!
!* Represented in normalized values in kbT
      epsLJ = epsLJ  *1.d-3/kbT
      eps_K = eps_K  *1.d-3/kbT
      eps_Cl= eps_Cl *1.d-3/kbT
!
      if(io_pe.eq.1) then
        OPEN (unit=11,file=praefixc//'.06'//suffix2,             &
              status='unknown',position='append',form='formatted')
!
        write(11,203) epsLJ,eps_K,eps_Cl
  203   format(/,'epsLJ, eps_K, eps_Cl /kT=',3f8.3)
        close(11)
      end if
!
!
!* Restart fixed - 7/15/2004
! -----------------------------
      if(kstart.ne.0) return
! -----------------------------
!
      if(np.ne.0) then
!     ++++++++++++++++
!**
      if(io_pe.eq.1) then
        OPEN (unit=11,file=praefixc//'.06'//suffix2,             &
              status='unknown',position='append',form='formatted')
!
        write(11,*) 'Load polyelectrolyte... np=',np 
        write(11,*) '    np= n_p*n_lp (READ_CONF)'
        close(11)
      end if
!*
      q_PE= 0.d0
!
      do i= 1,np
!   PO_4 and sugar ring, mod(i,2) <- 1,2, 3,4,... 
!
      if(mod(i,2).eq.1) then   
        ch(i)= -1.d0                 ! PO_4
        if(ifbase.eq.2) ch(i)= 0.d0  !  neutral chain 
!
        am(i)=  94.d0/wwat
        ag(i)=  a_phos/2     ! Diameter= 4.1 Ang
        ep(i)=  epsLJ
      else
!            ag(i)= 4.1/(2*1.5)= 1.37
!
        ch(i)= 0.d0                  ! Sugar ring
        am(i)= 218.d0/wwat
        ag(i)=  a_sugar/2    ! Diameter= 3.6 Ang
        ep(i)=  epsLJ
      end if
!
      q_PE= q_PE +ch(i)
!
      if(mod(i,2).eq.1) then
        rod_len(i)= -rod_leng/2 +ag(i)*(i-1)  !<- chain
      else
        rod_len(i)= rod_len(i-1)
      end if
!
      if(io_pe.eq.1) then
        OPEN (unit=11,file=praefixc//'.06'//suffix2,             &
              status='unknown',position='append',form='formatted')
!
        write(11,720) i,ch(i),am(i),ag(i),rod_len(i)
  720   format(' i=',i3,' ch,am,ag,rod_len=',4f8.2)
!
        close(11)
      end if
      end do
!*
      if(io_pe.eq.1) then
        OPEN (unit=11,file=praefixc//'.06'//suffix2,             &
              status='unknown',position='append',form='formatted')
!
        write(11,*)
        write(11,*) 'DNA (PO4 and sugar ring)  np=1 to ',np
!
        do i= 1,np
        write(11,210) i,ch(i),am(i),ag(i),ep(i)
  210   format('i=',i5,'  ch=',f6.1,'  am=',f6.1,'  ag,ep=',f6.1,1pd12.3)
        end do
!
        close(11)
      end if
!*
!
      do i= 1,np
      vx(i)= 0
      vy(i)= 0
      vz(i)= 0
      end do
!
!  Double or single stranded DNA
      if(mod(nseg,2).eq.0) then
        n_rodp= (np+1)/2  ! Double stranded DNA
      else
!
        n_rodp= np        ! Single stranded DNA
      end if
!
!  Quaternions
      if(kstart.eq.0) then
        xgc = 0.d0
        ygc = 0.d0
        zgc = 0.d0
!
        vxg = 0.d0
        vyg = 0.d0
        vzg = 0.d0
!
        q1 = 0.d0
        q2 = 0.d0
        q3 = 0.d0
        q4 = 1.d0
!
        angx= 0.d0
        angy= 0.d0
        angz= 0.d0
      end if
!
! -------------------------------------------
!*  Inertial moment: Get initial position
! -------------------------------------------
!   n_rodp= (np+1)/2 (double) or np (single) 
      airc= (n_rodp+1)/2.
!     Rhelix= Rmac    !<- Radius for monomer helix, READ_CONF 
!
!* Charge blobs of the rod: doxc-dozc are defined here !
!
      do i= 1,n_rodp
      th= pi2*(i-1)/(n_rodp -1)
!
      doxc(i)= Rhelix*sin(th)
      doyc(i)= Rhelix*cos(th)
      dozc(i)= rod_leng*(i -airc)/float(n_rodp-1)
!
      if(mod(nseg,2).eq.0) then
        doxc(n_rodp+i) = Rhelix*sin(th +pi)
        doyc(n_rodp+i) = Rhelix*cos(th +pi)
        dozc(n_rodp+i) = dozc(i)
      end if
      end do
!
! Also define the rod-axis (rod core)
!
      if(mod(nseg,2).eq.0) then
!
        do i= 1,n_rodp     ! double stranded
        dox(i)= 0.d0
        doy(i)= 0.d0
        doz(i)= (rod_leng -Rmac)*(i -airc)/float(n_rodp-1)
        end do
      else
!
        do i= 1,n_rodp     ! single stranded - no rod core
        dox(i)= doxc(i)    !    beads_on_the_chain
        doy(i)= doyc(i)
        doz(i)= dozc(i)
        end do
      end if
!
      if(io_pe.eq.1) then
        OPEN (unit=11,file=praefixc//'.06'//suffix2,             &
              status='unknown',position='append',form='formatted')
!
        write(11,*)
        write(11,*) 'dox and doxc'
        write(11,*) '  single(1)/double(2) helix at n_rodp/pi2...'
        write(11,*) 'nseg, n_rodp, Rhelix=',nseg,n_rodp,Rhelix
        write(11,*) ' Note: n_rodp= (np+1)/2 (double) or np (single)'
!
        do i= 1,n_rodp     !<-- n_rodp=24
        write(11,230) i,doxc(i),doyc(i),dozc(i),dox(i),doy(i),doz(i)
  230   format('i=',i4,1p3d12.3,2x,3d12.3)    !<-- only i=1-24 !
        end do                                !    i= 25-36 by Counterions
      end if                                  !    so, it is 0.0
!
      ipx= 0.d0
      ipy= 0.d0
      ipz= 0.d0
!
      do i= 1,np
      ipx= ipx + am(1)*(doy(i)**2+doz(i)**2)
      ipy= ipy + am(1)*(doz(i)**2+dox(i)**2)
      ipz= ipz + am(1)*(dox(i)**2+doy(i)**2)
      end do
!
! ip around the major axis
!
      ipz1= np*am(1)*Rmac**2

      if(io_pe.eq.1) then
        OPEN (unit=11,file=praefixc//'.06'//suffix2,             &
              status='unknown',position='append',form='formatted')
!
        write(11,240) ipx,ipy,ipz,ipz1
  240   format(' Inertial moments (ipx,ipy,ipz,ipz1): ',4f11.3,/)
        close(11)
      end if
!
!* Definition of X-Z: must be always performed.
! ----------------------------------------------
      call get_rod (xg,yg,zg,rxg,ryg,rzg,np)  !<-- xg-zg are defined
! ----------------------------------------------
!
!*  np,nseg: defined in /READ_CONF/
! ---------------------
      nps = np/nseg     !<-- nseg=1,nps=1
! ---------------------
      if(io_pe.eq.1) then
        OPEN (unit=11,file=praefixc//'.06'//suffix2,             &
              status='unknown',position='append',form='formatted')
!
        write(11,*)
        write(11,*) 'Polyelectrolyte chains ...'
        write(11,*) '  nseg, nops=',nseg,nps
        write(11,*)
        close(11) 
      end if
!
      if(nseg.gt.0) then
!*
        do k= 1,nseg
        nsg(k)= nps*(k-1)
        end do
!
        nsg(nseg+1)= nsg(nseg) + nps 
!
        if(io_pe.eq.1) then
          OPEN (unit=11,file=praefixc//'.06'//suffix2,             &
                status='unknown',position='append',form='formatted')
!
          write(11,*)
          write(11,*) 'Polyelectrolyte chains ...'
          write(11,*) '  nseg, nops=',nseg,nps
          write(11,*) ' chain=',k,' starts at i=',nsg(k)+1
!
          write(11,*) '  Last chain ends at nsg(nseg+1)=',nsg(nseg+1)
          write(11,*)                                   ! nsg()=24 
          close(11)
        end if
      end if
!**
      end if  !<- L.4490-4725
!
! --------------------------------------------------------
!*  Additional monomers attached to sugar ring - A,G,T,C
! --------------------------------------------------------
!   np is updated by attaching the base AGCT
!
      np00 = np  
      nnb = 0
!
      if(ifbase.eq.0) go to 290  !<-- Without sugar rings
!     *************************
!
      do i= 1,np00 
      if(abs(ch(i)).gt.0.1d0) go to 250  ! loop - only neutral
      if(mod(i,2).eq.1) go to 250        ! skip (was not neutral)   
!
      nnb = nnb +1
!
      ch(np00 +nnb)= 0.d0                ! base: AGCT,heavy,LJ
      am(np00 +nnb)= 130.d0/wwat 
      ag(np00 +nnb)= a_base1/2   ! Diameter
      ep(np00 +nnb)= epsLJ
!
      ist1(nnb)= i
      ist2(nnb)= np00 +nnb
!
  250 continue
      end do
!
!    ++++++++++++++++
      np = np00 +nnb  ! Now np= 24+12= 36, the bases are 12 monomers
!    ++++++++++++++++
!
      do i= np00,np 
      vx(i)= 0
      vy(i)= 0
      vz(i)= 0
      end do
!
        if(io_pe.eq.1) then
          OPEN (unit=11,file=praefixc//'.06'//suffix2,             &
                status='unknown',position='append',form='formatted')
          write(11,*) '273 continue'
          close(11)
        end if
!
! --------------------
      do jj= 1,nnb  !<-- do 270
!
      i= ist1(jj)   ! ist1(1),...
      k= ist2(jj)   ! ist2(1),...
!
      ss = ag(i) +ag(k) 
!     ll= 0
!
  273 continue  !ll= ll +1
      ph = 2*pi*ranff(0.d0)
      th =   pi*ranff(0.d0)
      ph3= pi2*(i-1)/(n_rodp -1)
!
      x0= xg(i) +ss*(cos(ph3 +ph))
      y0= yg(i) +ss*(sin(ph3 +ph))
!
      z0= zg(i)             !+ss*sin(th)
!         +++++
!
!* Avoid the separator plate.
!
      Hpore2= 0.5d0*Hpore +ag(k)
!
      if(abs(z0).le.Hpore2) then
        rr= sqrt(x0**2 +y0**2)
        if(rr.gt. Rpore-ag(k)) go to 273 
      end if 
!
        if(io_pe.eq.1) then
          OPEN (unit=11,file=praefixc//'.06'//suffix2,             &
                status='unknown',position='append',form='formatted')
          write(11,275) i,k,ph,x0,y0,z0
  275     format('i,k,ph,x0-z0=',2i4,4f8.2)
          close(11)
        end if
!
!* Rigid rod polyelectoryte         
      if(ifrodco.eq.1) then
        do j= 1,n_rodp
        rr= sqrt((x0-rxg(j))**2 +(y0-ryg(j))**2 +(z0-rzg(j))**2) ! rxg is used
        if(rr.le. ag(k)+Rmac) go to 273                     ! Rmac here...
        end do 
      end if
!
      do j= 1,i-1  !  j=1 to (i-1)
      rr= sqrt((x0-xg(j))**2 +(y0-yg(j))**2 +(z0-zg(j))**2) 
      if(rr.le. ag(k) +ag(j)) go to 273
      end do  
!
      xg(k)= x0    ! branch of AGCT
      yg(k)= y0
      zg(k)= z0
      end do       !<- do 270
  290 continue
!
!  Generate new np
      if(io_pe.eq.1) then
        OPEN (unit=11,file=praefixc//'.06'//suffix2,             &
              status='unknown',position='append',form='formatted')
!
        write(11,*)
        write(11,*) 'Generate i= 1,np'
!
        do i= 1,np  !<-- complete i=1,36
        write(11,295) i,xg(i),yg(i),zg(i),ch(i),am(i)
  295   format('i=',i5,' xg,yg,zg=',3f8.1,'  ch,am=',f8.1,1pd12.3)
        end do
! 
        write(11,*)
        close(11)
      end if
!
!
!**********************************************
!*  Give charges of counter/co-ions           *
!**********************************************
!* For (Zcp:Zcn) ions 
!          ch(i).lt.0: q_PE= -12, 24 monomers
!
!              ****
      Nzi= abs(q_PE)/Zcp +0.1  ! Polyelectrolyte q_PE, ch() < 0
      Nci= np + Nzi            ! positive Nzi=12 + np=24 -> Nci= 36
!
      Ncc_pore= 28 
      Ncc_wide= nq -Nzi -Ncc_pore  ! 300 -12 -28= 260
!     +++++++++++++++++++++++++++
!
      n1= 0
      n2= 0
      n3= 0
      n4= 0
! 
      do i= np+1,np+nq 
!     ----------------
      if(i.gt.np .and. i.le.np+Nzi) then  ! Counterions to charge neutrality
          n1= n1 +1                       ! they are Nci=12
!
!     ww1  = 38.d0  ! K,  2.27 Ang, K+ 2.98 Ang
!     ww2  = 34.d0  ! Cl, 1.75 Ang, Cl- 1.81 Ang
          ch(i)= Zcp                      ! due to PE=-12; that is i=37-48
          am(i)= 38.d0/wwat
          ag(i)= 2.98 ! acount
          ep(i)= 4.34 ! eV  eps_K
!
      else if(i.ge.np+Nzi+1 .and. i.le.np+Nzi+nq/2) then
          n3= n3 +1
!
          ch(i)=  Zcp
          am(i)= 38.d0/wwat
          ag(i)= 2.98 ! acount
          ep(i)= 4.34 ! eV
!
      else if(i.ge.np+Nzi+nq/2+1 .and. i.le.np+nq) then
          n4= n4 +1
!
          ch(i)=  Zcn
          am(i)= 34.d0/wwat
          ag(i)= 1.81 ! acoion
          ep(i)= eps_Cl
      end if
      end do
!
      cci= 0.
      do i= 1,np+nq
      cci= cci +ch(i)    !<- 0. ?
      end do
!
      if(io_pe.eq.1) then
        OPEN (unit=11,file=praefixc//'.06'//suffix2,             &
              status='unknown',position='append',form='formatted')
!
        write(11,350) q_PE,Rpore,Hpore,Zcp,Zcn,cci,n1,n2
        write(11,351) acount,acoion
        write(11,352) epsLJ,eps_K,eps_Cl
        write(11,*)
  350   format(' q_PE, Rpore(radius), Hpore(height) =',3f7.2,/,  &
               ' Zcp, Zcn, cci(PE+count/co)= ',3f7.2,/,          &
               ' Number of counterions/coions= ',2i5)
  351   format(' Size of counterion (Ang)=',f7.2, &
               '         coion      (Ang)=',f7.2)
  352   format(' LJ potential(water, Ka, Cl)=',3f7.2)
!
        close(11)
      end if
!
!     call shuffl_ions (ch,np,nq)
!
!--------------------------------------------------
!* Position of counterions and coions: np+1,np+nq.
!--------------------------------------------------
!* i= np+1,Nci are all neutralizing counterions
!  which must be located in the pore      7/15/2004
!
      n0= 0
!           ++++ ++++++
      do i= np+1,np+Nzi   !! Counterions only near PE
      j= 3*(i-np) -2
!                                 ! rr in Rpore-ag(i) and 0 
 330  rr = 0.7*(Rpore-ag(i))*ranff(0.d0) 
      th = pi2*ranff(0.d0)
!
      x0 = rr*cos(th)
      y0 = rr*sin(th) 
      z0 = Hpore*(ranff(0.d0) -0.5d0)
! 
      if( sqrt((x0-xg(j))**2 +(y0-yg(j))**2 +(z0-zg(j))**2) &
                                 .le. ag(i)+ag(j) ) go to 330 !<--
      n0= n0 +1
      ch(i)= Zcp
! 
      xg(i)= x0       ! Counterions near the polyelectroryte 
      yg(i)= y0
      zg(i)= z0
!
        if(io_pe.eq.1) then
          OPEN (unit=11,file=praefixc//'.06'//suffix2,             &
                status='unknown',position='append',form='formatted')
!
          write(11,*) 'Pore region:  i,n0(+)=',i,n0
          close(11)
        end if
      end do
!
!* Balance the number of ions
!
      n3= 0
      n4= 0
!           ++++++++ +++++++++++
      do i= np+Nzi+1,np+Nzi+nq/2
!
  460 z0= zmin +(zmax-zmin)*ranff(0.d0)
!
      if(abs(z0).lt.Hpore2) then   !<- z0 < Hpore/2
        th = pi2*ranff(0.d0)
        x0= (Rpore-ag(i))*cos(th)
        y0= (Rpore-ag(i))*sin(th)
      else
        x0 = xmin +(xmax-xmin)*ranff(0.d0) 
        y0 = ymin +(ymax-ymin)*ranff(0.d0) 
      end if
!
      do j= 1,i-1 
      rr= sqrt((x0-xg(j))**2 +(y0-yg(j))**2 +(z0-zg(j))**2)
      if(rr.le. ag(i)+ag(j)) go to 460 
      end do  
!
      n3= n3 +1
      ch(i)= Zcp
!
      xg(i)= x0 
      yg(i)= y0
      zg(i)= z0
!
        if(io_pe.eq.1) then
          OPEN (unit=11,file=praefixc//'.06'//suffix2,             &
                status='unknown',position='append',form='formatted')
!
          write(11,*) 'Counterion: i,n3(+), n4(-)=',i,n3,n4
          close(11)
        end if
      end do 
!
!           +++++++++++++++++
      do i= np+Nzi+nq/2,np+nq
!
  465 z0= zmin +(zmax-zmin)*ranff(0.d0) 
!
      if(abs(z0).lt.Hpore2) then
        th = pi2*ranff(0.d0)
        x0= (Rpore-ag(i))*cos(th)
        y0= (Rpore-ag(i))*sin(th)
      else
        x0 = xmin +(xmax-xmin)*ranff(0.d0) 
        y0 = ymin +(ymax-ymin)*ranff(0.d0) 
      end if
!
      do j= 1,i-1 
      rr= sqrt((x0-xg(j))**2 +(y0-yg(j))**2 +(z0-zg(j))**2)
      if(rr.le. ag(i)+ag(j)) go to 465 
      end do 
!
      n4= n4 +1
      ch(i)= Zcn
!
      xg(i)= x0 
      yg(i)= y0
      zg(i)= z0
!
        if(io_pe.eq.1) then
          OPEN (unit=11,file=praefixc//'.06'//suffix2,             &
                status='unknown',position='append',form='formatted')
!
          write(11,*) 'Coion:      i,n3(+), n4(-)=',i,n3,n4
          close(11)
        end if
      end do 
!
!
      if(io_pe.eq.1) then
        OPEN (unit=11,file=praefixc//'.06'//suffix2,             &
              status='unknown',position='append',form='formatted')
!
        write(11,*)
        write(11,*) 'Number of ions in the pore...'
        write(11,*) '  DNA: PO_4 =',q_PE
        write(11,*) '       Sugar and AGCT=',2.d0*abs(q_PE)
        write(11,*) '  Counterions and coions in pore and wide regions '
        write(11,*) '  (+)=',n0
        write(11,*) '  (+)=',n3
        write(11,*) '  (-)=',n4
        write(11,*) ' np+nq = ',np+nq 
        write(11,*)
!
        close(11)
      end if
!
! Swap Zcp and Zcn
!
      nn3= np+Nzi+nq/2  ! max number of counterions
      nn4= np+Nzi+nq
!
!     do k= 1,3
!     do k= 1,5  !<-- Look at numbers !
      do k= 1,6  !<-- Look at numbers !
!
      n3= n3 -1
      n4= n4 +1
!
! Border of Zcp and Zcn
      nn3= nn3 -1
      nn4= nn4 +1
!
!     ww1  = 38.d0  ! K,  2.27 Ang, K+ 2.98 Ang
!     ww2  = 34.d0  ! Cl, 1.75 Ang, Cl- 1.81 Ang
      ch(nn3)=  Zcp
      am(nn3)= 38.d0/wwat
      ag(nn3)= acount
      ep(nn3)= eps_K
!
      ch(nn4)=  Zcn
      am(nn4)= 34.d0/wwat
      ag(nn4)= acoion
      ep(nn4)= eps_Cl
      end do
!
      if(io_pe.eq.1) then
        OPEN (unit=11,file=praefixc//'.06'//suffix2,             &
              status='unknown',position='append',form='formatted')
!
        write(11,*)
        write(11,*) 'Swap Zcp and Zcn...'
        write(11,*) '  (+)=',n3
        write(11,*) '  (-)=',n4
        write(11,*) ' Zcp (last) =',nn3
        write(11,*) ' Zcn (first)=',nn4
        write(11,*) ' np+nq = ',np+nq 
        write(11,*)
!
        close(11)
      end if
!
!-----------------------------------------
!* 3. Solvent particles
!-----------------------------------------
!   Must avoid macroions. One molecule in every 3 Angstroms. 
!     if a= 2 Ang, then one H2O/(1.5a)**3.
!     for L= 32a, one has 21 H2O 
!
      ll= np + nq  !!<- Starting number of solvent, ll
!         **   **
!
      if(io_pe.eq.1) then
        OPEN (unit=11,file=praefixc//'.06'//suffix2,             &
              status='unknown',position='append',form='formatted')
!
        write(11,*)
        write(11,*) '< After swapping Zcp and Zcn >'
        write(11,*) 'Number of ions in the pore...'
        write(11,*) '  DNA: PO_4 (- charge)=',q_PE
        write(11,*) '       sugar and AGCT =',2.d0*abs(q_PE)
        write(11,*) 'Counterions and coions in pore and wide regions '
        write(11,*) '  (+)=',n0
        write(11,*) '  (+)=',n3
        write(11,*) '  (-)=',n4
        write(11,*) 'np+nq = ',np+nq 
        write(11,*)
!
        close(11)
      end if
!
      xcent= xleng/2.d0
      ycent= yleng/2.d0
      zcent= zleng/2.d0
!
      awat_diam= 3.10d0   !<- diameter
      kl= 24
      km= 24
      kn= 50
      istp= 0
!
      do k= 1,kn  ! 17,34     !<-- (i,j,k) are used for cells 
      do j= 1,km
      do i= 1,kl
      istp= istp +1
!
      x0= awat_diam*(i -0.5d0) -xcent  !<- middle center of water
      y0= awat_diam*(j -0.5d0) -ycent
      z0= awat_diam*(k -0.5d0) -zcent
!
      if(i.eq.1 .and. j.eq.1 .and. (k.le.5 .or. k.ge.kn-4)) then
      if(io_pe.eq.1) then
        OPEN (unit=11,file=praefixc//'.06'//suffix2,             &
              status='unknown',position='append',form='formatted')
!
        write(11,701) istp,i,j,k,x0,y0,z0
  701   format(' Solvent: istp,kl,km,kn,x0-z0=',i8,3i5,3f7.2)
        close(11)
      end if
      end if
!
!----------------------------
      Hpore2= 0.5d0*Hpore 
!
      if(z0.lt.zmin .or. z0.gt.zmax) go to 530
!
!* Over the pore
      if(abs(z0).gt.Hpore2) then   !<- z0 > Hpore/2
!
        do l= 1,ll  !! Nci+1,ll
        rr= sqrt((x0-xg(l))**2 +(y0-yg(l))**2 +(z0-zg(l))**2)
        if(rr.le.(ag(l)+awat)) go to 530
        end do 
!
!* Within the pore
      else if(abs(z0).le.Hpore2) then  !<- z0 < Hpore/2
!
        rr= sqrt(x0**2 +y0**2)
        if(rr.gt.Rpore-awat) go to 530 
!
        do l= 1,ll  !! Nci+1,ll   !<-- ions in existing PE or ions
        rr= sqrt((x0-xg(l))**2 +(y0-yg(l))**2 +(z0-zg(l))**2)
        if(rr.le.(ag(l)+awat)) go to 530
        end do 
      end if
!
!* Rod Polyelectrolyte
      if(ifrodco.eq.1) then
        do l= 1,n_rodp
        rr= sqrt((x0-rxg(l))**2 +(y0-ryg(l))**2 +(z0-rzg(l))**2) ! rxg is used
        if(rr.le.(Rmac +awat)) go to 530                         ! Rmac here...
        end do 
      end if
!
      ll= ll +1         !!<- next number of water, ll
!
      if(ll.gt.npqr0) then
        if(io_pe.eq.1) then
          OPEN (unit=11,file=praefixc//'.06'//suffix2,             &
                status='unknown',position='append',form='formatted')
!
          write(11,*) ' Solvent: the l is greater than npqr0..., l=',ll
          close(11)
        end if
!
        stop
      end if
!
      xg(ll)= x0   !<-- Solvent particles
      yg(ll)= y0
      zg(ll)= z0
!
      am(ll)= wwat
      ag(ll)= awat
      ch(ll)= 0.d0
!
  530 continue  !<- Occupied in (i,j,k) 
      end do
      end do
      end do
!
!***************************************************
      nCLp = np +nq      ! PE+AGCT, Counter/Coions
!                   
      npqr= ll           ! the total particles 
      nr  = npqr - nCLp  ! Water
!***************************************************
!
!* Fold back to (-L/2,L/2)
!
      do i= 1,npqr0
      xg(i)= xg(i) - NINT(xg(i)/xleng)*xleng
      yg(i)= yg(i) - NINT(yg(i)/yleng)*yleng
      zg(i)= zg(i) - NINT(zg(i)/zleng)*zleng
      end do
!
!  Generate nq
      if(io_pe.eq.1) then
        OPEN (unit=11,file=praefixc//'.06'//suffix2,             &
              status='unknown',position='append',form='formatted')
!
        write(11,*)
        write(11,*) 'Number of particles (Macro, Counter/Coions)'
        write(11,*) '  np, nq, nr=',np,nq,nr
        write(11,*)
!
        do i= np+1,np+nq
        write(11,540) i,xg(i),yg(i),zg(i),ch(i),am(i)
  540   format('i=',i5,' xg,yg,zg=',3f8.1,'  ch,am=',f8.1,1pd12.3)
        end do
!
        write(11,*) '--- Subroutine init (end) ------------------------'
        write(11,*)
        close(11)
      end if
!
!***************************
!* Dielectric constant     *
!***************************
!* used in /diel/ to calculate dielectric constant
!
!     call water_dens (xg,yg,zg,ipar,nCLp,npqr)
!
      return
      end subroutine init
!
!
!--------------------------------------------------------------------
      subroutine get_rod (xg,yg,zg,rxg,ryg,rzg,np)
!--------------------------------------------------------------------
      use, intrinsic :: iso_c_binding 
      implicit none
!*
      include   'paramAPGa.h'
!
      real(C_DOUBLE),dimension(npqr0) :: xg,yg,zg
      real(C_DOUBLE),dimension(np0) ::  rxg,ryg,rzg
      real(C_DOUBLE),dimension(np0) ::  dox,doy,doz,doxc,doyc,dozc
      common/macroin1/ dox,doy,doz,doxc,doyc,dozc
!
      real(C_DOUBLE) xgc,ygc,zgc,vxg,vyg,vzg,rod_leng,Rmac,Rhelix, &
                     q1,q2,q3,q4,angx,angy,angz,ipx,ipy,ipz,       &
                     rxx,rxy,rxz,ryx,ryy,ryz,rzx,rzy,rzz
!
      integer(C_INT) np,n_rodp,i
      common/macroion/ xgc,ygc,zgc,vxg,vyg,vzg,rod_leng,Rmac,  &
                       Rhelix,n_rodp
      common/macroiov/ q1,q2,q3,q4,angx,angy,angz,ipx,ipy,ipz
!
      integer(C_INT) io_pe
      common/sub_proc/ io_pe
!
!* First define the r-matrix by quaternion
!
      rxx= -q1**2 +q2**2 -q3**2 +q4**2
      rxy=  2.d0*(q3*q4 -q1*q2)
      rxz=  2.d0*(q2*q3 +q1*q4)
      ryx= -2.d0*(q1*q2 +q3*q4)
      ryy= q1**2 -q2**2 -q3**2 +q4**2
      ryz=  2.d0*(q2*q4 -q1*q3)
      rzx=  2.d0*(q2*q3 -q1*q4)
      rzy= -2.d0*(q1*q3 +q2*q4)
      rzz= -q1**2 -q2**2 +q3**2 +q4**2
!
      if(io_pe.eq.1) then
        OPEN (unit=11,file=praefixc//'.06'//suffix2,             &
              status='unknown',position='append',form='formatted')
!
        write(11,*) 'Get_rod: matrix'
        write(11,910) rxx,rxy,rxz,ryx,ryy,ryz,rzx,rzy,rzz
  910   format('rxx,rxy,rxz=',1p3d12.3,/, &
               'ryx,ryy,ryz=',3d12.3,/,   &
               'rzx,rzy,rzz=',3d12.3)
      end if
!
!* r-inverse is used for the Coulomb force.
!  Define xg,yp,zg by doxc - dozc
!
      do i= 1,np
      xg(i)= xgc +rxx*doxc(i) +ryx*doyc(i) +rzx*dozc(i)
      yg(i)= ygc +rxy*doxc(i) +ryy*doyc(i) +rzy*dozc(i)
      zg(i)= zgc +rxz*doxc(i) +ryz*doyc(i) +rzz*dozc(i)
      end do
!
      if(io_pe.eq.1) then
        OPEN (unit=11,file=praefixc//'.06'//suffix2,             &
              status='unknown',position='append',form='formatted')
        write(11,*) 'Get_rod  np=',np           !<-- np= 25-36 ???
        write(11,*) 'xgc,ygc,zgc=',xgc,ygc,zgc  !   Counterion, not a rod !
!
        do i=1,np
        write(11,920) i,doxc(i),doyc(i),dozc(i)  !<-- i= 25-36  ???
  920   format('i=',i4,' doxc,doyc,dozc=',1p3d12.3)
        end do
!
        do i= 1,np
        write(11,930) i,xg(i),yg(i),zg(i)
  930   format('i=',i4,' xg,yg,zg=',1p3d12.3)
        end do
        close(11)
      end if
!
!* Rod-axis center for the LJ forces
!  Defined here by dox-doz: rxg - rzg.
!
      do i= 1,n_rodp
      rxg(i)= xgc +rxx*dox(i) +ryx*doy(i) +rzx*doz(i)
      ryg(i)= ygc +rxy*dox(i) +ryy*doy(i) +rzy*doz(i)
      rzg(i)= zgc +rxz*dox(i) +ryz*doy(i) +rzz*doz(i)
      end do
!
      if(io_pe.eq.1) then
        OPEN (unit=11,file=praefixc//'.06'//suffix2,             &
              status='unknown',position='append',form='formatted')
        write(11,*) 'n_rodp=',n_rodp
!
        do i= 1,n_rodp
        write(11,940) i,xg(i),yg(i),zg(i)
  940   format('get_rod i=',i4,' rxg..=',1p3d12.3)
        end do
        close(11)
      end if
!
      return
      end subroutine get_rod
!
!                                        ***********
!---------------------------------------------------------------------
      subroutine get_torque (tqx,tqy,tqz,fcx,fcy,fcz,ftx,fty,ftz,np)
!---------------------------------------------------------------------
      use, intrinsic :: iso_c_binding 
      implicit none
!*
      include   'paramAPGa.h'
!
      real(C_DOUBLE),dimension(npq0) :: fcx,fcy,fcz
      real(C_DOUBLE),dimension(np0) ::  ftx,fty,ftz 
      real(C_DOUBLE),dimension(np0) ::  dox,doy,doz,doxc,doyc,dozc
      common/macroin1/ dox,doy,doz,doxc,doyc,dozc
!
      real(C_DOUBLE) tqx,tqy,tqz,                                &
                     q1,q2,q3,q4,angx,angy,angz,ipx,ipy,ipz,     &
                     rxx,rxy,rxz,ryx,ryy,ryz,rzx,rzy,rzz,        &
                     sxo,syo,szo,xmax,ymax,zmax,xmin,ymin,zmin , &
                     xgc,ygc,zgc,vxg,vyg,vzg,rod_leng,Rmac,Rhelix
      real(C_float)  phi,tht,dtwr1,dtwr2,dtwr3
      integer(C_INT) np,n_rodp,i
!
      common/parm3/ xmax,ymax,zmax,xmin,ymin,zmin
      common/parm4/ phi,tht,dtwr1,dtwr2,dtwr3
      common/macroion/ xgc,ygc,zgc,vxg,vyg,vzg,rod_leng,Rmac,  &
                       Rhelix,n_rodp
      common/macroiov/ q1,q2,q3,q4,angx,angy,angz,ipx,ipy,ipz
!
!--------------------------------------------------------------------
!*  Calculate site vectors with q, since the use of sx= xo(i) -x(i)
!   makes a problem because of the folding-back procedure.
!--------------------------------------------------------------------
!
      rxx= -q1**2 +q2**2 -q3**2 +q4**2
      rxy=  2.d0*(q3*q4 -q1*q2)
      rxz=  2.d0*(q2*q3 +q1*q4)
      ryx= -2.d0*(q1*q2 +q3*q4)
      ryy= q1**2 -q2**2 -q3**2 +q4**2
      ryz=  2.d0*(q2*q4 -q1*q3)
      rzx=  2.d0*(q2*q3 -q1*q4)
      rzy= -2.d0*(q1*q3 +q2*q4)
      rzz= -q1**2 -q2**2 +q3**2 +q4**2
!
!* The Coulomb force (site vector from the gc)
!
      tqx= 0.d0
      tqy= 0.d0
      tqz= 0.d0

      do i= 1,np
      sxo=  rxx*doxc(i) +ryx*doyc(i) +rzx*dozc(i)
      syo=  rxy*doxc(i) +ryy*doyc(i) +rzy*dozc(i)
      szo=  rxz*doxc(i) +ryz*doyc(i) +rzz*dozc(i)
!
      tqx=  tqx + syo*fcz(i) -szo*fcy(i)
      tqy=  tqy + szo*fcx(i) -sxo*fcz(i)
      tqz=  tqz + sxo*fcy(i) -syo*fcx(i)
      end do
!
!* The LJ force from the rod-axis
!
      do i= 1,n_rodp
      sxo=  rxx*dox(i) +ryx*doy(i) +rzx*doz(i)
      syo=  rxy*dox(i) +ryy*doy(i) +rzy*doz(i)
      szo=  rxz*dox(i) +ryz*doy(i) +rzz*doz(i)
!
      tqx=  tqx + syo*ftz(i) -szo*fty(i)
      tqy=  tqy + szo*ftx(i) -sxo*ftz(i)
      tqz=  tqz + sxo*fty(i) -syo*ftx(i)
      end do
!
      return
      end subroutine get_torque 
!
!
!---------------------------------------------------------------
      subroutine shuffl_ions (ch,np,nq)
!---------------------------------------------------------------
!*  Shuffle the order of ions (not chained beads) 
!
      use, intrinsic :: iso_c_binding 
      implicit none
!
      include   'paramAPGa.h'
!
      integer(C_INT) np,nq,kshuf,k,i1,i2
      real(C_DOUBLE) ranff,ch(npqr0),sch
!
      kshuf= 3*nq*ranff(0.d0) +1.001
!
      do k= 1,kshuf
      i1= np+1 +mod(int(nq*ranff(0.d0)),nq)
      i2= np+1 +mod(int(nq*ranff(0.d0)),nq)
      if(ch(i1)*ch(i2).gt.0.) go to 130
!
      sch   = ch(i1)
      ch(i1)= ch(i2)
      ch(i2)= sch
!
  130 continue
      end do
!
      return
      end subroutine shuffl_ions 
!
!
!---------------------------------------------------------------
      subroutine ggauss
!---------------------------------------------------------------
!*  Create the Gauss distribution
!
      use, intrinsic :: iso_c_binding 
      implicit none
!
      integer(C_int) i,j,k,k2,ns
      real(C_DOUBLE) fv,vv0,dv,vv,s,sdv,fun
      common/gaus1/  fv(101),vv0,dv
!
      fv(1)= 0
!
      vv0= -3.d0
      dv= 2.d0*abs(vv0)/100.d0
!
      vv= vv0
      do j= 1,100
      s= 0
      ns= 1000
      k2= ns/2
      sdv= dv/float(ns)
!
      do k= 1,k2
      vv= vv +2.d0*sdv
      s= s +4.d0*fun(vv-sdv) +2.d0*fun(vv)
      end do
!
      s= (s +4.d0*fun(vv+sdv) +fun(vv+2.d0*sdv))*sdv/3.d0
      fv(j+1)= fv(j)+s
      end do
!
      do i= 1,101
      fv(i)= fv(i)/fv(101)
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
      implicit none
!
      integer(C_INT) k,k2
      real(C_DOUBLE) dgaus2,vmax,eps,x2,y1,y2
      real(C_DOUBLE) fv,vv0,dv,ranff
      common/gaus1/ fv(101),vv0,dv
!
      fv(1)= 0
      vv0= -3.d0
      dv= 2.d0*abs(vv0)/100.d0
      eps= ranff(0.d0)
!
      do k= 1,101
      k2= k
      if(fv(k).gt.eps) go to 200
      end do
!
  200 y1= fv(k2-1)
      y2= fv(k2)
      x2= (eps-y2)/(y2-y1) +k2
      dgaus2= vmax*(vv0 +dv*(x2-1.d0))
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
!
      real(C_DOUBLE) fun,v
      fun= exp(-v**2/2.d0)
!
      return
      end function fun
!
!
!-----------------------------------------------------------
      subroutine vdistr (vx,vy,vz,np,nq,nCLp,nr,npqr)
!-----------------------------------------------------------
      use, intrinsic :: iso_c_binding 
      implicit none
!
      include   'paramAPGa.h'
!
      real(C_DOUBLE),dimension(npqr0) :: vx,vy,vz
      integer(C_INT) np,nq,nCLp,nr,npqr
!
!     character*1  char(8**3)
      real(C_DOUBLE) pi,dt,axi,Gamma,rbmax,vth,tmax, &
                     wwat,awat
      common/parm2/  pi,dt,axi,Gamma,rbmax,vth,tmax
      common/waterm/ wwat,awat
!
      integer(C_INT) i,ix,iy,iz,k,ILN
      real(C_float)  aiv,xsc(101),fvx(101),fvy(101),fvz(101), &
                     fmax1,fmax2,fmax3,fmin1,fmin2,fmin3
      real(C_DOUBLE) vmax1,vmax2
!
!*******************************
!*  Ions only.                 *
!*******************************
!     am(1)=  94.d0/wwat
!     ww1  = 38.d0  ! K,  2.27 Ang, K+ 2.98 Ang
!
!     vmax1= 10.d0*vth/sqrt(38.d0/wwat)  !<- mass of K(+)
      vmax1=  6.d0*vth/sqrt(18.d0/wwat)  !<- water
      vmax2=  6.d0*vth/sqrt(18.d0/wwat)  !<- water
!
!* (1) Coulomb particles.
      aiv= 50./vmax1
!
      do k= 1,101
      xsc(k)= (k -51)/aiv
      fvx(k)= 0
      fvy(k)= 0
      fvz(k)= 0
      end do
!
      do i= np+1,nCLp
      ix= aiv*vx(i) +50.5
      iy= aiv*vy(i) +50.5
      iz= aiv*vz(i) +50.5
!
      if(abs(ix-51).gt.50)  go to 230
        fvx(ix)= fvx(ix) +1.
!
  230 if(abs(iy-51).gt.50)  go to 250
        fvy(iy)= fvy(iy) +1.
!
  250 if(abs(iz-51).gt.50)  go to 200
        fvz(iz)= fvz(iz) +1.
  200 continue
      end do
!
      call lplmax (fvx,fmax1,fmin1,101)
      call lplmax (fvy,fmax2,fmin2,101)
      call lplmax (fvz,fmax3,fmin3,101)
!
      ILN= 1
      call hplot1 (2,4,101,xsc,fvx,fmax1,fmin1,iln,'fvx(K +)',8, &
                   '   vx   ',8,'        ',8)
      call hplot1 (2,5,101,xsc,fvy,fmax2,fmin2,iln,'fvy(K +)',8, &
                   '   vy   ',8,'        ',8)
      call hplot1 (2,6,101,xsc,fvz,fmax3,fmin3,iln,'fvz(K +)',8, &
                   '   vz   ',8,'        ',8)
!
!* (2) Water particles.
      aiv= 50.d0/vmax2
!
      do k= 1,101
      xsc(k)= (k -51)/aiv
      fvx(k)= 0
      fvy(k)= 0
      fvz(k)= 0
      end do
!
      do i= nCLp+1,npqr
      ix= aiv*vx(i) +50.5
      iy= aiv*vy(i) +50.5
      iz= aiv*vz(i) +50.5
!
      if(abs(ix-51).gt.50)  go to 430
        fvx(ix)= fvx(ix) +1.
!
  430 if(abs(iy-51).gt.50)  go to 450
        fvy(iy)= fvy(iy) +1.
!
  450 if(abs(iz-51).gt.50)  go to 400
        fvz(iz)= fvz(iz) +1.
  400 continue
      end do
!
      call lplmax (fvx,fmax1,fmin1,101)
      call lplmax (fvy,fmax2,fmin2,101)
      call lplmax (fvz,fmax3,fmin3,101)
!
      ILN= 1
      call hplot1 (3,4,101,xsc,fvx,fmax1,fmin1,iln,'fvx(wat)',8, &
                   '   vx   ',8,'        ',8)
      call hplot1 (3,5,101,xsc,fvy,fmax2,fmin2,iln,'fvy(wat)',8, &
                   '   vy   ',8,'        ',8)
      call hplot1 (3,6,101,xsc,fvz,fmax3,fmin3,iln,'fvz(wat)',8, &
                   '   vz   ',8,'        ',8)
!---------------------
      call chart
!---------------------
!
      return
      end subroutine vdistr 
!
!
!------------------------------------------------------
      subroutine averg1 (q,qav,is)
!------------------------------------------------------
      use, intrinsic :: iso_c_binding 
      implicit none
!
      integer(C_INT) is,i
      real(C_float)  q(is),qav
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
      subroutine clocks (walltime,size,cl_first)
!------------------------------------------------------
!*  Measure both cpu and elapsed times 
!
      use, intrinsic :: iso_c_binding 
      implicit none
!
      include   'mpif.h'
      include   'paramAPGa.h'
!
      integer(C_INT) io_pe
      common/sub_proc/ io_pe
!
      integer(c_int) size,cl_first,MPIerror 
      real(c_double) walltime,walltime0,buffer1,buffer2
      save       walltime0
!
      if(cl_first.eq.1) then
        walltime0= mpi_wtime() 
!
        buffer1= walltime0
        call mpi_allreduce (buffer1,buffer2,1,mpi_double_precision, &
                            mpi_sum,MPI_COMM_WORLD,MPIerror)
        walltime0 = buffer2/size
      end if
!
      walltime = mpi_wtime() - walltime0
!
      buffer1= walltime
      call mpi_allreduce (buffer1,buffer2,1,mpi_double_precision, &
                          mpi_sum,MPI_COMM_WORLD,MPIerror)
      walltime = buffer2/size
!
      return
      end subroutine clocks
!
!
!***********************************************************************
!*  The Poisson solver which uses /cressl/ or /bcgsts/ matrix solver   *
!***********************************************************************
!-----------------------------------------------------------------------
      subroutine poissn (rho8,pot8,ndim,itrmax,iterp,ipar)
!-----------------------------------------------------------------------
!***********************************************************************
!*    << Poisson solver in 3-D : cresmd method >>                      *
!*                                                                     *
!*        Laplacian*pot(x,y,z) = rho(x,y,z).                           *
!*                                                                     *
!*   pot(i,j,k)... solution in (x,y,z) space.                          *
!*              p comes with a guess.                  <----- input    *
!*                                                                     *
!*   rho(i,j,k)....... source term.                    <----- input    *
!*                                                                     *
!***********************************************************************
      use, intrinsic :: iso_c_binding 
      implicit none
!*
      include    'mpif.h'
      include    'paramAPGa.h'
!
      integer(C_INT) ndim,itrmax,iterp,ipar 
!     parameter     (nob3=13,iblk3=1)           ! paramAPGa.h
!
      real(C_DOUBLE),dimension(0:mxyz-1) :: rho8,pot8,xx8
!
      real(C_DOUBLE) Vtop0,Vbot0,Vtop,Vbot,t_recyc
      integer(C_INT) fixedbc_x,fixedbc_y,itermax,filtx,filty,filtz
      common/espot/  Vtop0,Vbot0,Vtop,Vbot,fixedbc_x,fixedbc_y,  &
                     filtx,filty,filtz
!*---------------------------------------------------------------------
      real(C_DOUBLE)  e_c_s,e_poisn,e_c_r,e_lj,e_elas
      common /energy/ e_c_s,e_poisn,e_c_r,e_lj,e_elas
!*----------------------------------------------------------------------
      real(C_DOUBLE) aa,ss,xx,rpr,eps
      integer(C_INT) na,nobx,ipr,itrm,ik,i,j,k,ijk
!
      common/cresm1/  aa(0:mxyz-1,nob3)
      common/cresma/  na(0:mxyz-1,nob3)
      common/cresmb/  nobx
      common/cresmc/  xx(0:mxyz-1)
      common/cresm2/  ss(0:mxyz-1)
      common/cresm3/  ipr(10),rpr(10)
!
      real(C_DOUBLE) ww(mxyz,6),wp(mxyz)     ! uses wp(1:mxyz)
      integer(C_INT) iw(mxyz),mxyz7,mxyz8,nob7
!
      integer(C_INT) iwrt1,iwrt2,iwrt3,iwrt4,ierror
      real(C_DOUBLE) sres,avex,rsdl,conv
      common/iotim/  iwrt1,iwrt2,iwrt3,iwrt4
      common/crconv/ sres,avex,rsdl,conv,ierror
!
      integer(C_INT) io_pe,cnt_recv,disp_recv,i0,i2,i3,i4,i00, &
                     cnt_send2,cnt_recv2,disp_recv2
      common/sub_proc/ io_pe
      common/dat_typ2/ i3(30),i4(30),cnt_recv2(30),disp_recv2(30)
!
      logical     first /.true./
!*
      do ik= i3(ipar),i4(ipar)  ! the index i3(1)= 0
      xx(ik)= pot8(ik)
      ss(ik)= rho8(ik)
      end do
!                                       Define the CFP matrix
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!  If dielectric constant remains unchanged, the matrix elements 
!  do not change, otherwise aa() must be defined each time. Boundary 
!  source terms need to be defined.
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      if(first) then
        call escof3 (aa,xx,ss,na,ndim)
        first= .false.
      else
        call bound_s (xx,ss)
      end if
!
      nobx = nob3
      itrm = itrmax
      eps  = 1.0d-5
!
      ipr(1) = itrm
      rpr(1) = eps
!
!                     Solve the equation by ipar=1,num_proc
!***********************************************************************
      call cresmd (xx,ss,ierror,ipar)
!
!     mxyz7= mxyz
!     nob7 = nob3
!
!     call bcgsts (aa,mxyz7,mxyz7,nob7,na,ss,xx,ipr,rpr, &
!                  iw,ww,wp,ierror)
!***********************************************************************
!* Define packed arrays for mpi_allgatherv.
!
      if(num_proc.eq.1) then
        do ik= 0,mxyz-1
        pot8(ik)= xx(ik)
        end do
!
      else
        i00= i3(ipar)
!
        do i= i3(ipar),i4(ipar)  !<- partial solution in i3,i4
        xx8(i-i00)= xx(i)
        end do
!
        cnt_send2= i4(ipar) -i3(ipar) +1
        call mpi_allgatherv (xx8, cnt_send2,           mpi_real8, &
                             pot8,cnt_recv2,disp_recv2,mpi_real8, &
                             mpi_comm_world,ierror)
      end if
! 
      iterp = ipr(2)
      rsdl  = rpr(2)
!
!     if(ierr.ne.0) then
!     if(io_pe.eq.1 .and. iwrt1.eq.0) then
!       OPEN (unit=11,file=praefixc//'.06'//suffix2,             &
!             status='unknown',position='append',form='formatted')
!
!       write(11,600) ierror,iterp,avex,rsdl
! 600   format('#(bcg-ps) ierr,iter,avex,rsdl= ',2i6,2f13.6)
!       close(11)
!     end if
!     end if
!
      return
      end subroutine poissn 
!
!
!-----------------------------------------------------------------------
      subroutine escof3 (aa,pot8,rho8,na,ndim)
!----------------------------------------------------------------------
      use, intrinsic :: iso_c_binding 
      implicit none
!
      include    'paramAPGa.h'
!
!     integer(C_INT) nob3     !<- paramAPGa.h
!     parameter     (nob3=13)
!
!*---------------------------------------------------------------------
      real(C_DOUBLE) aa(0:mxyz-1,nob3)
      integer(C_INT) na(0:mxyz-1,nob3)
!
      real(C_DOUBLE),dimension(0:mx-1,0:my-1,0:mz-1) :: pot8,rho8,dec2
      integer(C_INT),dimension(nob3) :: lai,laj,lak
!
      real(C_DOUBLE) Vtop0,Vbot0,Vtop,Vbot,t_recyc,ca(nob3)
      integer(C_INT) fixedbc_x,fixedbc_y,itermax,filtx,filty,filtz
      common/espot/  Vtop0,Vbot0,Vtop,Vbot,fixedbc_x,fixedbc_y,  &
                     itermax,filtx,filty,filtz
!*---------------------------------------------------------------------
!
      integer(C_INT) pxr,pxc,pxl,pyr,pyc,pyl,pzr,pzc,pzl
      common/ptable/ pxr(-1:mx+1),pxc(-1:mx+1),pxl(-1:mx+1), &
                     pyr(-1:my+1),pyc(-1:my+1),pyl(-1:my+1), &
                     pzr(-1:mz+1),pzc(-1:mz+1),pzl(-1:mz+1)
!
      integer(C_INT) io_pe
      common/sub_proc/ io_pe
!
      integer(C_INT) ndim,i,j,k,m,ii,jj,kk,ijk,lxyz
!
      real(C_DOUBLE) gx,gy,gz,ghx,ghy,ghz,ghx2,ghy2,ghz2, &
                     ghxsq,ghysq,ghzsq,hxi,hyi,hzi,diel
      integer(C_INT) ptx,pty,ptz
      common/xtable/ gx(0:mx),ghx(0:mx),ghx2(0:mx),ghxsq(0:mx), &
                     gy(0:my),ghy(0:my),ghy2(0:my),ghysq(0:my), &
                     gz(0:mz),ghz(0:mz),ghz2(0:mz),ghzsq(0:mz), &
                     hxi,hyi,hzi
      common/xtabl2/ ptx(-10:3000),pty(-10:3000),ptz(-10:4000)
!
      if(ndim.eq.2) then
        if(io_pe.eq.1) then
        OPEN (unit=11,file=praefixc//'.06'//suffix2,             &
              status='unknown',position='append',form='formatted')
!
        write(11,*) 'escof3: rewrite the coefficients for 2-D case'
        write(11,*) '   ca2 =  0.'
        write(11,*) '   ca4 =  2./hysq'
        close(11)
        end if
!
        stop
      end if
!
!*                       !!<-- escof3, i,j periodic/k bound
      do i= 0,mx-1 
!
      if(i.eq.0) then
        do j= 0,my-1 
        do k= 0,my-1 
!
        lai(1) =  0 
        laj(1) =  j  
        lak(1) =  k 
         ca(1) =  1.d0 
!
        lai(2)=   1
        laj(2)=   j
        lak(2)=   k
         ca(2)=   1.d0 
!*
        lxyz= 0 +mx*(j +my*k)
!
        aa(lxyz,1) = ca(1) 
        na(lxyz,1) = 0 +mx*(j +my*k)
! 
        aa(lxyz,2) = ca(2) 
        na(lxyz,2) = 1 +mx*(j +my*k) 
! 
        do m= 3,nob3
        aa(lxyz,m)= 0
        na(lxyz,m)= lxyz 
        end do
        end do
        end do
!
      else if(i.gt.0 .and. i.lt.mx-1) then
!
!***
      do j= 0,my-1
!
      if(j.eq.0) then
        do k= 0,mz-1
        lai(1) =  i 
        laj(1) =  0  
        lak(1) =  k 
         ca(1) =  1.d0 
!
        lai(2)=   i
        laj(2)=   1
        lak(2)=   k
         ca(2)=   1.d0 
!*
        lxyz= i +mx*(0 +my*k)
!
        aa(lxyz,1) = ca(1) 
        na(lxyz,1) = i +mx*(0 +my*k)
! 
        aa(lxyz,2) = ca(2) 
        na(lxyz,2) = i +mx*(1 +my*k) 
! 
        do m= 3,nob3
        aa(lxyz,m)= 0
        na(lxyz,m)= lxyz 
        end do
        end do
!
      else if(j.gt.0 .and. j.lt.my-1) then
!
!***
      do k= 0,mz-1
      if(k.eq.0) then
!
!* (i,j-1,k)
        lai(1) =  i      ! number i=0: lai(1,0)=0
        laj(1) =  j      !  na(,)=0    laj(1,0)=0
        lak(1) =  0      !             lak(1,0)=0
         ca(1) =  1.d0   !            1:nob2 0:mx-1
!
!* (i,j+1,k)
        lai(2)=   i
        laj(2)=   j
        lak(2)=   1
         ca(2)=   1.d0 
!*
        lxyz= i +mx*(j +my*0)
!
        aa(lxyz,1) = ca(1)  ! aa(1,m=1)
        na(lxyz,1) = i +mx*(j +my*0)  ! na(1,m=1)=0
!            l    m           +          +
        aa(lxyz,2) = ca(2)  ! aa(1,m=2) <- lxyz
        na(lxyz,2) = i +mx*(j +my*1) 
! 
        do m= 3,nob3
        aa(lxyz,m)= 0
        na(lxyz,m)= lxyz 
        end do
!***
      else if(k.gt.0 .and. k.lt.mz-1) then
!
!* (i,j,k-1)
      lai(1)=   i
      laj(1)=   j
      lak(1)=   k-1
       ca(1)=   1.d0/ghzsq(k)
!
!* (j-1,k)
      lai(2)=   i
      laj(2)=   j-1
      lak(2)=   k
       ca(2)=   1.d0/ghysq(j)
!
!* (i-1,j,k)
      lai(3)=   i-1
      laj(3)=   j
      lak(3)=   k
       ca(3)=   1.d0/ghxsq(i)
!
!* (i,j,k)
      lai(4)=   i
      laj(4)=   j
      lak(4)=   k
       ca(4)=   -2.d0*(1/ghxsq(i) +1/ghysq(j) +1/ghzsq(k))
!
!* (i+1,j,k)
      lai(5)=   i+1
      laj(5)=   j
      lak(5)=   k
       ca(5)=   1.d0/ghxsq(i)
!
!* (i,j+1,k)
      lai(6)=   i
      laj(6)=   j+1
      lak(6)=   k
       ca(6)=   1.d0/ghysq(j)
!
!* (i,j,k+1)
      lai(7)=   i
      laj(7)=   j
      lak(7)=   k+1
       ca(7)=   1.d0/ghzsq(k)
!
!* (deps/dx)/eps
      lai(8)=   i-1
      laj(8)=   j
      lak(8)=   k
       ca(8)=  -(dec2(i+1,j,k)-dec2(i-1,j,k))/(ghx2(i)**2*dec2(i,j,k))
!
      lai(9)=   i+1
      laj(9)=   j
      lak(9)=   k
       ca(9)=  (dec2(i+1,j,k)-dec2(i-1,j,k))/(ghx2(i)**2*dec2(i,j,k)) 
!
!* (deps/dy)/eps
      lai(10)=  i
      laj(10)=  j-1
      lak(10)=  k
       ca(10)= -(dec2(i,j+1,k)-dec2(i,j-1,k))/(ghy2(j)**2*dec2(i,j,k))
!
      lai(11)=  i
      laj(11)=  j+1
      lak(11)=  k
       ca(11)=  (dec2(i,j+1,k)-dec2(i,j-1,k))/(ghy2(j)**2*dec2(i,j,k))
!
!* (deps/dz)/eps
      lai(12)=  i
      laj(12)=  j
      lak(12)=  k-1
       ca(12)= -(dec2(i,j,k+1)-dec2(i,j,k-1))/(ghz2(k)**2*dec2(i,j,k))
!
      lai(13)=  i
      laj(13)=  j
      lak(13)=  k+1
       ca(13)=  (dec2(i,j,k+1)-dec2(i,j,k-1))/(ghz2(k)**2*dec2(i,j,k))
!
!
      lxyz= i +mx*(j +my*k)
!
      do m= 1,nob3
      ii= lai(m)      ! if periodic  pxc()
      jj= laj(m)
      kk= lak(m)      !<-- bound 
!
      aa(lxyz,m) = ca(m) 
      na(lxyz,m) = ii +mx*(jj +my*kk) 
      end do 
!***
!
      else if(k.eq.mz-1) then
!
!* (i,j-1,k)
        lai(1)=   i
        laj(1)=   j
        lak(1)=   mz-2
         ca(1)=   1.d0
!
!* (j+1,k)
        lai(2)=   i
        laj(2)=   j
        lak(2)=   mz-1
         ca(2)=   1.d0
!
        lxyz= i +mx*(j +my*(mz-1)) 
!
        aa(lxyz,1) = ca(1) 
        na(lxyz,1) = i +mx*(j +my*(mz-2)) 
! 
        aa(lxyz,2) = ca(2) 
        na(lxyz,2) = i +mx*(j +my*(mz-1)) 
! 
        do m= 3,nob3
        aa(lxyz,m)= 0
        na(lxyz,m)= lxyz 
        end do
!
      end if
      end do  !<- end of k
!**
!
      else if(j.eq.my-1) then
        do k= 0,mz-1
        lai(1) =  i 
        laj(1) =  my-2  
        lak(1) =  k 
         ca(1) =  1.d0 
!
        lai(2)=   i
        laj(2)=   my-1 
        lak(2)=   k
         ca(2)=   1.d0 
!*
        lxyz= i +mx*(my-1 +my*k)
!
        aa(lxyz,1) = ca(1) 
        na(lxyz,1) = i +mx*(my-2 +my*k)
! 
        aa(lxyz,2) = ca(2) 
        na(lxyz,2) = i +mx*(my-1 +my*k) 
! 
        do m= 3,nob3
        aa(lxyz,m)= 0
        na(lxyz,m)= lxyz 
        end do
        end do
      end if
      end do  ! end of j
!**
!
      else if(i.eq.mx-1) then
        do j= 0,my-1 
        do k= 0,my-1 
!
        lai(1) =  mx-2 
        laj(1) =  j  
        lak(1) =  k 
         ca(1) =  1.d0 
!
        lai(2)=   mx-1
        laj(2)=   j
        lak(2)=   k
         ca(2)=   1.d0 
!*
        lxyz= mx-1 +mx*(j +my*k)
!
        aa(lxyz,1) = ca(1) 
        na(lxyz,1) = mx-2 +mx*(j +my*k)
!            l    m           +          +
        aa(lxyz,2) = ca(2) 
        na(lxyz,2) = mx-1 +mx*(j +my*k) 
! 
        do m= 3,nob3
        aa(lxyz,m)= 0
        na(lxyz,m)= lxyz 
        end do
!
        end do
        end do
      end if 
      end do  !<- end of i
!
      return
      end subroutine escof3
!
!
!-----------------------------------------------------------------------
      subroutine bound_s (pot8,rho8)
!----------------------------------------------------------------------
! ******************************
! * Boundary values for source *
! ******************************
      use, intrinsic :: iso_c_binding 
      implicit none
!
      include    'paramAPGa.h'
!     integer(C_INT) nob3   !<- paramAPGa.h
!
      real(C_DOUBLE),dimension(0:mx-1,0:my-1,0:mz-1) :: pot8,rho8
      real(C_DOUBLE) Vtop0,Vbot0,Vtop,Vbot,t_recyc,ca(nob3)
      integer(C_INT) fixedbc_x,fixedbc_y,itermax,filtx,filty,filtz, &
                     i,j,k
      common/espot/  Vtop0,Vbot0,Vtop,Vbot,fixedbc_x,fixedbc_y,     &
                     itermax,filtx,filty,filtz
!*---------------------------------------------------------------------
!
      do k= 0,mz-1
      do j= 0,my-1
      do i= 0,mx-1
!
! if x= xmin or xmax
      if(i.eq.0 .or. i.eq.mx-1) then
        pot8(i,j,k)= 0.d0
        rho8(i,j,k)= 0.d0
      end if
!
! if y= ymin or ymax
      if(j.eq.0 .or. j.eq.my-1) then
        pot8(i,j,k)= 0.d0
        rho8(i,j,k)= 0.d0
      end if
!
! if z= zmin or zmax
      if(k.eq.0 .or. k.eq.mz-1) then
        if(k.eq.   0) pot8(i,j,k)= Vbot
        if(k.eq.mz-1) pot8(i,j,k)= Vtop
!
        rho8(i,j,k)= 0.d0
      end if
!
      end do
      end do
      end do
!
      return
      end subroutine bound_s 
!
!
!***********************************************************************
!*    Bi-conjugate gradient method: non-block type.                    *
!***********************************************************************
!----------------------------------------------------------------------
      subroutine bcgsts (aa,la,n,ma,ja,b,u,ipr,rpr,iw,w,wp,ierr)
!----------------------------------------------------------------------
!
!    subroutine name : bcgsts
!
!    purpose :
!       this subroutine solves non symmetric linear systems
!       by bi-cgstab method with scalar precondtioning.
!
!    format :
!                 -- -- - -- -- -   --- ---
!     call bcgsts(aa,la,n,ma,ja,b,u,ipr,rpr,iw,w,wp,ierr)
!                                 - --- ---      ----
!
!    arguments :
!
!    aa(la,ma) (real*8) ------------ coefficient matrix
!    ja(la,ma) (integer*4) --------- column no. table
!                                    ja(i,jj) is column no.of aa(i,jj)
!    n (integer*4) ----------------- vector size
!    b(n) (real*8) ----------------- right-hand side's vector
!    u(n) (real*8) ----------------- solution vector
!    ipr(10) (integer*4) ----------- integer parameter
!                   ipr(1) (in) .... maximum number of iteration
!                   ipr(2) (out) ... final iteration number
!                   ipr(3) (in) .... type of preconditioning
!    rpr(10) (real*8) -------------- real parameter
!                   rpr(1) (in) .... admissible residual norm
!                   rpr(2) (out) ... final residual norm
!                   rpr(3) (out) ... cpu time
!    iw(n) (integer*4) ------------- integer work area
!    w(n,6),wp(n) (real*8) --------- real work area
!                   w(*,1) ......... residual vector
!    ierr (integer*4) -------------- return code
!                   .eq. 0 ......... normal termination
!                   .ge. 1000 ...... warning: useing default values
!                   .ge. 3000 ...... fatal error: abnormal termination
!----------------------------------------------------------------------
      use, intrinsic :: iso_c_binding 
      implicit none 
!
      real(C_DOUBLE) aa(la,ma),b(n),u(n),rpr(10),w(n,6),wp(n)
      integer(C_INT) ja(la,ma),ipr(10),iw(n),la,ma,n,ierr
!
      integer(C_INT) itrm,iprc,iblk,i,iflg,jm
      real(C_DOUBLE) eps
!
!     cl_first= 2
!     call clocks(c0,size,cl_first)
      ierr = 0
!
!                                      Parameter check
!----------------------------------------------------------------------
      call wwstsp (la,n,ipr,rpr,ierr)
!----------------------------------------------------------------------
!
!                                      Parameter set
      itrm = ipr(1)
      iprc = ipr(3)
!     iblk = ipr(4)
      eps  = rpr(1)
!                                      Matrix check
!----------------------------------------------------------------------
      call wwstsc (aa,la,n,ma,ja,iw,jm,iflg,ierr)
!----------------------------------------------------------------------
!
      if ( ierr.ge.3000 ) goto 4000
!
!                                      Matrix sort if necessary
!
      if ( iprc.eq.1 .and. iflg.eq.1 ) then
         call wwstss (aa,la,ma,n,jm,ja,iw)
         ierr = 1200
      endif
!
!                                      Preconditioning
!----------------------------------------------------------------------
      call wwstsi (aa,la,ma,n,jm,ja,iw,iprc,wp,ierr)
!----------------------------------------------------------------------
!
!                                      Make initial value
!
      if ( ierr.eq.1000 .or. ierr.eq.1010 ) then
         do i=1,n
           u(i) = b(i)
         end do
      else
         do i=1,n
           u(i) = b(i)/aa(i,iw(i))
         end do
      end if
!
!                                      Bi-cgstab's iteration
!----------------------------------------------------------------------
      call wwstsk (aa,la,ma,n,jm,ja,b,u,itrm,iprc,eps,iw,wp,     &
                   w(1,1),w(1,2),w(1,3),w(1,4),w(1,5),w(1,6),ierr)
!----------------------------------------------------------------------
!
 4000 continue
!     cl_first= 2
!     call clocks (c1,size,cl_first)
!
      ipr(2) = itrm
      rpr(2) = eps
      rpr(3) = 0  ! c1 - c0
!
      return
      end subroutine bcgsts 
!
!
!----------------------------------------------------------------------
      subroutine wwstsk (aa,la,ma,n,jm,ja,b,u,itrm,iprc,eps,jd,wp,  &
                         r,p,q,s,t,v,ierr)
!----------------------------------------------------------------------
      use, intrinsic :: iso_c_binding 
      implicit none 
!
      include  'paramAPGa.h'
!
      real(C_DOUBLE) aa(la,ma),b(n),u(n),wp(n),      &
                     r(n),p(n),q(n),s(n),t(n),v(n),eps
      integer(C_INT) ja(la,ma),jd(n),la,ma,n,jm,prc,itrm,iprc,ierr
!
      integer(C_INT) iwrt1,iwrt2,iwrt3,iwrt4,io_pe
      common/iotim/  iwrt1,iwrt2,iwrt3,iwrt4
      common/sub_proc/ io_pe
!
      integer(C_INT) itr,i
      real(C_DOUBLE) bb,bnrm,alpa,beta,omeg,qrofn,qrofp, &
                     qv,tr,tt,rr,uu,rsdl
!
!                                  bi-cgstab method iteration
!                                       initialize
!
!                                            itr = 0
      itr = 0
!
!                                            bb = (b,b)
      bb = 0.d0
!
      do i=1,n
        bb = bb + b(i)*b(i)
      end do
!
      if ( bb.lt.1.d-72 ) then
         ierr = 1300
         do 110 i=1,n
            u(i) = 0.d0
  110    continue
         goto 2000
      elseif ( bb.gt.1.d+72 ) then
         ierr = 3500
         goto 2000
      endif
!
!                                            bnrm = sqrt(bb)
      bnrm = sqrt(bb)
!
!                                            r = a*u
      call wwstsm (aa,la,ma,n,jm,ja,u,r)
!
!                                            r = b - r
      do i=1,n
        r(i) = b(i) - r(i)
      end do
!                                            q = r
      do i=1,n
        q(i) = r(i)
      end do
!
!                                            p = v = 0
      do i=1,n
        p(i) = 0.d0
        v(i) = 0.d0
      end do
!
!                                            alpa = omeg = qrofn = 1
      alpa = 1.d0
      omeg = 1.d0
      qrofn = 1.d0
!
!                                       iteration
!*-----------------------------------------***--------------------------
 1000 continue
!                                            itr = itr + 1
      itr = itr + 1
!
!                                            qrofp = qrofn
      qrofp = qrofn
!
!                                            qrofn = (q,r)
      qrofn = 0.d0
      do i=1,n
        qrofn = qrofn + q(i)*r(i)
      end do
!
      if ( abs(qrofn).lt.1.d-72 ) then
         ierr = 3600
      endif
!                                            beta = (qrofn*alpa)
!                                                  /(qrofp*omgp)
      if ( abs(qrofp).lt.1.d-72 ) then
         ierr = 3601
         beta = alpa/omeg
      else
!
         beta = (qrofn*alpa)/(qrofp*omeg)
      end if
!
!                                            v = r - beta*omeg*v
      do i=1,n
        v(i) = r(i) - beta*omeg*v(i)
      end do
!
!                                            s = m^(-1)*v
      call wwstsj (aa,la,ma,n,jm,ja,jd,wp,iprc,v,s)
!
!                                            p = s + beta*p
      do i=1,n
         p(i) = s(i) + beta*p(i)
      end do
!                                            v = a*p
      call wwstsm (aa,la,ma,n,jm,ja,p,v)
!
!                                            qv = (q,v)
      qv = 0.d0
      do i=1,n
        qv = qv + q(i)*v(i)
      end do
!
      if ( abs(qv).lt.1.d-72 ) then
         ierr = 3700
         goto 2000
      endif
!                                            alpa = qrofn/qv
!
      alpa = qrofn/qv
!
!                                            r = r - alpa*v
      do i=1,n
        r(i) = r(i) - alpa*v(i)
      end do
!
!                                            s = m^(-1)*r
      call wwstsj (aa,la,ma,n,jm,ja,jd,wp,iprc,r,s)
!
!                                            t = a*s
      call wwstsm (aa,la,ma,n,jm,ja,s,t)
!
!                                            tr = (t,r) , tt = (t,t)
      tr = 0.d0
      tt = 0.d0
      do i=1,n
        tr = tr + t(i)*r(i)
        tt = tt + t(i)*t(i)
      end do
!
      if ( tt.lt.1.d-72 ) then
         ierr = 3700
         goto 2000
      endif
!
!                                            omeg = tr/tt
      omeg = tr/tt
!                                            u = u + alpa*p + omeg*s
!                                            r = r - omeg*t
      do i=1,n
        u(i) = u(i) + alpa*p(i) + omeg*s(i)
        r(i) = r(i) - omeg*t(i)
      end do
!                                            rr = (r,r)
      rr = 0.d0
      uu = 0.d0
!
      do i=1,n
        rr = rr + r(i)*r(i)
        uu = uu + u(i)*u(i)
      end do
!
      if ( rr.gt.1.d+72 ) then
        ierr = 3800
        goto 2000
      endif
!                                            rsdl = sqrt(rr)/bnrm
      rsdl = sqrt(rr)/(sqrt(uu) +1.e-15) ! bnrm
      uu   = sqrt(uu/float(n))
!
!     if (mod(itr,50).eq.0 .or. itr.lt.50) then
      if(io_pe.eq.1 .and. mod(itr,50).eq.1) then
        OPEN (unit=11,file=praefixc//'.06'//suffix2,             &
              status='unknown',position='append',form='formatted')
!
        write(11,900) itr,rsdl,uu
  900   format('(s) itr= ',i4,'  rsdl= ',1pe12.3,'  <x>= ',e12.3)
        close(11)
      end if
!     end if
!
      if ( rsdl.lt.eps ) goto 2000
      if ( itr.le.itrm ) then
         goto 1000
      else
         ierr = 2000
      endif
!*-----------------------------------------***--------------------------
!                                       end of iteration
 2000 continue
!
      if(iwrt1.eq.0 .and. io_pe.eq.1) then
        OPEN (unit=11,file=praefixc//'.06'//suffix2,             &
              status='unknown',position='append',form='formatted')
!
        write(11,600) itr,ierr,rsdl,uu
  600   format('#(s) itr, ierr= ',2i5,', rsdl, <x>= ',1p2e12.4)
        close(11)
      end if
!
      itrm = itr
      eps  = rsdl
!
      return
      end subroutine wwstsk 
!
!
!----------------------------------------------------------------------
      subroutine wwstsp (la,n,ipr,rpr,ierr)
!----------------------------------------------------------------------
      use, intrinsic :: iso_c_binding 
      implicit none 
!
      integer(C_INT) ipr(10),n,la,ierr
      real(C_DOUBLE) rpr(10)
!
      if ( n.le.0 ) then
         ierr = 3000
         return
      endif
!
      if ( la.lt.n ) then
         ierr = 3100
         return
      endif
!
      if ( ipr(1).le.0 ) then
         ipr(1) = 5*int(sqrt(real(n)))
         ierr   = 1100
      endif
!
      if ( ipr(3).lt.0 .or. ipr(3).gt.1 ) then
         ipr(3) = 0
         ierr   = 1110
      endif
!
      if ( rpr(1).lt.1.d-72 ) then
         rpr(1) = 1.d-12
         ierr   = 1150
      endif
!
      return
      end subroutine wwstsp 
!
!
!----------------------------------------------------------------------
      subroutine wwstsc (aa,la,n,ma,ja,jd,jm,iflg,ierr)
!----------------------------------------------------------------------
      use, intrinsic :: iso_c_binding 
      implicit none 
!
      real(C_DOUBLE) aa(la,ma)
      integer(C_INT) ja(la,ma),jd(n),la,ma,n,jm,iflg,ierr
!
      integer(C_INT) icnt,i,j,jj
      real(C_DOUBLE) sum
!
      iflg = 0
      jm = 0
      do 100 i=1,n
         icnt = 0
         jd(i) = 0
         sum = 0.d0
         do 110 jj=1,ma
            j = ja(i,jj)
            if ( j.lt.1 .or. j.gt.n ) then
               aa(i,jj) = 0.d0
               ja(i,jj) = i
               goto 110
            endif
            if ( j.eq.i .and. abs(aa(i,jj)).lt.1.d-72 ) goto 110
!
            jm = max(jm,jj)
!
            if ( j.eq.i ) then
               jd(i) = jj
               icnt  = icnt + 1
               goto 110
            endif
!
            sum = sum + abs(aa(i,jj))
            if ( jd(i).eq.0 .and. j.ge.i ) iflg = 1
            if ( jd(i).ne.0 .and. j.le.i ) iflg = 1
  110    continue
!
         if ( icnt.eq.0 ) then
            ierr = 1000
         elseif ( icnt.gt.1 ) then
            ierr = 4000
            return
         elseif ( abs(aa(i,jd(i))).lt.1.d-1*sum ) then
            ierr = 1010
         endif
  100 continue
!
      return
      end subroutine wwstsc 
!
!
!----------------------------------------------------------------------
      subroutine wwstss (aa,la,ma,n,jm,ja,jd)
!----------------------------------------------------------------------
      use, intrinsic :: iso_c_binding 
      implicit none 
!
      real(C_DOUBLE) aa(la,ma)
      integer(C_INT) ja(la,ma),jd(n)
!
      integer(C_INT) i,j,j0,j1,la,ma,n,jm,jj
      real(C_DOUBLE) a0
!
      do 100 i=1,n
         do 110 j0=1,jm
            j  = ja(i,j0)
            a0 = aa(i,j0)
            if ( j.eq.i .and. abs(a0).lt.1.d-72 ) then
               ja(i,j0) = n + 1
               goto 110
            endif
!
            do 120 j1=j0-1,1,-1
               if ( ja(i,j1).gt.j ) then
                  ja(i,j1+1) = ja(i,j1)
                  aa(i,j1+1) = aa(i,j1)
               else
                  ja(i,j1+1) = j
                  aa(i,j1+1) = a0
                  goto 110
               endif
  120       continue
            ja(i,1) = j
            aa(i,1) = a0
  110    continue
         do 130 jj=1,jm
            if ( ja(i,jj).eq.i ) then
               jd(i) = jj
            elseif ( ja(i,jj).eq.n+1 ) then
               ja(i,jj) = i
            endif
  130    continue
  100 continue
!
      return
      end subroutine wwstss
!
!
!----------------------------------------------------------------------
      subroutine wwstsi (aa,la,ma,n,jm,ja,jd,iprc,wp,ierr)
!----------------------------------------------------------------------
      use, intrinsic :: iso_c_binding 
      implicit none
!
      real(C_DOUBLE) aa(la,ma),wp(n)
      integer(C_INT) ja(la,ma),jd(n),la,ma,n,jm,iprc,ierr, &
                     i,j,j0,j1
!
      if ( ierr.eq.1000 .or. ierr.eq.1010 ) then
         do 100 i=1,n
            wp(i) = 1.d0
  100    continue
         return
      endif
!
      if ( iprc.eq.0 ) then
         do 110 i=1,n
            wp(i) = 1.d0/aa(i,jd(i))
  110    continue
         return
      endif
!
      do 120 i=1,n
         wp(i) = aa(i,jd(i))
         do 130 j0=1,jd(i)-1
            j = ja(i,j0)
            if ( j.ge.i ) goto 130
            do 140 j1=jd(j)+1,jm
               if ( ja(j,j1).eq.i ) then
                  wp(i) = wp(i) - aa(i,j0)*wp(j)*aa(j,j1)
                  goto 130
               endif
  140       continue
  130    continue
         if ( abs(wp(i)).lt.0.1d0*abs(aa(i,jd(i))) ) then
            wp(i) = 0.1d0*aa(i,jd(i))
            ierr = 1210
         endif
         wp(i) = 1.d0/wp(i)
  120 continue
!
      return
      end subroutine wwstsi 
!
!
!----------------------------------------------------------------------
      subroutine wwstsj (aa,la,ma,n,jm,ja,jd,wp,iprc,v,w)
!----------------------------------------------------------------------
      use, intrinsic :: iso_c_binding 
      implicit none
!
      include   'paramAPGa.h'
      real(C_DOUBLE) aa(la,ma),wp(n),v(n),w(n),wa(mx)
      integer(C_INT) ja(la,ma),jd(n),la,ma,n,jm,iprc
!
      integer(C_INT) i,jj,ik,k
!     integer(C_INT) io_pe
!     common/sub_proc/ io_pe
!
      if ( iprc.eq.0 ) then
         do 100 i=1,n
            w(i) = wp(i)*v(i)
  100    continue
         return
      endif
!
!***
!     if(jm.eq.9) go to 1000
!***
!
      do 110 i=1,n
         w(i) = v(i)
         do 120 jj=1,jd(i)-1
            w(i) = w(i) - aa(i,jj)*w(ja(i,jj))
  120    continue
         w(i) = w(i)*wp(i)
  110 continue
!
      do 130 i=n-1,1,-1
!
         do 140 jj=jd(i)+1,jm
            w(i) = w(i) - wp(i)*aa(i,jj)*w(ja(i,jj))
  140    continue
  130 continue
!
      return
!
!------------------------------------------------------------------
!   A shortcut version of ipc=1 case  (nob2=9).    date: 3/7/1993.
!------------------------------------------------------------------
!*    In a periodic case, k=1 and mz columns are special.
!*      because sorting has been made in /wwstss/.
!
 1000 continue
!     if(io_pe.eq.1) then
!     OPEN (unit=11,file=praefixc//'.06'//suffix2,             &
!           status='unknown',position='append',form='formatted')
!
!     write(11,*) 'this portion has not been corrected.......'
!     close(11)
!     end if
!
!********************************************************
!*  1. Going up the matrix /a/.                         *
!********************************************************
!
!*  for k=1.
!
      w(1)= v(1)*wp(1)
!
      do 1100 i= 1,mx
      ik= i
      w(ik)= (v(ik) -aa(ik,1)*w(ja(ik,1)))*wp(ik)
 1100 continue
!
!*  for general k.
!
      do 1200 k= 2,mz-1
      ik= 1 +mx*(k-1)
      w(ik)= v(ik)*wp(ik)
!
!
      do 1250 i= 1,mx
      ik= i +mx*(k-1)
      wa(i)=  aa(ik,1)*w(i-1+mx*(k-2))  &
             +aa(ik,2)*w(i  +mx*(k-2))  &
             +aa(ik,3)*w(i+1+mx*(k-2))
 1250 continue
!
      do 1300 i= 1,mx
      ik= i +mx*(k-1)
      w(ik)= (v(ik) -wa(i) -aa(ik,4)*w(i-1+mx*(k-1)) )*wp(ik)
 1300 continue
!
!
      ik= mx +mx*(k-1)
      w(ik)= (v(ik) -aa(ik,1)*w(mx-1+mx*(k-1)) )*wp(ik)
 1200 continue
!
!
!*  for k= mz.
!
      ik= 1 +mx*(mz-1)
      w(ik)= v(ik)*wp(ik)
!
      do 1400 i= 1,mx
      ik= i +mx*(mz-1)
      w(ik)= (v(ik) -aa(ik,1)*w(ja(ik,1))  &
                    -aa(ik,2)*w(ja(ik,2))  &
                    -aa(ik,3)*w(ja(ik,3))  &
                    -aa(ik,4)*w(ja(ik,4))  &
                    -aa(ik,5)*w(ja(ik,5))  &
                    -aa(ik,6)*w(ja(ik,6))  &
                    -aa(ik,7)*w(ja(ik,7)) )*wp(ik)
 1400 continue
!
      ik= mx +mx*(mz-1)
      w(ik)= ( v(ik) -aa(ik,1)*w(ja(ik,1)) )*wp(ik)
!
!********************************************************
!*  2. Going down the matrix /a/.                       *
!********************************************************
!
!*  for k= mz.
!
      do 1500 i= mx,1,-1
      ik= i +mx*(mz-1)
      w(ik)= w(ik) -wp(ik)*aa(ik,9)*w(ja(ik,9))
 1500 continue
!
      ik= 1 +mx*(mz-1)
      w(ik)= w(ik) -wp(ik)*aa(ik,2)*w(ja(ik,2))
!
!*  for general k.
!
      do 1600 k= mz-1,2,-1
!
!
      do 1650 i= mx,1,-1
      ik= i +mx*(k-1)
      wa(i)=  aa(ik,7)*w(i-1+mx*k)  &
             +aa(ik,8)*w(i  +mx*k)  &
             +aa(ik,9)*w(i+1+mx*k)
 1650 continue
!
      do 1700 i=   mx,1,-1
      ik= i +mx*(k-1)
      w(ik)= w(ik) -wp(ik)*( aa(ik,6)*w(i+1+mx*(k-1)) +wa(i) )
 1700 continue
!
!
      ik= 1 +mx*(k-1)
      w(ik)= w(ik) -wp(ik)*aa(ik,2)*w(2+mx*(k-1))
 1600 continue
!
!
!*  for k= 1.
!
      do 1800 i= mx,1,-1
      ik= i
      w(ik)= w(ik) -wp(ik)*( aa(ik,3)*w(ja(ik,3))  &
                            +aa(ik,4)*w(ja(ik,4))  &
                            +aa(ik,5)*w(ja(ik,5))  &
                            +aa(ik,6)*w(ja(ik,6))  &
                            +aa(ik,7)*w(ja(ik,7))  &
                            +aa(ik,8)*w(ja(ik,8))  &
                            +aa(ik,9)*w(ja(ik,9)))
 1800 continue
!
      ik= 1
      w(ik)= w(ik) -wp(ik)*aa(ik,2)*w(ja(ik,2))
!
!
      return
      end subroutine wwstsj 
!
!
!----------------------------------------------------------------------
      subroutine wwstsm (aa,la,ma,n,jm,ja,v,w)
!----------------------------------------------------------------------
      use, intrinsic :: iso_c_binding 
      implicit none 
!
      real(C_DOUBLE) aa(la,ma),v(n),w(n)
      integer(C_INT) ja(la,ma),la,ma,n,jm,i,j,jj
!
      do 100 i=1,n
        w(i) = 0.d0
  100 continue
!
      do 110 jj=1,jm
         do 120 i=1,n
            j = ja(i,jj)
            w(i) = w(i) + aa(i,jj)*v(j)
  120    continue
  110 continue
!
      return
      end subroutine wwstsm 
!
!
!******************************************************* 7/15/93 *******
!*                                                                     *
!*    Conjugate residual method  (for non-symmetric matrix).           *
!*                                                                     *
!**************************************** parallel version: 11/4/94 ****
!-----------------------------------------------------------------------
      subroutine cresmd (xx,ss,ierr,ipar)
!-----------------------------------------------------------------------
!
!       solve:  mat(aa)*x = b.
!
!     call cresmd (ierr)
!     common/cresm1/ aa(mxyz,nob3)
!     common/cresma/ na(mxyz,nob3)
!     common/cresm2/ x(mxyz)
!     common/cresmc/ b(mxyz)
!     common/cresmb/ nobx
!     common/cresm3/ ipr(10),rpr(10)
!
!    arguments :
!
!    aa(mxyz,nob3) (real*4, 8) --------- coefficient matrix
!    na(mxyz,nob3) (integer*4) --------- column no. table
!            na(i,jj) gives a location of the base vector for a(i,jj).
!    note: nob2= 9; but, any nobx .le. nob2 is acceptable.
!
!    b(mxyz)     (real*4, 8) -------------- source vector
!    x(mxyz)     (real*4, 8) -------------- solution (unknown)
!
!    ipr(10) (integer*4) ----------- integer parameter
!                   ipr(1) (in) .... maximum number of iteration
!                   ipr(2) (out) ... final iteration number
!    rpr(10) (real*8) -------------- real parameter
!                   rpr(1) (in) .... tolerance of convergence
!                   rpr(2) (out) ... residual norm
!                   rpr(3) (out) ... cpu time
!    ierr (integer*4) -------------- return code
!                     = 0    ......... normal termination
!                     = 1000 ......... source vector = 0.
!                     = 2000 ......... no convergence within itrm
!                     = 3000 ......... blow up (x > 1.d+50)
!
!                                   copyright : m.tanaka  7/15/1993.
!----------------------------------------------------------------------
      use, intrinsic :: iso_c_binding 
      implicit none
!
      include    'mpif.h'
      include    'paramAPGa.h'
!
!     integer(C_INT) nob3
!     parameter  (nob3=13)
!
      real(C_DOUBLE),dimension(0:mxyz-1) :: xx,ss
      integer(C_INT) ierr,ipar
!
      real(C_DOUBLE) sq0,sq1,rpr,eps
      integer(C_INT) ipr
      common/cresm3/ ipr(10),rpr(10)
      integer(C_INT) io_pe,i3,i4,cnt_rec2,disp_rec2,ncrs,itrm
      common/sub_proc/ io_pe
      common/dat_typ2/ i3(30),i4(30),cnt_rec2(30),disp_rec2(30)
!
!*-------------------------------------------------------------------
!*  Restart procedure is taken in /cressl/ iteration if convergence 
!*  slows down (ierr= 2000).
!*                             (also in /wwstbk/)  ** 04/22/1998 ** 
!*-------------------------------------------------------------------
!
      ncrs= 1
!
 2000 continue
      itrm = ipr(1)
      eps  = rpr(1)
!
      call cresin (xx,ss,sq0,sq1,ierr,ipar)
      if(ierr.ne.0) go to 1000
!
      call cressl (xx,ss,sq0,sq1,eps,ncrs,itrm,ierr,ipar)
!
!     if(ierr.eq.2000) then
!       if(ncrs.le.2) then
!          syme = -1.d0
!          syme2= -1.d0
!          call filt1L (xx,1,1,1,syme,syme2)
!          ierr= 1700 +ncrs
!
!         if(io_pe.eq.1) then
!         OPEN (unit=11,file=praefixc//'.06'//suffix2,             &
!               status='unknown',position='append',form='formatted')
!
!          write(11,*) ' reset and restart (cressl), ncrs=',ncrs
!          close(11)
!          end if
!
!          go to 2000
!       end if
!     end if
!
 1000 continue
!
      ipr(2) = itrm
      rpr(2) = eps
      rpr(3) = 0  ! c1 - c0
!
      return
      end subroutine cresmd 
!
!
!----------------------------------------------------------------------
      subroutine cresin (x,b,sq0,sq1,ierr,ipar)
!----------------------------------------------------------------------
!**********************************************************************
!*    Initialization.                                                 *
!**********************************************************************
      use, intrinsic :: iso_c_binding 
      implicit none
!
      include    'mpif.h' 
      include    'paramAPGa.h'
!
!     integer     nob3
!     parameter  (nob3=13)
!*--------------------------------------------------------------------
      real(C_DOUBLE),dimension(0:mxyz-1) :: x,b
      real(C_DOUBLE) sq0,sq1,swq,swq2
!
      real(C_DOUBLE),dimension(0:mxyz-1) :: p,q,r,s,p1,p2,q1,q2
      common/cresm4/  p,q,r,s,p1,p2,q1,q2
!***
      integer(C_INT)  ipar,mpierror,ierr,i
      integer(C_INT)  io_pe,i3,i4,cnt_rec2,disp_rec2, &
                      i2,cnt_recv,disp_recv
      common/sub_proc/ io_pe
      common/dat_typ2/ i3(30),i4(30),cnt_rec2(30),disp_rec2(30)
!*--------------------------------------------------------------------
!
!*  s = a*x
      call avmult (x,s,ipar)
!
      do i= i3(ipar),i4(ipar)
      p1(i)= 0
      q1(i)= 0
!
      r(i)= b(i) -s(i)
      p(i)= r(i)
      end do
!
!* q = a*r
      call avmult (r,q,ipar)
!
!
      swq= 0.d0
!
      do i= i3(ipar),i4(ipar)
      swq= swq + q(i)**2   ! mpi_sum
      end do
!
      if(num_proc.ne.1) then
        call mpi_allreduce (swq,swq2,1,mpi_real8,mpi_sum, &
                            mpi_comm_world,mpierror)
        swq= swq2
      end if
!
      if(swq.gt.1.d-30) then
         sq0= swq
         ierr= 0
      else
         sq0= 0.
         ierr= 1000
      end if
!
      sq1= 1.d0
!
      return
      end subroutine cresin 
!
!
!----------------------------------------------------------------------
      subroutine cressl (x,b,sq0,sq1,eps,ncrs,itrm,ierrr,ipar)
!----------------------------------------------------------------------
!**********************************************************************
!*    iteration step of conjugate residual method.                    *
!**********************************************************************
      use, intrinsic :: iso_c_binding 
      implicit none
!
      include    'mpif.h'
      include    'paramAPGa.h'
!
!     integer     nob3
!     parameter  (nob3=13)
!*--------------------------------------------------------------------
      real(C_DOUBLE),dimension(0:mxyz-1) :: x,b
      real(C_DOUBLE)  eps,srhs
      integer(C_INT)  ipar,ncrs
!
      real(C_DOUBLE),dimension(0:mxyz-1) :: p,q,r,s,p1,p2,q1,q2
      common/cresm4/  p,q,r,s,p1,p2,q1,q2
!
      integer(C_INT) io_pe,i3,i4,cnt_rec2,disp_rec2, &
                     i0,i2,cnt_recv,disp_recv
      common/sub_proc/ io_pe
      common/dat_typ2/ i3(30),i4(30),cnt_rec2(30),disp_rec2(30)
!*--------------------------------------------------------------------
      real(C_DOUBLE) sq0,sq1,sq2,alpha,beta1,beta2,wb1,wb2, &
                     wrq,wrq2,wr1,wxs,buf1(2),buf2(2),      &
                     sw1,sw2,wq,wq2,swr,ssw
!
      integer(C_INT) i,itr,itrm,mpierror,ierr,ierrr
      real(C_DOUBLE) sres,avex,rsdl,conv,conv0
      common/crconv/ sres,avex,rsdl,conv,ierr
!
!*
      wb1= 0.d0
!
      do i= i3(ipar),i4(ipar)
      wb1= wb1 +b(i)**2
      end do
!
      if(num_proc.ne.1) then
        call mpi_allreduce (wb1,wb2,1,mpi_real8,mpi_sum, &
                            mpi_comm_world,mpierror)
        wb1= wb2
      end if
!
      srhs=  sqrt(wb1/float(mxyz))
!
!                        ** Iteration starts **
      itr= 0
 1000 itr= itr +1
!
      sq2= sq1
      sq1= sq0
!
      wrq= 0.d0
!
      do i= i3(ipar),i4(ipar)
      p2(i)= p1(i)
      p1(i)=  p(i)
      q2(i)= q1(i)
      q1(i)=  q(i)
!
      wrq= wrq + r(i)*q1(i)  ! mpi_sum
      end do
!
      if(num_proc.ne.1) then
        call mpi_allreduce (wrq,wrq2,1,mpi_real8,mpi_sum, &
                            mpi_comm_world,mpierror)
        wrq= wrq2
      end if
!
      alpha= wrq/sq1
!                +++
      wr1= 0.d0
      wxs= 0.d0
!
      do i= i3(ipar),i4(ipar)
      x(i)= x(i) +alpha*p1(i)
      r(i)= r(i) -alpha*q1(i)
!
      wr1= wr1 +r(i)**2   ! mpi_sum
      wxs= wxs +x(i)**2
      end do
!
      if(num_proc.ne.1) then
        buf1(1)= wr1
        buf1(2)= wxs
        call mpi_allreduce (buf1,buf2,2,mpi_real8,mpi_sum, &
                            mpi_comm_world,mpierror)
        wr1= buf2(1)
        wxs= buf2(2)
      end if
!
!**********************************************************************
!*    Convergence check.  (compare residual with source)              *
!**********************************************************************
!----------------------------------------------------------------------
!
      sres= sqrt(wr1/float(mxyz))
      avex= sqrt(wxs/float(mxyz))
!
      conv0= conv
      conv= sres/(avex +1.e-15) ! (srhs +1.e-15)
!
!     if(io_pe.eq.1) then
!     if(mod(itr,100).eq.0) then
!       if(io_pe.eq.1) then
!       OPEN (unit=11,file=praefixc//'.06'//suffix2,             &
!             status='unknown',position='append',form='formatted')
!
!       write(11,900) itr,sres,conv,avex
! 900   format('#(crm) itr=',i3,' res=',1pe12.5,', res/<x>=', &
!            e12.5,', <x>=',e12.5)
!       close(11)
!       end if
!     end if
!     end if
!                                 ** Successfully converged **
      if(conv.lt.eps) then
         ierr= 0
         go to 7000
      end if
!
!                                 ** Source too small **
      if(srhs.lt.1.d-30) then
         ierr= 1000
         go to 7000
      end if
!
!                                 ** blow up **
      if(avex.gt.1.d+50) then
        ierr= 3000
        go to 7000
      end if
!----------------------------------------------------------------------
!                                       s = a*r
      call avmult (r,s,ipar)
!
      sw1= 0.d0
      sw2= 0.d0
!*
      do i= i3(ipar),i4(ipar)
      sw1= sw1 +s(i)*q1(i)  ! mpi_sum
      sw2= sw2 +s(i)*q2(i)
      end do
!
      if(num_proc.ne.1) then
        buf1(1)= sw1
        buf1(2)= sw2
        call mpi_allreduce (buf1,buf2,2,mpi_real8,mpi_sum, &
                            mpi_comm_world,mpierror)
        sw1= buf2(1)
        sw2= buf2(2)
      end if
!
      beta1= -sw1/sq1
      beta2= -sw2/sq2
!                 +++
!
      wq= 0.d0
!
      do i= i3(ipar),i4(ipar)
      p(i)= r(i) +beta1*p1(i) +beta2*p2(i)
      q(i)= s(i) +beta1*q1(i) +beta2*q2(i)
!
      wq= wq + q(i)**2   ! mpi_sum
      end do
!
      if(num_proc.ne.1) then
        call mpi_allreduce (wq,wq2,1,mpi_real8,mpi_sum, &
                            mpi_comm_world,mpierror)
        wq= wq2
      end if
!
      sq0= wq
!
!-------------------------------------------------------------------
!*  Round-off error is dominating, thus restart iterations.
!--------------------------------------------- added: 04/22/1998 ---
!
!     rcc= abs(conv -conv0)/(conv0 +1.e-10)
!     if(itr.gt.30.and.rcc.lt.1.e-4) go to 6000
      if(itr.lt.itrm) go to 1000
!
 6000 ierr= 2000
      ncrs= ncrs +1
!
 7000 ierrr= ierr
      eps=  conv
      itrm= itr
!
      return
      end subroutine cressl 
!
!
!----------------------------------------------------------------------
      subroutine avmult (v,w,ipar)
!----------------------------------------------------------------------
!**********************************************************************
!   Calculate: w(i) = sum(j)| a(i,j)*v(na(i,j)) |
!**********************************************************************
      use, intrinsic :: iso_c_binding 
      implicit  none
!
      include   'mpif.h'
      include   'paramAPGa.h'
!
      integer(C_INT) nobx 
      common/cresmb/ nobx
!     parameter  (nob3=13)
!***
      real(C_DOUBLE),dimension(0:mxyz-1) :: v,w,v8,w8,vv
      integer(C_INT) ipar,i,j,k,jj,i00,ierror
!
      real(C_DOUBLE) aa
      integer(C_INT) na
      common/cresm1/  aa(0:mxyz-1,nob3)
      common/cresma/  na(0:mxyz-1,nob3)
!
      integer(C_INT) io_pe,i3,i4,cnt_rec2,disp_rec2
      integer(C_INT) cnt_rec3,disp_rec3,disp_rec4,cnt_send
      common/sub_proc/ io_pe
      common/dat_typ2/ i3(30),i4(30),cnt_rec2(30),disp_rec2(30)
      common/dat_typ3/ cnt_rec3(30),disp_rec3(30),disp_rec4(30)
!
!* Parallelization is faster by num_proc times !
!
      do i= i3(ipar),i4(ipar)
      w8(i) = 0.d0
      v8(i) = v(i)
      end do
!
      if(num_proc.ne.1) then
!
        do k= 1,num_proc
        disp_rec3(k)= i3(k)
        disp_rec4(k)= i4(k)-mxy+1
!
        cnt_rec3(k)= mxy
        end do
!
        i00= i3(ipar)
        do i= i3(ipar),i3(ipar)+mxy-1
        vv(i-i00)= v(i)
        end do
!
        cnt_send= mxy
        call mpi_allgatherv (vv,cnt_send,          mpi_real8, &
                             v8,cnt_rec3,disp_rec3,mpi_real8, &
                             mpi_comm_world,ierror)
!
        i00= i4(ipar)-mxy+1
        do i= i4(ipar)-mxy+1,i4(ipar)
        vv(i-i00)= v(i)
        end do
!
        cnt_send= mxy
        call mpi_allgatherv (vv,cnt_send,          mpi_real8, &
                             v8,cnt_rec3,disp_rec4,mpi_real8, &
                             mpi_comm_world,ierror)
      end if
!
!* pgi - loop (nobx) should be an inner loop (faster)
!  vpp - nobx loop must be outside
!
      do i= i3(ipar),i4(ipar)
      do j= 1,nobx
      w8(i)= w8(i) + aa(i,j)*v8(na(i,j))
      end do
      end do
! 
!* --------------------------------------------------------------
!
      do i= i3(ipar),i4(ipar)
      w(i)= w8(i)
      end do
!
      return
      end subroutine avmult 
!
!
!------------------------------------------------------
      subroutine rehist
!------------------------------------------------------
      use, intrinsic :: iso_c_binding 
      implicit none
!
      include      'paramAPGa.h'
!
      real(C_float),dimension(ntmax) :: &
                    ekin,ppot,ekn2,etot,z_pe,vzpe,vzco,vzct, &
                    vdtm,vpor,time,ecr,elj,espr
      common/ehist/ ekin,ppot,ekn2,etot,z_pe,vzpe,vzco,vzct, &
                    vdtm,vpor,time,ecr,elj,espr
!
      integer(C_INT) io_pe
      common/sub_proc/ io_pe
      integer(C_INT) it,is
      common/parm1/ it,is
!
      real(C_float)  t,t00,xp_leng
      common/headr2/ t,t00,xp_leng
      integer(C_INT) iwa,iwb,iwc,iwd,j,k,k1
      common/imemo/  iwa,iwb,iwc,iwd
!
      real(C_DOUBLE) xmax,ymax,zmax,xmin,ymin,zmin
      real(C_float)  phi,tht,dtwr1,dtwr2,dtwr3
      common/parm3/ xmax,ymax,zmax,xmin,ymin,zmin
      common/parm4/ phi,tht,dtwr1,dtwr2,dtwr3
!
!
      do k= 1,is
      if(mod(k,2).eq.0) then
!
      k1= k/2
      ekin(k1)= (ekin(k-1) +ekin(k))/2
      ppot(k1)= (ppot(k-1) +ppot(k))/2
      ekn2(k1)= (ekn2(k-1) +ekn2(k))/2
      etot(k1)= (etot(k-1) +etot(k))/2
      z_pe(k1)= (z_pe(k-1) +z_pe(k))/2
      vzpe(k1)= (vzpe(k-1) +vzpe(k))/2
      vzco(k1)= (vzco(k-1) +vzco(k))/2
      vzct(k1)= (vzct(k-1) +vzct(k))/2
!
      vdtm(k1)= (vdtm(k-1) +vdtm(k))/2
       vpor(k1)= ( vpor(k-1) + vpor(k))/2
      time(k1)= (time(k-1) +time(k))/2
       ecr(k1)= ( ecr(k-1) + ecr(k))/2
       elj(k1)= ( elj(k-1) + elj(k))/2
      espr(k1)= (espr(k-1) +espr(k))/2
!
      end if
      end do
!
!* Reset iwa for proper timing
!
      is  = is/2
      dtwr1= 2*dtwr1
      iwa = t/dtwr1 
!
      if(io_pe.eq.1) then
        OPEN (unit=11,file=praefixc//'.06'//suffix2,             &
              status='unknown',position='append',form='formatted')
!
        write(11,*) 'rehist: new is, dtwr1=',is,dtwr1
        close(11)
      end if
!
      return
      end subroutine rehist
!
!
!------------------------------------------------------
      subroutine lplots
!------------------------------------------------------
      use, intrinsic :: iso_c_binding 
      implicit none
!
      include      'paramAPGa.h'
!
      real(C_float),dimension(ntmax) :: &
                    ekin,ppot,ekn2,etot,z_pe,vzpe,vzco,vzct, &
                    vdtm,vpor,time,ecr,elj,espr
      common/ehist/ ekin,ppot,ekn2,etot,z_pe,vzpe,vzco,vzct, &
                    vdtm,vpor,time,ecr,elj,espr
!
      real(C_float)  t,t00,xp_leng
      integer(C_INT) it,is,ILN,ILG,nsg,nseg
      character*8    label,cdate*10
!
      common/headr1/ label,cdate
      common/headr2/ t,t00,xp_leng
      common/parm1/ it,is
      common/psegm/ nsg(30),nseg
!
      real(C_float) emax1,emin1,emax2,emin2,pmax,pmin,    &
                    etmax,etmin,elmax,elmin,e3max,e3min,  &
                    vdtm1,vdtm2,emax,cjz1,cjz2,qz11,qz12, &
                    qz21,qz22,vzpe1,vzpe2,vzct1,vzct2,    &
                    vzco1,vzco2,qtop1a,qtop1b,qtop2a,qtop2b,   &
                    qpop1a,qpop1b,qpop2a,qpop2b,qbot1a,qbot1b, &
                    qbot2a, qbot2b,qqmax
!
      ekin(1)= ekin(3)
      ekin(2)= ekin(3)
      ppot(1)= ppot(3)
      ppot(2)= ppot(3)
      etot(1)= etot(3)
      etot(2)= etot(3)
!
      call lplmax (ekin,emax1,emin1,is)
      call lplmax (ekn2,emax2,emin2,is)
      call lplmax (ppot,pmax,pmin,is)
      call lplmax (etot,etmax,etmin,is)
      call lplmax ( elj,elmax,elmin,is)
      call lplmax (espr,e3max,e3min,is)
      call lplmax (vdtm,vdtm1,vdtm2,is)
      emax= max(emax1,emax2)
!
      ILN= 1
      ILG= 2
!
      call lplot1 (2,4,is,time,ekin,emax1,0.0,iln,'kin ions',8, &
                 '        ',8,'        ',8)
      call lplot1 (2,5,is,time,ekn2,emax2,0.0,iln,'kin solv',8, &
                 '        ',8,'        ',8)
      call lplot1 (2,6,is,time,ppot,pmax,pmin,iln,'es ener ',8, &
                 '  time  ',8,'        ',8)
!
!
      call lplot1 (3,4,is,time, elj,elmax,elmin,iln,'e_lj    ',8, &
                 '        ',8,'        ',8)
      call lplot1 (3,5,is,time,etot,etmax,etmin,iln,'e.total ',8, &
                 '        ',8,'        ',8)
      call lplot1 (3,6,is,time,vdtm,vdtm1,vdtm2,iln,'v*dt    ',8, &
                 '  time  ',8,'        ',8)
!     call lplot1 (3,6,is,time,espr,e3max,e3min,iln,'e_elas  ',8, &
!                '  time  ',8,'        ',8)
!------------------------
      call chart
!------------------------
!
      return
      end subroutine lplots
!
!
!------------------------------------------------------
      subroutine lplmax (f,fmax,fmin,is)
!------------------------------------------------------
      use, intrinsic :: iso_c_binding 
      implicit none
!
      integer(C_INT) is,i
      real(C_float)  f(is),fmax,fmin
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
!-------------------------------------------------------------
      subroutine ppl3da (xg,yg,zg,ch,ag,rod_leng,Rmac,  &
                         np,nq,nCLp,npqr,first_ppl)
!-------------------------------------------------------------
      use, intrinsic :: iso_c_binding 
      implicit none
!
      include   'paramAPGa.h'
!
      real(C_DOUBLE),dimension(npqr0) :: xg,yg,zg,ch,ag
      real(C_float),dimension(npqr0) ::  x,y,z
!
      real(C_DOUBLE) xmax,ymax,zmax,xmin,ymin,zmin,         &
                     pi,dt,axi,Gamma,rbmax,vth,tmax,        &
                     Bjerrum,qfrac,Rpore,Hpore,Zci,Zcp,Zcn, &
                     rod_leng,Rmac
      real(C_float)  phi,tht,dtwr1,dtwr2,dtwr3
      integer(C_INT) nsg,nseg,ifqq,nCLp,np,npqr,nq
!
      common/parm2/  pi,dt,axi,Gamma,rbmax,vth,tmax
      common/parm3/  xmax,ymax,zmax,xmin,ymin,zmin
      common/parm4/  phi,tht,dtwr1,dtwr2,dtwr3
      common /elsta/ Bjerrum
      common/cntion/ qfrac,Rpore,Hpore,Zci,Zcp,Zcn,ifqq
      common/psegm/  nsg(30),nseg
!
      integer(C_INT) n_rodp,jj,np00,nnb,ist1,ist2
      common/pbase/  np00,nnb,ist1(100),ist2(100)
!
      real(C_DOUBLE) Vtop0,Vbot0,Vtop,Vbot,t_recyc,rpr
      real(C_DOUBLE) diel2,dielpr,aa
      integer(C_INT) fixedbc_x,fixedbc_y,itermax,filtx,filty,filtz,ipr
      common/espot/  Vtop0,Vbot0,Vtop,Vbot,fixedbc_x,fixedbc_y, &
                     itermax,filtx,filty,filtz
      common/dielec/ diel2,dielpr,aa
!
      character*8    label,cdate*10,cax*1
      real(C_float)  t,t00,xp_leng
      common/headr1/ label,cdate
      common/headr2/ t,t00,xp_leng
!
      integer(C_INT) io_pe
      common/sub_proc/ io_pe
!
      integer(C_INT) i,k,kpl,ns,ia,ib,l
      real(C_float)  hh,Rpore4,Hpore4,Zcp4,Zcn4,Bjerrum4,xleng4,zleng4, &
                     Rmac4,Vtop4,Vbot4,Gamma4,rod_leng4,diel24,         &
                     fsize,hl,vd,pha,tha,cph,sph,cth,sth,xp,yp,zp,      &
                     rmax1,ps,x1,y1,z1,xx,yy,dd,xpp,ypp,zpp
      real(C_DOUBLE) xleng,yleng,zleng
      common/parm8/  xleng,yleng,zleng
!
      real(C_DOUBLE) awat,wwat
      common/parmd/  awat,wwat
      logical        first_ppl
!
      fsize= 8.
      hl=  11.
      vd=  6.
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
      xp= xleng*cph -yleng*sph
      yp= xleng*sph +yleng*cph
      zp= zleng
!
      ypp=  yp
      zpp= -xp*sth +zp*cth
!
      rmax1= sqrt(ypp**2 +zpp**2)
      ps= fsize/rmax1
!
!**********************
!*  Draw axes.        *
!**********************
!
      do i= 1,npqr
      x(i) = xg(i)
      y(i) = yg(i)
      z(i) = zg(i)
      end do
!
!
      do 1000 ns= 1,2
!*
      hh= 0.7
      call symbol (0.5,17.0,hh,label,0.,8)
      call symbol (4.5,17.0,hh,praefixc,0.,31)
!
      call symbol ( 0.5,1.0,hh,cdate, 0.,10)
      call symbol (15.9,1.0,hh,'t=', 0.,2)
      call number (999.,1.0,hh,t,0.,7)
!
      Rpore4= Rpore
      Hpore4= Hpore
      Zcp4  = Zcp
      Zcn4  = Zcn
      xleng4 = xleng
      zleng4 = zleng
!
      call symbol ( 0.5,16.0,hh,'Rpore=', 0.,6)
      call number ( 2.5,16.0,hh,Rpore4,0.,5)
      call symbol ( 5.5,16.0,hh,'Zcp=', 0.,4)
      call number ( 7.5,16.0,hh,Zcp4,0.,5)
      call symbol (10.5,16.0,hh,'np=', 0.,3)
      call number (13.0,16.0,hh,float(np),0.,5)
      call symbol (16.0,16.0,hh,'xleng=',0.,6)
      call number (18.0,16.0,hh,xleng4,0.,5)
!
      Rmac4 = Rmac
      call symbol ( 0.5,15.3,hh,'Hpore=', 0.,6)
      call number ( 2.5,15.3,hh,Hpore4,0.,5)
      call symbol ( 5.5,15.3,hh,'Zcn=', 0.,4)
      call number ( 7.5,15.3,hh,Zcn4,0.,5)
      call symbol (10.5,15.3,hh,'nq=', 0.,3)
      call number (13.0,15.3,hh,float(nq),0.,5)
      call symbol (16.0,15.3,hh,'zleng=', 0.,6)
      call number (18.0,15.3,hh,zleng4,0.,5)
!
      Vtop4= Vtop
      Vbot4= Vbot
      Gamma4= Gamma
      diel24= diel2
      rod_leng4= rod_leng
!
      call symbol ( 0.5,14.6,hh,'V_tb=', 0.,5)
      call number ( 2.5,14.6,hh,Vtop4-Vbot4,0.,5)
      call symbol ( 5.5,14.6,hh,'Gamma=', 0.,6)
      call number ( 7.5,14.6,hh,Gamma4,0.,5)
      call symbol (10.5,14.6,hh,'rod_l=', 0.,6)
      call number (13.0,14.6,hh,rod_leng4,0.,5)
      call symbol (16.0,14.6,hh,'diel2=', 0.,6)
      call number (18.0,14.6,hh,diel24,0.,5)
!*
!
      do i= 1,3
      if(i.eq.1) then
         x1= rmax1
         y1= 0.
         z1= 0.
         cax='x'
      else if(i.eq.2) then
         x1= 0.
         y1= rmax1
         z1= 0.
         cax='y'
      else if(i.eq.3) then
         x1= 0.
         y1= 0.
         z1= rmax1
         cax='z'
      end if
!
      xp= x1*cph -y1*sph
      yp= x1*sph +y1*cph
      zp= z1
!
      ypp= yp
      zpp= -xp*sth +zp*cth
!
      xx= ps*ypp  +hl
      yy= ps*zpp  +vd
      call plot (hl,vd,3)
      call plot (xx,yy,2)
!
      call symbol (xx-0.7,yy-0.5,hh,cax,0.,1)
      end do
!
!  Monomers
!-------------------
!     dd= 0.30
      dd= 1.0*ps
!-------------------
!
      do i= 1,nCLp
      xp= x(i)*cph -y(i)*sph
      yp= x(i)*sph +y(i)*cph
      zp= z(i)
!
      xpp=  xp*cth +zp*sth
      ypp= yp
      zpp= -xp*sth +zp*cth
!
      xx= ps*ypp +hl
      yy= ps*zpp +vd
!
      dd= ps
      if(i.le.np)   dd= 2.0*ps  ! 1.5*ps
      if(i.gt.nCLp) dd= 0.5*ps
!
      if(i.le.np) then
        if(ch(i).lt.0.d0) then
          call newcolor(3,1.,0.,0.)
          call circle (xx-0.15,yy-0.15,dd,2)  ! 2: filled circle, red
!
        else if(ch(i).ge.0.d0) then
          call newcolor(3,0.,0.,1.)
          call circle (xx-0.15,yy-0.15,dd,2)  ! 2: filled circle, blue
        end if
      else
!
        call newcolor(0,1.,0.,0.) 
        call circle (xx-0.15,yy-0.15,dd,1)  ! 1: open circle
      end if
      end do
!
      call newcolor(0,1.,0.,0.)  ! reset
!
!-----------
!*  Chain
!-----------
!     integer(C_INT) n_rodp,jj,np00,nnb,ist1,ist2
!     common/pbase/  np00,nnb,ist1(100),ist2(100)
!
      call newcolor(3,0.,1.,0.)
!
      do k= 1,nseg
      ia= nsg(k) +1 
      ib= nsg(k) +24  ! 12+12 monomers
!
      l= 0
      do i= ia,ib
      l= l +1
!
      xp= x(i)*cph -y(i)*sph
      yp= x(i)*sph +y(i)*cph
      zp= z(i)
!
      xpp=  xp*cth +zp*sth
      ypp= yp
      zpp= -xp*sth +zp*cth
!
      xx= ps*ypp +hl
      yy= ps*zpp +vd
!
      if(l.eq.1) then
        call plot (xx,yy,3)  ! line starts
      else
        call plot (xx,yy,2)  ! line continues
      end if
      end do
!
      call plot (xx,yy,3)    ! close the line 
      end do
!
      do k= 1,nnb
      ia= ist1(k) 
      ib= ist2(k)  !<- branch ia and ib 
!
      if(first_ppl) then
        OPEN (unit=11,file=praefixc//'.06'//suffix2,             &
              status='unknown',position='append',form='formatted')
!
        write(11,*) 'k,ist1(k),ist2(k)=',k,ist1(k),ist2(k)
        close(11)
      end if
!
      xp= x(ia)*cph -y(ia)*sph
      yp= x(ia)*sph +y(ia)*cph
      zp= z(ia)
!
      xpp=  xp*cth +zp*sth
      ypp= yp
      zpp= -xp*sth +zp*cth
!
      xx= ps*ypp +hl
      yy= ps*zpp +vd
!
      call plot (xx,yy,3)  ! line starts
!
      xp= x(ib)*cph -y(ib)*sph
      yp= x(ib)*sph +y(ib)*cph
      zp= z(ib)
!
      xpp=  xp*cth +zp*sth
      ypp= yp
      zpp= -xp*sth +zp*cth
!
      xx= ps*ypp +hl
      yy= ps*zpp +vd
      call plot (xx,yy,2)  ! line ends
      end do
      first_ppl= .false.
!
      call plot (xx,yy,3)
      call newcolor(0,1.,0.,0.)  ! Reset
!
!------------------------
!*  Solvent particles.
!------------------------
!
      if(ns.eq.2) then
        do i= nCLp+1,npqr,5
        xp= x(i)*cph -y(i)*sph
        yp= x(i)*sph +y(i)*cph
        zp= z(i)
!
        xpp=  xp*cth +zp*sth
        ypp= yp
        zpp= -xp*sth +zp*cth
!
        xx= ps*ypp +hl
        yy= ps*zpp +vd
!
        dd= 0.5*ps
        call circle (xx-0.15,yy-0.30,dd,2)
        end do
      end if
!
!   Make plots at 1) i=1,np+nq and 2) i=1,npqr 
!---------------------
      call chart
!---------------------
 1000 continue
!
      return
      end subroutine ppl3da 
!
!
!-------------------------------------------------
       subroutine circle (x,y,d,ic)
!-------------------------------------------------
!* Open circle centered at (x,y) /or outer edge.
      use, intrinsic :: iso_c_binding 
      implicit none
!
      integer(C_INT) ic,j,nc
      real(C_float) x,y,d,pi,dth,a,x0,y0,th,x1,y1,x2,y2
!*
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
!
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
!------------------------------------------------
      function iwrta (t,twr)
!------------------------------------------------
      use, intrinsic :: iso_c_binding 
      implicit none
!
      integer(C_INT) iwrta,iwa,iwb,iwc,iwd,iw
      real(C_float)  t,twr
      common/imemo/ iwa,iwb,iwc,iwd
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
      implicit none
!
      integer(C_INT) iwrtb,iwa,iwb,iwc,iwd,iw
      real(C_float)  t,twr
      common/imemo/ iwa,iwb,iwc,iwd
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
      implicit none
!
      integer(C_INT) iwrtc,iwa,iwb,iwc,iwd,iw
      real(C_float)  t,twr
      common/imemo/ iwa,iwb,iwc,iwd
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
      function iwrtd (t,twr)
!------------------------------------------------
      use, intrinsic :: iso_c_binding 
      implicit none
!
      integer(C_INT) iwrtd,iwa,iwb,iwc,iwd,iw
      real(C_float)  t,twr
      common/imemo/ iwa,iwb,iwc,iwd
!
      iw= t/twr 
      if(iw.gt.iwd) then
        iwd= iw
        iwrtd= 0
      else
        iwrtd= 1
      end if
!
      return
      end function iwrtd 
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
!
      use, intrinsic :: iso_c_binding 
      implicit none
!
      real(C_DOUBLE) ranff,invm,x
      integer(C_INT) ir,iq,mask,lambda
      common/ranfff/ ir,iq
!
      parameter  (mask=2**30+(2**30-1),invm= 0.5d0**31)
      parameter  (lambda=48828125)
!
      ir= iand( lambda*ir, mask)
      ranff= ir*invm
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
      subroutine lplot1 (ix,iy,npt1,x,y,ymax,ymin,IL,lab1,n1,lab2,n2, &
                         lab3,n3)
!-----------------------------------------------------------------------
!  <<warning>>  Order and number of arguments /lplot/ have been changed.
!               also, x (time) is defined for all range.
!               date: 5/18/96 at mit.
!***********************************************************************
!   il=1................ linear plot of (x,y)
!   il=2................ log10 plot of (x,log y)
!***********************************************************************
      use, intrinsic :: iso_c_binding 
      implicit none
!
      include   'paramAPGa.h'
!
      integer(C_INT) ix,iy,npt1,IL,n1,n2,n3
      real(C_float)  x(npt1),y(npt1),u(7000),v(7000),   &
                     ymax,ymin,time,t00,xp_leng,        &
                     nfine,pl1,pr1,ql1,qr1,xmin1,xmax1, &
                     ymin1,ymax1,xcm(6),ycm(6),         &
                     pl(6),pr(6),ql(6),qr(6)
      character*8    lab1,lab2,lab3,label,cdate*10
!
      common/headr1/ label,cdate
      common/headr2/ time,t00,xp_leng
      common/pplcom/ nfine,pl1(10),pr1(10),ql1(10),qr1(10), &
                     xmin1(10),xmax1(10),ymin1(10),ymax1(10)
!
!   for fujitsu.
!     data  xcm/18.46,2*9.867,3*6.18/,
!    *      ycm/16.85,2*7.435,3*4.381/,
!    *      pl/2*2.00,15.132,2.00,8.00,18.20/,
!    *      ql/1.95,10.885,1.95,13.832,7.891,1.95/
!
!   for nec.
      data  xcm/21.0, 2*9.00, 3*6.00/,        &
            ycm/15.0, 2*6.80, 3*3.90/,        &
            pl/2.0,  2.0,14.0, 2.0,9.0,16.0/, &
            ql/2.3, 10.5,2.3, 12.9,7.6,2.3/
!
      integer(C_INT) i,j,k,i1,j1,npt,iplot,isc
      real(C_float) hh,hhs,xmax,xmin,dx,dy,x0,y0,scx,scy, &
                    x1,x2,x3,x4,y1,y2,y3,y4,xc,xd,xl,xu,  &
                    yc,yl,yr
!
      iplot=1
      go to 1
!
!-----------------------------------------------------------------------
      entry hplot1 (ix,iy,npt1,x,y,ymax,ymin,il,lab1,n1,lab2,n2,lab3,n3)
!-----------------------------------------------------------------------
      iplot=2
!
    1 npt= npt1
      isc= 1
!
      do i=1,6
      pr(i)= pl(i) +xcm(i)
      end do
!
      do j=1,6
      qr(j)= ql(j) +ycm(j)
      end do
!
!                 ******************************************************
!*                **  make a copy before the top-left frame is drawn. **
!                 ******************************************************
      hh = 0.70
      hhs= 0.60
!
      i1= abs(ix)
      j1= abs(iy)
      if(i1.ge.3) go to 10
      if(j1.eq.3.or.j1.ge.5) go to 10
!                                              ************************
!                                              ** label of the page. **
!                                              ************************
      call symbol (0.1,18.7,hh,label,0.,8)
      call symbol (6.0,18.7,hh,'Date=',0.,5)
      call symbol (999.0,999.0,hh,cdate,0.,10)
!
      call symbol (15.9,0.1,hh,'t=',0.,2)
      call number (999.0,999.0,hh,time,0.,101)
!
   10 continue
!
      do i=1,npt
      u(i)= x(i)
      end do
!
      xmax= u(npt)
      xmin= u(1)
!                             ************************************
!                             ** three-point average if il > 0  **
!                             ************************************
      if(il.gt.0) then
        v(1)=   y(1)
        v(npt)= y(npt)
!
        do i=2,npt-1
        v(i)= 0.33333*(y(i-1)+y(i)+y(i+1))
        end do
      else
        do i=1,npt
        v(i)= y(i)
        end do
      end if
!                                                *****************
!                                                **  log. scale **
!                                                *****************
      if(abs(il).eq.2) then
        do i=1,npt
        if(v(i).gt.0.) then
          v(i)= alog10(v(i))
        else
          v(i)= -10.
        end if
        end do
      end if
!                                **************************************
!                                ** set a new scale and draw a frame.**
!                                **************************************
      if(iplot.eq.2) then
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
      end if
!
      if(ymax.le.ymin) ymax= ymin+1.0
      if(abs(il).eq.2) then
        if(ymax.gt.0.0) ymax= ymax+1.0
      end if
!
      dx= (xmax-xmin)/xcm(i1)
      dy= (ymax-ymin)/ycm(j1)
      x0= xmin
      y0= ymin
!
      call scalex (pl(i1),ql(j1),x0,y0,dx,dy,isc)
!
      pl1(isc)= pl(i1)
      pr1(isc)= pr(i1)
      ql1(isc)= ql(j1)
      qr1(isc)= qr(j1)
      xmin1(isc)= xmin
      xmax1(isc)= xmax
      ymax1(isc)= ymax
      ymin1(isc)= ymin
!                                                      *************
!                                                      **  frame. **
!                                                      *************
      call plot (pl(i1),ql(j1),3)
      call plot (pl(i1),qr(j1),2)
      call plot (pr(i1),qr(j1),2)
      call plot (pr(i1),ql(j1),2)
      call plot (pl(i1),ql(j1),2)
!                                                    ******************
!                                                    **  tick marks. **
!                                                    ******************
      scx= xcm(i1)/4.0
      scy= ycm(j1)/4.0
!
      x0= pl(i1)
      y1= ql(j1)
      y4= qr(j1)
      y2= y1 +0.25
      y3= y4 -0.25
!
      do k=1,3
      x0= x0 +scx
      call plot (x0,y1,3)
      call plot (x0,y2,2)
      call plot (x0,y3,3)
      call plot (x0,y4,2)
      end do
!
      y0= ql(j1)
      x1= pl(i1)
      x4= pr(i1)
      x2= x1 +0.25
      x3= x4 -0.25
!
      do k=1,3
      y0= y0 +scy
      call plot (x1,y0,3)
      call plot (x2,y0,2)
      call plot (x3,y0,3)
      call plot (x4,y0,2)
      end do
!                                                     **************
!                                                     ** numbers. **
!                                                     **************
!
      call number (pl(i1)-1.5,ql(j1)-0.45,hhs,xmin,0.,101) 
      call number (pr(i1)-1.5,ql(j1)-0.45,hhs,xmax,0.,101)
!
      call number (pl(i1)-1.30,ql(j1)     ,hhs,ymin,0.,101)
      call number (pl(i1)-1.30,qr(j1)-0.30,hhs,ymax,0.,101)
!
!                                                     **************
!                                                     **  labels. **
!                                                     **************
      xc= 0.5*(pl(i1)+pr(i1))
!
      xu= xc -1.60
      yr= qr(j1)+0.15
      call symbol (xu,yr,hh,lab1,0.,n1)
!
      xd= xc -0.3*n2/2
      yl= ql(j1)-0.50
      call symbol (xd,yl,hh,lab2,0.,n2)
!
      xl= pl(i1)-1.50
      yc= 0.5*(ql(j1)+qr(j1))
      call symbol (xl,yc,hh,lab3,0.,n3)
!
!                                     **********************************
!                                     **  no plot is made if npt1 < 0 **
!                                     **********************************
   70 if(npt1.lt.0) return
!
      call plotl (u(1),v(1),isc,3)
!**
      if(iplot.eq.1) then
         do i=1,npt
         call plotl (u(i),v(i),isc,2)
         end do
      else
         do i=1,npt-1
         call plotl (u(i+1),v(i)  ,isc,2)
         call plotl (u(i+1),v(i+1),isc,2)
         end do
      end if
!**
      call plotl (u(npt),v(npt),isc,3)
!
      return
      end subroutine lplot1 
!
!
!-----------------------------------------------------------------------
      function igz (iz,k0,wg)
!-----------------------------------------------------------------------
      use, intrinsic :: iso_c_binding 
      implicit none
      include   'paramAPGa.h'
!
      integer(C_INT) igz,iz,k0
      real(C_float)  wg,wgh
      real(C_DOUBLE) gx,gy,gz,ghx,ghy,ghz,ghx2,ghy2,ghz2,  &
                     ghxsq,ghysq,ghzsq,hxi,hyi,hzi
      integer(C_INT) ptx,pty,ptz
      common/xtable/ gx(0:mx),ghx(0:mx),ghx2(0:mx),ghxsq(0:mx), &
                     gy(0:my),ghy(0:my),ghy2(0:my),ghysq(0:my), &
                     gz(0:mz),ghz(0:mz),ghz2(0:mz),ghzsq(0:mz), &
                     hxi,hyi,hzi
      common/xtabl2/ ptx(-10:3000),pty(-10:3000),ptz(-10:4000)
!
!   //Note//
!*   To avoid truncation error, first think of the maximum of
!    gx(i)/hx0 = xmax, then add 1.01 so that the last 1 survive.
!
!*  Plot in every wgh interval.
!
      if(iz.eq.1) then
        wgh= 0.9999 * ghx(1)
        wg= wg + ghx(k0)
!
        if(wg.gt.wgh) then
          wg= amod(wg,wgh)
          igz= 0
        else
          igz= 1
        end if
      end if
!*
      if(iz.eq.2) then
        wgh= 0.9999 * ghy(1)
        wg= wg + ghy(k0)
!
        if(wg.gt.wgh) then
          wg= amod(wg,wgh)
          igz= 0
        else
          igz= 1
        end if
      end if
!*
      if(iz.eq.3) then
        wgh= 0.9999 * ghz(1)
        wg= wg + ghz(k0)
!
        if(wg.gt.wgh) then
          wg= amod(wg,wgh)
          igz= 0
        else
          igz= 1
        end if
      end if
!
      return
      end function igz 
!
!
!-----------------------------------------------------------------------
      subroutine cplot3 (q,xmax8,ymax8,zmax8,char,nc)
!-----------------------------------------------------------------------
!***********************************************************************
!*   Contour plots of scalar quantities.                               *
!***********************************************************************
!     Area of i= 1 to mx-2, etc.
!
      use, intrinsic :: iso_c_binding 
      implicit none

      include   'paramAPGa.h'
!
      real(C_DOUBLE) q(0:mx-1,0:my-1,0:mz-1),xmax8,ymax8,zmax8
      real(C_float)  a(17000),b(17000),ww(17000),cut(200,4),cut2(200), &
                     time,t00,xp_leng,xmax,ymax,zmax,xmin,ymin,zmin
      integer(C_INT) nc,i0,j0,k0
!
      character*8     char,label,cdate*10
      common/headr1/  label,cdate
      common/headr2/  time,t00,xp_leng
!
      real(C_DOUBLE) xleng,yleng,zleng
      common/parm8/  xleng,yleng,zleng
!
      integer(C_INT) pxr,pxc,pxl,pyr,pyc,pyl,pzr,pzc,pzl
      common/ptable/ pxr(-1:mx+1),pxc(-1:mx+1),pxl(-1:mx+1), &
                     pyr(-1:my+1),pyc(-1:my+1),pyl(-1:my+1), &
                     pzr(-1:mz+1),pzc(-1:mz+1),pzl(-1:mz+1)
!
      real(C_DOUBLE) gx,gy,gz,ghx,ghy,ghz,ghx2,ghy2,ghz2, &
                     ghxsq,ghysq,ghzsq,hxi,hyi,hzi
      integer(C_INT) ptx,pty,ptz
      common/xtable/ gx(0:mx),ghx(0:mx),ghx2(0:mx),ghxsq(0:mx), &
                     gy(0:my),ghy(0:my),ghy2(0:my),ghysq(0:my), &
                     gz(0:mz),ghz(0:mz),ghz2(0:mz),ghzsq(0:mz), &
                     hxi,hyi,hzi
      common/xtabl2/ ptx(-10:3000),pty(-10:3000),ptz(-10:4000)
!
      integer(C_INT) igz,i,j,ii,ij,ik,ir,il,npx,npx2,npy,npz,   &
                     jj,jr,jl,kk,k,kr,kl,nxz,ncontr,nxy
      real(C_float) qc,amax7,amin7,xc1,xc2,xl1,xr1,xl2,xr2,hh,      &
                    yc,yr,yl,yl2,yr2,zc,zr,zl,am21,am41,am22,am42,  &
                    ams1,ams2,ams,wamin,wamax,dd,dx,dy,x1,y1,xo,xu, &
                    wgx,wgy,wgz,gdx,gdy,gdz           
!
      real(C_DOUBLE) t8,cptot
      common/headr3/ t8,cptot
!
      i0= mx/2
      j0= my/2 
      k0= mz/2
!
      xmax=  xleng/2.d0 
      ymax=  yleng/2.d0
      zmax=  zleng/2.d0
!
      xmin= -xleng/2.d0
      ymin= -yleng/2.d0
      zmin= -zleng/2.d0
!
!* 1. Plot at constanat y: subscript i first.
!
      qc = 1./16.
      ik= 0
!
      npz= 0
      amax7= -1.e+10
      amin7=  1.e+10
!***
      do k= 1,mz-2
      npz= npz +1
      kr= k+1 
      kl= k-1
!
      cut2(npz)= q(i0,j0,k)
      amax7= max(cut2(npz),amax7)
      amin7= min(cut2(npz),amin7)
!
      npx= 0
!
      do i= 1,mx-2
      npx= npx +1
      ir= i+1
      il= i-1
!
      ik= ik +1
      a(ik)= qc*( q(ir,j0,kr) +2.*q(ir,j0,k)    +q(ir,j0,kl)  &
              +2.*q(i ,j0,kr) +4.*q(i ,j0,k) +2.*q(i ,j0,kl)  &
              +   q(il,j0,kr) +2.*q(il,j0,k)    +q(il,j0,kl) )
      end do
      end do
!
!* 2. Plot at constant z: subscript i first.
!
      qc = 1./16.
      ij= 0
!***
      npy= 0
!
      do j= 1,my-2
      npy= npy +1
      jr= j+1
      jl= j-1
!
      npx2= 0
!
      do i= 1,mx-2
      npx2= npx2 +1
      ir= i+1
      il= i-1
!
      ij= ij +1
      b(ij)= qc*( q(ir,jr,k0) +2.*q(ir,j,k0)    +q(ir,jl,k0)  &
              +2.*q(i ,jr,k0) +4.*q(i ,j,k0) +2.*q(i ,jl,k0)  &
              +   q(il,jr,k0) +2.*q(il,j,k0)    +q(il,jl,k0) )
      end do
      end do
!
!
      hh = 0.70
      call symbol (0.5,18.6,hh,'cplot3',0.,6)
      call symbol (1.8, 1.6,hh,char,0.,nc)
      call symbol (12.0,1.0,hh,label,0.,8)
      call symbol (16.0,1.0,hh,cdate,0.,10)
      call symbol (13.0,0.1,hh,'t=',0.,2)
      call number (999.0,999.0,hh,time,0.,7)
!
! X-Z
      xl1=  1.8
      xr1=  9.3
      zl =  3.0
      zr =  zl +(xr1 -xl1)*zleng/xleng 
      if(zr.gt.18.) zr= 18.
!
      xc1= 0.5*(xr1+xl1)
      zc = 0.5*(zr+zl)
!                          <--- limit elongated z-length.
! X-Y
      xl2= 10.0
      xr2= 17.5
      yr = zr 
      yl = yr -(xr2 -xl2)*yleng/xleng 
      if(yl.lt.3.0) yl= 3.0
!.
      xc2= 0.5*(xr2+xl2)
      yc = 0.5*(yr+yl)
!
!---------------------------------------------
!*  **Maximum of the vectors**
!---------------------------------------------
!
      am21= -777.
      am41= -777.
      am22=  777.
      am42=  777.
!
      do ij= 1,npx*npz
      am21= max(am21,a(ij))
      am22= min(am22,a(ij))
      end do
!
      do ij= 1,npx2*npy
      am41= max(am41,b(ij))
      am42= min(am42,b(ij))
      end do
!
      ams1= max(am21,am41)
      ams2= min(am22,am42)
      ams=  max(ams1,-ams2)
      if(ams.lt.1.e-10) ams=999.0
!
      call symbol  (xl1,0.80,hh,'max=',0.,4)
      call number2 (999.0,999.0,hh,ams1,0.,2)
      call symbol  (xl1,0.10,hh,'min=',0.,4)
      call number2 (999.0,999.0,hh,ams2,0.,2)
!
!---------------------------------------------
!*  (1): Contours in (x,z) plane.
!---------------------------------------------
!
      call setscl (xmin,zmin,xmax,zmax,xl1,zl,xr1,zr,gdx,gdz, &
                   3,'   ',3,'   ',0.4,                       &
                   1,'X',0.6,1,'Z',0.6,1)
!
      call number (xl1-1.3,zl-0.3, hh,zmin,0.,5)
      call number (xl1-1.3,zr-0.3, hh,zmax,0.,5)
      call number (xr1-1.3,zl-0.5, hh,xmax,0.,5)
!
      nxz= npx * npz
      call daisho (a,nxz,wamin,wamax)
!
      if(.false.) then
        OPEN (unit=11,file=praefixc//'.06'//suffix2,             &
              status='unknown',position='append',form='formatted')
        write(11,*)
        write(11,*) 'Contor at (x,z)  t8=',t8 
        write(11,*) 'npx,npz=',npx,npz
        write(11,*) 'xmin,zmin,xmax,zmax=',xmin,zmin,xmax,zmax
        write(11,*) 'xl1,zl,xr1,zr=',xl1,zl,xr1,zr
        write(11,*) 'a(nxz): wamax,wamin,=',wamax,wamin
!
!       do ik= 1,200,5
!       write(11,990) a(ik),a(ik+1),a(ik+2),a(ik+3),a(ik+4)
! 990   format(1p5d12.3)
!       end do
!
!       do ik= nxz-299,nxz,5
!       write(11,990) a(ik),a(ik+1),a(ik+2),a(ik+3),a(ik+4)
!       end do
!       write(11,*)
        close(11)
      end if
!
      if(wamax-wamin.gt.0.) then
        ncontr= 11
        call eqcntr (a,ww,npx,npz,xl1,zl,xr1,zr,wamin,0.0,wamax, &
                     ncontr,1)  !! xxxxxxxxxxxxxx
      end if
!
!---------------------------------------------
!*  (2): Contours in (x,y) plane.
!---------------------------------------------
!
      call setscl (xmin,ymin,xmax,ymax,xl2,yl,xr2,yr,gdx,gdy, &
                   3,'   ',3,'   ',0.4,                       &
                   1,'X',0.6,1,'Y',0.6,1)
!
      call number (xl2-0.45,yl-0.5,hh,ymin,0.,5)
      call number (xl2-1.3, yr-0.3,hh,ymax,0.,5)
      call number (xr2-1.3, yl-0.5,hh,xmax,0.,5)
!
      nxy= npx2*npy
      call daisho (b,nxy,wamin,wamax)
!
!       OPEN (unit=11,file=praefixc//'.06'//suffix2,             &
!             status='unknown',position='append',form='formatted')
!       write(11,*) 'Contor at (x,y)  t8=',t8 
!       write(11,*) 'xmin,ymin,xmax,ymax=',xmin,ymin,xmax,ymax
!       write(11,*) 'xl2,yl,xr2,yr=',xl2,yl,xr2,yr
!       write(11,*) 'b(nxy): wamax,wamin,=',wamax,wamin
!       write(11,*)
!       close(11)
!
      if(wamax-wamin.gt.0.) then
        ncontr= 11
        call eqcntr (b,ww,npx2,npy,xl2,yl,xr2,yr,wamin,0.0,wamax, &
                     ncontr,1)
      end if
!
!---------------------
      call chart
!---------------------
!
      return
      end subroutine cplot3 
!
!
!-----------------------------------------------------------------------
      subroutine filt1p (q,filtx,filty,filtz,sym,sym2)
!-----------------------------------------------------------------------
!***********************************************************************
!*  For duplicate local array  (e/m field).                            *
!***********************************************************************
!   The cell index is i= 0,mx-1, and q(0,mx-1,,)
!
      use, intrinsic :: iso_c_binding 
      implicit none
!
      include  'paramAPGa.h'
!
      real(C_DOUBLE) q(0:mx-1,0:my-1,0:mz-1),sym,sym2
      real(C_float)  a(-2:mx+1,-2:my+1,-2:mz+1)
      integer(C_INT) filtx,filty,filtz,ntx,nty,ntz,i,j,k,is,   &
                     il,ill,ir,irr,js,jl,jll,jr,jrr,ks,kl,kll, &
                     kr,krr  
!
      integer(C_INT)  pxr,pxc,pxl,pyr,pyc,pyl,pzr,pzc,pzl
      common/ptable/  pxr(-1:mx+1),pxc(-1:mx+1),pxl(-1:mx+1), &
                      pyr(-1:my+1),pyc(-1:my+1),pyl(-1:my+1), &
                      pzr(-1:mz+1),pzc(-1:mz+1),pzl(-1:mz+1)
!
      real(C_DOUBLE) gx,gy,gz,ghx,ghy,ghz,ghx2,ghy2,ghz2, &
                     ghxsq,ghysq,ghzsq,hxi,hyi,hzi
      integer(C_INT) ptx,pty,ptz
      common/xtable/ gx(0:mx),ghx(0:mx),ghx2(0:mx),ghxsq(0:mx), &
                     gy(0:my),ghy(0:my),ghy2(0:my),ghysq(0:my), &
                     gz(0:mz),ghz(0:mz),ghz2(0:mz),ghzsq(0:mz), &
                     hxi,hyi,hzi
      common/xtabl2/ ptx(-10:3000),pty(-10:3000),ptz(-10:4000)
!
!***********************************************************************
!*  Filtering in the x direction                                       *
!***********************************************************************
!
      ntx= 0
 1000 ntx= ntx +1
      if(ntx.gt.filtx) go to 2000
!
      do k= 0,mz-1
      do j= 0,my-1
      do i= 0,mx-1
      a(i,j,k)= q(i,j,k)
      end do
      end do
      end do
!                       sym= -1. : odd   -2 -1 (0) 1 2 
!                       sym=  1. : even  
      do is= -1,0
      do k= 0,mz-1
      do j= 0,my-1
      a(is-1,j,k)= sym * q(1-is,j,k)
      end do
      end do
      end do
!                            mx-3 mx-2 (mx-1) mx mx+1 
      do is= 0,1
      do j= 0,my-1
      do k= 0,mz-1
      a(mx+is,j,k)= sym * q(mx-2-is,j,k)
      end do
      end do
      end do
!
      do k= 0,mz-1
      do j= 0,my-1
      do i= 0,mx-1     ! odd for -2,-1,mx,mx+1
      ir = i+1
      il = i-1
      irr= i+2
      ill= i-2
!
      q(i,j,k)= -0.0625*a(ill,j,k) +0.25*a(il,j,k) +0.625*a(i,j,k) &
                  +0.25*a(ir,j,k) -0.0625*a(irr,j,k)
      end do
      end do
      end do
!
      go to 1000
 2000 continue
!
!
      nty= 0
 3000 nty= nty +1
      if(nty.gt.filty) go to 4000
!
      do k= 0,mz-1
      do j= 0,my-1
      do i= 0,mx-1
      a(i,j,k)= q(i,j,k)
      end do
      end do
      end do
!                   
!                    
      do js= -1,0
      do k= 0,mz-1
      do i= 0,mx-1
      a(i,js-1,k)= sym * q(i,1-js,k)
      end do
      end do
      end do
!                
      do js= 0,1
      do i= 0,mx-1
      do k= 0,mz-1
      a(i, my+js,k)= sym * q(i,my-2-js,k)
      end do
      end do
      end do
!
      do k= 0,mz-1
      do j= 0,my-1
      do i= 0,mx-1
      jr = j+1
      jl = j-1
      jrr= j+2
      jll= j-2
!
      q(i,j,k)= -0.0625*a(i,jll,k) +0.25*a(i,jl,k) +0.625*a(i,j,k) &
                 +0.25*a(i,jr,k) -0.0625*a(i,jrr,k)
      end do
      end do
      end do
!
      go to 3000
 4000 continue
!
!
      ntz= 0
 5000 ntz= ntz +1
      if(ntz.gt.filtz) go to 6000
!
      do k= 0,mz-1
      do j= 0,my-1
      do i= 0,mx-1
      a(i,j,k)= q(i,j,k)
      end do
      end do
      end do
!                       sym2= -1: odd   -2 -1 (0) 1 2 
!                       sym2=  1: even  
      do ks= -1,0
      do i= 0,mx-1
      do j= 0,my-1
      a(i,j,ks-1)= sym2 * q(i,j,1-ks)
      end do
      end do
      end do
!                        mz-3 mz-2 (mz-1) mz mz+1
      do ks= 0,1
      do j= 0,my-1
      do i= 0,mx-1
      a(i,j,mz+ks)= sym2 * q(i,j,mz-2-ks)
      end do
      end do
      end do
!
!
      do k= 1,mz-2
      do j= 0,my-1
      do i= 0,mx-1
      kr = k+1
      kl = k-1
      krr= k+2
      kll= k-2
!
      q(i,j,k)= -0.0625*a(i,j,kll) +0.25*a(i,j,kl) +0.625*a(i,j,k) &
                 +0.25*a(i,j,kr) -0.0625*a(i,j,krr)
      end do
      end do
      end do
!
      go to 5000
 6000 continue
!
      return
      end subroutine filt1p 
!
!
!-----------------------------------------------------------------------
      subroutine eqcntr (u,w,nx,ny,xl,yl,xr,yr,umin,ubund,umax,  &
                         lank,iwaku)
!-----------------------------------------------------------------------
!  << eqcntr >>
!         presented by kunihiko.watanabe 14.nov.1989
!         reviced   by hisanori.takamaru 16.mar.1990
!-----------------------------------------------------------------------
!     1. function
!        (1) to draw tokosen
!     2. arguments   (size)   (i/o)     (meaning)
!        (1) u       nx,ny     (i)       world value
!        (2) w       nx,ny     (i)       work array (real*4)
!        (3) xl,xr,yl,yr       (i)       absolute coordinate value
!        (4) umin,umax         (i)       height of max & min
!                                        umin>umax : automatic control
!        (5) ubund             (i)       draw dash line (u < ubund)
!        (6) lank              (i)       number of draw lines
!        (7) iwaku             (i)       =1 : draw frame
!     3. called by
!             (** nothing **)
!     4. calls
!             (** plot   **)
!-----------------------------------------------------------------------
      use, intrinsic :: iso_c_binding 
      real(C_float) u(nx),w(ny)  !<-- nx,ny
!
      if (nx.lt.2) return
      if (ny.lt.2) return
      if (xr.lt.xl) return
      if (yr.lt.yl) return
!
      nxy = nx*ny
      nxm1 = nx - 1
      nym1 = ny - 1
!
      dx = (xr-xl)/ float(nxm1)
      dy = (yr-yl)/ float(nym1)
!
      umax1 = umax
      umin1 = umin
!
      if(umax1.gt.(1.000001*umin1)) then
!
        do i = 1 , nxy
          w(i) = u(i) - umin1
          if(u(i).gt.umax1) w(i) = umax1 - umin1
          if(u(i).lt.umin1) w(i) = 0.
        end do
!
      else
!
        umax1=-1.e+30
        umin1= 1.e+30
        do i = 1 , nxy
          umax1=max(umax1,u(i))
          umin1=min(umin1,u(i))
        end do
!
        do i = 1 , nxy
          w(i) = u(i) - umin1
        end do
!
      endif
!
!------------------------------------------------
      if(umax1.le.(1.000001*umin1))  return
!------------------------------------------------
!
      if(iwaku.eq.1) then
        call plot(xl,yl,3)
        call plot(xr,yl,2)
        call plot(xr,yr,2)
        call plot(xl,yr,2)
        call plot(xl,yl,2)
        call plot(xl,yl,3)
      endif
!
      uld = float(lank+1) / (umax1-umin1)
      eps = 1.0e-8
!
      nxym1 = nxm1*nym1
      do 9000  ijnxy1 = 1,nxym1
        j = (ijnxy1-1)/nxm1 + 1
        i = ijnxy1 - (j-1)*nxm1
!
          i1 = i + nx * (j - 1)
          i2 = i1 + 1
          i3 = i1 + 1 + nx
          i4 = i1 + nx
!
          u1 =  w(i1) * uld
          u2 =  w(i2) * uld
          u3 =  w(i3) * uld
          u4 =  w(i4) * uld
!
          k1 = ifix(u1)
          k2 = ifix(u2)
          k3 = ifix(u3)
          k4 = ifix(u4)
!
!
          j1 = abs(k2-k1)
          j2 = abs(k3-k2)
          j3 = abs(k4-k3)
!
          if(j1.ne.0) then
            do 1000 ll = 1 , j1
              u0 = float(ll) + float(min0(k1,k2))
                ujouge = u0/uld + umin1
                if (ujouge.lt.ubund) then
                  jouge = 4
                else
                  jouge = 1
                end if
!
              if(abs(u2-u1).lt.eps)                 go to 1000
!
              x1 = xl + dx * ( (u0-u1)/(u2-u1) + float(i-1) )
              y1 = yl + dy * float(j-1)
!
              if( ((u3-u0)*(u2-u0)).gt.0. )         go to 1100
              if( ( (u0-u2).gt.0. ).and.( (u0-u4).gt.0. ) ) go to 1100
              if( abs(u3-u2).lt.eps )               go to 1100
!
                x2 = xl + dx * float(i)
                y2 = yl + dy * ( (u0-u2)/(u3-u2) + float(j-1) )
!
                call wdash(x1,y1,x2,y2,jouge)
!
 1100         continue
              if( ((u4-u0)*(u3-u0)).gt.0. )         go to 1200
              if( ((u1-u0)*(u3-u0)).gt.0. )         go to 1200
              if( ((u2-u0)*(u4-u0)).gt.0. )         go to 1200
              if( abs(u4-u3).lt.eps )               go to 1200
!
                x2 = xl + dx * ( (u0-u4)/(u3-u4) + float(i-1) )
                y2 = yl + dy * float(j)
                call wdash(x1,y1,x2,y2,jouge)
!
 1200         continue
              if( ((u1-u0)*(u4-u0)).gt.0. )         go to 1300
              if( ( (u0-u1).gt.0. ).and.( (u0-u3).gt.0. ) ) go to 1300
              if( abs(u1-u4).lt.eps )               go to 1300
!
                x2 = xl + dx * float(i-1)
                y2 = yl + dy*((u0-u1)/(u4-u1)+float(j-1))
                call wdash(x1,y1,x2,y2,jouge)
 1300         continue
 1000       continue
!
          endif
!
          if(j2.ne.0) then
!
            do 2000 ll = 1 , j2
              u0 = float(ll) + float(min0(k2,k3))
                ujouge = u0/uld + umin1
                if (ujouge.lt.ubund) then
                  jouge = 4
                else
                  jouge = 1
                end if
              if( abs(u3-u2).lt.eps )               go to 2000
!
              x1 = xl + dx * float(i)
              y1 = yl + dy * ( (u0-u2)/(u3-u2) + float(j-1) )
!
              if( ((u4-u0)*(u3-u0)).gt.0. )         go to 2100
              if( ( (u0-u1).gt.0. ).and.( (u0-u3).gt.0. ) ) go to 2100
              if( abs(u4-u3).lt.eps )               go to 2100
!
                x2 = xl + dx * ( (u0-u4)/(u3-u4) + float(i-1) )
                y2 = yl + dy * float(j)
!
                call wdash(x1,y1,x2,y2,jouge)
!
 2100         continue
              if( ((u1-u0)*(u4-u0)).gt.0. )         go to 2200
              if( ((u1-u0)*(u3-u0)).gt.0. )         go to 2200
              if( ((u2-u0)*(u4-u0)).gt.0. )         go to 2200
              if( abs(u1-u4).lt.eps )               go to 2200
!
                x2 = xl + dx * float(i-1)
                y2 = yl + dy * ( (u0-u1)/(u4-u1)+float(j-1) )
                call wdash(x1,y1,x2,y2,jouge)
 2200         continue
 2000       continue
!
          endif
!
          if(j3.ne.0) then
!
            do 3000 ll = 1 , j3
              u0 = float(ll) + float(min0(k3,k4))
                ujouge = u0/uld + umin1
                if (ujouge.lt.ubund) then
                  jouge = 4
                else
                  jouge = 1
                end if
              if( abs(u4-u3).lt.eps )               go to 3000
!
              x1 = xl + dx * ( (u0-u4)/(u3-u4) + float(i-1) )
              y1 = yl + dy * float(j)
!
              if( ((u1-u0)*(u4-u0)).gt.0. )         go to 3100
              if( ( (u0-u2).gt.0. ).and.( (u0-u4).gt.0. ) ) go to 3100
              if( abs(u1-u4).lt.eps )               go to 3100
!
                x2 = xl + dx * float(i-1)
                y2 = yl + dy * ( (u0-u1)/(u4-u1) + float(j-1) )
                call wdash(x1,y1,x2,y2,jouge)
 3100         continue
 3000       continue
          endif
 9000 continue
!
      return
      end subroutine eqcntr 
!
!
!------------------------------------------------------------------------
      subroutine setscl (wminx,wminy,wmaxx,wmaxy, xl,yl,xr,yr, gdx,gdy, &
                         n1,char1, n2,char2, hight1,                    &
                         nnx,charx, hightx, nny,chary, highty, iwaku)
!------------------------------------------------------------------------
!  << setscl >>                   /char1/
!                          wmaxy  +--------------------+  (xl,yl)
!                               y |             (xr,yr)|  (xr,yr) on 0
!                               r |                    |
!                               a |                    |
!    (wminx,wminy)              h |                    |
!    (wmaxx,wmaxy) on is        c |                    |
!                                 |(xl,yl)             |
!                          wminy  +--------+--+--------+
!                                 wminx  /charx/       wmaxx
!-----------------------------------------------------------------------
!
!     setscl
!
!     1. function
!        (1) to scale the graphics by calcomp specifications
!     2. arguments            (i/o)     (meaning)
!        (1) wminx,wmaxx,
!            wminy,wmaxy       (i)       world coordinate value
!        (2) xl,xr,yl,yr       (i)       absolute coordinate value
!        (3) gdx,gdy           (o)       scaling factor of coordinate
!                                        from world to absolute
!        (4) char1,charx,cahry (i)       title on graph,x-axis,y-axis
!        (5) iwaku             (i)       draw frame (0:off ; 1:on)
!                                         999 : write out only title,
!                                                    not draw otherwise
!     3. called by
!             (** nothing **)
!     4. calls
!             (** plot   **)
!             (** symbol **)
!             (** number **)
!-----------------------------------------------------------------------
      use, intrinsic :: iso_c_binding 
      implicit none
!
      include     'paramAPGa.h'
      character*1  char1(n1),char2(n2),charx(nnx),chary(nny)
!
      integer(C_INT) n1,n2,n3,iwaku,nnx,nny,nnx1,nny1
      real(C_float)  wmaxx,wminx,wmaxy,wminy,xr,xl,yr,yl, &
                     gdx,gdy,gdz,xc,yc,zc,hight1,xs1,xs2, &
                     hightx,highty,fctr
!
      integer(C_INT) io_pe
      common/sub_proc/ io_pe
!
      if (wmaxx.le.wminy) goto 9999
      if (wmaxx.le.wminy) goto 9999
      if (xr.le.xl)       goto 9999
      if (yr.le.yl)       goto 9999
!
      gdx= (xr-xl)/(wmaxx-wminx)
      gdy= (yr-yl)/(wmaxy-wminy)
!
      xc = 0.5*( xr + xl )
      yc = 0.5*( yr + yl )
!
      if (n1 .gt.0) then
        if (hight1.gt.0) then
          xs1= xc -0.5*n1*hight1
          xs2= xs1 +(n1+1)*hight1
          call symbol(xs1,yr+0.1,hight1,char1(1),0.,n1)
          call symbol(xs2,yr+0.1,hight1,char2(1),0.,n2)
        end if
      end if
!-----------------------------------------------------------------------
      if (iwaku.eq.999) return
!-----------------------------------------------------------------------
!
      if (iwaku.eq.1) then
        call plot (xl,yl,3)
        call plot (xl,yr,2)
        call plot (xr,yr,2)
        call plot (xr,yl,2)
        call plot (xl,yl,2)
        call plot (999.,999.0,3)
      end if
!
      if (nnx.gt.0) then
        if (hightx.gt.0) then
          call symbol(xc-0.5*hightx*nnx,yl-0.5,hightx,charx(1),0.,1)
          do nnx1=2,nnx
          call symbol(999.0,999.0,hightx,charx(nnx1),0.,1)
          end do
        end if
      end if
!
      if (nny.gt.0) then
        if (highty.gt.0) then
          call symbol(xl-0.5,yc-0.5*highty*nny,highty,chary(1),0.,1)
          do nny1=2,nny
          call symbol(999.0,999.0,highty,chary(nny1),0.,1)
          end do
        end if
      else if(nny.lt.0) then
        if (highty.gt.0) then
          call symbol(xc-0.5*highty*nny,yc,highty,chary(1),0.,1)
          do nny1=2,nny
          call symbol(999.0,999.0,highty,chary(nny1),0.,1)
          end do
        end if
      end if
!
      return
!
!-----------------------------------------------------------------------
!
 9999 continue
      if(io_pe.eq.1) then
      OPEN (unit=11,file=praefixc//'.06'//suffix2,             &
            status='unknown',position='append',form='formatted')
!
      write(11,*) '**********  abnormal world coordinate ********'
      write(11,*) '      '
      write(11,*) '    wmaxx =',wmaxx,' wminx = ',wminx
      write(11,*) '    wmaxy =',wmaxy,' wminy = ',wminy
      write(11,*) '    xl,yl,xr,yr =',xl,yl,xr,yr
      write(11,*) '    fctr  =',fctr
      write(11,*) '      '
      close(11)
      end if
!
      call chart
      call symbol(1.0,10.0,0.2,' abnormal world coordinate call',0.,31)
      call symbol(1.0,09.0,0.2,' wmaxx =',0.,8)
      call number(999.0,999.0,0.2,wmaxx,0.,2)
      call symbol(1.0,08.5,0.2,' wminx =',0.,8)
      call number(999.0,999.0,0.2,wminy,0.,2)
      call symbol(1.0,08.0,0.2,' wmaxy =',0.,8)
      call number(999.0,999.0,0.2,wmaxy,0.,2)
      call symbol(1.0,07.5,0.2,' wminy =',0.,8)
      call number(999.0,999.0,0.2,wminy,0.,2)
      call symbol(1.0,07.0,0.2,' fctr  =',0.,8)
      call number(999.0,999.0,0.2,fctr,0.,2)
      call symbol(1.0,06.5,0.2,' xleft =',0.,8)
      call number(999.0,999.0,0.2,xl,0.,2)
      call symbol(1.0,06.0,0.2,' yleft =',0.,8)
      call number(999.0,999.0,0.2,yl,0.,2)
      call symbol(1.0,05.5,0.2,' xright=',0.,8)
      call number(999.0,999.0,0.2,xr,0.,2)
      call symbol(1.0,05.0,0.2,' yright=',0.,8)
      call number(999.0,999.0,0.2,yr,0.,2)
      stop
!
      return
      end subroutine setscl 
!
!
!-----------------------------------------------------------------------
      subroutine scalex (xcm,ycm,x00,y00,dx,dy,isc)
!-----------------------------------------------------------------------
      use, intrinsic :: iso_c_binding 
      implicit none
!
      integer(C_INT) isc
      real(C_float)  xcm,ycm,x00,y00,dx,dy,x0,y0,xl,yl,dxi,dyi
      common/gscale/ x0(10),y0(10),xl(10),yl(10),dxi(10),dyi(10)
!
      x0(isc)= x00
      y0(isc)= y00
      dxi(isc)= 1./dx
      dyi(isc)= 1./dy
!
      xl(isc)= xcm
      yl(isc)= ycm
!
      return
      end subroutine scalex 
!
!
!-----------------------------------------------------------------------
      subroutine plotl (x,y,isc,ipl)
!-----------------------------------------------------------------------
      use, intrinsic :: iso_c_binding 
      implicit none
!
      integer(C_INT) isc,ipl
      real(C_float)  x,y,x0,y0,xl,yl,dxi,dyi,xcm,ycm
      common/gscale/ x0(10),y0(10),xl(10),yl(10),dxi(10),dyi(10)
!
      xcm= xl(isc) +dxi(isc)*(x -x0(isc))
      ycm= yl(isc) +dyi(isc)*(y -y0(isc))
!
      call plot (xcm,ycm,ipl)
!
      return
      end subroutine plotl 
!
!
!-----------------------------------------------------------------------
      subroutine values (x,y,height,val,theta,ifmat)
!-----------------------------------------------------------------------
!  << values >>
!     1. function
!        (1) to draw variable
!     2. arguments   (size)   (i/o)     (meaning)
!        (1) x,y               (i)       absolute coordinate value
!        (2) height            (i)       draw out size on paper
!        (3) val               (i)       variable
!        (4) theta             (i)       angle
!        (5) ifmat             (i)       format type
!     3. called by
!             (** nothing **)
!     4. calls
!             (** number **)
!             (** symbol **)
!-----------------------------------------------------------------------
!        ifmat = (n100)*100 + keta
!        n100 = 0 : integer format
!        n100 = 1 : f format ::  number(x,y,height,val,theta,keta)
!        n100 = 2 : e format ::
!        n100 = 3 : power of ten format
!        n100 = othewise : not write out
!-----------------------------------------------------------------------
      use, intrinsic :: iso_c_binding 
!
      integer(C_INT) ifmat
      real(C_float)  x,y,height,val,theta
!
      character   chr13*13,chr12*12,chr3*3
      character*1 minus,zero,blank
      parameter(ratio = 6./7. )
      data minus/'-'/,zero/'0'/,blank/' '/
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      if (ifmat.lt.0) return
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      n100 = ifmat/100
      keta = ifmat - n100*100
!
      if (n100.eq.0) then
        call number(x,y,height,val,theta,ifmat)
!       call number(x,y,height,val,theta,-1)
      else if (n100.eq.1) then
        call number(x,y,height,val,theta,ifmat)
!       call number(x,y,height,val,theta,keta)
      else if (n100.eq.2) then
        chr13 = '             '
        chr12 = '            '
        if (keta.eq.0) then
          write(chr13,'(1pe13.6)') val
          chr12(1:4) = chr13(1:3)//'e'
          numsym = 4
        else
          keta = keta + 1
          if (val.lt.0.) then
            chrval = val - 5.*10**float(-keta)
            write(chr13,'(1pe13.6)') chrval
            chr12(1:keta+3) = chr13(1:keta+2)//'e'
            numsym = keta + 3
          else if (val.eq.0) then
            chrval = val
            write(chr13,'(1pe13.6)') chrval
            chr12(1:keta+3) = chr13(1:keta+2)//'e'
            numsym = keta + 3
          else
            chrval = val + 5.*10**float(-keta)
            write(chr13,'(1pe13.6)') chrval
            chr12(1:keta+2) = chr13(2:keta+2)//'e'
            numsym = keta + 2
          end if
        end if
        chr3 = '   '
!
        if (chr13(11:11) .eq. minus) then
          if (chr13(12:12) .eq. zero  .or. &
              chr13(12:12) .eq. blank) then
            chr3(1:2) = '-'//chr13(13:13)
          else
            chr3(1:3) = '-'//chr13(12:13)
          end if
          numsy1 = 3
        else
          if (chr13(12:12) .eq. zero  .or. &
              chr13(12:12) .eq. blank) then
            chr3(1:1) = chr13(13:13)
            numsy1 = 1
          else
            chr3(1:2) = chr13(12:13)
            numsy1 = 2
          end if
        end if
        akaku = 2. * 3.1415927 / 360.
        cost = cos(theta*akaku)
        call symbol(x,y,height,chr12,theta,numsym)
        call symbol(999.,999.,height,chr3,theta,numsy1)
      else if (n100.eq.3) then
        chr13 = '             '
        chr12 = '            '
        if (keta.eq.0) then
          write(chr13,'(1pe13.6)') val
          chr12(1:6) = chr13(1:3)//'x10'
          numsym = 6
        else
          keta = keta + 1
          if (val.lt.0.) then
            chrval = val - 5.*10**float(-keta)
            write(chr13,'(1pe13.6)') chrval
            chr12(1:keta+5) = chr13(1:keta+2)//'x10'
            numsym = keta + 5
          else
            chrval = val + 5.*10**float(-keta)
            write(chr13,'(1pe13.6)') chrval
            chr12(1:keta+4) = chr13(2:keta+2)//'x10'
            numsym = keta + 4
          end if
        end if
        chr3 = '   '
!
        if (chr13(11:11) .eq. minus) then
          if (chr13(12:12) .eq. zero  .or. &
              chr13(12:12) .eq. blank) then
            chr3(1:2) = '-'//chr13(13:13)
          else
            chr3(1:3) = '-'//chr13(12:13)
          end if
          numsy1 = 3
        else
          if (chr13(12:12) .eq. zero  .or. &
              chr13(12:12) .eq. blank) then
            chr3(1:1) = chr13(13:13)
            numsy1 = 1
          else
            chr3(1:2) = chr13(12:13)
            numsy1 = 2
          end if
        end if
        akaku = 2. * 3.1415927 / 360.
        cost = cos(theta*akaku)
        sint = sin(theta*akaku)
        call symbol(x,y,height,chr12,theta,numsym)
!
!                                             *******************
!                                             ** exponent part **
!                                             *******************
!
        h2 = height * 5./7.
        x1 = (numsym+1)* height * ratio
        y1 = height * 4./7.
        if (abs(theta).lt.1e-04) then
          x1 = x + x1
          y1 = y + y1
        else
          x2 =     x1 * cost - y1 * sint
          y1 = y + x1 * sint + y1 * cost + h2*cost
          x1 = x + x2                    - h2*sint
        end if
        call symbol(x1,y1,h2,chr3,theta,numsy1)
      end if
      return
      end subroutine values 
!
!
!-----------------------------------------------------------------------
      subroutine wdash (x1,y1,x2,y2,ipen )
!-----------------------------------------------------------------------
!  << wdash  >>                      ver 2.00   16.mar.1990
!
!     1. function
!        (1) to draw line from (x1,y1) to (x2,y2) by wdash
!                            in absolute coordinate
!     2. arguments            (i/o)     (meaning)
!        (1) x1,x2,y1,y2       (i)       absolute coordinate value
!        (2) ipen              (i)       pen type of 'wdash'
!     3. called by
!             (** eqcntr  **)
!             (** wdashl  **)
!     4. calls
!             (** plot   **)
!-----------------------------------------------------------------------
!       ipen : meaning           - : 0.05 (cm)
!        1   :       line     -------------------
!        2   :  dash line     --- --- --- --- ---
!        3   :  dash line     -- -- -- -- -- -- --
!        4   :  dash line     - - - - - - - - - -
!        5   :  1 point dash  ---- - ---- - ---- -
!        6   :  2 point dash  --2.0-- - - --2.0--
!   otherwise:  line          ---------------------
!-----------------------------------------------------------------------
      use, intrinsic :: iso_c_binding 
!
      integer(C_INT) ipen
      real(C_float)  x1,y1,x2,y2
!
      h1  =  0.05
      h2  =  2.0 * h1
      h3  =  3.0 * h1
      h4  =  4.0 * h1
      h20 = 20.0 * h1
      call plot ( x1 , y1 , 3 )
      k = - 1
      if(ipen.lt.2) then
        go to 999
      else if(ipen.eq.2) then
        hh1 = h3
        hh2 = h1
      else if (ipen.eq.3) then
        hh1 = h2
        hh2 = h1
      else if (ipen.eq.4) then
        hh1 = h1
        hh2 = h1
      else if (ipen.eq.5) then
        hh1 = h4
        hh2 = h1
        hh3 = h1
        hh4 = h1
      else if (ipen.eq.6) then
        hh1 = h20
        hh2 = h1
        hh3 = h1
        hh4 = h1
        hh5 = h1
        hh6 = h1
      end if
      if(ipen.lt.5) then
        rleng = sqrt ( ( x2 - x1 ) **2 + ( y2 - y1 ) **2 )
        if(rleng.lt.1.0e-5) goto 999
        if(rleng.lt.hh1) goto 999
        costh = ( x2 - x1 ) / rleng
        sinth = ( y2 - y1 ) / rleng
        d = hh1
        x = x1 + d * costh
        y = y1 + d * sinth
        call plot ( x , y , ( 5 + k ) / 2 )
        k = - k
        d = d + hh2
        hhh = hh1
        hh1 = hh2
        hh2 = hhh
  200   if(d.le.rleng) then
          x = x1 + d * costh
          y = y1 + d * sinth
          call plot ( x , y , ( 5 + k ) / 2 )
          k = - k
          hhh = hh1
          hh1 = hh2
          hh2 = hhh
          d=d+hh1
          goto 200
        end if
      else if (ipen.eq.5) then
        rleng = sqrt ( ( x2 - x1 ) **2 + ( y2 - y1 ) **2 )
        if(rleng.lt.1.0e-5) goto 999
        if(rleng.lt.hh1) goto 999
        costh = ( x2 - x1 ) / rleng
        sinth = ( y2 - y1 ) / rleng
        d = hh1
        x = x1 + d * costh
        y = y1 + d * sinth
        call plot ( x , y , ( 5 + k ) / 2 )
        k = - k
        d = d + hh2
        hhh = hh1
        hh1 = hh2
        hh2 = hh3
        hh3 = hh4
        hh4 = hhh
  500   if(d.le.rleng) then
          x = x1 + d * costh
          y = y1 + d * sinth
          call plot ( x , y , ( 5 + k ) / 2 )
          k = - k
          hhh = hh1
          hh1 = hh2
          hh2 = hh3
          hh3 = hh4
          hh4 = hhh
          d=d+hh1
          goto 500
        end if
      else if (ipen.eq.6) then
        rleng = sqrt ( ( x2 - x1 ) **2 + ( y2 - y1 ) **2 )
        if(rleng.lt.1.0e-5) goto 999
        if(rleng.lt.hh1) goto 999
        costh = ( x2 - x1 ) / rleng
        sinth = ( y2 - y1 ) / rleng
        d = hh1
        x = x1 + d * costh
        y = y1 + d * sinth
        call plot ( x , y , ( 5 + k ) / 2 )
        k = - k
        d = d + hh2
        hhh = hh1
        hh1 = hh2
        hh2 = hh3
        hh3 = hh4
        hh4 = hh5
        hh5 = hh6
        hh6 = hhh
  600   if(d.le.rleng) then
          x = x1 + d * costh
          y = y1 + d * sinth
          call plot ( x , y , ( 5 + k ) / 2 )
          k = - k
          hhh = hh1
          hh1 = hh2
          hh2 = hh3
          hh3 = hh4
          hh4 = hh5
          hh5 = hh6
          hh6 = hhh
          d=d+hh1
          goto 600
        end if
      end if
  999 call plot ( x2 , y2 , ( 5 + k ) / 2 )
      call plot ( x2 , y2 , 3)
!
      return
      end subroutine wdash 
!
!
!-----------------------------------------------------------------------
      subroutine daisho (x,nx,xmin1,xmax1)
!-----------------------------------------------------------------------
      use, intrinsic :: iso_c_binding 
!
      integer(C_INT) nx
      real(C_float)  x(nx),xmin1,xmax1
!
      xmax1= -1.d30 
      xmin1=  1.d30 
!
      do i= 1,nx
      xmax1= max(xmax1,x(i))
      xmin1= min(xmin1,x(i))
      end do
!
      return
      end subroutine daisho
!
!
!***************************************************************
!*  This program package generates a UNIX postscript           *
!*  graphic file when called by calcomp-compatible /plot23.f/  *
!***************************************************************
!---------------------------------------------------------------
!    PostScript header by fortran
!      T. Ogino (Nagoya University) February 27, 1992
!    Modified to conform GSIPP commands
!      Motohiko Tanaka (NIFS)       November 23, 1993
!
!----------------------------------------------- 5/27/1996 -----
!   This PS-Adobe-2.0 header allows us full paging features in
!   the Ghostview.  To scroll up the page (backward), click the 
!   page number and press two buttons of mouse simultaneously.
!
!   Consult: A.Saitou (Kyoto U.)  The definition of /@eop  
!   needs stroke for line drawings (not in the TeX header).
!
!---------------------------------------------------------------
       subroutine gopen (nframe)
!----------------------------------------------------------
!*  This is an adobe-2.0 postscript file.
!
       use, intrinsic :: iso_c_binding 
       implicit none
!
       integer(C_INT) nframe,n0,ipage,nfrm
       real(C_float)  fmag,x0,y0,h0,xcm,xwc
!
       common/convsn/ fmag,x0,y0,h0,n0
       common/pages/  ipage,nfrm
!
       write(77,10)
   10  format('%!ps-adobe-2.0',/      &
              '%%pages: (atend)',/    &
              '%%pageorder: ascend',/ &
              '%%endcomments',/       &
              '%%begindocument')
!
!%%%%%%%%%%%%%%%%%%% Procedure defintions %%%%%%%%%%%%%%%%%%%%%%%%%%
!
!     write(77,11) 
!  11 format('%%boundingbox: 150. 400. 550. 600.')
!
      write(77,21) 
   21 format('/l {lineto} bind def  % x y l -- line to position',/ &
             '/m {moveto} bind def  % x y m -- move to position')
!
      write(77,23) 
   23 format('/tr {/times-roman findfont} bind def',/ &
             '/sf {scalefont} bind def',/ &            
             '/se {setfont} bind def',/   &
             '/ro {rotate}  bind def',/   &
             '/tl {translate} bind def',/ &
             '/sc {scale} bind def')
!
      write(77,24) 
   24 format('/@bop          % @bop -- begin the a new page',/ &
             '{erasepage newpath initgraphics',/ &
             '/saveimage save def',/             &
             '} bind def')
!
      write(77,25) 
   25 format('/@eop          % @eop -- end a page',/ &
             '{stroke showpage',/   &
             ' saveimage restore',/ &
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
             ' {gsave currentpoint translate 90 rotate 0 0 moveto',/ &
             ' show grestore}',/ &
             ' {show} ifelse',/  &
             '} bind def')
!
      write(77,31)
   31 format('%%enddocument',/  &
             '%%endprolog',/    &
             '%%beginsetup',/   &
             '/resolution 300 def',/ &
             '/#copies 1 def',/ &
             '%%endsetup')
!
!%%%%%%%%%%%%%%%%%%% End of the header %%%%%%%%%%%%%%%%%%%%%%%%%%
!
!*  Initiate the page one.
!
       nfrm = nframe
!
       ipage = 1
       write(77,12) ipage,ipage
   12  format('%%page:',1x,i2,1x,i2)
!
       write(77,30) 
   30  format('%%beginpagesetup',/ &
              '%%endpagesetup',/   &
              '@bop')
!
!*  Set magnifying factor (gsipp to sun coordinate).
!   rotate and translate to output on a4-L paper.
!      left corner ...... (  0.,  0.)
!      right corner ..... (600.,780.)
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
       end subroutine gopen
!
!
!-----------------------------
       subroutine gclose
!-----------------------------
       use, intrinsic :: iso_c_binding 
       implicit none
!
       call plote
!
       return
       end subroutine gclose
!
!
!-----------------------------
       subroutine plote
!-----------------------------
       use, intrinsic :: iso_c_binding 
       implicit none
!
       write(77,10) 
   10  format('@eop')
!
       return
       end subroutine plote
!
!
!-----------------------------------------
       subroutine chart
!-----------------------------------------
!*  four frames in a page (if nfrm=4).
       use, intrinsic :: iso_c_binding 
       implicit none
!
       integer(C_INT) ipage,nfrm,loc,lpage
       common/pages/  ipage,nfrm
!
       ipage = ipage +1
       loc= mod(ipage-1,nfrm)
!
!*  frame 1: open a new page.
!
       if(loc.eq.0) then
          call plote
!
          if(nfrm.eq.1) lpage= ipage
          if(nfrm.ne.1) lpage= (ipage+3)/4
!
          write(77,10) 
   10     format('%%pagetrailer    % need for the page count')
!
          write(77,20) lpage,lpage
   20     format('%%page:',1x,i2,1x,i2)
!
          write(77,30) 
   30     format('%%beginpagesetup',/ &
                 '%%endpagesetup',/   &
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
!   First cancel the previous translation, then
!   make a new translation (scale factor alive).
!-----------------------------------------------------
!*  Frames 2-4:
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
       end subroutine chart
!
!
!------------------------------------
       subroutine factor(fct)
!------------------------------------
       use, intrinsic :: iso_c_binding 
       implicit none
!
       real(C_float) fct
!
       write(77,10) fct,fct
   10  format(f6.2,1x,f6.2,' sc')
!
       return
       end subroutine factor
!
!
!------------------------------------
       subroutine newpen (ip)
!------------------------------------
       use, intrinsic :: iso_c_binding 
       implicit none
!
       integer(C_INT) ip,i1,i2
       real(C_float)  pi1
!
       i1=(ip-1)/2
       i2=ip-2*i1
!
       write(77,*) 'sn'
       pi1=0.40*float(i1-1)
       write(77,30) pi1
   30  format(f3.1,' sl')
!
       if(i2.ne.1) then
       write(77,*) '[2 2] 0 sd'
       endif
!
       return
       end subroutine newpen
!
!
!---------------------------------------
       subroutine newcolor (ic,r,g,b)
!---------------------------------------
!  ic= 3 tri-color
!  ic= 0 gray scale, r= 0. for black
!     black specified is... newcolor (0,1.,0.,0.)
!
       use, intrinsic :: iso_c_binding 
       implicit none
!
       integer(C_INT) ic
       real(C_float)  r,g,b
!
       write(77,*) 'stroke'
!
       if(ic.eq.0) then
         write(77,10) 1.-r  ! 0. for black
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
       use, intrinsic :: iso_c_binding 
       implicit none
!
       write(77,*) 'st'
!
       return
       end subroutine linee
!
!
!------------------------------------
       subroutine plot (x0,y0,ip)
!------------------------------------
       use, intrinsic :: iso_c_binding 
       implicit none
!
       integer(C_INT) ip,n
       real(C_float)  x0,y0,x,y,h
!
       x= x0
       y= y0
       h= 0.
       n= 777
       call sunscl (x,y,h,n)
!
       if(ip.eq.3)  write(77,10) x,y
   10  format(f5.1,1x,f5.1,' m')
!
       if(ip.eq.2)  write(77,20) x,y
   20  format(f5.1,1x,f5.1,' l')

       if(ip.eq.-3) write(77,30) x,y
   30  format(f5.1,1x,f5.1,' tl')
!
       if(ip.eq.-2) write(77,40) x,y,x,y
   40  format(f5.1,1x,f5.1,' l sn',1x,f5.1,1x,f5.1,' tl')
!       write(77,*) 'st'
!
       return
       end subroutine plot
!
!
!-------------------------------------------------
       subroutine symbol (x0,y0,h0,isymb,ang,n0)
!-------------------------------------------------
       use, intrinsic :: iso_c_binding 
       implicit none
!
       integer(C_INT) n0,n,i
       real(C_float)  x0,y0,h0,ang,x,y,z,h
!
       character(*)   isymb
       character ica*80,ich(80)*1
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
       end subroutine symbol
!
!
!-----------------------------------------------
       subroutine number (x0,y0,h0,anu,ang,n0)
!-----------------------------------------------
       use, intrinsic :: iso_c_binding 
       implicit none
!
       real(C_float)  x0,y0,h0,anu,ang,x,y,h
       integer(C_INT) n0,n
       character      isymb*9
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
       if(abs(anu).lt.10000.) then  ! 5 digits
         if(abs(anu).gt.0.1) then
           write(isymb,40) anu
   40      format(f6.1)
         else
           if(abs(anu).gt.0.001) then  ! f6.3
             write(isymb,41) anu
   41        format(f6.3)
           else
             if(abs(anu).gt.0.001) then  ! f6.3
               write(isymb,42) anu   ! e9.2
   42          format(1pe9.2)
             else
               write(isymb,40) anu   ! 0.0
             end if
           end if
         end if
!
       else
         if(abs(anu).lt.100000.) then
           write(isymb,51) anu     ! f7.1
   51      format(f7.1)
         else
           write(isymb,52) anu     ! e9.2
   52      format(1pe9.2)
         end if
       end if
!
       write(77,*) '(',isymb,') s'
!
       return
       end subroutine number
!
!-----------------------------------------------
       subroutine number2 (x0,y0,h0,anu,ang,n0)
!-----------------------------------------------
       use, intrinsic :: iso_c_binding 
       implicit none
!
       integer(C_INT) n0,n
       real(C_float)  x0,y0,h0,anu,ang,x,y,h
       character      isymb*9
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
      if(n0.eq.1) write(isymb,40) anu
      if(n0.eq.2) write(isymb,41) anu
   40  format(f6.1)
   41  format(f6.3)
!
       write(77,*) '(',isymb,') s'
!
       return
       end subroutine number2
!
!
!---------------------------------------------------
       subroutine sunscl (x,y,h,n)
!---------------------------------------------------
       use, intrinsic :: iso_c_binding 
       implicit none
!
       integer(C_INT) n,n0
       real(C_float)  x,y,h,fmag,x0,y0,h0 
       common/convsn/ fmag,x0,y0,h0,n0
!
       if(x.eq.999.) then
         x= x0 +abs(n0)*h0
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
       end subroutine sunscl
