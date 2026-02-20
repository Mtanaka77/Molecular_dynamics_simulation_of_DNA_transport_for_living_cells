!************************************************** 2025/07/13 ***
!*                                                               *
!*   ## Molecular Dynamics for ElectroStatic Living Cells ##     *
!*                                                               *
!*     @nanoWatPa.f03 - Short-range Coulomb and LJ forces, and   *
!*     long-range Poisson forces for electrostatic potentials.   *
!*     Constant electric field is applied in the z-direction.    *
!*                                                               *
!*   Author: Motohiko Tanaka, Ph.D.                              *
!*           Nature and Science Applications, Nagoya 464, Japan. *
!*   This code is permitted by GNU General Public License v3.0.  *
!*                                                               *
!*     The short-range Coulomb forces and the long-range         *
!*   electrostatic effects by the Poisson equation are treated   *
!*   in this simulation code.                                    *
!*                                                               *
!*     This molecular dynamics simulation of nanoscale pores     *
!*   in 3D electrostatic effects was updated in April 2025.      *
!*   The initial configuration of water and hydrate is made by   *
!*   TIP4 model in Dr.M.Matsumoto, https://github.com/vitroid.   *  
!*                                                               *
!*   References:                                                 *
!*   1) Y.Rabin and M.Tanaka, DNA in nanopores: Counterion cond- * 
!*    ensation and..., Phys.Rev.Letters, vol.94, 148103 (2005).  *
!*   2) M. Tanaka, https://github.com/Mtanaka77/                 *
!*    Electrostatic_molecular_dynamics_of_living_human_cells     *
!*                                                               *
!*---------------------------------------------------------------*
!*  >> Note: isize must be chosen such that the sub-box size     *
!*            Lx/isize > ag(counter) +ag(water)                  *
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
!*     mass........ m= 18*1.67 10^-24 g, mass of water molecule  *
!*     time........ t= 0.01 ps= 10^-14 s                         *
!*     charge...... 4.8 10^-10 esu= 1.60 10^-19 C                *
!*          with (1/2)M*(a/tau)**2= kT                           *
!*                                                               *
!*  Equation of motion:                                          *
!*      dv        Gamma*q'q grad r                               *
!*   m ---- = Sum --------- ------ - fgm*(2*r(i)-r(i+1)-r(i-1))  *
!*      dt           r^2      r                                  *
!*                                                               *
!*               epsLJ       sigma       sigma                   *
!*          + 48*----- grad[(-----)^12- (-----)^6] + q*E(i,j,k)  *
!*                 kT         r_ij        r_ij                   *
!*                                                               *
!*  Poisson equation:                                            *
!*   div(eps(i,j,k) [grad pot(i,j,k)]) = - 4*pi*Gamma*rho(i,j,k) *
!*   Gamma = Bjerrum/(a*kT) = e**2/(epsLJ*aLJ*kT)                *
!*                                                               *
!*****************************************************************
!*  Main program and subroutines:                                *
!*   Periodic (x,y) and bounded (z) boundaries                   *
!*                                                               *
!*   Program nanopore  MPI setup -> setups /Run_MD/ -> /moldyn/  *
!*    param_WatPa.h (parameter), PORW31_config.start3 (config)   *     
!*                                                               *
!*   /moldyn/     Time cycles, Coulomb and ES fields, L.675-     *
!*   /realteil/   Coulomb ans LJ forces, L.1665, 2270-           *
!*   /sprmul/     Spring forces, L.1670, 3350-                   * 
!*   /reflect_endpl/ Particles boundary, L.1895, 3760-           *
!*                                                               *
!*   /init/       Setups from /RUN_MD/, L.430,4430-              *
!*   /poissn_eq/  Poisson solver, L.1730,6710-                   *
!*   /escof3/     Large-scale electrostatic forces, L.6780,6840- *
!*     /bound_s/    for it > 1, L.6780,7060                      * 
!*   /cresmd/-/avmult/  Conjugate residual method, L.6790,7810-  *
!*    Graphics    /gopen/ (Adobe-2.0 postscript)                 *
!*****************************************************************
!*                                                               *
!*  To get a free format of Fortan f90 or f03, convert f77 by    *
!*    :%s/^c/!/  change 'c' of ALL column one to '!'             *
!*    :%s/^C/!/                                                  *
!*    tr 'A-Z' 'a-z' <@nanopor.f >@nanopor.f03, outside command  *
!*                                                               *
!*  For subroutines, write "use, intrinsic :: iso_c_binding",    *
!*  "implicit none" is recommended for minimizing typoerrors.    *
!*                                                               *
!*  Compilation by Linux, execution:                             *
!*  % mpif90 -mcmodel=medium -fpic -O2 -o a.out @nanoWatQa.f03 -I/opt/fftw3/include -L/opt/fftw3/lib -lfftw3 &> log       *
!*  % mpiexec -n <proc> a.out &  <proc>=6 or more                                       *
!*                                                               *
!*  First version.  2004/10/25                                   *
!*  Second version; 2006/12/18 (Complete F90)                    *
!*  Third version;  2025/07/03 (Fortran 2003, TIP5P model)       *
!*                                                               *
!*****************************************************************
!
      program nanopore
      use, intrinsic :: iso_c_binding 
!
      include    'mpif.h'
      include    'paramWatQa.h'
!
      integer(C_INT) size,rank,ierror,ipar,cl_first
      real(C_DOUBLE) wtime,walltime1,walltime2
      logical        first_recyc
!
      integer(C_INT) io_pe
      common/sub_proc/ io_pe
!
! Check results when it is necessary
      character(len=8) :: fort51
      common/fort_write/  fort51(6)
!
      fort51(1)='fortr.51'
      fort51(2)='fortr.52'
      fort51(3)='fortr.53'
      fort51(4)='fortr.54'
      fort51(5)='fortr.55'
      fort51(6)='fortr.56'
!
!  -----------------
!     num_proc = 6  !  <- size of job: paramWatQa.h
!  -----------------
!
      call mpi_init (ierror)
      call mpi_comm_rank (mpi_comm_world,rank,ierror)
      call mpi_comm_size (mpi_comm_world,size,ierror)
!
      ipar = 1 + rank             ! this is my rank number 
!
      open(50+ipar,file=fort51(ipar),status='replace',form='formatted')
      close(50+ipar)
!     +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
      io_pe= 0
      if(ipar.eq.1) io_pe= 1
!*
      if(kstart.eq.0) then        ! kstart in param.h
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
      call Run_MD (rank,size,ipar,first_recyc,cl_first) 
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
      subroutine Run_MD (rank,size,ipar,first_recyc,cl_first)
!--------------------------------------------------------------------
      use, intrinsic :: iso_c_binding 
      implicit none
!
      include    'mpif.h'
      include    'paramWatQa.h'
!
      integer(C_INT) rank,size,ipar,cl_first,ifbase,   &
                     ifrgrod,ifrodco,nframe,           &
                     np,nq,nr,nCLp,npqr,npqr3,npqr5,   &
                     io_pe,i,j
      logical        first_recyc,ende
      common/sub_proc/ io_pe
!
!  Three particle categories of np, nq, nr; npqr, npqr3, npqr5
!  Water in npqr3 is handled by xe-ze
! 
      real(C_DOUBLE),dimension(npqr50) :: xa,ya,za,ch,am        ! rotation
      real(C_DOUBLE),dimension(npqr0) ::  xg,yg,zg,amm,ag,ep    ! translation
      real(C_DOUBLE),dimension(npqr0) ::  vx,vy,vz
!
      real(C_DOUBLE),dimension(npqr0) ::  A11,A12,A13,A21,A22,A23,A31,  &
                                          A32,A33,e0,e1,e2,e3,Lgx,Lgy,Lgz
      real(C_DOUBLE),dimension(npqr50) :: xr,yr,zr
      real(C_DOUBLE),dimension(npqr0,3) :: Im 
!
      real(C_DOUBLE),dimension(np0) :: ftx,fty,ftz,fpx,fpy,fpz
      real(C_DOUBLE) xmax,ymax,zmax,xmin,ymin,zmin,Temp,  &
                     xleng,yleng,zleng,aLJ,wwat,awat
! 
!    ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      real(C_DOUBLE) doh,doL,dohL,dohcos,dohsin,doLcos,doLsin,phwat
      real(C_DOUBLE) q_O,q_H,q_L,masso,massh,epslj_w,epslj_A,epslj_B
      common/unitHL/ doh,doL,dohL,dohcos,dohsin,doLcos,doLsin,phwat
      common/unitOH/ q_O,q_H,q_L,masso,massh,epslj_w,epslj_A,epslj_B
!    ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
      integer(C_INT) ifqq
      real(C_DOUBLE) qfrac,Rpore,Hpore,Zci,Zcp,Zcn
      common/cntion/ qfrac,Rpore,Hpore,Zci,Zcp,Zcn,ifqq
!
      real(C_DOUBLE) xgc,ygc,zgc,vxg,vyg,vzg,rod_leng,Rmac,   &
                     Rhelix,q1,q2,q3,q4,angx,angy,angz,       &
                     ipx,ipy,ipz,tau_rot,gdump
      integer(C_INT) n_rodp,npqr5_0,npqr5_1,npqr5_2, &
                     npqr_0,npqr_1,npqr_2
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
      real(C_DOUBLE) pi,dt,dth,axi,Gamma,rbmax,vth,tmax,tmax7
      real(C_float) phi,tht,dtwr1,dtwr2,dtwr3
      common/parm1/ it,is
      common/parm2/ pi,dt,dth,axi,Gamma,rbmax,vth,tmax
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
!  
      real(C_DOUBLE)  t_pe,t_poisn
      common/timings/ t_pe,t_poisn
!
      real(C_DOUBLE)  Vtop,Vbot,Vtop0,Vbot0,rpr
      integer(C_INT)  fixedbc_x,fixedbc_y,itermax,filtx,filty,filtz,ipr
      common/espot/   Vtop,Vbot,Vtop0,Vbot0,fixedbc_x,fixedbc_y, &
                      itermax,filtx,filty,filtz
      common/cresm3/  ipr(10),rpr(10)
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
      real(C_DOUBLE) diel2,dielpr,daa,dbb
      common/dielec/ diel2,dielpr,daa,dbb
      common/ewald2/ nCLp     ! ip0, paramWatQa.h L.4320
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
      character(len=8) :: fort51
      common/fort_write/  fort51(6)
!
!**************************************************************
!*  np,nq,nseg :  defined in /READ_CONF/.
!
      label='+Recycle'
      call date_and_time_7 (cdate,ctime)
!
      call READ_CONF (np,nq,ifbase) 
!     +++++++++++++++++++++++++++++
!
      xp_leng= 7.
      nframe= 4
!
      call FLOPEN (nframe)
!
!--------------------------------------------------------------
!* Read basic data and constants in /READ_CONF/
!    >> Bjerrum, axi, 
!    >> dt, xmax, rcut_clf, rcutlj
!    >> kstart, cptot, dtwr1, dtwr2, dtwr3
!    >> ifqq, qfrac, Rpore, Hpore, Zcp, Zcn.
!
!   cptot (in min.)
!     cptot= 60.0  <--- 1 h
!     dtwr1=   10.
!     dtwr2=  500.
!     dtwr3=   50.
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
      rbmax= 1.5d0
      aLJ  = 1.4d0
!
      daa =   2.0d0   ! 1.5 Ang in function diel2
      dbb =   3.0d0   ! vertical depth
!
!  Water particles.
! ++++++++++++++++++++++++++++++++++++++++++++
      wwat = 18.d0         ! Water mass
      awat = 3.10d0 /2.d0  ! 3.166d0/2.d0  ! for hybrid LJ
! ++++++++++++++++++++++++++++++++++++++++++++
!
!* Bjerrum = Gamma aLJ*kT = e**2/epsLJ  
      Temp= 1.d0   ! T= 300 K
      Gamma= Bjerrum/(awat*Temp) 
!
      istop = 0    ! signal for termination: istop= 1
!**************************************************************
!
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
      do j= 1,npqr0
      A11(j)= 0
      A12(j)= 0
      A13(j)= 0
      A21(j)= 0
      A22(j)= 0
      A23(j)= 0
      A31(j)= 0
      A32(j)= 0
      A33(j)= 0
!
      e0(j)= 0
      e1(j)= 0
      e2(j)= 0
      e3(j)= 0
      end do
!
      call init_P (xg,yg,zg,vx,vy,vz,xa,ya,za,ch,am,ag,ep,amm,  &
                   A11,A12,A13,A21,A22,A23,A31,A32,A33,         & 
                   e0,e1,e2,e3,ifrgrod,ifrodco,ifbase,ipar,     & 
                   np,nq,nCLp,nr,npqr,npqr3,npqr5,              &
                   npqr5_0,npqr5_1,npqr5_2,npqr_0,npqr_1,npqr_2)
!                       np,nq are defined at /READ_CONF/
!                       nCLp, npqr, npqr5 are defined at the end of /init/
!
!  Check list
        if(kstart.eq.0) then
        open(50+ipar,file=fort51(ipar),position='append', &
                                             form='formatted')
        write(50+ipar,*) 'call init_P'
        write(50+ipar,*) 'nCLP=',nCLp

        do i= nCLp+1,nCLp+10,5 
        write(50+ipar,1903) i,xa(i),ya(i),za(i),ch(i)
        write(50+ipar,1903) i+1,xa(i+1),ya(i+1),za(i+1),ch(i+1)
        write(50+ipar,1903) i+2,xa(i+2),ya(i+2),za(i+2),ch(i+2)
 1903   format('1903 i,xa-za,ch=',i8,4f8.2)
        end do

        do j= nCLp+1,nCLp+7
        write(50+ipar,1904) j,xg(j),yg(j),zg(j),vx(j),amm(j)
 1904   format('j,xg-zg,vx,amm=',i8,3f8.2,2x,2f8.4)
        end do
!,
        do j= nCLp+1,nCLp+7 
        write(50+ipar,1908) j,A11(j),A12(j),A13(j),A21(j),A22(j),A23(j)
        write(50+ipar,1909) j,e0(j),e1(j),e2(j),e3(j)
 1908   format('1909 j,A11...=',i6,6f8.4) 
 1909   format('1909 j,e0,...=',i6,4f8.4)
        end do
        close(50+ipar)
        end if
!
      if(kstart.eq.0) then
      if(io_pe.eq.1) then
        OPEN (unit=11,file=praefixc//'.06'//suffix2,             &
              status='unknown',position='append',form='formatted')
        write(11,*) 'np,nq,nr=',np,nq,nr
        write(11,*) 'npqr.npqr3,npqr5=',npqr,npqr3,npqr5
        write(11,*) 'npqr5_2,npqr_2= ',npqr5_2,npqr_2
        write(11,*)
!
        do j= nCLp-3,nCLp+5
        write(11,1988) j,vx(j),vy(j),vz(j)
 1988   format('1988 j,vx...=',i6,1p3d12.3)
        end do
!
        do j= nCLp+1,nCLp+5 ! npqr  
        write(11,1989) j,A11(j),A12(j),A13(j),A21(j),A22(j),A23(j), &
                       A31(j),A32(j),A33(j)
 1989   format('1989 j,A11...=',i6,9f8.4)
        end do
!
        close(11)
      end if
      end if
!
!  ------------------------------------------
      call Initialisierung                 !<- in nCLp used
      call Interpol_charge_assign_function
!  ------------------------------------------
!
      if(kstart.eq.0) then
        if(io_pe.eq.1) then
          OPEN (unit=11,file=praefixc//'.06'//suffix2,             &
                status='unknown',position='append',form='formatted')
!
          write(11,*) '# Polyelectrolyte #'
          write(11,620) (ch(i),i=1,np)
!
          if(nq.ne.0) then
             write(11,*) '# cations and anions #'
             write(11,620) (ch(i),i=np+1,nCLp)
  620        format(10f5.1)
          end if
!
          close(11)
        end if
!
!------------------------------------------
!* FT12 must be mounted on all nodes
!------------------------------------------
!* Restart data: kstart= 1.
!    these data overwrite definitions made above.
!
      else
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
        OPEN (unit=12,file=praefixi//'.12'//suffix1,         &
                              status='old',form='unformatted')
!
        read(12) it,np,nq,nCLp,nr,npqr,npqr3,npqr5
        read(12) ifrgrod,ifrodco,ifbase
!                            ++++++++
        read(12) t8,xg,yg,zg,vx,vy,vz,ch,am,amm,ag,ep
        read(12) xa,ya,za
!
        read(12) A11,A12,A13,A21,A22,A23,A31,A32,A33
        read(12) e0,e1,e2,e3,xr,yr,zr,Im,Lgx,Lgy,Lgz
!
!                    + tmax is for write(12) only, is dummy here
        read(12) pi,dt,Gamma,rbmax,vth,tmax7 
        read(12) t,phi,tht,dtwr1,dtwr2,dtwr3,is
        read(12) ekin,ppot,ekn2,etot,z_pe,vzpe,vzco,vzct, &
                 vdtm,vpor,time,ecr,elj,espr
!
        read(12) iwa,iwb,iwc,iwd
        read(12) nsg,nseg
!       read(12) nbox,list
!
        read(12) np00,nnb,ist1,ist2
        read(12) diel2,daa,dbb
!                **** **** **** **** 
!       read(12) xmax,ymax,zmax,xmin,ymin,zmin
        read(12) xleng,yleng,zleng
        read(12) aLJ,epsLJ,eps_K,eps_Cl,Vtop0,Vbot0,ifLJ
        read(12) qfrac,Rpore,Hpore,Zci,Zcp,Zcn,ifqq
!
!   Rods are present
        if(ifrgrod.eq.1) then
          read(12) xgc,ygc,zgc,vxg,vyg,vzg,rod_leng,Rmac, &
                   dox,doy,doz,Rhelix,doxc,doyc,dozc,n_rodp
          read(12) q1,q2,q3,q4,angx,angy,angz,ipx,ipy,ipz
          read(12) tau_rot,gdump
        end if
!
        close(12)
!
!   These things are usually defined in the second half of /init/ !!
!   ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        epslj_w = 1.02d-14    ! 2.189d-14 ! (kjoule/mol) in erg
        epslj_A = 3.8538d-08  ! erg, A12, tip5p/e with Ewald sums
        epslj_B = 4.3676d-11  ! erg, A6
!
        q_O =  0.000d0
        q_H =  0.241d0
        q_L = -0.241d0  
!
        masso = 16.d0 /wwat   ! O
        massh =  1.d0 /wwat   ! H
!
        phwat= 104.52d0       ! tip4p  <- /moldyn/
!       awat = 3.166d0/2.d0   ! for hybrid LJ
! 
        doh  = 0.9572d0 
        doL  = 0.70d0
!
        dohcos = doh * cos(pi*phwat/(2*180.d0))
        dohsin = doh * sin(pi*phwat/(2*180.d0))
!
        doLcos = doL * cos(pi*phwat/(2*180.d0))
        doLsin = doL * sin(pi*phwat/(2*180.d0))
        dohL = dohcos +doLcos
!   ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      end if
!
! ------------------------
!     dtwr1=  10.  by /READ_CONF/
!     dtwr2= 500.
!     dtwr3=  50.
!
      if(first_recyc) then
        t8= -dt
        t = t8
!
        is= 0
        it= 0
        iwa= -1
        iwb= -1
        iwc= -1
        iwd= -1
      end if
!
!    nCLp= np +nq  !<-- defined in /init/ 
! ------------------------
      if(io_pe.eq.1) then
        OPEN (unit=11,file=praefixc//'.06'//suffix2,             &
              status='unknown',position='append',form='formatted')
!
        write(11,771) t8
        write(11,772) diel2,dielpr,daa,dbb,t00
  771   format(' t8=',f8.2)
  772   format(' diel2(pore),dielpr(water)=',2f8.2,/,'  daa,dbb,t00=',3f8.2)
!
        write(11,773) np,nq,nCLp,npqr,npqr3,npqr5
  773   format(' np,nq,nCLp=',3i6,'  npqr,npqr3,npqr5=',3i8)
!
        close(11)
      end if
!
!       call mpi_broadcast ()
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
!*****************************************************************
!*  Step 3: Subroutine /moldyn/                                  *
!*****************************************************************
!*----------------------------------------------------------------
      call moldyn (xa,ya,za,ch,am,xg,yg,zg,amm,ag,ep,     &
                   vx,vy,vz,                              &
                   A11,A12,A13,A21,A22,A23,A31,A32,A33,   &
                   e0,e1,e2,e3,Lgx,Lgy,Lgz,xr,yr,zr,Im,   &
                   first_recyc,ifbase,ifrgrod,ifrodco,    &
                   np,nq,nCLp,nr,npqr,npqr3,npqr5,        &
                   ipar,size,cl_first)
!*----------------------------------------------------------------
!
!*****************************************************************
!*  Step 4: Restart data.                                        *
!*****************************************************************
!
      if(io_pe.eq.1) then
!
        OPEN (unit=12,file=praefixc//'.12'//suffix2,     &
                      status='replace',form='unformatted')
!
        write(12) it,np,nq,nCLp,nr,npqr,npqr3,npqr5
        write(12) ifrgrod,ifrodco,ifbase
!
        write(12) t8,xg,yg,zg,vx,vy,vz,ch,am,amm,ag,ep
        write(12) xa,ya,za
!
        write(12) A11,A12,A13,A21,A22,A23,A31,A32,A33
        write(12) e0,e1,e2,e3,xr,yr,zr,Im,Lgx,Lgy,Lgz
!
        write(12) pi,dt,Gamma,rbmax,vth,tmax  !<- L.2160 IF
        write(12) t,phi,tht,dtwr1,dtwr2,dtwr3,is
        write(12) ekin,ppot,ekn2,etot,z_pe,vzpe,vzco,vzct, &
                  vdtm,vpor,time,ecr,elj,espr
!
        write(12) iwa,iwb,iwc,iwd
        write(12) nsg,nseg
!       write(12) nbox,list
!
        write(12) np00,nnb,ist1,ist2
        write(12) diel2,daa,dbb
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
      include    'paramWatQa.h'
!
      character   praefix_name*6
      common /confdatas/ praefix_name
!
      integer(C_INT) io_pe
      common/sub_proc/ io_pe
!
!   --------------------------
      if(io_pe.eq.1) then
!   --------------------------
!
        OPEN (unit=11,file=praefixc//'.06'//suffix2,             &
              status='unknown',position='append',form='formatted')
!
        write(11,*) ' Uses ',praefixs//'_config.start'//suffix0
        write(11,*) ' Uses ',praefix_name
        close(11)
!
        OPEN (unit=77,file=praefixc//'.77'//suffix2//'.ps', &
                                            form='formatted')
        call gopen (nframe)
        close (77)
!
      end if
!
!     OPEN (unit=12,file=praefixc//'.12'//suffix2,
!     OPEN (unit=13,file=praefixc//'.13'//suffix2,
!
      return
      end subroutine FLOPEN
!
!
!-----------------------------------------------------------------
      subroutine moldyn (xa,ya,za,ch,am,xg,yg,zg,amm,ag,ep,    &
                         vx,vy,vz,                             &
                         A11,A12,A13,A21,A22,A23,A31,A32,A33,  &
                         e0,e1,e2,e3,Lgx,Lgy,Lgz,xr,yr,zr,Im,  &
                         first_recyc,ifbase,ifrgrod,ifrodco,   &
                         np,nq,nCLp,nr,npqr,npqr3,npqr5,       &
                         ipar,size,cl_first)
!-----------------------------------------------------------------
!*  Double precision.
!
      use, intrinsic :: iso_c_binding 
      implicit none
!
      include    'mpif.h'
      include    'paramWatQa.h'
!
      real(C_DOUBLE),dimension(npqr50) :: xa,ya,za,ch,am,fcx,fcy,fcz
      real(C_DOUBLE),dimension(npqr30) :: xe,ye,ze,chg,fgx,fgy,fgz
!
      real(C_DOUBLE),dimension(npqr0) ::  xg,yg,zg,amm,ag,ep,vx,vy,vz
      real(C_DOUBLE),dimension(npqr0) ::  omg_x,omg_y,omg_z
!
      real(C_DOUBLE),dimension(npqr0) ::  A11,A12,A13,A21,A22,A23,A31,A32, &
                                          A33,e0,e1,e2,e3,Lgx,Lgy,Lgz
      real(C_DOUBLE),dimension(npqr50) :: xr,yr,zr
      real(C_DOUBLE),dimension(npqr0,3) :: Im 
!
      real(C_DOUBLE),dimension(np0) :: ftx,fty,ftz,fpx,fpy,fpz
!
      real(C_DOUBLE) xg1,yg1,zg1,xr1,yr1,zr1,xr2,yr2,zr2,        &
                     xr3,yr3,zr3,xr4,yr4,zr4,xr5,yr5,zr5,        &
                     Torqx1,Torqy1,Torqz1,LLgx1,LLgy1,LLgz1,     &
                     omg_x1,omg_y1,omg_z1,pe01,pe11,pe21,pe31,   &
                     xxa,yya,zza,xxb,yyb,zzb,xhh,yhh,zhh,vec1,   &
                     xpoint,ypoint,zpoint,xxc,yyc,zzc,vec3,      &
                     corr,fcxx,fcyy,fczz
!
      real(C_DOUBLE) doh,doL,dohL,dohcos,dohsin,doLcos,doLsin,phwat
      common/unitHL/ doh,doL,dohL,dohcos,dohsin,doLcos,doLsin,phwat
      integer(C_INT) ncorr
!
      integer(C_INT) ipar,size,np,nq,nCLp,nr,npqr,npqr3,npqr5,   &
                     io_pe,cl_first
      common/sub_proc/ io_pe
!
!  -----------------------------------------------------------------
      integer(C_INT) i0,i1,i2,i3,i4,cnt_recv3,disp_recv3,   &
                     cnt_recv5,disp_recv5,disp_recv6
      common/dat_typ1/ i0(30),i1(30),i2(30)
      common/dat_typ3/ i3(30),i4(30),cnt_recv3(30),disp_recv3(30)
      common/dat_typ5/ cnt_recv5(30),disp_recv5(30),disp_recv6(30)
!  -----------------------------------------------------------------
!
      real(C_DOUBLE) t8,cptot
      common/headr3/ t8,cptot
!
!  The space cell index (0,mx-1)
      real*8   t1,t2,t3,d1,xx,yy,zz
      real(C_DOUBLE),dimension(0:mx-1,0:my-1,0:mz-1) :: &
                                  rho,ex,ey,ez,dec2,pot,ez1
      real(C_DOUBLE),dimension(0:ip0**3-1,0:npqr30-1) :: ql
      integer(C_INT),dimension(0:npqr30-1,0:2) :: gg
!
      real(C_DOUBLE) diel
      common/potsav/ pot,dec2
      logical     first_recyc 
!
      integer(C_INT) ierr
      real(C_DOUBLE) sres,avex,rsdl,conv
      common/crconv/ sres,avex,rsdl,conv,ierr
!
      real(C_float),dimension(npqr50) :: x4,y4,z4,vx4,vy4,vz4,ch4
      real(C_float),dimension(npqr0) :: am4,ag4
      integer(C_INT) i,j,k,l,ll,m,m0,xpos,ypos,zpos, &
                     mx1,mx2,my1,my2,mz1,mz2
!
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
                     ndim,iterp,ncti,ncoi,ntimes,istep,ii
!
      integer(C_INT) np00,nnb,ist1,ist2,nc_DNA
      common/pbase/ np00,nnb,ist1(100),ist2(100)
      real(C_DOUBLE) diel2,dielpr,daa,dbb
      common/dielec/ diel2,dielpr,daa,dbb
!
      real(C_float) phi,tht,dtwr1,dtwr2,dtwr3,                &
                    qfrac4,Rpore4,Hpore4,Zcp4,Zcn4,Bjerrum4,  &
                    xmax4,ymax4,zmax4,walltime,               &
                    dtm,vm,s0,s1,s2,vsq,ekin0,ekin1,          &
                    ekin2,sz1,svz1,svz2,svz3,ef2,ef20,dx,dy
      real(C_DOUBLE) pi,dt,dth,axi0,Gamma,rbmax,vth,tmax,       &
                    xmax,ymax,zmax,xmin,ymin,zmin,symp,symp2,   &
                    xleng,yleng,zleng,dgaus2,vmax0,vmax1,vmax2, &
                    e_c_s,e_grid,e_c_r,e_lj,e_elas,q0,ss,dt0,   &
                    fv,vv,dv,Temp
      common/gaus1/ fv(101),vv,dv
!
      common/parm1/ it,is
      common/parm2/ pi,dt,dth,axi0,Gamma,rbmax,vth,tmax
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
      common/xtabl2/ ptx(-199:3000),pty(-199:3000),ptz(-100:4000)
!
      integer(C_INT) pxr,pxc,pxl,pyr,pyc,pyl,pzr,pzc,pzl
      common/ptable/ pxr(-10:mx+10),pxc(-10:mx+10),pxl(-10:mx+10), &
                     pyr(-10:my+10),pyc(-10:my+10),pyl(-10:my+10), &
                     pzr(-10:mz+10),pzc(-10:mz+10),pzl(-10:mz+10)
!
      real(C_DOUBLE) Vtop,Vbot,Vtop0,Vbot0,rpr
      integer(C_INT) fixedbc_x,fixedbc_y,itermax,        &
                     filtx,filty,filtz,ipr
      common/espot/  Vtop,Vbot,Vtop0,Vbot0,fixedbc_x,fixedbc_y, &
                     itermax,filtx,filty,filtz
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
      real(C_DOUBLE) t_recyc  !<- in /READ_CONF/
      common/trecyc/ t_recyc
      real(C_DOUBLE) wtime0,wtime1,wtime2,wtime3,wtime4,wtime5
!
      integer(C_INT) ifLJ
      real(C_float)  epsLJ4
      real(C_DOUBLE) epsLJ,eps_K,eps_Cl
      common/parmlj/ epsLJ,eps_K,eps_Cl
      common/parmlj2/ ifLJ
!
      logical    first_11/.true./,first_15/.true./, &
                 first_ppl/.true./,if_xyz3/.true./
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
!++++++++++++++++++++++++++++++++++
!* 1) Particles by i1(ipar)       +
!++++++++++++++++++++++++++++++++++
!   npqr = np+nq+nr/5    <- /init/ translation coordinate in oxygens 
!   npqr5= np+nq+nr      <- /init/ 5-point calculation in /realteil/
!
      do k= 1,8        
      if(num_proc.eq.1) then
        i0(k)= 1 
        i1(k)= nCLp+1
        i2(k)= npqr
      end if 
!
      if(num_proc.eq.6) then
        if(k.eq.1) then
          i0(k)= 1         ! i0(1)=1 for nCLp 
          i1(k)= nCLp+1    ! i1(1) starts here
          i2(k)= npqr/6 
        else if(k.eq.2) then
          i1(k)= npqr/6+1
          i2(k)= 2*npqr/6
        else if(k.eq.3) then
          i1(k)= 2*npqr/6+1
          i2(k)= 3*npqr/6 
        else if(k.eq.4) then
          i1(k)= 3*npqr/6+1
          i2(k)= 4*npqr/6
        else if(k.eq.5) then
          i1(k)= 4*npqr/6+1
          i2(k)= 5*npqr/6 
        else if(k.eq.6) then
          i1(k)= 5*npqr/6+1
          i2(k)= npqr
        end if
      end if
!
      if(num_proc.eq.8) then
        if(k.eq.1) then
          i0(k)= 1
          i1(k)= nCLp
          i2(k)= npqr/8
        else if(k.eq.2) then
          i1(k)= npqr/8+1
          i2(k)= 2*npqr/8
        else if(k.eq.3) then
          i1(k)= 2*npqr/8+1
          i2(k)= 3*npqr/8
        else if(k.eq.4) then
          i1(k)= 3*npqr/8+1
          i2(k)= 4*npqr/8
        else if(k.eq.5) then
          i1(k)= 4*npqr/8+1
          i2(k)= 5*npqr/8
        else if(k.eq.6) then
          i1(k)= 5*npqr/8+1
          i2(k)= 6*npqr/8
        else if(k.eq.7) then
          i1(k)= 6*npqr/8+1
          i2(k)= 7*npqr/8
        else if(k.eq.6) then
          i1(k)= 7*npqr/8+1
          i2(k)= npqr
        end if
      end if
!
!     disp_recv1(k)= i1(k) -1          ! i1(k)=1 if no displacement 
!     cnt_recv1(k) = i2(k) -i1(k) +1   ! i2(k), i1(k)
!
      if(io_pe.eq.1) then
        OPEN (unit=11,file=praefixc//'.06'//suffix2,             &
              status='unknown',position='append',form='formatted')
!
        write(11,115) k,i1(k),i2(k)
  115   format('# k=',i7,'  i1(k),i2(k),=',2i8)
!
        if(k.eq.8) write(11,*)
        close(11)
      end if
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
      disp_recv3(k)= i3(k)
      cnt_recv3(k) = i4(k) -i3(k) +1
!
      if(io_pe.eq.1) then
        OPEN (unit=11,file=praefixc//'.06'//suffix2,             &
              status='unknown',position='append',form='formatted')
!
        write(11,117) k,i3(k),i4(k)
  117   format('# k=',i7,'  i3(k),i4(k)=',2i8)
!
        if(k.eq.8) write(11,*)
        close(11)
      end if
      end do
!
!--------------------------
!*  Settings
!--------------------------
!
!* For ions.
!     ww1  = 38.d0 /wwat   ! K,  2.27 Ang, K+ 2.98 Ang
!     ww2  = 34.d0 /wwat   ! Cl, 1.75 Ang, Cl- 1.81 Ang
!
      vth = sqrt(epsLJ)
!     vmax0= vth/sqrt(216.d0/wwat)   !<-- mass of sugar,am(2)
      vmax0= vth/sqrt( 96.d0/wwat)   !<-- mass of PO_4, am(1)/wwat
      vmax1= vth/sqrt( 38.d0/wwat)   !<-- mass of K(+)
      vmax2= vth/sqrt( 18.d0/wwat)   !<-- water
!
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
        istep= 0  !<-- do for the Poisson solver
!
!* Solvent - water
        do j= 1,np
        vx(j)= dgaus2(vmax0)
        vy(j)= dgaus2(vmax0)
        vz(j)= dgaus2(vmax0)
        end do
!
        do j= np+1,nCLp
        vx(j)= dgaus2(vmax1)
        vy(j)= dgaus2(vmax1)
        vz(j)= dgaus2(vmax1)
        end do
!
        do j= nCLp+1,npqr
        vx(j)= dgaus2(vmax2)
        vy(j)= dgaus2(vmax2)
        vz(j)= dgaus2(vmax2)
        end do 
!
!  Initial value e0,e1,e2,e3 in /init/
        do j= nCLp+1,npqr
        Lgx(j)= 0
        Lgy(j)= 0
        Lgz(j)= 0
        end do
!
        do i= 1,npqr3
        fgx(i)= 0.d0
        fgy(i)= 0.d0
        fgz(i)= 0.d0
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
        do j= 1,nCLp
        ch4(j)= ch(j)
        am4(j)= am(j)
        ag4(j)= ag(j)
        end do
!
        ii= nCLp
        do i= nCLp+1,npqio
        if(mod(i-nCLp,5).ge.1 .and. mod(i-nCLp,5).le.3) then 
          ii= ii +1
!
          if(mod(i-nCLp,5).eq.1) ch4(ii)= -1.0d0
          if(mod(i-nCLp,5).eq.2) ch4(ii)=  0.5d0
          if(mod(i-nCLp,5).eq.3) ch4(ii)=  0.5d0
        end if
        end do
!
        j= nCLp
        do i= nCLp+1,npqr5,5 
        j= j +1
!
        am4(j)= am(i)+am(i+1)+am(i+2)
        ag4(j)= ag(j) 
        end do
!
        write(13) ch4,am4,ag4
        close(13)
      end if
!
!***********************************************
!* >> Reset ion velocities at t= t_vreset      *
!***********************************************
!*  Read data in /READ_CONF/
!   Values are suggested:
!     t_pe    =    3. ! PE motion starts at t_pe
!     t_poisn =    1. ! solve Poisson eq after t_poisn
!
!****************************
!* Step 1A: Initialization  *
!****************************
 1000 continue
!
!    ++++++++++++
      t8= t8 +dt
      t = t8
      it= it +1
!    ++++++++++++
!
      cl_first= 2
      call clocks (wtime0,size,cl_first)
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
        do i= 0,mx-1
        do j= 0,my-1
        do k= 0,mz-1
        ez1(i,j,k)= 0
        end do 
        end do 
        end do 
!
!  Present values of the top and bottom potential
      else if(t8.gt.t00) then
        Vtop = Vtop0 *(1.d0 -exp(-(t8-t00)/100.d0))
        Vbot = Vbot0 *(1.d0 -exp(-(t8-t00)/100.d0))
!
        do i= 0,mx-1
        do j= 0,my-1
        do k= 0,mz-1
        ez1(i,j,k)= - (Vtop-Vbot)/zleng   ! ez1= -grad pot
!                                              = -(Vtop-Vbot)/zleng 
        end do  
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
!* Count the number of ions in cells
!   Results are contained in common/ionpop/
!
!!      call cell_ions (xa,ya,za,vz,ch,np,nq,nCLp,is) 
      end if
!
      Temp= 1.d0
      Gamma= Bjerrum/(awat*Temp)
!
!***********************************************************
!*  Main Loop                                              *
!***********************************************************
!******************************
!*  Step 2A: Velocity update  *
!******************************
!
      e_lj   = 0.d0
      e_c_r  = 0.d0
!*
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
      do i = 1,npqr5
      fcx(i)= 0.d0
      fcy(i)= 0.d0
      fcz(i)= 0.d0
      end do
!
!  Separation of R and s_k {j=1; i=1,2,3}
!   (xr,yr,zr) = A*(xr,yr,zr)
!
!  Three sites O-H-H are the base water
!  Only the simulation starts to define Im(j,1-3) and save them
!
        cl_first= 2
        call clocks (wtime1,size,cl_first)
!
      if(it.eq.1) then 
!
        j= nCLp 
        do i= nCLp+1,npqr5,5
        j= j +1
!
        xg1= (16*xa(i) +xa(i+1) +xa(i+2))/18.d0
        yg1= (16*ya(i) +ya(i+1) +ya(i+2))/18.d0
        zg1= (16*za(i) +za(i+1) +za(i+2))/18.d0
!
!  A_ij is read from /init/ at the first step.
        xr1= xa(i) -xg1
        yr1= ya(i) -yg1
        zr1= za(i) -zg1
        xr(i)= A11(j)*xr1 +A12(j)*yr1 +A13(j)*zr1
        yr(i)= A21(j)*xr1 +A22(j)*yr1 +A23(j)*zr1
        zr(i)= A31(j)*xr1 +A32(j)*yr1 +A33(j)*zr1
!
        xr2= xa(i+1) -xg1
        yr2= ya(i+1) -yg1
        zr2= za(i+1) -zg1
        xr(i+1)= A11(j)*xr2 +A12(j)*yr2 +A13(j)*zr2
        yr(i+1)= A21(j)*xr2 +A22(j)*yr2 +A23(j)*zr2
        zr(i+1)= A31(j)*xr2 +A32(j)*yr2 +A33(j)*zr2
!
        xr3= xa(i+2) -xg1
        yr3= ya(i+2) -yg1
        zr3= za(i+2) -zg1
        xr(i+2)= A11(j)*xr3 +A12(j)*yr3 +A13(j)*zr3
        yr(i+2)= A21(j)*xr3 +A22(j)*yr3 +A23(j)*zr3
        zr(i+2)= A31(j)*xr3 +A32(j)*yr3 +A33(j)*zr3
!
!  Principal coordinate  L_P= R * L_R
!   Rotation matrix A_ij^n by (6.48); the base step 
!
        Im(j,1)=  am(i)*  (yr(i  )**2 +zr(i  )**2)  &
                 +am(i+1)*(yr(i+1)**2 +zr(i+1)**2)  &
                 +am(i+2)*(yr(i+2)**2 +zr(i+2)**2)   
!
        Im(j,2)=  am(i)*  (zr(i  )**2 +xr(i  )**2)  &
                 +am(i+1)*(zr(i+1)**2 +xr(i+1)**2)  &
                 +am(i+2)*(zr(i+2)**2 +xr(i+2)**2)  
!
        Im(j,3)=  am(i)*  (xr(i  )**2 +yr(i  )**2)  &
                 +am(i+1)*(xr(i+1)**2 +yr(i+1)**2)  &
                 +am(i+2)*(xr(i+2)**2 +yr(i+2)**2) 
        end do
      end if
!
!****************************************************************
!  All information in the next job is connected from xr-zr,     *
!  Im(), e0-e3, Lgx-Lgz, A11..A33 in the read(12)               *
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!  Translation: v(n-1/2) -> v(n+1/2), and x(n) -> x(n+1)
!  Electric field for water and ions, i >= nq+1 (fcx(i) +ch(i)*exc)
!
      j= nCLp
      do i= nCLp+1,npqr5,5   ! for i -> nCLp+1,npqr5,5
      j= j +1                ! for j
!
      xg1= (16*xa(i) +xa(i+1) +xa(i+2))/18.d0
      yg1= (16*ya(i) +ya(i+1) +ya(i+2))/18.d0
      zg1= (16*za(i) +za(i+1) +za(i+2))/18.d0
!
      xr1= xa(i) -xg1 
      yr1= ya(i) -yg1
      zr1= za(i) -zg1
!
      xr2= xa(i+1) -xg1
      yr2= ya(i+1) -yg1
      zr2= za(i+1) -zg1
!
      xr3= xa(i+2) -xg1
      yr3= ya(i+2) -yg1
      zr3= za(i+2) -zg1
!
      xr4= xa(i+3) -xg1
      yr4= ya(i+3) -yg1
      zr4= za(i+3) -zg1
!
      xr5= xa(i+4) -xg1
      yr5= ya(i+4) -yg1
      zr5= za(i+4) -zg1
!
      Torqx1 = & 
               (yr1*fcz(i)   -zr1*fcy(i)     & ! O
               +yr2*fcz(i+1) -zr2*fcy(i+1)   & ! H1
               +yr3*fcz(i+2) -zr3*fcy(i+2)   & ! H2
               +yr4*fcz(i+3) -zr4*fcy(i+3)   & ! L1
               +yr5*fcz(i+4) -zr5*fcy(i+4))    ! L2
!
      Torqy1 = & 
               (zr1*fcx(i)   -xr1*fcz(i)     & ! O
               +zr2*fcx(i+1) -xr2*fcz(i+1)   & ! H1
               +zr3*fcx(i+2) -xr3*fcz(i+2)   & ! H2
               +zr4*fcx(i+3) -xr4*fcz(i+3)   & ! L1
               +zr5*fcx(i+4) -xr5*fcz(i+4))    ! L2
!
      Torqz1 = & 
               (xr1*fcy(i)   -yr1*fcx(i)     & ! O
               +xr2*fcy(i+1) -yr2*fcx(i+1)   & ! H1
               +xr3*fcy(i+2) -yr3*fcx(i+2)   & ! H2
               +xr4*fcy(i+3) -yr4*fcx(i+3)   & ! L1
               +xr5*fcy(i+4) -yr5*fcx(i+4))    ! L2
!
!  Trial move on a half time steps 
!
      LLgx1= Lgx(j) +Torqx1*dth  ! Lgx(n-1/2),Torqx(n) -> LLgx(n)
      LLgy1= Lgy(j) +Torqy1*dth
      LLgz1= Lgz(j) +Torqz1*dth
!
      omg_x1= (A11(j)*LLgx1 +A12(j)*LLgy1 +A13(j)*LLgz1)/Im(j,1)
      omg_y1= (A21(j)*LLgx1 +A22(j)*LLgy1 +A23(j)*LLgz1)/Im(j,2)
      omg_z1= (A31(j)*LLgx1 +A32(j)*LLgy1 +A33(j)*LLgz1)/Im(j,3)
!
!  Prediction of e0(n+1/2): q(n+1/2)= q(n) +dth*Q(n)*omg(n)
!  Rotation matrix A_ij^(n+1/2) by e0(j),e1(j),...
!   in each pe01, mixes omg_x1, omg_y1, omg_z1
!
      pe01= e0(j) +(dth/2.d0)*(  &
                 -e1(j)*omg_x1 -e2(j)*omg_y1 -e3(j)*omg_z1 )
      pe11= e1(j) +(dth/2.d0)*(  &
                  e0(j)*omg_x1 -e3(j)*omg_y1 +e2(j)*omg_z1 )
      pe21= e2(j) +(dth/2.d0)*(  &
                  e3(j)*omg_x1 +e0(j)*omg_y1 -e1(j)*omg_z1 )
      pe31= e3(j) +(dth/2.d0)*(  &
                 -e2(j)*omg_x1 +e1(j)*omg_y1 +e0(j)*omg_z1 )
!
      A11(j)= pe01**2 +pe11**2 -pe21**2 -pe31**2 
      A12(j)= 2*(pe11*pe21 +pe01*pe31) 
      A13(j)= 2*(pe11*pe31 -pe01*pe21)
      A21(j)= 2*(pe11*pe21 -pe01*pe31) 
      A22(j)= pe01**2 -pe11**2 +pe21**2 -pe31**2
      A23(j)= 2*(pe21*pe31 +pe01*pe11)
      A31(j)= 2*(pe11*pe31 +pe01*pe21)
      A32(j)= 2*(pe21*pe31 -pe01*pe11)
      A33(j)= pe01**2 -pe11**2 -pe21**2 +pe31**2
!
!  All information of the next job is connected by xr-zr, Im(), 
!  The full time step
!
      Lgx(j)= Lgx(j) +Torqx1*dt  ! Lgx(n-1/2),Torq(n) -> Lgx(n+1/2)
      Lgy(j)= Lgy(j) +Torqy1*dt
      Lgz(j)= Lgz(j) +Torqz1*dt
!
!  omgx(n+1/2)
      omg_x1= (A11(j)*Lgx(j) +A12(j)*Lgy(j) +A13(j)*Lgz(j))/Im(j,1)
      omg_y1= (A21(j)*Lgx(j) +A22(j)*Lgy(j) +A23(j)*Lgz(j))/Im(j,2)
      omg_z1= (A31(j)*Lgx(j) +A32(j)*Lgy(j) +A33(j)*Lgz(j))/Im(j,3)
!
      omg_x(j)= omg_x1
      omg_y(j)= omg_y1
      omg_z(j)= omg_z1
!
!  Then, e0(n) -> e0(n+1)
      e0(j)= e0(j) +(dt/2.d0)*(  &
                -pe11*omg_x1 -pe21*omg_y1 -pe31*omg_z1 )
      e1(j)= e1(j) +(dt/2.d0)*(  &  
                 pe01*omg_x1 -pe31*omg_y1 +pe21*omg_z1 )
      e2(j)= e2(j) +(dt/2.d0)*(  &
                 pe31*omg_x1 +pe01*omg_y1 -pe11*omg_z1 )
      e3(j)= e3(j) +(dt/2.d0)*(  &
                -pe21*omg_x1 +pe11*omg_y1 +pe01*omg_z1 )
!
      A11(j)= e0(j)**2 +e1(j)**2 -e2(j)**2 -e3(j)**2 
      A12(j)= 2*(e1(j)*e2(j) +e0(j)*e3(j)) 
      A13(j)= 2*(e1(j)*e3(j) -e0(j)*e2(j))
      A21(j)= 2*(e1(j)*e2(j) -e0(j)*e3(j)) 
      A22(j)= e0(j)**2 -e1(j)**2 +e2(j)**2 -e3(j)**2
      A23(j)= 2*(e2(j)*e3(j) +e0(j)*e1(j))
      A31(j)= 2*(e1(j)*e3(j) +e0(j)*e2(j))
      A32(j)= 2*(e2(j)*e3(j) -e0(j)*e1(j))
      A33(j)= e0(j)**2 -e1(j)**2 -e2(j)**2 +e3(j)**2
!
!  All information of the next job is connected by xr-zr, Im(), 
!
!  Vector (xa,ya,za)= R_G(j) +R^-1*(xr,yr,zr) on A(n+1)
!   b= R^-1= (A11,A21,A31,...)  There are 5-point molecules, 
!   with zero-mass at the sites
!
      dohL = dohcos +doLcos
!
      xa(i  )= xg(j) +A11(j)*xr(i) +A21(j)*yr(i) +A31(j)*zr(i)
      ya(i  )= yg(j) +A12(j)*xr(i) +A22(j)*yr(i) +A32(j)*zr(i)
      za(i  )= zg(j) +A13(j)*xr(i) +A23(j)*yr(i) +A33(j)*zr(i)
!
      xa(i+1)= xg(j) +A11(j)*xr(i+1) +A21(j)*yr(i+1) +A31(j)*zr(i+1)
      ya(i+1)= yg(j) +A12(j)*xr(i+1) +A22(j)*yr(i+1) +A32(j)*zr(i+1)
      za(i+1)= zg(j) +A13(j)*xr(i+1) +A23(j)*yr(i+1) +A33(j)*zr(i+1)
! 
      xa(i+2)= xg(j) +A11(j)*xr(i+2) +A21(j)*yr(i+2) +A31(j)*zr(i+2)
      ya(i+2)= yg(j) +A12(j)*xr(i+2) +A22(j)*yr(i+2) +A32(j)*zr(i+2)
      za(i+2)= zg(j) +A13(j)*xr(i+2) +A23(j)*yr(i+2) +A33(j)*zr(i+2)
!  
!* Corrections
!  1) xr1= xa(i) -xg1
!  2) xr(i)= A11(j)*xr1 +A12(j)*yr1 +A13(j)*zr1
!  3) xa(i)= xg(j) +A11(j)*xr(i) +A21(j)*yr(i) +A31(j)*zr(i)
!  Outside the triangle plane on virtual sites
!
      xxa= xa(i+2) -xa(i+1)
      yya= ya(i+2) -ya(i+1)
      zza= za(i+2) -za(i+1)
!
      xxb= xa(i) -(xa(i+2)+xa(i+1))/2
      yyb= ya(i) -(ya(i+2)+ya(i+1))/2
      zzb= za(i) -(za(i+2)+za(i+1))/2
      vec1= sqrt(xxb*xxb +yyb*yyb +zzb*zzb)
!
      xhh= (xa(i+2)+xa(i+1))/2
      yhh= (ya(i+2)+ya(i+1))/2
      zhh= (za(i+2)+za(i+1))/2
!
      xpoint= xhh +dohL*xxb/vec1
      ypoint= yhh +dohL*yyb/vec1
      zpoint= zhh +dohL*zzb/vec1
!
      xxc= yya*zzb -zza*yyb
      yyc= zza*xxb -xxa*zzb
      zzc= xxa*yyb -yya*xxb
      vec3= sqrt(xxc*xxc +yyc*yyc +zzc*zzc)
!
      xa(i+3)= xpoint +doLsin*xxc/vec3 
      ya(i+3)= ypoint +doLsin*yyc/vec3
      za(i+3)= zpoint +doLsin*zzc/vec3
!
      xa(i+4)= xpoint -doLsin*xxc/vec3 
      ya(i+4)= ypoint -doLsin*yyc/vec3
      za(i+4)= zpoint -doLsin*zzc/vec3
      end do
!
!* Correction in every 10 time steps
!
      ncorr= 10 
      if(mod(it,ncorr).eq.0) then
!                        +++  
        do j= nCLp+1,npqr 
        corr= 1.d0/sqrt(e0(j)**2 +e1(j)**2 +e2(j)**2 +e3(j)**2)

        e0(j)= corr*e0(j)
        e1(j)= corr*e1(j)
        e2(j)= corr*e2(j)
        e3(j)= corr*e3(j)
        end do
      end if
!
!  All information of the next job is connected by xr-zr, Im(), 
!* End of the long i-loop
!
        cl_first= 2
        call clocks (wtime2,size,cl_first)
!
!--------------------
!*  Find forces
!--------------------
! 
      call realteil (xa,ya,za,ch,xg,yg,zg,ag,ep,fcx,fcy,fcz, & 
                     ipar,size,np,nq,nCLp,npqr,npqr5)
!
      call sprmul (xg,yg,zg,ag,fpx,fpy,fpz,np)
!
        cl_first= 2
        call clocks (wtime3,size,cl_first)
!
!**************************************************
!*  Step 3A: Coulomb forces on grids              *
!**************************************************
!   The space cell follows i= 0,mx-1 in /charges/.
!
!*  -4*pi*rho = div*[eps*grad pot]
!    pot = pot_q + pot_g
!        pot_q: direct sum of coulomb forces
!        pot_g: geometry correction (boundary)
!
! ------------------------------------
      ntimes= 10      !<-- Take average in every ntimes steps 
      if(t8.lt.t_poisn) go to 370 
! ------------------------------------
      if(mod(it,ntimes).eq.0) then 
!**   
        istep= istep +1  
!
        do i= 1,npqr3 
        fgx(i)= 0.d0
        fgy(i)= 0.d0
        fgz(i)= 0.d0
        end do
!                     *** 
        call charges (rho,xg,yg,zg,xa,ya,za,xe,ye,ze,ch,chg, &
                      A11,A12,A13,A21,A22,A23,A31,A32,A33,   &
                      xr,yr,zr,ql,gg,ipar,np,nq,nCLp,npqr,   &
                      npqr3,npqr5)
!
        do k= 0,mz-1
        do j= 0,my-1
        do i= 0,mx-1
        rho(i,j,k)= 0
        pot(i,j,k)= 0
!
        ex(i,j,k)= 0
        ey(i,j,k)= 0
        ez(i,j,k)= 0
        end do
        end do
        end do
!
        do k= 1,mz-2 
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
        call poissn_eq (rho,pot,ndim,itermax,iterp,ipar) 
!       ++++++++++++++++++++++++++++++++++++++++++++++++
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
        if(mod(istep,50).eq.1 .and. io_pe.eq.1) then   !  istep=50 for Dt=5 
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
        do i= 1,npqr3
        m0= -1
!
        do j = -1,1   ! symmetric
        do k = -1,1
        do m = -1,1
        xpos = pxc(gg(i,0) + j)   ! use fresh g()
        ypos = pyc(gg(i,1) + k)
        zpos = pzc(gg(i,2) + m)
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
!*
        end do
!
!* Coulomb forces are calculated under a large stride.
!
        if(t8.lt.t_pe) then
          do i= 1,np
          fgx(i)= 0.d0
          fgy(i)= 0.d0
          fgz(i)= 0.d0
          end do
        end if
!
      end if
  370 continue 
!
        cl_first= 2
        call clocks (wtime4,size,cl_first)
!
!* End of the long loop
!*  Do not fold back positions: (-l/2, l/2). 
!
!
      if(t8.ge.t_pe) then
        j= 0
!
        do i= 1,np
        j= j +1
        dtm= dt/amm(j)
!
!    ez1(i,j,k)= -grad pot= - (Vtop-Vbot)/zleng
!    charges -> poissn -> fgx-fgz, ql(m0,i)
!
        vx(j)= vx(j) +(fcx(i) +fpx(i) +fgx(i))*dtm 
        vy(j)= vy(j) +(fcy(i) +fpy(j) +fgy(i))*dtm
        vz(j)= vz(j) +(fcz(i) +fpz(j) +fgz(i))*dtm
!
        xg(j)= xg(j) +dt*vx(j) 
        yg(j)= yg(j) +dt*vy(j)
        zg(j)= zg(j) +dt*vz(j)
        end do
!
      else
        do i= 1,np
        vx(j)= 0
        vy(j)= 0
        vz(j)= 0
        end do  
      end if
! 
!     if(np.ne.0) then
!        nc_DNA= np00/2
!
!        vx(nc_DNA)= 0.d0
!        vy(nc_DNA)= 0.d0
!        vz(nc_DNA)= 0.d0
!     end if
!
!
      j= np
      do i= np+1,nCLp
      j= j +1
      dtm= dt/amm(j)
!
      vx(j)= vx(j) +(fcx(i) +fgx(i))*dtm
      vy(j)= vy(j) +(fcy(i) +fgy(i))*dtm
      vz(j)= vz(j) +(fcz(i) +fgz(i))*dtm
!
      xg(j)= xg(j) +dt*vx(j) 
      yg(j)= yg(j) +dt*vy(j)
      zg(j)= zg(j) +dt*vz(j)
      end do
!*
!
      j= nCLp
      k= nCLp-2
!
      do i= nCLp+1,npqr5,5
      j= j +1
      k= k +3
      dtm= dt/amm(j)
!                      Short-Coulomb+LJ  fgx:Long-range
      vx(j)= vx(j) +( fcx(i)+fcx(i+1)+fcx(i+2)+fcx(i+3)+fcx(i+4) & 
                     +fgx(k)+fgx(k+1)+fgx(k+2) )*dtm 
      vy(j)= vy(j) +( fcy(i)+fcy(i+1)+fcy(i+2)+fcy(i+3)+fcy(i+4) &
                     +fgy(k)+fgy(k+1)+fgy(k+2) )*dtm
      vz(j)= vz(j) +( fcz(i)+fcz(i+1)+fcz(i+2)+fcz(i+3)+fcz(i+4) &
                     +fgz(k)+fgz(k+1)+fgz(k+2) )*dtm
!
      xg1= (16*xa(i) +xa(i+1) +xa(i+2))/18.d0 
      yg1= (16*ya(i) +ya(i+1) +ya(i+2))/18.d0 
      zg1= (16*za(i) +za(i+1) +za(i+2))/18.d0 
!
      xg(j)= xg1 +dt*vx(j)
      yg(j)= yg1 +dt*vy(j)
      zg(j)= zg1 +dt*vz(j)
      end do
!
      if(t8.le.t_pe+10.d0) then   ! retain salt until DNA contracts
        do j= np+1,npqr 
        if(abs(zg(j)).lt.0.5*Hpore) then
          vz(j)= 0.d0
        end if
        end do
      end if 
!
!    --------------------------------------------------------
      call reflect_endpl_P (xa,ya,za,xg,yg,zg,vx,vy,vz,ag, &
                            np,nq,npqr,npqr5)
!    --------------------------------------------------------
!
        cl_first= 2
        call clocks (wtime5,size,cl_first)
!
      if(iwrt1.eq.0 .and. io_pe.eq.1) then
!     if(.false.) then
! 
        OPEN (unit=11,file=praefixc//'.06'//suffix2,             &
              status='unknown',position='append',form='formatted')
!
        write(11,612) t8,wtime5-wtime0,wtime1-wtime0,wtime2-wtime1, &
                      wtime3-wtime2,wtime4-wtime3,wtime5-wtime4,    &
                      wtime0/60.d0,it
  612   format('# t8,wtime=',f6.1,f8.2,' s, cpu=',5f7.4,f7.1,i7)
!*
        write(11,*) '##### t8,Vtop,Vbot=',t8,Vtop,Vbot
        write(11,*)
!
        close(11)
      end if
!
!********************************
!*  Step 4A: Diagnosis section  *
!********************************
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
!* amm(j) is defined in subroutine /init/
!
        do j= nCLp+1,npqr 
        s2= s2 +0.5*amm(j)*(vx(j)**2 +vy(j)**2 +vz(j)**2)
        end do
!
        Ekin0= s0/(np +1.e-5)
        Ekin1= s1/(nq +1.e-5)
        Ekin2= s2/(nr +1.e-5)
!
        time(is)= t
        vdtm(is)= vm*dt
!
        ekin(is) = Ekin1   ! s1/nq
        ekn2(is) = Ekin2   ! s2/nr
        ecr (is) = e_c_r   ! e_c_r/npqr5
        elj (is) = e_lj    ! e_lj/npqr
        espr(is) = e_elas  ! e_elas/np
        ppot(is) = e_grid  ! e_grid/(mx*my*mz)
        etot(is) = 0.5d0*( s0 +s1 +s2 ) +ecr(is) +elj(is) &  ! total
                                        +espr(is) +ppot(is)
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
          write(11,600)
  600     format(/,'   time:     E_kin0/np  E_kin1/nq  E_kin2/nr  ', &
                   'E_coulb    E_grid     E_lj       E_elas     ',   &
                   'E_tot    iterp rsdl      wall(sec)  it')
        end if
!
        write(11,610) t8,Ekin0,Ekin1,Ekin2,ecr(is),ppot(is),       &
                      elj(is),espr(is),etot(is),ipr(2),rpr(2),     &
                      wtime5/60.,it 
  610   format('#t8=',f8.2,1p8d11.3,i5,d9.2,0pf8.1,i6) 
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
!       do i= 1,npqr5
!       x4(i)= xg(i)
!       y4(i)= yg(i)
!       z4(i)= zg(i)
!
!       vx4(i)= vx(i)
!       vy4(i)= vy(i)
!       vz4(i)= vz(i)
!       end do
!
        write(13)  t,xg,yg,zg,vx,vy,vz
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
          call charges (rho,xg,yg,zg,xa,ya,za,xe,ye,ze,ch,chg, &
                        A11,A12,A13,A21,A22,A23,A31,A32,A33,   &
                        xr,yr,zr,ql,gg,ipar,np,nq,nCLp,npqr,   &
                        npqr3,npqr5)
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
!* All PE must call /poissn_eq/ (MPI wait !)
!
          ndim= 3
          call poissn_eq (rho,pot,ndim,itermax,iterp,ipar)
!
!
          if(io_pe.eq.1) then
!+
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
!           OPEN (unit=16,file=praefixc//'.16'//suffix2,             &
!               status='unknown',position='append',form='unformatted')
!
!           mx1= (mx-1)/2-10
!           mx2= (mx-1)/2+10
!           my1= (my-1)/2-10
!           my2= (my-1)/2+10
!           mz1= 0
!           mz2= mz-1
!
!           write(16) t,mx1,mx2,my1,my2,mz1,mz2
!           write(16) (((rho(i,j,k),i= mx1,mx2),j=my1,my2),k=mz1,mz2)
!           write(16) (((pot(i,j,k),i= mx1,mx2),j=my1,my2),k=mz1,mz2)
!
!           close(16)
          end if
        end if
!
!**
!  Convert from 5-point to 3-point molecules
!
        if(t8.gt.1.d0 .and. io_pe.eq.1) then
!
          OPEN (unit=77,file=praefixc//'.77'//suffix2//'.ps',      &
                status='unknown',position='append',form='formatted')
!
          do i= 1,nCLp
          xe(i)= xg(i)
          ye(i)= yg(i)
          ze(i)= zg(i)
!
          chg(i)= ch(i)
          end do
!
          do i= nCLp+1,npqr3
          if(mod(i-nCLp,3).eq.1) then
            chg(i)= -1.0d0
          else
            chg(i)=  0.5d0
          end if
          end do
!
!  For plots by ppl3da and vdistr subroutines
!
          j= nCLp
          do i= nCLp+1,npqr5,5
          j= j +1
!
          xa(i  )= xg(j) +A11(j)*xr(i) +A21(j)*yr(i) +A31(j)*zr(i)
          ya(i  )= yg(j) +A12(j)*xr(i) +A22(j)*yr(i) +A32(j)*zr(i)
          za(i  )= zg(j) +A13(j)*xr(i) +A23(j)*yr(i) +A33(j)*zr(i)
!  
          xa(i+1)= xg(j) +A11(j)*xr(i+1) +A21(j)*yr(i+1) +A31(j)*zr(i+1)
          ya(i+1)= yg(j) +A12(j)*xr(i+1) +A22(j)*yr(i+1) +A32(j)*zr(i+1)
          za(i+1)= zg(j) +A13(j)*xr(i+1) +A23(j)*yr(i+1) +A33(j)*zr(i+1)
! 
          xa(i+2)= xg(j) +A11(j)*xr(i+2) +A21(j)*yr(i+2) +A31(j)*zr(i+2)
          ya(i+2)= yg(j) +A12(j)*xr(i+2) +A22(j)*yr(i+2) +A32(j)*zr(i+2)
          za(i+2)= zg(j) +A13(j)*xr(i+2) +A23(j)*yr(i+2) +A33(j)*zr(i+2)
! 
!  Use three point O-H-H
          k= nCLp+1 +3*(i-nCLp-1)/5
!
          xe(k)  = xa(i)
          ye(k)  = ya(i)
          ze(k)  = za(i)
!
          xe(k+1)= (xa(i+1) +xa(i+3))/2
          ye(k+1)= (ya(i+1) +ya(i+3))/2
          ze(k+1)= (za(i+1) +za(i+3))/2
!
          xe(k+2)= (xa(i+2) +xa(i+4))/2
          ye(k+2)= (ya(i+2) +ya(i+4))/2
          ze(k+2)= (za(i+2) +za(i+4))/2
          end do
!
          call ppl3da (xe,ye,ze,chg,rod_leng,Rmac,np,nq,nCLp, &
                       npqr3,first_ppl)    
!
          call vdistr (vx,vy,vz,np,nq,nCLp,nr,npqr) 
!
          close(77)
!
        end if
!*
      end if
!
! -------------------
!*  Restart data.
! -------------------
!
      if(iwrt2.eq.0 .and. io_pe.eq.1) then
!
        OPEN (unit=12,file=praefixc//'.12'//suffix2,     &
                      status='replace',form='unformatted')
!
        write(12) it,np,nq,nCLp,nr,npqr,npqr3,npqr5
        write(12) ifrgrod,ifrodco,ifbase
!
        write(12) t8,xg,yg,zg,vx,vy,vz,ch,am,amm,ag,ep
        write(12) xa,ya,za
!
        write(12) A11,A12,A13,A21,A22,A23,A31,A32,A33
        write(12) e0,e1,e2,e3,xr,yr,zr,Im,Lgx,Lgy,Lgz
!
        write(12) pi,dt,Gamma,rbmax,vth,tmax  !<- L.2160 IF
        write(12) t,phi,tht,dtwr1,dtwr2,dtwr3,is
        write(12) ekin,ppot,ekn2,etot,z_pe,vzpe,vzco,vzct, &
                  vdtm,vpor,time,ecr,elj,espr
!
        write(12) iwa,iwb,iwc,iwd
        write(12) nsg,nseg
!       write(12) nbox,list
!
        write(12) np00,nnb,ist1,ist2
        write(12) diel2,daa,dbb
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
      end if
!
! --------------------------------------- On major node --------------
!
!     if(istop.ge.1) then
!       if(io_pe.eq.1) then
!       OPEN (unit=11,file=praefixc//'.06'//suffix2,             &
!             status='unknown',position='append',form='formatted')
!
!       write(11,*) ' Uses Abnormal termination (istop=1)... '
!       write(11,*) '  ipar,t8=',ipar,t8
!       close(11)
!       end if
!       
!       go to 2000
!     end if
!                         +++  it= 1 2 3 4 5 6 7 8 9 10
      if(mod(it,ntimes).eq.0) then 
        if(t8.gt.tmax) go to 2000
        if((wtime0/60.d0).gt.cptot) go to 2000
      end if
!
      go to 1000
!
 2000 continue
      if(io_pe.eq.1) then
        OPEN (unit=11,file=praefixc//'.06'//suffix2,             &
              status='unknown',position='append',form='formatted')
!
        write(11,*) ' Final: t, tmax=',t8,tmax
        write(11,*) '   wtime5/60., cptot(min)=',wtime5/60.,cptot
!
        close(11)
      end if
!
      return
      end subroutine moldyn
!
!
!--------------------------------------------------------------------
      subroutine realteil (xa,ya,za,ch,xg,yg,zg,ag,ep,fcx,fcy,fcz, & 
                           ipar,size,np,nq,nCLp,npqr,npqr5)
!--------------------------------------------------------------------
      use, intrinsic :: iso_c_binding 
      implicit  none
!*
      include   'mpif.h'
      include   'paramWatQa.h'
!
      real(C_DOUBLE),dimension(npqr50) :: xa,ya,za,ch,fcx,fcy,fcz, &
                                          fdx,fdy,fdz 
      real(C_DOUBLE),dimension(npqr0) ::  xg,yg,zg,ag,ep
      integer(C_INT),dimension(npqr0) ::  ibx
      integer(C_INT) ipar,size,np,nq,nCLp,npqr,npqr5
!
!  Labels_P in L.3700
      integer(C_INT) isize,isizeZ,isize2,isize4,nc3
      parameter  (isize=12,isizeZ=25,nc3=isize**2*isizeZ)
      parameter  (isize2=isize*isize,isize4=isize2+isize)
! 
      integer(C_INT) it,is
      common/parm1/ it,is
!
      integer(C_INT),dimension(nc3) :: ncel
      integer(C_INT),dimension(nbxs,nc3) :: lcel
      integer(C_INT),dimension(npqr0) :: nipl
      integer(C_INT),dimension(nbxs,npqr0) :: lipl
      integer(C_INT),dimension(27,nc3) :: ibind
      common /boxind/ ibind 
!
      integer(C_INT) io_pe,i0,i1,i2,i1a
      common/sub_proc/ io_pe
      common/dat_typ1/ i0(30),i1(30),i2(30)
!
      real(C_DOUBLE) pi,dt,dth,axi,Gamma,rbmax,vth,tmax, &
                     xmax,ymax,zmax,xmin,ymin,zmin,      &
                     xleng,yleng,zleng,aLJ
      common/parm2/  pi,dt,dth,axi,Gamma,rbmax,vth,tmax
      common/parm3/  xmax,ymax,zmax,xmin,ymin,zmin
      common/parm8/  xleng,yleng,zleng
      common/shrgf/  aLJ
!
      real(C_DOUBLE) rcut_clf,rcutlj,r_fene,  &
                     e_c_s,e_c_pme,e_c_r,e_lj,e_elas, &
                     xaa,xbb,yaa,ybb,zaa,zbb
      common /confdatar/ rcut_clf,rcutlj,r_fene
      common /energy/ e_c_s,e_c_pme,e_c_r,e_lj,e_elas
!
      real(C_DOUBLE) xmax3,ymax3,zmax3,xmin3,ymin3,zmin3 
      real(C_DOUBLE) t8,cptot
      common/headr3/ t8,cptot
!
      real(C_float) phi,tht,dtwr1,dtwr2,dtwr3
      common/parm4/ phi,tht,dtwr1,dtwr2,dtwr3
!
      real(C_float)  time,t00,xp_leng
      common/headr2/ time,t00,xp_leng
!
      integer(C_INT) i,j,k,l,ia,id,ja,jd,jx,jy,jz,ibox,neigh,ii,jj, &
                     ix,iy,iz,ierrer,istop,ierror
      real(C_DOUBLE) driwu,driwu2,rlj,rsi,snt,ccel,unif1(2),unif2(2)
!
      integer(C_INT) iwrt1,iwrt2,iwrt3,iwrt4
      common/iotim/  iwrt1,iwrt2,iwrt3,iwrt4
      common/abterm/ istop
!
      real(C_DOUBLE) rcutcl2,dx,dy,dz,dx1,dy1,dz1,dx2,dy2,dz2,    &
                     dx3,dy3,dz3,forceV,forceV1,forceV2,forceV3,  &
                     r2,r,f_cut,diel,e_c_r0,e_c_r1,e_c_r2,e_c_r3, &
                     asc2,pqq,eps_lj,addpot,rljcut,rljcut2,       &
                     Lx3,Ly3,Lz3,aix,aiy
!
      real(C_DOUBLE) epsLJ,eps_K,eps_Cl
      integer(C_INT) ifLJ
      common/parmlj/ epsLJ,eps_K,eps_Cl
      common/parmlj2/ ifLJ
!
      real(C_DOUBLE) xgc,ygc,zgc,vxg,vyg,vzg,rod_leng,Rmac,Rhelix
      integer(C_INT) n_rodp
      common/macroion/ xgc,ygc,zgc,vxg,vyg,vzg,rod_leng,Rmac,  &
                       Rhelix,n_rodp
!
      real(C_DOUBLE),dimension(np0) :: &
                       dox,doy,doz,doxc,doyc,dozc
      common/macroin1/ dox,doy,doz,doxc,doyc,dozc
!-----------
      integer(C_INT) pxr,pxc,pxl,pyr,pyc,pyl,pzr,pzc,pzl
      common/ptable/ pxr(-10:mx+10),pxc(-10:mx+10),pxl(-10:mx+10), &
                     pyr(-10:my+10),pyc(-10:my+10),pyl(-10:my+10), &
                     pzr(-10:mz+10),pzc(-10:mz+10),pzl(-10:mz+10)
!
      real(C_DOUBLE) gx,gy,gz,ghx,ghy,ghz,ghx2,ghy2,ghz2,  &
                     ghxsq,ghysq,ghzsq,hxi,hyi,hzi
      integer(C_INT) ptx,pty,ptz
      common/xtable/ gx(0:mx),ghx(0:mx),ghx2(0:mx),ghxsq(0:mx), &
                     gy(0:my),ghy(0:my),ghy2(0:my),ghysq(0:my), &
                     gz(0:mz),ghz(0:mz),ghz2(0:mz),ghzsq(0:mz), &
                     hxi,hyi,hzi
      common/xtabl2/ ptx(-199:3000),pty(-199:3000),ptz(-100:4000)
!*---------------------------------------------------------------
!
      character(len=8) :: fort51
      common/fort_write/  fort51(6)
!
!*---------------------------------------------------------
!*  Update interaction table in every << 5 steps.>> 
!*---------------------------------------------------------
!* Cutoff for table listing
! 
      rcutcl2= (2.d0* rcut_clf)**2    ! Cutoff of Coulomb forces
!     rcutlj =                        !<- common/confdatar/ rcut_clf,rcutlj
!
      driwu2 = 1.25992104989487316476721060728d0  ! 2**(1/3)
      driwu  = sqrt(driwu2)
!
      asc2 = (0.85d0*aLJ)**2 
!
      xmax3= xmax -0.01d0
      ymax3= ymax -0.01d0
      zmax3= zmax -0.5d0*zleng/isizeZ
!
      xmin3= xmin +0.01d0
      ymin3= ymin +0.01d0
      zmin3= zmin +0.5d0*zleng/isizeZ     !<-- (-0.5, 0.5)*zleng/isizeZ
!                                               index of 1 
      Lx3= xmax -xmin
      Ly3= ymax -ymin
      Lz3= zmax -zmin
!
!**************************************************************
!* Step 1R. Note: Large macroions must be treated separately  *
!**************************************************************
!
      if(mod(it,5).eq.1) then
!     +++++++++++++++++++++++
!
!* Step 1: Find particles 
! 
!    -------------------
        call Labels_P
!    -------------------
!
        do k= 1,nc3
        ncel(k)= 0     ! Clear the (ix,iy,iz) cell registers.
        end do
!
        do i = 1,npqr
        ix= pxc(int(isize* (xg(i)-xmin)/Lx3 +1.0001))  !<-- Periodic
        iy= pyc(int(isize* (yg(i)-ymin)/Ly3 +1.0001))
        iz=     int(isizeZ*(zg(i)-zmin)/Lz3 +1.5001)   !<-- Bounded
        if(iz.le.1 .or. iz.ge.isizeZ) go to 100
!
        ibox = ix + isize*(iy-1 + isize*(iz-1))
! 
        ncel(ibox)= ncel(ibox) +1     ! Number of targets is registered
        lcel(ncel(ibox),ibox)= i      ! Large bin is used
!
  100   end do
!
!
!* Step 2: Register particles around i-th particle
!
        if(ipar.eq.1) then
!*
          do i= 1,nCLp
          nipl(i)= 0
!
          ix= pxc(int(isize* (xg(i)-xmin)/Lx3 +1.0001))  !<-- Periodic
          iy= pyc(int(isize* (yg(i)-ymin)/Ly3 +1.0001)) 
          iz=     int(isizeZ*(zg(i)-zmin)/Lz3 +1.5001)   !<-- Bounded   
          if(iz.le.1 .or. iz.ge.isizeZ) go to 110
!
          ibox = ix + isize*(iy-1 + isize*(iz-1))
!
          do k = 1,27                  ! 130:
          neigh= ibind(k,ibox)
          if(neigh.le.0) go to 130     ! Far away
!
          do l= 1,ncel(neigh)          ! 140: Find ions in the boxes around i-th
          j= lcel(l,neigh)             ! j-th belongs to this box.
!
          if(i.eq.j) go to 140
! 
          dx= xg(i) -xg(j) 
          dy= yg(i) -yg(j)
          dz= zg(i) -zg(j)
!
          dx = dx - nint(dx/xleng -0.5d0)*xleng 
          dy = dy - nint(dy/yleng -0.5d0)*yleng 
!***
          r2 = dx**2 +dy**2 +dz**2
          if(r2.lt.rcutcl2) then
!
            nipl(i)= nipl(i) +1 
            lipl(nipl(i),i)= j
!
          end if
  140     end do
  130     end do
!
!  Occasionaly, nipl(i)= 0 happens resulting in null operation 
!  at this 5-step cycle !
!
!         if(nipl(i).eq.0) go to 110
!
  110     end do
!*
        end if
!
!  This step must execute with any ipar step 
!
        do i= i1(ipar),i2(ipar)
        nipl(i)= 0
!
        ix= pxc(int(isize* (xg(i)-xmin)/Lx3 +1.0001))  !<-- Periodic
        iy= pyc(int(isize* (yg(i)-ymin)/Ly3 +1.0001)) 
        iz=     int(isizeZ*(zg(i)-zmin)/Lz3 +1.5001)   !<-- Bounded
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
        dx = dx - nint(dx/xleng -0.5d0)*xleng 
        dy = dy - nint(dy/yleng -0.5d0)*yleng 
!***
        r2 = dx**2 +dy**2 +dz**2
        if(r2.lt.rcutcl2) then
!
          nipl(i)= nipl(i) +1 
          lipl(nipl(i),i)= j   !<<-- debugger
!
        end if
  240   end do
  230   end do
!
!       if(nipl(i).eq.0) go to 210
!
  210   end do
!*
        end if
!
!************************************
!*  Step 3R: The Coulomb forces.    *
!************************************
!   real(C_DOUBLE),dimension(npqr50) :: xa,ya,za,ch,fcx,fcy,fcz, &
!   real(C_DOUBLE),dimension(npqr0) ::  xg,yg,zg
!
      e_c_r= 0
!
      if(ipar.eq.1) then
!*
        do i = 1,nCLp
        if(zg(i).lt.zmin3 .or. zg(i).gt.zmax3) go to 50
        if(nipl(i).eq.0) go to 50
!
        do jd= 1,nipl(i)
        jj= lipl(jd,i)
!
        if(jj.eq.i) go to 300
!
        jz= isizeZ*(zg(jj)-zmin)/Lz3 +1.5001
        if(jz.le.1 .or. jz.ge.isizeZ) go to 300
!
        if(jj.le.nCLp) then  ! xg(j)
          j= jj 
          dx= xg(i) -xg(j)
          dy= yg(i) -yg(j)
          dz= zg(i) -zg(j)
!
          dx = dx - nint(dx/xleng -0.5d0)*xleng 
          dy = dy - nint(dy/yleng -0.5d0)*yleng 
        else
!                          <- H(2),H(3),L(4),L(5)
          j = nCLp +5*(jj-nCLp)-3  ! H(2)
          dx= xg(i) -xa(j) 
          dy= yg(i) -ya(j) 
          dz= zg(i) -za(j)
!
          dx = dx - nint(dx/xleng -0.5d0)*xleng 
          dy = dy - nint(dy/yleng -0.5d0)*yleng 
        end if
!
! (1a)
        r2 = max(asc2, dx**2 +dy**2 +dz**2)
        r  = sqrt(r2)
        if(r2.gt.rcutcl2) go to 300
!          +++++++++++++
!
        if(jj.le.nCLp) then                  !  ch(j)
!
          pqq  = f_cut(r,rcut_clf) *Gamma*ch(i)*ch(j)/r2         &
                  *diel((xg(i)+xg(j))/2.d0,(yg(i)+yg(j))/2.d0, &
                                           (zg(i)+zg(j))/2.d0)
          e_c_r0 = pqq/r
          forceV = pqq/(r2*r)
!
          fcx(i) = fcx(i) +forceV*dx
          fcy(i) = fcy(i) +forceV*dy
          fcz(i) = fcz(i) +forceV*dz
          e_c_r = e_c_r + e_c_r0
        else
!
          pqq  = f_cut(r,rcut_clf) *Gamma*ch(i)*ch(j)/r2       &
                  *diel((xg(i)+xa(j))/2.d0,(yg(i)+ya(j))/2.d0, &
                                           (zg(i)+za(j))/2.d0)
!
          e_c_r0 = pqq/r
          forceV = pqq/(r2*r)
!
! (1b)
          j= nCLp +5*(jj-nCLp)-2 
          dx1= xg(i) -xa(j)
          dy1= yg(i) -ya(j)
          dz1= zg(i) -za(j)
!
          dx1= dx1- nint(dx1/xleng -0.5d0)*xleng 
          dy1= dy1- nint(dy1/yleng -0.5d0)*yleng 
!
          r2 = max(asc2, dx1**2 +dy1**2 +dz1**2)
          r  = sqrt(r2)
!
          pqq  = f_cut(r,rcut_clf) *Gamma*ch(i)*ch(j)/r2       &
                  *diel((xg(i)+xa(j))/2.d0,(yg(i)+ya(j))/2.d0, &
                                           (zg(i)+za(j))/2.d0)
          e_c_r1 = pqq/r
          forceV1= pqq/(r2*r)
!
! (1c) 
          j= nCLp +5*(jj-nCLp)-1 
          dx2= xg(i) -xa(j)
          dy2= yg(i) -ya(j)
          dz2= zg(i) -za(j)
!
          dx2= dx2- nint(dx2/xleng -0.5d0)*xleng 
          dy2= dy2- nint(dy2/yleng -0.5d0)*yleng 
!
          r2 = max(asc2, dx2**2 +dy2**2 +dz2**2)
          r  = sqrt(r2)
!
          pqq  = f_cut(r,rcut_clf) *Gamma*ch(i)*ch(j)/r2       &
                  *diel((xg(i)+xa(j))/2.d0,(yg(i)+ya(j))/2.d0, &
                                           (zg(i)+za(j))/2.d0)
          e_c_r2 = pqq/r
          forceV2= pqq/(r2*r)
!
! (1d) 
          j= nCLp +5*(jj-nCLp)
          dx3= xg(i) -xa(j)
          dy3= yg(i) -ya(j)
          dz3= zg(i) -za(j)
!
          dx3= dx3- nint(dx3/xleng -0.5d0)*xleng 
          dy3= dy3- nint(dy3/yleng -0.5d0)*yleng 
!
          r2 = max(asc2, dx3**2 +dy3**2 +dz3**2)
          r  = sqrt(r2)
!
          pqq  = f_cut(r,rcut_clf) *Gamma*ch(i)*ch(j)/r2       &
                  *diel((xg(i)+xa(j))/2.d0,(yg(i)+ya(j))/2.d0, &
                                           (zg(i)+za(j))/2.d0)
          e_c_r3 = pqq/r
          forceV3= pqq/(r2*r)
!
          fcx(i) = fcx(i) +forceV*dx +forceV1*dx1 +forceV2*dx2 +forceV3*dx3
          fcy(i) = fcy(i) +forceV*dy +forceV1*dy1 +forceV2*dy2 +forceV3*dy3
          fcz(i) = fcz(i) +forceV*dz +forceV1*dz1 +forceV2*dz2 +forceV3*dz3
!
          e_c_r = e_c_r +e_c_r0 +e_c_r1 +e_c_r2 +e_c_r3
        end if
  300   end do
   50   continue
!
        end do
!*
      end if
!
!************************************
!*  Step 3S: The Coulomb forces (2) *
!************************************
!
      do l = i1(ipar),i2(ipar)
!                                <-- do reflect_endplP 
      if(zg(l).lt.zmin3 .or. zg(l).gt.zmax3) go to 60
      if(nipl(l).eq.0) go to 60
!
! -----------------------------
      i= nCLp +5*(l-nCLp)-3
! -----------------------------
!
      do jd= 1,nipl(l)
      jj= lipl(jd,l)
!
      if(jj.eq.l) go to 400
!
      jz= isizeZ*(zg(jj)-zmin)/Lz3 +1.5001
      if(jz.le.1 .or. jz.ge.isizeZ) go to 400
!
      if(jj.le.nCLp) then
        j= jj 
        dx= xa(i) -xg(j)
        dy= ya(i) -yg(j)
        dz= za(i) -zg(j)
!
        dx = dx - nint(dx/xleng -0.5d0)*xleng 
        dy = dy - nint(dy/yleng -0.5d0)*yleng 
      else
!                          <- H(2),H(3),L(4),L(5)
        j= nCLp +5*(jj-nCLp)-3
        dx= xa(i) -xa(j) 
        dy= ya(i) -ya(j)
        dz= za(i) -za(j)
!
        dx = dx - nint(dx/xleng -0.5d0)*xleng 
        dy = dy - nint(dy/yleng -0.5d0)*yleng 
      end if
!
      r2 = max(asc2, dx**2 +dy**2 +dz**2)
      r = sqrt(r2)
      if(r2.gt.rcutcl2) go to 400
!        ++++++++++++++++
!
! (1a) 
      if(jj.le.nCLp) then
!
        pqq  = f_cut(r,rcut_clf) *Gamma*ch(i)*ch(j)/r2         &
                *diel((xa(i)+xg(j))/2.d0,(ya(i)+yg(j))/2.d0, &
                                         (za(i)+zg(j))/2.d0)
        e_c_r0 = pqq/r
        forceV = pqq/(r2*r)
!
        fcx(i) = fcx(i) +forceV*dx
        fcy(i) = fcy(i) +forceV*dy
        fcz(i) = fcz(i) +forceV*dz
        e_c_r = e_c_r + e_c_r0
!
      else
!
        pqq  = f_cut(r,rcut_clf) *Gamma*ch(i)*ch(j)/r2       &
                *diel((xa(i)+xa(j))/2.d0,(ya(i)+ya(j))/2.d0, &
                                         (za(i)+za(j))/2.d0)
        e_c_r0 = pqq/r
        forceV = pqq/(r2*r)
!
! (1b)
        j= nCLp +5*(jj-nCLp) -2
        dx1= xa(i) -xa(j)
        dy1= ya(i) -ya(j)
        dz1= za(i) -za(j)
!
        dx1= dx1- nint(dx1/xleng -0.5d0)*xleng 
        dy1= dy1- nint(dy1/yleng -0.5d0)*yleng 
!
        r2 = max(asc2, dx1**2 +dy1**2 +dz1**2)
        r  = sqrt(r2)
!
        pqq  = f_cut(r,rcut_clf) *Gamma*ch(i)*ch(j)/r2       &
                *diel((xa(i)+xa(j))/2.d0,(ya(i)+ya(j))/2.d0, &
                                         (za(i)+za(j))/2.d0)
        e_c_r1 = pqq/r
        forceV1= pqq/(r2*r)
!
! (1c) 
        j= nCLp +5*(jj-nCLp) -1
        dx2= xa(i) -xa(j)
        dy2= ya(i) -ya(j)
        dz2= za(i) -za(j)
!
        dx2= dx2- nint(dx2/xleng -0.5d0)*xleng 
        dy2= dy2- nint(dy2/yleng -0.5d0)*yleng 
!
        r2 = max(asc2, dx2**2 +dy2**2 +dz2**2)
        r  = sqrt(r2)
!
        pqq  = f_cut(r,rcut_clf) *Gamma*ch(i)*ch(j)/r2       &
                *diel((xa(i)+xa(j))/2.d0,(ya(i)+ya(j))/2.d0, &
                                         (za(i)+za(j))/2.d0)
        e_c_r2 = pqq/r
        forceV2= pqq/(r2*r)
!
! (1d) 
        j= nCLp +5*(jj-nCLp)
        dx3= xa(i) -xa(j)
        dy3= ya(i) -ya(j)
        dz3= za(i) -za(j)
!
        dx3= dx3- nint(dx3/xleng -0.5d0)*xleng 
        dy3= dy3- nint(dy3/yleng -0.5d0)*yleng 
!
        r2 = max(asc2, dx3**2 +dy3**2 +dz3**2)
        r  = sqrt(r2)
!
        pqq  = f_cut(r,rcut_clf) *Gamma*ch(i)*ch(j)/r2       &
                *diel((xa(i)+xa(j))/2.d0,(ya(i)+ya(j))/2.d0, &
                                         (za(i)+za(j))/2.d0)
        e_c_r3 = pqq/r
        forceV3= pqq/(r2*r)
!
        fcx(i) = fcx(i) +forceV*dx +forceV1*dx1 +forceV2*dx2 +forceV3*dx3
        fcy(i) = fcy(i) +forceV*dy +forceV1*dy1 +forceV2*dy2 +forceV3*dy3
        fcz(i) = fcz(i) +forceV*dz +forceV1*dz1 +forceV2*dz2 +forceV3*dz3
!
        e_c_r = e_c_r +e_c_r0 +e_c_r1 +e_c_r2
      end if
!
! (2a)
! -----------------------------
      i= nCLp +5*(l-nCLp) -2 
! -----------------------------
!
      if(jj.le.nCLp) then
        j= jj
        dx= xa(i) -xg(j)
        dy= ya(i) -yg(j)
        dz= za(i) -zg(j)
!
        dx = dx - nint(dx/xleng -0.5d0)*xleng 
        dy = dy - nint(dy/yleng -0.5d0)*yleng 
!
        r2 = max(asc2, dx**2 +dy**2 +dz**2)
        r  = sqrt(r2)
        pqq  = f_cut(r,rcut_clf) *Gamma*ch(i)*ch(j)/r2         &
                *diel((xa(i)+xg(j))/2.d0,(ya(i)+yg(j))/2.d0, &
                                         (za(i)+zg(j))/2.d0)
        e_c_r0 = pqq/r
        forceV = pqq/(r2*r)
!
        fcx(i) = fcx(i) +forceV*dx
        fcy(i) = fcy(i) +forceV*dy
        fcz(i) = fcz(i) +forceV*dz
        e_c_r = e_c_r + e_c_r0
!
      else
        j= nCLp +5*(jj-nCLp)-3
        dx= xa(i) -xa(j)
        dy= ya(i) -ya(j)
        dz= za(i) -za(j)
!
        dx = dx - nint(dx/xleng -0.5d0)*xleng 
        dy = dy - nint(dy/yleng -0.5d0)*yleng 
!
        r2 = max(asc2, dx**2 +dy**2 +dz**2)
        r  = sqrt(r2)
        pqq  = f_cut(r,rcut_clf) *Gamma*ch(i)*ch(j)/r2       &
                *diel((xa(i)+xa(j))/2.d0,(ya(i)+ya(j))/2.d0, &
                                         (za(i)+za(j))/2.d0)
        e_c_r0 = pqq/r
        forceV = pqq/(r2*r)
!
! (2b) 
        j= nCLp +5*(jj-nCLp)-2
        dx1= xa(i) -xa(j)
        dy1= ya(i) -ya(j)
        dz1= za(i) -za(j)
!
        dx1= dx1- nint(dx1/xleng -0.5d0)*xleng 
        dy1= dy1- nint(dy1/yleng -0.5d0)*yleng 
!
        r2 = max(asc2, dx1**2 +dy1**2 +dz1**2)
        r  = sqrt(r2)
!
        pqq  = f_cut(r,rcut_clf) *Gamma*ch(i)*ch(j)/r2       &
                *diel((xa(i)+xa(j))/2.d0,(ya(i)+ya(j))/2.d0, &
                                         (za(i)+za(j))/2.d0)
        e_c_r1 = pqq/r
        forceV1= pqq/(r2*r)
!
! (2c)
        j= nCLp +5*(jj-nCLp)-1 
        dx2= xa(i) -xa(j)
        dy2= ya(i) -ya(j)
        dz2= za(i) -za(j)
!
        dx2= dx2- nint(dx2/xleng -0.5d0)*xleng 
        dy2= dy2- nint(dy2/yleng -0.5d0)*yleng 
!
        r2 = max(asc2, dx2**2 +dy2**2 +dz2**2)
        r  = sqrt(r2)
!
        pqq  = f_cut(r,rcut_clf) *Gamma*ch(i)*ch(j)/r2        &
                *diel((xa(i)+xa(j))/2.d0,(ya(i)+ya(j))/2.d0,  &
                                         (za(i)+za(j))/2.d0)
        e_c_r2 = pqq/r
        forceV2= pqq/(r2*r)
!
! (2d)
        j= nCLp +5*(jj-nCLp) 
        dx3= xa(i) -xa(j)
        dy3= ya(i) -ya(j)
        dz3= za(i) -za(j)
!
        dx3= dx3- nint(dx3/xleng -0.5d0)*xleng 
        dy3= dy3- nint(dy3/yleng -0.5d0)*yleng 
!
        r2 = max(asc2, dx3**2 +dy3**2 +dz3**2)
        r  = sqrt(r2)
!
        pqq  = f_cut(r,rcut_clf) *Gamma*ch(i)*ch(j)/r2        &
                *diel((xa(i)+xa(j))/2.d0,(ya(i)+ya(j))/2.d0,  &
                                         (za(i)+za(j))/2.d0)
        e_c_r3 = pqq/r
        forceV3= pqq/(r2*r)
!*
        fcx(i) = fcx(i) +forceV*dx +forceV1*dx1 +forceV2*dx2 +forceV3*dx3
        fcy(i) = fcy(i) +forceV*dy +forceV1*dy1 +forceV2*dy2 +forceV3*dy3
        fcz(i) = fcz(i) +forceV*dz +forceV1*dz1 +forceV2*dz2 +forceV3*dz3
!
        e_c_r = e_c_r +e_c_r0 +e_c_r1 +e_c_r2
      end if
!
! (3a) 
! ----------------------------
      i= nCLp +5*(l-nCLp) -1
! ----------------------------
!
      if(jj.le.nCLp) then
        j= jj
        dx= xa(i) -xg(j)
        dy= ya(i) -yg(j)
        dz= za(i) -zg(j)
!
        dx = dx - nint(dx/xleng -0.5d0)*xleng 
        dy = dy - nint(dy/yleng -0.5d0)*yleng 
!
        r2 = max(asc2, dx**2 +dy**2 +dz**2)
        r  = sqrt(r2)
        pqq  = f_cut(r,rcut_clf) *Gamma*ch(i)*ch(j)/r2       &
                *diel((xa(i)+xg(j))/2.d0,(ya(i)+yg(j))/2.d0, &
                                         (za(i)+zg(j))/2.d0)
        e_c_r0 = pqq/r
        forceV = pqq/(r2*r)
!
        fcx(i) = fcx(i) +forceV*dx
        fcy(i) = fcy(i) +forceV*dy
        fcz(i) = fcz(i) +forceV*dz
        e_c_r = e_c_r + e_c_r0
!
      else
        j= nCLp +5*(jj-nCLp)-3
        dx= xa(i) -xa(j)
        dy= ya(i) -ya(j)
        dz= za(i) -za(j)
!
        dx = dx - nint(dx/xleng -0.5d0)*xleng 
        dy = dy - nint(dy/yleng -0.5d0)*yleng 
!
        r2 = max(asc2, dx**2 +dy**2 +dz**2)
        r  = sqrt(r2)
        pqq  = f_cut(r,rcut_clf) *Gamma*ch(i)*ch(j)/r2       &
                *diel((xa(i)+xa(j))/2.d0,(ya(i)+ya(j))/2.d0, &
                                         (za(i)+za(j))/2.d0)
        e_c_r0 = pqq/r
        forceV = pqq/(r2*r)
!
! (3b) 
        j= nCLp +5*(jj-nCLp)-2
        dx1= xa(i) -xa(j)
        dy1= ya(i) -ya(j)
        dz1= za(i) -za(j)
!
        dx1= dx1- nint(dx1/xleng -0.5d0)*xleng 
        dy1= dy1- nint(dy1/yleng -0.5d0)*yleng 
!
        r2 = max(asc2, dx1**2 +dy1**2 +dz1**2)
        r  = sqrt(r2)
!
        pqq  = f_cut(r,rcut_clf) *Gamma*ch(i)*ch(j)/r2       &
                *diel((xa(i)+xa(j))/2.d0,(ya(i)+ya(j))/2.d0, &
                                         (za(i)+za(j))/2.d0)
        e_c_r1 = pqq/r
        forceV1= pqq/(r2*r)
!
! (3c) 
        j= nCLp +5*(jj-nCLp)-1 
        dx2= xa(i) -xa(j)
        dy2= ya(i) -ya(j)
        dz2= za(i) -za(j)
!
        dx2= dx2- nint(dx2/xleng -0.5d0)*xleng 
        dy2= dy2- nint(dy2/yleng -0.5d0)*yleng 
!
        r2 = max(asc2, dx2**2 +dy2**2 +dz2**2)
        r  = sqrt(r2)
!
        pqq  = f_cut(r,rcut_clf) *Gamma*ch(i)*ch(j)/r2       &
                *diel((xa(i)+xa(j))/2.d0,(ya(i)+ya(j))/2.d0, &
                                         (za(i)+za(j))/2.d0)
        e_c_r2 = pqq/r
        forceV2= pqq/(r2*r)
!
! (3d) 
        j= nCLp +5*(jj-nCLp)
        dx3= xa(i) -xa(j)
        dy3= ya(i) -ya(j)
        dz3= za(i) -za(j)
!
        dx3= dx3- nint(dx3/xleng -0.5d0)*xleng 
        dy3= dy3- nint(dy3/yleng -0.5d0)*yleng 
!
        r2 = max(asc2, dx3**2 +dy3**2 +dz3**2)
        r  = sqrt(r2)
!
        pqq  = f_cut(r,rcut_clf) *Gamma*ch(i)*ch(j)/r2       &
                *diel((xa(i)+xa(j))/2.d0,(ya(i)+ya(j))/2.d0, &
                                         (za(i)+za(j))/2.d0)
        e_c_r3 = pqq/r
        forceV3= pqq/(r2*r)
!*
        fcx(i) = fcx(i) +forceV*dx +forceV1*dx1 +forceV2*dx2 +forceV3*dx3
        fcy(i) = fcy(i) +forceV*dy +forceV1*dy1 +forceV2*dy2 +forceV3*dy3
        fcz(i) = fcz(i) +forceV*dz +forceV1*dz1 +forceV2*dz2 +forceV3*dz3
!
        e_c_r = e_c_r +e_c_r0 +e_c_r1 +e_c_r2
      end if
!
! (4a) 
! ----------------------------
      i= nCLp +5*(l-nCLp)
! ----------------------------
!
      if(jj.le.nCLp) then
        j= jj
        dx= xa(i) -xg(j)
        dy= ya(i) -yg(j)
        dz= za(i) -zg(j)
!
        dx = dx - nint(dx/xleng -0.5d0)*xleng 
        dy = dy - nint(dy/yleng -0.5d0)*yleng 
!
        r2 = max(asc2, dx**2 +dy**2 +dz**2)
        r  = sqrt(r2)
        pqq  = f_cut(r,rcut_clf) *Gamma*ch(i)*ch(j)/r2       &
                *diel((xa(i)+xg(j))/2.d0,(ya(i)+yg(j))/2.d0, &
                                         (za(i)+zg(j))/2.d0)
        e_c_r0 = pqq/r
        forceV = pqq/(r2*r)
!
        fcx(i) = fcx(i) +forceV*dx
        fcy(i) = fcy(i) +forceV*dy
        fcz(i) = fcz(i) +forceV*dz
        e_c_r = e_c_r + e_c_r0
!
      else
        j= nCLp +5*(jj-nCLp)-3
        dx= xa(i) -xa(j)
        dy= ya(i) -ya(j)
        dz= za(i) -za(j)
!
        dx = dx - nint(dx/xleng -0.5d0)*xleng 
        dy = dy - nint(dy/yleng -0.5d0)*yleng 
!
        r2 = max(asc2, dx**2 +dy**2 +dz**2)
        r  = sqrt(r2)
        pqq  = f_cut(r,rcut_clf) *Gamma*ch(i)*ch(j)/r2         &
                *diel((xa(i)+xa(j))/2.d0,(ya(i)+ya(j))/2.d0, &
                                         (za(i)+za(j))/2.d0)
        e_c_r0 = pqq/r
        forceV = pqq/(r2*r)
!
! (4b) 
        j= nCLp +5*(jj-nCLp)-2
        dx1= xa(i) -xa(j)
        dy1= ya(i) -ya(j)
        dz1= za(i) -za(j)
!
        dx1= dx1- nint(dx1/xleng -0.5d0)*xleng 
        dy1= dy1- nint(dy1/yleng -0.5d0)*yleng 
!
        r2 = max(asc2, dx1**2 +dy1**2 +dz1**2)
        r  = sqrt(r2)
!
        pqq  = f_cut(r,rcut_clf) *Gamma*ch(i)*ch(j)/r2       &
                *diel((xa(i)+xa(j))/2.d0,(ya(i)+ya(j))/2.d0, &
                                         (za(i)+za(j))/2.d0)
        e_c_r1 = pqq/r
        forceV1= pqq/(r2*r)
!
! (4c) 
        j= nCLp +5*(jj-nCLp)-1 
        dx2= xa(i) -xa(j)
        dy2= ya(i) -ya(j)
        dz2= za(i) -za(j)
!
        dx2= dx2- nint(dx2/xleng -0.5d0)*xleng 
        dy2= dy2- nint(dy2/yleng -0.5d0)*yleng 
!
        r2 = max(asc2, dx2**2 +dy2**2 +dz2**2)
        r  = sqrt(r2)
!
        pqq  = f_cut(r,rcut_clf) *Gamma*ch(i)*ch(j)/r2       &
                *diel((xa(i)+xa(j))/2.d0,(ya(i)+ya(j))/2.d0, &
                                         (za(i)+za(j))/2.d0)
        e_c_r2 = pqq/r
        forceV2= pqq/(r2*r)
!
! (4d) 
        j= nCLp +5*(jj-nCLp)
        dx3= xa(i) -xa(j)
        dy3= ya(i) -ya(j)
        dz3= za(i) -za(j)
!
        dx3= dx3- nint(dx3/xleng -0.5d0)*xleng 
        dy3= dy3- nint(dy3/yleng -0.5d0)*yleng 
!
        r2 = max(asc2, dx3**2 +dy3**2 +dz3**2)
        r  = sqrt(r2)
!
        pqq  = f_cut(r,rcut_clf) *Gamma*ch(i)*ch(j)/r2       &
                *diel((xa(i)+xa(j))/2.d0,(ya(i)+ya(j))/2.d0, &
                                         (za(i)+za(j))/2.d0)
        e_c_r3 = pqq/r
        forceV3= pqq/(r2*r)
!*
        fcx(i) = fcx(i) +forceV*dx +forceV1*dx1 +forceV2*dx2 +forceV3*dx3
        fcy(i) = fcy(i) +forceV*dy +forceV1*dy1 +forceV2*dy2 +forceV3*dy3
        fcz(i) = fcz(i) +forceV*dz +forceV1*dz1 +forceV2*dz2 +forceV3*dz3
!
        e_c_r = e_c_r +e_c_r0 +e_c_r1 +e_c_r2
      end if
  400 end do
   60 continue
!
      end do
!
!*************************************
!* Step 4: Lennard-Jones potential   *
!*************************************
!
      if(ipar.eq.1) then
        i1a= i0(ipar)
      else
        i1a= i1(ipar)
      end if
!
      do id = i1a,i2(ipar) 
      if(id.le.nCLp) then
        i= id
        xaa= xg(i)
        yaa= yg(i)
        zaa= zg(i)
      else
        i= nCLp +5*(id-nCLp) -4
        xaa= xa(i)
        yaa= ya(i)
        zaa= za(i) 
      end if
! 
      do jd= 1,npqr       ! np+nq+nr/5
      if(jd.le.nCLp) then
        j= jd
        xbb= xg(j)
        ybb= yg(j)
        zbb= zg(j)
      else
        j= nCLp +5*(jd-nCLp) -4 
        xbb= xa(j)
        ybb= ya(j)
        zbb= za(j)
      end if
!
      if(i.eq.j) go to 500
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
!
      dx = xaa -xbb
      dy = yaa -ybb
      dz = zaa -zbb
!
      dx = dx - nint(dx/xleng -0.5d0)*xleng 
      dy = dy - nint(dy/yleng -0.5d0)*yleng 
!
!*  Lennard-Jones force.
!   >> Energy unit: epsLJ (LJ energy).
!   >> rlj= 1.0 when i and j are touching.
!
      r2 = dx**2 + dy**2 + dz**2
      rlj = sqrt(r2)/(ag(id)+ag(jd))
!
      if(rlj.le.rljcut) then
        rsi = 1.d0/max(rlj**2,asc2)
        snt = rsi*rsi*rsi
!
        eps_lj= sqrt(ep(id)*ep(jd))
        ccel  = 48.d0*eps_lj*snt*(snt-0.5d0)/r2
!
        fcx(i) = fcx(i) + ccel*dx  !<-- fcx(i), i=1,6,... for LJ
        fcy(i) = fcy(i) + ccel*dy
        fcz(i) = fcz(i) + ccel*dz
!
        e_lj = e_lj  + 4.d0*eps_lj*(snt*(snt - 1.d0) - addpot)
      end if
!
  500 end do
      end do
!
! ---------------------
!*  Unify the force.
! ---------------------
!  fcx() spreads on all domains - must take "allreduce".....
!
      if(num_proc.ne.1) then
!
        call mpi_allreduce (fcx,fdx,npqr5,mpi_real8,mpi_sum, &
                            mpi_comm_world,ierror)
        call mpi_allreduce (fcy,fdy,npqr5,mpi_real8,mpi_sum, &
                            mpi_comm_world,ierror)
        call mpi_allreduce (fcz,fdz,npqr5,mpi_real8,mpi_sum, &
                            mpi_comm_world,ierror)
!               +++++ 
        do i= 1,npqr5 
        fcx(i)= fdx(i)
        fcy(i)= fdy(i)
        fcz(i)= fdz(i)
        end do
!
!
        unif1(1)= e_c_r
        unif1(2)= e_lj
        call mpi_allreduce (unif1,unif2,2,mpi_real8,mpi_sum, &
                            mpi_comm_world,ierror)
        e_c_r= unif2(1)
        e_lj = unif2(2)
!
        if(iwrt1.eq.0 .and. io_pe.eq.1) then
          OPEN (unit=11,file=praefixc//'.06'//suffix2,             &
                status='unknown',position='append',form='formatted')
          write(11,*)
          write(11,700) e_c_r/npqr5,e_lj/npqr
  700     format('# e_c_r/npqr5, e_lj/npqr=',1p2d12.3)
!
          close(11)
        end if
!
      end if
!
      return
      end subroutine realteil
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
      end if         ! rcut_clf=5 -> r/5
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
      include 'paramWatQa.h'
!
      integer(C_INT) ifqq
      real(C_DOUBLE) xg,yg,zg,rr,diel,dielv,diel_m,     &
                     qfrac,Rpore,Hpore,Hpore2,Zci,Zcp,Zcn
      real(C_DOUBLE) diel2,dielpr,dielp,tau,ddv,daa,dbb,r1,z1
!
      common/cntion/ qfrac,Rpore,Hpore,Zci,Zcp,Zcn,ifqq
      common/dielec/ diel2,dielpr,daa,dbb
!
      integer(C_INT) io_pe
      common/sub_proc/ io_pe
      logical        first
      data           first /.true./
!
! -------------------------
!     daa =   2.0d0   ! 1.5 Ang 
!     dbb =   3.0d0   ! vertical depth
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
          write(11,740) dielp
          write(11,750) daa,dbb
          write(11,760) dielv,dielpr/79.,diel_m
  740     format(' dielp=',f8.2)
  750     format(' daa,dbb=',2f8.2)
  760     format(' diel[v,p,m]=',3f8.2)
!
          close(11)
        end if
!
        first= .false.
      end if
!
      Hpore2= 0.5d0*Hpore
      r1 =  Rpore  -daa
      z1 =  Hpore2 +dbb
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
            diel= ddv -(ddv -diel_m)*min((rr-r1)/daa, 1.d0) !<-- pore
          end if
        end if
!
! 4) Topside of a membrane
        if(abs(zg).ge.Hpore2) then
          if(rr.gt.Rpore) then
            diel= dielv -(dielv -diel_m)*(z1-abs(zg))/dbb
!
! 5) in-between
          else if(rr.gt.r1) then
            diel= ( dielv -(dielv -diel_m)*(z1-abs(zg))/dbb     &
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
      subroutine sprmul (xg,yg,zg,ag,fpx,fpy,fpz,np)
!------------------------------------------------------------
!  >>> Inextensible chains + bond angle rigidity (ggg).
!
      use, intrinsic :: iso_c_binding 
      implicit none
!
      include   'paramWatQa.h'
!    
      real(C_DOUBLE),dimension(npqr0) :: xg,yg,zg,ag
!
      real(C_DOUBLE),dimension(np0) :: fpx,fpy,fpz
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
      real(C_DOUBLE) bond_ps,bond_ss,a_phos,a_sugar,a_baseA,   &
                     a_baseG,a_baseC,a_baseT
      common/pbase/ np00,nnb,ist1(100),ist2(100)
      common/pbond/ bond_ps,bond_ss,a_phos,a_sugar,a_baseA,   &
                    a_baseG,a_baseC,a_baseT 
!
      integer(C_INT) io_pe
      common/sub_proc/ io_pe
!
      integer(C_INT) iwrt1,iwrt2,iwrt3,iwrt4
      common/iotim/ iwrt1,iwrt2,iwrt3,iwrt4
!
!   fgm0= 1.5d0/2= 0.75d0 
      fgm0= 0.75d0
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
        write(11,700) e_elas/np
  700   format('# e_elas/np=',1pd12.3)
        close(11)
      end if
!
      return
      end subroutine sprmul 
!
!
!------------------------------------------------------------------
      subroutine charges (rho,xg,yg,zg,xa,ya,za,xe,ye,ze,ch,chg, &
                          A11,A12,A13,A21,A22,A23,A31,A32,A33,   &
                          xr,yr,zr,ql,gg,ipar,np,nq,nCLp,npqr,   &
                          npqr3,npqr5)
!------------------------------------------------------------------
!* Even meshes are always used.
      use, intrinsic :: iso_c_binding 
      implicit none
!
      include    'mpif.h'
      include    'paramWatQa.h'
!
      integer(C_INT) ipar,np,nq,nCLp,npqr,npqr3,npqr5 
!
      real(C_DOUBLE),dimension(0:mx-1,0:my-1,0:mz-1) :: rho,rho2 
      real(C_DOUBLE),dimension(0:npqr50-1) :: xa,ya,za,ch
      real(C_DOUBLE),dimension(0:npqr30-1) :: xe,ye,ze,chg
      real(C_DOUBLE),dimension(0:npqr0-1) ::  xg,yg,zg  
      real(C_DOUBLE),dimension(-1:1) :: ffx,ffy,ffz
!
      real(C_DOUBLE),dimension(0:ip0**3-1,0:npqr30-1) :: ql
      integer(C_INT),dimension(0:npqr30-1,0:2) :: gg
!
      real(C_DOUBLE),dimension(0:npqr0-1) ::  A11,A12,A13,A21,A22,A23,A31,A32,A33
      real(C_DOUBLE),dimension(0:npqr50-1) :: xr,yr,zr
!
      real(C_DOUBLE) xxa,yya,zza,xxb,yyb,zzb,vec1,xhh,yhh,zhh, &
                     xpoint,ypoint,zpoint,xxc,yyc,zzc,vec3,    &
                     dohL,doLsin,xg1,yg1,zg1
      real(C_DOUBLE) t1,t2,t3,d1,d3,xx,yy,zz
      integer(C_INT) i,j,k,m,m0,ik,xpos,ypos,zpos,xps,yps,zps
!                
      real(C_DOUBLE) xmax,ymax,zmax,xmin,ymin,zmin
      real(C_float)  phi,tht,dtwr1,dtwr2,dtwr3
      common/parm3/  xmax,ymax,zmax,xmin,ymin,zmin
      common/parm4/  phi,tht,dtwr1,dtwr2,dtwr3
!
      integer(C_INT) pxr,pxc,pxl,pyr,pyc,pyl,pzr,pzc,pzl
      common/ptable/ pxr(-10:mx+10),pxc(-10:mx+10),pxl(-10:mx+10), &
                     pyr(-10:my+10),pyc(-10:my+10),pyl(-10:my+10), &
                     pzr(-10:mz+10),pzc(-10:mz+10),pzl(-10:mz+10)
!
      real(C_DOUBLE) gx,gy,gz,ghx,ghy,ghz,ghx2,ghy2,ghz2, &
                     ghxsq,ghysq,ghzsq,hxi,hyi,hzi
      integer(C_INT) ptx,pty,ptz
      common/xtable/ gx(0:mx),ghx(0:mx),ghx2(0:mx),ghxsq(0:mx), &
                     gy(0:my),ghy(0:my),ghy2(0:my),ghysq(0:my), &
                     gz(0:mz),ghz(0:mz),ghz2(0:mz),ghzsq(0:mz), &
                     hxi,hyi,hzi
      common/xtabl2/ ptx(-199:3000),pty(-199:3000),ptz(-100:4000)
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
! Ions only 
      do i= 0,nCLp-1
      xe(i)= xg(i)
      ye(i)= yg(i)
      ze(i)= zg(i)
!
      chg(i)= ch(i)
      end do
!
! Convert from 5-point to 3-point molecules
!
      do i= nCLp,npqr3-1
      if(mod(i-nCLp,3).eq.0) then
        chg(i)= -1.0d0
      else
        chg(i)=  0.5d0
      end if
      end do
!
!  /init/ A11-A33 initialized, /charges/ A11-A33 changed
!
      j= nCLp-1
      k= nCLp-1
!
      do i= nCLp,npqr5-1,5
      j= j +1
      k= k +3
!
      xa(i  )= xg(j) +A11(j)*xr(i) +A21(j)*yr(i) +A31(j)*zr(i)
      ya(i  )= yg(j) +A12(j)*xr(i) +A22(j)*yr(i) +A32(j)*zr(i)
      za(i  )= zg(j) +A13(j)*xr(i) +A23(j)*yr(i) +A33(j)*zr(i)
!  
      xa(i+1)= xg(j) +A11(j)*xr(i+1) +A21(j)*yr(i+1) +A31(j)*zr(i+1)
      ya(i+1)= yg(j) +A12(j)*xr(i+1) +A22(j)*yr(i+1) +A32(j)*zr(i+1)
      za(i+1)= zg(j) +A13(j)*xr(i+1) +A23(j)*yr(i+1) +A33(j)*zr(i+1)
! 
      xa(i+2)= xg(j) +A11(j)*xr(i+2) +A21(j)*yr(i+2) +A31(j)*zr(i+2)
      ya(i+2)= yg(j) +A12(j)*xr(i+2) +A22(j)*yr(i+2) +A32(j)*zr(i+2)
      za(i+2)= zg(j) +A13(j)*xr(i+2) +A23(j)*yr(i+2) +A33(j)*zr(i+2)
! 
! Use three-point O-H-H
!     xg1= (16*xa(i) +xa(i+1) +xa(i+2))/18.d0 
!     yg1= (16*ya(i) +ya(i+1) +ya(i+2))/18.d0 
!     zg1= (16*za(i) +za(i+1) +za(i+2))/18.d0 
!
      xe(k)  = xa(i)
      ye(k)  = ya(i)
      ze(k)  = za(i)
!
      xe(k+1)= (xa(i+1) +xa(i+3))/2.d0
      ye(k+1)= (ya(i+1) +ya(i+3))/2.d0
      ze(k+1)= (za(i+1) +za(i+3))/2.d0
!
      xe(k+2)= (xa(i+2) +xa(i+4))/2.d0
      ye(k+2)= (ya(i+2) +ya(i+4))/2.d0
      ze(k+2)= (za(i+2) +za(i+4))/2.d0
      end do
!
!  ip0= 3, Non-periodic case in (-l/2,l/2).
!  definition of d1, gg() are modified.
!
!* Need all definition for ql() in /moldyn/
!
      do i= 0,npqr3-1 
!
      d1  = (xe(i) -xmin)*hxi
      xps = ptx(int(d1 + 0.5d0))
      gg(i,0) = pxc(xps)
!
      xx = (xe(i) -gx(xps))/ghx(xps)
      ffx(-1)= 0.50d0*(0.5d0-xx)*(0.5d0-xx)
      ffx( 0)= 0.75d0-xx*xx
      ffx( 1)= 0.50d0*(0.5d0+xx)*(0.5d0+xx)
!*      
      d1  = (ye(i) -ymin)*hyi
      yps = pty(int(d1 + 0.5d0))
      gg(i,1) = pyc(yps)
!
      yy = (ye(i) -gy(yps))/ghy(yps) 
      ffy(-1)= 0.50d0*(0.5d0-yy)*(0.5d0-yy)
      ffy( 0)= 0.75d0-yy*yy
      ffy( 1)= 0.50d0*(0.5d0+yy)*(0.5d0+yy)
!*      
      d3  = (ze(i) -zmin)*hzi
      zps = ptz(int(d3 + 0.5d0))
      gg(i,2) = pzc(zps)
!
      zz = (ze(i) -gz(zps))/ghz(zps) 
      ffz(-1)= 0.50d0*(0.5d0-zz)*(0.5d0-zz)
      ffz( 0)= 0.75d0-zz*zz
      ffz( 1)= 0.50d0*(0.5d0+zz)*(0.5d0+zz)
!
      m0= -1 
!
      do j = -1,1            ! must be centered at 0
      do k = -1,1
      do m = -1,1
      t1 = chg(i) *ffx(j)
      t2 = t1 *    ffy(k)
      t3 = t2 *    ffz(m)
!    
      m0= m0 + 1
      ql(m0,i) = t3          ! assignment factor.
      end do
      end do
      end do
!*
      end do
!      
!* The index of this loop must be partial (no double count).
!  round-robin partitioning is used for nclp particles.
!
!           ++++++ +++++++
      do i= ipar-1,npqr3-1,num_proc 
      m0= -1
!
      do j = -1,1
      do k = -1,1
      do m = -1,1
      xpos = pxc(gg(i,0) + j)
      ypos = pyc(gg(i,1) + k)
      zpos = pzc(gg(i,2) + m)
!    
      m0= m0 +1
      rho(xpos,ypos,zpos) = rho(xpos,ypos,zpos) + ql(m0,i)
      end do
      end do
      end do
!*
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
      subroutine Labels_P
!----------------------------------------------------------
!  Periodic (x,y) and bounded (z) boundaries
!   Indices of the neighboring (27) sub-boxes.
!
      use, intrinsic :: iso_c_binding 
      implicit none
!
      include  'paramWatQa.h'
      integer(C_INT) i,j,k,mm,nn,n,ip,im,jp,jm,kp,km,icmax
!
      integer(C_INT) isize,isizeZ,isize2,isize4,nc3
      integer(C_INT) ibind
      parameter (isize=12,isizeZ=25,nc3=isize**2*isizeZ, &
                 isize2=isize*isize,isize4=isize2+isize)
      common /boxind/ ibind(27,nc3)
!
!     icmax= nc3 + 1
!
      do nn= 1,nc3
      do mm= 1,27
      ibind(mm,nn)= 0
      end do
      end do
!
!*  Extension boxes
!
      do k = 1,isizeZ
      kp = k + 1
      km = k - 1
      if (k.eq.isizeZ) goto 1 
      if (k.eq.1) goto 1 
!+
      do j = 1,isize
      jp = j + 1
      jm = j - 1
      if(jp.gt.isize) jp= 1  !<-- periodic
      if(jm.lt.1) jm= isize 
!+
      do i = 1,isize
      ip = i + 1
      im = i - 1
      if(ip.gt.isize) ip= 1  !<-- periodic
      if(im.lt.1) im= isize
!     if (i.eq.isize) goto 1 !ip = isize !icmax ! 1     ! must be excluded 
!     if (i.eq.1) goto 1 !im=1         !icmax ! isize ! - otherwise, double counted
!
      n = i + isize*(j-1) + isize2*(k-1)
!
      ibind( 1,n) = i  +  isize*( j-1) +  k*isize2 - isize4 
      ibind( 2,n) = ip +  isize*( j-1) +  k*isize2 - isize4
      ibind( 3,n) = ip +  isize*(jm-1) +  k*isize2 - isize4
      ibind( 4,n) = ip +  isize*(jp-1) +  k*isize2 - isize4
      ibind( 5,n) = i  +  isize*(jp-1) +  k*isize2 - isize4 
!
      ibind( 6,n) = i  +  isize*( j-1) + kp*isize2 - isize4 
      ibind( 7,n) = im +  isize*( j-1) + kp*isize2 - isize4
      ibind( 8,n) = ip +  isize*( j-1) + kp*isize2 - isize4
      ibind( 9,n) = i  +  isize*(jm-1) + kp*isize2 - isize4
      ibind(10,n) = im +  isize*(jm-1) + kp*isize2 - isize4
      ibind(11,n) = ip +  isize*(jm-1) + kp*isize2 - isize4 
      ibind(12,n) = i  +  isize*(jp-1) + kp*isize2 - isize4
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
!
      end do
      end do
!
    1 continue
      end do
!
      return
      end subroutine Labels_P
!
!
!-------------------------------------------------------------------
      subroutine reflect_endpl_P (xa,ya,za,xg,yg,zg,vx,vy,vz,ag,  &
                                  np,nq,npqr,npqr5)
!-------------------------------------------------------------------
      use, intrinsic :: iso_c_binding 
      implicit none
!
      include   'paramWatQa.h'
!
      real(C_DOUBLE),dimension(npqr50) :: xa,ya,za
      real(C_DOUBLE),dimension(npqr0) ::  xg,yg,zg,vx,vy,vz,ag
      integer(C_INT) np,nq,npqr,npqr5
!
      integer(C_INT) isize,isizeZ,isize2,isize4,nc3
      parameter (isize=12,isizeZ=25,nc3=isize**2*isizeZ, &
                 isize2=isize*isize,isize4=isize2+isize)
!
      real(C_DOUBLE) doh,doL,dohL,dohcos,dohsin,doLcos,doLsin,phwat
      common/unitHL/ doh,doL,dohL,dohcos,dohsin,doLcos,doLsin,phwat
!
      integer(C_INT),dimension(npqr0) :: imark
      integer(C_INT) ifqq,nCLp,i,j 
!
      real(C_DOUBLE) xmax3,ymax3,zmax3,xmin3,ymin3,zmin3,    &
                     zm1_t,zm2_t,zm1_b,zm2_b,dr,xr1,xr2,xr3, &
                     yr1,yr2,yr3,zr1,zr2,zr3,xg1,yg1,zg1,    &
                     x1,y1,z1,pe01,pe11,pe21,pe31,           &
                     xxa,yya,zza,xxb,yyb,zzb,vec1,xhh,yhh,zhh, &
                     xpoint,ypoint,zpoint,xxc,yyc,zzc,vec3
!
      integer(C_INT) it,is
      common/parm1/ it,is
!
      integer(C_INT) io_pe
      common/sub_proc/ io_pe
!
      real(C_DOUBLE) t8,cptot
      common/headr3/ t8,cptot
!
      real(C_DOUBLE) xmax,ymax,zmax,xmin,ymin,zmin,      &
                     xleng,yleng,zleng,Hpore2,rmax,      & 
                     qfrac,Rpore,Hpore,Zci,Zcp,Zcn,      &
                     rr,rr2,vpara,zpor1,zpor2,ddz 
!
      real(C_float)  phi,tht,dtwr1,dtwr2,dtwr3
      common/parm3/  xmax,ymax,zmax,xmin,ymin,zmin
      common/parm4/  phi,tht,dtwr1,dtwr2,dtwr3
      common/parm8/  xleng,yleng,zleng
      common/cntion/ qfrac,Rpore,Hpore,Zci,Zcp,Zcn,ifqq
!
      real(C_DOUBLE) pi,dt,dth,axi,Gamma,rbmax,vth,tmax
      common/parm2/  pi,dt,dth,axi,Gamma,rbmax,vth,tmax
!
      real(C_float)  t,t00,xp_leng
      common/headr2/ t,t00,xp_leng
!
!-----------------------------------------
!* Center of the system is (0,0,0)
!-----------------------------------------
!* Particle box is smaller than the field box.
!* grids i= 0, ..., (mx-1)
!
      xmax3= xmax -0.01d0
      ymax3= ymax -0.01d0
      zmax3= zmax -0.5d0*zleng/isizeZ
!
      xmin3= xmin +0.01d0
      ymin3= ymin +0.01d0
      zmin3= zmin +0.5d0*zleng/isizeZ     !<-- (-0.5, 0.5)*zleng/isizeZ
!
      Hpore2 = 0.5d0*Hpore 
!
!  In end-plates
      do i= 1,npqr
      imark(i)= 0
!
! 1) Outside the pore region, hitting the box sides.
! 
      if(abs(zg(i)).gt.Hpore2) then  ! Outside region
!
!* The X,Y directions are periodic 
!
        if(xg(i).lt.xmin) then
          xg(i)= xg(i) +xleng
          imark(i)= 1
!
        else if(xg(i).gt.xmax) then
          xg(i)= xg(i) -xleng
          imark(i)= 1
        end if
!
        if(yg(i).lt.ymin) then
          yg(i)= yg(i) +yleng
          imark(i)= 1
!
        else if(yg(i).gt.ymax) then
          yg(i)= yg(i) -yleng
          imark(i)= 1
        end if
!
!* The Z boundaries are closed 
!
        if(zg(i).lt.zmin3) then
          zg(i) = 2.d0*zmin3 -zg(i)
          vz(i) = -vz(i)
          imark(i)= 1
!
        else if(zg(i).gt.zmax3) then
          zg(i) = 2.d0*zmax3 -zg(i)
          vz(i) = -vz(i)
          imark(i)= 1
        end if
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
          imark(i)= 1
!
!         else
!           if(zg(i).lt.0.d0) then    ! bottom side
!             zpor1 = -Hpore2 
!             zg(i)= 2.d0*zpor1 -zg(i)
!             vz(i)= -vz(i)
!             imark(i)= 1
!           else
!             zpor2 = Hpore2
!             zg(i)= 2.d0*zpor2 -zg(i) ! top side
!             vz(i)= -vz(i)
!             imark(i)= 1
!           end if
!         end if
       end if
!*
      end if
      end do
!
!  xg(j)= 2.d0*xmin1-xg(j), and similar with xa(i)-xa(i+2)
!  xa(), ya() and za() are processed with xg(), yg() and zg().
!
      j= np+nq
      do i= np+nq+1,npqr5,5
      j= j +1
! 
      if(imark(j).eq.0) go to 100
!
      if(xg(j).lt.xmin) then
        xa(i  )= xa(i  ) +xleng
        xa(i+1)= xa(i+1) +xleng
        xa(i+2)= xa(i+2) +xleng
!
      else if(xg(j).gt.xmax) then
        xa(i  )= xa(i  ) -xleng 
        xa(i+1)= xa(i+1) -xleng
        xa(i+2)= xa(i+2) -xleng
      end if
!
      if(yg(j).lt.ymin) then
        ya(i  )= ya(i  ) +yleng
        ya(i+1)= ya(i+1) +yleng
        ya(i+2)= ya(i+2) +yleng
!
      else if(yg(j).gt.ymax) then
        ya(i  )= ya(i  ) -yleng
        ya(i+1)= ya(i+1) -yleng
        ya(i+2)= ya(i+2) -yleng
      end if
!
      if(zg(j).lt.zmin3) then
        za(i  )= 2.d0*zmin3 -za(i  )
        za(i+1)= 2.d0*zmin3 -za(i+1)
        za(i+2)= 2.d0*zmin3 -za(i+2)
!
      else if(zg(j).gt.zmax3) then
        za(i  )= 2.d0*zmax3 -za(i  )
        za(i+1)= 2.d0*zmax3 -za(i+1)
        za(i+2)= 2.d0*zmax3 -za(i+2)
      end if
!
  100 continue
      end do
!
      return
      end subroutine reflect_endpl_P 
!
!
!------------------------------------------------------------------
      subroutine READ_CONF (np,nq,ifbase) 
!------------------------------------------------------------------
      use, intrinsic :: iso_c_binding 
      implicit none
!
      include  'paramWatQa.h'
!
      integer(C_INT) np,nq,ifbase
      integer(C_INT) io_pe
      common/sub_proc/ io_pe
!
      character praefix_name*6,text1*40
      integer(C_INT) nsg,nseg,mpc,ladabst,n_p,n_lp,v_g,n_smol, &
                     v_sp,v_sm,seed,Nzions,npartcls,     &
                     SubBox  
!
      real(C_DOUBLE) rcut_clf,rcutlj,r_fene,     &
                     skin,skin2,k_fene,Bjerrum,r_fene2
      common /poly/    n_p,mpc,ladabst
      common /fenepara/ k_fene,r_fene2
      common /elsta/ Bjerrum
      common /einles/ skin,skin2
!
      common /confdatai/ n_lp,v_g,n_smol,v_sp,v_sm,seed
      common /confdatar/ rcut_clf,rcutlj,r_fene
      common /confdatas/ praefix_name
!----------------------------------------------------------------
      real(C_DOUBLE) pi,dt,dth,axi,Gamma,rbmax,vth,tmax,tmin, &
                     xmax,ymax,zmax,xmin,ymin,zmin,       &
                     xleng,yleng,zleng,qfrac,Rpore,Hpore, &
                     Zci,Zcp,Zcn,acount,acoion
      real(C_float)  phi,tht,dtwr1,dtwr2,dtwr3
      integer(C_INT) ifqq
!
      common/parm2/  pi,dt,dth,axi,Gamma,rbmax,vth,tmax
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
      real(C_DOUBLE) Vtop,Vbot,Vtop0,Vbot0,t_recyc
      integer(C_INT) fixedbc_x,fixedbc_y,itermax,filtx,filty,filtz
      common/espot/  Vtop,Vbot,Vtop0,Vbot0,fixedbc_x,fixedbc_y, &
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
      real(C_DOUBLE) diel2,dielpr,daa,dbb
      common/dielec/ diel2,dielpr,daa,dbb
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
        write(11,*)
        write(11,*) 'READ_CONF: parameter read... start'
        close(11)
      end if
!*
      read(08,'(a40,a6)') text1,praefix_name ! String der Simulationserkennung
!
      read(08,'(a40,i12)') text1,ifqq       !*Non-neutral Pl (if=1)   1 
      ifqq = 1 !!
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
      if(kstart.eq.0) then
        np  = n_p * n_lp
        nseg= n_p
!
        npartcls = np + Nzions
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
!noporWat.f03
      xmin= -0.5d0*xleng
      ymin= -0.5d0*yleng
      zmin= -0.5d0*zleng
!  -----------------------
!
      read(08,'(a40,f20.0)') text1,Vtop0    ! Potential of the top plate 0.
      read(08,'(a40,f20.0)') text1,Vbot0    ! Potential of the bot plate -0.
!
      if(Vtop0.lt.Vbot0) then
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
        write(11,*) 'dielpr,diel2=',dielpr,diel2
        write(11,*) 'Rmac,rod_leng (Ang)=',Rmac,rod_leng
        write(11,*) 'acount,acoion (Ang)=',acount,acoion
        write(11,*) 'rcut_clf=',rcut_clf
        write(11,*) 'Bjerrum (Ang)=',Bjerrum
!
        write(11,*) 'READ_CONF: parameter read... end'
        write(11,*)
!
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
      include    'paramWatQa.h'
      integer(C_INT) ifbase,np,nq,npqr
      integer(C_INT) ifqq,mpc,ladabst,n_p,n_lp,v_g,n_smol, &
                     v_sp,v_sm,seed
      real(C_DOUBLE) rcut_clf,rcutlj,r_fene,   &
                     k_fene,skin,Bjerrum,r_fene2,skin2
      character   praefix_name*6
      logical     ende
!
      common /poly/     n_p,mpc,ladabst
      common /fenepara/ k_fene,r_fene2
      common /elsta/ Bjerrum
      common /einles/ skin,skin2
!
      common /confdatai/ n_lp,v_g,n_smol,v_sp,v_sm,seed
      common /confdatar/ rcut_clf,rcutlj,r_fene
      common /confdatas/ praefix_name
!
      real(C_DOUBLE) qfrac,Rpore,Hpore,Zci,Zcp,Zcn
      common/cntion/ qfrac,Rpore,Hpore,Zci,Zcp,Zcn,ifqq
!
      real(C_DOUBLE) diel2,dielpr,daa,dbb
      common/dielec/ diel2,dielpr,daa,dbb
!
      real(C_DOUBLE) Vtop,Vbot,Vtop0,Vbot0
      integer(C_INT) fixedbc_x,fixedbc_y,itermax, &
                     filtx,filty,filtz,ipr
      common/espot/  Vtop,Vbot,Vtop0,Vbot0,fixedbc_x,fixedbc_y, &
                     itermax,filtx,filty,filtz
!----------------------------------------------------------------
      real(C_DOUBLE) pi,dt,dth,axi,Gamma,rbmax,vth,tmax, &
                     xmax,ymax,zmax,xmin,ymin,zmin,      &
                     xleng,yleng,zleng,acount,acoion   
      common/parm8/  xleng,yleng,zleng
      real(C_float)  phi,tht,dtwr1,dtwr2,dtwr3
!
      common/parm2/ pi,dt,dth,axi,Gamma,rbmax,vth,tmax
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
      write(07,'("# praefix-string................: ",a6)')praefix_name
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
      include   'paramWatQa.h'
!
      integer(C_INT) n_g,n_sp,n_sm,n_si,mpc,n_p,ladabst, &
                     n_lp,v_g,n_smol,v_sp,v_sm,seed,     &
                     verletfix,qlj,i,iverg,i_steps,      &
                     measstep,confstep,nCLp
      real(C_DOUBLE) qwert,rcutlj,rcut_clf,skin,         &
                     rcps,fmesh,skin2,r_fene2,rcutlj2,   &
                     k_fene,r_fene,Bjerrum,wupi,         &
                     e_c_s,e_c_pme,e_c_r,e_lj,e_elas
!-----------
!                    **** 
      common/ewald2/ nCLp   ! ,ip0 paramWatQa.h
!
      common /inidati/ n_g,n_sp,n_sm,n_si
      common /inidatr/ qwert,fmesh
      common /chargelj/ qlj(npqr50)
      common /poly/ n_p,mpc,ladabst
      common /einles/ skin,skin2
      common /fenepara/ k_fene,r_fene2

      common /confdatai/ n_lp,v_g,n_smol,v_sp,v_sm,seed
      common /confdatar/ rcut_clf,rcutlj,r_fene
      common /energy/ e_c_s,e_c_pme,e_c_r,e_lj,e_elas
      common /elsta/ Bjerrum
!----------------------------------------------------------------
      real(C_DOUBLE) pi,dt,dth,axi,Gamma,rbmax,vth,tmax, &
                    xmax,ymax,zmax,xmin,ymin,zmin,       &
                    xleng,yleng,zleng,acount,acoion   
      real(C_float) phi,tht,dtwr1,dtwr2,dtwr3
      common/parm2/ pi,dt,dth,axi,Gamma,rbmax,vth,tmax
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
!     rcut_clf2  = rcut_clf* rcut_clf
      r_fene2   = r_fene * r_fene
      rcutlj2   = rcutlj * rcutlj
!     rcps2     = rcps   * rcps
  
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
      include    'paramWatQa.h'
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
!--------------------------------------------------------------------
      subroutine init_P (xg,yg,zg,vx,vy,vz,xa,ya,za,ch,am,ag,ep,amm, &
                         A11,A12,A13,A21,A22,A23,A31,A32,A33,        &
                         e0,e1,e2,e3,ifrgrod,ifrodco,ifbase,ipar,    &
                         np,nq,nCLp,nr,npqr,npqr3,npqr5,             &
                         npqr5_0,npqr5_1,npqr5_2,npqr_0,npqr_1,npqr_2)
!--------------------------------------------------------------------
!  Periodic (x,y) and bounded (z) boundaries
!
      use, intrinsic :: iso_c_binding 
      implicit none
!
      include   'mpif.h'
      include   'paramWatQa.h'
!
      real(C_DOUBLE),dimension(npqr50) :: xa,ya,za,ch,am
      real(C_DOUBLE),dimension(npqr0) ::  xg,yg,zg,vx,vy,vz,amm,ag,ep
!
      real(C_DOUBLE),dimension(npqr0) ::  A11,A12,A13,A21,A22,A23,A31,A32,A33
      real(C_DOUBLE),dimension(npqr0) ::  e0,e1,e2,e3
!
      integer(C_INT) ifrgrod,ifrodco,ifbase,ipar,np,nq,nCLp,nr,  &
                     npqr,npqr3,npqr5
!
      integer(C_INT) io_pe
      common/sub_proc/ io_pe
! 
      integer(C_INT) it,is,ifqq,nsg,nseg,i,j,k,l,nps,k1,       &
                     Nzi,Nci,Ncc_pore,Ncc_wide,n0,n1,n2,n3,n4, &
                     nn3,nn4,ntried,i1,i2,ll,kl,km,kn,i0,      &
                     natom1,npqr100,npqr200,npqr5_0,npqr5_1,   &
                     npqr5_2,npqr_0,npqr_1,npqr_2
      character(len=4) line82*82,line79*79,line5*5,line4*4,    &
                       dummy1*1
      character(len=4),dimension(npqr50) :: tip5
!
      real(C_DOUBLE) pi,dt,dth,axi,Gamma,rbmax,vth,tmax, &
                     xmax,ymax,zmax,xmin,ymin,zmin,      &
                     xleng,yleng,zleng,                  &
                     ddz,KJ,KCal,mol,kbT,ranff,ch8,abin, & 
                     dgaus2,vmax1,vmax2,a_cos,a_sin
      real(C_float)  phi,tht,dtwr1,dtwr2,dtwr3
!
      common/parm1/ it,is
      common/parm2/ pi,dt,dth,axi,Gamma,rbmax,vth,tmax
      common/parm3/ xmax,ymax,zmax,xmin,ymin,zmin
      common/parm4/ phi,tht,dtwr1,dtwr2,dtwr3
      common/parm8/ xleng,yleng,zleng
      common/psegm/ nsg(30),nseg
!
      real(C_DOUBLE) qfrac,Rpore,Hpore,Zci,Zcp,Zcn,Hpore2,   &
                     wwat,awat,awat_diam,acount,acoion,ag0,am0
      real(C_DOUBLE) doh,doL,dohL,dohcos,dohsin,doLcos,doLsin,phwat
! 
      common/cntion/ qfrac,Rpore,Hpore,Zci,Zcp,Zcn,ifqq
      common/waterm/ wwat,awat
      common/ionsiz/ acount,acoion   
      common/unitHL/ doh,doL,dohL,dohcos,dohsin,doLcos,doLsin,phwat
!
      real(C_DOUBLE) pi2,q_PE,aLJ,cci,dnint,rr,x0,y0,z0, &
                     rcut_clf,rcutlj,r_fene
      common/shrgf/ aLJ
      common/confdatar/ rcut_clf,rcutlj,r_fene
!
      real(C_DOUBLE) fv,vv,dv
      common/gaus1/ fv(101),vv,dv
!
      integer(C_INT) pxr,pxc,pxl,pyr,pyc,pyl,pzr,pzc,pzl
      common/ptable/ pxr(-10:mx+10),pxc(-10:mx+10),pxl(-10:mx+10), &
                     pyr(-10:my+10),pyc(-10:my+10),pyl(-10:my+10), &
                     pzr(-10:mz+10),pzc(-10:mz+10),pzl(-10:mz+10)
!
      real(C_DOUBLE) gx,gy,gz,ghx,ghy,ghz,ghx2,ghy2,ghz2,  &
                     ghxsq,ghysq,ghzsq,hxi,hyi,hzi,        &
                     ghx8(0:mx),ghy8(0:my),ghz8(0:mz),     &
                     ghx0,ghy0,ghz0,hx8,hy8,hz8,hz0,hhz0,  &
                     zz0,zm,hx0,hhx0,xx0,xm,hy0,hhy0,wg
      common/xtable/ gx(0:mx),ghx(0:mx),ghx2(0:mx),ghxsq(0:mx), &
                     gy(0:my),ghy(0:my),ghy2(0:my),ghysq(0:my), &
                     gz(0:mz),ghz(0:mz),ghz2(0:mz),ghzsq(0:mz), &
                     hxi,hyi,hzi
!   Ranges: ptx: -199->i1+200, ptz: -100->k1+100
      integer(C_INT) ptx,pty,ptz,icentx,icentz,ismall,k0
      common/xtabl2/ ptx(-199:3000),pty(-199:3000),ptz(-100:4000)
!
      real(C_DOUBLE) xx,yy,zz,diel,xgc,ygc,zgc,vxg,vyg,vzg,    &
                     rod_leng,rod_len(100),Rmac,Rhelix,        &
                     th,ph,ph3,ss,airc,rxg(np0),ryg(np0),rzg(np0), &
                     q1,q2,q3,q4,angx,angy,angz,ipx,ipy,ipz,   &
                     ipz1,xyz1,xyz3,rr0,                       &
                     xxa,yya,zza,xxb,yyb,zzb,vec1,xhh,yhh,zhh, &
                     xpoint,ypoint,zpoint,xxc,yyc,zzc,vec3,    &
                     xnum,ynum,znum,xcent,ycent,zcent,         &
                     z_top,z_bottom,z_bot  
!
      integer(C_INT) n_rodp,jj,np00,nnb,ist1,ist2
      common/pbase/  np00,nnb,ist1(100),ist2(100)
!
      real(C_DOUBLE) bond_ps,bond_ss,a_phos,a_sugar,a_baseA,   &
                     a_baseG,a_baseC,a_baseT,ww1,ww2
      common/pbond/  bond_ps,bond_ss,a_phos,a_sugar,a_baseA,   &
                     a_baseG,a_baseC,a_baseT 
!
      common/macroion/ xgc,ygc,zgc,vxg,vyg,vzg,rod_leng,Rmac,  &
                       Rhelix,n_rodp
      real(C_DOUBLE),dimension(np0) :: &
                                dox,doy,doz,doxc,doyc,dozc
      common/macroin1/ dox,doy,doz,doxc,doyc,dozc
      common/macroiov/ q1,q2,q3,q4,angx,angy,angz,ipx,ipy,ipz
!
      real(C_DOUBLE) Vtop,Vbot,Vtop0,Vbot0
      integer(C_INT) fixedbc_x,fixedbc_y,itermax,filtx,filty,filtz
      common/espot/  Vtop,Vbot,Vtop0,Vbot0,fixedbc_x,fixedbc_y, &
                     itermax,filtx,filty,filtz
!
      real(C_DOUBLE) epsLJ,eps_K,eps_Cl,rhow
      real(C_DOUBLE) q_O,q_H,q_L,masso,massh,epslj_w,epslj_A,epslj_B
      integer(C_INT) ifLJ,Ncc_1,Ncc_2,Ncc_3
      common/unitOH/ q_O,q_H,q_L,masso,massh,epslj_w,epslj_A,epslj_B
!
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
!  Periodic
!
      do i= 0,mx-1
      pxr(i)= i +1
      pxc(i)= i 
      pxl(i)= i -1
      end do
!
      pxc(-3)= mx-3
      pxc(-2)= mx-2
      pxc(-1)= mx-1
!
      pxc(mx  )= 0
      pxc(mx+1)= 1
      pxc(mx+2)= 2
!
      pxr(mx-1)=  0
      pxl(0)   = mx-1 
!
!
      do j= 0,my-1
      pyr(j)= j +1
      pyc(j)= j 
      pyl(j)= j -1
      end do
!
      pyc(-3)= my-3
      pyc(-2)= my-2
      pyc(-1)= my-1
!
      pyc(mx  )= 0
      pyc(mx+1)= 1
      pyc(mx+2)= 2
!
      pyr(my-1)=  0
      pyl(0)   = my-1 
!
!  Non-periodic
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
!
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
      if(kstart.eq.0) then
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
      end if
!
!--------------------------------
!* Find the nearest grid.
!--------------------------------
! The X direction
!
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
! Ranges: ptx: -199->i1+200, ptz: -100->k1+100
!
      do i= 1,200
      ptx(i1+i)= ptx(i)         ! 1-> i1+1, 200-> i1+200
!
      if(kstart.eq.0) then
        if(io_pe.eq.1) then
          OPEN (unit=11,file=praefixc//'.06'//suffix2,             &
                status='unknown',position='append',form='formatted')
          write(11,751) i,i1+i,i
  751     format('i: i1+i,i=',i5,2i8)
          close(11)
        end if
      end if
      end do
!
!
      do i= 1,199
      ptx(i-200)= ptx(i+i1-200) ! i+i1-200-> -199, i1-> -1
!
      if(kstart.eq.0) then
        if(io_pe.eq.1) then
          OPEN (unit=11,file=praefixc//'.06'//suffix2,             &
                status='unknown',position='append',form='formatted')
          write(11,753) i,i-200,i+i1-200
  753     format('i: i-200,i+i1-200=',i5,2i8)
          close(11)
        end if
        end if
      end do
!
! Y: from j= -10 to j= i1+10
!
      hhy0 = 0.3d0 * hy0
!
      do j= -199,i1+200
      pty(j)= ptx(j)
      end do
!
! -----------------------------------------
!  The Z direction
!
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
      do k= 1,100
      ptz(k1+k)= mz-1
      ptz(-k  )= 0
!
      if(kstart.eq.0) then
        if(io_pe.eq.1) then
          OPEN (unit=11,file=praefixc//'.06'//suffix2,             &
                status='unknown',position='append',form='formatted')
!
          write(11,757) k,k1+k,-k
  757     format('k: k1+k,-k=',i5,2i8)
          close(11)
        end if
       end if
      end do
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
        write(11,160) hxi,hyi,hzi
  160   format(/,' hxi, hyi, hzi=',3f12.7)
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
      a_baseA= 4.3d0   ! A, diameter
      a_baseG= 4.3d0   ! G 
      a_baseC= 3.3d0   ! C
      a_baseT= 3.3d0   ! T
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
! -----------------------------
      if(kstart.ge.1) then
        if(io_pe.eq.1) then
        OPEN (unit=11,file=praefixc//'.06'//suffix2,             &
              status='unknown',position='append',form='formatted')
!
        write(11,*) '# kstart >= 1, then subroutine /init/ retruns...'
        end if
!
        return
      end if
! -----------------------------
!
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
      end if 
!
! --------------------------------------------------------
!*  Additional monomers attached to sugar rings - A,G,T,C
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
      if(mod(i,2).eq.1) go to 250        ! skip (is not neutral)   
!
      nnb = nnb +1
      ch(np00 +nnb)= 0.d0                ! base: AGCT,heavy,LJ
!
      if(mod(nnb,4).eq.1) then
        ag0= a_baseA /2.d0     ! diameter/2
        am0= 134.d0/wwat
! 
      else if(mod(nnb,4).eq.2) then
        ag0= a_baseG /2.d0
        am0= 150.d0/wwat
! 
      else if(mod(nnb,4).eq.3) then
        ag0= a_baseC /2.d0
        am0= 125.d0/wwat
! 
      else if(mod(nnb,4).eq.0) then
        ag0= a_baseT /2.d0
        am0= 110.d0/wwat
      end if
!
      ag(np00 +nnb)= ag0
      am(np00 +nnb)= am0
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
!
  273 continue 
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
          write(11,275) i,k,ag(k),x0,y0,z0
  275     format('i,k,ag(k),x0-z0=',2i4,4f8.2)
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
      cci= cci +ch(i)    !<- 0
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
!* 3. Water particles
!-----------------------------------------
!  Must avoid macroions. One molecule in every 3.1 Angstrom
!
      if(io_pe.eq.1) then
        OPEN (unit=11,file=praefixc//'.06'//suffix2,             &
              status='unknown',position='append',form='formatted')
!
        write(11,*) 'xleng,yleng,zleng=',xleng,yleng,zleng
        close(11)
      end if
!
      do j= np+nq+1,npqr0
      ag(j)= 3.10d0 /2.d0 
      end do
!
      ll= np + nq   !!<- ions, rotation
      jj= np + nq   !<-- translation, quaternion, A11-A33
!  
      xcent= xleng/2.d0
      ycent= yleng/2.d0
      zcent= zleng/2.d0
!
      awat_diam= 3.10d0   !<- diameter
!
      do k= 17,34     !<-- (i,j,k) are used for cells 
      do j= 1,km
      do i= 1,kl
!
      x0= awat_diam*(i -0.5d0) -xcent  !<- middle center of water
      y0= awat_diam*(j -0.5d0) -ycent
      z0= awat_diam*(k -0.5d0) -zcent
!
      Hpore2= 0.5d0*Hpore 
!
      if(z0.lt.zmin .or. z0.gt.zmax) go to 530
!
!* Over the pore
!     if(abs(z0).gt.Hpore2) then   !<- z0 > Hpore/2
!
!       do l= 1,np+nq
!       rr= sqrt((x0-xg(l))**2 +(y0-yg(l))**2 +(z0-zg(l))**2)
!       if(rr.le.(ag(l)+awat)) go to 530
!       end do 
!
!       do l= np+nq+1,jj
!       rr= sqrt((x0-xg(l))**2 +(y0-yg(l))**2 +(z0-zg(l))**2)
!       if(rr.le.(ag(l)+awat)) go to 530
!       end do 
!
!* Within the pore
      if(abs(z0).ge.Hpore2) go to 530  !<- z0 < Hpore/2
!
      rr= sqrt(x0**2 +y0**2)
      if(rr.gt.Rpore-awat) go to 530 
!
      do l= 1,np+nq 
      rr= sqrt((x0-xg(l))**2 +(y0-yg(l))**2 +(z0-zg(l))**2)
      if(rr.le.(ag(l) +awat)) go to 530
      end do 
!
      do l= np+nq+1,jj
      rr= sqrt((x0-xg(l))**2 +(y0-yg(l))**2 +(z0-zg(l))**2)
      if(rr.le.(awat +awat)) go to 530
      end do 
!
! ll=(null)    np+nq+1(5)
! ll=np+nq+5   np+nq+10
! jj=np+nq+1   jj=np+nq+2
!
      ll= ll +5         !!<- ions, rotation
      jj= jj +1         !!<- translation, quaternion, A11-A33
!     +++++++++
!
      if(ll.gt.npqr50) then
        if(io_pe.eq.1) then
          OPEN (unit=11,file=praefixc//'.06'//suffix2,             &
                status='unknown',position='append',form='formatted')
!
          write(11,*) ' Water: the l is greater than npqr50.., ll=',ll
          close(11)
        end if
!
        stop
      end if
!
!
      phwat= 104.52d0       ! tip4p  <- /moldyn/
! 
      doh  = 0.9572d0 
      doL  = 0.70d0
!
      dohcos = doh * cos(pi*phwat/(2*180.d0))
      dohsin = doh * sin(pi*phwat/(2*180.d0))
!
      doLcos = doL * cos(pi*phwat/(2*180.d0))
      doLsin = doL * sin(pi*phwat/(2*180.d0))
      dohL = dohcos +doLcos
!
!  1:O, 2:H, 3:H, dummy 4:L, 5:L. Centered at x0,y0,z0
      xa(ll-4)= x0 +0.d0 
      ya(ll-4)= y0 +0.d0
      za(ll-4)= z0
!
      xa(ll-3)= x0 +dohsin 
      ya(ll-3)= y0 +dohcos
      za(ll-3)= z0 
!
      xa(ll-2)= x0 -dohsin
      ya(ll-2)= y0 +dohcos
      za(ll-2)= z0 
! 
!  Outside the triangle plane by virtual sites
!             H         H
      xxa= xa(ll-2) -xa(ll-3)
      yya= ya(ll-2) -ya(ll-3)
      zza= za(ll-2) -za(ll-3)
!             O          H        H
      xxb= xa(ll-4) -(xa(ll-3)+xa(ll-2))/2.d0
      yyb= ya(ll-4) -(ya(ll-3)+ya(ll-2))/2.d0
      zzb= za(ll-4) -(za(ll-3)+za(ll-2))/2.d0
      vec1= sqrt(xxb*xxb +yyb*yyb +zzb*zzb)
!
      xhh= (xa(ll-3)+xa(ll-2))/2.d0
      yhh= (ya(ll-3)+ya(ll-2))/2.d0
      zhh= (za(ll-3)+za(ll-2))/2.d0
!
      xpoint= xhh +dohL*xxb/vec1
      ypoint= yhh +dohL*yyb/vec1
      zpoint= zhh +dohL*zzb/vec1
!
      xxc= yya*zzb -zza*yyb
      yyc= zza*xxb -xxa*zzb
      zzc= xxa*yyb -yya*xxb
      vec3= sqrt(xxc*xxc +yyc*yyc +zzc*zzc)
!
      xa(ll-1)= xpoint +doLsin*xxc/vec3 
      ya(ll-1)= ypoint +doLsin*yyc/vec3
      za(ll-1)= zpoint +doLsin*zzc/vec3
!
      xa(ll)  = xpoint -doLsin*xxc/vec3 
      ya(ll)  = ypoint -doLsin*yyc/vec3
      za(ll)  = zpoint -doLsin*zzc/vec3
!
!   e0-e3 and A11-A33 
!
      xg(jj)= (16*xa(ll-4) +xa(ll-3) +xa(ll-2))/18.d0 
      yg(jj)= (16*ya(ll-4) +ya(ll-3) +ya(ll-2))/18.d0 
      zg(jj)= (16*za(ll-4) +za(ll-3) +za(ll-2))/18.d0 
!
      a_cos= cos(pi*phwat/(2*180.d0))
      a_sin= sin(pi*phwat/(2*180.d0))
!
      e0(jj)= a_cos !cos(theta/2)*cos((phii+psii)/2)
      e1(jj)= 0     !sin(theta/2)*cos((phii-psii)/2)
      e2(jj)= 0     !sin(theta/2)*sin((phii-psii)/2)
      e3(jj)= a_sin !cos(theta/2)*sin((phii+psii)/2)
!
!   Rotation matrix R(3 x 3) in terms of Quarternion
      A11(jj)= e0(jj)**2 +e1(jj)**2 -e2(jj)**2 -e3(jj)**2 
      A12(jj)= 2*(e1(jj)*e2(jj) +e0(jj)*e3(jj)) 
      A13(jj)= 2*(e1(jj)*e3(jj) -e0(jj)*e2(jj))
      A21(jj)= 2*(e1(jj)*e2(jj) -e0(jj)*e3(jj)) 
      A22(jj)= e0(jj)**2 -e1(jj)**2 +e2(jj)**2 -e3(jj)**2
      A23(jj)= 2*(e2(jj)*e3(jj) +e0(jj)*e1(jj))
      A31(jj)= 2*(e1(jj)*e3(jj) +e0(jj)*e2(jj))
      A32(jj)= 2*(e2(jj)*e3(jj) -e0(jj)*e1(jj))
      A33(jj)= e0(jj)**2 -e1(jj)**2 -e2(jj)**2 +e3(jj)**2
!
  530 continue     !<- Occupied in (i,j,k) 
      end do
      end do
      end do
!
!     +++++++++++
      npqr5_0= ll  !<- np+nq+(in the pore)
      npqr_0 = jj
!     +++++++++++
!
      if(io_pe.eq.1) then
        OPEN (unit=11,file=praefixc//'.06'//suffix2,             &
              status='unknown',position='append',form='formatted')
!
        write(11,*)
        write(11,*) 'npqr5_0,npqr_0=',npqr5_0,npqr_0
!
        do j= np+nq+1,npqr_0  ! jj
        write(11,430) j,xg(j),yg(j),zg(j)
  430   format('j,xg,yg,zg=',i7,3f8.2)
        end do
!
        close(11)
      end if
!
!**********************************************************************
!* Top and bottom parts except for the middle one                     *
!**********************************************************************
!  Periodic (x,y) and bounded (z) boundaries
!  Original coordinate at positive shapes > 0 
!
      z_top    =  27.90d0    !<- the pore top
      z_bottom = -27.90d0    !<- the pore bottom
!
      z_bot    = z_bottom -49.66d0
!                +++++++++++++++++  
!
      open (unit=12,file='1cx12128z.exyz',status='old',form='formatted')
!                      from 1cx12128z.exyz
!***
      read(12,710) line79
      read(12,713) natom1    ! number of atoms
      read(12,715) line4
  710 format(a79)
  713 format(i5)
  715 format(a4)
!
      npqr100= 36864         !<- 4-atoms by natom1, if 5-atoms is 46080
!
      do i= 1,npqr100,4
      ll= ll +5
!
      read(12,700) tip5(ll-4),xa(ll-4),ya(ll-4),za(ll-4)
      read(12,700) tip5(ll-3),xa(ll-3),ya(ll-3),za(ll-3)
      read(12,700) tip5(ll-2),xa(ll-2),ya(ll-2),za(ll-2)
      read(12,700) tip5(ll-1),xa(ll-1),ya(ll-1),za(ll-1)
  700 format(a4,f15.5,f15.5,f15.5)
      end do
! 
!     read(12,300) vtext1,vtextx1,vtexty1,vtextz1
!     read(12,300) vtext2,vtextx2,vtexty2,vtextz2
!     read(12,300) vtext3,vtextx3,vtexty3,vtextz3
!     close(12)
!
!     +++++++++++
      npqr5_1= ll            !<- include top shape
!     +++++++++++
!
!
      jj= npqr5_0 
      do i= npqr5_0+1,npqr5_1,5   !<- i > npqr5_0, jj > npqr5_0
      jj= jj +1                   !<-  jj= npqr5_0,npqr_1
!
      xa(i)= xa(i) -xcent    !<- half positive, half negative
      ya(i)= ya(i) -ycent
      za(i)= za(i) +z_top    !<- z > z_top
!
      xa(i+1)=  xa(i+1) -xcent
      ya(i+1)=  ya(i+1) -ycent
      za(i+1)=  za(i+1) +z_top
!
      xa(i+2)=  xa(i+2) -xcent
      ya(i+2)=  ya(i+2) -ycent
      za(i+2)=  za(i+2) +z_top
!
!  Outside the triangle plane by virtual sites
!             H         H
      xxa= xa(i+2) -xa(i+1)
      yya= ya(i+2) -ya(i+1)
      zza= za(i+2) -za(i+1)
!             O          H        H
      xxb= xa(i) -(xa(i+1)+xa(i+2))/2.d0
      yyb= ya(i) -(ya(i+1)+ya(i+2))/2.d0
      zzb= za(i) -(za(i+1)+za(i+2))/2.d0
      vec1= sqrt(xxb*xxb +yyb*yyb +zzb*zzb)
!
      xhh= (xa(i+1)+xa(i+2))/2.d0
      yhh= (ya(i+1)+ya(i+2))/2.d0
      zhh= (za(i+1)+za(i+2))/2.d0
!
      xpoint= xhh +dohL*xxb/vec1
      ypoint= yhh +dohL*yyb/vec1
      zpoint= zhh +dohL*zzb/vec1
!
      xxc= yya*zzb -zza*yyb
      yyc= zza*xxb -xxa*zzb
      zzc= xxa*yyb -yya*xxb
      vec3= sqrt(xxc*xxc +yyc*yyc +zzc*zzc)
!
      xa(i+3)= xpoint +doLsin*xxc/vec3 
      ya(i+3)= ypoint +doLsin*yyc/vec3
      za(i+3)= zpoint +doLsin*zzc/vec3
!
      xa(i+4)  = xpoint -doLsin*xxc/vec3 
      ya(i+4)  = ypoint -doLsin*yyc/vec3
      za(i+4)  = zpoint -doLsin*zzc/vec3
!
      xg(jj)= (16.d0*xa(i) +xa(i+1)+xa(i+2))/18.d0 
      yg(jj)= (16.d0*ya(i) +ya(i+1)+ya(i+2))/18.d0 
      zg(jj)= (16.d0*za(i) +za(i+1)+za(i+2))/18.d0 
      end do
!***
!     ++++++++++
      npqr_1= jj             !<- include the top
!     +++++++++++
!
!*
      rewind (12)
!     open (unit=12,file='1cx12128z.exyz',status='old',form='formatted')
!
      read(12,710) line79
      read(12,713) natom1  
      read(12,715) line4
!
!
      do i= 1,npqr100,4
      ll= ll +5
!
      read(12,700) tip5(ll-4),xa(ll-4),ya(ll-4),za(ll-4)
      read(12,700) tip5(ll-3),xa(ll-3),ya(ll-3),za(ll-3)
      read(12,700) tip5(ll-2),xa(ll-2),ya(ll-2),za(ll-2)
      read(12,700) tip5(ll-1),xa(ll-1),ya(ll-1),za(ll-1)
      end do
! 
!     read(12,300) vtext1,vtextx1,vtexty1,vtextz1
!     read(12,300) vtext2,vtextx2,vtexty2,vtextz2
!     read(12,300) vtext3,vtextx3,vtexty3,vtextz3
      close(12)
!
!     +++++++++++
      npqr5_2= ll            !<- include the bottom
!     +++++++++++
!
      if(io_pe.eq.1) then
        OPEN (unit=11,file=praefixc//'.06'//suffix2,             &
              status='unknown',position='append',form='formatted')
        write(11,992) npqr5_0,npqr5_1,npqr5_2
  992   format('992  npqr5_0,npqr5_1,npqr5_2=',3i8)
        close(11)
      end if
!
!
!   jj= npqr_1
      do i= npqr5_1+1,npqr5_2,5  !<- i > npqr5_1, jj > npqr5_1
      jj= jj +1                  !<- jj= npqr5_1,npqr_2
!
      xa(i)=  xa(i) -xcent
      ya(i)=  ya(i) -ycent
      za(i)=  za(i) +z_bot       !<- z < z_bottom
!
      xa(i+1)=  xa(i+1) -xcent
      ya(i+1)=  ya(i+1) -ycent
      za(i+1)=  za(i+1) +z_bot
!
      xa(i+2)=  xa(i+2) -xcent
      ya(i+2)=  ya(i+2) -ycent
      za(i+2)=  za(i+2) +z_bot
! 
!  Outside the triangle plane by virtual sites
!             H         H
      xxa= xa(i+2) -xa(i+1)
      yya= ya(i+2) -ya(i+1)
      zza= za(i+2) -za(i+1)
!             O          H        H
      xxb= xa(i) -(xa(i+1)+xa(i+2))/2.d0
      yyb= ya(i) -(ya(i+1)+ya(i+2))/2.d0
      zzb= za(i) -(za(i+1)+za(i+2))/2.d0
      vec1= sqrt(xxb*xxb +yyb*yyb +zzb*zzb)
!
      xhh= (xa(i+1)+xa(i+2))/2.d0
      yhh= (ya(i+1)+ya(i+2))/2.d0
      zhh= (za(i+1)+za(i+2))/2.d0
!
      xpoint= xhh +dohL*xxb/vec1
      ypoint= yhh +dohL*yyb/vec1
      zpoint= zhh +dohL*zzb/vec1
!
      xxc= yya*zzb -zza*yyb
      yyc= zza*xxb -xxa*zzb
      zzc= xxa*yyb -yya*xxb
      vec3= sqrt(xxc*xxc +yyc*yyc +zzc*zzc)
!
      xa(i+3)= xpoint +doLsin*xxc/vec3 
      ya(i+3)= ypoint +doLsin*yyc/vec3
      za(i+3)= zpoint +doLsin*zzc/vec3
!
      xa(i+4)  = xpoint -doLsin*xxc/vec3 
      ya(i+4)  = ypoint -doLsin*yyc/vec3
      za(i+4)  = zpoint -doLsin*zzc/vec3
!
      xg(jj)= (16.d0*xa(i) +xa(i+1)+xa(i+2))/18.d0 
      yg(jj)= (16.d0*ya(i) +ya(i+1)+ya(i+2))/18.d0 
      zg(jj)= (16.d0*za(i) +za(i+1)+za(i+2))/18.d0 
      end do
!***
!     ++++++++++
      npqr_2= jj
!     ++++++++++
!
      if(io_pe.eq.1) then
        OPEN (unit=11,file=praefixc//'.06'//suffix2,             &
              status='unknown',position='append',form='formatted')
        write(11,993) npqr_0,npqr_1,npqr_2
  993   format('993  npqr_0,npqr_1,npqr_2=',3i8)
        close(11)
      end if
!
!
!****************
!* Quaternion   *
!****************
!Command line: /usr/local/bin/analice2 1cx12128z.exyz -O OW -H HW[12] -w tip4p -f q
!@BOX3
!74.48265999999998 74.48265999999998 49.65509999999999
!@NX4A
!9216
!   0.0000    0.0000    0.0000     0.0000    0.9239   -0.3827    0.0000
!
      open (unit=12,file='1cx12128z.q',status='old',form='formatted')
!***
      read(12,*)
      read(12,730) line82
      read(12,732) line5
      read(12,735) xnum,ynum,znum
      read(12,732) line5
      read(12,737) natom1    !<- 1-atom molecules 
  730 format(a82)
  732 format(a5)
  735 format(f15.5,f15.5,f15.5)
  737 format(i4)
!
      npqr200= 9216          !<- natom1(i4) at top and bottom parts
!
      jj= npqr_0
      do i= 1,npqr200
      jj= jj +1              !<- npqr_0,npqr_1
!
      read(12,738) xnum,ynum,znum,dummy1,e0(jj),e1(jj),e2(jj),e3(jj)
  738 format(3f10.4,a1,4f10.4)         ! +++++++++++++++++++++++++++
      end do
!
!*
      rewind(12)
!
      read(12,*)
      read(12,730) line82
      read(12,732) line5
      read(12,735) xnum,ynum,znum
      read(12,732) line5
      read(12,737) natom1     !<- 1-point molecules
!
      jj= npqr_1
      do i= 1,npqr200
      jj= jj +1               !<- npqr_1,npqr_2
!
      read(12,738) xnum,ynum,znum,dummy1,e0(jj),e1(jj),e2(jj),e3(jj)
      end do                           ! +++++++++++++++++++++++++++ 
!
!***
      close(12)
!
!  Rotation matrix R(3 x 3) in terms of Quarternion
!
      j= npqr_0
      do i= npqr5_0+1,npqr5_2,5
      j= j +1
!
      A11(j)= e0(j)**2 +e1(j)**2 -e2(j)**2 -e3(j)**2 
      A12(j)= 2*(e1(j)*e2(j) +e0(j)*e3(j)) 
      A13(j)= 2*(e1(j)*e3(j) -e0(j)*e2(j))
      A21(j)= 2*(e1(j)*e2(j) -e0(j)*e3(j)) 
      A22(j)= e0(j)**2 -e1(j)**2 +e2(j)**2 -e3(j)**2
      A23(j)= 2*(e2(j)*e3(j) +e0(j)*e1(j))
      A31(j)= 2*(e1(j)*e3(j) +e0(j)*e2(j))
      A32(j)= 2*(e2(j)*e3(j) -e0(j)*e1(j))
      A33(j)= e0(j)**2 -e1(j)**2 -e2(j)**2 +e3(j)**2
      end do
!
!*******************************************************
      nCLp = np +nq          ! PE+AGCT, Counter/Coions, nCLp is first defined
      npqr5= ll              ! ll: the total particles, 46080 atoms 
      nr   = npqr5 -(np+nq)  ! Water
!
      npqr3= np +nq +3*nr/5 
      npqr = np +nq +  nr/5 
!*******************************************************
!
      if(io_pe.eq.1) then
        OPEN (unit=11,file=praefixc//'.06'//suffix2,             &
              status='unknown',position='append',form='formatted')
!
        write(11,*)
        write(11,*) 'Number of particles (Macro, Counter/Coions): init'
        write(11,*) '  np, nq, nr=',np,nq,nr
!
        write(11,780) npqr5_0,npqr5_1,npqr5_2
        write(11,783) npqr_0,npqr_1,npqr_2
  780   format(' npqr5_0, npqr5_1,npqr5_2=',3i8)
  783   format(' npqr_0, npqr_1, npqr_2  =    ',3i8)
        write(11,*)
!
        j= npqr-6 ! nCLp 
        do i= npqr5-29,npqr5,5  !nCLp+25,5 
        j= j +1
!
        write(11,786) i,xa(i),ya(i),za(i),j,xg(j),yg(j),zg(j)
        write(11,787) i+1,xa(i+1),ya(i+1),za(i+1)
        write(11,787) i+2,xa(i+2),ya(i+2),za(i+2)
        write(11,787) i+3,xa(i+3),ya(i+3),za(i+3)
        write(11,787) i+4,xa(i+4),ya(i+4),za(i+4)
  786   format('786 i,xa,ya,za.xg,.=',i6,3f10.2,2x,i6,3f10.2)
  787   format('787 i,xa,ya,za.....=',i6,3f10.2)
        end do
        write(11,*)
!
        close(11)
      end if
!
!* Fold back to (-L/2,L/2)
!
      do j= 1,npqr
      xg(j)= xg(j) - NINT(xg(j)/xleng)*xleng
      yg(j)= yg(j) - NINT(yg(j)/yleng)*yleng
      zg(j)= zg(j) - NINT(zg(j)/zleng)*zleng
!     xg(j) = xg(j) - nint(xg(j)/xleng -0.5d0)*xleng  NG
!     yg(j) = yg(j) - nint(yg(j)/yleng -0.5d0)*yleng 
!     zg(j) = zg(j) - nint(zg(j)/zleng -0.5d0)*zleng 
      end do
!
!
      if(io_pe.eq.1) then
        OPEN (unit=11,file=praefixc//'.06'//suffix2,             &
              status='unknown',position='append',form='formatted')
!
        j= npqr-6 ! nCLp 
        do i= npqr5-29,npqr5,5  !nCLp+25,5 
        j= j +1
!
        write(11,789) i,xa(i),ya(i),za(i),j,xg(j),yg(j),zg(j)
  789   format('789 i,xa,ya,za.xg,.=',i6,3f10.2,2x,i6,3f10.2)
        end do
        write(11,*)
!
        close(11)
      end if
!
!
!***************************************************************
!* The water molecules except top and bottom parts are defined *
!* in the DNA and pore region above                            *
!***************************************************************
!
      epslj_w = 1.02d-14 ! 2.189d-14 ! (kjoule/mol) in erg
!!    awat = 3.166d0 /2.d0 ! for hybrid LJ
!
!  phi= A/r^12 -B/r^6 -> d.phi/dr= -12A/r^13 +6B/r^7= 0
!       solve 12A/6B=r_eq^6 -> r_eq= (2A/B)^(1/6)= (1.7647^3)^(1/6) 
!               -> rcutlj(tip5p-Ewald)= r_eq= 3.4763 Ang       
!                  given in read_conf
!
      epslj_A = 3.8538d-08 ! erg, A12, tip5p/e with Ewald sums
      epslj_B = 4.3676d-11 ! erg, A6
!
!     epslj_A = 3.7856d-08 ! tip5p plain
!     epslj_B = 4.1041d-11
!     epslj_A = 4.1715d-08 ! tip4p
!     epslj_B = 4.2410d-11 
!

      q_O =  0.000d0
      q_H =  0.241d0
      q_L = -0.241d0  
!
      masso = 16.d0 /18 ! O, wwat=18
      massh =  1.d0 /18 ! H
!
      l= nCLp
      do i= nCLp+1,npqr5,5   ! O-H(1)-H(2)-L(1)-L(2)
      l= l +5
!
!  For Exc             TIP4P        SPC
!     ch(l-3)=  0    ! 0         ! 
!     ch(l-2)=  q_H  ! 0.4238d0  !  0.42 e
!     ch(l-1)=  q_O  !-0.8476d0  ! -0.85 e
!     ch(l  )=  q_H  ! 0.4238d0  !  0.42 e
!
      ch(l-4)=  q_O 
      ch(l-3)=  q_H
      ch(l-2)=  q_H
      ch(l-1)=  q_L
      ch(l  )=  q_L
!
      am(l-4)= masso 
      am(l-3)= massh 
      am(l-2)= massh 
      am(l-1)= 0
      am(l  )= 0
      end do
!
!
      j= 0
      do i= 1,nCLp
      j= j +1
      amm(j)= am(i)
      end do
!
      do i= nCLp+1,npqr5,5   ! O-H(1)-H(2)-L(1)-L(2)
      j= j +1
!              O        H        H
      amm(j)= am(i) +am(i+1) +am(i+2)
      ag (j)= 3.166d0/2.d0 
      ep (j)= epslj_w        ! LJ-potential = 4*ep(o)
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
        do j= np+1,np+5 !nCLp
        write(11,567) j,xg(j),yg(j),zg(j),amm(j),ag(j)
  567   format('567 j=',i8,' xg,yg,zg=',3f8.1,'  am,ag=',2f8.3)
        end do
!
        do j= nCLp+1,nCLp+100 !npqr
        write(11,567) j,xg(j),yg(j),zg(j),amm(j),ag(j)
        end do
!
        write(11,*) '--- Subroutine init (end) ------------------------'
        write(11,*)
        close(11)
      end if
!
!*******************************
!*  Velocity for i= np+1,npqr  *
!*******************************
!
      vmax1= vth/sqrt(38.d0/wwat)  !<-- mass of K(+)
      vmax2= vth/sqrt(18.d0/wwat)
!
      do j= np+1,nCLp
      vx(j)= 0 !dgaus2(vmax1) later at L.1100
      vy(j)= 0 
      vz(j)= 0
      end do
!
      do j= nCLp+1,npqr
      vx(j)= 0          !<- TIP5P water must be cold
      vy(j)= 0
      vz(j)= 0
      end do
!
!***************************
!* Dielectric constant     *
!***************************
!* used in /diel/ to calculate dielectric constant
!
!     call water_dens (xg,yg,zg,ipar,nCLp,npqr5)
!
      return
      end subroutine init_P
!
!
!--------------------------------------------------------------------
      subroutine get_rod (xg,yg,zg,rxg,ryg,rzg,np)
!--------------------------------------------------------------------
      use, intrinsic :: iso_c_binding 
      implicit none
!*
      include   'paramWatQa.h'
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
        write(11,*) 'Get_rod  np=',np           !<-- np= 25-36 
        write(11,*) 'xgc,ygc,zgc=',xgc,ygc,zgc  !   Counterion, not a rod !
!
        do i=1,np
        write(11,920) i,doxc(i),doyc(i),dozc(i)  !<-- i= 25-36 
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
! 
!---------------------------------------------------------------------
      subroutine get_torque (tqx,tqy,tqz,fcx,fcy,fcz,ftx,fty,ftz,np)
!---------------------------------------------------------------------
      use, intrinsic :: iso_c_binding 
      implicit none
!*
      include   'paramWatQa.h'
!
      real(C_DOUBLE),dimension(npqr50) :: fcx,fcy,fcz
      real(C_DOUBLE),dimension(np0)   :: ftx,fty,ftz 
      real(C_DOUBLE),dimension(np0)   :: dox,doy,doz,doxc,doyc,dozc
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
      include   'paramWatQa.h'
!
      integer(C_INT) np,nq,kshuf,k,i1,i2
      real(C_DOUBLE) ranff,ch(npqr50),sch
!
      kshuf= 3*nq*ranff(0.d0) +1.001
!
      do k= 1,kshuf
      i1= np+1 +mod(int(nq*ranff(0.d0)),nq)
      i2= np+1 +mod(int(nq*ranff(0.d0)),nq)
      if(ch(i1)*ch(i2).gt.0.d0) go to 130
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
      include   'paramWatQa.h'
!
      real(C_DOUBLE),dimension(npqr0) :: vx,vy,vz
      integer(C_INT) np,nq,nCLp,nr,npqr
!
!     character*1  char(8**3)
      real(C_DOUBLE) pi,dt,dth,axi,Gamma,rbmax,vth,tmax, &
                     wwat,awat
      common/parm2/  pi,dt,dth,axi,Gamma,rbmax,vth,tmax
      common/waterm/ wwat,awat
!
      integer(C_INT) i,ix,iy,iz,k,ILN
      real(C_float)  aiv,xsc(101),fvx(101),fvy(101),fvz(101), &
                     fmax1,fmax2,fmax3,fmin1,fmin2,fmin3
      real(C_DOUBLE) vmax0,vmax2
!
!*******************************
!*  Ions only.                 *
!*******************************
!
!       if(ifbase.eq.2) ch(i)= 0.d0  !  neutral chain 
!       am(i)=  94.d0/wwat
!     vmax0= 10.d0*vth/sqrt(218.d0/wwat)  !<- mass of sugar ring
!     vmax0= 10.d0*vth/sqrt( 38.d0/wwat)  !<- mass of K
      vmax0= 10.d0*vth/sqrt( 18.d0/wwat)  !<- 
      vmax2= 10.d0*vth/sqrt( 18.d0/wwat)  !<- water
!
!* (1) Coulomb particles.
      aiv= 50./vmax0
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
      integer(C_INT) is,is1,i
      real(C_float)  q(is),qav
!
      qav= 0
!
      is1= is
      if(is.ge.10) is1= 10
!
      do i= 1,is1
      qav= qav +q(i)
      end do
!
      qav= qav/is1
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
      include   'paramWatQa.h'
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
      subroutine poissn_eq (rho8,pot8,ndim,itrmax,iterp,ipar)
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
      include    'paramWatQa.h'
!
      integer(C_INT) ndim,itrmax,iterp,ipar 
!     parameter     (nob3=13)           ! paramWatQa.h
!
      real(C_DOUBLE),dimension(0:mxyz-1) :: rho8,pot8,xx8
!
      real(C_DOUBLE) Vtop,Vbot,Vtop0,Vbot0
      integer(C_INT) fixedbc_x,fixedbc_y,itermax,filtx,filty,filtz
      common/espot/  Vtop,Vbot,Vtop0,Vbot0,fixedbc_x,fixedbc_y, &
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
      common/cresm3/  ipr(10),rpr(10)  !<-- values of ipr(10), rpr(10)
!
      real(C_DOUBLE) ww(mxyz,6),wp(mxyz)     ! uses wp(1:mxyz)
      integer(C_INT) iw(mxyz),mxyz7,mxyz8,nob7
!
      integer(C_INT) iwrt1,iwrt2,iwrt3,iwrt4,ierror
      real(C_DOUBLE) sres,avex,rsdl,conv
      common/iotim/  iwrt1,iwrt2,iwrt3,iwrt4
      common/crconv/ sres,avex,rsdl,conv,ierror
!
      integer(C_INT) io_pe,i3,i4,cnt_recv3,disp_recv3, &
                     cnt_send3,i00
      common/sub_proc/ io_pe
      common/dat_typ3/ i3(30),i4(30),cnt_recv3(30),disp_recv3(30)
!
      logical     first /.true./
!*
      do ik= i3(ipar),i4(ipar)    ! the index i3(1)= 0
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
        call escof3_P (aa,xx,ss,na,ndim)
        first= .false.
      else
        call bound_s (xx,ss,Vtop,Vbot)
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
        cnt_send3= i4(ipar) -i3(ipar) +1
        call mpi_allgatherv (xx8, cnt_send3,           mpi_real8, &
                             pot8,cnt_recv3,disp_recv3,mpi_real8, &
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
      end subroutine poissn_eq 
!
!
!-----------------------------------------------------------------------
      subroutine escof3_P (aa,pot8,rho8,na,ndim)
!----------------------------------------------------------------------
!  Periodic (x,y) and bounded (z) boundaries
!
      use, intrinsic :: iso_c_binding 
      implicit none
!
      include    'paramWatQa.h'
!
!     integer(C_INT) nob3     !<- paramWatQa.h
!     parameter     (nob3=13)
!
!*---------------------------------------------------------------------
      real(C_DOUBLE) aa(0:mxyz-1,nob3)
      integer(C_INT) na(0:mxyz-1,nob3)
!
      real(C_DOUBLE),dimension(0:mx-1,0:my-1,0:mz-1) :: pot8,rho8,dec2
      integer(C_INT),dimension(nob3) :: lai,laj,lak
      real(C_DOUBLE) ca(nob3)
!
      real(C_DOUBLE) Vtop,Vbot,Vtop0,Vbot0
      integer(C_INT) fixedbc_x,fixedbc_y,itermax,filtx,filty,filtz
      common/espot/  Vtop,Vbot,Vtop0,Vbot0,fixedbc_x,fixedbc_y, &
                     itermax,filtx,filty,filtz
!*---------------------------------------------------------------------
!
      integer(C_INT) pxr,pxc,pxl,pyr,pyc,pyl,pzr,pzc,pzl
      common/ptable/ pxr(-10:mx+10),pxc(-10:mx+10),pxl(-10:mx+10), &
                     pyr(-10:my+10),pyc(-10:my+10),pyl(-10:my+10), &
                     pzr(-10:mz+10),pzc(-10:mz+10),pzl(-10:mz+10)
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
      common/xtabl2/ ptx(-199:3000),pty(-199:3000),ptz(-100:4000)
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
!*                   !!<-- escof3P, i,j: periodic, k: bounded
!                          periodic means pxc() and pyc() defined in /init/
      do i= 0,mx-1 
      do j= 0,my-1
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
      laj(2)=  pyc(j-1)
      lak(2)=   k
       ca(2)=   1.d0/ghysq(j)
!
!* (i-1,j,k)
      lai(3)=  pxc(i-1)
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
      lai(5)=  pxc(i+1)
      laj(5)=   j
      lak(5)=   k
       ca(5)=   1.d0/ghxsq(i)
!
!* (i,j+1,k)
      lai(6)=   i
      laj(6)=  pyc(j+1)
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
      lai(8)=  pxc(i-1)
      laj(8)=   j
      lak(8)=   k
       ca(8)=  -(dec2(pxc(i+1),j,k)-dec2(pxc(i-1),j,k))/(ghx2(i)**2*dec2(i,j,k))
!
      lai(9)=  pxc(i+1)
      laj(9)=   j
      lak(9)=   k
       ca(9)=  (dec2(pxc(i+1),j,k)-dec2(pxc(i-1),j,k))/(ghx2(i)**2*dec2(i,j,k)) 
!
!* (deps/dy)/eps
      lai(10)=  i
      laj(10)= pyc(j-1)
      lak(10)=  k
       ca(10)= -(dec2(i,pyc(j+1),k)-dec2(i,pyc(j-1),k))/(ghy2(j)**2*dec2(i,j,k))
!
      lai(11)=  i
      laj(11)= pyc(j+1)
      lak(11)=  k
       ca(11)=  (dec2(i,pyc(j+1),k)-dec2(i,pyc(j-1),k))/(ghy2(j)**2*dec2(i,j,k))
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
      ii= pxc(lai(m))      ! periodic, pxc()
      jj= pyc(laj(m))
      kk=     lak(m)       !<-- bounded 
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
      end do  !<- end of j
      end do  !<- end of i
!
      return
      end subroutine escof3_P
!
!
!-----------------------------------------------------------------------
      subroutine bound_s (pot8,rho8,Vtop,Vbot)
!----------------------------------------------------------------------
! ******************************
! * Boundary values for source *
! ******************************
      use, intrinsic :: iso_c_binding 
      implicit none
!
      include    'paramWatQa.h'
!     integer(C_INT) nob3   !<- paramWatQa.h
!
      integer(C_INT) i,j,k
      real(C_DOUBLE),dimension(0:mx-1,0:my-1,0:mz-1) :: pot8,rho8
      real(C_DOUBLE) Vtop,Vbot
!*---------------------------------------------------------------------
!     ez1(i,j,k)= - (Vtop-Vbot)/zleng   ! ez1= -grad pot
!
      do k= 0,mz-1
      do j= 0,my-1
      do i= 0,mx-1
!
      if(k.eq.0 .or. k.eq.mz-1) then
        if(k.eq.   0) pot8(i,j,k)= Vbot  !<-- present value of potentals
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
      include  'paramWatQa.h'
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
      if(mod(itr,50).eq.1 .and. io_pe.eq.1) then
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
      include   'paramWatQa.h'
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
      include    'paramWatQa.h'
!
!     integer(C_INT) nob3
!     parameter  (nob3=13)
!
      real(C_DOUBLE),dimension(0:mxyz-1) :: xx,ss
      integer(C_INT) ierr,ipar
!
      real(C_DOUBLE) rpr,eps,sq0,sq1
      integer(C_INT) ipr
      common/cresm3/ ipr(10),rpr(10)
!
      integer(C_INT) io_pe,i3,i4,cnt_recv3,disp_recv3,ncrs,itrm
      common/sub_proc/ io_pe
      common/dat_typ3/ i3(30),i4(30),cnt_recv3(30),disp_recv3(30)
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
      include    'paramWatQa.h'
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
      integer(C_INT)  io_pe,i3,i4,cnt_recv3,disp_recv3
      common/sub_proc/ io_pe
      common/dat_typ3/ i3(30),i4(30),cnt_recv3(30),disp_recv3(30)
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
      include    'paramWatQa.h'
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
      integer(C_INT) io_pe,i3,i4,cnt_recv3,disp_recv3
      common/sub_proc/ io_pe
      common/dat_typ3/ i3(30),i4(30),cnt_recv3(30),disp_recv3(30)
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
      include   'paramWatQa.h'
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
      integer(C_INT) io_pe
      integer(C_INT) i3,i4,cnt_recv3,disp_recv3
      integer(C_INT) cnt_recv5,disp_recv5,disp_recv6,cnt_send5
      common/sub_proc/ io_pe
      common/dat_typ3/ i3(30),i4(30),cnt_recv3(30),disp_recv3(30)
      common/dat_typ5/ cnt_recv5(30),disp_recv5(30),disp_recv6(30)
!
!  Bounded in z=0 and z= zmax
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
        disp_recv5(k)= i3(k)
        disp_recv6(k)= i4(k)-mxy+1
!
        cnt_recv5(k)= mxy
        end do
!
        i00= i3(ipar)
        do i= i3(ipar),i3(ipar)+mxy-1
        vv(i-i00)= v(i)
        end do
!
        cnt_send5= mxy
        call mpi_allgatherv (vv,cnt_send5,           mpi_real8, &
                             v8,cnt_recv5,disp_recv5,mpi_real8, &
                             mpi_comm_world,ierror)
!
        i00= i4(ipar)-mxy+1
        do i= i4(ipar)-mxy+1,i4(ipar)
        vv(i-i00)= v(i)
        end do
!
        cnt_send5= mxy
        call mpi_allgatherv (vv,cnt_send5,           mpi_real8, &
                             v8,cnt_recv5,disp_recv6,mpi_real8, &
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
      include      'paramWatQa.h'
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
      include      'paramWatQa.h'
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
                    vzco1,vzco2
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
      call lplot1 (2,4,is,time,ekin,emax1,0.0,iln,'Kin ions',8, & !ekin
                 '        ',8,'        ',8)
      call lplot1 (2,5,is,time,ekn2,emax2,0.0,iln,'Kin watr',8, & !ekn2
                 '        ',8,'        ',8)
      call lplot1 (2,6,is,time,ppot,pmax,pmin,iln,'es ener ',8, & !ppot
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
!-------------------------------------------------------------------------
      subroutine ppl3da (xe,ye,ze,chg,rod_leng,Rmac,np,nq,nCLp, &
                         npqr3,first_ppl)
!-------------------------------------------------------------------------
      use, intrinsic :: iso_c_binding 
      implicit none
!
      include   'paramWatQa.h'
!
      real(C_DOUBLE),dimension(npqr30) :: xe,ye,ze,chg,x,y,z
!
      real(C_DOUBLE) pi,dt,dth,axi,Gamma,rbmax,vth,tmax,    &
                     xmax,ymax,zmax,xmin,ymin,zmin,         &
                     Bjerrum,qfrac,Rpore,Hpore,Zci,Zcp,Zcn, &
                     rod_leng,Rmac
      real(C_float)  phi,tht,dtwr1,dtwr2,dtwr3
      integer(C_INT) nsg,nseg,ifqq,np,nq,nCLp,npqr3,j 
!
      common/parm2/  pi,dt,dth,axi,Gamma,rbmax,vth,tmax
      common/parm3/  xmax,ymax,zmax,xmin,ymin,zmin
      common/parm4/  phi,tht,dtwr1,dtwr2,dtwr3
      common /elsta/ Bjerrum
      common/cntion/ qfrac,Rpore,Hpore,Zci,Zcp,Zcn,ifqq
      common/psegm/  nsg(30),nseg
!
      integer(C_INT) np00,nnb,ist1,ist2,n_rodp,jj
      common/pbase/  np00,nnb,ist1(100),ist2(100)
!
      real(C_DOUBLE) Vtop,Vbot,Vtop0,Vbot0,rpr
      real(C_DOUBLE) diel2,dielpr,daa,dbb
      integer(C_INT) fixedbc_x,fixedbc_y,itermax,filtx,filty,filtz
      common/espot/  Vtop,Vbot,Vtop0,Vbot0,fixedbc_x,fixedbc_y, &
                     itermax,filtx,filty,filtz
      common/dielec/ diel2,dielpr,daa,dbb
!
      character*8    label,cdate*10,cax*1
      real(C_float)  t,t00,xp_leng
      common/headr1/ label,cdate
      common/headr2/ t,t00,xp_leng
!
      integer(C_INT) io_pe
      common/sub_proc/ io_pe
!
      integer(C_INT) i,k,kpl,ns,ia,ib,l,nml
      real(C_float)  hh,Rpore4,Hpore4,Zcp4,Zcn4,Bjerrum4,xleng4,zleng4, &
                     Rmac4,Vtop4,Vbot4,Gamma4,rod_leng4,diel24,         &
                     fsize,hl,vd,pha,tha,cph,sph,cth,sth,xp,yp,zp,      &
                     rmax1,ps,x1,y1,z1,xx,yy,dd,xpp,ypp,zpp
      real(C_DOUBLE) xleng,yleng,zleng
      common/parm8/  xleng,yleng,zleng
!
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
      do i= 1,npqr3
      x(i) = xe(i)
      y(i) = ye(i)
      z(i) = ze(i)
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
      if(i.gt.nCLp) dd= 1.0*ps  ! 0.5*ps
!
      if(i.le.np) then
        if(chg(i).lt.0.d0) then
          call newcolor(3,1.,0.,0.)
          call circle (xx-0.15,yy-0.15,dd,2)  ! 2: filled circle, red
!
        else if(chg(i).ge.0.d0) then
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
      call newcolor(0,1.,0.,0.)  ! Reset
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
        write(11,300) k,ist1(k),ist2(k)
  300   format('k=',i4,'  ist1(k),ist2(k)=',2i6)
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
!
      first_ppl= .false.
!
      call plot (xx,yy,3)
      call newcolor(0,1.,0.,0.)  ! Reset
!
!------------------------
!*  Solvent water
!------------------------
!
      if(ns.eq.2) then
        nml= 3 * 5     ! in 5 times by 3
!
        do i= nCLp+1,npqr3
        if(mod(i-nCLp,nml).ge.1 .and. mod(i-nCLp,nml).le.3) then
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
        dd= 0.8*ps  ! 0.5*ps
        call circle (xx-0.15,yy-0.30,dd,2)
!
        end if
        end do
      end if
!
      call plot (xx,yy,3)
      call newcolor(0,1.,0.,0.)  ! Reset
!
!   Make plots at 1) i=1,np+nq, and 2) i=1,npqr3 
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
      include   'paramWatQa.h'
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
      include   'paramWatQa.h'
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
      common/xtabl2/ ptx(-199:3000),pty(-199:3000),ptz(-100:4000)
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

      include   'paramWatQa.h'
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
      common/ptable/ pxr(-10:mx+10),pxc(-10:mx+10),pxl(-10:mx+10), &
                     pyr(-10:my+10),pyc(-10:my+10),pyl(-10:my+10), &
                     pzr(-10:mz+10),pzc(-10:mz+10),pzl(-10:mz+10)
!
      real(C_DOUBLE) gx,gy,gz,ghx,ghy,ghz,ghx2,ghy2,ghz2, &
                     ghxsq,ghysq,ghzsq,hxi,hyi,hzi
      integer(C_INT) ptx,pty,ptz
      common/xtable/ gx(0:mx),ghx(0:mx),ghx2(0:mx),ghxsq(0:mx), &
                     gy(0:my),ghy(0:my),ghy2(0:my),ghysq(0:my), &
                     gz(0:mz),ghz(0:mz),ghz2(0:mz),ghzsq(0:mz), &
                     hxi,hyi,hzi
      common/xtabl2/ ptx(-199:3000),pty(-199:3000),ptz(-100:4000)
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
!     do i= 1,mx-2
      do i= 0,mx-1
      npx= npx +1
      ir= i+1
      il= i-1
      if(i.eq.0) il= pxl(0) 
      if(i.eq.mx-1) ir= pxr(mx-1)
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
!     do j= 1,my-2
      do j= 0,my-1
      npy= npy +1
      jr= j+1
      jl= j-1
      if(j.eq.0) jl= pyl(0)
      if(j.eq.my-1) jr= pyr(my-1) 
!
      npx2= 0
!
      do i= 0,mx-1
      npx2= npx2 +1
      ir= i+1
      il= i-1
      if(i.eq.0) il= pxl(0) 
      if(i.eq.mx-1) ir= pxr(mx-1)
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
!
          write(11,*)
          write(11,*) 'Contor at (x,z)  t8=',t8 
          write(11,*) 'npx,npz=',npx,npz
          write(11,*) 'xmin,zmin,xmax,zmax=',xmin,zmin,xmax,zmax
          write(11,*) 'xl1,zl,xr1,zr=',xl1,zl,xr1,zr
          write(11,*) 'a(nxz): wamax,wamin,=',wamax,wamin
!
!         do ik= 1,200,5
!         write(11,990) a(ik),a(ik+1),a(ik+2),a(ik+3),a(ik+4)
! 990     format(1p5d12.3)
!         end do
!
!         do ik= nxz-299,nxz,5
!         write(11,990) a(ik),a(ik+1),a(ik+2),a(ik+3),a(ik+4)
!         end do
!         write(11,*)
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
      include  'paramWatQa.h'
!
      real(C_DOUBLE) q(0:mx-1,0:my-1,0:mz-1),sym,sym2
      real(C_float)  a(-2:mx+1,-2:my+1,-2:mz+1)
      integer(C_INT) filtx,filty,filtz,ntx,nty,ntz,i,j,k,is,   &
                     il,ill,ir,irr,js,jl,jll,jr,jrr,ks,kl,kll, &
                     kr,krr  
!
      integer(C_INT)  pxr,pxc,pxl,pyr,pyc,pyl,pzr,pzc,pzl
      common/ptable/ pxr(-10:mx+10),pxc(-10:mx+10),pxl(-10:mx+10), &
                     pyr(-10:my+10),pyc(-10:my+10),pyl(-10:my+10), &
                     pzr(-10:mz+10),pzc(-10:mz+10),pzl(-10:mz+10)
!
      real(C_DOUBLE) gx,gy,gz,ghx,ghy,ghz,ghx2,ghy2,ghz2, &
                     ghxsq,ghysq,ghzsq,hxi,hyi,hzi
      integer(C_INT) ptx,pty,ptz
      common/xtable/ gx(0:mx),ghx(0:mx),ghx2(0:mx),ghxsq(0:mx), &
                     gy(0:my),ghy(0:my),ghy2(0:my),ghysq(0:my), &
                     gz(0:mz),ghz(0:mz),ghz2(0:mz),ghzsq(0:mz), &
                     hxi,hyi,hzi
      common/xtabl2/ ptx(-199:3000),pty(-199:3000),ptz(-100:4000)
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
      real(C_float) u(17000),w(17000)  !<-- nx,ny
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
      include     'paramWatQa.h'
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
