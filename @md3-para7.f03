!*-------------------------------------------------- 2025/4/01 ----*
!*                                                                 *
!* << Simple Program of "3D Electrostatic Simulation Code" >>      *
!*                                                                 *
!*   This program explains a very easy code showing Fortran 2003   *
!*  compilation and then parallel execution, using Linux OS of     *
!*  6 cores of 3 GMz. It uses MPI v.3 and the Ghostviewer script.  *
!*                                                                 *
!*  Code: @md3-para7.f03                                           *
!*        paramE7.h                                                *
!*                                                                 *
!*  Author: Motohiko Tanaka, Ph.D., Chikusa, Nagoya 464, Japan     *
!*                                                                 *
!*  Released by GPL-3.0 License, https://github.com/Mtanaka77/     *
!*  Copyright(C) 2006-2025. All rights reserved.                   *
!*                                                                 *
!*-----------------------------------------------------------------*
!*$ mpif90 -mcmodel=medium -fPIC -O2 -o a.out @md3-para7.f03       *
!* For one cpu execution, % mpiexec -n 1 a.out &  (very slow !)    *
!*-----------------------------------------------------------------*  
!
      program es3d_CL
!
      use, intrinsic :: iso_c_binding
      implicit none
!
      include    'paramE7.h'
      include    'mpif.h'
!
      integer(C_INT) ipar,igrp,size,rank,key,ierror
      real(C_DOUBLE) cptot,tmax,wtime,walltime1,walltime2
!
      integer(C_INT) kstart,irefl,io_pe,nframe,it_write06, &
                     it_write23
      common/sub_proc/ io_pe
!
      real(C_float) phi,tht,rgmax
      common/parm4/ phi,tht,rgmax
!
!
! Setup of MPI
      call mpi_init (ierror)
      call mpi_comm_rank (MPI_COMM_WORLD,rank,ierror)
      call mpi_comm_size (MPI_COMM_WORLD,size,ierror)
!
      ipar = 1 + rank
      igrp = 1
!
      io_pe= 0
      if(ipar.eq.1) io_pe= 1  !<- For major node
!
      if(io_pe.eq.1) then     !<- Initial time with unit=11 and 77
        open (unit=11,file=praefixc//'.06'//cname,form='formatted')
!       
        write(11,*) 'Molecular Dynamics (parallel)'
        write(11,*) ' Number of processors=',size
        write(11,*) '  num_proc=',num_proc
        write(11,*)
        close(11)
!
        open (unit=77,file=praefixc//'.77'//cname//'.ps',    &
                            status='replace',form='formatted')
!                                   +++++++
        nframe= 4
        call gopen (nframe)
        close(77)
      end if
!
! ++++++++++++++++++++++++++++++++++++++++++++++
!  Parameters of this run to start /moldyn/
! ++++++++++++++++++++++++++++++++++++++++++++++
!   The parameter file paramE7.h for basic ones
!
      irefl = 0         ! =0 for an open system
      kstart= 0         ! =0 for an initial time
!
      tmax = 7.0d+3     ! physical time
      cptot= 70.3d0     ! minutes   3.d0
!
      phi= -60.
      tht=  15.
!
      it_write06=  500 ! steps   200 !  2000 
      it_write23= 5000 ! steps   400 ! 10000 ! 50000 
!
! ---------------------------------------------------------
      call moldyn (ipar,igrp,irefl,kstart,tmax,cptot,  &
                   it_write06,it_write23)
! ---------------------------------------------------------
!
!* Closing of FT77 and /mpi_finalize/
!
      if(io_pe.eq.1) then
        open (unit=77,file=praefixc//'.77'//cname//'.ps',        &
              status='unknown',position='append',form='formatted')
!       
        call gclose
        close(77)
      end if
!
      call mpi_finalize (ierror)
!
      stop
      end program es3d_CL
!
!
!------------------------------------------------------------------------
      subroutine date_and_time_7 (date_now,time_now)
!------------------------------------------------------------------------
      use, intrinsic :: iso_c_binding  ! <-
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
!---------------------------------------------------------------
      subroutine moldyn (ipar,igrp,irefl,kstart,tmax,cptot,  &
                         it_write06,it_write23)
!---------------------------------------------------------------
      use, intrinsic :: iso_c_binding
      implicit none
!
      include    'paramE7.h'
      include    'mpif.h'
!
      integer(C_INT) ipar,igrp,irefl,kstart,np,nq,nCLp,io_pe, &
                     it_write06,it_write23
      real(C_DOUBLE) tmax,cptot
      common/sub_proc/ io_pe
!
      real(C_DOUBLE),dimension(npq0) :: &
                           x,y,z,px,py,pz,vx,vy,vz,ch,am,ag
      real(C_DOUBLE),dimension(npq0,3) :: ffr
!
      integer(C_INT) i,k,it
      real(C_DOUBLE) xmax,ymax,zmax,xmin,ymin,zmin,        &
                     xleng,yleng,zleng,                    &
                     E_C_s,E_grid,E_C_r,E_LJ,E_elas,q0,    &
                     a_unit,m_unit,c_unit,e_unit,          &
                     Wrest,kJoule,kcal,mol,Temp,Temp_erg,  &
                     stime,t,t8
      real(C_float)  dt0,s0,s1,etot
!
      real(C_DOUBLE) pi,dt,pref_CL,pref_LJ,pthe
      common/parm2/  pi,dt,pref_CL,pref_LJ,pthe
      common/parm3/  xmax,ymax,zmax,xmin,ymin,zmin
      common/parm5/  xleng,yleng,zleng
      COMMON/ENERGY/ E_C_s,E_grid,E_C_r,E_LJ,E_elas
      common/tempr/  Temp
!
      real(C_DOUBLE) epsLJ,gamma,W_1p,Pot0,               &
                     wtime,wtime1,wtime2,                 &
                     R_sp,D_sp,N_sp,W_sp,massi,Lambda_D,  &
                     ch_ion,wt_ion,rd_ion,ch_el,wt_el,rd_el
      common/ionsiz/ R_sp,D_sp,N_sp,massi,Lambda_D,       &
                     ch_ion,wt_ion,rd_ion,ch_el,wt_el,rd_el 
!
      character(len=8) label,date_now*10,time_now*8
      real(c_float)  xp_leng
      common/HEADR1/ label,date_now
      common/HEADR2/ t8,xp_leng
!
      real(C_DOUBLE) energ1,energ2
!
      real(C_float)  tcpu,cputime
      real(C_DOUBLE) walltime1,walltime2,walltime3
      logical        first /.true./
!    
      walltime1= mpi_wtime()
!
      call date_and_time_7 (date_now,time_now)
!
      if(io_pe.eq.1) then
        open (unit=11,file=praefixc//'.06'//cname,               &
              status='unknown',position='append',form='formatted')
!
        write(11,100) date_now,time_now
  100   format(' date: ',a10,'  time: ',a8,/)
        close(11)
      end if
!  --------------------------------------------
!
      pi = 4.d0*atan(1.d0)
!
      dt = 0.005d0        ! time t
      np = np0
      nq = nq0
!
      a_unit= 1.0000d-4   ! a_unit= 1.d-4 cm= 1 micron
      m_unit= 0.9109d-27
      c_unit= 2.9979d+10  ! t_unit= a_unit/c_unit= 1.0d-4/3.0d10= 3.33d-15 s
      e_unit= 4.8032d-10
      Wrest = m_unit * c_unit**2
!
!     parameter  (R_sp0=1.0d-6,D_sp0=1.0d21) <- parameters
      massi= 1836.d0 ! 100.d0
      R_sp= R_sp0
      D_sp= D_sp0
!
      if(io_pe.eq.1) then
        open (unit=11,file=praefixc//'.06'//cname,               &
              status='unknown',position='append',form='formatted')
!       
        if(irefl.eq.0) then
          write(11,*) 'Open system'
!
        else if(irefl.eq.1) then
          write(11,*) 'Closed system'
        end if
!
        close(11)
      end if
!
      if(irefl.eq.1) then
        xleng= 1.d-4/a_unit
        yleng= 1.d-4/a_unit
        zleng= 1.d-4/a_unit
      else
        xleng= 1.d-5/a_unit
        yleng= 1.d-5/a_unit
        zleng= 1.d-5/a_unit
      end if
!
      xmax=  0.5d0*xleng
      xmin= -0.5d0*xleng
      ymax=  0.5d0*yleng
      ymin= -0.5d0*yleng
      zmax=  0.5d0*zleng
      zmin= -0.5d0*zleng
!
      epsLJ= 0.0148d0
      Temp= 4.5d0     !<- cold temperature in K
!
      N_sp = (4.d0*pi/3.d0)* R_sp**3 * D_sp
      W_sp = N_sp * massi 
!
      ch_ion = 1.d0 
      wt_ion = massi 
      rd_ion = (R_sp/1.d-5)/np**(1.d0/3.d0)
!     rd_ion = (R_sp/1.d-4)/np**(1.d0/3.d0)
!
      ch_el = -1.d0 
      wt_el =  1.d0 
      rd_el =  rd_ion * (wt_el/wt_ion)**(1.d0/3.d0)
!
      pref_CL = (e_unit**2/a_unit) /Wrest
      Temp_erg= Temp*1.6d-12
      pthe    = sqrt( (Temp_erg/Wrest +1.d0)**2 -1.d0)
!
      kJoule = 1.d10
      kcal   = 4.1868d0 *kJoule
      mol    = 6.0220d23
!
      pref_LJ = 48.0d0*epsLJ*(kcal/mol) /Wrest
      Lambda_D= sqrt(Temp_erg /(4*pi*D_sp*e_unit**2))
!
!
      if(io_pe.eq.1) then
        open (unit=11,file=praefixc//'.06'//cname,               &
              status='unknown',position='append',form='formatted')
!
        write(11,*) 'a_unit= 1.0000d-04 cm = 1 micron'
        write(11,*) 'm_unit= 0.9109d-27 g'
        write(11,*) 'c_unit= 2.998d+10 cm/s'
        write(11,*) 'e_unit= 4.803d-10 esu'
        write(11,*) 'dt=',dt,',  t_unit= a_unit/c_unit= 3.33x10^-15 s' 
        write(11,*)
        write(11,*) 'tmax =',tmax,' physical time'
        write(11,*) 'cptot=',cptot,' minutes'
        write(11,*)
        write(11,*) 'np0,nq0=',np0,nq0
        write(11,*) 'massi=',massi
        write(11,*) '  real mass for hydrogen !'
        write(11,*) 'Open (0)/Closed (1) boundary=',irefl
        write(11,*) 'R_sp=',R_sp
        write(11,*) 'D_sp=',D_sp
        write(11,*) 'n_sp=',n_sp
        write(11,*) 'w_sp=',w_sp
        write(11,*) 'temperature in K,erg=',Temp,Temp_erg
        write(11,*)
!
        close(11)
      end if
!
      call ggauss 
!
!  Restart in parameter paramE7.h
!
      if(kstart.eq.0) then  
!
        call init (x,y,z,px,py,pz,ch,am,ag,ipar,igrp, &
                   np,nq,nCLp,irefl)
!
      else if(kstart.eq.1) then
!
        open (unit=12,file=praefixc//'.12'//cname, &
                    status='old',form='unformatted')
!
        read(12) irefl,np,nq,nCLp,it_write06,it_write23
        read(12) x,y,z,px,py,pz
        read(12) ch,am,ag
!
        close(12)
!
        if(io_pe.eq.1) then
          open (unit=11,file=praefixc//'.06'//cname,               &
                status='unknown',position='append',form='formatted')
!
          write(11,*) 'Start/restart: it_write06: it=',it
          close(11)
        end if
!*
      end if
!
!
      if(io_pe.eq.1) then
      if(first) then
!
        open (unit=23,file=praefixc//'.23'//cname,form='formatted')
!
        first= .false.
!
        write(23,*) 'A cube is R_sp= 1 micron and D_sp= 10^21/cm^3'
        write(23,*) 'charge is defined as ch_i= e/np and ch_e= -e/nq'
        write(23,230) np,nq,nCLp
  230   format('np,nq,nCLp=',3i10,/)
!
        do i= 1,np+nq
        write(23,231) x(i),y(i),z(i),px(i),py(i),pz(i),  &
                        ch(i),am(i),ag(i)
  231   format(3d12.4,2x,3d12.4,2x,3d12.4)
        end do
!
        close(23)
!
      end if
      end if
!
!
      t8= 0.d0
      t = 0.d0
      it= 0
!
      wtime2= mpi_wtime()
!
      if(it.eq.0) then
      if(io_pe.eq.1) then
        open (unit=11,file=praefixc//'.06'//cname,               &
              status='unknown',position='append',form='formatted')
!
        write(11,*) 'tmax,cptot=',tmax,cptot
        write(11,*)
        close(11)
!
      end if
      end if
!
      walltime2= mpi_wtime()
!------------------------------------------
!** Loop starts here **********************
 1000 continue
      wtime=  mpi_wtime() -walltime1
      tcpu =  wtime/60.d0
!
      if(t8.gt.tmax) then
        if(io_pe.eq.1) then
          open (unit=11,file=praefixc//'.06'//cname,               &
                status='unknown',position='append',form='formatted')
          write(11,*) 'Go to 2000'
          close(11)
        end if
!
        go to 2000
      end if
!
      if(tcpu.gt.cptot) go to 2000
!
!
      t8 = t8 +dt
      t  = t8
      it = it +1
!
      E_c_r  = 0.d0
      E_LJ   = 0.d0
!
      do i = 1,nCLp
      ffr(i,1)= 0.d0
      ffr(i,2)= 0.d0
      ffr(i,3)= 0.d0
      end do
!
      call cl_forces (x,y,z,ch,am,ag,ffr,ipar,nCLp)
!
      do i= 1,nCLp
      px(i)= px(i) + ffr(i,1)*dt
      py(i)= py(i) + ffr(i,2)*dt
      pz(i)= pz(i) + ffr(i,3)*dt
!
      gamma= am(i) *sqrt(1.d0 +(px(i)**2+py(i)**2+pz(i)**2)/am(i)**2)
      vx(i)= px(i)/gamma
      vy(i)= py(i)/gamma
      vz(i)= pz(i)/gamma
!
      x(i)= x(i) + vx(i)*dt
      y(i)= y(i) + vy(i)*dt
      z(i)= z(i) + vz(i)*dt
      end do
!
      if(irefl.eq.1) then
        call reflect (x,y,z,vx,vy,vz,ch,ag,nCLp,irefl)
      end if
!
!* Energy history
      if(mod(it,it_write06).eq.1) then
      if(io_pe.eq.1) then 
!
        energ1= 0
        energ2= 0
!
        do i= 1,np0
        energ1= energ1 +am(i)*(vx(i)**2 +vy(i)**2 +vz(i)**2)
        end do
!
        do i= np0+1,np0+nq0
        energ2= energ2 +am(i)*(vx(i)**2 +vy(i)**2 +vz(i)**2)
        end do
!
        open (unit=11,file=praefixc//'.06'//cname,               &
              status='unknown',position='append',form='formatted')
!
        write(11,650) t8,tcpu,energ1,energ2,E_c_r,E_LJ
  650   format('# t8,tcpu=',f8.1,f8.2,'  energy1,2=',1p2d11.3, &
               '  E_c_r,E_LJ=',2d11.3)
        close(11)
!
      end if
      end if
!
!*
      if(mod(it,it_write23).eq.1) then
      if(io_pe.eq.1) then 
!
        open (unit=23,file=praefixc//'.23'//cname,               &
              status='unknown',position='append',form='formatted')
!
        write(23,240) np,nq,nCLp
  240   format(3i10)
!
        write(23,241) x,y,z,px,py,pz,ch,am,ag
  241   format(3d12.4,2x,3d12.4,2x,3d12.4)
        close(23)
!
!
        open (unit=12,file=praefixc//'.12'//cname,    &
                   status='replace',form='unformatted')
!
        write(12) irefl,np,nq,nCLp,it_write06,it_write23
        write(12) x,y,z,px,py,pz
        write(12) ch,am,ag
        close(12)
!
!
        open (unit=11,file=praefixc//'.06'//cname,               &
              status='unknown',position='append',form='formatted')
!
        write(11,700) it,t8
  700   format('% /ppl3da/ is plotted, it, t8=',i8,f8.2)
        close(11)
!
        open (unit=77,file=praefixc//'.77'//cname//'.ps',        &
              status='unknown',position='append',form='formatted')
!
        call ppl3da (x,y,z,ch,ag,np,nq,nCLp,irefl)
        close(77)
!
      end if 
      end if 
!*
      go to 1000
!** Loop ends here ************************
!
 2000 continue
      if(io_pe.eq.1) then
!*
        open (unit=11,file=praefixc//'.06'//cname,               &
              status='unknown',position='append',form='formatted')
!
        walltime3= mpi_wtime()
        write(11,730) walltime3-walltime1,    &
                      walltime2-walltime1,walltime3-walltime2
  730   format(/,'# Walltime total=',f10.3,'  wallt2,wallt3=',2f10.3,/)
!
        call date_and_time_7 (date_now,time_now)
!
        write(11,760) date_now,time_now
  760   format('# Write to FT12 and closing gclose of this job',/, &
               '# Final date: ',a10,'  time: ',a8,/)
        close(11)
!
!
        open (unit=12,file=praefixc//'.12'//cname,    &
                   status='replace',form='unformatted')
!
        write(12) irefl,np,nq,nCLp,it_write06,it_write23
        write(12) x,y,z,px,py,pz
        write(12) ch,am,ag
!
        close(12)
!*
      end if
!
      return
      end subroutine moldyn
!
!
!------------------------------------------------------------------
      subroutine cl_forces (x,y,z,ch,am,ag,ffr,ipar,nCLp)
!------------------------------------------------------------------
      use, intrinsic :: iso_c_binding
      implicit none
!*
      include    'paramE7.h'
      include    'mpif.h'
!
      integer(C_INT) ipar,nCLp
      real(C_DOUBLE),dimension(npq0) :: x,y,z,ch,am,ag
      real(C_DOUBLE),dimension(npq0,3) :: ffr,ffc
!
      integer(C_INT) i,j,k,l,itabl,ierror
      REAL(C_DOUBLE) dx,dy,dz,r2,r,forceV,ccel,driwu,driwu2, &
                     addpot,asc2,rlj,rsi,snt,rcl2,asc3,      &
                     pi,dt,pref_CL,pref_LJ,pthe,             &
                     rcutlj,rcutlj2,unif1(2),unif2(2),       &
                     E_C_s,E_C_PME,E_c_r,E_LJ,E_elas
      common/parm2/  pi,dt,pref_CL,pref_LJ,pthe
      COMMON /ENERGY/ E_C_s,E_C_PME,E_c_r,E_LJ,E_elas
!
      real(C_DOUBLE) R_sp,D_sp,N_sp,W_sp,Lambda_D,massi,   &
                     ch_ion,wt_ion,rd_ion,ch_el,wt_el,rd_el
      common/ionsiz/ R_sp,D_sp,N_sp,massi,Lambda_D,        &
                     ch_ion,wt_ion,rd_ion,ch_el,wt_el,rd_el 
!*---------------------------------------------------------------
!
      rcutlj = 1.122462048d0
      rcutlj2= rcutlj**2
!
      driwu2 = 1.25992104989487316476721060728d0
      driwu  = sqrt(driwu2)
!
      addpot =  1.d0/rcutlj2**6 - 1.d0/rcutlj2**3
      asc2 = (0.85d0)**2
      asc3 = (0.001d0)**2
!
      do i= ipar,nCLp,num_proc
      do j= 1,nCLp
      if(j.eq.i) go to 300
!
      dx = x(i) -x(j)
      dy = y(i) -y(j)
      dz = z(i) -z(j)
!
      r2 = dx**2 + dy**2 + dz**2
      r  = sqrt(r2)
      rcl2= max(r**2,asc3)
!
      forceV = pref_CL * ch(i)*ch(j)/(rcl2 *r)
      E_C_r = E_C_r + pref_CL*ch(i)*ch(j)/r
!
      rlj = r/(ag(i)+ag(j))
      ccel= 0.d0
!
      if(rlj.le.driwu) then
        rsi = 1.d0/max(rlj**2,asc2)
        snt = rsi*rsi*rsi
!
        ccel  = pref_LJ*am(i) * snt*(snt-0.5d0)/r2
        E_LJ  = E_LJ  + (pref_LJ*am(i)/12.d0)*(snt*(snt -1.d0) -addpot)
      end if
!
      ffr(i,1) = ffr(i,1) + (forceV +ccel)*dx
      ffr(i,2) = ffr(i,2) + (forceV +ccel)*dy
      ffr(i,3) = ffr(i,3) + (forceV +ccel)*dz
  300 end do
      end do
!
! =====================================================
!   All reduce: fec(nCLp,3), ffr(nCLp,3)
! =====================================================
!
      if(num_proc.gt.1) then
!*
        call mpi_allreduce (ffr(1,1),ffc(1,1),nCLp,mpi_real8,mpi_sum, &
                            mpi_comm_world,ierror)
        call mpi_allreduce (ffr(1,2),ffc(1,2),nCLp,mpi_real8,mpi_sum, &
                            mpi_comm_world,ierror)
        call mpi_allreduce (ffr(1,3),ffc(1,3),nCLp,mpi_real8,mpi_sum, &
                            mpi_comm_world,ierror)
!
        do i= 1,nCLp
        ffr(i,1)= ffc(i,1)
        ffr(i,2)= ffc(i,2)
        ffr(i,3)= ffc(i,3)
        end do
!
        unif1(1)= E_c_r
        unif1(2)= E_LJ
        call mpi_allreduce (unif1,unif2,2,mpi_real8,mpi_sum, &
                            MPI_COMM_WORLD,ierror)
        E_c_r= unif2(1)
        E_LJ = unif2(2)
!*
      end if
!
      return
      end subroutine cl_forces
!
!
!--------------------------------------------------------------
      subroutine init (x,y,z,px,py,pz,ch,am,ag,ipar,igrp, &
                       np,nq,nCLp,irefl)
!--------------------------------------------------------------
      use, intrinsic :: iso_c_binding
      implicit none
!
      include     'paramE7.h'
!
      integer(C_INT) ipar,igrp,np,nq,nCLp,irefl,io_pe
      common/sub_proc/ io_pe
!
      real(C_DOUBLE),dimension(npq0) :: x,y,z,px,py,pz,ch,am,ag
      real(C_DOUBLE) pi,dt,pref_CL,pref_LJ,pthe,           &
                     xmax,ymax,zmax,xmin,ymin,zmin,        &
                     xleng,yleng,zleng,gamma,vv,s1,ranff,  &
                     vmax2,dgaus2,x0,y0,z0,rr,th,ph
      integer(C_INT) i,j,k,l,ntried,k1,n1,n2
      common/parm5/  xleng,yleng,zleng
!
      common/parm2/ pi,dt,pref_CL,pref_LJ,pthe
      common/parm3/ xmax,ymax,zmax,xmin,ymin,zmin
!
      real(C_DOUBLE) R_sp,D_sp,N_sp,massi,W_sp,Lambda_D,     &
                     ch_ion,wt_ion,rd_ion,ch_el,wt_el,rd_el
      common/ionsiz/ R_sp,D_sp,N_sp,massi,Lambda_D,          &
                     ch_ion,wt_ion,rd_ion,ch_el,wt_el,rd_el 
!
      do i= 1,np
      ch(i)= ch_ion
      am(i)= wt_ion
      ag(i)= rd_ion
      end do
!
      do i= np+1,np+nq
      ch(i)= ch_el
      am(i)= wt_el
      ag(i)= rd_el
      end do
!
      if(io_pe.eq.1) then
        open (unit=11,file=praefixc//'.06'//cname,               &
              status='unknown',position='append',form='formatted')
!
        write(11,100) ch_ion,wt_ion,rd_ion
  100   format('ch_ion,wt_ion,rd_ion=',2f10.3,1pd12.3)
!
        write(11,130) ch_el,wt_el,rd_el
  130   format('ch_el, wt_el, rd_el =',2f10.3,1pd12.3)
!
        close(11)
      end if
!
!
      do i= 1,np+nq
      ntried= 0
!
  460 ntried= ntried +1
      if(ntried.gt.100000) stop
!
      x0= xmin +xleng*ranff(0.d0)
      y0= ymin +yleng*ranff(0.d0)
      z0= zmin +zleng*ranff(0.d0)
!
      x(i)= x0
      y(i)= y0
      z(i)= z0
      end do
!
      nCLp = np +nq
      vmax2= 0.7*pthe
!
      do i= 1,np
      px(i)= 0.d0
      py(i)= 0.d0
      pz(i)= 0.d0
      end do
!
      do i= np+1,np+nq
      px(i)=  am(i)*dgaus2(vmax2)
      py(i)=  am(i)*dgaus2(vmax2)
      pz(i)=  am(i)*dgaus2(vmax2)
      end do
!
      return
      end subroutine init
!
!
!-------------------------------------------------------------------
      subroutine reflect (x,y,z,vx,vy,vz,ch,ag,nCLp,irefl)
!-------------------------------------------------------------------
      use, intrinsic :: iso_c_binding
      implicit none
!
      include    'paramE7.h'
      real(C_DOUBLE),dimension(npq0) :: x,y,z,ch,vx,vy,vz,ag
      real(C_DOUBLE) xmax1,ymax1,zmax1,xmin1,ymin1,zmin1
!
      integer(C_INT) np,nq,nCLp,irefl,i,is
      real(C_DOUBLE) xmax,ymax,zmax,xmin,ymin,zmin,     &
                     pi,dt,pref_CL,pref_LJ,pthe
      common/parm2/  pi,dt,pref_CL,pref_LJ,pthe
      common/parm3/  xmax,ymax,zmax,xmin,ymin,zmin
!
      do i= 1,nCLp
      xmax1 =  xmax -ag(i)
      ymax1 =  ymax -ag(i)
      zmax1 =  zmax -ag(i)
!
      xmin1 =  xmin +ag(i)
      ymin1 =  ymin +ag(i)
      zmin1 =  zmin +ag(i)
!
      if(x(i).lt.xmin1) then
        x(i) = 2.d0*xmin1 -x(i)
        vx(i)= -vx(i)
!
      else if(x(i).gt.xmax1) then
        x(i) = 2.d0*xmax1 -x(i)
        vx(i)= -vx(i)
      end if
!***
      if(y(i).lt.ymin1) then
        y(i) = 2.d0*ymin1 -y(i)
        vy(i)= -vy(i)
!
      else if(y(i).gt.ymax1) then
        y(i) = 2.d0*ymax1 -y(i)
        vy(i)= -vy(i)
      end if
!***
      if(z(i).lt.zmin1) then
        z(i) = 2.d0*zmin1 -z(i)
        vz(i)= -vz(i)
!
      else if(z(i).gt.zmax1) then
        z(i) = 2.d0*zmax1 -z(i)
        vz(i)= -vz(i)
      end if
      end do
!
      return
      end subroutine reflect
!
!
!------------------------------------------------------------------
      subroutine ppl3da (x,y,z,ch,ag,np,nq,nCLp,irefl)
!------------------------------------------------------------------
      use, intrinsic :: iso_c_binding 
      implicit none
      include    'paramE7.h'
!
      integer(C_INT) np,nq,nCLp,irefl
!  
      real(C_DOUBLE),dimension(npq0) :: x,y,z,ch,ag
      real(C_float) ypp,zpp,t4,xp_leng,phi,tht,rgmax
      common/parm4/ phi,tht,rgmax
!
      real(C_DOUBLE) t8
      character(len=8) label,date_now*10,cax*1
      common/HEADR1/ label,date_now
      common/HEADR2/ t8,xp_leng
!
      real(C_DOUBLE) xmax,ymax,zmax,xmin,ymin,zmin
      real(C_DOUBLE) Temp,                                &
                     R_sp,D_sp,N_sp,W_sp,massi,Lambda_D,  &
                     ch_ion,wt_ion,rd_ion,ch_el,wt_el,rd_el
      common/parm3/  xmax,ymax,zmax,xmin,ymin,zmin
      common/tempr/  Temp
      common/ionsiz/ R_sp,D_sp,N_sp,massi,Lambda_D,       &
                     ch_ion,wt_ion,rd_ion,ch_el,wt_el,rd_el 
!
      integer(C_INT) io_pe,i,nskip,Nsp4
      common/sub_proc/ io_pe
      real(C_float)  cph,cth,dd,Dens4,fsize,hh,hl,pha,pi, &
                     Rsp4,sph,sth,Temp4,tha,vd,           &
                     x1,xp,xpp,xx,y1,yp,yy,z1,zp
!
      real(C_float),save :: rmax1,ps,ps2
      logical,save :: if_refl
      data  if_refl/.true./
!
!
      HH= 0.7
      CALL SYMBOL (0.5,18.0,HH,praefixc,0.,6)     !<- paramE7.h
      CALL SYMBOL (4.5,18.0,HH,praefixc//'.77'//cname,0.,10)
!
      t4= t8
      CALL SYMBOL ( 0.5,1.0,HH,date_now, 0.,10)
      CALL SYMBOL (17.0,1.0,HH,'t4=', 0.,3)
      call number (20.0,1.0,HH,t4,0.,101)
!
      Temp4 = Temp
      Dens4 = D_sp
      Rsp4  = R_sp  ! in micron
      Nsp4  = N_sp
!     CALL SYMBOL (16.4,16.0,HH,'Temp=', 0.,5)
!     call number (18.5,16.0,HH,Temp4,0.,5)
!     CALL SYMBOL (16.4,15.2,HH,'Dens=', 0.,5)
!     call number (18.5,15.2,HH,Dens4,0.,5)
!     CALL SYMBOL (16.4,14.4,HH,'Rsp=', 0.,4)
!     call number (18.5,14.4,HH,Rsp4,0.,5)
!
!     CALL SYMBOL (16.4,13.4,HH,'Nsp=', 0.,4)
!     call number (18.5,13.4,HH,float(Nsp4),0.,5)
!     CALL SYMBOL (16.4,12.6,HH,'np=', 0.,3)
!     call number (18.5,12.6,HH,float(np),0.,5)
!
      HL=  10.
      VD=   8.
!
!     phi= -60. !<-- main program
!     tht=  15.
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
      if(irefl.eq.0 .and. if_refl) then
!*
        if_refl= .false.
        fsize= 13. ! 7. ! 13.
!
        rmax1= 0.
        do i= 1,nq+np 
        rmax1= max(rmax1,real(sqrt(x(i)**2 +y(i)**2 +z(i)**2)))
        end do
!
!       ps  = 1.5*fsize/rmax1
        ps  = 0.1*fsize/rmax1 !<- may be the best
        ps2 = 1.3*fsize/rmax1
!       ps2 = 1.5*fsize/rmax1
!
        if(io_pe.eq.1) then
          open (unit=11,file=praefixc//'.06'//cname,               &
                status='unknown',position='append',form='formatted')
!
          write(11,*) 'rmax1, ps, ps2=',rmax1,ps,ps2
          close(11)
        end if
!
      else if(irefl.eq.1) then
!
        fsize= 13.
!
        rmax1= 0.
        do i= 1,nq+np 
        rmax1= max(rmax1,real(sqrt(x(i)**2 +y(i)**2 +z(i)**2)))
!       rmax1= 2.8d-2   ! fixed in time
        end do
!
        ps  = fsize/rmax1
        ps2 = 1.5*fsize/rmax1
!
        if(io_pe.eq.1) then
          open (unit=11,file=praefixc//'.06'//cname,               &
                status='unknown',position='append',form='formatted')
          write(11,*) 'rmax1=',rmax1
          close(11)
        end if
!*
      end if
!
!**********************
!*  Draw Axes.        *
!**********************
!
      CALL SYMBOL (HL+6.5,VD-3.5,HH,'Rmax=', 0.,5)
      call number (HL+9.0,VD-3.5,HH,rmax1,0.,101)
      CALL SYMBOL (HL+6.5,VD-4.3,HH,'massi=', 0.,6)
      call number (HL+9.2,VD-4.3,HH,real(massi),0.,101)
!
      do i= 1,3
      if(i.eq.1) then
         x1= 0.7*rmax1
         y1= 0.
         z1= 0.
         cax='X'
      else if(i.eq.2) then
         x1= 0.
         y1= 0.7*rmax1
         z1= 0.
         cax='Y'
      else if(i.eq.3) then
         x1= 0.
         y1= 0.
         z1= 0.7*rmax1
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
      xx= ps2*ypp  +HL
      yy= ps2*zpp  +VD
      CALL PLOT (HL,VD,3)
      CALL PLOT (xx,yy,2)
!
      CALL SYMBOL (xx-0.7,yy-0.5,HH,cax,0.,1)
      end do
!
!---    
!
      do i= 1,np,10
      xp= x(i)*cph -y(i)*sph
      yp= x(i)*sph +y(i)*cph
      zp= z(i)
!
      xpp=  xp*cth +zp*sth
      ypp= yp
      zpp= -xp*sth +zp*cth
!
      xx= ps*ypp +HL
      yy= ps*zpp +VD
!
      if(xx.lt.0.1 .or. xx.gt.23.) go to 300
      if(yy.lt.0.1 .or. yy.gt.18.) go to 300
!
      dd= 1*ps*ag(i)
!     dd= 7*ps*ag(i)
      call newcolor (3,0.,0.,1.)  ! blue
      call circle (xx-0.12,yy-0.12,dd,2)
!
  300 continue
      end do
!---
!
      do i= np+1,np+nq,10
      xp= x(i)*cph -y(i)*sph
      yp= x(i)*sph +y(i)*cph
      zp= z(i)
!
      xpp=  xp*cth +zp*sth
      ypp= yp
      zpp= -xp*sth +zp*cth
!
      xx= ps*ypp +HL
      yy= ps*zpp +VD
!
      if(xx.lt.0.1 .or. xx.gt.23.) go to 500
      if(yy.lt.0.1 .or. yy.gt.18.) go to 500
!
!     dd= 4*ps*ag(i)
      dd= 7*ps*ag(i)
      call newcolor (3,1.,0.,0.)  ! red
      call circle (xx-0.12,yy-0.12,dd,1)
!
  500 continue
      end do
!
      call newcolor (0,1.,0.,0.)  ! Reset colors
!---------------------
      CALL CHART
!---------------------
!
      return
      end subroutine ppl3da
!
!
!---------------------------------------------------------------
      subroutine ggauss
!---------------------------------------------------------------
      use, intrinsic :: iso_c_binding
      implicit none
!
      real(C_DOUBLE) fun,ranff
      real(C_DOUBLE) fv,vv0,vv,dv,s,sdv
      integer(C_INT) NS,K2,i,j,k
!
      common/gaus1/ fv(51),vv0,dv
!
      FV(1)=0.0
      VV0= -3.
      DV= 2.*ABS(VV0)/50.0
!
      VV= VV0
      DO 100 J=1,50
      S=0.0
      NS=1000
      K2=NS/2
      SDV=DV/FLOAT(NS)
!
      DO 130 K=1,K2
      VV=VV +2.0*SDV
      S=S +4.0*FUN(VV-SDV) +2.0*FUN(VV)
  130 CONTINUE
      S= (S +4.0*FUN(VV+SDV) +FUN(VV+2.0*SDV))*SDV/3.0
      FV(J+1)= FV(J)+S
  100 CONTINUE
!
      DO 200 I=1,51
      FV(I)=FV(I)/FV(51)
  200 CONTINUE
!
      return
      end subroutine ggauss
!
!
!---------------------------------------------------------------
      function dgaus2(vmax)
!---------------------------------------------------------------
      use, intrinsic :: iso_c_binding
      implicit none
!
      integer(C_INT) k,k2
      real(C_DOUBLE) dgaus2,vmax,ranff,    &
                     fv,vv0,dv,eps,x2,y1,y2
      common/gaus1/ fv(51),vv0,dv
!
      eps= ranff(0.d0)
!
      do k= 1,51
      K2=K
      IF(FV(K).GT.eps) GO TO 200
      end do
!
  200 y1= FV(K2-1)
      y2= FV(K2)
      x2= (EPS-Y2)/(Y2-Y1)+K2
!
      dgaus2= vmax*(vv0 +dv*(x2-1.0))
!
      return
      end function dgaus2
!
!
!---------------------------------------------------------------
      function fun(v)
!---------------------------------------------------------------
      use, intrinsic :: iso_c_binding
      implicit none
!
      real(C_DOUBLE) fun,v
      fun= exp(-v**2/2.)
!
      return
      end function fun
!
!
!------------------------------------------------
      block data
!------------------------------------------------
      common/ranfff/ ir,iq
      data  ir/3021/,iq/7331/
      end
!
!
!------------------------------------------------
      function ranff(x)
!------------------------------------------------
      use, intrinsic :: iso_c_binding
      implicit none
!
      real(C_DOUBLE) ranff,x,INVM
      integer(C_INT) ir,iq,MASK,LAMBDA
!
      common/ranfff/ ir,iq
      PARAMETER  (MASK=2**30+(2**30-1),INVM= 0.5D0**31)
      PARAMETER  (LAMBDA=48828125)
!
      IR= IAND( LAMBDA*IR, MASK)
      ranff= IR*INVM
!
      return
      end function ranff
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
        write(77,10) 1.-r  ! 0. for black
   10   format(f4.1,' setgray')
      end if
!
      if(ic.eq.3) then
        write(77,30) r,g,b
   30   format(3f4.1,' setrgbcolor')
      end if
!
      return
      end subroutine newcolor 
!
!
!-------------------------------------------------
      subroutine circle (x,y,d,ic)
!-------------------------------------------------
!*  Open circle centered at (x,y) /or outer edge.
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
      do 100 j= 1,nc
      th= dth*j
!
      x1= x +a*cos(th)
      y1= y +a*sin(th)
!
      call plot (x1,y1,2)
  100 continue
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
      do 300 j= -nc,nc
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
  300 continue
!
      call plot (x2,y2,3)
      write(77,*) " 1.0 setlinewidth"
!
      return
      end subroutine circle 
!
!
!***************************************************************
!*  This program package generates a UNIX postscript graphics  *
!***************************************************************
!---------------------------------------------------------------
!  PostScript header by fortran
!   Dr. T.Ogino (Nagoya University) February 27, 1992
!   Modified to conform GSIPP commands
!   Motohiko Tanaka (NIFS)       November 23, 1993
!----------------------------------------------- 5/27/1996 -----
!   This PS-Adobe-2.0 header allows us full paging features in
!  the Ghostview.  To scroll up the page (backward), click the 
!  page number and press two buttons of mouse simultaneously.
!
!---------------------------------------------------------------
       subroutine gopen (nframe)
!----------------------------------------------------------
       use, intrinsic :: iso_c_binding 

       common/convsn/ fmag,x0,y0,h0,n0
       common/pages/  ipage,nfrm
!
!*  This is an Adobe-2.0 postscript file.
!
       write(77,10)
   10  format('%!PS-Adobe-2.0',/       &
              '%%Pages: (atend)',/     &
              '%%PageOrder: Ascend',/  &
              '%%EndComments',/        &
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
             '{erasepage newpath initgraphics',/               &
             '/SaveImage save def',/                           &
             '} bind def')
!
      write(77,25) 
   25 format('/@eop          % @eop -- end a page',/ &
             '{stroke showpage',/                    &
             ' SaveImage restore',/                  &
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
             ' show grestore}',/                              &
             ' {show} ifelse',/                               &
             '} bind def')
!
      write(77,31)
   31 format('%%EndDocument',/  &
             '%%EndProlog',/    &
             '%%BeginSetup',/   &
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
!
       return
       end
!
!
!-----------------------------
       subroutine plote
!-----------------------------
       write(77,10) 
   10  format('@eop')
!
       return
       end
!
!
!-----------------------------------------
       subroutine chart
!-----------------------------------------
!*  Four frames in a page (if nfrm=4).
       common/pages/ ipage,nfrm
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
!
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
!
       return
       end
!
!
!-----------------------------
       subroutine linee
!-----------------------------
       write(77,*) 'st'
!
       return
       end
!
!
!------------------------------------
       subroutine plot (x0,y0,ip)
!------------------------------------
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
!
       return
       end
!
!
!-------------------------------------------------
       subroutine symbol (x0,y0,h0,isymb,ang,n0)
!-------------------------------------------------
       use, intrinsic :: iso_c_binding 
       implicit none
!
       character    ica*80,ich(80)*1
       character(*) isymb 
       equivalence (ica,ich(1))
!
       real(c_float)  x0,y0,h0,ang,x,y,h
       integer(c_int) n0,n,i
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
       use, intrinsic :: iso_c_binding 
       implicit none
!
       real(c_float) x0,y0,h0,anu,ang,x,y,h
       integer(c_int) n0,n,i
!      character(*)  isymb 
       character(len=9) isymb
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
       write(isymb,101) anu
  101  format(1pe9.2)
!      
       if(.true.) go to 300
       if(abs(anu).gt.1.e+1 .or.  &
          abs(anu).lt.1.e-1) then
        write(isymb,31) anu
   31   format(1pe9.2)
       else
        write(isymb,32) anu
   32   format(f7.2)
       end if
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
  300  continue
!
       write(77,*) '(',isymb,') s'
!
       return
       end
!
!
!-----------------------------------------------
       subroutine number2 (x0,y0,h0,anu,ang,n0)
!-----------------------------------------------
       use, intrinsic :: iso_c_binding  ! <-
       implicit none
!
       real(c_float)  x0,y0,h0,anu,ang,x,y,h
       integer(c_int) n0,n
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
      if(n0.eq.1) write(isymb,41) anu
      if(n0.eq.2) write(isymb,42) anu
   41  format(f6.1)
   42  format(f6.2)
!
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
