module percroutines
INTEGER(KIND=2),    PARAMETER     :: DBL = 8, LONG = 4
!--2D
character(*),       parameter     :: name = "perc2D_f40_dyn1e5"
integer(kind=long), parameter     :: dim  = 2
!--3D
!integer(kind=long), parameter     :: dim  = 3
!character(*),       parameter     :: name = "perc3D_f68_dyn1e4"

integer(kind=long), parameter     :: &
     N    = 10,&                       !--> # tracers, fixed
     ze   = 5d6,&                      !--> trajectory length, old 5d4
     Mtra = 1000,&                     !--> # of trajectories in total
     Ntra = ze/100,&                   !--> time steps recorded in trajectory
     Nxy  = 300+(3-dim)*100,&          !--> lattice size x,y; 400 in 2D; 300 in 3D
     Nz   = 1+(dim-2)*(Nxy-1),&        !--> lattice size z
     Ntot = Nxy*Nxy*Nz                 !--> # lattice points

real(kind=dbl),     parameter     :: &
     dt        = 4*6.25d-4/(1d0*dim),& !--> time increment, 100% step probab.
     dx        = 0.1d0,&               !--> spatial discretization 100nm
     Dconst    = 2*dim*dt*2d0/dx**2,&  !--> diffusion coefficient D=2
     f_obst    = 0.4d0+(dim-2)*0.28d0  !--> fraction of obstacles ca. 40.7% in 2D, 68.8% in 3D

!------------------------------------------------------------------------------

integer(kind=long), dimension(N)                  :: p
integer(kind=long), dimension(Ntot)               :: g
integer(kind=long), dimension(Ntot,6)             :: neigh
integer(kind=long), dimension(6,3)                :: inc
integer(kind=long), dimension(Nxy,Nxy,Nz)         :: liste
integer(kind=long), dimension(:), allocatable     :: orb
real(kind=dbl),     dimension(N,ze,3)             :: pos
integer(kind=long)                                :: Nobst,seed

contains
!##################################################################

!-----------------------------------------------------------------------
! initialize system
!-----------------------------------------------------------------------

subroutine initialize
implicit none
integer(kind=long)               :: i,j,k,nn,nr,xl,yl,zl,xr,yr,zr

call initrnd
if (Dconst.le.0.999d0) then
   print*,"probability failure",Dconst,dt
   stop
endif

Nobst      = Ntot*f_obst
inc(1,1:3) = (/ -1, 0, 0 /)
inc(2,1:3) = (/ +1, 0, 0 /)
inc(3,1:3) = (/  0,-1, 0 /)
inc(4,1:3) = (/  0,+1, 0 /)
inc(5,1:3) = (/  0, 0,-1 /)
inc(6,1:3) = (/  0, 0,+1 /)

!---> ijk<->nr 
nr=0
do i=1,Nxy
   do j=1,Nxy
      do k=1,Nz
         nr            = nr+1
         liste(i,j,k)  = nr
      enddo
   enddo
enddo

!---> NN
nr=0
do i=1,Nxy
   xl = i-1
   xr = i+1
   if (i.eq.1)   xl=Nxy
   if (i.eq.Nxy) xr=1
   do j=1,Nxy
      yl = j-1
      yr = j+1
      if (j.eq.1)   yl=Nxy
      if (j.eq.Nxy) yr=1
      do k=1,Nz
         zl = k-1
         zr = k+1
         if (k.eq.1)  zl=Nz
         if (k.eq.Nz) zr=1
         nr            = nr+1
         neigh(nr,1:6) = (/ liste(xl,j,k),liste(xr,j,k),&
                            liste(i,yl,k),liste(i,yr,k),&
                            liste(i,j,zl),liste(i,j,zr)/)
      enddo
   enddo
enddo
allocate(orb(Nobst))
end subroutine initialize

!-----------------------------------------------------------------------
! create percolation cluster and distribute tracers
!-----------------------------------------------------------------------

subroutine cluster
implicit none
integer(kind=long)                     :: i,j,k,nn,nr
real(kind=dbl)                         :: mranf

!--- position obstacles
g(:) = 0
nr   = 0
do 
   i  = Nxy*mranf()+1
   j  = Nxy*mranf()+1
   k  = Nz*mranf() +1
   nn = liste(i,j,k)
   if (g(nn).eq.0) then
      nr      = nr+1
      orb(nr) = nn
      g(nn)   = 1      
   endif
   if (nr.eq.Nobst) exit
enddo

!---> position tracer particles, invisible to each other!
pos(:,:,:) = 0
nr         = 0
do 
   i  = Nxy*mranf()+1
   j  = Nxy*mranf()+1
   k  = Nz*mranf() +1
   nn = liste(i,j,k)
   if (g(nn).eq.0) then
      nr    = nr+1
      p(nr) = nn
   endif
   if (nr.eq.N) exit
enddo
end subroutine cluster

!-----------------------------------------------------------------------
! diffusion of tracer, blind ant
!-----------------------------------------------------------------------

subroutine diffuse(t)
implicit none
integer(kind=long), intent(in)             :: t
integer(kind=long)                         :: i,j,k
real(kind=dbl)                             :: mranf

do i=1,N
   j            = 2*dim*mranf()+1
   k            = neigh(p(i),j)
   pos(i,t,1:3) = pos(i,t-1,1:3)
   if (g(k).eq.0) then
      p(i) = k
      pos(i,t,1:3) = pos(i,t-1,1:3)+inc(j,1:3)
   endif
enddo
end subroutine diffuse

!-----------------------------------------------------------------------
! diffusion of obstacles (blind ant, no tracer hit)
!-----------------------------------------------------------------------

subroutine dynobst
implicit none
integer(kind=long)                         :: i,j,k
real(kind=dbl)                             :: mranf

do i=1,Nobst 
   j            = 2*dim*mranf()+1
   k            = neigh(orb(i),j)
   if ((g(k).eq.0).and.(minval(abs(p(1:10)-k)).gt.0)) then
      g(orb(i)) = 0
      g(k)      = 1      
      orb(i)    = k
   endif
enddo
end subroutine dynobst

!--------------------------------------------------------------------
! initialize rnd-generator
!--------------------------------------------------------------------

subroutine initrnd
implicit none
integer(kind=long),dimension(8)               :: t_init

call date_and_time(values=t_init)
seed =t_init(1)+70*(t_init(2)+12*(t_init(3)+31*(t_init(5)+23*(t_init(6)+59*t_init(7)))))
call init_genrand(seed)
end subroutine initrnd

!##################################################################

end module percroutines

