      implicit none
      integer*4 i,n,j,t,natom,iseed,nhis,ig
      parameter(n=11,natom=n*n*n,iseed=123456,nhis=750)
      character*2,ELEM(natom)
      real*8 l,x(natom),zero,y(natom),z(natom),a
      real*8 vx(natom),vy(natom),vz(natom),ranf,ekin,ekinp
      real*8 vsq,deltagv,gv(nhis),v
      real*8 xx,yy,zz,d0,epot,r,sigma,r6,eps,epotp
      real*8 deltagr,pi,gr(natom),vol,rsq,rdis
      parameter(zero=0.d0,eps=0.1d0)
      pi=4.d0*dtan(1.d0)
      l=1
C     Calculation of the positions 
      x(1)=0.d0
      do i=2, (n+1)/2
      x(i)=x(i-1)+l
      x((n+1)/2+i-1)=-x(i)
      enddo
      y(1)=0.d0
      do i=2, (n+1)/2
      y(i)=x(i-1)+l
      y((n+1)/2+i-1)=-y(i)
      enddo
      z(1)=0.d0
      do i=2, (n+1)/2
      z(i)=z(i-1)+l
      z((n+1)/2+i-1)=-z(i)
      enddo
C     Positions saved to a file 
      open(99,file='temp.dat',status='unknown')
      do i=1, n
      do j=1, n
      do t=1, n
      write(99,'(3f20.10)')x(i),y(j),z(t)
      enddo
      enddo
      enddo
      close(99)
C     Change of index i to natom 
      write(3,'(i4,/)')n**3
      open(99,file='temp.dat',status='unknown')
      do i=1,natom
      read(99,'(3f20.10)')x(i),y(i),z(i)
      if(((i/2)*2).eq.i) then
      ELEM(i)='N'
      write(3,'(a2,3f20.10)')ELEM(i),x(i),y(i),z(i)
      else
      ELEM(i)='O'
      write(3,'(a2,3f20.10)')ELEM(i),x(i),y(i),z(i)
      endif
      enddo
      close(99)
C     Assignation of random initial velocities
      call rantest(iseed)
      do i=1,natom
      vx(i)=ranf(iseed)-0.5d0
      vy(i)=ranf(iseed)-0.5d0
      vz(i)=ranf(iseed)-0.5d0
      enddo
C     Calculation of the total kinetic energy
      ekin=0.d0
      do i=1,natom
      ekin=ekin+(vx(i)**2+vy(i)**2+vz(i)**2)
      enddo
      ekin=ekin/2.d0
      write(*,*)'The total kinetic energy is', ekin, 'a.u.'
C     Kinetic energy per particle
      ekinp=ekin/natom
      write(*,*)'The kinetic energy per particle is', ekinp, 'a.u'
C     Velocity radial distribution function
      deltagv=dsqrt(3.d0)/dble(2*nhis)
      do i=1, nhis
      gv(i)=0
      enddo
      do i=1, natom
      vsq=dsqrt(vx(i)**2+vy(i)**2+vz(i)**2)
      ig=vsq/deltagv
      gv(ig)=gv(ig)+1
      enddo
      do i=1, nhis
      v=deltagv*i
      write(4,'(2f20.10)')v,gv(i)
      enddo
C     Introduction of the Lennard-Jones potential
      d0=l
      sigma=d0*(1.d0/2.d0)**(1.d0/6.d0)
      epot=0.d0
      do i=1, natom
      do j=i+1, natom-1
      xx=x(i)-x(j)
      yy=y(i)-y(j)
      zz=z(i)-z(j)
      r=dsqrt(xx**2+yy**2+zz**2)
      r6=(sigma/r)**6
      epot=epot+4.d0*eps*r6*(r6-1.d0)
      enddo
      enddo
      epotp=epot/natom
      write(*,*)'The L-J potential energy is', epot
      write(*,*)'The L-J potential energy per particle is', epotp
C     Radial distribution function of distances
      deltagr=dsqrt(3.d0)/2.d0+l
      do i=1,nhis
      gr(i)=0.d0
      enddo
      do i=1,natom
      do j=1+1,natom-1
      xx=x(i)-x(j)
      yy=y(i)-y(j)
      zz=z(i)-z(j)
      rsq=dsqrt(xx**2+yy**2+zz**2)
      ig=idint(rsq/deltagr)
      gr(ig)=gr(ig)+2
      enddo
      enddo
      do i=1,nhis
      rdis=deltagr*(i+0.5d0)
      vol=(4.d0/3.d0)*pi*((i+1)**3-i**3)*deltagr
      gr(i)=gr(i)/(natom*vol)
      enddo
      write(7,'(2f20.10)')rdis,gr 
      stop
      end
