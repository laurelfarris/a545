!Programmer(s):         Laurel Farris
!Course:                ASTR 545; hw03(1)
!Files:                 2partfunc.txt (removed first column for baby fortran)
!Functions:             polint.f90
!Last modified:         26 November 2014
!Description:           This program calculates the fraction of neutral 
!                       Hydrogen (H) at the first (ground), second, and 
!                       third excitation states to all H 
!                       in a gas structure in thermal equilibrium

program excitation

real*8                  :: k=1.3806488e-16 !Boltzmann const (erg)
real*8                  :: R=13.598*1.6e-12 !Rydberg const (erg)
real*8                  :: C=4.83e15  !Constant from fraction
real*8                  :: P,T,Tmin=3000.0,Tmax=20000.0,dT=1
real*8                  :: frac,pf11,pf21,e_pf11,e_pf21,num_columns=11.0
real,dimension(11)      :: Hminus,Hneutral,Hplus,theta,othertheta
real*8                  :: Xion_H=13.5984340*1.6e-12 !ionization pot (erg)
real*8                  :: g, Xexc_H

!Read in partfunc values
open(unit=1,file='2partfunc.txt',status='old',action='read')
read(1,*) theta
read(1,*) Hminus
read(1,*) Hneutral
read(1,*) Hplus
close(unit=1)

!Open files for plotting
open(unit=2,file='temp.txt')
open(unit=3,file='exc1pres1.txt')
open(unit=4,file='exc1pres10.txt')
open(unit=5,file='exc1pres100.txt')
open(unit=6,file='exc2pres1.txt')
open(unit=7,file='exc2pres10.txt')
open(unit=8,file='exc2pres100.txt')
open(unit=9,file='exc3pres1.txt')
open(unit=10,file='exc3pres10.txt')
open(unit=11,file='exc3pres100.txt')

!Write temperature increments to a file for plotting
T=Tmin
do while (T<Tmax)
 write(2,*) T
 T = T+dt
end do

!Compute fraction over the range T=3000K to T=20000K

do i=1,3         
  g=2*i**2      
  Xexc_H=R*(1.0-(1.0/(i**2)))

  do j=1,3      
    P=10.0**(j-1)
    T=Tmin

    do while (T<Tmax)
      othertheta=5040/T
      CALL POLINT(theta,Hneutral,num_columns,othertheta,pf11,e_pf11)
      CALL POLINT(theta,Hplus,num_columns,othertheta,pf21,e_pf21)

       frac = ((g/(10**pf11))*exp(-Xion_H/(k*T)))/ &
          (1+((k*T*C*T**(3/2)*(10**pf21))/(P*(10**pf11)))* &
                       exp(-Xexc_H/(k*T)))
      frac=log(frac)
      if (i==1 .and. P==1) then
      write(3,*) frac
      else if (i==1 .and. P==10) then
      write(4,*) frac
      else if (i==1 .and. P==100) then
      write(5,*) frac
      else if (i==2 .and. P==1) then
      write(6,*) frac
      else if (i==2 .and. P==10) then
      write(7,*) frac
      else if (i==2 .and. P==100) then
      write(8,*) frac
      else if (i==3 .and. P==1) then
      write(9,*) frac
      else if (i==3 .and. P==10) then
      write(10,*) frac
      else if (i==3 .and. P==100) then
      write(11,*) frac
      end if
     T = T + dT
    end do
  end do
end do

close(unit=2)
close(unit=3)
end program excitation


