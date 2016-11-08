!Programmer(s):         Laurel Farris
!Course:                ASTR 545 hw03_2
!Files:                 partfunc2.txt; ionpots.txt; amu.txt
!Functions:             zbrent.f90 (zbrent.f); polint.f90 (ionpots.tab)
!Last modified:         11 December 2014
!Description:           Particle/Charge conservation

!program ionization

!--------------------------------Variables------------------------------------!
module stuff
real,dimension(11) :: HM,H,HP,He,HeP,HePP,Ca,CaP,CaPP,theta
real*8             :: othertheta,dt=10.0

real*8,parameter    :: constant=4.83d15
real*8,parameter    :: boltzmann=1.3806488d-16  !Boltzmann's const
real*8,parameter    :: boltzmannEV=8.6173324d-5
real*8,parameter    :: Mamu=1.66054d-24         !grams per amu
real*8,parameter    :: mass_electron = 9.1d-28  !electron mass (grams)

real*8,parameter    :: species=3,ionization=3   !3 elements and 3 ion. stages
real*8,parameter    :: num_columns=11.0         !# columns in pf table
real*8,parameter    :: Tmin=1500.0,Tmax=25000.0 !temperature range (part 3) 
real*8              :: temp=6000.0              !******** part 2            
real*8,parameter    :: P_gas=100.0              !pressure
real*8,parameter    :: tol=1d-5                 !electron density tolerance
real*8              :: ne_min=0.0001,ne_max=1.0d15,ne_calc,totalAbundance
real*8              :: summation, summation2

real,dimension(ionization,species)  :: Xion     !ionization potential
real,dimension(species)             :: A,X,alpha  !abundance/mass fraction
real,dimension(ionization,species)  :: pf,phi   !part. func. and phi
real*8                              :: error      !part. func. error

real*8              :: nTotal,nNuclear,nElectron !number densities
real*8              :: uTotal,uNuclear,uElectron !mean molecular weight
real*8              :: pTotal,pNuclear,pElectron !mass densities
real*8              :: P_electron,P_nuclear,P_ratio !Pressures
real,dimension(ionization,species) :: Y,njk,fjk,fne   !Y and fractions
real,dimension(species)   :: nk,ionStage        !number density of each species

real*8 :: num,num2                           !counters for program   
end module stuff
!-----------------------------------------------------------------------------!
program ionization

use stuff
implicit none
external func
real*8 :: func,zbrent

!Open files for writing data
open(unit=20,file='hw03_2_table2.txt')
open(unit=21,file='hw03_3_table.txt')

!--------Begin temperature loop for part 03-----------!

!temp=Tmin                    !********* part 3
!do while (temp<Tmax)         !********* part 3

!Number densities
nElectron = zbrent(func,ne_min,ne_max,tol)
print*, nElectron
nTotal = P_gas/(boltzmann*temp)
nNuclear = nTotal - nElectron

!Molecular weights
uNuclear = 1/((X(1)/A(1))+(X(2)/A(2))+(X(3)/A(3)))
uElectron = 1/summation2
uTotal = 1/((1/uNuclear)+(1/uElectron))

!Electron pressure and pressure ratio
P_electron = nElectron*boltzmann*temp
P_nuclear = P_gas - P_electron
P_ratio = P_electron/P_gas

!Mass densities
pElectron = mass_electron*nElectron
pNuclear = P_nuclear*uNuclear*Mamu/(boltzmann*temp)
pTotal = pElectron+pnuclear

!Fortran doesn't like k's and j's used in multiple modules
do num=1,species
 nk(num)=alpha(num)*nNuclear
 do num2=1,ionStage(num) 
  njk(num2,num) = fjk(num2,num)*nk(num)
  fne(num2,num) = njk(num2,num)*(num2-1)/nElectron
 end do
end do

!Write all data to the same file and hope it works
write(20,*) temp,nElectron,nNuclear,nTotal,nk(1),nk(2),nk(3),&
            njk(1,1),njk(1,2),njk(2,1),njk(2,2),njk(2,3),njk(1,3),njk(2,3),njk(3,3),&
            pElectron,pNuclear,pTotal,P_ratio,&
            fjk(1,1),fjk(2,1),fjk(1,2),fjk(2,2),fjk(3,2),fjk(1,3),fjk(2,3),fjk(3,3),&
            fne(1,1),fne(1,2),fne(2,1),fne(2,2),fne(3,2),fne(1,3),fne(2,3),fne(3,3)
             
!temp=temp+dt   !******* part 3
!end do         !******* part 3

close(unit=21)
close(unit=22)

end program ionization
!------------------------------------------------------------------------------!
double precision function func(nElec)
use stuff
real*8 :: nElec

A(1)=1.00780;A(2)=4.00260;A(3)=40.08000 !atomic weight A_k
X(1)=0.71;X(2)=0.27;X(3)=0.02
Xion(1,1)=13.5984340*1.6e-12 !Hydrogen (erg)
Xion(1,2)=24.5873876*1.6e-12 !Helium (erg)
Xion(2,2)=54.4177630*1.6e-12
Xion(1,3)=6.113158*1.6e-12 !Calcium (erg)
Xion(2,3)=11.87172*1.6e-12 

ionStage(1)=2; ionStage(2)=3; ionStage(3)=3 !Assign ion stage for each species

open(unit=1,file='partfunc2.txt',status='old',action='read')
read(1,*) theta
read(1,*) HM
read(1,*) H
read(1,*) HP
read(1,*) He
read(1,*) HeP
read(1,*) HePP
read(1,*) Ca
read(1,*) CaP
read(1,*) CaPP
close(unit=1)

othertheta=(5040.0/temp)
CALL POLINT(theta,H,num_columns,othertheta,pf(1,1),error)
CALL POLINT(theta,HP,num_columns,othertheta,pf(2,1),error)
CALL POLINT(theta,He,num_columns,othertheta,pf(1,2),error)
CALL POLINT(theta,HeP,num_columns,othertheta,pf(2,2),error)
CALL POLINT(theta,HePP,num_columns,othertheta,pf(3,2),error)
CALL POLINT(theta,Ca,num_columns,othertheta,pf(1,3),error)
CALL POLINT(theta,CaP,num_columns,othertheta,pf(2,3),error)
CALL POLINT(theta,CaPP,num_columns,othertheta,pf(3,3),error)

! Actual values of partition functions
do k=1,species
  do j=1,ionStage(k)
     pf(j,k)=10**pf(j,k)
  end do
end do

! Phi
do k=1,species
 do j=1,ionStage(k)
   phi(j,k) = constant*(temp**(3/2))*(pf((j+1),k)/pf(j,k))*exp(-(Xion_H)/(boltzmann*temp))
 end do
end do

! Y
!do k=1,species
! do j=1,(k+1)
!  if (array(j,k)/=0) then
!   Y(j,k)=(nElec**(-1))*phi(j,k)
!  end if 
! end do
!end do

! Total Abundance
totalAbundance=0
do k=1,species
 totalAbundance=totalAbundance+(X(k)/A(k))
end do

! ionization fraction
do k=1,species
 do j=1,ionStage(k)
   Y(j,k)=(nElec**(-1))*phi(j,k)
   if(j==1) then
    fjk(j,k)=1/(1+Y(j,k)+Y(j+1,k)*Y(j,k))
   else
    fjk(j,k)=fjk(j-1,k)*Y(j-1,k)
   end if 
 end do
 alpha(k)=(X(k)/A(k))/totalAbundance
end do
 
summation=0
do k=1,species
 do j=1,ionStage(k)
   summation = summation + alpha(k)*(j-1)*fjk(j,k) 
 end do
end do
summation2 = summation*totalAbundance ! for uElectron

ne_calc=(P_gas/(boltzmann*temp))*summation
func=nElec-ne_calc

return
end
!---------------------------------------------------------------------------!
