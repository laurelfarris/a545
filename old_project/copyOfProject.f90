!Programmer(s):         Laurel Farris
!Course:                ASTR 545
!Files:                 partfunc2.txt; ionpots.txt; amu.txt; opacity.txt;
!                       mihalas.txt (from notes-to calculate q)
!Functions:             zbrent.f90; interp.f90
!Last modified:         17 December 2014
!Description:           Stellar Atmosphere Model for Sirius A type star


!--------------------------------Variables------------------------------------!
module stuff
double precision,dimension(11) :: HM,H,HP,He,HeP,HePP,Ca,CaP,CaPP,theta
double precision               :: othertheta,d_opd
integer                        :: i,j,k,bool,counter

double precision,parameter :: constant=4.83d15
double precision,parameter :: boltzmann=1.3806488d-16  !Boltzmann's const
double precision,parameter :: Mamu=1.66054d-24         !grams per amu
double precision,parameter :: mass_electron=9.1d-28    !electron mass (grams)
double precision,parameter :: grav_const=6.67d-08      !Gravitational constant

integer,parameter                 :: array=45
double precision                  :: Teff=9940.0            !temperature
double precision                  :: radius=1.1905138d12    !Stellar radius
double precision                  :: mass=4.04d33           !Stellar mass
double precision,dimension(array) :: P_gas,opd,temp,X_final,depth
double precision                  :: e=0.001,P_est1,P_est2

integer,parameter                :: species=3           !#elements
integer,parameter                :: ionization=3        !#ionization levels
double precision,parameter       :: pf_num_columns=11.0 !#columns in pf table
double precision,parameter       :: opd_min=1.0d-3 !minimum optical depth (opd)
double precision,parameter       :: opd_max=1.0d+2 !maximum optical depth (opd) 
double precision,dimension(20)   :: opd_in,q_in,dummy
double precision                 :: q_out,R_calc
double precision,dimension(20,71):: opacity_data,X_calc
double precision,dimension(19)   :: R_table
double precision,dimension(70)   :: T_table
double precision,dimension(19,70):: X_table
double precision,parameter       :: tol=1.0d-5      !electron density tolerance
double precision    :: ne_min=1.0d-4,ne_max=1.0d15,ne_calc
double precision    :: totalAbundance
double precision    :: summation, summation2

double precision                                :: g
double precision,dimension(ionization,species)  :: Xion      !ionization potential
double precision,dimension(species)             :: A,X,alpha !abundance/mass fraction
double precision,dimension(ionization,species)  :: pf,phi    !part. func. and phi
double precision                                :: error     !part. func. error

double precision    :: nTotal,nNuclear,nElectron !number densities
double precision    :: uTotal,uNuclear,uElectron !mean molecular weight
double precision    :: pTotal,pNuclear,pElectron !mass densities
double precision    :: P_electron,P_nuclear,P_ratio !Pressures
double precision    :: exc1,exc2,exc3,partial
double precision,dimension(ionization,species) :: Y,njk,fjk,fne !Y and fractions
double precision,dimension(species)   :: nk,ionStage  !number density of each species

end module stuff
!-----------------------------------------------------------------------------!
program ionization

use stuff
implicit none
external func
!double precision :: func, zbrent
real*8 :: func,zbrent

!Open files for writing data (one for each required plot)
open(unit=11,file='plot1.txt')
open(unit=12,file='plot2.txt')
open(unit=13,file='plot3.txt')
open(unit=14,file='plot4.txt')
open(unit=15,file='plot5.txt')

open(unit=1,file='mihalas.txt',status='old',action='read')
read(1,*) opd_in
read(1,*) q_in
close(unit=1)

!Read in values from opacity table
open(unit=1,file='opacity.txt',status='old',action='read')
read(1,*) opacity_data
R_table=opacity_data(2:20,1)
T_table=opacity_data(1,2:71)
X_table=opacity_data(2:20,2:71)
close(unit=1)

!Set up tau array
opd(1)=opd_min
d_opd=1.0d-3
counter=1
do while (counter <= array)
counter=counter+1
opd(counter)=opd(counter-1)+d_opd 
 if (mod(counter-1,10)==0) then
  d_opd=d_opd*10
 end if
end do

g=grav_const*mass/(radius**2)
P_gas(1)=10**(1.6+(2/3)*(log(g)-4.0))
temp(1)=Teff
depth(1)=0.0
counter=1
bool=1
!----------------------Begin optical depth loop-------------------------------!
!do while (opd(counter)<opd_max)

100 do while (counter <= array)
200 FORMAT (30ES12.3E3)


!interpolate tau_in and q_in from table to get q_out to use for 
!temperature as a function of our tau

  CALL SPLINE(opd_in,q_in,20,2.0E30,2.0E30,dummy)
  CALL SPLINT(opd_in,q_in,dummy,20,opd(counter),q_out)

!Use q output to calculate temp as a function of tau
!That thing drifted
temp(counter)=(0.75*(Teff**4.0)*(opd(counter)+q_out))**0.25

!Solve for ne(opd) and p(opd) (#density and mass density)
nElectron = zbrent(func,ne_min,ne_max,tol)
nTotal = P_gas(counter)/(boltzmann*temp(counter))
nNuclear = nTotal - nElectron

!Electron pressure and pressure ratio
P_electron = nElectron*boltzmann*temp(counter)
P_nuclear = P_gas(counter) - P_electron
P_ratio = P_electron/P_gas(counter)

!Molecular weights
uNuclear = 1/((X(1)/A(1))+(X(2)/A(2))+(X(3)/A(3)))
uElectron = 1/summation2
uTotal = 1/((1/uNuclear)+(1/uElectron))

!Mass densities
pElectron = mass_electron*nElectron
pNuclear = (P_nuclear*uNuclear*Mamu)/(boltzmann*temp(counter))
pTotal = pElectron+pNuclear

!Calculate number densities and ionization fractions
do k=1,species
 nk(k)=alpha(k)*nNuclear
 do j=1,ionStage(k) 
  njk(j,k) = fjk(j,k)*nk(k)
  fne(j,k) = njk(j,k)*(j-1)/nElectron
 end do
end do

!Calculate R for opacity interpolation
R_calc=pTotal/((temp(counter)/(10**6))**3)

!Solve for opacity by interpolating table
CALL SPLIE2(R_table,T_table,X_table,19,70,X_calc)
CALL SPLIN2(R_table,T_table,X_table,X_calc,19,70,R_calc,temp(counter),X_final(counter))

!calculate physical depth
 if (counter/=1) then
  depth(counter)=depth(counter-1)+ &
                ((opd(counter)-opd(counter-1))/((X_final(counter))*pTotal))
 end if

!------ layer has been solved here for new pressure-------!

!Write data to files, increment loop, and estimate new pressure
if (bool==1) then
 write(11,200) opd(counter),depth(counter),temp(counter),X_final(counter),pTotal
 write(12,200) opd(counter),P_gas(counter),P_electron,P_ratio
 write(13,200) opd(counter),fjk(1,1),fjk(2,1),fjk(1,2),fjk(2,2), &
                 fjk(3,2),fjk(1,3),fjk(2,3),fjk(3,3)
 write(14,200) opd(counter),pElectron!,partial
 !write(15,200) opd,exc1,exc2,exc3
 counter=counter+1

 P_est1=(P_gas(counter-1))+g* &
      ((opd(counter)-opd(counter-1))/(X_final(counter-1)))
 P_gas(counter)=P_est1
 bool=0
 goto 100
end if

! pressure estimates
 P_est2=(P_gas(counter-1))+2*g*( (opd(counter)-opd(counter-1))/ &
                              (X_final(counter-1)+X_final(counter)))
 if (abs((P_est2-P_est1)/P_est1)<=e) then
  bool=1
 end if
 P_gas(counter)=P_est2
 goto 100

end do

end program ionization

!---------------Function for calculating electron density------------------!

double precision function func(nElec)
use stuff
real*8 :: nElec

A(1)=1.00780;A(2)=4.00260;A(3)=40.08000 !atomic weight A_k
X(1)=0.70;X(2)=0.28;X(3)=0.02
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

othertheta=(5040.0/temp(counter))
CALL POLINT(theta,H,pf_num_columns,othertheta,pf(1,1),error)
CALL POLINT(theta,HP,pf_num_columns,othertheta,pf(2,1),error)
CALL POLINT(theta,He,pf_num_columns,othertheta,pf(1,2),error)
CALL POLINT(theta,HeP,pf_num_columns,othertheta,pf(2,2),error)
CALL POLINT(theta,HePP,pf_num_columns,othertheta,pf(3,2),error)
CALL POLINT(theta,Ca,pf_num_columns,othertheta,pf(1,3),error)
CALL POLINT(theta,CaP,pf_num_columns,othertheta,pf(2,3),error)
CALL POLINT(theta,CaPP,pf_num_columns,othertheta,pf(3,3),error)

! Actual values of partition functions
do k=1,species
  do j=1,ionStage(k)
     pf(j,k)=10**pf(j,k)
  end do
end do

! Phi
do k=1,species
 do j=1,ionStage(k)
   phi(j,k) = constant*(temp(counter)**(3/2))*(pf(j+1,k)/pf(j,k))* &
                   exp(-(Xion(j,k))/(boltzmann*temp(counter)))
 end do
end do

! Y
do k=1,species
 do j=1,ionStage(k)
   Y(j,k)=(nElec**(-1))*phi(j,k)
 end do
end do

! Total Abundance
totalAbundance=0
do k=1,species
 totalAbundance=totalAbundance+(X(k)/A(k))
end do

! ionization fraction
do k=1,species
 do j=1,ionStage(k)
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

ne_calc=(P_gas(counter)/(boltzmann*temp(counter)))*summation
func=nElec-ne_calc

return
end
!---------------------------------------------------------------------------!
