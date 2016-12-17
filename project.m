% Laurel Farris
% ASTR 545 project
% 14 December 2014

close all; clear all; clc


fid = fopen('Physical_Depth.dat','r');    
A = fscanf(fid,'%f');          
fclose(fid);
fid = fopen('Optical_Depth.dat','r');
B = fscanf(fid,'%f');
fclose(fid);
fid = fopen('Temperature.dat','r');
C = fscanf(fid,'%f');
fclose(fid);
fid = fopen('Opacity.dat','r');
D = fscanf(fid,'%f');
fclose(fid);
fid = fopen('Total_Mass_Density.dat','r');
E = fscanf(fid,'%f');
fclose(fid);
fid = fopen('Gas_Pressure.dat','r');
F = fscanf(fid,'%f');
fclose(fid);
fid = fopen('Electron_Pressure.dat','r');
G = fscanf(fid,'%f');
fclose(fid);
fid = fopen('Pressure_Ratio','r');
H = fscanf(fid,'%f');
fclose(fid);
fid = fopen('fjk11.dat','r');
I = fscanf(fid,'%f');
fclose(fid);
fid = fopen('fjk21.dat','r');
J = fscanf(fid,'%f');
fclose(fid);
fid = fopen('fjk12.dat','r');
K = fscanf(fid,'%f');
fclose(fid);
fid = fopen('fjk22.dat','r');
L = fscanf(fid,'%f');
fclose(fid);
fid = fopen('fjk32.dat','r');
M = fscanf(fid,'%f');
fclose(fid);
fid = fopen('fjk13.dat','r');
N = fscanf(fid,'%f');
fclose(fid);
fid = fopen('fjk23.dat','r');
O = fscanf(fid,'%f');
fclose(fid);
fid = fopen('fjk33.dat','r');
P = fscanf(fid,'%f');
fclose(fid);
fid = fopen('Electron_Density.dat','r');
Q = fscanf(fid,'%f');
fclose(fid);
fid = fopen('Partial.dat','r');
R = fscanf(fid,'%f');
fclose(fid);
fid = fopen('Excitation1.dat','r');
S = fscanf(fid,'%f');
fclose(fid);
fid = fopen('Excitation2.dat','r');
T = fscanf(fid,'%f');
fclose(fid);
fid = fopen('Excitation3.dat','r');
U = fscanf(fid,'%f');
fclose(fid);

 hold on;
  plot(A,B)
      
% axis([3000 27500 -70 -10]);
% 
% legend('i=1,P=1','i=1,P=10','i=1,P=100','i=2,P=1','i=2,P=10',...
%     'i=2,P=100','i=3,P=1','i=3,P=10','i=3,P=100');
% 
 xlabel('Temperature (Kelvin)');
 ylabel('log(number density (cm^{-3}))');