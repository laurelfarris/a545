% Laurel Farris
% ASTR 545 HW03 (2,3)
% 12 December 2014

close all; clear all; clc

Temp = [];
A=[];B=[];C=[];D=[];E=[];F=[];G=[];H=[];I=[];J=[];K=[];

fid = fopen('temperature.txt','r');    
Temp = fscanf(fid,'%f');          
fclose(fid);
fid = fopen('n11.txt','r');
A = fscanf(fid,'%f');
fclose(fid);
fid = fopen('n21.txt','r');
B = fscanf(fid,'%f');
fclose(fid);
fid = fopen('n12.txt','r');
C = fscanf(fid,'%f');
fclose(fid);
fid = fopen('n22.txt','r');
D = fscanf(fid,'%f');
fclose(fid);
fid = fopen('n32.txt','r');
E = fscanf(fid,'%f');
fclose(fid);
fid = fopen('n13.txt','r');
F = fscanf(fid,'%f');
fclose(fid);
fid = fopen('n23.txt','r');
G = fscanf(fid,'%f');
fclose(fid);
fid = fopen('n33.txt','r');
H = fscanf(fid,'%f');
fclose(fid);
fid = fopen('n1.txt','r');
I = fscanf(fid,'%f');
fclose(fid);
fid = fopen('n2.txt','r');
J = fscanf(fid,'%f');
fclose(fid);
fid = fopen('n3.txt','r');
K = fscanf(fid,'%f');
fclose(fid);

 hold on;
  plot(Temp,log(A),Temp,log(B),Temp,log(C),Temp,log(D),Temp,log(E),...
      Temp,log(F),Temp,log(G),Temp,log(H),Temp,...
      log(I),Temp,log(J),Temp,log(K));
%       Temp,log(G),...
%       Temp,log(H))%,...
 %    Temp,log(D),...
 %    Temp,log(E));
     %Temp,F,...
%     Temp,G,...
%     Temp,H,...
%     Temp,I)
% axis([3000 27500 -70 -10]);
% 
 legend('n_{11}','n_{21}','n_{12}','n_{22}',...
          'n_{32}','n_{13}','n_{23}','n_{33}',...
           'n_{1}','n_{2}','n_{3}',-1);
% 
 xlabel('Temperature (Kelvin)');
 ylabel('Number density (cm^{-3}');
