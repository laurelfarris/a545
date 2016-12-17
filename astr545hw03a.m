% Laurel Farris
% ASTR 545 HW03 (1)
% 24 November 2014

close all; clear all; clc

Temp = [];
A=[];B=[];C=[];D=[];E=[];F=[];G=[];H=[];I=[];

fid = fopen('temp.txt','r');    
Temp = fscanf(fid,'%f');          
fclose(fid);
fid = fopen('exc1pres1.txt','r');
A = fscanf(fid,'%f');
fclose(fid);
fid = fopen('exc1pres10.txt','r');
B = fscanf(fid,'%f');
fclose(fid);
fid = fopen('exc1pres100.txt','r');
C = fscanf(fid,'%f');
fclose(fid);
fid = fopen('exc2pres1.txt','r');
D = fscanf(fid,'%f');
fclose(fid);
fid = fopen('exc2pres10.txt','r');
E = fscanf(fid,'%f');
fclose(fid);
fid = fopen('exc2pres100.txt','r');
F = fscanf(fid,'%f');
fclose(fid);
fid = fopen('exc3pres1.txt','r');
G = fscanf(fid,'%f');
fclose(fid);
fid = fopen('exc3pres10.txt','r');
H = fscanf(fid,'%f');
fclose(fid);
fid = fopen('exc3pres100.txt','r');
I = fscanf(fid,'%f');
fclose(fid);

max(A)
max(B)
max(C)
max(D)
max(E)
max(F)
max(G)
max(H)
max(I)

hold on;
 plot(Temp,A,...
     Temp,B,...
     Temp,C,...
    Temp,D,...
    Temp,E,...
    Temp,F,...
    Temp,G,...
    Temp,H,...
    Temp,I)
axis([3000 27500 -70 -10]);

legend('i=1,P=1','i=1,P=10','i=1,P=100','i=2,P=1','i=2,P=10',...
    'i=2,P=100','i=3,P=1','i=3,P=10','i=3,P=100');

xlabel('Temperature (Kelvin)');
ylabel('log(n_{i11}/n_{1})');


