clear all;
close all;
clc;
%% Step 1: Intial Data
nd=1; %dielectric constant
lambda=1; %wavelength
e=1; %initial thickness
M=100; %discretiziation in azimuth
N=100; %discretization in elevation
Ef=zeros(2*M+1,N+1); %Far field
g=zeros(2*M+1,N+1); %Primary Feed radiation intensity
h=zeros(2*M+1,N+1); %Desired intensity
r=e*ones(2*M+1,N+1); %thickness matrix
alpham=zeros(2*M+1,1); %Maximum angle of refraction along azimuth
%% Step 2: Initialization
k=0; %counter
th=1; %threshold
e=10; %error
n=zeros(2*M+1,N+1,3); %normal vector
t=zeros(2*M+1,N+1,3); %transverse vector
ki=zeros(2*M+1,N+1,3); %incident wave vector
kt=zeros(2*M+1,N+1,3); %transmitted wave vector
theta=0:(pi/2*N):pi/2; %discretized elevation angles 
phi=-pi:(pi/M):pi; %discretized azimuthal angles
thetai=zeros(2*M+1,N+1); %discretized incident angles 
thetat=zeros(2*M+1,N+1); %discretized transmitted angles 
while(e>th)
%% Step 3: Characterization of Plane of Incidence
for i=1:N+1
    for j=1:2*M+1
        %Initialise n for each loop
        n(i,j,1)=-cos(theta(i))*cos(phi(j));
        n(i,j,2)=-cos(theta(i))*sin(phi(j));
        n(i,j,3)=-sin(theta(i));
        %Initialise ki for each loop
        ki(i,j,1)=r*cos(theta(i))*cos(phi(j));
        ki(i,j,2)=r*cos(theta(i))*sin(phi(j));
        ki(i,j,3)=r*sin(theta(i));
        
        t(i,j,:)=cross(-n(i,j,:),(cross(n(i,j,:),ki(i,j,:)))); %transverse vector
        t(i,j,:)=t(i,j,:)/norm(t(i,j,:));
        
        thetai(i,j)=atan(dot(ki(i,j,:),t(i,j,:))/dot(ki(i,j,:),n(i,j,:)));
        thetat(i,j)=asin(nd*sin(thetai(i,j,:)));
        
        kt(i,j,:)=cos(thetat(i,j))*n(i,j,:)+sin(thetat(i,j))*t(i,j,:);
    end
end
%%
end