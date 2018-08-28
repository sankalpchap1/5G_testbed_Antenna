clear all;
close all;
clc;
%% Step 1: Intial Data
nd=1; %dielectric constant
lambda=1; %wavelength
e=1; %initial thickness
M=100; %discretiziation in azimuth
N=100; %discretization in elevation
% needs to be initialised
Ef=zeros(N+1,2*M+1,3); %Far field
g=zeros(N+1,2*M+1); %Primary Feed radiation intensity
h=zeros(N+1,2*M+1); %Desired intensity
alpham=zeros(2*M+1,1); %Maximum angle of refraction along azimuth
%initialized
r=e*ones(N+1,2*M+1); %thickness matrix
delta=1e-3; %Loss angle of dielectric
kd=(2*pi*nd/lambda)*(1-(1i*delta/2)); %wave number 
Ei=zeros(N+1,2*M+1,3); %Incident field
Ein=zeros(N+1,2*M+1); %Magnitude of Perpendicular component of Incident field
Eip=zeros(N+1,2*M+1); %Magnitude of Parallel component of Incident field
%% Step 2: Initialization
k=0; %counter
th=1; %threshold
e=10; %error
n=zeros(N+1,2*M+1,3); %normal vector
t=zeros(N+1,2*M+1,3); %transverse vector
ki=zeros(N+1,2*M+1,3); %incident wave vector
kt=zeros(N+1,2*M+1,3); %transmitted wave vector
deltheta=pi/(2*N); %Each discretisation in theta
delphi=pi/M; %Each discretisation in phi
theta=0:deltheta:pi/2; %discretized elevation angles 
phi=-pi:delphi:pi; %discretized azimuthal angles
thetai=zeros(N+1,2*M+1); %discretized incident angles 
thetat=zeros(N+1,2*M+1); %discretized transmitted angles 
Rn=zeros(N+1,2*M+1); %Perpendicular Fresnel Coefficients
Rp=zeros(N+1,2*M+1); %Parallel Fresnel Coefficients
Tk=zeros(N+1,2*M+1); %Power Transmission Coefficient
alphak=zeros(N+1,2*M+1); %direction of refracted rays in kth iteration
betak=zeros(N+1,2*M+1);
% Declaring partial derivatives of rk
rt=zeros(N+1,2*M+1); %along theta
rp=zeros(N+1,2*M+1); %along phi
rtt=zeros(N+1,2*M+1); %along theta twice
rpp=zeros(N+1,2*M+1); %along phi twice
rtp=zeros(N+1,2*M+1); %along theta and then phi

while(e>th)
%% Step 3: Characterization of Plane of Incidence
%Defining partial derivatives for edge cases

    for j=2:2*M
        % i=1
        rt(1,j)=(r(2,j)-r(1,j))/(2*deltheta);
        rp(1,j)=(r(1,j+1)-r(1,j-1))/(2*delphi);
        rtt(1,j)=(r(2,j)-2*r(1,j)+r(1,j))/(deltheta^2);
        rpp(1,j)=(r(1,j+1)-2*r(1,j)+r(1,j-1))/(delphi^2);
        rtp(1,j)=(r(2,j+1)-r(2,j-1)-r(1,j+1)+r(1,j-1))/(4*deltheta*delphi);
        
        % i=N+1
        rt(N+1,j)=(r(N+1,j)-r(N,j))/(2*deltheta);
        rp(N+1,j)=(r(N+1,j+1)-r(N+1,j-1))/(2*delphi);
        rtt(N+1,j)=(r(N+1,j)-2*r(N+1,j)+r(N,j))/(deltheta^2);
        rpp(N+1,j)=(r(N+1,j+1)-2*r(N+1,j)+r(N+1,j-1))/(delphi^2);
        rtp(N+1,j)=(r(N+1,j+1)-r(N+1,j-1)-r(N,j+1)+r(N,j-1))/(4*deltheta*delphi);
    end
    
    for i=2:N
        %j=1
        rt(i,1)=(r(i+1,1)-r(i-1,1))/(2*deltheta);
        rp(i,j)=(r(i,2)-r(i,1))/(2*delphi);
        rtt(i,j)=(r(i+1,1)-2*r(i,1)+r(i-1,1))/(deltheta^2);
        rpp(i,j)=(r(i,2)-2*r(i,1)+r(i,1))/(delphi^2);
        rtp(i,j)=(r(i+1,2)-r(i+1,1)-r(i-1,2)+r(i-1,1))/(4*deltheta*delphi);
        
        %j=2*M+1
        rt(i,2*M+1)=(r(i+1,2*M+1)-r(i-1,2*M+1))/(2*deltheta);
        rp(i,2*M+1)=(r(i,2*M+1)-r(i,2*M))/(2*delphi);
        rtt(i,2*M+1)=(r(i+1,2*M+1)-2*r(i,2*M+1)+r(i-1,2*M+1))/(deltheta^2);
        rpp(i,2*M+1)=(r(i,2*M+1)-2*r(i,2*M+1)+r(i,2*M))/(delphi^2);
        rtp(i,2*M+1)=(r(i+1,2*M+1)-r(i+1,2*M)-r(i-1,2*M+1)+r(i-1,2*M))/(4*deltheta*delphi);
    end
%corner points
        rt(1,1)=(r(2,1)-r(1,1))/(2*deltheta);
        rp(1,1)=(r(1,2)-r(1,1))/(2*delphi);
        rtt(1,1)=(r(2,1)-2*r(1,1)+r(1,1))/(deltheta^2);
        rpp(1,1)=(r(1,2)-2*r(1,1)+r(1,1))/(delphi^2);
        rtp(1,1)=(r(2,2)-r(2,1)-r(1,2)+r(1,1))/(4*deltheta*delphi);

        rt(1,2*M+1)=(r(2,2*M+1)-r(1,2*M+1))/(2*deltheta);
        rp(1,2*M+1)=(r(1,2*M+1)-r(1,2*M))/(2*delphi);
        rtt(1,2*M+1)=(r(2,2*M+1)-2*r(1,2*M+1)+r(1,2*M+1))/(deltheta^2);
        rpp(1,2*M+1)=(r(1,2*M+1)-2*r(1,2*M+1)+r(1,2*M))/(delphi^2);
        rtp(1,2*M+1)=(r(2,2*M+1)-r(2,2*M)-r(1,2*M+1)+r(1,2*M))/(4*deltheta*delphi);

        rt(N+1,1)=(r(N+1,1)-r(N,1))/(2*deltheta);
        rp(N+1,1)=(r(N+1,2)-r(N+1,1))/(2*delphi);
        rtt(N+1,1)=(r(N+1,1)-2*r(N+1,1)+r(N,1))/(deltheta^2);
        rpp(N+1,1)=(r(N+1,2)-2*r(N+1,1)+r(N+1,1))/(delphi^2);
        rtp(N+1,1)=(r(N+1,2)-r(N+1,1)-r(N,2)+r(N,1))/(4*deltheta*delphi);

        rt(N+1,2*M+1)=(r(N+1,2*M+1)-r(N,2*M+1))/(2*deltheta);
        rp(N+1,2*M+1)=(r(N+1,2*M+1)-r(N+1,2*M))/(2*delphi);
        rtt(N+1,2*M+1)=(r(N+1,2*M+1)-2*r(N+1,2*M+1)+r(N,2*M+1))/(deltheta^2);
        rpp(N+1,2*M+1)=(r(N+1,2*M+1)-2*r(N+1,2*M+1)+r(N+1,2*M))/(delphi^2);
        rtp(N+1,2*M+1)=(r(N+1,2*M+1)-r(N+1,2*M)-r(N,2*M+1)+r(N,2*M))/(4*deltheta*delphi);
        
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
%% Step 4: Calculation of power transmission coefficient
        Ei(i,j,:)=Ef(i,j,:)*(exp(-1i*kd*r(i,j))/r(i,j));
        Rn(i,j)=abs((nd*cos(thetai(i,j))-cos(thetat(i,j)))/(nd*cos(thetai(i,j))+cos(thetat(i,j))))^2;
        Rp(i,j)=abs((nd*cos(thetat(i,j))-cos(thetai(i,j)))/(nd*cos(thetat(i,j))+cos(thetai(i,j))))^2;
        Ein(i,j)=dot(Ei(i,j,:),(cross(n(i,j,:),ki(i,j,:))/norm(cross(n(i,j,:),ki(i,j,:)))));
        Eip(i,j)=norm(E(i,j,:)-(Ein(i,j)*(cross(n(i,j,:),ki(i,j,:))/norm(cross(n(i,j,:),ki(i,j,:))))));
        Ein(i,j)=norm(Ein(i,j));
        Tk(i,j)=1-(((Rp(i,j)^2*Eip(i,j)^2)+(Rn(i,j)^2*Ein(i,j)^2))/(Ein(i,j)^2+Eip(i,j)^2));
%% Step 5: Calculation of \alpha_k and \beta_k
        alphak(i,j)=acos(kt(i,j,3));
        betak(i,j)=atan(kt(i,j,2)/kt(i,j,1));
%% Step 6: Calculation of partial derivatives of rk
        if((i~=1)&&(i~=N+1)&&(j~=1)&&(j~=2*M+1))
            rt(i,j)=(r(i+1,j)-r(i-1,j))/(2*deltheta);
            rp(i,j)=(r(i,j+1)-r(i,j-1))/(2*delphi);
            rtt(i,j)=(r(i+1,j)-2*r(i,j)+r(i-1,j))/(deltheta^2);
            rpp(i,j)=(r(i,j+1)-2*r(i,j)+r(i,j-1))/(delphi^2);
            rtp(i,j)=(r(i+1,j+1)-r(i+1,j-1)-r(i-1,j+1)+r(i-1,j-1))/(4*deltheta*delphi);
        end
    end
end
%%
end