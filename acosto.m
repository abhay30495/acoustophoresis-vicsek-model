%Prediction model acoustic forces in O/W emulsion
%%
clear all
close all
clc

%% Initialization
Tmax=2000;          %duration
dt=0.1;             %time between two frames [s]
L=1320*10^-6;       %channel length [m]
w0=600*10^-6;       %channel width [m]
N=350;              %amount of oil droplets [-]
R=3*10^-6;          %average diameter of oil droplets [m]
r=4*R;              %interaction radius [m]
noise=0.05;         %noise

%% Boundary condition
PERIODIC=1;                   %1: periodic boundary condition, 0: unlimited

%% Initial condition
x=L*rand(1,N);                %random initial x-position
y=w0*rand(1,N);               %random initial y-position
w=w0*rand(1,N);
vel=ones(1,350)*0.5*10^-6;    %initial velocity [m/s]
theta=2*pi*(rand(1,N)-0.5);   %random distribution of particle directions
lambda_h=w0*2;                %half-wavelength
c=1497;                       %speed of sound in water [m/s]
f_p=c/(2*w0);                 %periodic frequency
Ek_e = zeros(1,N);
t = 1000*ones(1,N);

 
xP = [0.15*lambda_h    0.4*lambda_h     0.65*lambda_h    0.9*lambda_h]; %x-position antinodes
yP = [w0    0      w0   0];                                             %y-position antinodes


tmp_x = zeros(1,N);                   
tmp_y = zeros(1,N);
ave_theta=2*pi*(rand(1,N)-0.5);   
tmp = zeros(1,N);
tsep = 0.05*ones(1,N);
%% Video

vidObj = VideoWriter('acousto.avi');
open(vidObj);
for time=1:Tmax
    
    % Periodic boundary %
    if PERIODIC==1
        tmp_x(x<r) = L + x(x<r);
        tmp_x(x>L-r) = x(x>L-r)-L;
        tmp_x(r<=x & x<=L-r) = x(r<=x & x<=L-r);
        
        tmp_y(y<r) = w0 + y(y<r);
        tmp_y(y>w0-r) = y(y>w0-r)-w0;
        tmp_y(r<=y & y<=w0-r) = y(r<=y & y<=w0-r);
        
      %  tmp_D = pdist([tmp_x' tmp_y'],'euclidean');
      %  D = min([D; tmp_D]);
    end
    
    
%     tapp = 0.05;
%     dt = tsep + tapp;

%     
%     M = squareform(D); % Matrix representation for the disaatan22ce between particles
%     [l1,l2]=find(0<M & M<r);
        
%     for i = 1:N
%         list = l1(l2==i);
%         if ~isempty(list)
%             ave_theta(i) = atan222(mean(sin(theta(list))),mean(cos(theta(list))));
%         else
%             ave_theta(i) = theta(i);
%         end
%     end
    
    for i = 1:N
       list1 = sqrt((tmp_x(i)-xP(1))^2+(tmp_y(i)-yP(1))^2);
       list2 = sqrt((tmp_x(i)-xP(2))^2+(tmp_y(i)-yP(2))^2);
       list3 = sqrt((tmp_x(i)-xP(3))^2+(tmp_y(i)-yP(3))^2);
       list4 = sqrt((tmp_x(i)-xP(4))^2+(tmp_y(i)-yP(4))^2);
       
        a = min([list1 list2 list3 list4]);
        w(i) = a;
        
     
       if list1 == min([list1 list2 list3 list4])
        ave_theta(i) = (atan((tmp_y(i)-yP(1))/(tmp_x(i)-xP(1))));
    elseif list2 == min([list1 list2 list3 list4])
        ave_theta(i) = (atan((tmp_y(i)-yP(2))/(tmp_x(i)-xP(2))));
    elseif list3 == min([list1 list2 list3 list4])
        ave_theta(i) = (atan((tmp_y(i)-yP(3))/(tmp_x(i)-xP(3))));
       else 
        ave_theta(i) = (atan((tmp_y(i)-yP(4))/(tmp_x(i)-xP(4))));
       end
   
    end
      %%
   %  if position droplet is ~equal to antinode -> velocity = 0
%      for i=1:N
%          if round(min([list1 list2 list3 list4]),6)==0
%              vel=0;
%          end
%      end 

       
    
    %% Update
     c2 = 1274;                                                  %speed of sound in SO [m/s] http://www.precisionflow.co.uk/speed_of__sound.htm
     eta = 0.89*10^-3;                                           %viscosity water [Pa*s]=[kg*m-1*s-1]
     rho_p = 770;                                                %density of SO AR20 [kg/m3]
     rho_m = 997;                                                %density water [kg/m3]
     beta_p = 1/(rho_p*c2^2);                                    %compressibility oil [m2s2kg-1]
     beta_m = 1/(rho_m*c^2);                                     %compressibility water [m2s2kg-1]
     phi = (5*rho_p-2*rho_m)/(2*rho_p+rho_m)-(beta_p)/(beta_m);  %acoustic constrast factor SO AR20 [-]
     k = 2*pi/(1200*10^-6);                                      %wavenumber [m-1]
     p0 = 25;
     p = p0.*cos(k*w);
     Eac = 0.25*beta_m*(p.*p);
%    Eac = 0.5*10^-6*3*eta./(2*phi*k*R^2.*sin(2*k.*w));          %acoustic energy density [J/m3]=[kgs-2m-1]
     m = rho_p*(4/3)*pi*R^3;                                     %mass droplet [kg]
     
     
     k_v = 0.8*10^-6;                                            %kinematic viscosity of water m^2/sec
     mu = 8.90*10^-4;                                            %dynamic viscosity of water [Pa.s]
     gamma = 1.4;                                                %specific heat ratio for air
     stk = 1;                                                    %Stokes number assumed 1 for now 
     dels = sqrt(k_v/(pi*f_p));                                  %Thickness of viscous boundary layer
     wt = 2*pi*stk;
     b = 2*dels/R;
     mu_p = sqrt(1/(1+(wt)^2));                                  %Entrainment coefficient for the particle
     mu_g = sqrt(1-mu_p);                                        %flow around coefficient
     tht = pi/2;
     omega = 2*pi*f_p;
     d = 16*R;
     
     Fac = (-4*pi/3).*(R^3*k*Eac*phi.*sin(2*k.*w));                                          %acoustic radiation force [N]
    %Fdrag = 6*pi*eta*R*vel;                                                                 %stokes drag force [N]
     Fassym = -pi/24*Eac*k*(R^3)*mu_p*(4.5*((b^2)+b)*mu_g + (3+4.5*b)*mu_p).*sin(2*k.*w);    %assymeteric force [N] 
     Fvisc = 3*pi*mu/(rho_m*c)*(gamma-3)*R*(mu_g^2)*Eac.*sin(2*k.*w);                        %Viscous force [N]            
     Fsec = (pi*R^6/16)*((3.*cos(tht)-1)*(rho_p - rho_m)^2*(0.5*10^-6)^2/(6*rho_m*d^4)-(beta_p-beta_m)^2*rho_m*omega^2*(p.*p)/(9*d^2));  
     
     
     
     F_net = Fac + Fassym - Fvisc + Fsec ;                                                   %resulaatan22t force
    
     dv=(dt/m)*(F_net);                                                                      %acceleration
     
     vel=vel+dv; 
                                                                                %velocity = old velocity + delta velocity
    
     aa(time,:) = vel;
      

    %% Theta and position update
      theta = ave_theta;                   %direction of movement (+ noise*(rand(1,N) - 0.5);)

      x = x + vel.*cos(2*pi.*theta)*dt;    %new x-position
      y = y + vel.*sin(2*pi.*theta)*dt;    %new y-position

     if PERIODIC==1
         x(x<0) = L + x(x<0);
         x(L<x) = x(L<x) - L;
         y(y<0) = w0 + y(y<0);
         y(w0<y) = y(w0<y) - w0;
     end
     
  %   Ek_e = Ek_e + 0.5*m*(vel.*vel)*N; 
%% Figure
    figure(1)
    plot(x,y,'.','MarkerSize',7)
    %hold on
    % quiver(x,y,x.*cos(theta),y.*sin(theta),'LineWidth',1,...
    %'Color',[0 0 0],'AutoScaleFactor',0.3);
    xlim([0 L]); 
    ylim([0 w0]);
    axis auto
   % pause(0.01);    
 %      hold on
 %      plot(xP,yP,'ko')
     currFrame = getframe;
     writeVideo(vidObj,currFrame);
end
close(vidObj);
%% PLotting K.E. of the system 
% figure(2)
% %Ek_e = Ek_e + 0.5*m*(aa.*aa)*N; 
% for i = 1:Tmax
% bb(i) = sum(0.5*m*(aa(:,i).*aa(:,i))); 
% end
% 
% plot(1:Tmax,bb)