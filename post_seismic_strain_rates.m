% Postseismic deformation in a power law shear zone from Montesi (2004)
clear all
close all

% Vs(t) = V0[1+(1-1/n)t/tau]^(-1/(1-1/n))
%let g=(1-1/n)
% Vs(t) = V0[1+g*t/tau]^(-1/g)
% differentiate to fine the strain rate
% epsilon_dot=V0[-1/g - t/tau]^(-1/g)-1

%scenario - initial velocity of 1 mm/yr with a maxwell decay time of 50
%years 

V0=0.0005; %1 mm/yr
V0=0.5; 
n=3;
tau=0.09; %in years
t=0:0.1:10000;

vs = V0*(1+(1-1/n).*t./tau).^(-1/(1-1/n));
n=5;
vsn5=V0*(1+(1-1/n).*t./tau).^(-1/(1-1/n));
n=1.000000001;
vsn1=V0*(1+(1-1/n).*t./tau).^(-1/(1-1/n));

%strain rate - divide by a distance varying from 100 m to 1m
epsilon_dot_100=vs/100;
epsilon_dot_10 = vs/10;
epsilon_dot_1 = vs/1;

epsilon_dot_n5_100=vsn5/100;
epsilon_dot_n5_10 = vsn5/10;
epsilon_dot_n5_1 = vsn5/1;

epsilon_dot_n1_100=vsn1/100;
epsilon_dot_n1_10 = vsn1/10;
epsilon_dot_n1_1 = vsn1/1;

%%
z = [0:100:100000];
%Calculate temperature
%Q0 = 63e-3; %mW m^-2 from Nyblade 1990
Q0 = 68e-3; %the mean value from the southwestern and western rift branches
H=2.1e-6;
K=2.8;
T0=273.15; %zero degrees celsius in kelvin

T_rift=zeros(1,401);
T_rift(z<=9000) = T0 + ((Q0.*z(z<=9000))/K) - ((H.*z(z<=9000).^2)/(2*K));

%update values for lower crust
z_up = z(z>9000)-max(z(z<=9000));
T0=max(T_rift);
Q0=Q0-H*9000;
H=0.4e-6;
K=2.0;

T_rift(z>9000) = T0 + ((Q0.*z_up)./K) - ((H.*(z_up.^2))/(2*K));

%%Calculate differential stress

%Frictional regime
mu = 0.60;
%alpha = (sqrt(1+mu^2))^-2; %original equation given to me by Ake who
%claimed it's adapted from Sibson 1974 (Nature) - but not seen replicated
%elsewhere and weirdly results in an increase in fault strength with a
%reduction in friction
alpha = 2*mu/((1+mu^2)^0.5 + mu); %From Brace and Kohlstead, 1980 (JGR) and replicated elsewhere
rho_uc = 2650;
rho_lc = 2800;
g = 9.8;
lambda_v = 0.4;

frictional_strength = nan(size(z));

%upper crust
frictional_strength(z<=16000) = alpha * rho_uc * g .* z(z<=16000) * (1-lambda_v);

%lower crust
fs_update = max(frictional_strength);
frictional_strength(z>16000) = fs_update + (alpha * rho_lc * g .* (z(z>16000)-16000) * (1-lambda_v));

%viscous regime

%--------------------------------
%Dislocation Creep
%--------------------------------
strain_rate_background=1e-16;
strain_rate_1_week = epsilon_dot_10(t==0.1)/31557600;
strain_rate_1_yr = epsilon_dot_10(t==1)/31557600;
strain_rate_10_yrs = epsilon_dot_10(t==10)/31557600;
strain_rate_100_yrs = epsilon_dot_10(t==100)/31557600;
strain_rate_1000_yrs = epsilon_dot_10(t==1000)/31557600;

epsilon_dot_mantle_lithosphere = 1e-15;

%universal gas constant
R =  8.3145;

crustal_thickness=41000; %median crustal thickness from receiver function measurements of 3 transects


%-------------------
%Quartz upper crust
%from Rutter and Brodie 2004
n_uc = 3;
Q_uc = 242e3;
A_uc = 10^-4.9 * 1e6^-n_uc;
     
%-------------------
%FELSIC LOWER CRUST
%Values for synthetic anorthite (feldspar) from Rybacki and Dresen, 2000;
%JGR, Table 2.
%dry
n_flc = 3.0; %±0.4
Q_flc = 648e3; %± 20; kJ mol^-1; convert to J.
A_flc = 10^12.7 * 1e6^-n_flc; %LogA = 12.7 ± 0.8 MPa^-n µm^m s^-1; convert to Pa s-1 (div by 1e6^-n, m=0).

%wet
n_flc_w = 3.0; %±0.2
Q_flc_w = 356e3; %± 9; kJ mol^-1; convert to J.
A_flc_w = 10^2.6 * 1e6^-n_flc_w; %LogA = 2.6 ± 0.3 MPa^-n µm^m s^-1; convert to m (div by 1e6).
%-------------------
%DRY OLIVINE MANTLE LITHOSPHERE
%from Hirth and Kohlstedt, 2003 - taken from Burgmann and Dresen 2009
n_doml_disc = 3.5;
Q_doml_disc = 530e3; %±4; kJ mol^-1; convert to J.
A_doml_disc = 10^5 * 1e6^-n_doml_disc; %LogA in MPa^-n µm^m s^-1; convert to m (div by 1e6).

%--------------------------------
%Diffusion Creep
%--------------------------------
%for use in the mantle lithosphere below the rifts.

%-------------------
%DRY OLIVINE MANTLE LITHOSPHERE
%from Faul and Jackson, 2007, JGR. see paragraph 32.
Q_doml_difc = 484e3; %±30
n_doml_difc = 1.37; %±0.06
m_doml_difc = 3;
d_doml_difc = 50/1e6; %grain size range given as 2.7-5.6µm
A_doml_difc = (10^10.3 * 1e6^-n_doml_difc) / 1e6^m_doml_difc; % ± 4 MPa^-n µm^m s^-1; convert to m (div by 1e6).

%--------------------------------
%calculate viscous strength curves
%--------------------------------

viscous_stress_background = NaN(size(z));
viscous_stress_1_week = NaN(size(z));
viscous_stress_1_yr = NaN(size(z));
viscous_stress_10_yrs = NaN(size(z));
viscous_stress_100_yrs = NaN(size(z));
viscous_stress_1000_yrs = NaN(size(z));

%upper crust
%dislocation creep of Quartz
viscous_stress_background(z<=9000) = (strain_rate_background./(A_uc*exp(-Q_uc./(R.*T_rift(z<=9000))))).^(1/n_uc);
viscous_stress_1_week(z<=9000) = (strain_rate_1_week./(A_uc*exp(-Q_uc./(R.*T_rift(z<=9000))))).^(1/n_uc);
viscous_stress_1_yr(z<=9000) = (strain_rate_1_yr./(A_uc*exp(-Q_uc./(R.*T_rift(z<=9000))))).^(1/n_uc);
viscous_stress_10_yrs(z<=9000) = (strain_rate_10_yrs./(A_uc*exp(-Q_uc./(R.*T_rift(z<=9000))))).^(1/n_uc);
viscous_stress_100_yrs(z<=9000) = (strain_rate_100_yrs./(A_uc*exp(-Q_uc./(R.*T_rift(z<=9000))))).^(1/n_uc);
viscous_stress_1000_yrs(z<=9000) = (strain_rate_1000_yrs./(A_uc*exp(-Q_uc./(R.*T_rift(z<=9000))))).^(1/n_uc);

%mid- to lower-curst
%dislocation creep of feldspar (anorthite)
viscous_stress_background(z>9000 & z<=crustal_thickness) = (strain_rate_background./(A_flc*exp(-Q_flc./(R.*T_rift(z>9000 & z<=crustal_thickness))))).^(1/n_flc);
viscous_stress_1_week(z>9000 & z<=crustal_thickness) = (strain_rate_1_week./(A_flc*exp(-Q_flc./(R.*T_rift(z>9000 & z<=crustal_thickness))))).^(1/n_flc);
viscous_stress_1_yr(z>9000 & z<=crustal_thickness) = (strain_rate_1_yr./(A_flc*exp(-Q_flc./(R.*T_rift(z>9000 & z<=crustal_thickness))))).^(1/n_flc);
viscous_stress_10_yrs(z>9000 & z<=crustal_thickness) = (strain_rate_10_yrs./(A_flc*exp(-Q_flc./(R.*T_rift(z>9000 & z<=crustal_thickness))))).^(1/n_flc);
viscous_stress_100_yrs(z>9000 & z<=crustal_thickness) = (strain_rate_100_yrs./(A_flc*exp(-Q_flc./(R.*T_rift(z>9000 & z<=crustal_thickness))))).^(1/n_flc);
viscous_stress_1000_yrs(z>9000 & z<=crustal_thickness) = (strain_rate_1000_yrs./(A_flc*exp(-Q_flc./(R.*T_rift(z>9000 & z<=crustal_thickness))))).^(1/n_flc);

%mantle lithosphere
%diffusion creep of dry olivine (rift)
viscous_stress_background(z>crustal_thickness) = (epsilon_dot_mantle_lithosphere./(A_doml_difc*(d_doml_difc^-m_doml_difc)*exp(-Q_doml_difc./(R.*T_rift(z>crustal_thickness))))).^(1/n_doml_difc);
viscous_stress_1_week(z>crustal_thickness) = (epsilon_dot_mantle_lithosphere./(A_doml_difc*(d_doml_difc^-m_doml_difc)*exp(-Q_doml_difc./(R.*T_rift(z>crustal_thickness))))).^(1/n_doml_difc);
viscous_stress_1_yr(z>crustal_thickness) = (epsilon_dot_mantle_lithosphere./(A_doml_difc*(d_doml_difc^-m_doml_difc)*exp(-Q_doml_difc./(R.*T_rift(z>crustal_thickness))))).^(1/n_doml_difc);
viscous_stress_10_yrs(z>crustal_thickness) = (epsilon_dot_mantle_lithosphere./(A_doml_difc*(d_doml_difc^-m_doml_difc)*exp(-Q_doml_difc./(R.*T_rift(z>crustal_thickness))))).^(1/n_doml_difc);
viscous_stress_100_yrs(z>crustal_thickness) = (epsilon_dot_mantle_lithosphere./(A_doml_difc*(d_doml_difc^-m_doml_difc)*exp(-Q_doml_difc./(R.*T_rift(z>crustal_thickness))))).^(1/n_doml_difc);
viscous_stress_1000_yrs(z>crustal_thickness) = (epsilon_dot_mantle_lithosphere./(A_doml_difc*(d_doml_difc^-m_doml_difc)*exp(-Q_doml_difc./(R.*T_rift(z>crustal_thickness))))).^(1/n_doml_difc);

%combine different curves and convert to MPa
stress_background = min(viscous_stress_background./1e6,frictional_strength./1e6);
stress_1_week = min(viscous_stress_1_week./1e6,frictional_strength./1e6);
stress_1_yr = min(viscous_stress_1_yr./1e6,frictional_strength./1e6);
stress_10_yrs = min(viscous_stress_10_yrs./1e6,frictional_strength./1e6);
stress_100_yrs = min(viscous_stress_100_yrs./1e6,frictional_strength./1e6);
stress_1000_yrs = min(viscous_stress_1000_yrs./1e6,frictional_strength./1e6);

%--------------------------------------------------------------------------
%%Plot
figure('Position',[200,200,700,500])
subplot(2,2,1)
loglog(t,vsn1.*1e3,'black')
hold on
loglog(t,vs.*1e3,'blue')
loglog(t,vsn5.*1e3,'red')

% plot data from Mozambique earthquake (data in Ingleby and Wright, 2017).
time = [0.25,1.85,3.74];
velocity = [104.29,17.36,2.83]; %in mm/yr
loglog(time,velocity,'+')
%loglog(t,vt100.*1e3)

txt = '(a)';
text(0,0,txt,'Position',[0.05 0.95 0],'Units','normalized','FontSize',12,'FontName','helvetica')


xlabel('time (years)')
ylabel('velocity (mm/yr)')
title('Postseismic deformation in shear zones')
txt = {'V_0 = 500 mm/yr','\tau = 0.09 yrs^{-1}'};
text(300,0.08,txt)
txt = {'$$ V_s(t) = V_{0} \left[ 1 + \left( 1-\frac{1}{n} \right)\frac{1}{\tau} \right]^{ \left( \frac{-1}{1-1/n} \right)} $$'};
text(2,100,txt,'Interpreter','latex')

legend('n=1','n=3','n=5','Mozambique','Location','east')
ylim([0.01 600])
xlim([0.1 10000])

subplot(2,2,3)
loglog(t,epsilon_dot_100./31557600,'-.blue'); %divide to get strain rate per second
hold on
loglog(t,epsilon_dot_10./31557600,'-blue');
loglog(t,epsilon_dot_1./31557600,'--blue');

loglog(t,epsilon_dot_n5_100./31557600,'-.red');
loglog(t,epsilon_dot_n5_10./31557600,'-red');
loglog(t,epsilon_dot_n5_1./31557600,'--red');

loglog(t,epsilon_dot_n1_100./31557600,'-.black');
loglog(t,epsilon_dot_n1_10./31557600,'-black');
loglog(t,epsilon_dot_n1_1./31557600,'--black');

xlabel('time (years)')
ylabel('strain rate (s^{-1})')
legend('100 m wide shear zone','10 m wide shear zone','1 m wide shear zone','Location','northeast')
title('Postseismic strain rate with time')
%text(0.1,3e-11,'1 m wide shear zone')
%text(0.1,3e-12,'10 m wide shear zone')
%text(0.1,3e-13,'100 m wide shear zone')
txt = '(b)';
text(0,0,txt,'Position',[0.05 0.95 0],'Units','normalized','FontSize',12,'FontName','helvetica')


xlim([0.1 10000])
ylim([1e-16 1e-7])

subplot(2,2,[2,4])
plot(stress_background,z./1e3,'k','LineWidth',2)
hold on
plot(stress_1_yr,z./1e3,'LineWidth',1)
plot(stress_10_yrs,z./1e3,'LineWidth',1)
plot(stress_100_yrs,z./1e3,'LineWidth',1)
plot(stress_1000_yrs,z./1e3,'LineWidth',1)
set(gca,'YDir','reverse')
xlabel('\sigma_{1} - \sigma_{3} (MPa)')
ylabel('Depth (km)')
xlim([0 500])
ylim([0 50])
title('Crustal strength changes during postseismic deformation')
legend('background strain rate','1 yr','10 yrs','100 yrs','1000 yrs')

txt = '(c)';
text(0,0,txt,'Position',[0.05 0.98 0],'Units','normalized','FontSize',12,'FontName','helvetica')

time_ps=[1,10,100,1000];

%calculate BDT
[~,y1yr] = max(stress_1_yr(z<=crustal_thickness));
[~,y10yr] = max(stress_10_yrs(z<=crustal_thickness));
[~,y100yr] = max(stress_100_yrs(z<=crustal_thickness));
[~,y1000yr] = max(stress_1000_yrs(z<=crustal_thickness));
bdt(1) = z(y1yr);
bdt(2) = z(y10yr);
bdt(3) = z(y100yr);
bdt(4) = z(y1000yr);

% figure
% semilogx(time_ps,bdt./1e3,'+-')
% xlabel('time (years)')
% ylabel('brittle ductile transition (km)')
% 
% ylim([30 40])
