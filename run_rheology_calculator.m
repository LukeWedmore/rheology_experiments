% run scripts
close all
clear all

z = [0:100:100000];

%temperature parameters
Q0=68e-3;
T0 = 273.15;

H_uc = 2.1e-6;
K_uc = 2.8;

H_lc = 0.4e-6;
K_lc = 2.0;

uc_thickness = 9000;
crustal_thickness = 41000;

%rheological parameters
mu = 0.4;
rho_uc = 2650;
rho_lc = 2800;
pore_fluid_factor = 0.4;

%strain rate parameters
crustal_strain_rate = 1e-15;
mantle_strain_rate = 1e-15;

%dislocation creep parameters
%for Quartz Upper crust from Rutter and Brodio 2004
n_uc = 3;
Q_uc = 242e3;
log_A_uc = -4.9;

%for dry Felsic lower crust
%Values for synthetic anorthite (feldspar) from Rybacki and Dresen, 2000;
%JGR, Table 2.
n_lc = 3.0;
Q_lc = 648e3;
log_A_lc = 12.7;

%Diffusion Creep Parameters
%DRY OLIVINE MANTLE LITHOSPHERE
%from Faul and Jackson, 2007, JGR. see paragraph 32.
Q_ml = 484e3;
n_ml = 1.37;
m_ml = 3;
d_ml = 50/1e6; %grain size range given as 2.7-5.6µm
log_A_ml = 10.3;




T = calc_temperature(z,Q0,H_uc,K_uc,T0,H_lc,K_lc,uc_thickness,crustal_thickness);
diff_stress = calc_rheology(mu,rho_uc,rho_lc,pore_fluid_factor,z,uc_thickness,crustal_thickness,crustal_strain_rate,mantle_strain_rate,n_uc,Q_uc,log_A_uc,n_lc,Q_lc,log_A_lc,Q_ml,n_ml,m_ml,d_ml,log_A_ml,T);

figure
subplot(1,2,1)
plot(T-273.15,z./1e3)
hold on
set(gca,'YDir','reverse')
ylim([0 50])

subplot(1,2,2)
plot(diff_stress,z./1e3)
hold on
%plot(fric_strength./1e6,z./1e3)
set(gca,'YDir','reverse')
ylim([0 50])

%vary strain rate and pore fliud factor and analyse depth of Brittle Ductil
%Transition

figure('Position',[200,200,600,1200])
subplot(3,2,1)
plot(T-273.15,z./1e3,'LineWidth',2) 
hold on
plot([0 1200],[9 9],'k--')
plot([0 1200],[41 41],'k--')
hold on
set(gca,'YDir','reverse')
ylabel('Depth (km)')
xlabel('Temperature ({\circ}C)')
ylim([0 50])
title('Temperature profile')

txt='upper crust: 0-9 km';
text(225,7,txt)
txt='mid-lower crust: 9-41 km';
text(100,39,txt)
txt='lithospheric mantle: 41+ km';
text(100,48,txt)

text(0,0,'(a)','Position',[0.93 0.96 0],'Units','normalized','FontSize',12,'FontName','helvetica')


subplot(3,2,2)
clear T
T_2 = calc_temperature(z,Q0,H_uc,K_uc,T0,H_lc,K_lc,uc_thickness,25000);
diff_stress_25 = calc_rheology(mu,rho_uc,rho_lc,pore_fluid_factor,z,uc_thickness,25000,crustal_strain_rate,mantle_strain_rate,n_uc,Q_uc,log_A_uc,n_lc,Q_lc,log_A_lc,Q_ml,n_ml,m_ml,d_ml,log_A_ml,T_2);
plot(diff_stress_25,z./1e3,'--','LineWidth',2)
hold on
T = calc_temperature(z,Q0,H_uc,K_uc,T0,H_lc,K_lc,uc_thickness,crustal_thickness);
diff_stress_41 = calc_rheology(mu,rho_uc,rho_lc,pore_fluid_factor,z,uc_thickness,crustal_thickness,crustal_strain_rate,mantle_strain_rate,n_uc,Q_uc,log_A_uc,n_lc,Q_lc,log_A_lc,Q_ml,n_ml,m_ml,d_ml,log_A_ml,T);
plot(diff_stress_41,z./1e3,'r','LineWidth',2)
set(gca,'YDir','reverse')
xlabel('\sigma_{1} - \sigma_{3} (MPa)')
ylabel('Depth (km)')
xlim([0 550])
ylim([0 50])
title('Crustal thickness and lithospheric strength')
l1 = sprintf('25 km: strength = %3.1f GPa', (trapz(diff_stress_25)/1e3));
l2 = sprintf('41 km: strength = %3.1f GPa', (trapz(diff_stress_41)/1e3));

legend(l1,l2);
legend('boxoff')
text(0,0,'(c)','Position',[0.93 0.96 0],'Units','normalized','FontSize',12,'FontName','helvetica')



subplot(3,2,3)
T = calc_temperature(z,Q0,H_uc,K_uc,T0,H_lc,K_lc,uc_thickness,crustal_thickness);
hold on
set(gca,'YDir','reverse')
xlabel('\sigma_{1} - \sigma_{3} (MPa)')
ylabel('Depth (km)')
xlim([0 550])
ylim([0 50])
title('Change in crustal strain rate')
text(0,0,'(c)','Position',[0.93 0.96 0],'Units','normalized','FontSize',12,'FontName','helvetica')


%vary strain rate
count = 1;
for i=[logspace(-15,-10,6)]
    crustal_strain_rate = i;
    diff_stress = calc_rheology(mu,rho_uc,rho_lc,pore_fluid_factor,z,uc_thickness,crustal_thickness,crustal_strain_rate,mantle_strain_rate,n_uc,Q_uc,log_A_uc,n_lc,Q_lc,log_A_lc,Q_ml,n_ml,m_ml,d_ml,log_A_ml,T);
    plot(diff_stress,z./1e3);
    
    %calculate the BDT
    [~,y] = max(diff_stress(z<=crustal_thickness));
    bdt_strain_rate(count) = z(y);
    strain_rate(count) = i;
    
    count = count+1;
end

lgd = legend('1e-15 s^{-1}','1e-14 s^{-1}','1e-13 s^{-1}','1e-12 s^{-1}','1e-11 s^{-1}','1e-10 s^{-1}','Location','northeast');
title(lgd,'strain rate')
legend('boxoff')

subplot(3,2,4)
hold on
set(gca,'YDir','reverse')
xlabel('\sigma_{1} - \sigma_{3} (MPa)')
ylabel('Depth (km)')
xlim([0 550])
ylim([0 50])
title('Change in pore fluild factor')
text(0,0,'(d)','Position',[0.93 0.96 0],'Units','normalized','FontSize',12,'FontName','helvetica')


crustal_strain_rate = 1e-15;
count = 1;
for i = [0.4:0.1:0.9]
    pore_fluid_factor = i;
    diff_stress = calc_rheology(mu,rho_uc,rho_lc,pore_fluid_factor,z,uc_thickness,crustal_thickness,crustal_strain_rate,mantle_strain_rate,n_uc,Q_uc,log_A_uc,n_lc,Q_lc,log_A_lc,Q_ml,n_ml,m_ml,d_ml,log_A_ml,T);
    plot(diff_stress,z./1e3);
    
    %calculate the BDT
    [~,y] = max(diff_stress(z<=crustal_thickness));
    bdt_pore_fluid(count) = z(y);
    pore_fluid(count) = i;
    
    count = count+1;    
end

lgd=legend('0.4','0.5','0.6','0.7','0.8','0.9','Location','northeast');
title(lgd,'pore fluid factor')
legend('boxoff')

subplot(3,2,5)
semilogx(strain_rate,bdt_strain_rate./1e3,'-+')
ylim([30 40])
set(gca,'YDir','reverse')
xlabel('crustal strain rate (s^{-1})')
ylabel('brittle ductile transition (km)')
text(0,0,'(e)','Position',[0.93 0.96 0],'Units','normalized','FontSize',12,'FontName','helvetica')


subplot(3,2,6)
plot(pore_fluid,bdt_pore_fluid./1e3,'-+')
ylim([30 40])
set(gca,'YDir','reverse')
xlabel('pore fluild factor')
ylabel('brittle ductile transition (km)')
text(0,0,'(f)','Position',[0.93 0.96 0],'Units','normalized','FontSize',12,'FontName','helvetica')

exportgraphics(gcf,'Figure_4.png','Resolution',300)
