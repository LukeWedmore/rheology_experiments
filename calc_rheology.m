function [differential_stress] = calc_rheology(mu,rho_uc,rho_lc,pore_fluid_factor,z,uc_thickness,crustal_thickness,crustal_strain_rate,mantle_strain_rate,n_uc,Q_uc,log_A_uc,n_lc,Q_lc,log_A_lc,Q_ml,n_ml,m_ml,d_ml,log_A_ml,T)
% mu = coefficient of friction
% rho_uc = density of material in upper crust
% rho_lc = density of material in lower crust
% pore fluid factor = 
% z = depth in m
% uc_thickness = upper crustal thickness (in m)
% crustal_thickness = crustal thickness in m
% crustal_strain_rate = in s^-1

%% Constants
g = 9.8; %accelleration due to gravity
R = 8.3145; %universal gas constant

%% SETUP PARAMETERS
frictional_strength = NaN(size(z));
viscous_strength = NaN(size(z));

%% Upper Crust
alpha = 2 * mu/((1 + mu^2)^0.5 + mu); %From Brace and Kohlstead, 1980 (JGR) and replicated elsewhere

frictional_strength(z<=uc_thickness) = alpha * rho_uc * g .* z(z<=uc_thickness) * (1 - pore_fluid_factor);

%% Lower Crust
fs_update = max(frictional_strength);
frictional_strength(z>uc_thickness) = fs_update + (alpha * rho_lc * g .* (z(z>uc_thickness) - uc_thickness) * (1 - pore_fluid_factor));

%% Viscous Regime (Tempterature dependent)

%crust - dislocation creep
A_uc = (10^log_A_uc) * 1e6^-n_uc; %convert to Pa (*1e6) - original unit MPa^-n µm^m s^-1 but m=0
A_lc = (10^log_A_lc) * 1e6^-n_lc;

%mantle lithosphere - diffusion creep
A_ml = (10^log_A_ml * 1e6^-n_ml) / 1e6^m_ml;

%upper crust - dislocation creep
viscous_strength(z<=uc_thickness) = (crustal_strain_rate./(A_uc * exp(-Q_uc./(R.*T(z<=uc_thickness))))).^(1/n_uc);

%lower crust - dislocation creep
viscous_strength(z>uc_thickness & z<=crustal_thickness) = (crustal_strain_rate./(A_lc * exp(-Q_lc./(R.*T(z>uc_thickness & z<=crustal_thickness))))).^(1/n_lc);

%mantle lithosphere - diffusion creep
viscous_strength(z>crustal_thickness) = (mantle_strain_rate./(A_ml*(d_ml^-m_ml) * exp(-Q_ml./(R.*T(z>crustal_thickness))))).^(1/n_ml);

%% Differential Stress
%calculate minimum by combining frictional and viscous curves and convert
%to MPa
differential_stress=min(frictional_strength,viscous_strength)./1e6;

% plot(frictional_strength./1e6,z./1e3)
% hold on
% plot(viscous_strength./1e6,z./1e3)
% set(gca,'YDir','reverse')
% xlabel('\sigma_{1} - \sigma_{3} (MPa)')
% ylabel('Depth (km)')
% xlim([0 550])
% ylim([0 50])


