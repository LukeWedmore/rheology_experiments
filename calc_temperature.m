function [T] = calc_temperature(z,Q0,H_uc,K_uc,T0,H_lc,K_lc,uc_thickness,lc_thickness)
% z = depth in m
% Q0 = surface heat flow value (W m^-2)
% H_uc = surface heat production in the upper crust (Wm^-3)
% K_uc = thermal conductivity in the upper crust (Wm^-1 K^-1)
% T0 = initial temperature
% H_lc = surface heat production in the lower crust (Wm^-3)
% K_lc = thermal conductivity in the lower crust (Wm^-1 K^-1)
% uc_thickness = thickness of the upper crust (in m)
% lc_thickness = thickness of the lower crust (in m)

T = zeros(size(z));

%% Calculate Upper crust temperature

T(z<=uc_thickness) = T0 + (Q0.*z(z<=uc_thickness)/K_uc) - (H_uc.*z(z<=uc_thickness).^2)/(2*K_uc);

%update values for the lower crust
z_up = z(z>uc_thickness) - max(z(z<=uc_thickness));
T0 = max(T);
Q0 = Q0 - H_uc*uc_thickness;

T(z>uc_thickness) = T0 + (Q0.*z_up/K_lc) - (H_lc.*z_up.^2)/(2*K_lc);

% 