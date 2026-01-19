function Greenp_grd = generate_greenp_grd(Xc, Yc, Zc, ll, ww, DIP, STRIKE, ...
    xP, yP, tpP, looks, nu, fault_type, dat_ph_grd, dx, dy)
% GENERATE_GREENP_GRD Compute Green's function for phase gradients data.
%
%   Inputs:
%      Xc, Yc, Zc, ll, ww, DIP, STRIKE : Fault patch geometry arrays
%      xP, yP, tpP, looks              : UTM coordinates, topography, Looks
%      nu, fault_type                  : Poisson ratio, Fault slip tpye 
%      dat_ph_grd                      : the num of data for each 
%      dx, dy                          : the interval for gradient computation
%
%   Output:
%       Greenp_grd :Green’s matrix

    % --- Split counts for x- and y-gradient subsets ---
    n = numel(dat_ph_grd);
    n_half = n / 2;
    n_x = sum(dat_ph_grd(1:n_half));  % number of gradients along W - E 
    n_y = sum(dat_ph_grd(n_half+1:n)); % number of  gradients along S - N 

    % --- W - E direction finite difference ---
    idx_x = 1:n_x;
    Gp_xp = generate_green_p(Xc, Yc, Zc, ll, ww, DIP, STRIKE, ...
        xP(idx_x) + dx, yP(idx_x), tpP(idx_x), looks(idx_x,:), nu, fault_type);
    Gp_xm = generate_green_p(Xc, Yc, Zc, ll, ww, DIP, STRIKE, ...
        xP(idx_x) - dx, yP(idx_x), tpP(idx_x), looks(idx_x,:), nu, fault_type);
    Gx = (Gp_xp - Gp_xm) / (2 * dx);

    % --- S - N direction finite difference ---
    idx_y = n_x + 1 : numel(xP);
    Gp_yp = generate_green_p(Xc, Yc, Zc, ll, ww, DIP, STRIKE, ...
        xP(idx_y), yP(idx_y) + dy, tpP(idx_y), looks(idx_y,:), nu, fault_type);
    Gp_ym = generate_green_p(Xc, Yc, Zc, ll, ww, DIP, STRIKE, ...
        xP(idx_y), yP(idx_y) - dy, tpP(idx_y), looks(idx_y,:), nu, fault_type);
    Gy = (Gp_yp - Gp_ym) / (2 * dy);

    % --- Combine results ---
    Greenp_grd = [Gx; Gy];
end
