%% Pin, Hertz contact, and diagonal parameter sweeps
% Keep this file beside DP1_calcs.m and the three beam helper functions.
% Physics and material inputs come from that original calculation script.
% Only sweep/plot settings belong here, avoiding a second copy of the model.
% Model formulas and existing corner/midpoint torsion approximation are preserved.
% Required MATLAB features: script-local functions, tiledlayout, exportgraphics.

close all;

%% Sweep and plotting inputs
L_d_values = (400:800)/1000;    % Horizontal reach, m: previous 1 mm sweep
d_pin_values = 5:50;           % Pin diameters, mm: previous 1 mm sweep
plot_length = 'L_d';        % Plot horizontal reach L_d; Euler still uses actual member length
image_resolution = 200;
% b/h for the support mean outer width w and outer height h_support.
% Diagonal b/h mean solid-section w_diag and h_diag, respectively.
support_min_b = 4*25.4/1000;          % Minimum support width: 4 in
support_b_values = unique([support_min_b (105:5:300)/1000]); % m
support_h_values = (50:5:250)/1000;   % m
support_thickness_values = [5:1:25, 25.4]/1000; % m, broad sweep starting at 5 mm
support_thickness_plot_values = [5 6.35 9.525 12.7 19.05 25.4]/1000; % m, panels
thickness_panels_per_figure = 2;     % Keep comparisons readable on screen
diagonal_b_values = (6:1:40)/1000;    % m
diagonal_h_values = (10:1:60)/1000;   % m
support_plot_Sy_MPa = 35e3*0.006894757293; % 35 ksi = 241.3165 MPa
geometry_FOS = 2;
deflection_limit_mm = 0.7;   % Maximum allowable total CG displacement
% All fixed dimensions come from final_design_inputs() at the VERY END.
% Each selected point is rerun exactly, rather than snapped to the grid.
figure_outputs = struct('handle', {}, 'name', {});
fos_diameter_max = 16;        % Boundary-map zoom; only feasible diameters retained
reference_plane_alpha = 0.07; % Transparent stress-limit planes
hertz_plot_part = 'diagonal';  % 'bracket', 'support', or 'diagonal' HOLE material
% Hertz plots use this part's computed contacts. FOS references use the
% LOWER yield strength of the actual pin and selected hole material.
% The pressure calculation still uses BOTH materials' elastic properties.

script_dir = fileparts(mfilename('fullpath'));
if isempty(script_dir), script_dir = pwd; end
original_file = fullfile(script_dir, 'DP1_calcs.m');
assert(isfile(original_file), 'Keep the two MATLAB files in the same folder.');
% Keep helper functions beside this file; add the folder for standalone runs.
addpath(script_dir);
base = evaluate_original(original_file, [], []);
final_pick = final_design_inputs();
edge_distance_factor = 2;  % User requirement: hole CENTER to edge >= 2*d_hole
required_section_factor = 2*edge_distance_factor; % Centered hole: height >= 4*d_hole
% Include selected inputs exactly in every relevant grid.
L_d_values = unique([L_d_values base.L_d]);
d_pin_values = unique([d_pin_values base.d_pin]);
assert(base.w >= support_min_b, 'Final support width must be at least 4 inches.');
assert(required_section_factor*(base.d_pin+base.hole_tol) <= 1000*min(base.h_diag,base.h_support), ...
    'Final section heights must provide at least 2 hole diameters from center to each edge.');
% The pin runs through the member thickness; that thickness need NOT exceed
% the hole diameter. The in-plane HEIGHT must satisfy the 2d edge rule.
d_pin_values = d_pin_values(required_section_factor*(d_pin_values+base.hole_tol) ...
    <= 1000*min(base.h_diag,base.h_support));
assert(~isempty(d_pin_values), 'No pin diameter meets the selected 2d center-to-edge requirement.');
support_b_values = unique([support_b_values base.w]);
support_b_values = support_b_values(support_b_values >= support_min_b);
support_h_values = unique([support_h_values base.h_support]);
support_thickness_values = unique([support_thickness_values base.th_support]);
support_thickness_plot_values = unique([support_thickness_plot_values base.th_support]);
support_thickness_values = unique([support_thickness_values support_thickness_plot_values]);
[~, support_thickness_plot_indices] = ismember(support_thickness_plot_values, support_thickness_values);
assert(all(isfinite(support_thickness_values)) && all(support_thickness_values > 0), ...
    'Thickness sweep values must be positive and finite.');
diagonal_b_values = unique([diagonal_b_values base.w_diag]);
diagonal_h_values = unique([diagonal_h_values base.h_diag]);
assert(2*base.th_support < min(base.w,base.h_support), ...
    'Final support thickness must leave a hollow interior.');
FOS = base.FOS_target;         % Required FOS = 2; not an extra strength divisor
FOS_levels = [FOS base.FOS_contact_target]; % Show BOTH 2 and conservative 3
reference_Sy_MPa = base.reference_Sy_MPa;
reference_material_names = base.reference_material_names;
reference_material_short = base.reference_material_short;
reference_colors = base.reference_colors;
hertz_selected_case_mask = strcmp(base.contact_hole_part, hertz_plot_part);
assert(any(hertz_selected_case_mask), 'No contact cases match hertz_plot_part.');
hertz_first_case = find(hertz_selected_case_mask, 1);
hertz_selected_material_name = base.contact_material_names{hertz_first_case};
hertz_selected_Sy_MPa = base.hole_Sy_MPa_cases(hertz_first_case);
hertz_selected_E_Pa = base.hole_E_cases(hertz_first_case);
hertz_selected_v = base.hole_v_cases(hertz_first_case);
assert(all(base.hole_E_cases(hertz_selected_case_mask) == hertz_selected_E_Pa) ...
    && all(base.hole_v_cases(hertz_selected_case_mask) == hertz_selected_v), ...
    'Selected Hertz cases must have one consistent material.');
assert(isequaln(base.hole_Sy_MPa_cases(hertz_selected_case_mask), ...
    repmat(hertz_selected_Sy_MPa, size(base.hole_Sy_MPa_cases(hertz_selected_case_mask)))), ...
    'Selected Hertz cases must have one consistent yield strength.');
% The weaker material controls BOTH the subsurface-shear and pressure limits.
% Keep unknown strength explicit instead of allowing min() to omit NaN.
hertz_governing_Sy_MPa = NaN;
hertz_governing_material_name = 'Unspecified pin or hole yield';
if isfinite(base.Sy_pin_MPa) && isfinite(hertz_selected_Sy_MPa)
    hertz_governing_Sy_MPa = min(base.Sy_pin_MPa, hertz_selected_Sy_MPa);
    if hertz_selected_Sy_MPa <= base.Sy_pin_MPa
        hertz_governing_material_name = hertz_selected_material_name;
    else
        hertz_governing_material_name = 'Pin material';
    end
end
pi_value = base.pi_value;
pin_names = base.pin_names;
pin_shear_max_FOS = base.Sy_pin_MPa/(2*FOS);
n_lengths = numel(L_d_values);
n_diameters = numel(d_pin_values);
n_contacts = numel(base.contact_pin_index);
L_diag_values = sqrt(L_d_values.^2 + base.L_y^2);
[L_d_grid, d_pin_grid] = meshgrid(L_d_values, d_pin_values);
[L_diag_grid, ~] = meshgrid(L_diag_values, d_pin_values);

length_values = L_d_values;
length_grid = L_d_grid;
length_label = 'Horizontal diagonal reach L_d (m)';

%% Storage
reaction_solution = zeros(9, n_lengths);
pin_reactions_arr = zeros(6, n_lengths);
max_pin_reaction_mag = zeros(1, n_lengths);
governing_pin_index = zeros(1, n_lengths);
equilibrium_residual = zeros(1, n_lengths);
pin_shear_by_pin_MPa = zeros(6, n_diameters, n_lengths);
contact_p_max_by_case_MPa = zeros(n_contacts, n_diameters, n_lengths);
contact_tau_by_case_MPa = zeros(n_contacts, n_diameters, n_lengths);
contact_b_by_case_m = zeros(n_contacts, n_diameters, n_lengths);
contact_FOS_pin_by_case = zeros(n_contacts, n_diameters, n_lengths);
contact_FOS_hole_by_case = zeros(n_contacts, n_diameters, n_lengths);
contact_FOS_by_case = zeros(n_contacts, n_diameters, n_lengths);
axial_end_forces = zeros(2, 2, n_lengths);
axial_stress_by_diag_MPa = zeros(2, n_lengths);
sigma_bending_MPa = zeros(1, n_lengths);
sigma_diag_by_member_MPa = zeros(2, n_lengths);
sigma_max_diag_MPa = zeros(1, n_lengths);
FOS_diag_yield_by_member = zeros(2, n_lengths);
FOS_diag_yield = zeros(1, n_lengths);
P_cr_diag = zeros(1, n_lengths);
compression_by_diag = zeros(2, n_lengths);
FOS_diag_buckling_by_member = zeros(2, n_lengths);
FOS_diag_buckling = zeros(1, n_lengths);

%% Diameter scaling (exact for the supplied fixed-clearance model)
% Solve the original matrix once per length at the original pin diameter.
% F/A scales as 1/d^2. With fixed hole-minus-pin clearance, Hertz pressure
% scales as sqrt(1/d_pin - 1/d_hole), and b scales by the inverse factor.
d_pin_m_values = d_pin_values/1000;
d_hole_m_values = (d_pin_values + base.hole_tol)/1000;
curvature_values = 1./d_pin_m_values - 1./d_hole_m_values;
curvature_base = 1/(base.d_pin/1000) ...
    - 1/((base.d_pin + base.hole_tol)/1000);
pressure_scale = sqrt(curvature_values/curvature_base);
width_scale = 1./pressure_scale;
direct_shear_scale = (base.d_pin./d_pin_values).^2;

for i = 1:n_lengths
    c = evaluate_original(original_file, L_d_values(i), base.d_pin);
    reaction_solution(:, i) = c.x;
    pin_reactions_arr(:, i) = c.pin_reactions(:);
    max_pin_reaction_mag(i) = c.max_pin_reaction_mag;
    governing_pin_index(i) = c.max_pin_index;
    equilibrium_residual(i) = c.equilibrium_residual;
    pin_shear_by_pin_MPa(:, :, i) = c.pin_shear_stress_MPa(:)*direct_shear_scale;
    contact_p_max_by_case_MPa(:, :, i) = c.p_max_cases_MPa(:)*pressure_scale;
    contact_tau_by_case_MPa(:, :, i) = c.tau_contact_cases_MPa(:)*pressure_scale;
    contact_b_by_case_m(:, :, i) = c.b_cases_m(:)*width_scale;
    contact_FOS_pin_by_case(:, :, i) = c.FOS_contact_pin_cases(:)*width_scale;
    contact_FOS_hole_by_case(:, :, i) = c.FOS_contact_hole_cases(:)*width_scale;
    contact_FOS_by_case(:, :, i) = c.FOS_contact_cases(:)*width_scale;
    axial_end_forces(:, :, i) = c.axial_end_forces;
    axial_stress_by_diag_MPa(:, i) = c.axial_stress_by_diag_MPa(:);
    sigma_bending_MPa(i) = c.sigma_max_bend_MPa;
    sigma_diag_by_member_MPa(:, i) = c.sigma_diag_by_member_MPa(:);
    sigma_max_diag_MPa(i) = c.sigma_max_diag_MPa;
    FOS_diag_yield_by_member(:, i) = c.FOS_diag_yield_by_member(:);
    FOS_diag_yield(i) = c.FOS_diag_yield;
    P_cr_diag(i) = c.P_cr_diag;
    compression_by_diag(:, i) = c.compression_by_diag(:);
    FOS_diag_buckling_by_member(:, i) = c.FOS_diag_buckling_by_member(:);
    FOS_diag_buckling(i) = c.FOS_diag_buckling;
end

% Matrix orientation: rows = d_pin; columns = length.
[max_pin_reaction_stress, shear_governing] = max(pin_shear_by_pin_MPa, [], 1);
max_pin_reaction_stress = reshape(max_pin_reaction_stress, n_diameters, n_lengths);
shear_governing_pin_index = reshape(shear_governing, n_diameters, n_lengths);
FOS_pin_shear = (base.Sy_pin_MPa/2)./max_pin_reaction_stress;
[hertz_p_max_MPa, contact_pressure_governing] = max(contact_p_max_by_case_MPa, [], 1);
hertz_p_max_MPa = reshape(hertz_p_max_MPa, n_diameters, n_lengths);
hertz_pressure_case_index = reshape(contact_pressure_governing, n_diameters, n_lengths);
tau_contact_MPa = reshape(max(contact_tau_by_case_MPa, [], 1), n_diameters, n_lengths);
FOS_contact = NaN(n_diameters, n_lengths);
contact_FOS_case_index = NaN(n_diameters, n_lengths);
if all(base.contact_strength_known)
    [minimum_contact_FOS, contact_FOS_governing] = min(contact_FOS_by_case, [], 1);
    FOS_contact = reshape(minimum_contact_FOS, n_diameters, n_lengths);
    contact_FOS_case_index = reshape(contact_FOS_governing, n_diameters, n_lengths);
end
bracket_cases = strcmp(base.contact_hole_part, 'bracket');
FOS_contact_brackets = reshape(min(contact_FOS_by_case(bracket_cases, :, :), [], 1), ...
    n_diameters, n_lengths);
FOS_contact_pin = reshape(min(contact_FOS_pin_by_case, [], 1), n_diameters, n_lengths);
hertz_p_selected_MPa = reshape(max(contact_p_max_by_case_MPa(hertz_selected_case_mask, :, :), [], 1), ...
    n_diameters, n_lengths);
hertz_tau_selected_MPa = base.hertz_shear_factor*hertz_p_selected_MPa;
FOS_hertz_selected = (hertz_governing_Sy_MPa/2)./hertz_tau_selected_MPa;
passes_FOS = FOS_pin_shear >= FOS;
% NaN flags preserve the distinction between unknown and failed ratings.
passes_contact_FOS = NaN(size(FOS_contact));
known_contact = isfinite(FOS_contact);
passes_contact_FOS(known_contact) = FOS_contact(known_contact) >= base.FOS_contact_target;
passes_diag_yield_FOS = NaN(size(FOS_diag_yield));
known_diagonal = isfinite(FOS_diag_yield);
passes_diag_yield_FOS(known_diagonal) = FOS_diag_yield(known_diagonal) >= FOS;
passes_diag_buckling_FOS = FOS_diag_buckling >= FOS;
hertz_halfwidth_over_pin_radius = zeros(n_contacts, n_diameters, n_lengths);
for i = 1:n_lengths
    hertz_halfwidth_over_pin_radius(:, :, i) = contact_b_by_case_m(:, :, i) ...
        ./ repmat(d_pin_m_values/2, n_contacts, 1);
end
% Centered-hole height screen: center-to-edge >= 2*d_hole.
% Longitudinal end distances require the actual lug/member end geometry.
hole_fits_diag_height = required_section_factor*(d_pin_grid + base.hole_tol) <= 1000*base.h_diag;

%% Material reference checks at the SAME computed stress/load field
% These compare strengths; they do not switch the actual beam material or
% recalculate self-weight. The 36 ksi reference belongs to the brackets.
n_reference_materials = numel(reference_Sy_MPa);
FOS_pin_shear_ref = zeros(n_reference_materials, n_diameters, n_lengths);
FOS_diag_yield_ref = zeros(n_reference_materials, n_lengths);
FOS_diag_buckling_ref = zeros(n_reference_materials, n_lengths);
for material = 1:n_reference_materials
    FOS_pin_shear_ref(material, :, :) = ...
        reshape((reference_Sy_MPa(material)/2)./max_pin_reaction_stress, ...
        1, n_diameters, n_lengths);
    FOS_diag_yield_ref(material, :) = reference_Sy_MPa(material)./sigma_max_diag_MPa;
    FOS_diag_buckling_ref(material, :) = ...
        FOS_diag_buckling*(base.reference_E_Pa(material)/base.E_diag);
end
% Rows = reference material; columns = target FOS (2, 3).
allowable_shear_ref_MPa = (reference_Sy_MPa(:)/2)./FOS_levels;
% Hertz references use the WEAKER actual material in the selected contact pair.
% tau_allow = Sy/(2*FOS); p0_allow = tau_allow/(tau_max/p0).
allowable_hertz_shear_selected_MPa = (hertz_governing_Sy_MPa/2)./FOS_levels;
allowable_hertz_pressure_selected_MPa = ...
    allowable_hertz_shear_selected_MPa/base.hertz_shear_factor;
allowable_normal_ref_MPa = reference_Sy_MPa(:)./FOS_levels;

%% Checks
assert(all(isfinite(max_pin_reaction_stress(:))) && all(max_pin_reaction_stress(:) > 0));
assert(all(isfinite(hertz_p_selected_MPa(:))) && all(hertz_p_selected_MPa(:) > 0));
assert(max(equilibrium_residual) < 1e-7, 'Equilibrium residual exceeds tolerance.');
% Compare a direct rerun of the original script with the scaling calculation.
check_i = round(n_lengths/2);
check_j = round(n_diameters/2);
check = evaluate_original(original_file, L_d_values(check_i), d_pin_values(check_j));
assert(abs(max(check.p_max_cases_MPa) - hertz_p_max_MPa(check_j, check_i)) < 1e-8);
assert(abs(check.FOS_pin_shear - FOS_pin_shear(check_j, check_i)) < 1e-8);
assert(isequaln(check.FOS_contact, FOS_contact(check_j, check_i)) || ...
    abs(check.FOS_contact - FOS_contact(check_j, check_i)) < 1e-8);

%% Support outer width/height/thickness sweep: same calculation script
% Self-weight, stiffness, torque lever arms, and reactions are recomputed.
% Arrays below are indexed (height, width, thickness).
[support_b_grid, support_h_grid] = meshgrid(support_b_values, support_h_values);
assert(ismember(thickness_panels_per_figure, [1 2]), ...
    'Choose one or two thickness panels per figure for readable windows.');
support_sweep_size = [size(support_b_grid) numel(support_thickness_values)];
support_geometry_valid_by_thickness = false(support_sweep_size);
support_VM_by_thickness_MPa = NaN(support_sweep_size);
support_deflection_by_thickness_mm = NaN(support_sweep_size);
support_twist_by_thickness_rad = NaN(support_sweep_size);
for thickness_index = 1:numel(support_thickness_values)
    current_thickness = support_thickness_values(thickness_index);
    valid = support_b_grid >= support_min_b ...
        & support_b_grid > 2*current_thickness ...
        & support_h_grid > 2*current_thickness ...
        & support_h_grid >= required_section_factor*(base.d_pin+base.hole_tol)/1000;
    support_geometry_valid_by_thickness(:,:,thickness_index) = valid;
    for row = 1:numel(support_h_values)
        for column = 1:numel(support_b_values)
            if ~valid(row,column), continue; end
            geometry = struct('w', support_b_values(column), ...
                'h_support', support_h_values(row), 'th_support', current_thickness);
            c = evaluate_original(original_file, base.L_d, base.d_pin, geometry);
            support_VM_by_thickness_MPa(row,column,thickness_index) = c.support_max_VM_MPa;
            support_deflection_by_thickness_mm(row,column,thickness_index) = c.def_magnitude_mm;
            support_twist_by_thickness_rad(row,column,thickness_index) = c.alpha_sh_torsion_rad;
        end
    end
    fprintf('Support thickness %.4g mm: slice %d of %d complete.\n', ...
        current_thickness*1000, thickness_index, numel(support_thickness_values));
end
support_twist_by_thickness_deg = abs(support_twist_by_thickness_rad)*180/pi;
support_FOS_by_thickness = support_plot_Sy_MPa./support_VM_by_thickness_MPa;
support_deflection_passes_by_thickness = support_geometry_valid_by_thickness ...
    & support_deflection_by_thickness_mm <= deflection_limit_mm;
[~, support_selected_thickness_index] = min(abs(support_thickness_values-base.th_support));

% Preserve the earlier 2-D outputs as the selected-thickness slice.
support_geometry_valid = support_geometry_valid_by_thickness(:,:,support_selected_thickness_index);
support_VM_grid_MPa = support_VM_by_thickness_MPa(:,:,support_selected_thickness_index);
support_deflection_grid_mm = support_deflection_by_thickness_mm(:,:,support_selected_thickness_index);
support_twist_grid_rad = support_twist_by_thickness_rad(:,:,support_selected_thickness_index);
support_twist_grid_deg = support_twist_by_thickness_deg(:,:,support_selected_thickness_index);
support_FOS_grid = support_FOS_by_thickness(:,:,support_selected_thickness_index);
support_deflection_passes = support_deflection_passes_by_thickness(:,:,support_selected_thickness_index);
support_VM_allowable_MPa = support_plot_Sy_MPa/geometry_FOS;

% Common scales make the thickness panels directly comparable.
support_stress_color_max = max(support_VM_by_thickness_MPa(support_geometry_valid_by_thickness));
support_stress_zmax = 1.08*max(support_stress_color_max, support_VM_allowable_MPa);
support_deflection_zmax = 1.08*max([deflection_limit_mm; ...
    support_deflection_by_thickness_mm(support_geometry_valid_by_thickness)]);
support_twist_color_max = max([eps; support_twist_by_thickness_deg(support_geometry_valid_by_thickness)]);
support_twist_zmax = 1.08*support_twist_color_max;

%% Axial diagonal width/height sweep: Euler buckling about the weak axis
% Axes use b/h, with final L_d fixed. Actual L_diag enters Euler's formula.
% Demand is recomputed because changing section size changes self-weight.
[diagonal_b_grid, diagonal_h_grid] = meshgrid(diagonal_b_values, diagonal_h_values);
diagonal_geometry_valid = diagonal_h_grid >= required_section_factor*(base.d_pin+base.hole_tol)/1000;
% Centered holes: minimum height is 4*d_hole for 2d center-to-edge.
% w_diag is along the pin axis, so it is plate thickness, not hole height.
diagonal_Pcr_grid_N = NaN(size(diagonal_b_grid));
diagonal_compression_grid_N = NaN(size(diagonal_b_grid));
diagonal_buckling_FOS_grid = NaN(size(diagonal_b_grid));
for row = 1:numel(diagonal_h_values)
    for column = 1:numel(diagonal_b_values)
        if ~diagonal_geometry_valid(row,column), continue; end
        geometry = struct('w_diag', diagonal_b_values(column), ...
            'h_diag', diagonal_h_values(row));
        c = evaluate_original(original_file, base.L_d, base.d_pin, geometry);
        diagonal_Pcr_grid_N(row,column) = c.P_cr_diag;
        diagonal_compression_grid_N(row,column) = max(c.compression_by_diag);
        diagonal_buckling_FOS_grid(row,column) = c.FOS_diag_buckling;
    end
end
diagonal_required_Pcr_grid_N = geometry_FOS*diagonal_compression_grid_N;

%% Contact half-width / pin DIAMETER (not radius)
contact_halfwidth_over_diameter = zeros(size(contact_b_by_case_m));
for i = 1:n_lengths
    contact_halfwidth_over_diameter(:,:,i) = contact_b_by_case_m(:,:,i) ...
        ./repmat(d_pin_m_values, n_contacts, 1);
end
hertz_bd_selected = reshape(max( ...
    contact_halfwidth_over_diameter(hertz_selected_case_mask,:,:), [], 1), ...
    n_diameters, n_lengths);
hertz_bd_all = reshape(max(contact_halfwidth_over_diameter, [], 1), ...
    n_diameters, n_lengths);

%% Validate the new grids and the exact selected-design intersections
assert(all(isfinite(support_VM_by_thickness_MPa(support_geometry_valid_by_thickness))));
assert(all(isfinite(support_deflection_by_thickness_mm(support_geometry_valid_by_thickness))));
assert(all(isfinite(support_twist_by_thickness_rad(support_geometry_valid_by_thickness))));
assert(all(isfinite(support_VM_grid_MPa(support_geometry_valid))) && ...
    all(support_VM_grid_MPa(support_geometry_valid) >= 0));
assert(all(isfinite(support_deflection_grid_mm(support_geometry_valid))) && ...
    all(support_deflection_grid_mm(support_geometry_valid) >= 0));
assert(all(isfinite(support_twist_grid_rad(support_geometry_valid))));
assert(all(isfinite(diagonal_Pcr_grid_N(diagonal_geometry_valid))) && ...
    all(diagonal_Pcr_grid_N(diagonal_geometry_valid) > 0));
assert(all(isfinite(hertz_bd_selected(:))) && all(hertz_bd_selected(:) >= 0));
[~, picked_b] = min(abs(support_b_values-base.w));
[~, picked_h] = min(abs(support_h_values-base.h_support));
assert(abs(support_VM_grid_MPa(picked_h,picked_b)-base.support_max_VM_MPa) ...
    < 1e-8*max(1,base.support_max_VM_MPa));
assert(abs(support_deflection_grid_mm(picked_h,picked_b)-base.def_magnitude_mm) ...
    < 1e-8*max(1,base.def_magnitude_mm));
assert(abs(support_twist_grid_rad(picked_h,picked_b)-base.alpha_sh_torsion_rad) ...
    < 1e-8*max(1,abs(base.alpha_sh_torsion_rad)));
[~, picked_b_diag] = min(abs(diagonal_b_values-base.w_diag));
[~, picked_h_diag] = min(abs(diagonal_h_values-base.h_diag));
assert(abs(diagonal_Pcr_grid_N(picked_h_diag,picked_b_diag)-base.P_cr_diag) ...
    < 1e-8*max(1,base.P_cr_diag));

%% Save all numerical data (actual FOS is NEVER capped in storage)
model_inputs = base;
save(fullfile(script_dir, 'pin_stress_results.mat'), ...
    'L_d_values', 'L_diag_values', 'd_pin_values', 'L_d_grid', 'L_diag_grid', ...
    'd_pin_grid', 'max_pin_reaction_stress', 'max_pin_reaction_mag', ...
    'pin_reactions_arr', 'reaction_solution', 'governing_pin_index', ...
    'shear_governing_pin_index', 'pin_names', 'FOS', 'pin_shear_max_FOS', ...
    'pi_value', 'equilibrium_residual', 'pin_shear_by_pin_MPa', ...
    'FOS_pin_shear', 'passes_FOS', 'contact_p_max_by_case_MPa', ...
    'contact_tau_by_case_MPa', 'contact_b_by_case_m', ...
    'contact_FOS_pin_by_case', 'contact_FOS_hole_by_case', 'contact_FOS_by_case', ...
    'hertz_p_max_MPa', 'tau_contact_MPa', 'FOS_contact', ...
    'hertz_pressure_case_index', 'contact_FOS_case_index', 'passes_contact_FOS', ...
    'axial_end_forces', 'axial_stress_by_diag_MPa', 'sigma_bending_MPa', ...
    'sigma_diag_by_member_MPa', 'sigma_max_diag_MPa', ...
    'FOS_diag_yield_by_member', 'FOS_diag_yield', 'passes_diag_yield_FOS', ...
    'P_cr_diag', 'compression_by_diag', 'FOS_diag_buckling_by_member', ...
    'FOS_diag_buckling', 'passes_diag_buckling_FOS', ...
    'hertz_halfwidth_over_pin_radius', 'hole_fits_diag_height', 'model_inputs', ...
    'FOS_levels', 'reference_Sy_MPa', 'reference_material_names', ...
    'reference_material_short', 'reference_colors', 'FOS_pin_shear_ref', ...
    'FOS_diag_yield_ref', 'FOS_diag_buckling_ref', 'allowable_shear_ref_MPa', ...
    'allowable_normal_ref_MPa', 'FOS_contact_pin', 'FOS_contact_brackets', ...
    'hertz_plot_part', 'hertz_selected_case_mask', 'hertz_selected_material_name', ...
    'hertz_selected_Sy_MPa', 'hertz_selected_E_Pa', 'hertz_selected_v', ...
    'hertz_governing_Sy_MPa', 'hertz_governing_material_name', ...
    'hertz_p_selected_MPa', 'hertz_tau_selected_MPa', 'FOS_hertz_selected', ...
    'allowable_hertz_shear_selected_MPa', 'allowable_hertz_pressure_selected_MPa', ...
    'support_b_values', 'support_h_values', 'support_b_grid', 'support_h_grid', ...
    'support_thickness_values', 'support_selected_thickness_index', ...
    'support_thickness_plot_values', 'support_thickness_plot_indices', ...
    'support_geometry_valid_by_thickness', 'support_VM_by_thickness_MPa', ...
    'support_deflection_by_thickness_mm', 'support_twist_by_thickness_rad', ...
    'support_twist_by_thickness_deg', 'support_FOS_by_thickness', ...
    'support_deflection_passes_by_thickness', ...
    'support_geometry_valid', 'diagonal_geometry_valid', 'support_min_b', ...
    'edge_distance_factor', 'required_section_factor', ...
    'support_VM_grid_MPa', 'support_deflection_grid_mm', ...
    'deflection_limit_mm', 'support_deflection_passes', ...
    'support_twist_grid_rad', 'support_twist_grid_deg', 'support_FOS_grid', ...
    'support_plot_Sy_MPa', 'support_VM_allowable_MPa', 'geometry_FOS', ...
    'diagonal_b_values', 'diagonal_h_values', 'diagonal_b_grid', 'diagonal_h_grid', ...
    'diagonal_Pcr_grid_N', 'diagonal_compression_grid_N', ...
    'diagonal_required_Pcr_grid_N', 'diagonal_buckling_FOS_grid', ...
    'contact_halfwidth_over_diameter', 'hertz_bd_selected', 'hertz_bd_all', ...
    'final_pick');

%% Direct pin shear: material-specific FOS 2 and 3 limit planes
fig = fitted_figure();
ax = axes(fig);
stress_surface(ax, length_grid, d_pin_grid, max_pin_reaction_stress, ...
    length_label, 'Maximum nominal pin shear stress (MPa)', 'Direct pin shear');
stress_limit_planes(ax, length_grid, d_pin_grid, max_pin_reaction_stress, ...
    allowable_shear_ref_MPa, reference_material_short, reference_colors, ...
    FOS_levels, reference_plane_alpha);
mark_final_3d(ax, base.L_d, base.d_pin, base.max_pin_reaction_stress_MPa, ...
    sprintf('Final: shear = %.4g MPa', base.max_pin_reaction_stress_MPa));
figure_outputs(end+1) = struct('handle', fig, 'name', 'pin_stress_surface.png');

%% Hertz: lower pin/hole strength controls both FOS reference plots
fig = fitted_figure();
tl = tiledlayout(fig, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
hertz_color = [0.06 0.39 0.70];
ax = nexttile(tl);
stress_surface(ax, length_grid, d_pin_grid, hertz_p_selected_MPa, length_label, ...
    'Maximum pressure p_0 (MPa)', 'Hertz peak contact pressure');
stress_limit_planes(ax, length_grid, d_pin_grid, hertz_p_selected_MPa, ...
    allowable_hertz_pressure_selected_MPa, {'Weaker material'}, ...
    hertz_color, FOS_levels, reference_plane_alpha);
mark_final_3d(ax, base.L_d, base.d_pin, max(base.p_max_cases_MPa(hertz_selected_case_mask)), ...
    sprintf('Final: p0 = %.4g MPa', max(base.p_max_cases_MPa(hertz_selected_case_mask))));
ax = nexttile(tl);
stress_surface(ax, length_grid, d_pin_grid, hertz_tau_selected_MPa, length_label, ...
    'Maximum contact-induced shear (MPa)', 'Hertz subsurface shear');
stress_limit_planes(ax, length_grid, d_pin_grid, hertz_tau_selected_MPa, ...
    allowable_hertz_shear_selected_MPa, {'Weaker material'}, ...
    hertz_color, FOS_levels, reference_plane_alpha);
mark_final_3d(ax, base.L_d, base.d_pin, max(base.tau_contact_cases_MPa(hertz_selected_case_mask)), ...
    sprintf('Final: shear = %.4g MPa', max(base.tau_contact_cases_MPa(hertz_selected_case_mask))));
title(tl, sprintf('Hertz limits: weaker pin/hole yield = %.1f MPa', ...
    hertz_governing_Sy_MPa), ...
    'Interpreter', 'none');
xlabel(tl, sprintf(['Pin: E = %.1f GPa, nu = %.2f | tau = %.2f p_0 | ' ...
    'tau_allow = Sy/(2 FOS); p0_allow = tau_allow/%.2f'], ...
    base.pin_E/1e9, base.pin_v, base.hertz_shear_factor, base.hertz_shear_factor));
figure_outputs(end+1) = struct('handle', fig, 'name', 'hertz_contact_stress.png');

%% Diagonal normal stress: section size and shape stay constant
fig = fitted_figure();
ax = axes(fig);
plot(ax, length_values, sigma_max_diag_MPa, 'k-', 'LineWidth', 2.4, ...
    'DisplayName', 'Maximum axial + bending envelope');
hold(ax, 'on');
plot(ax, length_values, max(axial_stress_by_diag_MPa, [], 1), ...
    '--', 'Color', [0.45 0.45 0.45], 'LineWidth', 1.4, ...
    'DisplayName', 'Maximum axial magnitude');
plot(ax, length_values, sigma_bending_MPa, ':', ...
    'Color', [0.45 0.45 0.45], 'LineWidth', 1.4, ...
    'DisplayName', 'Bending contribution');
line_styles = {'-', '--'};
for material = 1:n_reference_materials
    for target = 1:numel(FOS_levels)
        yline(ax, allowable_normal_ref_MPa(material, target), line_styles{target}, ...
            'Color', reference_colors(material, :), 'LineWidth', 1.6, ...
            'DisplayName', sprintf('%s FOS %g: %.1f MPa', ...
            reference_material_short{material}, FOS_levels(target), ...
            allowable_normal_ref_MPa(material, target)));
    end
end
hold(ax, 'off');
grid(ax, 'on');
xlabel(ax, length_label);
ylabel(ax, 'Normal stress (MPa)');
title(ax, sprintf(['Diagonal stress | fixed %.0f x %.0f mm solid rectangle\n' ...
    'Reference limits compare strengths at the current computed stress'], ...
    1000*base.w_diag, 1000*base.h_diag));
legend(ax, 'Location', 'southoutside', 'NumColumns', 2, 'Interpreter', 'none');
xlim(ax, [min(length_values) max(length_values)]);
ylim(ax, [0 1.10*max([sigma_max_diag_MPa allowable_normal_ref_MPa(:).'])]);
if ~isfinite(base.Sy_diag_MPa)
    text(ax, 0.02, 0.96, 'Actual diagonal yield/FOS unset: enter Sy_diag_MPa', ...
        'Units', 'normalized', 'VerticalAlignment', 'top', ...
        'BackgroundColor', 'w', 'Margin', 3);
end
mark_final_2d(ax, base.L_d, base.sigma_max_diag_MPa, ...
    sprintf('Final: stress = %.4g MPa', base.sigma_max_diag_MPa));
figure_outputs(end+1) = struct('handle', fig, 'name', 'diagonal_stress.png');

%% FOS plots: explicit 2 and 3 boundaries, with uncapped material checks
fig = fitted_figure();
tl = tiledlayout(fig, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
ax = nexttile(tl);
fos_boundary_map(ax, length_grid, d_pin_grid, FOS_pin_shear_ref, FOS_levels, ...
    reference_material_short, reference_colors, length_label, ...
    'Direct pin shear: strength references', fos_diameter_max);
mark_final_2d(ax, base.L_d, base.d_pin, sprintf('Final: FOS = %.4g', base.FOS_pin_shear));
ax = nexttile(tl);
fos_boundary_map(ax, length_grid, d_pin_grid, ...
    reshape(FOS_hertz_selected, 1, n_diameters, n_lengths), FOS_levels, ...
    {'Weaker material'}, hertz_color, length_label, ...
    sprintf('Hertz: weaker pin / %s hole', hertz_plot_part), ...
    fos_diameter_max);
mark_final_2d(ax, base.L_d, base.d_pin, sprintf('Final: FOS = %.4g', ...
    (hertz_governing_Sy_MPa/2)/max(base.tau_contact_cases_MPa(hertz_selected_case_mask))));
title(tl, 'Pin/contact FOS: solid = 2, dashed = 3');
figure_outputs(end+1) = struct('handle', fig, 'name', 'pin_contact_fos.png');
fig = fitted_figure();
tl = tiledlayout(fig, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
ax = nexttile(tl);
fos_reference_curves(ax, length_values, FOS_diag_yield_ref, FOS_levels, ...
    reference_material_short, reference_colors, length_label, ...
    'Diagonal yield: strength references', FOS_diag_yield);
mark_final_2d(ax, base.L_d, base.FOS_diag_yield, ...
    sprintf('Final: FOS = %.4g', base.FOS_diag_yield));
if ~isfinite(base.Sy_diag_MPa)
    text(ax, 0.02, 0.96, 'Actual diagonal FOS unset; curves are reference checks', ...
        'Units', 'normalized', 'VerticalAlignment', 'top', ...
        'BackgroundColor', 'w', 'Margin', 3);
end
ax = nexttile(tl);
fos_reference_curves(ax, length_values, FOS_diag_buckling_ref, FOS_levels, ...
    reference_material_short, reference_colors, length_label, ...
    sprintf('Euler buckling: E references (K = %g)', base.K_diag), []);
mark_final_2d(ax, base.L_d, base.FOS_diag_buckling, ...
    sprintf('Final: FOS = %.4g', base.FOS_diag_buckling));
title(tl, sprintf('Diagonal FOS | fixed %.1f x %.1f mm', ...
    1000*base.w_diag, 1000*base.h_diag));
figure_outputs(end+1) = struct('handle', fig, 'name', 'fos_plots.png');

%% Support comparisons: multiple fixed-thickness surfaces per metric
% Full thickness data are saved above. Panel thicknesses are set near the top.
% Matching axes/color limits are used across every thickness of each metric.
support_metrics = {'von_mises', 'total_deflection', 'torsional_angle'};
for metric_index = 1:numel(support_metrics)
    metric = support_metrics{metric_index};
    switch metric
        case 'von_mises'
            group_title = 'Support von Mises: thickness comparison';
        case 'total_deflection'
            group_title = 'Total CG deflection: gray exceeds 0.7 mm';
        case 'torsional_angle'
            group_title = 'Support torsional angle: thickness comparison';
    end
    page = 0;
    for first_panel = 1:thickness_panels_per_figure:numel(support_thickness_plot_indices)
        page = page+1;
        last_panel = min(first_panel+thickness_panels_per_figure-1, ...
            numel(support_thickness_plot_indices));
        page_indices = support_thickness_plot_indices(first_panel:last_panel);
        fig = fitted_figure();
        tl = tiledlayout(fig, 1, numel(page_indices), ...
            'TileSpacing', 'compact', 'Padding', 'compact');
        for panel = 1:numel(page_indices)
            thickness_index = page_indices(panel);
            current_thickness = support_thickness_values(thickness_index);
            ax = nexttile(tl);
            valid = support_geometry_valid_by_thickness(:,:,thickness_index);
            if ~any(valid(:))
                text(ax, 0.5, 0.5, 'No feasible hollow section for this thickness', ...
                    'Units', 'normalized', 'HorizontalAlignment', 'center');
                title(ax, sprintf('t = %.4g mm', current_thickness*1000));
                continue;
            end
            switch metric
                case 'von_mises'
                    Z = support_VM_by_thickness_MPa(:,:,thickness_index);
                    geometry_surface(ax, support_b_grid*1000, support_h_grid*1000, Z, ...
                        'Support outer b (mm)', 'Support outer h (mm)', ...
                        'Maximum von Mises (MPa)', '');
                    horizontal_plane(ax, support_b_grid*1000, support_h_grid*1000, ...
                        support_VM_allowable_MPa, ...
                        sprintf('35 ksi / FOS 2 = %.2f MPa', support_VM_allowable_MPa));
                    caxis(ax, [0 support_stress_color_max]);
                    zlim(ax, [0 support_stress_zmax]);
                    final_value = base.support_max_VM_MPa;
                    final_label = sprintf('Final: von Mises = %.4g MPa', final_value);
                case 'total_deflection'
                    Z = support_deflection_by_thickness_mm(:,:,thickness_index);
                    deflection_surface(ax, support_b_grid*1000, support_h_grid*1000, ...
                        Z, deflection_limit_mm);
                    zlim(ax, [0 support_deflection_zmax]);
                    final_value = base.def_magnitude_mm;
                    final_label = sprintf('Final: deflection = %.4g mm', final_value);
                case 'torsional_angle'
                    Z = support_twist_by_thickness_deg(:,:,thickness_index);
                    geometry_surface(ax, support_b_grid*1000, support_h_grid*1000, Z, ...
                        'Support outer b (mm)', 'Support outer h (mm)', ...
                        '|Twist at shaft| (deg)', '');
                    caxis(ax, [0 support_twist_color_max]);
                    zlim(ax, [0 support_twist_zmax]);
                    final_value = abs(base.alpha_sh_torsion_rad)*180/pi;
                    final_label = sprintf('Final: |twist| = %.4g deg', final_value);
            end
            title(ax, sprintf('t = %.4g mm', current_thickness*1000));
            if thickness_index == support_selected_thickness_index
                mark_final_3d(ax, base.w*1000, base.h_support*1000, ...
                    final_value, final_label);
                title(ax, sprintf('t = %.4g mm (selected)', current_thickness*1000));
            end
        end
        title(tl, sprintf('%s | L_d = %.3f m', group_title, base.L_d));
        if page == 1
            output_name = sprintf('support_%s.png', metric);
        else
            output_name = sprintf('support_%s_%02d.png', metric, page);
        end
        figure_outputs(end+1) = struct('handle', fig, 'name', output_name);
    end
end

%% Selected-thickness von Mises map retained for reading the design boundary
fig = fitted_figure(); ax = axes(fig);
contourf(ax, support_b_grid*1000, support_h_grid*1000, support_VM_grid_MPa, ...
    24, 'LineStyle', 'none');
hold(ax, 'on');
if min(support_VM_grid_MPa(support_geometry_valid)) < support_VM_allowable_MPa && ...
        max(support_VM_grid_MPa(support_geometry_valid)) > support_VM_allowable_MPa
    contour(ax, support_b_grid*1000, support_h_grid*1000, ...
        support_VM_grid_MPa, [support_VM_allowable_MPa support_VM_allowable_MPa], ...
        'k-', 'LineWidth', 2);
end
hold(ax, 'off');
xlabel(ax, 'Support outer b (mm)'); ylabel(ax, 'Support outer h (mm)');
title(ax, sprintf('Von Mises map | selected t = %.4g mm | black boundary = FOS 2', ...
    base.th_support*1000));
cb = colorbar(ax); cb.Label.String = 'Maximum von Mises (MPa)';
grid(ax, 'on');
mark_final_2d(ax, base.w*1000, base.h_support*1000, ...
    sprintf('Final: von Mises = %.4g MPa', base.support_max_VM_MPa));
figure_outputs(end+1) = struct('handle', fig, 'name', 'support_von_mises_map.png');

%% Euler buckling versus axial diagonal solid-section b and h
fig = fitted_figure();
tl = tiledlayout(fig, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
ax = nexttile(tl);
geometry_surface(ax, diagonal_b_grid*1000, diagonal_h_grid*1000, ...
    diagonal_Pcr_grid_N/1000, 'Diagonal b (mm)', 'Diagonal h (mm)', ...
    'Euler P_{cr} (kN)', 'Weak-axis Euler capacity');
hold(ax, 'on');
surf(ax, diagonal_b_grid*1000, diagonal_h_grid*1000, ...
    diagonal_required_Pcr_grid_N/1000, 'FaceColor', [0.85 0.22 0.12], ...
    'FaceAlpha', 0.28, 'EdgeColor', 'none', ...
    'DisplayName', 'Required capacity = 2 x compression');
hold(ax, 'off');
legend(ax, 'show', 'Location', 'southoutside', 'FontSize', 9);
mark_final_3d(ax, base.w_diag*1000, base.h_diag*1000, base.P_cr_diag/1000, ...
    sprintf('Final: Pcr = %.4g kN', base.P_cr_diag/1000));
ax = nexttile(tl);
fos_plot = diagonal_buckling_FOS_grid;
fos_plot(~isfinite(fos_plot)) = NaN; % Stored values remain Inf when unloaded.
if any(isfinite(fos_plot(:)))
    geometry_surface(ax, diagonal_b_grid*1000, diagonal_h_grid*1000, ...
        fos_plot, 'Diagonal b (mm)', 'Diagonal h (mm)', ...
        'P_{cr} / compression', 'Achieved buckling FOS');
    horizontal_plane(ax, diagonal_b_grid*1000, diagonal_h_grid*1000, ...
        geometry_FOS, 'Required FOS = 2');
    mark_final_3d(ax, base.w_diag*1000, base.h_diag*1000, base.FOS_diag_buckling, ...
        sprintf('Final: FOS = %.4g', base.FOS_diag_buckling));
else
    text(ax, 0.5, 0.5, 'No compressive demand: buckling FOS is infinite', ...
        'Units', 'normalized', 'HorizontalAlignment', 'center');
    mark_final_2d(ax, base.w_diag*1000, base.h_diag*1000, 'Final: FOS = Inf (no compression)');
end
title(tl, sprintf('Axial diagonal (blank = 2d edge distance fails) | L_d = %.3f m | K = %g | E = %.1f GPa', ...
    base.L_d, base.K_diag, base.E_diag/1e9));
figure_outputs(end+1) = struct('handle', fig, 'name', 'diagonal_geometry_buckling.png');

%% Contact half-width / pin diameter and the requested 0.5 reference
fig = fitted_figure();
tl = tiledlayout(fig, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
ax = nexttile(tl);
geometry_surface(ax, L_d_grid, d_pin_grid, hertz_bd_selected, ...
    length_label, 'd_{pin} (mm)', 'Contact half-width / d_{pin}', ...
    sprintf('b/d: %s contacts', hertz_plot_part));
horizontal_plane(ax, L_d_grid, d_pin_grid, 0.5, 'b/d = 0.5');
final_bd = max(base.b_cases_m(hertz_selected_case_mask))/(base.d_pin/1000);
mark_final_3d(ax, base.L_d, base.d_pin, final_bd, sprintf('Final: b/d = %.4g', final_bd));
ax = nexttile(tl);
[~, selected_length_index] = min(abs(L_d_values-base.L_d));
plot(ax, d_pin_values, hertz_bd_selected(:,selected_length_index), ...
    'k-', 'LineWidth', 2, 'DisplayName', 'Selected hole material contacts');
hold(ax, 'on');
plot(ax, d_pin_values, hertz_bd_all(:,selected_length_index), ...
    '--', 'Color', [0.45 0.45 0.45], 'LineWidth', 1.5, ...
    'DisplayName', 'Envelope of all contacts');
yline(ax, 0.5, 'b-', 'LineWidth', 1.5, 'DisplayName', 'b/d = 0.5');
hold(ax, 'off');
xlabel(ax, 'd_{pin} (mm)'); ylabel(ax, 'Contact half-width / d_{pin}');
title(ax, sprintf('Diameter sweep at selected L_d = %.3f m', base.L_d));
grid(ax, 'on'); legend(ax, 'Location', 'southoutside', 'FontSize', 9);
mark_final_2d(ax, base.d_pin, final_bd, sprintf('Final: b/d = %.4g', final_bd));
title(tl, 'Contact geometry diagnostic: b is contact half-width, d is pin diameter');
figure_outputs(end+1) = struct('handle', fig, 'name', 'contact_halfwidth_over_diameter.png');

%% Export the finished figures, INCLUDING the red final-design markers
for figure_index = 1:numel(figure_outputs)
    exportgraphics(figure_outputs(figure_index).handle, ...
        fullfile(script_dir, figure_outputs(figure_index).name), ...
        'Resolution', image_resolution);
end

%% Summary
fprintf('Evaluated %d length/diameter combinations.\n', n_lengths*n_diameters);
fprintf('L_d: %.3f to %.3f m; actual L_diag: %.3f to %.3f m.\n', ...
    min(L_d_values), max(L_d_values), min(L_diag_values), max(L_diag_values));
fprintf('Constant diagonal section: %.1f x %.1f mm, solid rectangle.\n', ...
    1000*base.w_diag, 1000*base.h_diag);
fprintf('Bracket yield = %.3f MPa (36 ksi); diagonal yield input = %g MPa.\n', ...
    base.Sy_bracket_MPa, base.Sy_diag_MPa);
fprintf('Selected Hertz material: %s; E = %.1f GPa; yield = %g MPa.\n', ...
    hertz_selected_material_name, hertz_selected_E_Pa/1e9, hertz_selected_Sy_MPa);
if isfinite(base.Sy_diag_MPa)
    fprintf('Actual diagonal yield FOS range: %.3f to %.3f.\n', ...
        min(FOS_diag_yield), max(FOS_diag_yield));
else
    fprintf('Actual diagonal yield FOS is unrated; set Sy_diag_MPa.\n');
end
fprintf('Current diagonal Euler buckling FOS range: %.3f to %.3f.\n', ...
    min(FOS_diag_buckling), max(FOS_diag_buckling));
fprintf('Largest Hertz contact half-width / pin radius: %.3f.\n', ...
    max(hertz_halfwidth_over_pin_radius(:)));
fprintf('Hertz limits: tau_allow = Sy/(2*FOS), p0_allow = tau_allow/%.2f.\n', ...
    base.hertz_shear_factor);
fprintf('Controlling Hertz yield = %.3f MPa: %s.\n', ...
    hertz_governing_Sy_MPa, hertz_governing_material_name);
fprintf('Total-deflection limit = %.3f mm; gray surface exceeds it.\n', deflection_limit_mm);
fprintf('All numerical FOS values are saved uncapped; unknown ratings are NaN.\n');
fprintf('Support plot yield = %.3f MPa (35 ksi); FOS 2 limit = %.3f MPa.\n', ...
    support_plot_Sy_MPa, support_VM_allowable_MPa);
fprintf('Support sweep: %d thicknesses from %.4g to %.4g mm; %d plotted.\n', ...
    numel(support_thickness_values), min(support_thickness_values)*1000, ...
    max(support_thickness_values)*1000, numel(support_thickness_plot_indices));
fprintf('Final support: b = %.2f, h = %.2f, t = %.2f mm.\n', ...
    base.w*1000, base.h_support*1000, base.th_support*1000);
fprintf('Final VM = %.4f MPa; total CG displacement = %.4f mm; twist = %.5f deg.\n', ...
    base.support_max_VM_MPa, base.def_magnitude_mm, base.alpha_sh_torsion_rad*180/pi);

%% Local helpers
function result = evaluate_original(original_path, reach, diameter, geometry)
    sweep_override = final_design_inputs();
    sweep_override.verbose = false;
    if ~isempty(reach), sweep_override.L_d = reach; end
    if ~isempty(diameter), sweep_override.d_pin = diameter; end
    if nargin >= 4
        names = fieldnames(geometry);
        for k = 1:numel(names)
            sweep_override.(names{k}) = geometry.(names{k});
        end
    end
    % Capture the uploaded script's unsuppressed displays during large sweeps.
    evalc('run(original_path);');
    result = calc;
    result.max_pin_reaction_stress_MPa = max_pin_reaction_stress;
end

function stress_surface(ax, X, Y, Z, x_label, z_label, plot_title)
    surf(ax, X, Y, Z, 'EdgeColor', 'none', 'FaceColor', 'interp', ...
        'HandleVisibility', 'off');
    xlabel(ax, x_label);
    ylabel(ax, 'd_{pin} (mm)');
    zlabel(ax, z_label);
    title(ax, plot_title);
    shades = linspace(0.94, 0.40, 256).';
    colormap(ax, repmat(shades, 1, 3));
    caxis(ax, [min(Z(:)) max(Z(:))]);
    cb = colorbar(ax);
    cb.Label.String = 'Computed stress (MPa)';
    grid(ax, 'on');
    axis(ax, 'tight');
    view(ax, 45, 28);
end

function stress_limit_planes(ax, X, Y, Z, limits, names, colors, targets, alpha)
    % The stress-axis levels are converted from the applicable failure check:
    % direct/contact shear Sy/(2*N), Hertz p0 Sy/(2*k*N), normal Sy/N.
    styles = {'-', '--'};
    xmin = min(X(:)); xmax = max(X(:));
    ymin = min(Y(:)); ymax = max(Y(:));
    [PX, PY] = meshgrid([xmin xmax], [ymin ymax]);
    handles = gobjects(0);
    labels = {};
    hold(ax, 'on');
    for material = 1:size(limits, 1)
        for target = 1:numel(targets)
            limit = limits(material, target);
            if ~isfinite(limit), continue; end
            surf(ax, PX, PY, limit*ones(2), 'FaceColor', colors(material, :), ...
                'FaceAlpha', alpha, 'EdgeColor', 'none', 'HandleVisibility', 'off');
            plot3(ax, [xmin xmax xmax xmin xmin], [ymin ymin ymax ymax ymin], ...
                limit*ones(1, 5), styles{target}, 'Color', colors(material, :), ...
                'LineWidth', 1.1, 'HandleVisibility', 'off');
            if min(Z(:)) < limit && max(Z(:)) > limit
                contour3(ax, X, Y, Z, [limit limit], styles{target}, ...
                    'LineColor', colors(material, :), 'LineWidth', 2, ...
                    'HandleVisibility', 'off');
            end
            handles(end+1) = plot3(ax, NaN, NaN, NaN, styles{target}, ...
                'Color', colors(material, :), 'LineWidth', 1.8); %#ok<AGROW>
            labels{end+1} = sprintf('%s FOS %g: %.1f MPa', ...
                names{material}, targets(target), limit); %#ok<AGROW>
        end
    end
    hold(ax, 'off');
    xlim(ax, [xmin xmax]); ylim(ax, [ymin ymax]);
    finite_limits = limits(isfinite(limits));
    zlim(ax, [0 1.08*max([Z(:); finite_limits(:)])]);
    if ~isempty(handles)
        legend(ax, handles, labels, 'Location', 'southoutside', ...
            'NumColumns', min(2, numel(handles)), 'Interpreter', 'none', 'FontSize', 9);
    else
        text(ax, 0.02, 0.96, 'Selected yield unspecified: no FOS limits', ...
            'Units', 'normalized', 'VerticalAlignment', 'top');
    end
end

function fos_boundary_map(ax, X, Y, material_values, targets, names, colors, ...
        x_label, plot_title, diameter_max)
    styles = {'-', '--'};
    handles = gobjects(0);
    labels = {};
    hold(ax, 'on');
    for material = 1:size(material_values, 1)
        values = reshape(material_values(material, :, :), size(X));
        if ~any(isfinite(values(:))), continue; end
        for target = 1:numel(targets)
            if min(values(:)) < targets(target) && max(values(:)) > targets(target)
                contour(ax, X, Y, values, [targets(target) targets(target)], ...
                    styles{target}, 'LineColor', colors(material, :), 'LineWidth', 2);
            end
            handles(end+1) = plot(ax, NaN, NaN, styles{target}, ...
                'Color', colors(material, :), 'LineWidth', 2); %#ok<AGROW>
            labels{end+1} = sprintf('%s: FOS %g', names{material}, targets(target)); %#ok<AGROW>
            if min(values(:)) >= targets(target)
                labels{end} = [labels{end} ' (whole sweep meets)'];
            end
        end
    end
    hold(ax, 'off');
    xlabel(ax, x_label);
    ylabel(ax, 'd_{pin} (mm)');
    title(ax, plot_title);
    xlim(ax, [min(X(:)) max(X(:))]);
    ylim(ax, [min(Y(:)) min(max(Y(:)), diameter_max)]);
    grid(ax, 'on');
    if ~isempty(handles)
        legend(ax, handles, labels, 'Location', 'southoutside', ...
            'NumColumns', min(2, numel(handles)), 'Interpreter', 'none', 'FontSize', 9);
        text(ax, 0.02, 0.96, 'Above a boundary meets that material/target FOS', ...
            'Units', 'normalized', 'VerticalAlignment', 'top', ...
            'BackgroundColor', 'w', 'Margin', 3);
    else
        text(ax, 0.5, 0.5, 'Selected yield unspecified: FOS boundaries are unrated', ...
            'Units', 'normalized', 'HorizontalAlignment', 'center');
    end
end

function fos_reference_curves(ax, X, material_values, targets, names, colors, ...
        x_label, plot_title, actual_values)
    hold(ax, 'on');
    for material = 1:size(material_values, 1)
        plot(ax, X, material_values(material, :), '-', ...
            'Color', colors(material, :), 'LineWidth', 2, ...
            'DisplayName', names{material});
    end
    if ~isempty(actual_values) && any(isfinite(actual_values))
        plot(ax, X, actual_values, 'k-.', 'LineWidth', 1.4, ...
            'DisplayName', 'Actual diagonal: specified yield');
    end
    styles = {'-', '--'};
    for target = 1:numel(targets)
        yline(ax, targets(target), styles{target}, 'Color', [0.25 0.25 0.25], ...
            'LineWidth', 1.3, 'DisplayName', sprintf('FOS %g', targets(target)));
    end
    hold(ax, 'off');
    xlabel(ax, x_label);
    ylabel(ax, 'Achieved FOS (uncapped)');
    title(ax, plot_title);
    xlim(ax, [min(X) max(X)]);
    finite_values = material_values(isfinite(material_values));
    if ~isempty(actual_values)
        finite_values = [finite_values(:); actual_values(isfinite(actual_values)).'];
    end
    ylim(ax, [0 1.12*max([targets(:); finite_values(:)])]);
    grid(ax, 'on');
    legend(ax, 'Location', 'southoutside', 'NumColumns', 2, ...
        'Interpreter', 'none', 'FontSize', 9);
end

function fig = fitted_figure()
    % Normalized bounds fit the current monitor, including smaller displays.
    fig = figure('Color', 'w', 'Units', 'normalized', ...
        'Position', [0.04 0.10 0.90 0.78]);
    set(fig, 'DefaultAxesFontSize', 10, 'DefaultTextFontSize', 10);
end

function geometry_surface(ax, X, Y, Z, x_label, y_label, z_label, plot_title)
    surf(ax, X, Y, Z, 'EdgeColor', 'none', 'FaceColor', 'interp', ...
        'DisplayName', 'Computed response', 'HandleVisibility', 'off');
    xlabel(ax, x_label); ylabel(ax, y_label); zlabel(ax, z_label);
    title(ax, plot_title); grid(ax, 'on'); axis(ax, 'tight');
    colormap(ax, parula(256));
    cb = colorbar(ax); cb.Label.String = z_label;
    view(ax, 45, 28);
end

function horizontal_plane(ax, X, Y, value, label)
    [PX, PY] = meshgrid([min(X(:)) max(X(:))], [min(Y(:)) max(Y(:))]);
    hold(ax, 'on');
    surf(ax, PX, PY, value*ones(2), 'FaceColor', [0.15 0.30 0.75], ...
        'FaceAlpha', 0.15, 'EdgeColor', [0.15 0.30 0.75], ...
        'DisplayName', label);
    hold(ax, 'off');
    limits = zlim(ax);
    zlim(ax, [min([0 limits(1) value]) max(limits(2), 1.08*value)]);
    legend(ax, 'show', 'Location', 'southoutside', 'FontSize', 9);
end

function mark_final_3d(ax, xvalue, yvalue, zvalue, label_text)
    if nargin < 5, label_text = sprintf('Final: %.4g', zvalue); end
    if ~all(isfinite([xvalue yvalue zvalue])), return; end
    hold(ax, 'on');
    plot3(ax, xvalue, yvalue, zvalue, 'rx', 'MarkerSize', 12, ...
        'LineWidth', 2.5, 'HandleVisibility', 'off');
    text(ax, xvalue, yvalue, zvalue, label_text, 'Color', [0.80 0 0], ...
        'FontWeight', 'bold', 'FontSize', 9, 'Clipping', 'off', ...
        'Interpreter', 'none', 'VerticalAlignment', 'bottom', ...
        'HorizontalAlignment', 'right', 'BackgroundColor', 'w', 'Margin', 2);
    hold(ax, 'off');
end

function mark_final_2d(ax, xvalue, yvalue, label_text)
    if nargin < 4, label_text = sprintf('Final: %.4g', yvalue); end
    if ~all(isfinite([xvalue yvalue])), return; end
    hold(ax, 'on');
    plot(ax, xvalue, yvalue, 'rx', 'MarkerSize', 11, ...
        'LineWidth', 2.5, 'HandleVisibility', 'off');
    text(ax, xvalue, yvalue, label_text, 'Color', [0.80 0 0], ...
        'FontWeight', 'bold', 'FontSize', 9, 'Clipping', 'off', ...
        'Interpreter', 'none', 'VerticalAlignment', 'bottom', ...
        'HorizontalAlignment', 'right', 'BackgroundColor', 'w', 'Margin', 2);
    hold(ax, 'off');
    limits = xlim(ax);
    if xvalue < limits(1) || xvalue > limits(2)
        pad = 0.03*max(diff(limits), 1);
        xlim(ax, [min(limits(1),xvalue-pad), max(limits(2),xvalue+pad)]);
    end
    limits = ylim(ax);
    if yvalue < limits(1) || yvalue > limits(2)
        pad = 0.03*max(diff(limits), 1);
        ylim(ax, [min(limits(1),yvalue-pad), max(limits(2),yvalue+pad)]);
    end
end

function deflection_surface(ax, X, Y, Z, limit_mm)
    % Split each surface triangle at the limit. This keeps ALL above-limit
    % portions gray, including cells that cross the 0.7 mm boundary.
    passing_vertices = zeros(0,3); passing_faces = zeros(0,3);
    failing_vertices = zeros(0,3); failing_faces = zeros(0,3);
    for row = 1:size(Z,1)-1
        for column = 1:size(Z,2)-1
            corners = [X(row,column), Y(row,column), Z(row,column); ...
                X(row+1,column), Y(row+1,column), Z(row+1,column); ...
                X(row+1,column+1), Y(row+1,column+1), Z(row+1,column+1); ...
                X(row,column+1), Y(row,column+1), Z(row,column+1)];
            triangles = [1 2 3; 1 3 4];
            for k = 1:2
                triangle = corners(triangles(k,:),:);
                if any(~isfinite(triangle(:))), continue; end
                passing = clip_height(triangle, limit_mm, true);
                failing = clip_height(triangle, limit_mm, false);
                [passing_vertices, passing_faces] = ...
                    append_polygon(passing_vertices, passing_faces, passing);
                [failing_vertices, failing_faces] = ...
                    append_polygon(failing_vertices, failing_faces, failing);
            end
        end
    end
    hold(ax, 'on');
    if ~isempty(passing_faces)
        patch(ax, 'Vertices', passing_vertices, 'Faces', passing_faces, ...
            'FaceVertexCData', passing_vertices(:,3), 'FaceColor', 'interp', ...
            'CDataMapping', 'scaled', 'EdgeColor', 'none', ...
            'DisplayName', sprintf('Within %.1f mm: blue to red', limit_mm));
    end
    if ~isempty(failing_faces)
        patch(ax, 'Vertices', failing_vertices, 'Faces', failing_faces, ...
            'FaceColor', [0.65 0.65 0.65], 'EdgeColor', 'none', ...
            'DisplayName', sprintf('Above %.1f mm: gray', limit_mm));
    end
    [PX, PY] = meshgrid([min(X(:)) max(X(:))], [min(Y(:)) max(Y(:))]);
    surf(ax, PX, PY, limit_mm*ones(2), 'FaceColor', [0.2 0.2 0.2], ...
        'FaceAlpha', 0.06, 'EdgeColor', [0.15 0.15 0.15], ...
        'LineWidth', 1.2, 'DisplayName', sprintf('Allowable = %.1f mm', limit_mm));
    finite_Z = Z(isfinite(Z));
    if ~isempty(finite_Z) && min(finite_Z) < limit_mm && max(finite_Z) > limit_mm
        contour3(ax, X, Y, Z, [limit_mm limit_mm], 'k-', ...
            'LineWidth', 1.7, 'HandleVisibility', 'off');
    end
    hold(ax, 'off');
    colormap(ax, jet(256));
    caxis(ax, [0 limit_mm]);
    cb = colorbar(ax); cb.Label.String = 'Within-limit displacement (mm)';
    xlabel(ax, 'Support outer b (mm)'); ylabel(ax, 'Support outer h (mm)');
    zlabel(ax, 'Total CG displacement (mm)');
    title(ax, sprintf('Total CG deflection | allowable %.1f mm', limit_mm));
    grid(ax, 'on'); axis(ax, 'tight'); view(ax, 45, 28);
    legend(ax, 'show', 'Location', 'southoutside', 'FontSize', 9);
end

function polygon = clip_height(triangle, limit_mm, keep_below)
    % Clip a triangle to one side of z = limit_mm, preserving the boundary.
    polygon = zeros(0,3);
    previous = triangle(end,:);
    previous_inside = previous(3) <= limit_mm;
    if ~keep_below, previous_inside = ~previous_inside; end
    for k = 1:size(triangle,1)
        current = triangle(k,:);
        current_inside = current(3) <= limit_mm;
        if ~keep_below, current_inside = ~current_inside; end
        if current_inside ~= previous_inside
            fraction = (limit_mm-previous(3))/(current(3)-previous(3));
            intersection = previous + fraction*(current-previous);
            intersection(3) = limit_mm;
            polygon(end+1,:) = intersection; %#ok<AGROW>
        end
        if current_inside
            polygon(end+1,:) = current; %#ok<AGROW>
        end
        previous = current;
        previous_inside = current_inside;
    end
end

function [vertices, faces] = append_polygon(vertices, faces, polygon)
    if size(polygon,1) < 3, return; end
    offset = size(vertices,1);
    vertices = [vertices; polygon];
    for k = 2:size(polygon,1)-1
        faces(end+1,:) = offset + [1 k k+1]; %#ok<AGROW>
    end
end

%% FINAL DESIGN INPUTS -- EDIT ONLY THIS BOTTOM BLOCK TO MOVE RED MARKERS
function picked = final_design_inputs()
    % These selections also set the fixed inputs in all other sweeps.
    % b_support maps to w; h_support is the outer tube height.
    % The two diagonal dimensions describe the solid axial diagonal.
    b_support_mm = 101.6;
    h_support_mm = 152.4;
    th_support_mm = 12.7;      % Matches the uploaded script: 6.35 mm x 2
    b_diagonal_mm = 12.7/2;
    h_diagonal_mm = 50.8;    % Original 25.4 fails 2d center-to-edge for a 10.1 mm hole
    L_d_m = 0.800;
    d_pin_mm = 10;

    picked = struct('w', b_support_mm/1000, ...
        'h_support', h_support_mm/1000, 'th_support', th_support_mm/1000, ...
        'w_diag', b_diagonal_mm/1000, 'h_diag', h_diagonal_mm/1000, ...
        'L_d', L_d_m, 'd_pin', d_pin_mm);
end
