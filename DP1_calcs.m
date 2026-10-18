%% Calculations
% Full matrix solution to get the required forces.
% Part 1: Calculate loading on pins, size according to worst-case loading.
% This keeps the original calculation order and working notes.
% Geometry is in m, except pin/hole diameters and clearance are in mm.
% Loads are in N; E is in Pa; specified yield strengths are in MPa.

%% Input variables
% Units / constants
pi_value = 3.1415;             % Original numerical approximation retained
g = 9.81;
inch_to_m = 25.4/1000;
psi_to_MPa = 0.006894757293;
FOS = 2;                      % Required FOS, given by assignment. DO NOT CHANGE.
FOS_contact_target = 3;       % Additional conservative contact target
verbose = true;

% Geometry
L_d = 800/1000;                % Horizontal reach to diagonal attachment, m
L_y = 0.5;                    % Vertical separation of diagonal ends, m
d_pin = 10;                   % Nominal pin diameter, mm
L_sh = 1.2;                   % Shaft location from support origin, m
L = 1.3;                      % Support beam length, m
shaft_drag_lever = 0.5;        % Drag-force height before beam half-height, m
h_support = 152.4/1000;           % Support beam section height, m
w = 8*inch_to_m;               % Support beam section width, m
th_support = 6.35/1000*2;        % Support beam thickness
h_diag = 50.8/1000;              % Diagonal rectangular-section height, m
w_diag = 12.7/1000;              % Diagonal rectangular-section width, m
pulley_large_in = 9.554;
pulley_small_in = 3.038;
pulley_center_in = 39.49;
ratio = 60/19;
h_cg = 350/1000;

% Standalone geometry defaults above are hardcoded to the current final
% design in pin_stress_sweep.m. The two input blocks stay independently editable.

% Plotting-only geometry overrides, before dependent geometry/contact inputs.
% The default calculation and all existing formulas remain unchanged.
if exist('sweep_override', 'var') && isstruct(sweep_override)
    if isfield(sweep_override, 'w'), w = sweep_override.w; end
    if isfield(sweep_override, 'h_support'), h_support = sweep_override.h_support; end
    if isfield(sweep_override, 'th_support'), th_support = sweep_override.th_support; end
    if isfield(sweep_override, 'w_diag'), w_diag = sweep_override.w_diag; end
    if isfield(sweep_override, 'h_diag'), h_diag = sweep_override.h_diag; end
end


% Applied loading and the original fluid/mechanical inputs
fluid_density = 1000;
fluid_scale = 0.5;
fluid_speed = 2.5;
mechanical_efficiency = 0.3;
rotor_rpm = 285;
% T_gen = 0;                    
% F_tr = 0;                   
F_d = 560;
W_sh = 150;
W_gen = 120;
% F_d = 0;
% W_sh = 0;
% W_gen = 0;

% Material inputs: keep each physical part independent.
% The supplied 36 ksi yield belongs to the BRACKETS, not the diagonals.
% NaN means a yield strength is not specified; its actual FOS stays unrated.
reference_steel_Sy_MPa = 413;  % User's assumed carbon-steel comparison
reference_6061_Sy_MPa = 36e3*psi_to_MPa; % User's supplied 6061 bracket reference
reference_steel_E_Pa = 482e9;
reference_6061_E_Pa = 70e9;
density_sup = 2700;           % kg/m^3, assumed aluminum support beam
G_support = 26 * 10^9;         % Pa, shear modulus
E_support = 70e9;             % Pa
v_support = 0.33;
Sy_support_MPa = reference_6061_Sy_MPa;         % Set from the support beam's grade/temper
support_material_name = 'Aluminum support: grade/temper not specified';
density_diag = 2700;          % kg/m^3, assumed aluminum diagonal
E_diag = 70e9;                % Pa; independent of bracket properties
v_diag = 0.33;
diag_material_name = 'Aluminum diagonal: grade/temper not specified';
Sy_diag_MPa = reference_6061_Sy_MPa;            % 6061 T6
E_bracket = reference_6061_E_Pa; % Pa; edit for the selected bracket material
v_bracket = 0.33;
Sy_bracket_MPa = reference_6061_Sy_MPa; % Supplied BRACKET yield: 36 ksi
bracket_material_name = 'Aluminum bracket (supplied 36 ksi)';
pin_E = reference_steel_E_Pa;  % Pa, assumed carbon steel
pin_v = 0.29;
Sy_pin_MPa = reference_steel_Sy_MPa; % Assumed pin yield; confirm grade/condition

% Plot reference materials. The 6061 comparison uses the supplied bracket
% rating, not a generic temper-dependent catalog yield for the diagonal.
reference_material_names = {'Carbon steel (assumed 413 MPa)', ...
                            '6061 aluminum (36 ksi bracket reference)'};
reference_material_short = {'Carbon steel', '6061 aluminum'};
reference_Sy_MPa = [reference_steel_Sy_MPa reference_6061_Sy_MPa];
reference_E_Pa = [reference_steel_E_Pa reference_6061_E_Pa];
reference_colors = [0.82 0.30 0.10; 0.06 0.39 0.70];

% Direct shear: preserve the original single-shear model by default.
% Set an entry to 2 only for a confirmed symmetric double-shear pin.
pin_shear_planes = ones(1, 6); % Order: R1_pin, R2, D1, D2, R3, R4

% Hertz contact inputs. Each entry represents one specific hole contact.
% Both connected parts are checked at each joint: R1/R2 connect the support
% beam to brackets; D1/D2 connect support beam to diagonals; R3/R4 connect
% diagonals to brackets. Review this explicit topology against your assembly.
% Add repeated pin indices to check additional plates at the same joint.
% Example outer lug: force fraction = 0.5, length = that lug's thickness.
% The contact length is the loaded plate thickness along the pin axis.
hole_tol = 0.1;               % DIAMETRAL clearance, mm; hole minus pin diameter
hertz_shear_factor = 0.30;    % Frictionless Hertz line contact: tau_max / p_0
contact_pin_index = [1 2 1 2 3 4 3 4 5 6 5 6];
contact_hole_part = {'bracket', 'bracket', 'support', 'support', ...
    'support', 'support', 'diagonal', 'diagonal', ...
    'diagonal', 'diagonal', 'bracket', 'bracket'};
contact_case_names = {'R1 bracket', 'R2 bracket', 'R1 support', 'R2 support', ...
    'D1 support', 'D2 support', 'D1 diagonal', 'D2 diagonal', ...
    'R3 diagonal', 'R4 diagonal', 'R3 bracket', 'R4 bracket'};
c_length_single_m = max([w_diag, th_support]);
contact_length_m = c_length_single_m*ones(size(contact_pin_index)); % Provisional 10 mm
contact_force_fraction = ones(size(contact_pin_index)); % Full reaction per part
% If a bracket has two equal outer lugs, replace its case with two half-load
% cases and use each lug's actual thickness for its contact length.
hole_E_cases = [E_bracket E_bracket E_support E_support E_support E_support ...
    E_diag E_diag E_diag E_diag E_bracket E_bracket];
hole_v_cases = [v_bracket v_bracket v_support v_support v_support v_support ...
    v_diag v_diag v_diag v_diag v_bracket v_bracket];
hole_Sy_MPa_cases = [Sy_bracket_MPa Sy_bracket_MPa Sy_support_MPa Sy_support_MPa ...
    Sy_support_MPa Sy_support_MPa Sy_diag_MPa Sy_diag_MPa Sy_diag_MPa Sy_diag_MPa ...
    Sy_bracket_MPa Sy_bracket_MPa];
contact_material_names = {bracket_material_name, bracket_material_name, ...
    support_material_name, support_material_name, support_material_name, support_material_name, ...
    diag_material_name, diag_material_name, diag_material_name, diag_material_name, ...
    bracket_material_name, bracket_material_name};

% Separate optional compression-buckling diagnostic for the diagonals.
K_diag = 1;                   % Effective-length factor: ideal pin-ended member

% The plotting file can override ONLY the swept variables and console output.
% Ordinary runs use the input values above. No shared workspace is required.
if exist('sweep_override', 'var') && isstruct(sweep_override)
    if isfield(sweep_override, 'L_d'), L_d = sweep_override.L_d; end
    if isfield(sweep_override, 'd_pin'), d_pin = sweep_override.d_pin; end
    if isfield(sweep_override, 'verbose'), verbose = sweep_override.verbose; end
end

assert(isscalar(L_d) && L_d > 0 && isscalar(d_pin) && d_pin > 0, ...
    'L_d and d_pin must be positive scalars.');
assert(hole_tol > 0, 'The Hertz approximation requires positive clearance.');
assert(all(pin_shear_planes >= 1) && numel(pin_shear_planes) == 6, ...
    'Specify a shear-plane count for each of the six pins.');
n_contacts = numel(contact_pin_index);
assert(all(contact_pin_index >= 1 & contact_pin_index <= 6), ...
    'Contact pin indices must be between 1 and 6.');
assert(all([numel(contact_length_m), numel(contact_force_fraction), ...
    numel(hole_E_cases), numel(hole_v_cases), numel(hole_Sy_MPa_cases), ...
    numel(contact_case_names), numel(contact_hole_part)] == n_contacts), ...
    'Every contact-case input must have the same number of entries.');
assert(all(contact_length_m > 0) && all(contact_force_fraction > 0) ...
    && all(hole_E_cases > 0), ...
    'Contact lengths, load fractions, and moduli must be positive.');
assert(all(isnan(hole_Sy_MPa_cases) | ...
    (isfinite(hole_Sy_MPa_cases) & hole_Sy_MPa_cases > 0)), ...
    'Each hole yield must be positive and finite, or NaN if unspecified.');
assert(isnan(Sy_diag_MPa) || (isfinite(Sy_diag_MPa) && Sy_diag_MPa > 0));
assert(all(reference_Sy_MPa > 0) && all(reference_E_Pa > 0));

%% Constant / dependent variables
P_fluid = 1/2*fluid_density*(fluid_scale^2)*fluid_speed^3;
% Original fluid-power expression retained; it does not drive the current
% equilibrium loads while T_gen and F_tr are explicitly set to zero.
P_mech = P_fluid*mechanical_efficiency;
w_rot = rotor_rpm*pi_value/30;  % rad/s
T_sh = P_mech/w_rot;            % N*m
T_gen = T_sh/ratio;
F_tr = T_gen/((pulley_small_in/2)*inch_to_m);

% OVERRIDE HERE
% T_gen = 0;
% F_tr = 0;
% density_diag = 0;
% density_sup = 0;

theta = atan((pulley_large_in/2 - pulley_small_in/2)/pulley_center_in);
L_g = L_sh - pulley_center_in*inch_to_m;
alpha = atan(L_y/L_d);
a = sqrt(L_d^2 + L_y^2);
L_diag = a;                    % Actual member length; L_d is horizontal reach

T_x = F_d*(shaft_drag_lever + h_support/2);
% W_sup = density_sup*w*h_support*L*g;
W_sup = density_sup*(w*h_support - (w-2*th_support)*(h_support-2*th_support))*L*g;
W_diag = density_diag*w_diag*h_diag*a*g;
% W_sup = 0;
% W_diag = 0;

% Solve for Ax + B = 0, with x composed of:
% [R1x, R2x, R1y, R2y, D1x, D2x, D1y, D2y, R1z]
A = [0 w 0 0 0 w 0 0 0];                         % M_y
A = [A; 0 0 0 0 0 0 -L_d -L_d 0];               % M_z
A = [A; 0 0 -w/2 w/2 0 0 -w/2 w/2 0];           % M_x
A = [A; 1 1 0 0 1 1 0 0 0];                     % F_x
A = [A; 0 0 1 1 0 0 1 1 0];                     % F_y
A = [A; 0 0 0 0 0 0 0 0 -1];                    % F_z
A = [A; 0 0 0 0 0 sin(alpha) 0 -cos(alpha) 0];    % Left diagonal
A = [A; 0 0 0 0 sin(alpha) 0 -cos(alpha) 0 0];     % Right diagonal
A = [A; 0 0 1 -1 0 0 0 0 0];                    % Torsion BC assumption

% Define constant matrix
% Opposing F_tr*cos(theta)*w/2 terms cancel in the original first row.
B = [-T_gen - F_tr*sin(theta)*(L_sh - L_g) - F_d*L_sh];
B = [B; W_sh*L_sh + W_sup*(L/2) + W_gen*L_g];     % M_z
B = [B; -T_x];                                  % M_x
B = [B; 0];                                     % F_x
B = [B; -W_sup - W_sh - W_gen];                  % F_y
B = [B; F_d];                                   % F_z
B = [B; -W_diag*cos(alpha)/2];                   % Left diagonal
B = [B; -W_diag*cos(alpha)/2];                   % Right diagonal
B = [B; 0];

% Solve for x using matrix: Ax = -B
format shortG
x = -(A\B);
equilibrium_residual = norm(A*x + B, inf);

R1_pin = sqrt(x(1)^2 + x(3)^2);
R2 = sqrt(x(2)^2 + x(4)^2);
D1 = sqrt(x(5)^2 + x(7)^2);
D2 = sqrt(x(6)^2 + x(8)^2);
magnitude = [R1_pin R2 D1 D2];

% This is only for pins in the support beam; separate equations are needed
% for pins in the diagonals. Size to the worst case as well.
% R3 ('front' diagonal); R4 ('rear' diagonal)
R3_y = x(7) + W_diag;
R4_y = x(8) + W_diag;
R3_x = x(5);
R4_x = x(6);
R3 = sqrt(R3_y^2 + R3_x^2);
R4 = sqrt(R4_y^2 + R4_x^2);

pin_names = {'R1_pin', 'R2', 'D1', 'D2', 'R3', 'R4'};
pin_reactions_arr = [R1_pin, R2, D1, D2, R3, R4];
[max_pin_reaction_mag, max_pin_index] = max(pin_reactions_arr);

% Tresca criterion: its added conservatism is helpful for shear.
pin_shear_yield_MPa = Sy_pin_MPa/2;
pin_shear_max_FOS = pin_shear_yield_MPa/FOS;       % Required allowable, MPa
pin_shear_stress_arr = pin_reactions_arr ...
    ./ ((pi_value*d_pin^2/4)*pin_shear_planes);    % Nominal shear, MPa
[max_pin_reaction_stress, max_shear_pin_index] = max(pin_shear_stress_arr);
pin_shear_FOS_arr = pin_shear_yield_MPa./pin_shear_stress_arr;
FOS_pin_shear = min(pin_shear_FOS_arr);            % Achieved FOS, not allowable

%% Hertzian contact stress
% Assuming a carbon-steel pin and aluminum hole as specified in inputs.
% Evaluate every contact case; the largest pin reaction alone may not govern
% contact stress when plate thickness, material, or load split differs.
% This replaces the original always-true (~=5 || ~=6) material branch.
d_hole = d_pin + hole_tol;                        % mm
d_pin_m = d_pin/1000;                            % Hertz calculations use SI
d_hole_m = d_hole/1000;
curvature_difference = 1/d_pin_m - 1/d_hole_m;    % Convex pin, concave hole
F_contact_cases = pin_reactions_arr(contact_pin_index).*contact_force_fraction;
contact_compliance = (1-hole_v_cases.^2)./hole_E_cases ...
    + (1-pin_v^2)/pin_E;
K_b_cases = sqrt((2./(pi_value*contact_length_m)) ...
    .* contact_compliance/curvature_difference);
b_cases = K_b_cases.*sqrt(F_contact_cases);       % Contact half-width, m
p_max_cases = 2*F_contact_cases ...
    ./ (pi_value*b_cases.*contact_length_m);     % Peak contact pressure, Pa
p_max_cases_MPa = p_max_cases/1e6;
tau_contact_cases = hertz_shear_factor*p_max_cases; % Hertz line-contact shear, Pa
tau_contact_cases_MPa = tau_contact_cases/1e6;

% Compare maximum subsurface shear with Tresca shear yield in BOTH bodies.
FOS_contact_pin_cases = (Sy_pin_MPa/2)./tau_contact_cases_MPa;
FOS_contact_hole_cases = (hole_Sy_MPa_cases/2)./tau_contact_cases_MPa;
FOS_contact_cases = min(FOS_contact_pin_cases, FOS_contact_hole_cases);
contact_strength_known = isfinite(hole_Sy_MPa_cases);
FOS_contact_cases(~contact_strength_known) = NaN;
% Never allow MATLAB's missing-value behavior to report a complete joint FOS
% when any loaded hole material has an unspecified yield strength.
FOS_contact = NaN;
contact_FOS_case_index = NaN;
if all(contact_strength_known)
    [FOS_contact, contact_FOS_case_index] = min(FOS_contact_cases);
end
FOS_contact_pin = min(FOS_contact_pin_cases);
FOS_contact_brackets = min(FOS_contact_cases(strcmp(contact_hole_part, 'bracket')));
[p_max, max_contact_case_index] = max(p_max_cases);
tau_contact = max(tau_contact_cases);
p_max_MPa = p_max/1e6;
tau_contact_MPa = tau_contact/1e6;
% Convenient scalar aliases for the highest-pressure contact case.
b = b_cases(max_contact_case_index);
K_b = K_b_cases(max_contact_case_index);
l = contact_length_m(max_contact_case_index);
hole_E = hole_E_cases(max_contact_case_index);
hole_v = hole_v_cases(max_contact_case_index);

% Diagnostic only: Hertz assumes contact dimensions are small relative to
% body radii. A large ratio means a conformal-contact model may be needed.
hertz_halfwidth_over_pin_radius = b_cases/(d_pin_m/2);

%% Diagonal Sizing
% Plan: decouple into bending and axial forces.
% Weight becomes a distributed transverse load. Its axial component changes
% the axial force between the ends, so record both ends for each diagonal.
% Size and shape are CONSTANT for this analysis: solid rectangular 10 x 10 mm.
A_diag = w_diag*h_diag;

%% Diagonal Axial Loading / Deflection
W_diag_axial = W_diag*sin(alpha);
F_axial_support_1 = x(5)*cos(alpha) + x(7)*sin(alpha);
F_axial_support_2 = x(6)*cos(alpha) + x(8)*sin(alpha);
F_axial_1 = W_diag_axial + F_axial_support_1;
F_axial_2 = W_diag_axial + F_axial_support_2;
axial_end_forces = [F_axial_support_1 F_axial_1; ...
                    F_axial_support_2 F_axial_2];
max_axial_by_diag = max(abs(axial_end_forces), [], 2).';
[max_axial, max_ax_index] = max(max_axial_by_diag);
max_ax_stress = max_axial/A_diag;                 % Maximum magnitude, Pa
axial_stress_by_diag = max_axial_by_diag/A_diag;
delta_ax_diag = max_ax_stress*a/E_diag*1000;     % mm

%% Diagonal Bending Loading / Deflection
q_diag = W_diag*cos(alpha)/a;
I_diag = (1/12)*w_diag*h_diag^3;

% Diagonal bending and bending stress
diag_1_int_force = x(5)*sin(alpha) - x(7)*cos(alpha);
diag_1_bracket_force = -R3_x*sin(alpha) + R3_y*cos(alpha);
diag_2_int_force = x(6)*sin(alpha) - x(8)*cos(alpha);
diag_2_bracket_force = -R4_x*sin(alpha) + R4_y*cos(alpha);
% These transverse end forces are equal: both diagonals have the same
% self-weight bending; their axial loads differ.
M_max = q_diag*(a/4)*(a-a/2);                     % q*a^2/8, pinned uniform load
sigma_max_bend = M_max*(h_diag/2)/I_diag;
% delta_diag_bend = (5/384)*q_diag*a^4/(E_diag*I_diag)*1000; % mm

%% Combined Stress, Diagonal
% Only bending and axial normal stresses are combined here.
% Adding each member's worst axial magnitude and peak bending gives a
% conservative gross-section envelope; the peaks need not coincide spatially.
sigma_max_diag_by_member = sigma_max_bend + axial_stress_by_diag;
sigma_max_diag = max(sigma_max_diag_by_member);
sigma_max_diag_MPA = sigma_max_diag/1e6;
FOS_diag_yield_by_member = Sy_diag_MPa./(sigma_max_diag_by_member/1e6);
FOS_diag_yield = NaN;
if isfinite(Sy_diag_MPa), FOS_diag_yield = min(FOS_diag_yield_by_member); end

% Separate Euler buckling diagnostic: high yield FOS does not check buckling.
% Positive axial force is taken as tension; negative is compression.
compression_by_diag = max(0, -min(axial_end_forces, [], 2)).';
I_diag_weak = min(I_diag, h_diag*w_diag^3/12);
P_cr_diag = pi_value^2*E_diag*I_diag_weak/(K_diag*a)^2;
FOS_diag_buckling_by_member = inf(1, 2);
is_compressive = compression_by_diag > 0;
FOS_diag_buckling_by_member(is_compressive) = ...
    P_cr_diag./compression_by_diag(is_compressive);
FOS_diag_buckling = min(FOS_diag_buckling_by_member);

%% Support Beam Calculations

%% Rotor CG deflection
%% ONLY TORSION
x_locs = 0:10/1000:L;
% for vals, CCW about global x axis is positive, because this is positive z
T_x_locs = [L_g , L_d, L_sh, L_sh];
M_x_chain = F_tr*sin(theta)*(50/1000 + h_support/2);
M_Dx = x(7)*w/2 - x(8)*w/2;
T_x_vals = [-M_x_chain, -M_Dx, M_x_chain, T_x];

[alpha_sh_torsion, ~, M_x_int, ~] = beamTorsion(x_locs, L_sh, T_x_locs, T_x_vals, w, h_support, G_support, th_support);
alpha_sh_torsion
M_x_int;

def_torsion = alpha_sh_torsion * 350  % mm
pin_def_torsion = alpha_sh_torsion * w/2 * 1000 %mm
%% Deflection due to MY
% [v_sh, M] = beamResponse(x, L_sh, L, EI, forcePos, forceVal, momentPos, momentVal, q_sup)
% MY deflects in z direction
% MZ deflects in y direction
% I_support_MY = (1/12)*(h_support*w^3);
I_support_MY = (1/12)*(h_support*w^3 - (h_support-2*th_support)*(w-2*th_support)^3);
My_forces_beam = [F_d];
My_forces_loc = [L_sh];
M_Dy = (x(5)-x(6))*(w/2);
My_moment_loc = [L_sh, L_d];
My_moments_beam = [T_sh, M_Dy];
[def_sh_z, M_y_x] = beamResponse(x_locs, L_sh, L, E_support*I_support_MY, ...
    My_forces_loc, My_forces_beam, My_moment_loc, My_moments_beam, 0);
def_sh_z = def_sh_z * 1000 + def_torsion; % mm, torsion also acts in z direction

%% Deflection due to MZ
% I_support_MZ = (1/12)*(w*h_support^3);
I_support_MZ = (1/12)*(w*h_support^3 - (w-2*th_support)*(h_support-2*th_support)^3);
Mz_forces_beam = [-W_sh, (x(7)+x(8)), -W_gen];
Mz_forces_loc = [L_sh, L_d, L_g];
q_sup = W_sup / L;
[def_sh_y, M_z_x] = beamResponseMZ(x_locs, L_sh, L, E_support*I_support_MZ, ...
    Mz_forces_loc, Mz_forces_beam, 0, 0, -q_sup, L_d, 0);
% def_lg_y = beamResponse(x_locs)
def_sh_y = def_sh_y * 1000; % mm

%% SUPPORT BEAM Axial Deflection / STRESS
% There are only 2 axial forces, D12x and R12x, which both balance out
% sup_cs_area = w * h_support
sup_cs_area = w * h_support - (w-th_support*2)*(h_support-th_support*2);
def_x = (x(7)+x(8))*L_g / (E_support * sup_cs_area) * 1000 % mm

sup_axial_str = zeros(1, numel(x_locs));

for i = 1:length(x_locs)
    if (x_locs(i) > L_d)
        sup_axial_str(i) = (x(7)+x(8))/sup_cs_area;
    end
end

sup_axial_str

%% TOTAL DEFLECTION OF CG
% def_vector = [def_sh_y, def_sh_z, def_x]
def_vector = [def_x, def_sh_y, def_sh_z]
def_magnitude = norm(def_vector)

%% STRESSES IN THE SUPPORT BEAM
% Doesn't take into account peak stresses, that will be a separate section
% TO-DO HERE
support_point_names = { ...
    'Top z+', 'Top z-', 'Bottom z+', 'Bottom z-', ...
    'Top midpoint', 'Bottom midpoint', 'z+ midpoint', 'z- midpoint'};

support_y = [h_support/2, h_support/2, -h_support/2, -h_support/2, ...
             h_support/2, -h_support/2, 0, 0];

support_z = [w/2, -w/2, w/2, -w/2, ...
             0, 0, w/2, -w/2];

%% Signed normal stresses [Pa]
support_sigma_MY = -(M_y_x(:)/I_support_MY) * support_z;
support_sigma_MZ = -(M_z_x(:)/I_support_MZ) * support_y;

support_sigma_bending = support_sigma_MY + support_sigma_MZ;

% Add axial normal stress at the SAME longitudinal location.
support_sigma = support_sigma_bending ...
    + repmat(sup_axial_str(:), 1, numel(support_point_names));

% Maximum absolute bending stress at each sampled section.
support_max_bending = max(abs(support_sigma_bending(:,1:4)), [], 2);

%% Torsional shear magnitude [Pa]
% Nominal shear from the existing closed rectangular tube model.
support_A_mid = (w-th_support)*(h_support-th_support);
support_tau_nominal = abs(M_x_int(:))/(2*support_A_mid*th_support);

% Select midpoints of the longer sides.
support_torsion_mask = zeros(1, numel(support_point_names));

if w > h_support
    support_torsion_mask(5:6) = 1; % Top and bottom sides are longer
elseif h_support > w
    support_torsion_mask(7:8) = 1; % z+ and z- sides are longer
else
    support_torsion_mask(5:8) = 1; % Square: all side midpoints
end

support_tau = support_tau_nominal * support_torsion_mask;

%% Combined stress at each point and longitudinal location
% One normal-stress component plus torsional shear.
support_vonMises = sqrt(support_sigma.^2 + 3*support_tau.^2);

% Tresca equivalent stress = largest principal-stress difference.
support_Tresca = sqrt(support_sigma.^2 + 4*support_tau.^2);

%% Maximum equivalent stresses at each sampled section
[support_VM_by_x, support_VM_point_by_x] = ...
    max(support_vonMises, [], 2);

[support_Tresca_by_x, support_Tresca_point_by_x] = ...
    max(support_Tresca, [], 2);

%% Overall maxima across the sampled sections
[support_max_VM, support_VM_x_index] = max(support_VM_by_x);
support_VM_point_index = support_VM_point_by_x(support_VM_x_index);

[support_max_Tresca, support_Tresca_x_index] = max(support_Tresca_by_x);
support_Tresca_point_index = ...
    support_Tresca_point_by_x(support_Tresca_x_index);

%% Results
fprintf('Max absolute bending stress: %.3f MPa\n', ...
    max(support_max_bending)/1e6);

fprintf('Max von Mises: %.3f MPa at x = %.4f m, %s\n', ...
    support_max_VM/1e6, x_locs(support_VM_x_index), ...
    support_point_names{support_VM_point_index});

fprintf('Max Tresca equivalent: %.3f MPa at x = %.4f m, %s\n', ...
    support_max_Tresca/1e6, x_locs(support_Tresca_x_index), ...
    support_point_names{support_Tresca_point_index});

% Achieved yield factors of safety; NaN if yield is unspecified.
FOS_support_VM = Sy_support_MPa/(support_max_VM/1e6);
FOS_support_Tresca = Sy_support_MPa/(support_max_Tresca/1e6);


%% Results for the companion plotting file
% This packaging does not replace the original workspace variables.
calc = struct();
calc.L_d = L_d;
calc.L_diag = L_diag;
calc.d_pin = d_pin;
calc.FOS_target = FOS;
calc.FOS_contact_target = FOS_contact_target;
calc.pi_value = pi_value;
calc.L_y = L_y;
calc.w_diag = w_diag;
calc.h_diag = h_diag;
calc.A_diag = A_diag;
calc.I_diag = I_diag;
calc.I_diag_weak = I_diag_weak;
calc.K_diag = K_diag;
calc.Sy_pin_MPa = Sy_pin_MPa;
calc.Sy_diag_MPa = Sy_diag_MPa;
calc.diag_material_name = diag_material_name;
calc.support_material_name = support_material_name;
calc.bracket_material_name = bracket_material_name;
calc.Sy_bracket_MPa = Sy_bracket_MPa;
calc.Sy_support_MPa = Sy_support_MPa;
calc.density_diag = density_diag;
calc.E_diag = E_diag;
calc.v_diag = v_diag;
calc.density_sup = density_sup;
calc.E_support = E_support;
calc.v_support = v_support;
calc.E_bracket = E_bracket;
calc.v_bracket = v_bracket;
calc.reference_material_names = reference_material_names;
calc.reference_material_short = reference_material_short;
calc.reference_Sy_MPa = reference_Sy_MPa;
calc.reference_E_Pa = reference_E_Pa;
calc.reference_colors = reference_colors;
calc.pin_E = pin_E;
calc.pin_v = pin_v;
calc.hole_tol = hole_tol;
calc.hertz_shear_factor = hertz_shear_factor;
calc.pin_shear_planes = pin_shear_planes;
calc.pin_names = pin_names;
calc.x = x;
calc.pin_reactions = pin_reactions_arr;
calc.pin_shear_stress_MPa = pin_shear_stress_arr;
calc.max_pin_reaction_mag = max_pin_reaction_mag;
calc.max_pin_index = max_pin_index;
calc.FOS_pin_shear = FOS_pin_shear;
calc.equilibrium_residual = equilibrium_residual;
calc.contact_pin_index = contact_pin_index;
calc.contact_case_names = contact_case_names;
calc.contact_hole_part = contact_hole_part;
calc.contact_material_names = contact_material_names;
calc.contact_strength_known = contact_strength_known;
calc.contact_length_m = contact_length_m;
calc.contact_force_fraction = contact_force_fraction;
calc.hole_E_cases = hole_E_cases;
calc.hole_v_cases = hole_v_cases;
calc.hole_Sy_MPa_cases = hole_Sy_MPa_cases;
calc.contact_compliance = contact_compliance;
calc.b_cases_m = b_cases;
calc.p_max_cases_MPa = p_max_cases_MPa;
calc.tau_contact_cases_MPa = tau_contact_cases_MPa;
calc.FOS_contact_pin_cases = FOS_contact_pin_cases;
calc.FOS_contact_hole_cases = FOS_contact_hole_cases;
calc.FOS_contact_cases = FOS_contact_cases;
calc.FOS_contact = FOS_contact;
calc.FOS_contact_pin = FOS_contact_pin;
calc.FOS_contact_brackets = FOS_contact_brackets;
calc.hertz_halfwidth_over_pin_radius = hertz_halfwidth_over_pin_radius;
calc.axial_end_forces = axial_end_forces;
calc.axial_stress_by_diag_MPa = axial_stress_by_diag/1e6;
calc.sigma_max_bend_MPa = sigma_max_bend/1e6;
calc.sigma_diag_by_member_MPa = sigma_max_diag_by_member/1e6;
calc.sigma_max_diag_MPa = sigma_max_diag_MPA;
calc.FOS_diag_yield_by_member = FOS_diag_yield_by_member;
calc.FOS_diag_yield = FOS_diag_yield;
calc.compression_by_diag = compression_by_diag;
calc.P_cr_diag = P_cr_diag;
calc.FOS_diag_buckling_by_member = FOS_diag_buckling_by_member;
calc.FOS_diag_buckling = FOS_diag_buckling;

% Additional outputs for plotting; no calculation formulas changed.
calc.w = w;
calc.h_support = h_support;
calc.th_support = th_support;
calc.L = L;
calc.L_sh = L_sh;
calc.G_support = G_support;
calc.support_max_VM_MPa = support_max_VM/1e6;
calc.support_max_Tresca_MPa = support_max_Tresca/1e6;
calc.def_vector_mm = def_vector;
calc.def_magnitude_mm = def_magnitude;
calc.alpha_sh_torsion_rad = alpha_sh_torsion;
calc.x_locs = x_locs;
calc.M_y_x = M_y_x;
calc.M_z_x = M_z_x;
calc.M_x_int = M_x_int;
calc.support_VM_by_x_MPa = support_VM_by_x/1e6;


if verbose
    disp('Reaction solution [R1x R2x R1y R2y D1x D2x D1y D2y R1z]:');
    disp(x.');
    disp('Pin reactions [R1_pin R2 D1 D2 R3 R4], N:');
    disp(pin_reactions_arr);
    fprintf('L_d = %.3f m; actual L_diag = %.3f m; d_pin = %.1f mm.\n', ...
        L_d, L_diag, d_pin);
    fprintf('Max nominal pin shear = %.3f MPa; achieved FOS = %.3f.\n', ...
        max_pin_reaction_stress, FOS_pin_shear);
    fprintf('Max Hertz pressure = %.3f MPa; contact shear = %.3f MPa.\n', ...
        p_max_MPa, tau_contact_MPa);
    fprintf('Pin Hertz FOS = %.3f; bracket-contact FOS = %.3f (target %.0f).\n', ...
        FOS_contact_pin, FOS_contact_brackets, FOS_contact_target);
    if all(contact_strength_known)
        fprintf('Complete pin/hole contact FOS = %.3f.\n', FOS_contact);
    else
        fprintf('Complete contact FOS: unrated; set the unknown hole material yields.\n');
    end
    fprintf('Max diagonal stress envelope = %.3f MPa.\n', sigma_max_diag_MPA);
    if isfinite(Sy_diag_MPa)
        fprintf('Diagonal yield FOS = %.3f using its own yield %.3f MPa.\n', ...
            FOS_diag_yield, Sy_diag_MPa);
    else
        fprintf('Diagonal yield FOS: unrated; set Sy_diag_MPa for its grade/temper.\n');
    end
    fprintf('Separate Euler buckling FOS (K = %g) = %.3f.\n', K_diag, FOS_diag_buckling);
    fprintf('Max Hertz half-width / pin radius = %.3f (must be small for Hertz).\n', ...
        max(hertz_halfwidth_over_pin_radius));
    fprintf('Cross section fixed: %.1f x %.1f mm, solid rectangle.\n', ...
        1000*w_diag, 1000*h_diag);
end
