#=
p3_analysis.jl:
- Date: 2022-10-18
=#
using Unitful

function k(b_v, r_1)
    # implements what's in Figure 3.75 functionally
    indep_var = b_v / (2 * r_1);
    if indep_var <= 2
        return 0.75;
    elseif indep_var <= 3.5
        m = (1-0.75) / (3.5 - 2);
        return m * indep_var + (0.75 - m * 2);
    else # indep_var > 4
        return 1.0;
    end
end

# Define Environmental Parameters
g = 9.81u"m/s^2";
ρ_air = 0.0752u"lb/ft^3";

# Define Aircraft Flight Dynamics constants
a0 = 2 * pi;

# Define balsa wood densities (https://www.modelaviation.com/balsa)
ρ_balsa_min = 4u"lb/ft^3"; # ultra-light density of balsa wood (lb/ft^3)
ρ_balsa_max = 14u"lb/ft^3"; # hard density of balsa wood (lb/ft^3)

# Define wing dimensions
wing_span = 36u"inch";
wing_thickness = (1/8)u"inch";
wing_chord = 2u"inch";
wing_incidence_angle = 0u"°"; # just assign wing incidence to be equal to 0

# Define fuselage dimensions (modifiable)
fuselage_length = 18u"inch"; # 9 not 12 b/c 18 / 2 = 9
fuselage_width = (1/4)u"inch"; # shouldn't really change, but hypothetically changeable
fuselage_thickness = (1/4)u"inch"; # shouldn't really change, but hypothetically changeable

# Define horizontal tail parameters (modifiable)
htail_span = 8u"inch"; # 36 inch -> 12 -> try to make htail_ar really really small so CL_tail is insignificant.
htail_thickness = (1/8)u"inch";
htail_chord = 3u"inch"; # up to 3 inch

# Define vertical tail parameters (modifiable)
vtail_span = 4u"inch"; # up to 36 - htail_span
vtail_thickness = (1/8)u"inch";
vtail_chord = 3u"inch"; # up to 3 inch

# Define tail placement w.r.t. wing (modifiable)
lt = (0.75 * fuselage_length); # lt = dist from aero center of wing to aero center of tail

# Define vtail placement w.r.t. wing (modifiable)
lv = (0.75 * fuselage_length);

# Implicitly define c.g. ("xa" is modifiable)
xa = 1; # number of chords fore of the center of gravity that the aero. center is for wing
x_cg = lt # compute the center of gravity
@assert(x_cg < fuselage_length, "System c.g. exceeds fuselage length")
@assert(htail_span + vtail_span < 36u"inch", "This design requires using more tail wood than allowed")


####################################################################################################
####################################################################################################
# At this point, we have defined all the parameters necessary to determine build our aircraft and  #
# to determine whether the aircraft is statically stable. Before building it, let's check if it's  #
# statically stable, then update the parameters above incrementally.                               #
####################################################################################################
####################################################################################################


# Compute wing parameters
wing_surface_area = wing_span * wing_chord;
wing_ar = uconvert(NoUnits, wing_span / wing_chord); # dimensionless

# Compute the coefficient of lift for a flat plate wing
CL_wing = a0 / (1 + a0 / (pi * wing_ar)); # dimensionless (function of wing span and wing chord)

# This restriction enforces Cm(alpha = 0) > 0, which is necessary condition for stability.
# Proof:
#   1. Cm(alpha) = (CL)_w * x_bar - (CL)_t Vbar eta
#   2. <=> Cm(alpha) = (a0 / (1 + a0 / (pi * wing_ar))) * (alpha + i_w) * x_bar - (a0 / (1 + a0 / (pi * tail_ar))) * (alpha + i_t - epsilon) Vbar eta
#   3. <=> Cm(alpha=0) = -(a0 / (1 + a0 / (pi * wing_ar))) * (i_t - epsilon) Vbar eta > 0
#   4. <=> Cm(alpha=0) = (i_t - epsilon) > 0
#   5. <=> Cm(alpha=0) = i_t > epsilon
epsilon = (1.62 * CL_wing / (pi * wing_ar))u"rad"; # depends on the CL_wing we chose above
tail_incidence_angle = uconvert(u"°", epsilon) - 10u"°" # adding -15 degree to create some margin -> this choice basically chooses the aoa trim

# Compute horizontal tail parameters
htail_surface_area = htail_span * htail_chord;
htail_ar = uconvert(NoUnits, htail_span / htail_chord); # dimensionless
htail_volume = htail_span * htail_thickness * htail_chord;

# Compute the coefficient of lift for a flat plate tail
CL_tail = uconvert(NoUnits, a0 / (1 + a0 / (pi * htail_ar)))

# Compute the overall aircraft partial of the coefficient of moment w.r.t. angle of attack
# assume eta_t = 0...
V1 = (htail_surface_area * lt) / (wing_surface_area * wing_chord);
Cm_alpha = CL_wing * xa - CL_tail * V1

# We can compute the trim angle of attack from these calculations
alpha = CL_tail * deg2rad(tail_incidence_angle - epsilon) * V1 / (CL_wing * xa - CL_tail * V1) # this value must be positive by construction

@assert(Cm_alpha < 0, "System is statically longintudinally unstable")
@assert(alpha > 0, "System must have positive angle-of-attack trim to be fly-able")

# Overall goals of changing xa, lt, htail_span, htail_chord, fuselage_length is:
# 1. make sure Cm_alpha is negative for stability
# 2. Make sure x_cg < fuselage_length (and not too far forward b/c we'd end up needing to weight it
#    down more in order to achieve that x_cg, which would mean more weight)
# 3. Make sure CL_tail is minimized. This is because, CL_tail operates downwards, so we lose overall
#    lift if CL_tail is super high (even though that would make it more stable) and we won't go as
#    far.
# Summary: Minimize Cm_alpha, minimize x_cg, minimize CL_tail.

# Directional Stability
# - (d C_n / d beta) must be positive to be stable
# - only affected by vertical tail (pg 276 in book)
#   - (d C_n / d beta)_tail,fixed = k * a_v * (1 + dsigma/dbeta) * eta_v * V_2
#      - V2 = Sv * lv / (S * c), surface area of vertical and distance to ac from aircraft cg
#      - (1 + dsigma/dbeta) * eta_v = 0.724 + (3.06 * Sv / S) / 2 + 0.4 * zw / d_fmax + 0.009 A
#           - df_max is max fuselage depth = (1/4) inch
#           - zw is vertical dist btw wing root quarter chord to fuselage center line, positive if wing below fuselage, i.e. (1/8) / 2 + (1/4) / 2 = 1/16 + 1/8 = 3/16 inch
#           - a_v is just the a0/(1+a0/(pi * AR_v)) of the vertical tail -> depends on AR of VTail
#           - k (pg 274) "parameter 'k' is given by Fig. 3.75 as a function of b_v / 2r_1", where b_v is the vertical tail span measured up to the fuselage centerline and r_1 is the average radius of the fuselage section underneath tail"
#                - b_v = (1/4) / 2 + vertical tail span
#                - r_1 = (for rectangular planform) vertical tail chord

# Compute horizontal tail parameters
vtail_surface_area = vtail_span * vtail_chord;
vtail_ar = uconvert(NoUnits, vtail_span / vtail_chord); # dimensionless
vtail_volume = vtail_span * vtail_thickness * vtail_chord;

# Cn_beta = k * a_v * (1 + dsigma/dbeta) * eta_ * V_2
k_variable = k(fuselage_width / 2 + vtail_span, vtail_chord);
a_v = a0 / (1 + a0 / (pi * vtail_ar));
zw = fuselage_thickness / 2 + wing_thickness / 2;
dsigma_dbeta_eta_v = 0.724 + (3.06 * vtail_surface_area / wing_surface_area) / 2 + 0.4 * zw / fuselage_thickness + 0.009 * wing_ar;
V_2 = vtail_surface_area * lv / (wing_surface_area * wing_span);
Cn_beta = k_variable * a_v * dsigma_dbeta_eta_v * V_2

@assert(Cn_beta > 0, "System is statically directionally stable")


# Lateral Stability
# stable if Cl_beta < 0
# - fuselage contribution... insignificant
# - wing contribution... insignificant
# - horizontal tail contribution... insignificant
# - vertical tail contribution... significant

# vertical tail contribution basically exactly same as directional static stability coeff, but with
# minor modification to the last variable
# If we are directionally stable, then we are guaranteed to be laterally stable if zv cos(alpha) > lv sin(alpha);
# - We must choose angle of attack to trim the vehicle at, then check that the above inequality is true...
#   - zv is the vertical distance from the a.c. of vtail to centerline of fuselage (which we will assume is the c.g. of aircraft)
#   - lv is the longitudinal distance from the a.c. of vtail to c.g. in longitudinal plane
zv = fuselage_thickness / 2 + vtail_span / 2;
# lv = (whatever we picked above)
Cl_beta = -k_variable * a_v * dsigma_dbeta_eta_v * (zv * cos(alpha) - lv * sin(alpha)) / wing_span;
@assert(zv * cos(alpha) - lv * sin(alpha) > 0u"inch" && Cn_beta > 0, "System is laterally unstable")

# Estimated Center of Gravity of Balsa Wood Glider (without added mass)

# Compute Wing mass
wing_volume = wing_span * wing_thickness * wing_chord;
wing_mass_min = wing_volume * ρ_balsa_min; # min estimated wing mass in (lb)
wing_mass_max = wing_volume * ρ_balsa_max; # max estimated wing mass in (lb)

# Compute Fuselage mass
fuselage_volume = fuselage_length * fuselage_width * fuselage_thickness;
fuselage_mass_min = fuselage_volume * ρ_balsa_min; # min estimated wing mass in (lb)
fuselage_mass_max = fuselage_volume * ρ_balsa_max; # max estimated wing mass in (lb)

# Compute Horizontal tail mass
htail_volume = htail_span * htail_chord * htail_thickness;
htail_mass_min = htail_volume * ρ_balsa_min; # min estimated wing mass in (lb)
htail_mass_max = htail_volume * ρ_balsa_max; # max estimated wing mass in (lb)

# Compute Vertical tail mass
vtail_volume = vtail_span * vtail_chord * vtail_thickness;
vtail_mass_min = vtail_volume * ρ_balsa_min; # min estimated wing mass in (lb)
vtail_mass_max = vtail_volume * ρ_balsa_max; # max estimated wing mass in (lb)

# Compute overall weight bounds (out of curiosity)
weight_min = uconvert(u"lbf", (wing_mass_min + fuselage_mass_min + htail_mass_min + vtail_mass_min) * g);
weight_max = uconvert(u"lbf", (wing_mass_max + fuselage_mass_max + htail_mass_max + vtail_mass_max) * g);

# Compute estimated c.g.
# - we assumed a.c. of all surfaces is "middle" of surfaces. Assume uniform distribution of mass in
#   surface. Therefore, c.g. of each surface will be in "middle" of surface. Therefore, a.c. of
#   surface is same as c.g. of surface.
# - we defined position of each surface in terms of it's relationship btw the surface's a.c. to the
#   desired c.g.
wing_cg = -xa * wing_chord; # wing's c.g. must be 'xa' chord lengths fore of desired c.g.
htail_cg = lt; # htail's c.g. must be 'lt' aft of the desired c.g.
vtail_cg = lv; # vtail's c.g. must be 'lv' aft of the desired c.g.

# assume we attach the vtail at the "end" of the fuselage.... geometry blah ...
fuselage_cg = lv - fuselage_length / 2

#  m1 * d1 + m2 * d2 = 0 (@ cg)
# (Find m1) m2 * d2 / (-d1) = m1
offset_weight_min = uconvert(u"lb", (wing_cg * wing_mass_min + htail_cg * htail_mass_min + vtail_cg * vtail_mass_min + fuselage_cg * fuselage_mass_min) / (fuselage_length - x_cg));
offset_weight_max = uconvert(u"lb", (wing_cg * wing_mass_max + htail_cg * htail_mass_max + vtail_cg * vtail_mass_max + fuselage_cg * fuselage_mass_max) / (fuselage_length - x_cg));

# Compute overall weight bounds (out of curiosity)
weight_min = uconvert(u"lbf", (wing_mass_min + fuselage_mass_min + htail_mass_min + vtail_mass_min + offset_weight_min) * g);
weight_max = uconvert(u"lbf", (wing_mass_max + fuselage_mass_max + htail_mass_max + vtail_mass_max + offset_weight_max) * g);

# Compute overall required throw speed for steady-flight (assuming CL_overall = CL_wing)
launch_speed_min = uconvert(u"mi/hr", 1 / sqrt(CL_wing * (1/2 * ρ_air * wing_surface_area) / weight_min))
launch_speed_max = uconvert(u"mi/hr", 1 / sqrt(CL_wing * (1/2 * ρ_air * wing_surface_area) / weight_max))

# Report Results:
using PrettyTables
print("Aircraft Parameters (we choose these): ")
aircraft_parameters_table = [
    "Wing Span" wing_span;
    "Wing Chord" wing_chord;
    "Wing Thickness" wing_thickness;
    "Wing Incidence Angle" wing_incidence_angle;
    "Fuselage Length" fuselage_length;
    "Fuselage Width" fuselage_width;
    "Fuselage Thickness" fuselage_thickness;
    "Horizontal Tail Span" htail_span;
    "Horizontal Tail Chord" htail_chord;
    "Horizontal Tail Thickness" htail_thickness;
    "Horizontal Tail Incidence Angle" tail_incidence_angle;
    "Vertical Tail Span" vtail_span;
    "Vertical Tail Chord" vtail_chord;
    "Vertical Tail Thickness" vtail_thickness;
    "Wing Tail Placement (xa in chord lengths)" xa;
    "Horizontal Tail Placement (lt)" lt;
    "Vertical Tail Placement (lv)" lv;
    "Desired CG from Backend of Fuselage" x_cg;
];

pretty_table(aircraft_parameters_table; header=["Parameter Name", "Parameter Value"], header_crayon=crayon"green bold")


print("Stability Parameters (we compute these):")
stability_parameters_table = [
    "Trim alpha" alpha "-";
    "Cm_alpha" Cm_alpha true;
    "Cn_beta" Cn_beta true;
    "Cl_beta" "-" true;
    "Min. Mass Offset @ Nose" offset_weight_min "-";
    "Max. Mass Offset @ Nose" offset_weight_max "-";
];
pretty_table(stability_parameters_table; header=["Parameter Name", "Parameter Value", "Is Statically Stable"], header_crayon=crayon"green bold")