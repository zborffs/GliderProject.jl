#=
p2_analysis.jl:
- Julia version: 
- Author: zachbortoff
- Date: 2022-10-18
=#
using Unitful
using IntervalArithmetic

# Define System Parameters
g = 9.81u"m/s^2";
ρ_air = 0.0752u"lb/ft^3";

# Balsa wood densities (https://www.modelaviation.com/balsa)
ρ_balsa_min = 4u"lb/ft^3"; # ultra-light density of balsa wood (lb/ft^3)
ρ_balsa_max = 14u"lb/ft^3"; # hard density of balsa wood (lb/ft^3)

# Balsa wood wing dimensions
wing_span = 36u"inch";
wing_thickness = (1/8)u"inch";
wing_chord = 2u"inch";

# Relevant balsa wood wing equations
wing_surface_area = wing_span * wing_chord;
wing_ar = uconvert(NoUnits, wing_span / wing_chord); # dimensionless
wing_volume = wing_span * wing_thickness * wing_chord;
wing_mass_min = wing_volume * ρ_balsa_min; # min estimated wing mass in (lb)
wing_mass_max = wing_volume * ρ_balsa_max; # max estimated wing mass in (lb)
weight_min = wing_mass_min * g; # min estimated weight of wing in (ft * lb / s^2)
weight_max = wing_mass_max * g; # max estimated weight of wing in (ft * lb / s^2)

# Lift must equal weight for level flight.
L_min = weight_min; # assuming L = W
L_max = weight_max; # assuming L = W

# Free airstream velocity must be in this range to achieve L = W
# L = CL * V^2 * rho / 2 * S
# S = wing surface area
# V = free-stream velocity
# rho = air density
V = 10u"mi/hr"; # (estimated launch speed... roughly running pace of human)
CL_min = uconvert(NoUnits, L_min / (1/2 * ρ_air * V^2 * wing_surface_area));
CL_max = uconvert(NoUnits, L_max / (1/2 * ρ_air * V^2 * wing_surface_area));

# We need to be able to compute the following:
# 1. aerodynamic center of the wing
# 2. center of gravity of the whole glider
# 3. aerodynmic center of the horizontal tail
# 4. coefficient of lift of the wing
# 5. coefficient of lift of the horizontal tail
# 6. change in coefficient of moment of the wing w.r.t. angle-of-attack
# 7. change in coefficient of moment of the fuselage w.r.t. angle-of-attack
# 8. change in coefficient of moment of the horizontal tail w.r.t. angle-of-attack


# We get to choose:
wing_incidence_angle = 0u"rad";
# 1. dimensions of the fuselage
fuselage_length = 36u"inch";
fuselage_width = (1/4)u"inch";
fuselage_depth = (1/4)u"inch";
# 2. indirectly the mass of the fuselage
# 3. dimensions of the horizontal tail
horizontal_tail_length =
horizontal_tail_ar = 3.97; # tail_span / tail_chord
horizontal_tail_surface_area = 38.72u"cm^2"; # (tail_span / 1.16) * (tail_chord / 3)
horizontal_tail_incidence_angle = -2u"°";

# (area of tail) * (dist. btw aero center of tail and c.g.) / [(area of wing) x (mean chord)]
# horizontal_tail_length = 36u"inch";
# horizontal_tail_width = (1/8)u"inch";
# horizontal_tail_depth = (1/8)u"inch";
# 4. indirectly the mass of the horizontal tail


# Idea 1: directly compute all of these
# Idea 2: choose arbitrarily the dimensions -> masses -> coefficients and moments
#           --> choose the one with the best lift/drag ratio and with static stability
wing_ratio = uconvert(NoUnits, wing_span / 297u"mm")
wing_chord_ratio = uconvert(NoUnits, wing_chord / 4.37u"cm")