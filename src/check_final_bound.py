"""Independent high-precision arithmetic check for the 67.32059784% deduction."""
from fractions import Fraction
from mpmath import iv

iv.dps = 80
m = 250
q = 6
epsilon = Fraction(51_063, 10_000_000)
pressure_total = Fraction(3, 1_150)
h_cert = Fraction(672_457_041_414_544_284, 10**18)
claimed = Fraction(6_732_059_784, 10_000_000_000)

def I(x: Fraction):
    return iv.mpf(x.numerator) / x.denominator

A = I(epsilon) * (m - q)
R = A / m + 2 * (((m - 1) * A / m) ** iv.mpf('0.5')) - 1
eta = R / A
bound = (m * I(h_cert) - eta * I(pressure_total) * (m - 1)) / (m - R)

print(f"m={m}")
print(f"epsilon={epsilon} = {float(epsilon):.10f}")
print(f"pressure_total={pressure_total}")
print(f"A={A}")
print(f"R={R}")
print(f"eta={eta}")
print(f"bound={bound}")
print(f"percent={100 * bound}")
print(f"claimed_lower={claimed}")
print(f"verified={bound.a > I(claimed).b}")
