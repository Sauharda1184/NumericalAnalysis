from main import find_all_roots, horner

# Each polynomial's coefficients, ascending by power (matches main.py's
# convention), chosen to exercise every path through find_all_roots:
#   - a plain quadratic: hits the degree==2 base case directly, no
#     Newton/Muller/deflation loop at all.
#   - three real roots: exercises Newton + deflation twice in a row.
#   - two complex conjugate pairs: Newton can't start here (p'(0)=0 on
#     x^2+1), so this exercises the Muller fallback.
#   - three real roots + one complex pair: the general degree-5 case,
#     mixing both root-finding paths and multiple deflation steps.
TEST_CASES = [
    [6, -5, 1],                  # (x-2)(x-3):              roots 2, 3
    [-6, 11, -6, 1],             # (x-1)(x-2)(x-3):         roots 1, 2, 3
    [4, 0, 5, 0, 1],             # (x^2+1)(x^2+4):          roots +/-i, +/-2i
    [-24, 44, -30, 15, -6, 1],   # (x-1)(x-2)(x-3)(x^2+4):  roots 1, 2, 3, +/-2i
]

RESIDUAL_TOL = 1e-4


def test_find_all_roots():
    all_passed = True

    for coefficients in TEST_CASES:
        degree = len(coefficients) - 1
        roots = find_all_roots(coefficients)

        # A degree-n polynomial must have exactly n roots (counted with
        # multiplicity); losing or duplicating one here would mean
        # synth_div deflated incorrectly.
        count_ok = len(roots) == degree

        # Each reported root must actually satisfy p(root) ~ 0.
        residuals = [abs(horner(coefficients, r)[0]) for r in roots]
        residuals_ok = all(res < RESIDUAL_TOL for res in residuals)

        passed = count_ok and residuals_ok
        all_passed = all_passed and passed

        status = "PASS" if passed else "FAIL"
        print(f"{status}: coefficients={coefficients}")
        print(f"      roots={roots}")
        print(f"      residuals={residuals}")

    return all_passed


if __name__ == "__main__":
    if test_find_all_roots():
        print("\nAll tests passed.")
    else:
        print("\nSome tests FAILED.")
