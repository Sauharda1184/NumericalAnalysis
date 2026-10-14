
import math
import cmath

def g(coefficients, x): 
    value, derivative = horner(coefficients, x)

    return x - value / derivative 

def newton(coefficients, p_0):
    TOL = 1e-6
    i = 1 

    while i <= 100: 
        p = g(coefficients, p_0) 

        if (abs(p - p_0) < TOL): 
            print(f"Newton: I reached x={p} in {i} iterations. f(E)={g(coefficients, p)}") 
            return
        else: 
            i += 1 
            p_0 = p
    print("Newton's method failed. Did not reach TOL in 100 iterations.")

def muller(coefficients, p0, p1, p2):
    '''Fits a parabola through (p0,f(p0)), (p1,f(p1)), (p2,f(p2)) and takes
     the parabola's root closest to p2 as the next estimate, then slides
     the window forward. Using cmath (not math) lets the discriminant go
     negative, which is what lets this converge to complex roots even
     when p0, p1, p2 all start out real.'''
    TOL = 1e-6
    f = lambda x: horner(coefficients, x)[0]

    # h1, h2: step sizes between the three points.
    # delta1, delta2: slopes of the secants p0->p1 and p1->p2.
    # d: the parabola's leading (quadratic) coefficient, from the second
    # divided difference of the three points.
    h1 = p1 - p0
    h2 = p2 - p1
    delta1 = (f(p1) - f(p0)) / h1
    delta2 = (f(p2) - f(p1)) / h2
    d = (delta2 - delta1) / (h2 + h1)

    i = 3
    while i <= 100:
        # Parabola through the 3 points, in the local variable h = x - p2:
        #   f(p2) + b*h + d*h^2 = 0
        # b is the parabola's slope at p2; D is its discriminant's root.
        b = delta2 + h2 * d
        D = cmath.sqrt(b**2 - 4 * f(p2) * d)

        # Of the quadratic formula's two forms (b+D)/2d and (b-D)/2d,
        # pick whichever denominator has larger magnitude. That denominator
        # is E here; it avoids catastrophic cancellation when b and D are
        # close in value.
        if abs(b - D) < abs(b + D):
            E = b + D
        else:
            E = b - D

        # h is the step from p2 to the new estimate p (this is the
        # quadratic formula solved for the root nearest p2).
        h = -2 * f(p2) / E
        p = p2 + h

        if abs(h) < TOL:
            print(f"Muller: I reached x={p} in {i} iterations. f(p)={f(p)}")
            return p

        # Slide the window forward: drop p0, keep p1/p2, add the new estimate.
        p0, p1, p2 = p1, p2, p
        h1 = p1 - p0
        h2 = p2 - p1
        delta1 = (f(p1) - f(p0)) / h1
        delta2 = (f(p2) - f(p1)) / h2
        d = (delta2 - delta1) / (h2 + h1)
        i += 1

    print("Muller's method failed. Did not reach TOL in 100 iterations.")

def horner(coefficients, x):
    n = len(coefficients)-1
    y = coefficients[n]
    z = coefficients[n]
    for j in range (n-1, 0 ,-1): 
        y = x * y + coefficients[j] 
        z = x * z + y 
    y = x * y + coefficients[0] 
    return y, z

def synth_div(coefficients, x):
    n = len(coefficients)-1
    new_coefficients = [coefficients[n]]

    for j in range (n-1, 0 ,-1): 
        new_coefficients.append(x * new_coefficients[-1] + coefficients[j])
        new_coefficients.append(x * new_coefficients[-1] + coefficients[0]) 

    return new_coefficients

a = [-1, -1, 0, 1]

newton(a, 3)
