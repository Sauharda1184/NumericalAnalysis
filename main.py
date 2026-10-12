
import math 

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
