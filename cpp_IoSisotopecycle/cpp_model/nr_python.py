
import time
import pandas as pd
import numpy as np   


# newton raphson solver
def newton_raphson(x0,constants,e1,step,eqs,deriv,maxiter=1000):

    def dx(x,eqs):
        f_ = eqs(x,constants)
        result =(abs(0-f_))
        return result
    
    def nr(x0,eqs,deriv):
        f_ = eqs(x0,constants)
        df_ = deriv
        x0 = x0 - step*(f_/df_)
        return x0
    
    # create results table        
    results = pd.DataFrame([["guessx","diff","step","f","df","f/df"]])  
    
    delta1 = dx(x0,eqs)
    n=1
    results1 = pd.DataFrame([[x0,delta1,step]]) 
    results = pd.concat([results, results1], ignore_index=True)
    
    while delta1 > e1:
        f_ = eqs(x0,constants)
        df_ = deriv(x0,constants)
        x0 = x0 - step*(f_/df_)
        n=n+1.
        #while x0 < 0.:
        #    step = step/10.
        #    x0 = x0 - step*(f_/df_)
        delta1 = dx(x0,eqs)      
        results1 = pd.DataFrame([[x0,delta1,step,f_,df_,f_/df_]])
        results = pd.concat([results, results1], ignore_index=True)  
        if n > 10:
            results.to_csv('results_nr.csv', index=False, header=False) 
    return x0    


start_time = time.perf_counter()

print("\n=== Test newton_raphson (module version) ===")

# Solve x^2 - A = 0 with A=2  (root = sqrt(2))
constants = {"A": 2.0}

def eqs(x, c):
    return x*x - c["A"]

def deriv(x, c):
    return 2.0*x

x0   = 1.0    # initial guess
e1   = 1e-12  # tolerance
step = 1.0    # NR step scaling

root = newton_raphson(x0, constants, e1, step, eqs, deriv)
err  = abs(root - 2.0**0.5)

print("root  =", root)
print("error =", err)

end_time = time.perf_counter()
print(f"\nExecution time: {end_time - start_time} seconds")

