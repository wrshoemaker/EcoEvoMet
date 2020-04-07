from __future__ import division
import numpy as np
import math

# manipulate S to control excreted concentraiton, manipulating beta
def per_cell_metab_excret(X_new, X_ac = 0.0001, S_ac = 0.25):
    if X_new >= X_ac:
        return S_ac * (X_new - X_ac)
    else:
        return

def run_sim():

    iter = 10
    timesteps = 100
    N = 10000
    n1=5000
    n2=5000

    alpa_1_1 = 0.8
    alpa_1_2 = 1 - alpa_1_1

    alpa_2_1 = 0.2
    alpa_2_2 = 1 - alpa_2_1

    X_1 = 1
    for x in range
    X_2 = 1

    X_2_delta = 0.1 * X_2

    X_2 = X_2 - X_2_delta

    print(X_1 - X_2)

    # assume that theyre substituitable
    c_1 = 1000

    c_2 = (n1*per_cell_metab_excret(X_1)) + (n2*per_cell_metab_excret(X_2))

    beta_1 = c_1 / (c_1 + c_2)
    beta_2 = c_2 / (c_1 + c_2)

    for x in range(timesteps):

        mean_X_1 = np.log((alpa_1_1 * math.exp(X_1)/beta_1 * (n1/(n1+n2))) + \
                    (alpa_2_1 * math.exp(X_2)/beta_1 * (n2/(n1+n2))))
        mean_X_2 = np.log((alpa_1_2 * math.exp(X_1)/beta_2 * (n1/(n1+n2))) + \
                    (alpa_2_2 * math.exp(X_2)/beta_2 * (n2/(n1+n2))))

        lambda_1 = (alpa_1_1 * math.exp(X_1-mean_X_1)) + (alpa_1_2 * math.exp(X_1-mean_X_2))
        lambda_2 = (alpa_2_1 * math.exp(X_2-mean_X_1)) + (alpa_2_2 * math.exp(X_2-mean_X_2))

        lambda_1_new = lambda_1 * N/(lambda_1+lambda_2)
        lambda_2_new = lambda_2 * N/(lambda_1+lambda_2)

        n1 = np.random.poisson(lam=lambda_1_new)
        n2 = np.random.poisson(lam=lambda_2_new)


    print(n1 / (n1+n2))



run_sim()

#




#np.random.normal(loc=0.0, scale=s0)
