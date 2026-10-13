import numpy as np
import scipy as sp
import scipy.stats as spstats

from scipy.special import polygamma
from scipy.special import gammaln


class BetaMixture(object):

    def __init__(self, alpha=None, beta=None, step=0.01):
        self.prev_alpha = alpha
        self.prev_beta = beta
        self.prev_L = -np.inf
        self.iteration = 0
        self.tolerance = 0.01
        self.step_size = step

        self.prev_H = np.array([self.prev_alpha, self.prev_beta])
        
    def initialize_parameters(self, X):
        """
        ssize = 100000 #10000
        if len(X) >= ssize:
            D = np.random.choice(X, ssize, replace=False)
        else:
            D = X
        """
        D = X
        if self.prev_alpha == None and self.prev_beta == None:
            theta = spstats.beta.fit(D)
            self.prev_alpha = theta[0]
            self.prev_beta = theta[1]
            self.prev_H = np.array([self.prev_alpha, self.prev_beta])
            
    def initialize_states(self):
        s = np.count_nonzero(self.prev_P[0] >= 0.5)/float(len(self.prev_P[0]))
        self.prev_S = np.array([s, 1-s])

    def initialize_posterior(self, X):
        self.prev_P = np.random.rand(2, len(X))  # posterior

    def infer_stats(self):
        self.means_ = np.array([0.0, 0.0])
        self.covariances_ = np.array([0.0, 0.0])

        beta_params = spstats.beta.stats(self.prev_alpha, self.prev_beta)
        unif_params = spstats.beta.stats(1.0, 1.0)

        """
        self.means_[0] = float(beta_params[0])
        self.means_[1] = float(unif_params[0])

        self.covariances_[0] = float(beta_params[1])
        self.covariances_[1] = float(unif_params[1])
        """

        self.means_[0] = self.prev_alpha
        self.means_[1] = 1.0

        self.covariances_[0] = self.prev_beta
        self.covariances_[1] = 1.0

    
    def fit(self, X):

        X[X==0.0] = 0.0001
        X[X==1.0] = 1-0.0001

        self.initialize_parameters(X)
        self.initialize_posterior(X)
        self.initialize_states()
        
        while True:

            curr_N0 = self.prev_S[0] * spstats.beta.pdf(X, self.prev_alpha, self.prev_beta)
            curr_N1 = self.prev_S[1] * spstats.beta.pdf(X, 1, 1)

            curr_P0 = (curr_N0)/(curr_N0 + curr_N1)
            curr_P1 = 1 - curr_P0
            curr_P = np.array([curr_P0, curr_P1])

            curr_S0 = curr_P0.mean()
            curr_S1 = curr_P1.mean()

            curr_S = np.array([curr_S0, curr_S1])
            
            L0 = (curr_P[0] *
                  (gammaln(self.prev_alpha) + gammaln(self.prev_beta) - gammaln(self.prev_alpha+self.prev_beta)
                   + ((self.prev_alpha-1) * np.log(X)) + ((self.prev_beta-1) * np.log(1-X)))).sum()
            L1 = (curr_P[1] * 0).sum()
            
            curr_L = L0 + L1


            J0 = (curr_P[0] * (polygamma(0, self.prev_alpha)
                               - polygamma(0, self.prev_alpha+self.prev_beta) + np.log(X))).sum()
            J1 = (curr_P[0] * (polygamma(0, self.prev_beta)
                               - polygamma(0, self.prev_alpha+self.prev_beta) + np.log(1-X))).sum()

            J = np.array([J0, J1])

            H00 = (curr_P[0] * (polygamma(1, self.prev_alpha)
                                - polygamma(1, self.prev_alpha+self.prev_beta))).sum()
            H01 = (curr_P[0] * (-polygamma(1, self.prev_alpha+self.prev_beta))).sum()
            H11 = (curr_P[0] * (polygamma(1, self.prev_beta)
                                - polygamma(1, self.prev_alpha+self.prev_beta))).sum()

            H = np.array([[H00, H01],
                          [H01, H11]])

            R = np.matmul(np.linalg.inv(H), J)
            
            curr_H = self.prev_H - (self.step_size * R)

            change = curr_L-self.prev_L

            print (curr_H, self.prev_H, self.step_size, R, curr_L, self.prev_L, change)

            self.iteration += 1

            if self.iteration > 100:
                break
            
            if change < self.tolerance:
                break

            # Update
            self.prev_L = curr_L
            self.prev_S = curr_S
            self.prev_P = curr_P
            self.prev_H = curr_H

            self.prev_alpha = self.prev_H[0]
            self.prev_beta = self.prev_H[1]

            
    def predict_proba(self):
        return self.prev_P.T
    

if __name__ == "__main__":
    
    
    #X = np.random.rand(10)
    #X = np.array([0.15, 0.2, 0.13, 0.9, 0.94, 0.81, 0.4, 0.2, 0.12, 0.74])
    X = np.array([0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.95, 0.91, 0.93, 0.92, 0.91])
    X = np.concatenate([X, X, X])
    
    model = BetaMixture()
    model.fit(X)
    Z = model.predict_proba()
    
    for i, j in zip(X, Z):
        print (i, j)

    
    
    
    
