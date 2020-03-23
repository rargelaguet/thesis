#!/scratch/ricard/software/anaconda3/bin/python
#SBATCH -N 1                        # number of nodes
#SBATCH -n 1                        # number of cores
#SBATCH --mem 15G                   # memory pool for all cores
#SBATCH -o slurm.%N.%j.out          # STDOUT
#SBATCH -e slurm.%N.%j.err          # STDERR
#SBATCH --mail-type=END,FAIL        # notifications for job done & fail
#SBATCH --mail-user=ricard@ebi.ac.uk

# sbatch -p gpu --gres gpu:1 simulate.py

from __future__ import division
import scipy as s
import pandas as pd
from scipy.stats import bernoulli, norm, gamma, uniform, poisson, binom
from random import sample
import pickle
import os
from re import search

def sigmoid(X):
    return s.divide(1.,1.+s.exp(-X))

class Simulate(object):
    def __init__(self, M, N, D, K, G):
        """General method to Simulate from the generative model

        PARAMETERS
        ----------
        M (int): number of features groups (views)
        N (int): number of samples
        D (list/tuple of length M): dimensionality of each view
        K (int): number of factors
        G (int): number of samples groups
        """

        # Sanity checks
        assert len(D) == M
        assert K < min(D)

        # Sanity checks
        assert len(N) == G
        assert K < min(N)

        self.M = M
        self.N = N
        self.K = K
        self.D = D
        self.G = G

    def sampleAlphaW_discrete(self, min_value=1., max_value=100, proportion_active=0.5):
        """ Initialisation of ARD on the weights"""
        alpha = [ s.zeros(self.K,) for m in range(self.M) ]
        for m in range(self.M):
            tmp = bernoulli.rvs(p=proportion_active, size=self.K)
            tmp[tmp==1] = min_value
            tmp[tmp==0] = max_value
            alpha[m] = tmp
        return alpha

    def sampleAlphaW_gamma(self, a=1, b=1):
        """ Initialisation of ARD on the weights"""
        alpha = [ s.random.gamma(shape=a, scale=b, size=self.K) for m in range(self.M)]
        return alpha


    def sampleAlphaZ_discrete(self, min_value=1., max_value=100., proportion_active=0.5):
        """ Initialisation of ARD on the weights"""
        alpha = [ s.zeros(self.K,) for g in range(self.G) ]
        for g in range(self.G):
            tmp = bernoulli.rvs(p=proportion_active, size=self.K)
            tmp[tmp==1] = min_value
            tmp[tmp==0] = max_value
            alpha[g] = tmp
        return alpha

    def sampleAlphaZ_gamma(self, shape=1, scale=1):
        """ Initialisation of ARD on the weights"""
        alpha = [ s.random.gamma(shape=shape, scale=scale, size=self.K) for g in range(self.G)]
        return alpha

    def sampleThetaW(self, a=1, b=1):
        theta = [ s.random.beta(a, b, size=self.K) for m in range(self.M)]
        return theta

    def sampleThetaZ(self, a=1, b=1):
        theta = [ s.random.beta(a, b, size=self.K) for g in range(self.G)]
        return theta

    def sampleW(self, theta, alpha):
        """ Initialisation of weights with a spike and slab prior"""

        # Simulate bernoulli variable S
        S = [ s.zeros((self.D[m],self.K)) for m in range(self.M) ]
        for m in range(self.M):
            for k in range(self.K):
                S[m][:,k] = bernoulli.rvs(p=theta[m][:,k], size=self.D[m])

        # Simulate gaussian weights W
        W = [ s.empty((self.D[m],self.K)) for m in range(self.M) ]
        for m in range(self.M):
            for k in range(self.K):
                W[m][:,k] = norm.rvs(loc=0, scale=s.sqrt(1./alpha[m][k]), size=self.D[m])

            # Take the product S*W
            W[m] *= S[m]

        return W

    def sampleZ(self, theta, alpha):
        """ Initialisation of factors with a spike and slab prior"""

        # Simulate bernoulli variable S
        S = [ s.zeros((self.N[g],self.K)) for g in range(self.G) ]
        for g in range(self.G):
            for k in range(self.K):
                S[g][:,k] = bernoulli.rvs(p=theta[g][:,k], size=self.N[g])

        # Simulate gaussian factors Z
        Z = [ s.empty((self.N[g],self.K)) for g in range(self.G) ]
        for g in range(self.G):
            for k in range(self.K):
                Z[g][:,k] = norm.rvs(loc=0, scale=s.sqrt(1./alpha[g][k]), size=self.N[g])

            # Take the product S*Z
            Z[g] *= S[g]

        return Z

    def sampleTau_uniform(self, min_value=1., max_value=3.):
        """ Initialisation of noise precision"""
        Tau = [ [uniform.rvs(loc=min_value, scale=max_value, size=self.D[m]) for g in range(self.G)] for m in range(self.M) ]
        return Tau

    def sampleTau_gamma(self, shape=1, scale=1):
        """ Initialisation of noise precision"""
        Tau = [ [s.random.gamma(shape=shape, scale=scale, size=self.D[m]) for g in range(self.G)] for m in range(self.M) ]
        return Tau
     

    def sampleY(self, W, Z, Tau, likelihood="gaussian", missingness=0.0, zeros=0.0, missing_view=False):
        """ Initialisation of observations 

        PARAMETERS
        ----------
        W (list of length M where each element is a np array with shape (Dm,K)): weights
        Z (np array with shape (N,K): latent variables
        Tau (list of length M where each element is a np array with shape (Dm,)): precision of the normally-distributed noise
        likelihood (str): type of likelihood
        missingness (float): percentage of missing values
        """
        Y = [ [None for g in range(self.G)] for m in range(self.M) ]

        # F = [ s.zeros((self.N,self.D[m])) for m in range(self.M) ]

        if likelihood == "gaussian":
            for m in range(self.M):
                for g in range(self.G):
                    # Y[m][g] = s.dot(Z[g],W[m].T) + norm.rvs(loc=0, scale=1/s.sqrt(Tau[m][g]), size=self.D[m])
                    Y[m][g] = s.dot(Z[g],W[m].T) + norm.rvs(loc=0, scale=1/s.sqrt(Tau[m][g]), size=[self.N[g], self.D[m]])

        # Sample observations using a poisson likelihood
        elif likelihood == "poisson":
            for m in range(self.M):
                for g in range(self.G):
                    F = s.dot(Z[g],W[m].T)
                    Y[m][g] = s.special.round(s.log(1.+s.exp(F))) # Without noise
                    # Y[m] = poisson.rvs(rate).astype(float) # With noise, sample from the Poisson distribution

        # Sample observations using a bernoulli likelihood
        elif likelihood == "bernoulli":
            for m in range(self.M):
                for g in range(self.G):
                    F = s.dot(Z[g],W[m].T)
                    Y[m][g] = s.special.round(sigmoid(F)) # without noise
            #     # Y[m] = bernoulli.rvs(f).astype(float) # with noise

        # Introduce zeros
        if zeros > 0.0:
            for m in range(self.M):
                for g in range(self.G):
                    nas = s.random.choice(range(self.N[g]*self.D[m]), size=int(zeros*self.N[g]*self.D[m]), replace=False)
                    tmp = Y[m][g].flatten()
                    tmp[nas] = 0.
                    Y[m][g] = tmp.reshape((self.N[g],self.D[m]))

        # Introduce missing values into the data
        if missingness > 0.0:
            for m in range(self.M):
                for g in range(self.G):
                    nas = s.random.choice(range(self.N[g]*self.D[m]), size=int(missingness*self.N[g]*self.D[m]), replace=False)
                    tmp = Y[m][g].flatten()
                    tmp[nas] = s.nan
                    Y[m][g] = tmp.reshape((self.N[g],self.D[m]))


        if missing_view > 0.0:   # percentage of samples missing a view
            return
            # # select samples missing one view
            # n_missing = s.random.choice(range(self.N), int(missing_view * self.N), replace=False)
            # Y[0][n_missing,:] = s.nan

        # return F,Y
        return Y



if __name__ == "__main__":

    # Define I/O
    host = os.uname()[1]
    if search("ricard", host):
        outdir = '/Users/ricard/data/mofa2_simulations/test'
    elif search("embl", host):
        outdir = '/g/stegle/ricard/mofa2_simulations/test'
    else:
        print("Computer not recognised"); exit()

    D = [int(1e3),int(1e3),int(1e3)]
    M = len(D)
    K = 15
    N = [int(5e2)]
    G = len(N)

    tmp = Simulate(M=M, N=N, D=D, K=K, G=G)

    data = {}

    # Sample Z
    thetaZ = [ s.ones((N[g],K)) for g in range(G) ]
    AlphaZ = [ s.ones((K,)) for g in range(G) ]
    data['Z'] = tmp.sampleZ(theta=thetaZ, alpha=AlphaZ)

    # Sample W
    thetaW = [ s.ones((D[m],K))*0.75 for m in range(M) ]
    # AlphaW = [ s.ones((K,)) for m in range(M) ]
    AlphaW = tmp.sampleAlphaW_discrete(min_value=1., max_value=1000, proportion_active=0.35)
    # AlphaW = tmp.sampleAlphaW_gamma()
    data['W'] = tmp.sampleW(theta=thetaW, alpha=AlphaW)

    # Sample Tau
    data['Tau'] = tmp.sampleTau_uniform(min_value=1., max_value=3.)
    

    # Sample 
    data['Y'] = tmp.sampleY(W=data["W"], Z=data["Z"], Tau=data["Tau"], likelihood="gaussian")

    # Save output as a long data.table format
    dt = list()
    for m in range(M):
        for g in range(G):
            foo = pd.DataFrame(data["Y"][m][g])
            foo["group"],foo["view"] = ["group_"+str(g), "view_"+str(m)]
            foo["sample"] = s.char.mod('sample_%d', foo.index) + ("_"+foo["group"].values)
            bar = pd.melt(foo, id_vars=["sample","group","view"], var_name="feature", value_name='value')
            bar["feature"] = s.char.mod('feature_%d', bar["feature"].values) + ("_"+bar["view"].values)
            dt.append(bar)
    dt = pd.concat(dt)[["sample","group","feature","view","value"]]

    outfile = "%s/data.txt.gz" % (outdir)
    dt.to_csv(outfile, float_format='%.2f', sep='\t', index=False)
