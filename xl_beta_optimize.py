import sys
import gzip
import pandas as pd
import numpy as np
import scipy.stats as spstats
from sklearn.mixture import GaussianMixture
from statsmodels.stats.multitest import multipletests

from xl_betamm import BetaMixture

class XLBeta(object):

    def __init__(self):
        
        self.n_components = 2
        
        self.strand = None
        self.df = {}
        self.Z = {}
        self.use_col = {}
        self.bed = {}
        self.norm = {}
        self.params = {}
        self.posterior = {}

        self.chunk_size = 100000
    
    def set_df(self, df):
        self.df = df
        return self
    
    def fit(self):

        """ !!!Restore once confirmed!!!
        norm = self.df["norm"]
        zero = np.zeros(len(norm))
        data = np.stack((norm, zero), axis=1)
        """

        eclip = self.df["eclip"]
        ctrl = self.df["ctrl"]

        p_vals = spstats.poisson.cdf(eclip, ctrl)
        
        """ Use P values instead of Norms from now onwards"""
        self.df["norm"] = p_vals
        
        zero = np.zeros(len(p_vals))
        data = np.stack((p_vals, zero), axis=1)
        
        model = BetaMixture()
        model.fit(data[:,0])
        model.infer_stats()
        
        self.Z = model.predict_proba()

        self.use_col = 0
        
        self.params = {
            0: (model.means_[1], model.covariances_[1]),
            1: (model.means_[0],model.covariances_[0])}
        
        """
        sys.stderr.write("Converged: %s\n" % model.converged_)
        sys.stderr.write("Weights:\n%s\n" % model.weights_)
        sys.stderr.write("Means:\n%s\n" % model.means_)
        sys.stderr.write("Covariances:\n%s\n" % model.covariances_)
        """
        
        return self
    
    def output_posterior(self, output_prefix, posterior_threshold):
        f_posterior = output_prefix.rstrip("/") + ".txt"

        with open(f_posterior, "w") as o_posterior:
             # Write means and variance to header.
            o_posterior.write("# 0,%s,%s\n" % self.params[0])
            o_posterior.write("# 1,%s,%s\n" % self.params[1])
            
            colnames  = ["chrom", "start", "end", "eclip", "ctrl", "call",
                    "norm","posterior", "strand"]
            line = "\t".join(colnames) + "\n"
            o_posterior.write(line)

            for d, z in zip(self.df.values, self.Z):
                chrom, start, end, sign = (
                    d[0], int(d[1]), int(d[2]), d[3]
                )
                eclip, ctrl, norm = (
                    int(d[4]), int(d[5]), d[6]
                )
                posterior = z[self.use_col]

                call = 1 if posterior >= posterior_threshold else 0

                out = [
                    chrom, start, end, eclip, ctrl, call,
                    norm, posterior, sign
                ]

                o_posterior.write(
                    "%s\n" % "\t".join(map(str, out))
                )

    def load_posterior(self, f_posterior):

        # =========================
        # Load Posterior File
        # =========================
        colnames = ["chrom", "start", "end", "eclip", "ctrl", "call",
                    "norm", "posterior", "strand"]
        
        df = pd.read_table(f_posterior, sep="\t", skiprows=3,
                           names=colnames)
        
        self.df = df[["chrom", "start", "end", "strand",
                      "eclip", "ctrl", "norm"]]
        self.Z = df[["posterior"]].values
        self.use_col = 0

        # ==========================
        # Load Gaussian Parameters
        # ==========================
        with open(f_posterior, "r") as f:
            self.params = {}
            for _ in range(2):
                row = next(f).strip("\r\n").split()[1].split(",")
                self.params[int(row[0])] = (float(row[1]), float(row[2]))
        
        return self


    def call_bed(self, posterior_threshold, motif_size, st="both"):
        
        self.bed = []
        self.norm = []
        self.posterior = []
        for strand, sign in [(0, "+"), (1, "-")]:

            current_chrom = None
            current_start = None
            current_end = None
            current_freq = None
            current_posterior = None
            counter = 0
            
            for d, z in zip(self.df[self.df["strand"] == sign].values,
                            self.Z[self.df["strand"] == sign]):

                chrom, start, end = (d[0], int(d[1]), int(d[2]))
                eclip, ctrl,norm = (int(d[4]), int(d[5]), d[6])
                posterior = z[self.use_col]

                if posterior >= posterior_threshold:
                    call = 1
                else:
                    call = 0

                if call == 1:

                    # initialize current locus
                    if current_chrom == None:
                        current_chrom = chrom
                        current_start = start
                        current_end = start
                        current_norm = [norm]
                        current_posterior = [posterior]
                        continue
                    
                    # change in chrom: new region
                    if current_chrom != chrom:

                        # check size and store bed region
                        if (current_end - current_start + 1) >= motif_size:
                            name = "RBPS_%s_%i" % (strand, counter)
                            self.bed.append([current_chrom, current_start,
                                             current_end+1, name, 0, sign])
                            self.norm.append(current_norm)
                            self.posterior.append(current_posterior)
                            counter += 1

                        current_chrom = chrom
                        current_start = start
                        current_end = start
                        current_norm = [norm]
                        current_posterior = [posterior]

                    elif start-1 != current_end:

                        # check size and store bed region
                        if (current_end - current_start + 1) >= motif_size:
                            name = "RBPS_%s_%i" % (strand, counter)
                            self.bed.append([current_chrom, current_start,
                                             current_end+1, name, 0, sign])
                            self.norm.append(current_norm)
                            self.posterior.append(current_posterior)
                            counter += 1
                        
                        current_chrom = chrom
                        current_start = start
                        current_end = start
                        current_norm = [norm]
                        current_posterior = [posterior]

                    else:
                        current_end = start
                        current_norm.append(norm)
                        current_posterior.append(posterior)

            
            # check size and store bed region
            if current_chrom and (current_end - current_start + 1) >= motif_size:
                name = "RBPS_%s_%i" % (strand, counter)
                self.bed.append([current_chrom, current_start,
                                 current_end+1, name, 0, sign])
                self.norm.append(current_norm)
                self.posterior.append(current_posterior)
                

    def output_bed(self, f_name=None, prefix=None, st="both"):
        output_name = None

        if f_name:
            output_name = f_name
        if prefix:
            output_name = prefix.rstrip("/") + ".bed"

        if output_name is None:
            raise ValueError("Provide either f_name or prefix.")
            
        with open(output_name, "w") as o_bed:
            for out in sorted(
                self.bed, key=lambda x: (x[0], x[1], x[2])
            ):
                o_bed.write(
                    "%s\n" % "\t".join(map(str, out))
                )
            

    # ============================
    # Calc Significance LRT
    # ============================
    def calc_significance(self, posterior_threshold, distribution="beta"):
        
        pvals = []
        dofps = []
        for x, z, out in zip(self.norm, self.posterior,
                             self.bed):
            
            if distribution == "norm":
                L0 = spstats.norm.logpdf(x,loc=self.params[0][0],
                                         scale=self.params[0][1]).sum()
                L1 = spstats.norm.logpdf(x,loc=self.params[1][0],
                                         scale=self.params[1][1]).sum()
                dof = 0
            elif distribution == "beta":
                L0 = spstats.beta.logpdf(x,self.params[0][0],
                                         self.params[0][1]).sum()
                L1 = spstats.beta.logpdf(x,self.params[1][0],
                                         self.params[1][1]).sum()
            else:
                raise ValueError(
                    "distribution must be 'norm' or 'beta'"
                )
                
                #dof = 1
                #dof = 2 * len(x)
            
            d = 2 * (L1 - L0)
            
            pval = spstats.chi2.sf(d, 1)
            if pval == 0.0:
                pval = 1e-313
            elif np.isnan(pval):
                pval = 1e-313
            
            logp = -np.log10(pval)
            score = round((1-pval)*1000)
            signal = "%.5f" % np.array(z).mean()

            dofp = spstats.chi2.sf(d, len(x))
            if dofp == 0.0:
                dofp = 1e-313
            elif np.isnan(dofp):
                dofp = 1e-313
            
            pvals.append(pval)
            dofps.append(dofp)
            
            out[4] = int(score)
            out.append(signal)
            out.append(logp)

        if len(pvals) != 0:
            rejection, qvals, sidak, bonf = multipletests(pvals, method="fdr_bh")
            logq = -np.log10(qvals)
            logd = -np.log10(dofps)
            for qval, logdofp, out in zip(logq, logd, self.bed):
                out.append(qval)
                out.append(logdofp)


    
    # ======================================
    # TO-REMOVE Calc Significance in Chunks
    # ======================================
    def generate_rows(self, d):
        for (chrom, start), counts in sorted(d.items()):
            yield (chrom, start, counts)
    
    def generate_beds(self, bed):
        for row in bed:
            chrom, start, end, name, sign = (row[0], row[1], row[2], row[3], row[4])
            yield (chrom, start, end, name, sign)
    
    def calc_significance_old(self, raw_data, local_window):
        fr = []
        chroms = self.chroms
        chrom_sizes = self.chrom_sizes
        chunk_size = self.chunk_size
        local_flank = local_window//2
        local_window = local_window
        
        for (raw_datum, bed, chrom) in zip(raw_data, self.bed.values(), chroms):
            g = self.generate_rows(raw_datum)
            b = self.generate_beds(bed)
            sig = self.calc_significance_from_chunks(g, b, chrom,
                                                     chrom_sizes, chunk_size,
                                                     local_flank, local_window)
            fr.append(sig)
        self.bed = {}
        self.bed["forward"] = fr[0]
        self.bed["reverse"] = fr[1]
        
    def calc_significance_from_chunks(self, g, b, chroms,
                                      chrom_sizes, chunk_size,
                                      local_flank, local_window):
        sig = []
        pval = []
        g_chrom, g_start, g_counts = next(g)
        b_chrom, b_start, b_end, b_name, b_sign = next(b)
        b_size = (b_end - b_start) + local_window
        b_block = np.zeros(b_size, dtype=int)
        while True:
            try:
                
                if ((g_chrom, g_start) < (b_chrom, b_start-local_flank+1)):
                    g_chrom, g_start, g_counts = next(g)
                elif (((g_chrom, g_start) >= (b_chrom, b_start-local_flank+1)) and
                      ((g_chrom, g_start) < (b_chrom, b_end+local_flank+1))):
                    r = g_start - (b_start-local_flank+1)
                    b_block[r] = g_counts
                    g_chrom, g_start, g_counts = next(g)
                else:
                    """
                    print b_block
                    """

                    bed_size = (b_end - b_start)
                    raw = b_block[local_flank:local_flank+(b_end-b_start)].sum()
                    mean, std = self.permute(b_block, bed_size, times=100)
                    if mean != 0.0 and std != 0.0:
                        cdf = spstats.norm.cdf((raw-mean)/std)
                        score = int(round(cdf * 1000))
                        fc = raw/mean
                        p_val = 1.0 - cdf
                        log_p_val = -np.log10(p_val)
                        q_val = -1

                        pval.append(p_val)
                        sig.append([b_chrom, b_start, b_end, b_name,
                                    score, b_sign, "%.5f" % fc, "%.5f" % log_p_val])
                    
                    b_chrom, b_start, b_end, b_name, b_sign = next(b)
                    b_size = (b_end - b_start) + local_window
                    b_block = np.zeros(b_size, dtype=int)
            
            except StopIteration:
                break

        # Adjust P-value using Benjamini-Hochberg for positive correlation
        rejection, qvals, sidak, bonf = multipletests(pval, method="fdr_bh")
        for qval, s in zip(qvals, sig):
            log_q_val = -np.log10(qval)
            summit = -1
            s.append(log_q_val)
            s.append(summit)
        
        return sig

    def permute(self, b_block, bed_size, times=100):
        block_size = len(b_block)
        permutation = []
        for i in np.random.randint(0, block_size, times):
            permutation.append(np.roll(b_block, i)[:bed_size].sum())
        permutation = np.array(permutation)
        return permutation.mean(), permutation.std()
            
            
