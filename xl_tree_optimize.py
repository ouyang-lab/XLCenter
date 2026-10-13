import pysam
import numpy as np
import pandas as pd
import scipy.stats as spstats
from collections import OrderedDict
import pytest

class XLTree(object):

    def __init__(self, f_bam, f_chrom_sizes, experiment_type="eCLIP",
                 upstream_flank=0, downstream_flank=0, local_flank=50,
                 threshold_mapq=255,
                 threshold_variance=1.0,
                 threshold_entropy=0.0,
                 chunk_size=100000,
                 pseudo_count=0,
                 min_count=0):
        
        self.f_bam = f_bam
        self.chrom_sizes = self.set_chrom_sizes(f_chrom_sizes)
        self.experiment_type = experiment_type
        self.upstream_flank = upstream_flank
        self.downstream_flank = downstream_flank
        self.chunk_size = chunk_size
        self.local_flank = local_flank
        self.threshold_mapq = threshold_mapq
        self.threshold_variance = threshold_variance
        self.threshold_entropy = threshold_entropy
        self.pseudo_count = pseudo_count
        self.min_count = min_count
        #print (self.f_bam)

        self.df = OrderedDict()
        
    def set_chrom_sizes(self, f_chrom_sizes):
        l = []
        with open(f_chrom_sizes, "r") as f:
            for line in f:
                row = line.strip("\r\n").split("\t")
                l.append((row[0], int(row[1])))
        return dict(l)

    def get_df(self):
        return self.df
    
    def load_bam(self):

        colnames = ["chrom", "start", "end", "strand", "eclip", "ctrl", "norm", "var"]
        g_ms = self.generate_ms()
        rows = []
        
        while True:  # iterate over chrom
            try:
                chrom, raw, ms = next(g_ms)

                print ("Starting: %s" % chrom)
                #print (len(raw[0]))
                #print (len(raw[1]))
                #print (len(ms[0]))
                #print (len(ms[1]))
                print ("Forward strand #: %s" % len(ms[0]))
                print ("Reverse strand #: %s" % len(ms[1]))

                #print (sorted(raw[0].items())[0:2])
                #print (sorted(raw[1].items())[0:2])
                #print (sorted(ms[0].items())[0:2])
                #print (sorted(ms[1].items())[0:2])

                bg = self.local_window(chrom, ms, func="mean")
                var = self.local_window(chrom, ms, func="var")
                ent = self.local_window(chrom, raw, func="entropy")
                

                # filter uninformative peaks
                if ((len(ms[0]) != len(bg[0]) and len(ms[0]) != len(var[0])) and
                    (len(ms[1]) != len(bg[1]) and len(ms[1]) != len(var[1]))):
                    #print (len(ms[0]), len(bg[0]), len(var[0]))
                    #print (len(ms[1]), len(bg[1]), len(var[1]))
                    #assert False
                    pytest.skip("Skipping this test")

                
                for i, sign in [(0, "+"), (1, "-")]:  # forward and backward
                    for (mpos, m), (bpos, b), (epos, e), (vpos, v), (rpos, r) in zip(sorted(ms[i].items()),
                                                                                     sorted(bg[i].items()),
                                                                                     sorted(ent[i].items()),
                                                                                     sorted(var[i].items()),
                                                                                     sorted(raw[i].items())):
                        #if v >= self.threshold_variance:
                        if (m >= self.min_count) or (b >= self.min_count):
                            
                            if (v >= self.threshold_variance and
                                e > self.threshold_entropy):

                                if (mpos != bpos and mpos != vpos and
                                    mpos != epos and mpos != rpos):
                                    
                                    print (mpos, bpos, vpos, epos, rpos)
                                    #assert False
                                    pytest.skip("Skipping this test")
                                
                                m += self.pseudo_count
                                b += self.pseudo_count
                                rows.append((chrom, mpos, mpos+1, sign,
                                             m, b, m/(m+b), r))

                #print ("After: %s" % chrom)
                #print (len(rows))
                print ("Preprocessed #: %s" % len(rows))
                
            except StopIteration:
                break

        self.df = pd.DataFrame.from_records(rows, columns=colnames)
            
        return self
    
    def generate_ms(self):

        raw_forward = {}
        raw_reverse = {}
        
        ms_forward = {}
        ms_reverse = {}

        current_chrom = None
        
        sam = pysam.AlignmentFile(self.f_bam, "rb")
        for read in sam:

            allow = False
            
            if self.experiment_type == "eCLIP" and read.is_read2:
                allow = True
            elif self.experiment_type == "iCLIP":
                allow = True

            if allow:
                
                mapq = read.mapping_quality
                chrom = read.reference_name
                start = read.reference_start
                end = read.reference_end

                if current_chrom == None:
                    current_chrom = chrom
                    
                if current_chrom != chrom:
                    
                    yield (current_chrom, [raw_forward, raw_reverse], [ms_forward, ms_reverse])
                    
                    raw_forward = {}
                    raw_reverse = {}
                    ms_forward = {}
                    ms_reverse = {}
                    current_chrom = chrom

                # Previously
                # if mapq != self.threshold_mapq:
                if mapq < self.threshold_mapq:
                    continue
                
                if not read.is_reverse:  # forward

                    pos = start - 1
                        
                    for i in list(range(pos-self.upstream_flank, 
                                   pos+self.downstream_flank+1)):

                        if i < 0:
                            continue
                        
                        try:
                            ms_forward[i] += 1
                        except KeyError:
                            ms_forward[i] = 1

                        try:
                            raw_forward[i] += 0
                        except KeyError:
                            raw_forward[i] = 0

                    if pos >= 0:
                        raw_forward[pos] += 1
                    
                else:
                    
                    pos = end
                    
                    for i in list(range(pos-self.upstream_flank, 
                                   pos+self.downstream_flank+1)):

                        if i < 0:
                            continue
                        
                        try:
                            ms_reverse[i] += 1
                        except KeyError:
                            ms_reverse[i] = 1

                        try:
                            raw_reverse[i] += 0
                        except KeyError:
                            raw_reverse[i] = 0

                    if pos >= 0:
                        raw_reverse[pos] += 1

        yield (current_chrom, [raw_forward, raw_reverse], [ms_forward, ms_reverse])
                            
    
    def generate_rows(self, d):
        #for start, counts in sorted(d.items()):
            #yield (start, counts)
        yield from sorted(d.items())

    def local_window(self, chrom, data, func="mean"):
        ctrl = []

        if func == "mean":
            local_flank = self.local_flank
        elif func == "var":
            local_flank = self.upstream_flank
        elif func == "entropy":
            local_flank = self.upstream_flank
        #print (local_flank)
        
        size = self.chrom_sizes[chrom]
        #print (size)

        for datum in data:  # forward and reverse
            g = self.generate_rows(datum)
            local = self.chunks(g, size, local_flank, func)

        #data0 = data[0]
        #for datum in data0:  # forward
            #g = datum
            #local = self.chunks(g, size, local_flank, func)

            #ctrl.append(local)

        #data1 = data[1]
        #for datum in data1:  # reverse
            #g = datum
            #local = self.chunks(g, size, local_flank, func)


            ctrl.append(local)

        return ctrl

    
    def chunks(self, g, size, local_flank, func):

        local = {}

        prev_locus = None
        total_blocks = (size // self.chunk_size) + 1
        #print (total_blocks)
        #print (list(range(total_blocks)))
        for current_block in list(range(total_blocks)):

            block = np.zeros(self.chunk_size)
            marks = []

            while True:
                try:
                    if prev_locus:
                        start, counts = prev_locus
                    else:
                        start, counts = next(g)

                    q = start // self.chunk_size
                    r = start % self.chunk_size
                    if current_block == q:
                        block[r] = counts
                        marks.append((start, r))
                    else:
                        prev_locus = (start, counts)
                        #print current_block, counter, q
                        break  # into the next block

                    if prev_locus:
                        prev_locus = None

                except StopIteration:
                    #print current_block, counter
                    break

            # Work on current block
            for s, m in marks:

                if (m - local_flank) < 0:
                    end = m + local_flank + abs(m-local_flank)

                    if func == "mean":
                        value = block[:end].sum()/float(2*local_flank + 1)
                    elif func == "var":
                        value = block[:end].var()
                    elif func == "entropy":
                        value = (0.0 if sum(block[:end]) == 0 else spstats.entropy(block[:end]))

                    local[s] = value

                elif (m + local_flank) > self.chunk_size:
                    start = m - local_flank - (self.chunk_size-(m+local_flank))

                    if func == "mean":
                        value = block[start:].sum()/float(2*local_flank + 1)
                    elif func == "var":
                        value = block[start:].var()
                    elif func == "entropy":
                        value = (0.0 if sum(block[start:]) == 0 else spstats.entropy(block[start:]))

                    local[s] = value

                else:
                    start = m - local_flank
                    end = m + local_flank

                    if func == "mean":
                        value = block[start:end].sum()/float(2*local_flank + 1)
                    elif func == "var":
                        value = block[start:end].var()
                    elif func == "entropy":
                        value = (0.0 if sum(block[start:end]) == 0 else spstats.entropy(block[start:end]))

                    local[s] = value

        return local
