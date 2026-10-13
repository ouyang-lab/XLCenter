import sys
import argparse

from xl_tree_optimize import XLTree
from xl_beta_optimize import XLBeta

def main():

    parser = argparse.ArgumentParser(description=__doc__)
    subparser = parser.add_subparsers(help="Type <sub-command> -h for help ",
                                      dest="sub")

    # =====================================
    # Threshold - Description of Arguments
    # =====================================
    sub_threshold = subparser.add_parser("threshold",
                                         help="use a diff posterior threshold")
    sub_threshold.add_argument("-p", required=True,
                               help="Posterior File")
    sub_threshold.add_argument("-b", required=True,
                               help="Bed File")
    sub_threshold.add_argument("-t", default=0.5, type=float,
                               help="Posterior Threshold")
    sub_threshold.add_argument("-v", default=7, type=int,
                               help="Motif Size Window")
    sub_threshold.add_argument("-d", default="norm",
                               help="Distribution Type")


    # ================================
    # Dump - Description of Arguments
    # ================================
    sub_dump = subparser.add_parser("dump",
                                    help="dump normalized moving sum as bedGraph")
    sub_dump.add_argument("-i", required=True,
                          help="Input BAM (required)")
    sub_dump.add_argument("-o", required=True,
                          help="Filename of bedGraph output (required)")
    sub_dump.add_argument("-f", required=True,
                          help="Chrom Sizes file (Required)")
    sub_dump.add_argument("-e", required=True,
                          help="Type of CLIP Experiment [eCLIP|iCLIP] (Required)")
    sub_dump.add_argument("-n", default=0.0, type=float,
                          help="Pseudo count")
    sub_dump.add_argument("-m", default=3.0, type=float,
                          help="Minimum count")
    sub_dump.add_argument("-v", default=7, type=int,
                          help="Motif Size Window")
    sub_dump.add_argument("-w", default=100, type=int,
                          help="Local Window")
    sub_dump.add_argument("-q", default=255, type=int,
                          help="MAPQ Threshold")
    sub_dump.add_argument("-a", default=0.0, type=float,
                          help="Variance Threshold")

    # =====================================
    # callsites - Description of Arguments
    # =====================================
    sub_callsites = subparser.add_parser("callsites",
                                         help="Call RBP binding sites")
    sub_callsites.add_argument("-i", required=True,
                               help="CLIP BAM (Required)")
    sub_callsites.add_argument("-o", required=True,
                               help="Output Prefix (Required)")
    sub_callsites.add_argument("-f", required=True,
                               help="Chrom Sizes file (Required)")
    sub_callsites.add_argument("-e", required=True,
                               help="Type of CLIP Experiment [eCLIP|iCLIP] (Required)")
    sub_callsites.add_argument("-n", default=0.0, type=float,
                               help="Pseudo count")
    sub_callsites.add_argument("-m", default=3.0, type=float,
                               help="Minimum count")
    sub_callsites.add_argument("-q", default=255, type=int,
                               help="MAPQ Threshold")
    sub_callsites.add_argument("-t", default=0.5, type=float,
                               help="Posterior Threshold")
    sub_callsites.add_argument("-a", default=0.0, type=float,
                               help="Variance Threshold")
    sub_callsites.add_argument("-v", default=7, type=int,
                               help="Motif Size Window")
    sub_callsites.add_argument("-w", default=100, type=int,
                               help="Local Window")
    sub_callsites.add_argument("-d", default="bum",
                               help="Model [bum]")
    
    args = parser.parse_args()
    
    if args.sub == "threshold":
        f_posterior = args.p
        f_bed = args.b
        posterior_threshold = args.t
        motif_size = args.v
        distribution = args.d
        
        model = XLBeta()
        model.load_posterior(f_posterior)
        model.call_bed(posterior_threshold, motif_size)
        model.calc_significance(posterior_threshold, distribution)
        model.output_bed(f_name=f_bed)
        print ("Thresholding completed.")

    if args.sub == "dump":
        f_eclip = args.i
        f_bedGraph = args.o
        min_count = args.m
        pseudo_count = args.n
        experiment = args.e
        motif_size = args.v
        local_window = args.w
        local_flank = int(local_window/2)
        f_chrom_sizes = args.f
        threshold_variance = args.a
        threshold_mapq = args.q
        
        tree = XLTree(f_eclip, f_chrom_sizes, experiment,
                      upstream_flank=3, downstream_flank=3,
                      local_flank=local_flank, threshold_mapq=threshold_mapq,
                      threshold_variance=threshold_variance,
                      pseudo_count=pseudo_count,
                      min_count=min_count).load_bam()
        df = tree.get_df()
        with open(f_bedGraph, "w") as out:
            for row in df.itertuples():
                if row[4] == "+":
                    out.write("\t".join(map(str, [row[1], row[2], row[3], row[5]])) + "\n")
                else:
                    out.write("\t".join(map(str, [row[1], row[2], row[3], -1.0*row[5]])) + "\n")
        
    elif args.sub == "callsites":
        
        posterior_threshold = args.t
        motif_size = args.v
        min_count = args.m
        pseudo_count = args.n
        experiment = args.e
        local_window = args.w
        local_flank = int(local_window/2)
        f_chrom_sizes = args.f
        threshold_variance = args.a
        threshold_mapq = args.q
        model = args.d
        
        f_eclip = args.i
        output_prefix = args.o

        # ==========================
        # Stage 1 - Load data
        # ==========================
        print ("Stage 1 - loading data...")
        tree = XLTree(f_eclip, f_chrom_sizes, experiment,
                      upstream_flank=3, downstream_flank=3,
                      local_flank=local_flank, threshold_mapq=threshold_mapq,
                      threshold_variance=threshold_variance,
                      pseudo_count=pseudo_count,
                      min_count=min_count).load_bam()
        df = tree.get_df()
                    
        # =========================
        # Stage 2 - Call boundaries
        # =========================
        print ("Stage 2 - calling boundaries...")
        model = XLBeta()

        model.set_df(df)
        model.fit()
        model.call_bed(posterior_threshold, motif_size)

        # =======================================
        # Stage 3 - Calc Significance
        # =======================================
        print ("Stage 3 - calculating significance...")
        model.calc_significance(posterior_threshold)
        
        # =======================================
        # Stage 4 - Output
        # =======================================
        print ("Stage 4 - output files...")
        model.output_posterior(output_prefix, posterior_threshold)
        model.output_bed(prefix=output_prefix)
        print ("Task completed.")
        
if __name__ == "__main__":

    main()
