import pysam
import csv
import edlib
import subprocess
import os
import logging
import shutil
import tempfile
import argparse
import sys
from collections import Counter
from tqdm import tqdm
import pandas as pd


pd.set_option('display.max_rows', None)


def featureCounts(bam, gtf, summary, threads=12):
    p = subprocess.run(
        ["featureCounts",
            "-T", str(threads),
            "-L", "-s", "1",
            "-t", "exon", "-g", "gene_id", "-a", gtf,
            "-R", "CORE",
            "-o", summary, bam],
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT)
    
    return p.stdout


def pysam_open(alignment_file, in_format='BAM'):
    """Open SAM/BAM file using pysam.

    :param alignment_file: Input file.
    :param in_format: Format (SAM or BAM).
    :returns: pysam.AlignmentFile
    :rtype: pysam.AlignmentFile
    """
    if in_format == 'BAM':
        mode = "rb"
    elif in_format == 'SAM':
        mode = "r"
    else:
        raise Exception("Invalid format: {}".format(in_format))

    aln_iter = pysam.AlignmentFile(alignment_file, mode)
    return aln_iter


def align(query, pattern, max_ed, normalise=False):
        
    # https://github.com/nanoporetech/pipeline-umi-amplicon/blob/master/lib/umi_amplicon_tools/extract_umis.py
    seq = pattern
    for c in 'actgACTG':
        seq = seq.replace(c, "")
    wildcard = set(''.join(seq))

    equalities=[("M", "A"), ("M", "C"), ("R", "A"), ("R", "A"), ("W", "A"), ("W", "A"), ("S", "C"), ("S", "C"), ("Y", "C"), ("Y", "C"), ("K", "G"), ("K", "G"), ("V", "A"), ("V", "C"), ("V", "G"), ("H", "A"), ("H", "C"), ("H", "T"), ("D", "A"), ("D", "G"), ("D", "T"), ("B", "C"), ("B", "G"), ("B", "T"), ("N", "G"), ("N", "A"), ("N", "T"), ("N", "C")]
    
    result = edlib.align(
        pattern,
        query,
        task="path",
        mode="HW",
        k=max_ed,
        additionalEqualities=equalities,
    )
    if result["editDistance"] == -1:
        return None, None

    ed = result["editDistance"]
    if not normalise:
        locs = result["locations"][0]
        umi = query[locs[0]:locs[1]+1]
        return ed, umi

    # Extract and normalise UMI
    umi = ""
    align = edlib.getNiceAlignment(result, pattern, query)
    for q, t in zip(align["query_aligned"], align["target_aligned"]):
        if q not in wildcard:
            continue
        if t == "-":
            umi += "N"
        else:
            umi += t

    if len(umi) != 16:
        raise RuntimeError("UMI length incorrect: {}".format(umi))

    return ed, umi


# run vsearch
def cluster_umis(in_fasta, clusters_folder, centroids, consensus, threads=6):
    if not os.path.exists(clusters_folder):
        os.mkdir(clusters_folder)

    subprocess.run(
        ["vsearch", "--clusterout_id", "--clusters", clusters_folder, 
            "--threads", str(threads), "--cluster_fast", in_fasta, "--clusterout_sort",
            "--minseqlength", "16", "--maxseqlength", "16",
            "--centroids", centroids, "--consout", consensus,
            "--gapopen", "0E/5I", "--gapext", "0E/2I", "--mismatch", "-8", "--match", "6",
            "--iddef", "0", "--minwordmatches", "0", "-id", "0.85"],
        stdout=subprocess.DEVNULL,
        stderr=subprocess.STDOUT)


# append to fasta file
def write_fasta(header, seq, filename):
    with open(filename, 'a') as f:
        f.write('>' + header + "\n")
        f.write(seq + "\n")


# handle sys args
def parse_args(argv):
    """Read arguments from command line."""

    parser = argparse.ArgumentParser()

    parser.add_argument("-b", "--bam", type=str)
    parser.add_argument("-g", "--gtf", type=str)
    parser.add_argument("-o", "--output", type=str)
    parser.add_argument("-l", "--log", type=str, default = None)
    parser.add_argument("-p", "--pattern", type=str, default="TTTVVVVTTVVVVTTVVVVTTVVVVTTT")
    parser.add_argument("-c", "--cores", type=int, default=os.cpu_count() - 2,
        help="Number of threads used by featureCounts, vsearch and medaka.")

    args = parser.parse_args()

    return args


def main(argv=sys.argv[1:]):

    args = parse_args(argv=argv)

    gtffile = args.gtf
    bamfile = args.bam
    outfile = args.output
    logfile = args.log

    pattern = args.pattern

    threads = args.cores


    logging.basicConfig(
        format='%(asctime)s - %(levelname)s - %(message)s',
        level=logging.INFO
    )


    # make sure working dir exists
    os.makedirs(os.path.dirname(outfile), exist_ok=True)



    # tempdir
    tempfile.tempdir = os.path.dirname(outfile)
    tempfolder = tempfile.TemporaryDirectory()


    # assign reads to genes using featureCounts
    logging.info("Assigning reads to genes using featureCounts...")

    fcfile = os.path.join(tempfolder.name, os.path.basename(bamfile))
    os.makedirs(os.path.dirname(fcfile), exist_ok=True)
    fc_stdout = featureCounts(bamfile, gtffile, fcfile, threads)
    print(fc_stdout.decode())

    # create one folder per gene with a fa file holding all umis
    logging.info("Extracting umis and grouping by gene...")

    # count why reads are removed
    reads_total = 0
    good_reads = 0
    no_rx_tag = 0
    no_gene_attached = 0
    normalisation_failed = 0

    with open(fcfile + ".featureCounts", newline='') as fc, pysam_open(bamfile) as bam:

        fcin = csv.reader(fc, delimiter='\t')

        for (read, row) in tqdm(zip(bam, fcin)):

            reads_total += 1

            if read.query_name != row[0]:
                raise RuntimeError('Read names in bam and fc file are different: {} {} \n'.format(read.query_name, row[0]))

            try:
                umi = read.get_tag(tag="RX")
            
            except:
                no_rx_tag += 1
                continue

            gene_id = row[3]

            if umi != "None" and gene_id != "NA":
                ed, norm_umi = align(umi, pattern, 3, True)

                if norm_umi is None:
                    normalisation_failed += 1
                    continue

                good_reads += 1
            
                gene_folder = os.path.join(tempfolder.name, gene_id)

                if not os.path.exists(gene_folder):
                    os.mkdir(gene_folder)

                fasta_header = "|".join([read.query_name, read.get_tag(tag="TS"), read.seq]) 
                write_fasta(fasta_header, norm_umi, os.path.join(gene_folder, "norm_umi.fa"))
        
            elif gene_id == "NA":
                no_gene_attached += 1


    if logfile is None:
        print("Total reads analysed:", reads_total)
        print("Successfully recovered umi:", good_reads, "{0:.0%}".format(good_reads / reads_total))
        print("Reads without RX tag::", no_rx_tag, "{0:.0%}".format(no_rx_tag / reads_total))
        print("Reads not mapping to gene:", no_gene_attached, "{0:.0%}".format(no_gene_attached / reads_total))
        print("Failed to normalise umi:", normalisation_failed, "{0:.0%}".format(normalisation_failed / reads_total))

    else:
        with open(logfile, "a") as f:   
            print("Total reads analysed:", reads_total, file=f)
            print("Successfully recovered umi:", good_reads, "{0:.0%}".format(good_reads / reads_total), file=f)
            print("Reads without RX tag::", no_rx_tag, "{0:.0%}".format(no_rx_tag / reads_total), file=f)
            print("Reads not mapping to gene:", no_gene_attached, "{0:.0%}".format(no_gene_attached / reads_total), file=f)
            print("Failed to normalise umi:", normalisation_failed, "{0:.0%}".format(normalisation_failed / reads_total), file=f)


    # cluster umis
    logging.info("Clustering umis...")

    gene_folders = [f.path for f in os.scandir(tempfolder.name) if f.is_dir()]

    for gene_folder in tqdm(gene_folders):
        cluster_umis(
            os.path.join(gene_folder, "norm_umi.fa"),
            os.path.join(gene_folder, "vsearch_clusters/"),
            os.path.join(gene_folder, "clusters_centroids.fasta"),
            os.path.join(gene_folder, "clusters_consensus.fasta"),
            threads=threads)


    # write fa cluster files with actual sequences instead of umi seq
    logging.info("Write cluster fasta files with query sequences...")

    reads_per_cluster = []

    for gene_folder in tqdm(gene_folders):
        vsearch_folder = os.path.join(gene_folder, "vsearch_clusters/")
        cluster_folder = os.path.join(gene_folder, "clusters_fa/")

        if not os.path.exists(cluster_folder):
            os.mkdir(cluster_folder)

        for cluster in os.scandir(vsearch_folder):
            with pysam.FastxFile(cluster) as fh:
                out_file = os.path.join(cluster_folder, os.path.basename(gene_folder) + "_" + os.path.basename(cluster) + ".fa")

                reads = 0

                for entry in fh:

                    reads += 1

                    cols = entry.name.split("|")

                    if len(cols) != 3:
                        raise RuntimeError("Fasta header does not have 4 entries")

                    read_id = '|'.join([cols[0], cols[1], entry.sequence])
                    seq = cols[2]

                    write_fasta(read_id, seq, out_file)
                
                reads_per_cluster.append(reads)


    num_clusters = len(reads_per_cluster)
    num_reads = sum(reads_per_cluster)
    umi_distribution = Counter(reads_per_cluster)
    umi_df = pd.DataFrame.from_dict(umi_distribution, orient='index').reset_index()
    umi_df = umi_df.rename(columns={'index':'readsPerCluster', 0:'count'})

    if logfile is None:
        print("Clustered", num_reads, "reads into", num_clusters, "clusters, removing", "{0:.0%}".format((num_reads - num_clusters) / num_reads), "of the reads due to PCR duplication")
        print(umi_df.sort_values("readsPerCluster"))
    else:
        with open(logfile, "a") as f:
            print("Clustered", num_reads, "reads into", num_clusters, "clusters, removing", "{0:.0%}".format((num_reads - num_clusters) / num_reads), "of the reads due to PCR duplication", file=f)
            print(umi_df.sort_values("readsPerCluster"), file=f)

    logging.info("Extracting longest sequences...")
    for gene_folder in tqdm(gene_folders):
        cluster_folder = os.path.join(gene_folder, "clusters_fa/")
        consensus_file = os.path.join(gene_folder, "longest.fa")

        clusters = [f.path for f in os.scandir(cluster_folder)]

        for cluster in clusters:
            id = ""
            seq = ""

            with pysam.FastxFile(cluster) as fh:
                for entry in fh:
                    if len(entry.sequence) > len(seq):
                        id = entry.name
                        seq = entry.sequence
                       
            write_fasta(id, seq, consensus_file)
                        

    # concatenate
    logging.info("Concatenating consensus files...")
    consensus_files = [os.path.join(g, "longest.fa") for g in gene_folders]

    with open(outfile,'wb') as wfd:
        for f in consensus_files:
            with open(f, 'rb') as fd:
                shutil.copyfileobj(fd, wfd)


    logging.info('Done.')
    return 0


if __name__ == "__main__":
    exit(main())
