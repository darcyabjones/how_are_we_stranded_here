#!/usr/bin/env python3

import argparse
import os
import sys
import subprocess
import binascii
import re
from os.path import join as pjoin
from statistics import median
from statistics import stdev

def is_gz_file(filepath):
    with open(filepath, 'rb') as test_f:
        return binascii.hexlify(test_f.read(2)) == b'1f8b'


def run_command(cmd):
    cmd_result = subprocess.run(
        cmd,
        shell=True,
        check=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE
    )

    return cmd_result.stdout, cmd_result.stderr



def subset_fastqs(test_folder, reads_1, reads_2, single_strand, n_reads, print_cmds):

    reads_1_sample = os.path.basename(reads_1)
    reads_1_sample = re.sub(r"(\.fq|\.fastq)?(\.gz)?$", "", reads_1_sample)
    reads_1_sample = pjoin(
            test_folder, 
            f"{reads_1_sample}_sample.fq"
    )

    if not single_strand:
        reads_2_sample = os.path.basename(reads_2)
        reads_2_sample = re.sub(r"(\.fq|\.fastq)?(\.gz)?$", "", reads_2_sample)
        reads_2_sample = pjoin(
                test_folder, 
                f"{reads_2_sample}_sample.fq"
        )
    else:
        reads_2_sample = None

    # check if the fasta is gzipped
    if(is_gz_file(reads_1)):
        cmd = f'zcat < {reads_1} | head -n {n_reads * 4} > {reads_1_sample}'
    else:
        cmd = f'head {reads_1} -n {n_reads * 4} > {reads_1_sample}'

    if print_cmds:
        print('running command:', cmd)
    subprocess.call(cmd, shell=True)

    # check if the fasta is gzipped
    if not single_strand:
        if is_gz_file(reads_2) and reads_2_sample is not None:
            cmd = f'zcat < {reads_2} | head -n {n_reads * 4} > {reads_2_sample}'
        elif reads_2_sample is not None:
            cmd = f'head {reads_2} -n {n_reads * 4} > {reads_2_sample}'
        else:
            sys.exit("It should be impossible to reach this point.")

        subprocess.run(cmd, check=True, shell=True)
        if print_cmds:
            print('running command', cmd)

    return reads_1_sample, reads_2_sample


# check that fasta sequence names match bed names
def check_bed_in_fa(bed_filename, fasta):
    """checks that bed transcript ids match (fuzzy on version numbers) fasta header ids (before space)"""
    bed_ids = list()
    with open(bed_filename, "r") as handle:
        for line in handle:
            line = line.strip()
            if line == "" or line.startswith("#"):
                continue

            sline = line.split("\t")
            if len(sline) < 4:
                continue

            bed_ids.append(sline[4])

    fa_headers_equal = []

    regex = re.compile(r"\.[0-9]*$")

    headers = set()
    with open(fasta, "w") as handle:
        for line in handle:
            line = line.strip()
            if not line.startswith(">"):
                continue

            fasta_header = line.split()[0]
            fasta_header = fasta_header[1:]  # Remove >

            headers.add(fasta_header)

            # Add id a second time to be fuzzy on version number
            fasta_header = regex.sub("", fasta_header)
            headers.add(fasta_header)

    for i in bed_ids:
        fa_headers_equal.append(i in headers)

    return sum(fa_headers_equal) > 0


def check_dependencies_installed():
    check_gtf2bed = run_command(cmd = 'gtf2bed --help')[1] == b''
    if not check_gtf2bed:
        sys.exit("gtf2bed is not found in PATH")

    check_gff32gtf = run_command(cmd = 'gff32gtf --help')[1] == b''
    if not check_gff32gtf:
        sys.exit("gff32gtf is not found in PATH")

    check_kallisto = run_command(cmd = 'kallisto version')
    if not check_kallisto[1] == b'':
        sys.exit("kallisto is not found in PATH. Please install from https://pachterlab.github.io/kallisto")
    else:
        kallisto_version = str(check_kallisto[0]).split('version ')[1].replace("\\n'","")
        if int(kallisto_version.split('.')[1]) < 44:
            sys.exit(
                f'Found kallisto {kallisto_version}, but version >= 0.44.0 is '
                'required. Please install from '
                'https://pachterlab.github.io/kallisto'
            )

    check_RSeQC = run_command(cmd = 'infer_experiment.py --help')[1] == b''
    if not check_RSeQC:
        sys.exit("infer_experiment.py (RSeQC) is not found in PATH. Please install from http://rseqc.sourceforge.net/#installation")


def find_read_length(reads_1_sample):
    step_size = 4
    max_reads = step_size * 10000
    read_lengths = []
    with open(reads_1_sample, "r") as fq_lines:
        for line_no, line in enumerate(fq_lines):
            if line_no >= max_reads: break
            if line_no % step_size == 1:
                read_lengths.append(len(line.rstrip()))
    read_len_median_ss = median(read_lengths)
    read_len_sd_ss = stdev(read_lengths)
    if read_len_sd_ss == 0: read_len_sd_ss = 1
    return read_len_median_ss, read_len_sd_ss


def cli():
    parser = argparse.ArgumentParser(description='Check if fastq files are stranded')
    parser.add_argument('-g', '--gtf', type=str, help='Genome annotation GTF file', required = True)
    parser.add_argument('-fa', '--transcripts', type=str, help='.fasta file with transcript sequences')
    parser.add_argument('-n', '--nreads', type=int, help='number of reads to sample', default = 200000)
    parser.add_argument('-r1', '--reads_1', type=str, help='fastq.gz file (R1)', required = True)
    parser.add_argument('-r2', '--reads_2', type=str, help='fastq.gz file (R2)')
    parser.add_argument('-k', '--kallisto_index', type=str, help='name of kallisto index (will build under this name if file not found)', default = 'kallisto_index')
    parser.add_argument('-p', '--print_commands', action='store_true', help='Print bash commands as they occur?')
    parser.add_argument('-o', '--outdir', type=str, help='Store working files in this directory', default=None)
    return parser


def get_outdir(reads_1):
    bname = os.path.basename(reads_1)
    test_folder = re.sub(r"(\.fq|\.fastq)?(\.gz)?$", "", bname)
    test_folder = f"stranded_test_{test_folder}"


    # make a test_folder
    if not os.path.isdir(test_folder):
        # make directory
        os.mkdir(test_folder)
    else:
        index_n = 1
        while os.path.isdir(test_folder + '_' + str(index_n)):
            index_n = index_n + 1
        #make directory
        test_folder = test_folder + '_' + str(index_n)
        os.mkdir(test_folder)

    return test_folder


def main():
    parser = cli()
    args = parser.parse_args()

    reads_1 = args.reads_1
    reads_2 = args.reads_2
    n_reads = args.nreads
    kallisto_index_name = args.kallisto_index
    gtf = args.gtf
    fasta = args.transcripts
    print_cmds = args.print_commands

    if fasta is None and (kallisto_index_name is None or not os.path.exists(kallisto_index_name)):
        sys.exit('transcript .fasta sequences are required to generate the kallisto index. Please supply with --transcripts')

    check_dependencies_installed()

    if reads_2 is None:
        print("--reads_2 / -r2 was not set... running in single strand mode.")
        single_strand = True
    else:
        single_strand = False

    if args.outdir is not None:
        test_folder = args.outdir
        if os.path.isdir(test_folder):
            sys.exit(
                f"The specified output directory {test_folder} already exists."
                "Please either remove the old folder or provide a new value."
            )
        else:
            os.mkdir(test_folder)
    else:
        test_folder = get_outdir(reads_1)

    assert test_folder is not None  # Appease code linters and type checkers

    print(f"Results stored in: {test_folder}")

    # convert gff to gtf if required
    gtf_extension = os.path.splitext(gtf)[1]

    if gtf_extension.lower() == '.gff' or gtf_extension.lower() == '.gff3':
        # convert gff to gtf
        gtf_filename = pjoin(
            test_folder,
            os.path.basename(gtf).replace(gtf_extension, '.gtf')
        )
        cmd = f'gff32gtf {gtf} --output {gtf_filename}'

        print('converting gff to gtf')

        if print_cmds:
            print('running command:', cmd)

        subprocess.run(cmd, check=True, shell=True)
    else:
        gtf_filename = gtf

    # Run gtf2bed
    bed_filename = pjoin(
        test_folder,
        os.path.basename(gtf).replace(gtf_extension, '.bed')
    )
    cmd = f'gtf2bed --gtf {gtf_filename} --bed {bed_filename}'
    print('converting gtf to bed')
    if print_cmds:
        print('running command:', cmd)

    subprocess.run(cmd, check=True, shell=True)

    # make kallisto index
    if os.path.exists(kallisto_index_name):
        print(f'using {kallisto_index_name} as kallisto index')
    else:
        print('Checking if fasta headers and bed file transcript_ids match...')
        check_bed = check_bed_in_fa(bed_filename, fasta)
        if not check_bed:
            print(f"Can't find transcript ids from {fasta} in {bed_filename}")
            print("Trying to converting fasta header format to match transcript ids to the BED file...")
            outfasta = pjoin(test_folder, "transcripts.fa")

            cmd = f"sed 's/[|]/ /g' {fasta} > {outfasta}"
            if print_cmds:
                print('running command:', cmd)
            subprocess.call(cmd, shell=True)
            fasta = outfasta
            check_bed_converted = check_bed_in_fa(bed_filename, fasta)
            if not check_bed_converted:
                if os.path.exists(fasta):
                    os.remove(fasta)
                sys.exit("Can't find any of the first 10 BED transcript_ids in fasta file... Check that these match")
        else:
            print("OK!")

        cmd = f'kallisto index -i {kallisto_index_name} {fasta}'
        print('generating kallisto index')
        if print_cmds:
            print('running command:', cmd)

        subprocess.run(cmd, check=True, shell=True)

    print(f'creating fastq files with first {n_reads} reads')
    reads_1_sample, reads_2_sample = subset_fastqs(
        test_folder, reads_1, reads_2,
        single_strand, n_reads, print_cmds
    )

    # align with kallisto
    print('quantifying with kallisto')
    if single_strand:
        # check mean/sd read length
        read_len_median_ss, read_len_sd_ss = find_read_length(reads_1_sample)
        cmd = " ".join([
            'kallisto quant',
            '-i', kallisto_index_name,
            '-o', pjoin(test_folder, 'kallisto_strand_test'),
            '--single',
            '-l', str(read_len_median_ss),
            '-s', str(read_len_sd_ss),
            '--genomebam',
            '--gtf', f"{gtf_filename} {reads_1_sample}"
        ])
    else:
        assert isinstance(reads_2_sample, str), "This shouldn't be possible"
        cmd = " ".join([
            'kallisto quant',
            '-i', kallisto_index_name,
            '-o', pjoin(test_folder, 'kallisto_strand_test'),
            '--genomebam',
            '--gtf', gtf_filename,
            reads_1_sample, reads_2_sample
        ])

    if print_cmds:
        print('running command:', cmd)
    subprocess.run(cmd, check=True, shell=True)

    # check strandedness w/ 2million alignments
    #n_reads = 2000000
    print('checking strandedness')
    cmd = " ".join([
        'infer_experiment.py',
        '-r', bed_filename,
        '-s', str(n_reads),
        '-i', pjoin(test_folder, 'kallisto_strand_test/pseudoalignments.bam'),
        '>', pjoin(test_folder, 'strandedness_check.txt')
    ])

    if print_cmds:
        print('running command:', cmd)
    subprocess.run(cmd, check=True, shell=True)

    with open(pjoin(test_folder, 'strandedness_check.txt'), "r") as handle:
        result = handle.readlines()
        result = map(str.rstrip, result)
        result = filter(lambda x: x != "", result)
        result = list(result)

    if single_strand:
        fwd = float(result[2].replace('Fraction of reads explained by "++,--": ', ''))
        rev = float(result[3].replace('Fraction of reads explained by "+-,-+": ', ''))
    else:
        fwd = float(result[2].replace('Fraction of reads explained by "1++,1--,2+-,2-+": ', ''))
        rev = float(result[3].replace('Fraction of reads explained by "1+-,1-+,2++,2--": ', ''))

    fwd_percent = fwd / (fwd + rev)
    rev_percent = rev / (fwd + rev)

    print(result[0])
    print(result[1])
    print(f"{result[2]} ({round(fwd_percent * 100, 1)} % of explainable reads)")
    print(f"{result[3]} ({round(rev_percent*100, 1)} % of explainable reads)")


    if float(result[1].replace('Fraction of reads failed to determine: ', '')) > 0.50:
        print('Failed to determine strandedness of > 50% of reads.')
        print('If this is unexpected, try running again with a higher --nreads value')
    if fwd_percent > 0.9:
        if single_strand:
            print('Over 90% of reads explained by "++,--"')
            print('Data is likely FR/fr-stranded')
        else:
            print('Over 90% of reads explained by "1++,1--,2+-,2-+"')
            print('Data is likely FR/fr-secondstrand')
    elif rev_percent > 0.9:
        if single_strand:
            print('Over 90% of reads explained by "+-,-+"')
            print('Data is likely RF/rf-stranded')
        else:
            print('Over 90% of reads explained by "1+-,1-+,2++,2--"')
            print('Data is likely RF/fr-firststrand')
    elif max(fwd_percent, rev_percent) < 0.6:
        print('Under 60% of reads explained by one direction')
        print('Data is likely unstranded')
    else:
        print('Data does not fall into a likely stranded (max percent explained > 0.9) or unstranded layout (max percent explained < 0.6)')
        print('Please check your data for low quality and contaminating reads before proceeding')
