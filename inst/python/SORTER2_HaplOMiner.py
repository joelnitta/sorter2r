import sys
import subprocess
import csv
import argparse
import os
import shutil


parser = argparse.ArgumentParser()
parser.add_argument("-wd", "--workingdir")
parser.add_argument("-outdir", "--outputdir")
parser.add_argument("-cpref", "--cpref")
parser.add_argument("-c", "--coverage", type=float)
parser.add_argument("-d", "--depth", type=float)
# Optional: directory containing raw FASTQ files when Stage 1A was run
# without Trim Galore (trim=F).  When omitted, reads are expected inside
# each *_assembly/ subdirectory as *_R1_val_1.fq / *_R2_val_2.fq.
parser.add_argument("-reads", "--readsdir", default=None)
parser.add_argument(
    "-v", "--verbose", action="store_true", default=False,
    help="Print debug information during processing"
)
args = parser.parse_args()

quiet = {} if args.verbose else {
    "stdout": subprocess.DEVNULL, "stderr": subprocess.DEVNULL
}

workingdir = (
    args.workingdir if args.workingdir.endswith('/')
    else args.workingdir + '/'
)
outdir = (
    args.outputdir if args.outputdir.endswith('/')
    else args.outputdir + '/'
)
os.makedirs(outdir, exist_ok=True)

readsdir = None
if args.readsdir is not None:
    readsdir = (
        args.readsdir if args.readsdir.endswith('/')
        else args.readsdir + '/'
    )

# All HaplOMiner writes go under outdir — never into workingdir.
cpdir = outdir + 'all_chloroplasts/'
workfilesdir = outdir + 'haplominer_workfiles/'
os.makedirs(cpdir, exist_ok=True)
os.makedirs(workfilesdir, exist_ok=True)

cpref = args.cpref

print('HaplOMiner: indexing organellar reference')
subprocess.call(["bwa index %s" % cpref], shell=True, **quiet)
subprocess.call(["samtools faidx %s" % cpref], shell=True, **quiet)

# Map reads and generate per-sample consensus.
for folder in sorted(os.listdir(workingdir)):
    if not folder.endswith('_assembly'):
        continue
    asm_folder = workingdir + folder + '/'
    sample_base = folder[:-8]  # e.g. 'Iimura12_cthysanostomum_'
    raw_base = sample_base.rstrip('_')  # e.g. 'Iimura12_cthysanostomum'
    wf_folder = workfilesdir + folder + '/'
    os.makedirs(wf_folder, exist_ok=True)
    prefix = wf_folder + sample_base

    # Find reads: prefer Trim Galore output, fall back to raw FASTQ.
    trimmed_R1 = asm_folder + sample_base + 'R1_val_1.fq'
    trimmed_R2 = asm_folder + sample_base + 'R2_val_2.fq'
    if os.path.exists(trimmed_R1):
        read_path = trimmed_R1
        R2_path = trimmed_R2
    elif readsdir is not None:
        read_path = readsdir + raw_base + '_R1.fastq'
        R2_path = readsdir + raw_base + '_R2.fastq'
        if not os.path.exists(read_path):
            raise RuntimeError(
                'reads not found for %s. Checked %s and %s'
                % (folder, trimmed_R1, read_path)
            )
    else:
        raise RuntimeError(
            'reads not found: %s. Pass -reads to specify a raw FASTQ '
            'directory when Stage 1A was run without Trim Galore.'
            % trimmed_R1
        )

    print('HaplOMiner: mapping %s to organellar reference' % raw_base)
    if args.verbose:
        print('  R1: %s' % read_path)
        print('  R2: %s' % R2_path)

    subprocess.call(
        ["bwa mem -V %s %s %s > %smapreads.sam" % (
            cpref, read_path, R2_path, prefix)],
        shell=True, **quiet
    )
    subprocess.call(
        ["samtools sort %smapreads.sam -o %smapreads.bam" % (
            prefix, prefix)],
        shell=True, **quiet
    )
    subprocess.call(
        ["samtools index %smapreads.bam" % prefix],
        shell=True, **quiet
    )

    stat_path = wf_folder + raw_base + '_readstats.txt'
    subprocess.call(
        ["samtools flagstat %smapreads.bam > %s" % (prefix, stat_path)],
        shell=True, **quiet
    )
    subprocess.call(
        ["samtools depth -a %smapreads.bam "
         "| awk '{c++;s+=$3}END{if(c>0)print s/c; else print 0}'"
         " >> %s" % (prefix, stat_path)],
        shell=True
    )
    subprocess.call(
        ["samtools depth -a %smapreads.bam "
         "| awk '{c++; if($3>0) total+=1}"
         "END{if(c>0)print (total/c)*100; else print 0}'"
         " >> %s" % (prefix, stat_path)],
        shell=True
    )

    # Count mapped reads before attempting consensus generation.
    mapped_result = subprocess.run(
        ['samtools', 'view', '-c', '-F', '4',
         prefix + 'mapreads.bam'],
        capture_output=True, text=True
    )
    mapped_count = int(
        mapped_result.stdout.strip()
    ) if mapped_result.stdout.strip() else 0

    if mapped_count > 0:
        raw_consensus = wf_folder + raw_base + '_raw_consensus.fasta'
        subprocess.call(
            ["samtools consensus -f fasta -o %s %smapreads.bam" % (
                raw_consensus, prefix)],
            shell=True, **quiet
        )
        # Rename consensus headers to sample name.
        cp_fasta = cpdir + raw_base + '_chloroplast.fasta'
        seq_idx = 0
        with open(raw_consensus, 'r') as fin, \
                open(cp_fasta, 'w') as fout:
            for line in fin:
                if line.startswith('>'):
                    seq_idx += 1
                    if seq_idx > 1:
                        fout.write('>%s_%d\n' % (raw_base, seq_idx))
                    else:
                        fout.write('>%s\n' % raw_base)
                else:
                    fout.write(line)
        print(
            'HaplOMiner: consensus written for %s (%d mapped reads)'
            % (raw_base, mapped_count)
        )
    else:
        print(
            'HaplOMiner: no reads mapped for %s, skipping consensus'
            % raw_base
        )

    if os.path.exists(prefix + 'mapreads.sam'):
        os.remove(prefix + 'mapreads.sam')

# Compile read statistics into readstats_cp.csv.
STATDICT = {}

for folder in os.listdir(workingdir):
    if folder.endswith('_assembly'):
        raw_base = folder[:-8].rstrip('_')
        if raw_base not in STATDICT:
            STATDICT[raw_base] = {}

for folder in os.listdir(workingdir):
    if not folder.endswith('_assembly'):
        continue
    raw_base = folder[:-8].rstrip('_')
    wf_folder = workfilesdir + folder + '/'
    stat_path = wf_folder + raw_base + '_readstats.txt'
    if not os.path.exists(stat_path):
        continue
    with open(stat_path, 'r') as statfile:
        for position, line in enumerate(statfile):
            if '+' in line:
                parts = line.split(' ')
                statlabel = parts[3].strip('\n')
                statint = parts[0].strip('\n')
                STATDICT[raw_base][statlabel] = []
                STATDICT[raw_base][statlabel].append(int(statint))
            else:
                if position == 16:
                    statint = line.strip('\n')
                    val = int(float(statint)) if statint else 0
                    STATDICT[raw_base]['readdepth'] = [val]
                elif position == 17:
                    statint = line.strip('\n')
                    val = int(float(statint)) if statint else 0
                    STATDICT[raw_base]['coverage'] = [val]

if STATDICT:
    het_values_list = list(STATDICT.values())
    columns = list(het_values_list[0].keys())
    header = ['Individual'] + columns
    tmp_csv = outdir + 'readstats_cp_tmp.csv'

    with open(tmp_csv, 'w') as outfile:
        writer = csv.writer(outfile)
        writer.writerow(header)
        for ind in STATDICT.keys():
            writer.writerow(
                [ind] + [STATDICT[ind].get(col, '') for col in columns]
            )

    with open(tmp_csv, 'r') as infile:
        with open(outdir + 'readstats_cp.csv', 'w') as outfile:
            data = infile.read()
            data = data.replace(']', '')
            data = data.replace('[', '')
            outfile.write(data)

    os.remove(tmp_csv)

# Filter samples by coverage and depth thresholds.
min_coverage = float(args.coverage)
min_depth = float(args.depth)
cov_str = str(int(min_coverage))
dep_str = str(int(min_depth))
filtdir = outdir + 'coverage%s_depth%s_filtered/' % (cov_str, dep_str)
os.makedirs(filtdir, exist_ok=True)

passing_samples = []
for ind in STATDICT:
    cov_val = STATDICT[ind].get('coverage', [0])[0]
    dep_val = STATDICT[ind].get('readdepth', [0])[0]
    if cov_val >= min_coverage and dep_val >= min_depth:
        passing_samples.append(ind)
        cp_fasta = cpdir + ind + '_chloroplast.fasta'
        if os.path.exists(cp_fasta):
            shutil.copy(cp_fasta, filtdir)

print(
    'HaplOMiner: %d sample(s) passed coverage >= %s%% and depth >= %s: %s'
    % (len(passing_samples), cov_str, dep_str,
       ', '.join(passing_samples) if passing_samples else 'none')
)

# Concatenate filtered samples and align.
if passing_samples:
    all_fasta = filtdir + 'all_filtered.fasta'
    with open(all_fasta, 'w') as outf:
        for ind in passing_samples:
            cp_fasta = filtdir + ind + '_chloroplast.fasta'
            if os.path.exists(cp_fasta):
                with open(cp_fasta, 'r') as inf:
                    outf.write(inf.read())

    al_fasta = (
        outdir
        + 'HaplOMiner_coverage%s_depth%s_filtered_al.fasta'
        % (cov_str, dep_str)
    )
    print('HaplOMiner: aligning filtered sequences with mafft')
    subprocess.call(
        ['mafft --auto %s > %s' % (all_fasta, al_fasta)],
        shell=True, **quiet
    )

print('HaplOMiner complete. Results in: %s' % outdir)
