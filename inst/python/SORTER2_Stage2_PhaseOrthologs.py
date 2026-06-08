import sys
import subprocess
import fileinput
import csv
import itertools
import argparse
import re
import glob
import os
from os import path
import shutil
from shutil import move
from shutil import copyfile
import Bio
from Bio import SeqIO
from concurrent.futures import ProcessPoolExecutor


def run_usearch(cmd):
    cwd = os.getcwd()
    if shutil.which('usearch'):
        subprocess.call('usearch ' + cmd, shell=True, **quiet)
    elif shutil.which('docker'):
        docker_cmd = (
            'docker run --rm -v %s:%s -w %s '
            'joelnitta/usearch:latest %s' % (cwd, cwd, cwd, cmd)
        )
        subprocess.call(docker_cmd, shell=True, **quiet)
    else:
        sys.exit(
            'Error: usearch not found. Install usearch locally or install '
            'Docker and build the joelnitta/usearch image.'
        )


parser = argparse.ArgumentParser()
parser.add_argument("-wa", "--assemblydir")
parser.add_argument("-wc", "--clusterdir")
parser.add_argument("-outdir", "--outputdir")
parser.add_argument("-pq", "--phasequal")
parser.add_argument("-al", "--aliter")
parser.add_argument("-indel", "--indelrep")
parser.add_argument("-idformat", "--idformat")
# Optional: directory containing raw FASTQ files when Stage 1A was run
# without Trim Galore (trim=F).  When omitted, reads are expected inside
# each *_assembly/ subdirectory as *_R1_val_1.fq / *_R2_val_2.fq.
parser.add_argument("-reads", "--readsdir", default=None)
parser.add_argument(
    "-t", "--threads", type=int, default=1,
    help="Number of parallel worker processes for per-sample phasing"
)
parser.add_argument(
    "-v", "--verbose", action="store_true", default=False,
    help="Print debug information during processing"
)
args = parser.parse_args()
quiet = {} if args.verbose else {
    "stdout": subprocess.DEVNULL, "stderr": subprocess.DEVNULL
}

assemblydir = args.assemblydir if args.assemblydir.endswith('/') else args.assemblydir + '/'
clusterdir = args.clusterdir if args.clusterdir.endswith('/') else args.clusterdir + '/'
readsdir = None
if args.readsdir is not None:
    readsdir = args.readsdir if args.readsdir.endswith('/') else args.readsdir + '/'
outdir = args.outputdir if args.outputdir.endswith('/') else args.outputdir + '/'
os.makedirs(outdir, exist_ok=True)

diploidclusters = clusterdir + 'diploidclusters/'
contigdir = clusterdir + 'diploids/'

# All Stage 2 writes go under outdir — never into assemblydir or clusterdir.
workfilesdir = outdir + 'assembly_workfiles/'
os.makedirs(workfilesdir, exist_ok=True)

# Intermediate diploidcluster files (degap, deinterleaved, udb) written here.
workdipclusters = outdir + 'diploidclusters_work/'
os.makedirs(workdipclusters, exist_ok=True)


#Make diploids_phased folder if needed
def check_folder(directory, folder_name):
	folder_path = os.path.join(directory, folder_name)
	if not os.path.exists(folder_path):
		os.makedirs(folder_path)
		if args.verbose:
			print("Created folder: " + folder_path)
	else:
		if args.verbose:
			print("Folder already exists: " + folder_path)

folders_to_check = ["diploids_phased"]

for folder_name in folders_to_check:
	check_folder(outdir, folder_name)

#Set Directory Variables
phasedir = outdir + 'diploids_phased/'
direc = os.listdir(assemblydir)
map_contigs_to_baits_dir = sorted(os.listdir(contigdir))

def remove_if_exists(path):
    try:
        os.remove(path)
    except FileNotFoundError:
        pass

#Define function to change sequence IDs
def replaceAll(file, searchExp, replaceExp):
	for line in fileinput.input(file, inplace=1):
		if searchExp in line:
			line = line.replace(searchExp, replaceExp)
		sys.stdout.write(line)

#Functions to get longest contigs in fasta file
def get_longest_sequences(input_file, output_file, X):
	X = int(X)
	sequences = list(SeqIO.parse(input_file, "fasta"))
	longest_sequences = sorted(sequences, key=lambda seq: len(seq.seq), reverse=True)[:X]
	SeqIO.write(longest_sequences, output_file, "fasta")

def extract_longest_sequences(input_file, X):
	X = int(X)
	output_file = input_file + "longest.fa"
	get_longest_sequences(input_file, output_file, X)

#Function to remove duplicate sequence IDs and keep the longest sequence
def remove_dup(input_file, output_file):
	sequences = list(SeqIO.parse(input_file, "fasta"))
	sequence_dict = {}
	for seq in sequences:
		seq_id = seq.id
		if seq_id not in sequence_dict:
			sequence_dict[seq_id] = seq
		else:
			if len(seq.seq) > len(sequence_dict[seq_id].seq):
				sequence_dict[seq_id] = seq
	SeqIO.write(sequence_dict.values(), output_file, "fasta")

#Function to deinterleave fasta files
def deinterleave_fasta(input_file, output_file):
	with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
		sequences = {}
		current_id = None
		for line in infile:
			line = line.strip()
			if line.startswith(">"):
				current_id = line
				if current_id not in sequences:
					sequences[current_id] = []
			else:
				if current_id:
					sequences[current_id].append(line)
		for seq_id, seq_lines in sequences.items():
			outfile.write(seq_id + "\n")
			outfile.write("".join(seq_lines) + "\n")


#Degap locus-clusters; read from diploidclusters (Stage 1B, read-only),
#write to workdipclusters
print('Degapping locus-clusters for phasing')

for file in os.listdir(diploidclusters):
	if file.endswith('_'):
		newfilename = workdipclusters + file + 'degap.fasta'
		inpath = diploidclusters + file
		if args.verbose:
			print('degapping ' + file)
		with open(newfilename, "w") as o:
			for record in SeqIO.parse(inpath, "fasta"):
				record.seq = record.seq.replace("-", "")
				SeqIO.write(record, o, "fasta")

#deinterleave locus-clusters
print('deinterleaving locus-clusters for phasing')
for file in os.listdir(workdipclusters):
	if file.endswith('degap.fasta'):
		src = workdipclusters + file
		dst = workdipclusters + file[:-11] + 'deinterleaved_degap.fasta'
		if args.verbose:
			print("deinterleaving " + file)
		deinterleave_fasta(src, dst)
		os.remove(src)

#add new line in cluster files
for file in os.listdir(workdipclusters):
	if 'deinterleaved' in file:
		with open(workdipclusters + file, 'a+') as cluster:
			cluster.write('\n')

#make ublast data
os.chdir(workdipclusters)
subprocess.call("cat *degap.fasta > ALLsamples_allcontigs_allbaits_clusterannotated.fasta", shell=True, **quiet)
run_usearch("-makeudb_usearch ALLsamples_allcontigs_allbaits_clusterannotated.fasta -output diploid_master.udb")

print('Preparing Phase Directories...')

# Recompile sample-specific clusterbaits into per-sample annotated FASTAs,
# writing into workfilesdir instead of assemblydir.
for file in os.listdir(workdipclusters):
	if 'deinterleaved' in file:
		with open(workdipclusters + file, 'r') as allsamplefile:
			for line in allsamplefile:
				if '>' in line:
					linspl = line.split('_')
					folder_name = (linspl[2] + '_' + linspl[3].strip('\n')
						+ '_assembly')
					sample = (linspl[2] + '_' + linspl[3].strip('\n')
						+ '_allcontigs_allclusterbaits_annotated.fasta')
					if args.verbose:
						print(linspl)
					wf_folder = workfilesdir + folder_name + '/'
					os.makedirs(wf_folder, exist_ok=True)
					with open(wf_folder + sample, 'a+') as idx:
						if not line.endswith('\n'):
							idx.write('\n')
						else:
							while True:
								try:
									idx.write(line)
									seq = next(allsamplefile)
									idx.write(seq)
									if args.verbose:
										print(line)
									break
								except StopIteration as e:
									if args.verbose:
										print(e)
									break

os.chdir(assemblydir)


def _map_sample(folder):
    """Per-sample read mapping, phasing, and variant calling worker."""
    wf_folder = workfilesdir + folder + '/'
    asm_folder = assemblydir + folder + '/'
    for baits in os.listdir(wf_folder):
        if baits.endswith('allcontigs_allclusterbaits_annotated.fasta'):
            baits_path = wf_folder + baits
            sample_base = folder[:-8]  # e.g. 'Iimura18_cthysanostomum_'
            # Prefer Trim Galore output inside the assembly dir;
            # fall back to raw FASTQ files in readsdir when trim=F.
            trimmed_R1 = asm_folder + sample_base + 'R1_val_1.fq'
            trimmed_R2 = asm_folder + sample_base + 'R2_val_2.fq'
            if os.path.exists(trimmed_R1):
                read_path = trimmed_R1
                R2_path = trimmed_R2
            elif readsdir is not None:
                raw_base = sample_base.rstrip('_')
                read_path = readsdir + raw_base + '_R1.fastq'
                R2_path = readsdir + raw_base + '_R2.fastq'
                if not os.path.exists(read_path):
                    raise RuntimeError('reads not found for %s. '
                        'Checked %s and %s' % (folder, trimmed_R1, read_path))
            else:
                raise RuntimeError('reads not found: %s. '
                    'Pass -reads to specify a raw FASTQ directory '
                    'when Stage 1A was run without Trim Galore.' % trimmed_R1)
            prefix = wf_folder + sample_base
            if args.verbose:
                print(read_path)
                print(R2_path)
            subprocess.call(["bwa index %s" % (baits_path)], shell=True, **quiet)
            subprocess.call(["bwa mem -V %s %s %s > %smapreads.sam" % (
                baits_path, read_path, R2_path, prefix)], shell=True, **quiet)
            subprocess.call(["samtools sort %smapreads.sam -o %smapreads.bam" % (
                prefix, prefix)], shell=True, **quiet)
            subprocess.call(["samtools index %s" % (prefix + 'mapreads.bam')], shell=True, **quiet)
            subprocess.call(["samtools phase -A -Q %s -b %s %s" % (
                args.phasequal, prefix, prefix + 'mapreads.bam')], shell=True, **quiet)
            subprocess.call(["samtools sort %s -o %s" % (
                prefix + '.0.bam', prefix + '0srt.bam')], shell=True, **quiet)
            subprocess.call(["samtools sort %s -o %s" % (
                prefix + '.1.bam', prefix + '1srt.bam')], shell=True, **quiet)
            subprocess.call(["samtools sort %s -o %s" % (
                prefix + '.chimera.bam', prefix + 'chimerasrt.bam')], shell=True, **quiet)
            subprocess.call([
                "bcftools mpileup -B --min-BQ 20 -Ou -d 500 -f %s %s "
                "| bcftools call -mv --ploidy 1 -Oz -p .0001 -o %s" % (
                baits_path, prefix + '0srt.bam', prefix + '0.vcf.gz')], shell=True, **quiet)
            subprocess.call([
                "bcftools mpileup -B --min-BQ 20 -Ou -d 500 -f %s %s "
                "| bcftools call -mv --ploidy 1 -Oz -p .0001 -o %s" % (
                baits_path, prefix + '1srt.bam', prefix + '1.vcf.gz')], shell=True, **quiet)
            subprocess.call([
                "bcftools mpileup -B --min-BQ 20 -Ou -d 500 -f %s %s "
                "| bcftools call -mv --ploidy 1 -Oz -p .0001 -o %s" % (
                baits_path, prefix + 'chimerasrt.bam', prefix + 'chimera.vcf.gz')], shell=True, **quiet)
            subprocess.call(["bcftools index %s" % (prefix + '0.vcf.gz')], shell=True, **quiet)
            subprocess.call(["bcftools index %s" % (prefix + '1.vcf.gz')], shell=True, **quiet)
            subprocess.call(["bcftools index %s" % (prefix + 'chimera.vcf.gz')], shell=True, **quiet)
            subprocess.call(["cat %s | bcftools consensus %s > %s0_Final.fasta" % (
                baits_path, prefix + '0.vcf.gz', prefix)], shell=True, **quiet)
            subprocess.call(["cat %s | bcftools consensus %s > %s1_Final.fasta" % (
                baits_path, prefix + '1.vcf.gz', prefix)], shell=True, **quiet)
            remove_if_exists(prefix + 'mapreads.sam')
            remove_if_exists(prefix + '.0.bam')
            remove_if_exists(prefix + '0srt.bam')
            remove_if_exists(prefix + '0.vcf.gz')
            remove_if_exists(prefix + '0.vcf.gz.csi')
            remove_if_exists(prefix + '.1.bam')
            remove_if_exists(prefix + '1srt.bam')
            remove_if_exists(prefix + '1.vcf.gz')
            remove_if_exists(prefix + '1.vcf.gz.csi')
            remove_if_exists(prefix + '.chimera.bam')
            remove_if_exists(prefix + 'chimerasrt.bam')
            bam = prefix + 'mapreads.bam'
            stat_path = wf_folder + folder[:-8] + 'readstats.txt'
            with open(stat_path, 'a+') as statfile:
                statfile.write(folder[:-8] + " Read Statistics\n")
                subprocess.call(["samtools flagstat %s >> %s" % (bam, stat_path)], shell=True, **quiet)
                subprocess.call([
                    "samtools depth -a %s | awk '{c++;s+=$3}END{print s/c}' >> %s" % (
                    bam, stat_path)], shell=True, **quiet)
                subprocess.call([
                    "samtools depth -a %s | awk '{c++; if($3>0) total+=1}END{print (total/c)*100}' >> %s" % (
                    bam, stat_path)], shell=True, **quiet)
                statfile.close()


#Map reads to consensus allele references, phase bi-allelic haplotypes.
#Reads (FASTQ, bwa index) come from assemblydir (Stage 1A, read-only).
#All outputs (BAM, VCF, Final.fasta, readstats) go to workfilesdir.
sample_folders = [f for f in direc if 'assembly' in f]
with ProcessPoolExecutor(max_workers=args.threads) as pool:
    list(pool.map(_map_sample, sample_folders))


os.chdir(assemblydir)

# annotate phased fasta files with alternative phase
for folder in direc:
	if 'assembly' in folder:
		wf_folder = workfilesdir + folder + '/'
		for filename in os.listdir(wf_folder):
			if filename.endswith("0_Final.fasta"):
				filepath = wf_folder + filename
				with open(filepath, 'r') as infile:
					for line in infile:
						if '>' in line:
							if args.verbose:
								print(line)
							name = line.rstrip('\n') + '_ph0' + '\n'
							if args.verbose:
								print(name)
							replaceAll(filepath, line, name)

for folder in direc:
	if 'assembly' in folder:
		wf_folder = workfilesdir + folder + '/'
		for filename in os.listdir(wf_folder):
			if filename.endswith("1_Final.fasta"):
				filepath = wf_folder + filename
				with open(filepath, 'r') as infile:
					for line in infile:
						if '>' in line:
							if args.verbose:
								print(line)
							name = line.rstrip('\n') + '_ph1' + '\n'
							if args.verbose:
								print(name)
							replaceAll(filepath, line, name)

##Concatenate Phased seqs files
for folder in direc:
	if 'assembly' in folder:
		wf_folder = workfilesdir + folder + '/'
		finals = glob.glob(wf_folder + '*_Final.fasta')
		subprocess.call(
			"cat %s > %s" % (
				' '.join(finals),
				wf_folder + folder[:-8] + 'allcontigs_allclusterbaits_contigs_phased.fasta'),
			shell=True, **quiet)

#Move Phased allcontigs_allbaits files to phasedir
for folder in direc:
	if 'assembly' in folder:
		wf_folder = workfilesdir + folder + '/'
		for filename in os.listdir(wf_folder):
			if filename.endswith('_allcontigs_allclusterbaits_contigs_phased.fasta'):
				src = wf_folder + filename
				dst = phasedir + filename
				os.rename(src, dst)

os.chdir(diploidclusters)

DICT2 = {}

for baitcluster in os.listdir(diploidclusters):
	if baitcluster.endswith('_'):
		if baitcluster not in DICT2:
			DICT2[baitcluster] = {}

for folder in map_contigs_to_baits_dir:
		if folder.endswith('_'):
			for baitcluster in DICT2.keys():
				DICT2[baitcluster][folder + 'ph0'] = []

for folder in map_contigs_to_baits_dir:
		if folder.endswith('_'):
			for baitcluster in DICT2.keys():
				DICT2[baitcluster][folder + 'ph1'] = []

os.chdir(phasedir)

#concatenated phased sequences for all samples
subprocess.call(["cat *_allcontigs_allclusterbaits_contigs_phased.fasta  > ALLsamples_allcontigs_allbaitclusters_contigs_phased.fasta"], shell=True, **quiet)

#filling in the dictionary
input_fasta = SeqIO.parse("ALLsamples_allcontigs_allbaitclusters_contigs_phased.fasta", "fasta")
for folder in map_contigs_to_baits_dir:
	if folder.endswith('_'):
		for baitcluster in DICT2.keys():
			for record in input_fasta:
				bait = record.id.split('_', 1)[0]
				baitcluster = 'L' + bait.split('L', 1)[1] + '_' + record.id.split('_', 3)[1] + '_'
				if args.verbose:
					print(baitcluster)
				phase = (record.id.split('_', 4)[4]).rstrip('\n')
				if args.verbose:
					print(phase)
				folder = record.id.split('_', 4)[2] + '_' + record.id.split('_', 4)[3] + '_' + phase
				if args.verbose:
					print(record.id.split('_', 4)[3])
					if args.verbose:
						print(folder)
				seq = record.seq
				DICT2[baitcluster][folder].append(seq)


#write fasta output summary files by bait
for baitcluster in DICT2.keys():
	if len(DICT2[baitcluster]) > 0:
		outfile = open(baitcluster + "_allsamples_allcontigs.fasta", 'w+')
		for folder in DICT2[baitcluster].keys():
			if len(DICT2[baitcluster][folder]) > 0:
				seq_list = DICT2[baitcluster][folder]
				sorted_seq_list = sorted(seq_list, key=lambda id: int(len(seq)), reverse=True)
				for seq in sorted_seq_list:
					index = str(sorted_seq_list.index(seq))
					outfile.write(str('>' + baitcluster + folder + '_' + index + '\n' + seq + '\n'))

#output the nested dictory to a csv file
columns = [x for x in map_contigs_to_baits_dir if x.endswith('_')]

ph0 = ["{}{}".format(i, 'ph0') for i in columns]
ph1 = ["{}{}".format(i, 'ph1') for i in columns]

header = ['baitcluster'] + ph0 + ph1

with open('ALLsamples_allcontigs_allbaits_SUMMARY_TABLE_phased.csv', 'w') as outfile:
	writer = csv.writer(outfile)
	writer.writerow(header)
	first_value = list(DICT2.values())[0]
	samples = sorted(first_value.keys())
	for baitcluster in DICT2.keys():
		writer.writerow([baitcluster] + [len(DICT2[baitcluster][sample]) for sample in samples])

os.chdir(phasedir)

with open('ALLsamples_allcontigs_allbaits_SUMMARY_TABLE_phased.csv', 'w') as outfile:
	columns = [x for x in map_contigs_to_baits_dir if x.endswith('_')]
	ph0 = ["{}{}".format(i, 'ph0') for i in columns]
	ph1 = ["{}{}".format(i, 'ph1') for i in columns]
	header = ['baitcluster'] + ph0 + ph1
	writer = csv.writer(outfile)
	writer.writerow(header)

#Align and trim clusterbaits.
print('Only keeping the longest sequence for sequences with duplicate name IDs prior to alignment/trimming')
for file in os.listdir(phasedir):
	if file.endswith('allsamples_allcontigs.fasta'):
		with open(file, 'r') as infile:
					for line in infile:
						if '>' in line:
							linspl = line.split('_')
							name = linspl[0] + '_' + linspl[1] + '_' + linspl[2] + '_' + linspl[3] + '_' + linspl[4] + '\n'
							replaceAll(file, line, name)

for file in os.listdir(phasedir):
	if file.endswith('allsamples_allcontigs.fasta'):
		remove_dup(file, file[:-6] + '_duprem.fasta')

def _align_locus(file_path):
    """Per-locus alignment and trimming worker."""
    base = file_path[:-6]  # strip .fasta
    subprocess.call(
        ["mafft --globalpair --maxiterate %s %s > %s_al.fasta" % (
            args.aliter, file_path, base)],
        shell=True, **quiet)
    subprocess.call(
        ["trimal -in %s -out %s -gt %s" % (
            base + "_al.fasta", base + '_trimmed.fasta', args.indelrep)],
        shell=True, **quiet)

duprem_files = [
    os.path.join(phasedir, f)
    for f in os.listdir(phasedir)
    if f.endswith('_duprem.fasta')
]
with ProcessPoolExecutor(max_workers=args.threads) as pool:
    list(pool.map(_align_locus, duprem_files))

#deinterleave
for file in os.listdir(phasedir):
	if 'trimmed' in file:
		if args.verbose:
			print("deinterleaving " + file)
		deinterleave_fasta(file, file[:-6] + '_final.fasta')
		os.remove(file)

#compile stats; read readstats.txt from workfilesdir, write output to outdir
HETDICT = {}

for folder in os.listdir(workfilesdir):
	if folder.endswith("assembly"):
		ind = folder[:-9]
		if args.verbose:
			print(ind)
		if ind not in HETDICT:
			HETDICT[ind] = {}

for samp in os.listdir(workfilesdir):
	if samp.endswith('assembly'):
		samp_dir = workfilesdir + samp + '/'
		for readstat in os.listdir(samp_dir):
			if readstat.endswith('readstats.txt'):
				for ind in HETDICT:
					if ind in readstat:
						if args.verbose:
							print(readstat)
						with open(samp_dir + readstat, "r") as statfile:
							lines_to_read = [16, 17]
							for position, line in enumerate(statfile):
								if '+' in line:
									statlabela = line.split(" ")[3]
									statlabel = statlabela.strip('\n')
									statinta = line.split(" ")[0]
									statint = statinta.strip('\n')
									if args.verbose:
										print(statlabel)
										if args.verbose:
											print(statint)
									HETDICT[ind][statlabel] = []
									HETDICT[ind][statlabel].append(int(statint))
								else:
									if position in lines_to_read:
										if position == 16:
											statlabel = 'readdepth'
											statint = line.strip('\n')
											if args.verbose:
												print(statlabel + ' = ' + statint)
											HETDICT[ind][statlabel] = []
											HETDICT[ind][statlabel].append(int(float(statint)))
										elif position == 17:
											statlabel = 'coverage'
											statint = line.strip('\n')
											if args.verbose:
												print(statlabel + ' = ' + statint)
											HETDICT[ind][statlabel] = []
											HETDICT[ind][statlabel].append(int(float(statint)))

#get readstat values
het_values_list = list(HETDICT.values())
columns = [x for x in het_values_list[0]]
header = ['Individual'] + columns

with open(outdir + 'readstats.csv', 'w') as outfile:
	writer = csv.writer(outfile)
	writer.writerow(header)
	statlist = list(HETDICT.values())
	stats = statlist[0].keys()
	for ind in HETDICT.keys():
		writer.writerow([ind] + [HETDICT[ind].get(stat, '') for stat in stats])

with open(outdir + 'readstats.csv', 'r') as infile:
	with open(outdir + 'readstats_fin.csv', 'w') as outfile:
		data = infile.read()
		data = data.replace("]", "")
		data = data.replace("[", "")
		outfile.write(data)

os.remove(outdir + 'readstats.csv')

os.chdir(phasedir)

if 'full' in args.idformat:
	#Reformat as >L100_cl0_@@##_sampleid_phase_index ; Keeps full annotation.
	for file in os.listdir(phasedir):
		if file.endswith('final.fasta'):
			with open(file, 'r') as infile:
				for line in infile:
					if '>' in line:
						if args.verbose:
							print(line)
						linspl = line.split(' ')[0]
						linspl2 = linspl.split('_')
						if args.verbose:
							print(linspl2)
						name = linspl2[0] + '_' + linspl2[1] + '_' + linspl2[2] + '_' + linspl2[3] + '_' + linspl2[4]
						if args.verbose:
							print(name)
						replaceAll(file, line, name)
	print('Kept full sequence ID annotations; e.g. >L100_cl0_WA10_sampleid_0')
else:
	if 'phase' in args.idformat:
		for file in os.listdir(phasedir):
			if file.endswith('final.fasta'):
				with open(file, 'r') as infile:
					for line in infile:
						if '>' in line:
							if args.verbose:
								print(line)
							linspl = line.split(' ')[0]
							linspl2 = linspl.split('_')
							if args.verbose:
								print(linspl2)
							name = '>' + linspl2[2] + '_' + linspl2[3] + '_' + linspl2[4]
							if args.verbose:
								print(name)
							replaceAll(file, line, name)
		print('Annotated alignments as: >@@##_sampleid_0/1 (annotated with phase)')
	else:
		if 'onlysample' in args.idformat:
			#Reformat as >@@##_sampleid ; simplest format for concatenation.
			for file in os.listdir(phasedir):
				if file.endswith('final.fasta'):
					with open(file, 'r') as infile:
						for line in infile:
							if '>' in line:
								if args.verbose:
									print(line)
								linspl = line.split(' ')[0]
								linspl2 = linspl.split('_')
								if args.verbose:
									print(linspl2)
								name = '>' + linspl2[2] + '_' + linspl2[3] + '\n'
								if args.verbose:
									print(name)
								replaceAll(file, line, name)
			print('Annotated alignments as: >@@##_sampleid (no phase annotations)')
		else:
			sys.exit("-idformat flag not set or did not correspond to 'full', 'copies', or 'onlysample' keeping default trimal headers; e.g. >L100_cl0_WA10_sampleid_0 1230 bp ")
