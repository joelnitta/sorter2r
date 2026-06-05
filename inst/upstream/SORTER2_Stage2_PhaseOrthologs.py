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


def run_usearch(cmd):
    cwd = os.getcwd()
    if shutil.which('usearch'):
        subprocess.call('usearch ' + cmd, shell=True)
    elif shutil.which('docker'):
        docker_cmd = (
            'docker run --rm -v %s:%s -w %s '
            'joelnitta/usearch:latest %s' % (cwd, cwd, cwd, cmd)
        )
        subprocess.call(docker_cmd, shell=True)
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
args = parser.parse_args()

assemblydir = args.assemblydir if args.assemblydir.endswith('/') else args.assemblydir + '/'
clusterdir = args.clusterdir if args.clusterdir.endswith('/') else args.clusterdir + '/'
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
		print("Created folder: " + folder_path)
	else:
		print("Folder already exists: " + folder_path)

folders_to_check = ["diploids_phased"]

for folder_name in folders_to_check:
	check_folder(outdir, folder_name)

#Set Directory Variables
phasedir = outdir + 'diploids_phased/'
direc = os.listdir(assemblydir)
map_contigs_to_baits_dir = sorted(os.listdir(contigdir))

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
subprocess.call("cat *degap.fasta > ALLsamples_allcontigs_allbaits_clusterannotated.fasta", shell=True)
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
									print(line)
									break
								except StopIteration as e:
									print(e)
									break

os.chdir(assemblydir)


#Map reads to consensus allele references, phase bi-allelic haplotypes.
#Reads (FASTQ, bwa index) come from assemblydir (Stage 1A, read-only).
#All outputs (BAM, VCF, Final.fasta, readstats) go to workfilesdir.
for folder in direc:
	if 'assembly' in folder:
		wf_folder = workfilesdir + folder + '/'
		asm_folder = assemblydir + folder + '/'
		for baits in os.listdir(wf_folder):
			if baits.endswith('allcontigs_allclusterbaits_annotated.fasta'):
				baits_path = wf_folder + baits
				read = folder[:-8] + 'R1_val_1.fq'
				R2 = folder[:-8] + 'R2_val_2.fq'
				read_path = asm_folder + read
				R2_path = asm_folder + R2
				prefix = wf_folder + folder[:-8]
				print(read)
				print(R2)
				subprocess.call(["bwa index %s" % (baits_path)], shell=True)
				subprocess.call(["bwa mem -V %s %s %s > %smapreads.sam" % (
					baits_path, read_path, R2_path, prefix)], shell=True)
				subprocess.call(["samtools sort %smapreads.sam -o %smapreads.bam" % (
					prefix, prefix)], shell=True)
				subprocess.call(["samtools index %s" % (prefix + 'mapreads.bam')], shell=True)
				subprocess.call(["samtools phase -A -Q %s -b %s %s" % (
					args.phasequal, prefix, prefix + 'mapreads.bam')], shell=True)
				subprocess.call(["samtools sort %s -o %s" % (
					prefix + '.0.bam', prefix + '0srt.bam')], shell=True)
				subprocess.call(["samtools sort %s -o %s" % (
					prefix + '.1.bam', prefix + '1srt.bam')], shell=True)
				subprocess.call(["samtools sort %s -o %s" % (
					prefix + '.chimera.bam', prefix + 'chimerasrt.bam')], shell=True)
				subprocess.call([
					"bcftools mpileup -B --min-BQ 20 -Ou -d 500 -f %s %s "
					"| bcftools call -mv --ploidy 1 -Oz -p .0001 -o %s" % (
					baits_path, prefix + '0srt.bam', prefix + '0.vcf.gz')], shell=True)
				subprocess.call([
					"bcftools mpileup -B --min-BQ 20 -Ou -d 500 -f %s %s "
					"| bcftools call -mv --ploidy 1 -Oz -p .0001 -o %s" % (
					baits_path, prefix + '1srt.bam', prefix + '1.vcf.gz')], shell=True)
				subprocess.call([
					"bcftools mpileup -B --min-BQ 20 -Ou -d 500 -f %s %s "
					"| bcftools call -mv --ploidy 1 -Oz -p .0001 -o %s" % (
					baits_path, prefix + 'chimerasrt.bam', prefix + 'chimera.vcf.gz')], shell=True)
				subprocess.call(["bcftools index %s" % (prefix + '0.vcf.gz')], shell=True)
				subprocess.call(["bcftools index %s" % (prefix + '1.vcf.gz')], shell=True)
				subprocess.call(["bcftools index %s" % (prefix + 'chimera.vcf.gz')], shell=True)
				subprocess.call(["cat %s | bcftools consensus %s > %s0_Final.fasta" % (
					baits_path, prefix + '0.vcf.gz', prefix)], shell=True)
				subprocess.call(["cat %s | bcftools consensus %s > %s1_Final.fasta" % (
					baits_path, prefix + '1.vcf.gz', prefix)], shell=True)
				os.remove(prefix + 'mapreads.sam')
				os.remove(prefix + '.0.bam')
				os.remove(prefix + '0srt.bam')
				os.remove(prefix + '0.vcf.gz')
				os.remove(prefix + '0.vcf.gz.csi')
				os.remove(prefix + '.1.bam')
				os.remove(prefix + '1srt.bam')
				os.remove(prefix + '1.vcf.gz')
				os.remove(prefix + '1.vcf.gz.csi')
				os.remove(prefix + '.chimera.bam')
				os.remove(prefix + 'chimerasrt.bam')
				#generate read statistics text file
				bam = prefix + 'mapreads.bam'
				stat_path = wf_folder + folder[:-8] + 'readstats.txt'
				with open(stat_path, 'a+') as statfile:
					statfile.write(folder[:-8] + " Read Statistics\n")
					subprocess.call(["samtools flagstat %s >> %s" % (bam, stat_path)], shell=True)
					subprocess.call([
						"samtools depth -a %s | awk '{c++;s+=$3}END{print s/c}' >> %s" % (
						bam, stat_path)], shell=True)
					subprocess.call([
						"samtools depth -a %s | awk '{c++; if($3>0) total+=1}END{print (total/c)*100}' >> %s" % (
						bam, stat_path)], shell=True)
					statfile.close()


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
							print(line)
							name = line.rstrip('\n') + '_ph0' + '\n'
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
							print(line)
							name = line.rstrip('\n') + '_ph1' + '\n'
							print(name)
							replaceAll(filepath, line, name)

##Concatenate Phased seqs files
for folder in direc:
	if 'assembly' in folder:
		wf_folder = workfilesdir + folder + '/'
		os.chdir(wf_folder)
		subprocess.call(["cat *_Final.fasta > %sallcontigs_allclusterbaits_contigs_phased.fasta" % (
			folder[:-8])], shell=True)

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
subprocess.call(["cat *_allcontigs_allclusterbaits_contigs_phased.fasta  > ALLsamples_allcontigs_allbaitclusters_contigs_phased.fasta"], shell=True)

#filling in the dictionary
input_fasta = SeqIO.parse("ALLsamples_allcontigs_allbaitclusters_contigs_phased.fasta", "fasta")
for folder in map_contigs_to_baits_dir:
	if folder.endswith('_'):
		for baitcluster in DICT2.keys():
			for record in input_fasta:
				bait = record.id.split('_', 1)[0]
				baitcluster = 'L' + bait.split('L', 1)[1] + '_' + record.id.split('_', 3)[1] + '_'
				print(baitcluster)
				phase = (record.id.split('_', 4)[4]).rstrip('\n')
				print(phase)
				folder = record.id.split('_', 4)[2] + '_' + record.id.split('_', 4)[3] + '_' + phase
				print(record.id.split('_', 4)[3])
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

for file in os.listdir(phasedir):
	if file.endswith('_duprem.fasta'):
		subprocess.call(["mafft --globalpair --maxiterate %s %s > %s_al.fasta" % (args.aliter, file, file[:-6])], shell=True)
		subprocess.call(["trimal -in %s -out %s -gt %s" % (file[:-6] + "_al.fasta", file[:-6] + '_trimmed.fasta', args.indelrep)], shell=True)

#deinterleave
for file in os.listdir(phasedir):
	if 'trimmed' in file:
		print("deinterleaving " + file)
		deinterleave_fasta(file, file[:-6] + '_final.fasta')
		os.remove(file)

#compile stats; read readstats.txt from workfilesdir, write output to outdir
HETDICT = {}

for folder in os.listdir(workfilesdir):
	if folder.endswith("assembly"):
		ind = folder[:-9]
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
						print(readstat)
						with open(samp_dir + readstat, "r") as statfile:
							lines_to_read = [16, 17]
							for position, line in enumerate(statfile):
								if '+' in line:
									statlabela = line.split(" ")[3]
									statlabel = statlabela.strip('\n')
									statinta = line.split(" ")[0]
									statint = statinta.strip('\n')
									print(statlabel)
									print(statint)
									HETDICT[ind][statlabel] = []
									HETDICT[ind][statlabel].append(int(statint))
								else:
									if position in lines_to_read:
										if position == 16:
											statlabel = 'readdepth'
											statint = line.strip('\n')
											print(statlabel + ' = ' + statint)
											HETDICT[ind][statlabel] = []
											HETDICT[ind][statlabel].append(int(float(statint)))
										elif position == 17:
											statlabel = 'coverage'
											statint = line.strip('\n')
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
		writer.writerow([ind] + [HETDICT[ind][stat] for stat in stats])

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
						print(line)
						linspl = line.split(' ')[0]
						linspl2 = linspl.split('_')
						print(linspl2)
						name = linspl2[0] + '_' + linspl2[1] + '_' + linspl2[2] + '_' + linspl2[3] + '_' + linspl2[4]
						print(name)
						replaceAll(file, line, name)
	sys.exit('Kept full sequence ID annotations; e.g. >L100_cl0_WA10_sampleid_0')
else:
	if 'phase' in args.idformat:
		for file in os.listdir(phasedir):
			if file.endswith('final.fasta'):
				with open(file, 'r') as infile:
					for line in infile:
						if '>' in line:
							print(line)
							linspl = line.split(' ')[0]
							linspl2 = linspl.split('_')
							print(linspl2)
							name = '>' + linspl2[2] + '_' + linspl2[3] + '_' + linspl2[4]
							print(name)
							replaceAll(file, line, name)
		sys.exit('Annotated alignments as: >@@##_sampleid_0/1 (annotated with phase)')
	else:
		if 'onlysample' in args.idformat:
			#Reformat as >@@##_sampleid ; simplest format for concatenation.
			for file in os.listdir(phasedir):
				if file.endswith('final.fasta'):
					with open(file, 'r') as infile:
						for line in infile:
							if '>' in line:
								print(line)
								linspl = line.split(' ')[0]
								linspl2 = linspl.split('_')
								print(linspl2)
								name = '>' + linspl2[2] + '_' + linspl2[3] + '\n'
								print(name)
								replaceAll(file, line, name)
			sys.exit('Annotated alignments as: >@@##_sampleid (no phase annotations)')
		else:
			sys.exit("-idformat flag not set or did not correspond to 'full', 'copies', or 'onlysample' keeping default trimal headers; e.g. >L100_cl0_WA10_sampleid_0 1230 bp ")
