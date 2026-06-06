import os
import re
import sys
import shutil
import argparse
import subprocess

parser = argparse.ArgumentParser()
parser.add_argument("-o", "--outdir", required=True, help='Output directory for assembled results')
parser.add_argument("-spades", "--spades", required=True, help='Run Spades Assembly? (T/F)')
parser.add_argument("-trim", "--trim", required=True, help='Run Trim Galore? (T/F)')
parser.add_argument("-clean_tmp", "--clean_tmp", default='F',
                    help='Delete SPAdes K-mer graph dirs after each sample? (T/F)')
args = parser.parse_args()


def _clean_kmer_dirs(spades_out):
    """Remove K[0-9]+ subdirectories from a SPAdes output directory."""
    if not os.path.isdir(spades_out):
        return
    for item in os.listdir(spades_out):
        if re.match(r'^K\d+$', item):
            shutil.rmtree(os.path.join(spades_out, item), ignore_errors=True)


#working directory iterator
rootwd=os.getcwd()+'/'

#Output directory
dst = args.outdir if args.outdir.endswith('/') else args.outdir + '/'
os.makedirs(dst, exist_ok=True)
wdlist = os.listdir(rootwd)
print('Reads will be processed in:\n' + dst)

if args.trim == 'T':
#Trim reads, put paired reads in folders corresponding to each sample

	for file in wdlist:
		if 'R1' in file:
			filefolder=file.split('_')[0]+'_'+file.split('_')[1]+'_assembly'
			readst= dst + filefolder + '/'
			os.makedirs(os.path.join(readst))
			print('Moving ' + file + ' to: \n' + readst)
			R1_old_path = os.path.join(rootwd, file)
			R1_new_path = os.path.join(readst, file)
			R2_old_path = os.path.join(rootwd, file.replace('_R1.','_R2.'))
			R2_new_path = os.path.join(readst, file.replace('_R1.','_R2.'))
			os.renames(R1_old_path, R1_new_path)
			os.renames(R2_old_path, R2_new_path)
			os.chdir(readst)
			subprocess.run(["trim_galore --quality 20 --length 30 --paired --fastqc %s %s" % (file, file.replace('_R1.','_R2.'))], shell=True)
			os.remove(file)
			os.remove(file.replace('_R1.','_R2.'))
			os.chdir(rootwd)
		else:
			continue

#Assemble Contigs with Spades
if args.spades == 'T':

	if args.trim == 'T':
		for file in os.listdir(dst):
			if 'assembly' in file:
				print("Running Spades on " + file)
				os.chdir(dst + file)
				for read in os.listdir(dst + file):
					if 'R1_val_1.fq' in read:
						R2 = read[:-11] + 'R2_val_2.fq'
						subprocess.call(["spades.py --only-assembler -1 %s -2 %s -o spades_hybrid_assembly" % (read, R2)], shell=True)
						if args.clean_tmp == 'T':
							_clean_kmer_dirs('spades_hybrid_assembly')
						os.chdir(dst)
					else:
						continue
			else:
				continue
	else:
		for file in wdlist:
			if 'R1' in file and file.endswith('.fastq'):
				filefolder=file.split('_')[0]+'_'+file.split('_')[1]+'_assembly'
				readst= dst + filefolder + '/'
				os.makedirs(readst, exist_ok=True)
				print("Running Spades on " + file)
				R1 = os.path.join(rootwd, file)
				R2 = os.path.join(rootwd, file.replace('_R1.','_R2.'))
				os.chdir(readst)
				subprocess.call(["spades.py --only-assembler -1 %s -2 %s -o spades_hybrid_assembly" % (R1, R2)], shell=True)
				if args.clean_tmp == 'T':
					_clean_kmer_dirs('spades_hybrid_assembly')
				os.chdir(rootwd)
			else:
				continue


print("Trimgalore and SPADES processing has finished, exiting script")
sys.exit(0)
