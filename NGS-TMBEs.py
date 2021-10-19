import os
import re
from Bio import SeqIO
from Bio.Seq import Seq

# Place fastq and reference files into a folder
# Modify these according to your situation
folder='/Users/hekezhang/Desktop/MTs-ACBE-U-G-1_FKDL210026767-1a-N710-N503'  # Input the folder path
f1='/Users/hekezhang/Desktop/MTs-ACBE-U-G-1_FKDL210026767-1a-N710-N503/MTs-ACBE-U-G-1_FKDL210026767-1a-N710-N503_1.fq.gz'  # Input the read1 path
f2='/Users/hekezhang/Desktop/MTs-ACBE-U-G-1_FKDL210026767-1a-N710-N503/MTs-ACBE-U-G-1_FKDL210026767-1a-N710-N503_2.fq.gz'  # Input the read2 path
ref="/Users/hekezhang/Desktop/MTs-ACBE-U-G-1_FKDL210026767-1a-N710-N503/EGFP1.fa"  # Input the reference file path
FwdPrimer = "gagggcgatgccacctacgg"  # Input the forward primer sequence
RevPrimer = "gtcggccatgatatagacgttg"  # Input the reverse primer sequence
sitenumbers=366  # Input the amplicon length


ali1=[]
ali2=[]
# Alignment1 for variant analysis, mutation numbers per read
def align1():
    global ali1
    os.system("/opt/anaconda3/bin/bwa mem -t 6 " + ref+ " " + folder + "'/temp.fastq' > " + folder + "'/temp.sam'")
    os.system("/opt/anaconda3/bin/samtools fixmate -O bam -@ 5 " + folder + "'/temp.sam' " + folder + "'/fixmate_temp.bam'")
    os.system("/opt/anaconda3/bin/samtools sort -O bam -@ 5 -o " + folder + "'/sorted_temp.bam' " + folder + "'/fixmate_temp.bam'")
    os.system("/opt/anaconda3/bin/samtools index " + folder + "'/sorted_temp.bam'")
    os.system("/opt/anaconda3/bin/samtools tview -d T -w 500 " + folder + "'/sorted_temp.bam' " +ref+ " > " + folder +"'/alignment.txt'")
    with open(folder + '/alignment.txt', "r") as f:
        lines = f.readlines()
    for line in lines[3:]:
        ali1.append(line.strip().replace('*',''))

# Alignment2 for variant analysis, mutation numbers per site
def align2():
    global ali2
    os.system("/opt/anaconda3/bin/bwa mem -t 6 " + ref + " " + folder + "'/Equallength_reads.fastq' > " + folder + "'/sam.sam'")
    os.system("/opt/anaconda3/bin/samtools fixmate -O bam -@ 5 " + folder + "'/sam.sam' " + folder + "'/fixmate.bam'")
    os.system("/opt/anaconda3/bin/samtools sort -O bam -@ 5 -o " + folder + "'/sorted.bam' " + folder + "'/fixmate.bam'")
    os.system("/opt/anaconda3/bin/samtools index " + folder + "'/sorted.bam'")
    os.system("/opt/anaconda3/bin/samtools mpileup -R -d 0 -f " + ref + " " + folder + "'/sorted.bam' > " + folder + "'/alignment.mpileup.txt'")
    with open(folder + '/alignment.mpileup.txt', "r") as f:
        lines = f.readlines()
    for line in lines:
        temp = line.strip('\n').split('\t')
        indeltemp = re.findall("\d+", temp[4])
        for j in indeltemp:
            temp[4] = re.sub('[+-][0-9]+[ATCGatcg]{}'.format({int(j)}), '', temp[4], count=1, flags=0)
        ali2.append(temp)

# Variant analysis , mutation numbers per site, mutation numbers per read and the number of mutation combinations
def analy():
    rmtnumber = []
    diffs = list(set(ali1))
    for read in ali1:
        mprN = read.count('A') + read.count('T') + read.count('C') + read.count('G')+read.count('a') + read.count('t') + read.count('c') + read.count('g')
        rmtnumber.append(mprN)
    f = open(folder + '/mutation_statistics.txt', 'w')
    f.write('Equallength reads number(Alignment depth): ' + '\t' + str(count3) + '\n')
    f.write('Indel reads number: ' + '\t' + str(count4) + '\n')
    f.write('Mutation combinations: ' + '\t' + str(len(diffs)) + '\n\n')
    f.write('Mutation per Read' + '\t' + 'Count ' + '\n')
    for i in range(max(rmtnumber)+1):
        f.write(str(i) + '\t' + str(rmtnumber.count(i)) + '\n')
    f.write('\nMutation numbers per site' + '\n')
    f.write('Site' + '\t' + 'Base' + '\t' + 'A' + '\t' + 'T' + '\t' + 'C' + '\t' + 'G' + '\n')
    for temp in ali2:
        mpsA = temp[4].count('A')+temp[4].count('a')
        mpsT = temp[4].count('T')+temp[4].count('t')
        mpsC = temp[4].count('C')+temp[4].count('c')
        mpsG = temp[4].count('G')+temp[4].count('g')
        f.write(temp[1] + '\t' + temp[2] + '\t' + str(mpsA) + '\t' + str(mpsT) + '\t' + str(mpsC) + '\t' + str(mpsG) + '\n')
    f.close()

# Quality filtering & merging (merge each pair of reads into a single read if they are overlapped)
os.system("/opt/anaconda3/bin/fastp -w 6 -5 -3 -M 30 -q 30 -u 30 -i "+f1+" -I "+f2+ " -h "+folder+"'/fastp.html'" " -m --merged_out "+folder+"'/merged.fastq'" )

# Make a file containing reads that can be aligned with the reference sequence and have the forward and reverse primer (to remove non-target reads)
os.system("/opt/anaconda3/bin/bwa index " + ref)
os.system("/opt/anaconda3/bin/bwa mem -t 6 " + ref + " " + folder + "'/merged.fastq' > " + folder + "'/temp.sam'")
os.system("/opt/anaconda3/bin/samtools fixmate -O bam -r -@ 5 " + folder + "'/temp.sam' " + folder + "'/fixmate_temp.bam'")
os.system("/opt/anaconda3/bin/samtools sort -n -O bam -@ 5 -o " + folder + "'/sorted_temp.bam' " + folder + "'/fixmate_temp.bam'")
os.system("/opt/anaconda3/bin/samtools fastq -@ 5 -0 " + folder + "'/aligned.fastq' " + folder + "'/sorted_temp.bam'")
primer_match = []
unprimer_match = []
for rec in SeqIO.parse(folder + '/aligned.fastq', "fastq"):
    if (FwdPrimer.upper() == rec[0:len(FwdPrimer)].seq) & (Seq(RevPrimer.upper()).reverse_complement() == rec.seq[-len(RevPrimer):]):
        primer_match.append(rec)
    else:
        unprimer_match.append(rec)
count1 = SeqIO.write(primer_match, folder + "/primer_match.fastq", "fastq")
count2 = SeqIO.write(unprimer_match, folder + "/unprimer_match.fastq", "fastq")
print("Found %i reads matched with primers, %i reads unmatched with primers" % (count1, count2))

# Detect Indel reads
Equallength_reads = []
Indel_reads = []
for rec in SeqIO.parse(folder + "/primer_match.fastq", "fastq"):
    if len(rec.seq) == sitenumbers:
        Equallength_reads.append(rec)
    else:
        Indel_reads.append(rec)
count3 = SeqIO.write(Equallength_reads, folder + "/Equallength_reads.fastq","fastq")
count4 = SeqIO.write(Indel_reads, folder + "/Indel_reads.fastq", "fastq")

# fastqc
os.system("/opt/anaconda3/bin/fastqc -t 2 -q " + folder + "'/Equallength_reads.fastq' " + folder + "'/Indel_reads.fastq'" )
print("Found %i Equal length reads, %i Indel reads" % (count3, count4))

# Split sequences, because samtools tview only supports viewing 8000 reads
tempfq=[]
n=1
for rec in SeqIO.parse(folder+"/Equallength_reads.fastq", "fastq"):
    tempfq.append(rec)
    n+=1
    if n % 8000 == 0:
        SeqIO.write(tempfq, folder + "/temp.fastq", "fastq")
        tempfq = []
        align1()
if n% 8000 != 0:
    align1()

align2()
analy()

print('finish')