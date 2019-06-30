#!python
#Author: Shenqi

import argparse
import os.path
#parse the arguments
parser=argparse.ArgumentParser(description='generate shell scripts for each sample to detect SNP')
parser.add_argument('-R',help='Reference genome sequence',required=True)
parser.add_argument('-t',help='Input file type',required=True,choices=['fq','bam','sam'])
parser.add_argument('-i',help='Input file directory',default='.')
parser.add_argument('-o',help='Output file directory',default='.')
parser.add_argument('-f',help='Used when "-t=fq", paired end rads is seperated by ":" and samples by ",". all *.fq files in input directory are read as single end',default=None)
parser.add_argument('-b',help='used when "-t=bam/sam", seperated by ",". all *.bam or *.sam files in input directory are read',default=None)
parser.add_argument('-n',help='Used when -f/-b/-s defined, sample names in output, the same order as -f/-b/-s. prefix of the input files is used',default=None)
parser.add_argument('-k',help='know vcf file',default=None)
parser.add_argument('-p',help='Used when calling non-diploid organisms',default=2)
parser.add_argument('-gff',help='Used when SNP annotation is performed, only a GTF/GFF format file with CDS annotation is supported',default=None)
parser.add_argument('-c',help='Used when "-t=fq", thread number for mapping',default=10)
parser.add_argument('-Scrip','--scriptsdir',help='scripts dir')
parser.add_argument('-Tools','--toolsdir',help='tools dir')
parser.add_argument('-bulksize',help='size of bulk')
parser.add_argument('-popstruc',help='F2/RIL')
#parser.add_argument('-gff')

argv=vars(parser.parse_args())
ref=argv['R'].strip()
type=argv['t'].strip()
indir=argv['i'].strip()
outdir=argv['o'].strip()
fq=argv['f']
bam=argv['b'].strip()
name=argv['n'].strip()
known=argv['k']
ploidy=argv['p']
gff=argv['gff'].strip()
core=argv['c']
bulksize=int(argv['bulksize'].strip())
popstruc=argv['popstruc'].strip()
if argv['scriptsdir'] == None:
    raise Exception('You should provide the Script dir !')
else:
    scripts_dir = argv['scriptsdir']

if argv['toolsdir'] == None:
    raise Exception('You should provide the Tools dir !')
else:
    tools_dir =argv['toolsdir']

#gff=argv['gff'].strip()

#Configure
#gatk='/PUBLIC/source/RNA/RefRNA/TransRef/GATKsnp/bin/GenomeAnalysisTK.jar'

gatk = scripts_dir+'/Script/SNP_TR/GATKsnp/GenomeAnalysisTK.jar'
bwa = tools_dir+'/bwa/bwa-0.7.6a/bwa'
samtools = tools_dir+'/samtools/samtools-0.1.18/samtools'
picardtools = tools_dir+'/picard/picard-tools-1.96'
java = tools_dir+'/jre/jre1.7.0_25/bin/java'
scriptsdir = scripts_dir+'/Script/SNP_TR/GATKsnp/bin'

dirfa=os.path.dirname(ref)
fa=os.path.basename(ref)
if fa[-3:] == '.fa':
    prefa=fa[:-3]
if fa[-6:] == '.fasta':
    prefa=fa[:-6]
if 'indir' not in vars():
    indir=os.getcwd()
if 'outdir' not in vars():
    outdir=os.getcwd()
if 'ploidy' not in vars():
    ploidy=2
if 'core' not in vars():
    core=10
if not os.path.exists(outdir):
    os.system('mkdir %s' %(outdir))

# Reading input files and sample names
fa_list=[]
name_list=[]
if type=='fq':
    if 'fq' in vars():
        fq_list=fq.split(',')
        if 'name' in vars():
            name_list=name.split(',')
        else:
            for eachfq in fq_list:
                name=eachfq.split(':')[0][:-3]
                name_list.append(name)
    else:
        for eachfile in os.listdir(indir):
            if eachfile[-3:]=='.fq':
                fq_list.append(eachfile)
        for eachfq in fa_list:
            name=eachfq[:-3]
            name_list.append(name)

bam_list=[]
name_list=[]
if type=='bam' or type=='sam':
    if 'bam' in vars():
        bam_list=bam.split(',')
        if 'name' in vars():
            name_list=name.split(',')
        else:
            for eachbam in bam_list:
                name=eachbam[:-4]
                name_list.append(name)
    else:
        for eachfile in os.listdir(indir):
            if eachfile[-4:]=='.bam' or eachfile[-4:]=='sam':
                bam_list.append(eachfile)
            for eachbam in bam_list:
                name_list.append(eachbam[:-4])
# generate a shell script to create directories and prepare reference index files
os.chdir(outdir)
os.chdir('..')
workflow1=open('workflow1.sh','w')
workflow1.write('mkdir %s/ref %s/tmp\n' %(outdir,outdir))
workflow1.write('ln -s %s/%s %s/ref/%s\n' %(dirfa,fa,outdir,fa))
workflow1.write('%s faidx %s/ref/%s\n' %(samtools,outdir,fa))
workflow1.write('%s/CreateSequenceDictionary.jar R=%s/ref/%s o=%s/ref/%s.dict TMP_DIR=%s/tmp\n' %(picardtools,outdir,fa,outdir,prefa,outdir))
workflow1.close()

if type=='fq' or type=='sam':
    suffix='sam'
if type=='bam':
    suffix='bam'

# generate a shell script for each sample to precess data
os.system('mkdir %s/dedup' %(outdir))
for file in name_list:
    eachfile=open('workflow2_'+file+'.sh','w')
    #eachfile.write('mkdir %s/dedup\n' %(outdir))
    # add RG and sort
    eachfile.write('%s -jar %s/AddOrReplaceReadGroups.jar I=%s/%s.%s O=%s/dedup/%s.sorted.bam SO=coordinate ID=%s PL=illumina LB=LB PU=PU SM=%s TMP_DIR=%s/tmp VALIDATION_STRINGENCY=SILENT\n' %(java,picardtools,indir,file,suffix,outdir,file,file,file,outdir))
    # remove secondary alignment reads in some bam files (e.g. Tophats)
    eachfile.write('%s view -bF 0x100 %s/dedup/%s.sorted.bam > %s/dedup/%s.rmdup.bam\n' %(samtools,outdir,file,outdir,file))
    # duplicate marking
    eachfile.write('%s -jar %s/MarkDuplicates.jar INPUT=%s/dedup/%s.rmdup.bam OUTPUT=%s/dedup/%s.dedup.bam METRICS_FILE=%s/dedup/%s.metricsfile MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=250 ASSUME_SORTED=true TMP_DIR=%s/tmp VALIDATION_STRINGENCY=SILENT REMOVE_DUPLICATES=True\n' %(java,picardtools,outdir,file,outdir,file,outdir,file,outdir))
    # reorder
    eachfile.write ('%s -jar %s/ReorderSam.jar I=%s/dedup/%s.dedup.bam O=%s/dedup/%s.reorder.bam R=%s/ref/%s TMP_DIR=%s/tmp VALIDATION_STRINGENCY=SILENT\n' %(java,picardtools,outdir,file,outdir,file,outdir,fa,outdir))
    # index reorder.bam
    eachfile.write('%s index %s/dedup/%s.reorder.bam\n' %(samtools,outdir,file))
    eachfile.write('rm %s/dedup/%s.sorted.bam %s/dedup/%s.rmdup.bam %s/dedup/%s.dedup.bam\n' %(outdir,file,outdir,file,outdir,file))
    eachfile.close()

# generate a shell script to train SNP detection and variant analysis
workflow3=open('workflow3.sh','w')
if 'know' not in vars():
    workflow3.write('echo ==========Training SNP Detection============\ndate\n')
    workflow3.write('mkdir %s/SNPQ30\n' %(outdir))
    # variant calling
    workflow3.write('%s -jar %s -T UnifiedGenotyper --filter_reads_with_N_cigar -R %s/ref/%s -o %s/SNPQ30/rawVariantsQ30.vcf -stand_call_conf 30 -stand_emit_conf 10 -glm BOTH -ploidy %s' %(java,gatk,outdir,fa,outdir,ploidy))
    for eachfile in name_list:
        workflow3.write(' -I %s/dedup/%s.reorder.bam' %(outdir,eachfile))
    workflow3.write('\n')
    # extract SNPs from the call set
    workflow3.write('%s -jar %s -T SelectVariants -R %s/ref/%s -V %s/SNPQ30/rawVariantsQ30.vcf -o %s/SNPQ30/rawSNPQ30.vcf -selectType SNP\n' %(java,gatk,outdir,fa,outdir,outdir))
    # SNP filtration
    workflow3.write('%s -jar %s -T VariantFiltration -R %s/ref/%s -V %s/SNPQ30/rawSNPQ30.vcf -o %s/SNPQ30/filtSNPQ30.vcf --clusterWindowSize 10 --filterExpression "MQ0 >= 4 && ((MQ0/(1.0*DP))>0.1)" --filterName "HARD_TO_VALIDATE" --filterExpression "QUAL < 30.0 || QD < 5.0 || FS > 60.0" --filterName "GATKStandard"\n' %(java,gatk,outdir,fa,outdir,outdir))
    # extract InDels from the call set
    workflow3.write('%s -jar %s -T SelectVariants -R %s/ref/%s -V %s/SNPQ30/rawVariantsQ30.vcf -o %s/SNPQ30/rawInDelQ30.vcf -selectType INDEL\n' %(java,gatk,outdir,fa,outdir,outdir))
    workflow3.write('%s -jar %s -T VariantFiltration -R %s/ref/%s -V %s/SNPQ30/rawInDelQ30.vcf -o %s/SNPQ30/filtInDelQ30.vcf --filterExpression "QD < 2.0 || FS > 200.0" --filterName "GATKStandard"\n' %(java,gatk,outdir,fa,outdir,outdir))
    # define know vcf file
    known='%s/SNPQ30/filtSNPQ30.vcf' %(outdir)
    # Variant Analysis
    workflow3.write('echo ===============Variant Analysis=================\n')
    workflow3.write('mkdir %s/ResultsQ30\n' %(outdir))
    # extract pass SNPs
    workflow3.write('grep "PASS\\|^#" %s/SNPQ30/filtSNPQ30.vcf > %s/ResultsQ30/SNPpassQ30.vcf\n' %(outdir,outdir))
    workflow3.write('grep "PASS\\|^#" %s/SNPQ30/filtInDelQ30.vcf > %s/ResultsQ30/InDelpassQ30.vcf\n' %(outdir,outdir))
    # vcf format conversion
    workflow3.write('perl %s/formatVCF.pl %s/ResultsQ30/SNPpassQ30.vcf > %s/ResultsQ30/SNPpassQ30.tmp\n' %(scriptsdir,outdir,outdir))
    workflow3.write('perl %s/formatVCF.pl %s/ResultsQ30/InDelpassQ30.vcf > %s/ResultsQ30/InDels.xls\n' %(scriptsdir,outdir,outdir))
    # additional distance filter
    workflow3.write('perl %s/distanceFilt.pl -i %s/ResultsQ30/SNPpassQ30.tmp -o %s/ResultsQ30/SNPs.xls\n' %(scriptsdir,outdir,outdir))
    # variant statistics and plot
    # locate SNPs to gene
    workflow3.write('python %s/snp2gene.py --gtf %s --snp %s/ResultsQ30/SNPs.xls --out %s/ResultsQ30/SNPs2gene.res --filter n\n' %(scriptsdir,gff,outdir,outdir))
    workflow3.write('mv %s/ResultsQ30/SNPs2gene.res %s/ResultsQ30/SNPs.xls\n' %(outdir,outdir))
    workflow3.write('python %s/snp2gene.py --gtf %s --snp %s/ResultsQ30/InDels.xls --out %s/ResultsQ30/InDels2gene.res --filter n\n' %(scriptsdir,gff,outdir,outdir))
    workflow3.write('mv %s/ResultsQ30/InDels2gene.res %s/ResultsQ30/InDels.xls\n' %(outdir,outdir))
    workflow3.write('perl %s/variantPlot_v1.1.pl -S %s/ResultsQ30/SNPs.xls -I %s/ResultsQ30/InDels.xls -o %s/ResultsQ30\n' %(scriptsdir,outdir,outdir,outdir))
workflow3.close()

# generate a shell script to convert VCF file to table format
os.chdir(outdir)
os.chdir('..')
workflow4 = open('workflow4.sh','w')
workflow4.write('%s -Xmx8g -jar %s -T VariantsToTable -R %s -V %s/ResultsQ30/SNPpassQ30.vcf -F CHROM -F POS -F REF -F ALT -GF AD -GF DP -GF GQ -GF PL -o %s/SNP.table\n' %(java,gatk,ref,outdir,outdir))

# generate Chromosome number flle
chromlist = open('chromlist.txt','w')
with open(ref, 'r') as f:
    for line in f:
        if '>' in line:
            chromlist.write(line[1:].strip()+'\n')
chromlist.close()


# generate a shell script to run QTLseqr

os.system('cp /TJPROJ2/RESEQ/Project_Crop/Crop6/shenqi/BSRScripts/QTLseqr.R %s/../QTLseqr.R' %(outdir))

highbulk =  str(bam.split(',')[0])
lowbulk = str(bam.split(',')[1])

workflow5 = open('workflow5.sh','w')
workflow5.write('export PATH=/PUBLIC/software/RESEQ/software/anaconda3/setup/bin/:$PATH\n')
workflow5.write('unset PYTHONPATH\n')
workflow5.write('export PATH=/PUBLIC/software/RESEQ/software/anaconda3/setup/envs/perl/bin/:$PATH\n')
workflow5.write('source activate R350\n')
workflow5.write('''Rscript %s/../QTLseqr.R --data  %s/SNP.table --highbulk %s --lowbulk %s --q_value 0.01 --windowsize 1e6 --bulksize  %d --popstruc  %s\n''' %(outdir,outdir,highbulk,lowbulk,bulksize,popstruc))

# generate a python script to qsub the shell scripts generated above

code='''#!python

import misopy.cluster_utils as cluster_utils
import re
import os

node = 'crop1.q,crop2.q'
def generate_qsub(sh,vf,p):
    cluster_script='qsub -terse -V -cwd -l vf=%dG,p=%d -q %s %s' %(int(vf),int(p),node,sh)
    return cluster_script

# qsub the workflow1.sh
run_workflow1=generate_qsub("workflow1.sh",1,1)
job_id=cluster_utils.launch_job(run_workflow1,cmd_name='qsub')
cluster_utils.wait_on_jobs([job_id], cluster_cmd='qsub')

print "1done"

# qsub the workflow_sample.sh of each sample parallel
samples=[]
job_ids=[]
for file in os.listdir("."):
    if re.search(r'^workflow2_.*\.sh',file):
        samples.append(re.findall(r'^workflow2_.*\.sh',file)[0])
#print samples
for sample in samples:
    run_workflow2=generate_qsub(sample,10,6)
    job_ids.append(cluster_utils.launch_job(run_workflow2,cmd_name='qsub'))
cluster_utils.wait_on_jobs(job_ids, cluster_cmd='qsub')

print "2done"
# qsub the workflow3.sh
run_workflow3=generate_qsub("workflow3.sh",10,6)
job_id3=cluster_utils.launch_job(run_workflow3,cmd_name='qsub')
cluster_utils.wait_on_jobs([job_id3], cluster_cmd='qsub')

# qsub the workflow4.sh
run_workflow4=generate_qsub("workflow4.sh",1,1)
job_id4=cluster_utils.launch_job(run_workflow4,cmd_name='qsub')
cluster_utils.wait_on_jobs([job_id4], cluster_cmd='qsub')

# qsub the workflow5.sh
os.system('qsub -V -cwd -l vf=5G,p=2 -q crop1.q,crop2.q workflow5.sh')

'''
open('qsub_workflow.py','w').write(code)
qsub = open('qsub.sh','w')
qsub.write('unset PYTHONPATH\n')
qsub.write('nohup /PUBLIC/software/RESEQ/software/anaconda2/setup/bin/python qsub_workflow.py &\n')
