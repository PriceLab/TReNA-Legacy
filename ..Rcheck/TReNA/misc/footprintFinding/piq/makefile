# environment variables for pshannon
# PIQDIR=~/github/piq-single
# WORKDIR=~/s/work/priceLab/cory/piq
# BAMFILE=~/s/work/priceLab/cory/piq/ENCSR000DBY.chr20.bam
# OUTPUTDIR=~/s/work/priceLab/cory/piq/piq.out.hg38/
# TMPDIR=~/s/work/priceLab/cory/piq/tmp

# environment variables for cory on aws
PIQDIR=/home/piq-single
WORKDIR=/scratch/piq
BAMFILE=/scratch/piq/ENCSR000DBY.chr20.bam
OUTPUTDIR=/scratch/piq/results/
TMPDIR=/scratch/piq/tmp

default:
	echo "testInstall targets: pwmHits, bamConvert, call"
        
testInstall:
	(cd $(PIQDIR); R -f common.r)

pwmHits:
	(cd $(PIQDIR); Rscript pwmmatch.exact.r common.r pwms/jasparfix.txt 139 $(OUTPUTDIR))

bamConvert:
	(cd $(PIQDIR); Rscript bam2rdata.r common.r $(OUTPUTDIR)chr20.RData $(BAMFILE))

call:
	(cd $(PIQDIR); Rscript pertf.r common.r $(OUTPUTDIR) $(TMPDIR) $(OUTPUTDIR) $(OUTPUTDIR)chr20.RData 139)


