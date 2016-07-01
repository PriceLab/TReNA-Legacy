PIQDIR=~/github/piq-single
WORKDIR=~/s/work/priceLab/cory/piq
BAMFILE=~/s/work/priceLab/cory/piq/ENCSR000DBY.chr20.bam
OUTPUTDIR=~/s/work/priceLab/cory/piq/piq.out.hg38/
TMPDIR=~/s/work/priceLab/cory/piq/tmp

default:
	echo "testInstall targets: pwmHits, bamConvert, call"
        
testInstall:
	(cd $(PIQDIR); R -f common.r)

pwmHits:
	(cd $(PIQDIR); Rscript pwmmatch.exact.r common.r pwms/jasparfix.txt 139 $(OUTPUTDIR))

bamConvert:
	(cd $(PIQDIR); Rscript bam2rdata.r common.r $(OUTPUTDIR)d0.RData $(BAMFILE))

call:
	(cd $(PIQDIR); Rscript pertf.r common.r $(OUTPUTDIR) $(TMPDIR) $(OUTPUTDIR) $(OUTPUTDIR)chr20.RData 139)
#	(cd $(PIQDIR); Rscript pertf.r common.r $(OUTPUTDIR) $(TMPDIR) $(OUTPUTDIR) $(OUTPUTDIR)d0.RData 139)


