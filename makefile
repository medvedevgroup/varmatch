all: vcf vt
.PHONY: all vcf vt clean

vcf:
	$(MAKE) -C concurrent/concurrent all
	-chmod +x vcfcompare
	
vt:
	$(MAKE) -C vt all

clean:
	$(MAKE) -C concurrent/concurrent clean
	$(MAKE) -C vt clean

