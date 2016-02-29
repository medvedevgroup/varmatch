all: vcf vt
.PHONY: all vcf vt clean

vcf:
	$(MAKE) -C src all
	-chmod +x vcfcompare
	
vt:
	$(MAKE) -C vt all

clean:
	$(MAKE) -C src clean
	$(MAKE) -C vt clean

