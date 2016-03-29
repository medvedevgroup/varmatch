all: vm vt
.PHONY: all vm vt clean

vm:
	$(MAKE) -C src all
	chmod +x varmatch
	chmod +x purify
	chmod +x filter
	
vt:
	$(MAKE) -C vt 

clean:
	$(MAKE) -C src clean
	$(MAKE) -C vt clean

