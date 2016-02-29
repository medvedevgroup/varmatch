all: vm vt
.PHONY: all vm vt clean

vm:
	$(MAKE) -C vm all
	chmod +x varmatch
	
vt:
	$(MAKE) -C vt 

clean:
	$(MAKE) -C vm clean
	$(MAKE) -C vt clean

