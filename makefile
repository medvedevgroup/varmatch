all: vm
.PHONY: all vm clean

vm:
	$(MAKE) -C src all
	chmod +x varmatch
	chmod +x purify
	chmod +x filter

clean:
	$(MAKE) -C src clean

