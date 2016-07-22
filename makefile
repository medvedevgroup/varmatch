all: vm
.PHONY: all vm clean

vm:
	$(MAKE) -C src all
	chmod +x varmatch

clean:
	$(MAKE) -C src clean

