.PHONY: all
all: data.h

node_modules:
	mkdir -p node_modules

node_modules/webppl-json: node_modules
	cd node_modules; git clone https://github.com/stuhlmueller/webppl-json

data.json: data.wppl node_modules/webppl-json
	webppl data.wppl --require webppl-json

data.h: data.json
	./data.py >$@

.PHONY: clean
clean:
	rm -f data.h data.json
