# Latex Makefile using latexmk
# Modified by Dogukan Cagatay <dcagatay@gmail.com>
# Modified by Philipp Jund
# Originally from : http://tex.stackexchange.com/a/40759


PROJNAME=report_main
OUT_DIR=out
.PHONY: $(PROJNAME).pdf all clean

all: dir $(PROJNAME).pdf

dir:
	mkdir -p $(OUT_DIR)
	mkdir -p $(OUT_DIR)/chapters

$(PROJNAME).pdf: $(PROJNAME).tex
	latexmk -outdir=$(OUT_DIR) -pdf -use-make -file-line-error $<
	cp $(OUT_DIR)/$(PROJNAME).pdf $(PROJNAME).pdf

cleanall:
	rm -rf $(OUT_DIR)/*

clean:
	latexmk -outdir=$(OUT_DIR)/ -c

# Build report whenever a source file is saved.
#
# The timeout is necessary in case your source contains invalid latex syntax.
# Just fix the error, wait until the timeout terminates the stuck build process
# and save again. For large projects/slow hardware you might need to increase
# the timeout of 10s.
watch:
	evince $(PROJNAME).pdf &
	while true; do \
		timeout 10s make; \
		EXIT_CODE=$$?; \
		if [ $$EXIT_CODE -ne 0 ]; then \
			make clean; \
		fi; \
		inotifywait -qre close_write --exclude "out/.*" .; \
	done
