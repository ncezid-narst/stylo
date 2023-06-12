#!/usr/bin/env make

# Manual makefile
# Major targets are: install, test, clean

.DELETE_ON_ERROR:

.default: install

install: t/.data-uncompress.done
	@echo "Please edit t/environment.bash before running make test"

# Make all large files depend on git-lfs
t/data/.git-lfs.done:
	git lfs pull origin
	touch $@

# Uncompress large file(s)
t/.data-uncompress.done: t/data/.git-lfs.done t/data/SalmonellaLitchfield.FAST5.tar.xz
	tar --directory t/data -Jxvf t/data/SalmonellaLitchfield.FAST5.tar.xz
	zcat t/data/polished.fasta.gz > t/data/polished.fasta
	zcat t/data/unpolished.fasta.gz > t/data/unpolished.fasta
	makeblastdb -in t/data/polished.fasta -dbtype nucl
	touch $@

# Wrapper for cleaning all the things
clean: clean-large-files clean-tests

# Clean only the large files
clean-large-files:
	rm -vf t/data/SalmonellaLitchfield.FAST5.tar.xz
	rm -rvf t/data/SalmonellaLitchfield.FAST5
	rm -vf t/data/*.gz
	rm -vf t/data/*.fasta
	rm -vf t/.data-uncompress.done
	rm -vf t/data/.git-lfs.done

# Clean test results
clean-tests:
	rm -f t/*.done
	rm -rvf t/vanilla.project

# For any test, run it with bats and then touch the file
t/%.done:
	exe=t/$$(basename $@ .done) && \
		echo bats $$exe &&\
		bats $$exe
	touch $@

test: t/np_basecall-demux_guppy.sh.done t/np_prepSample_readLengths.sh.done t/np_assemble_wtdbg2.sh.done t/np_polish_medaka.sh.done

# The subset of tests to be done on travis, limited by prerequisites
test-travis: t/np_basecall-demux_guppy.sh.done t/np_prepSample_readLengths.sh.done t/np_assemble_wtdbg2.sh.done t/np_polish_medaka.sh.done

