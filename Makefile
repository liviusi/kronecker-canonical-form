SAGE = $(shell which sage)

TARGETS = tests build

.DEFAULT_GOAL := build

build:
	$(SAGE) ./kcf/kcf.sage
	@mv ./kcf/kcf.sage.py ./kcf/kcf_sage.py

tests: build
	pytest ./tests/