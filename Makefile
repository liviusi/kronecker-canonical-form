SAGE = $(shell which sage)

FLAGS = -preparse

TARGETS = tests build

.DEFAULT_GOAL := build

build:
	$(SAGE) $(FLAGS) ./kcf/kcf.sage
	@mv ./kcf/kcf.sage.py ./kcf/kcf_sage.py

tests: build
	pytest ./tests/