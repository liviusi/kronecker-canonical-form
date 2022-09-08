SAGE = $(shell which sage)

.PHONY: tests

default: tests

tests:
	$(SAGE) ./kcf/kcf.sage
	@mv ./kcf/kcf.sage.py ./kcf/kcf_sage.py
	pytest ./tests/