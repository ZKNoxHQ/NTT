PY = python3
VENV = myenv
PIP = $(VENV)/bin/pip
PYTHON = $(VENV)/bin/python
AUX = *.pyc *.cprof */*.pyc

install:
	$(PY) -m venv $(VENV)
	$(PIP) install pycryptodome

test:
	$(PYTHON) -m polyntt.generate_test_vectors
	$(PYTHON) -m unittest  $(if $(TEST),polyntt.test_$(TEST),discover) -v

bench:
	$(PYTHON) polyntt/bench_poly_mul.py
	$(PYTHON) polyntt/bench_iterative_recursive.py
	
clean:
	rm -f $(AUX)
	rm -rf __pycache__ */__pycache__
	rm -rf scripts/*.sage.py
	@echo "Clean done"
