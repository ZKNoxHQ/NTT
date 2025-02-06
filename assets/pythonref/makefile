PY = python3
VENV = myenv
PIP = $(VENV)/bin/pip
PYTHON = $(VENV)/bin/python
AUX = *.pyc *.cprof */*.pyc

install:
	$(PY) -m venv $(VENV)
	$(PIP) install eth_abi pandas pycryptodome scipy matplotlib

test:
	$(PYTHON) polyntt/generate_test_vectors.py
	$(PYTHON) -m unittest  $(if $(TEST),test_$(TEST).py,discover) -v -s polyntt
	
profile:
	rm -f $(AUX)
	rm -rf __pycache__
	touch profile_action.cprof
	$(PY) -m cProfile -o profile_action.cprof profile_action.py
	pyprof2calltree -k -i profile_action.cprof &

clean:
	rm -f $(AUX)
	rm -rf __pycache__ */__pycache__
	rm -rf scripts/*.sage.py
	@echo "Clean done"
