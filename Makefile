.PHONY: setup format lint check test clean

setup:
	python3 -m venv .venv
	.venv/bin/pip install --upgrade pip
	.venv/bin/pip install -e .
	.venv/bin/pip install black isort ruff flake8 mypy

install:
	.venv/bin/pip install -e .

format:
	.venv/bin/black src/
	.venv/bin/isort src/

lint:
	.venv/bin/ruff check src/

check: format lint

test:
	.venv/bin/asap-to-kite -ff tests/data1 -sp test1 -of tests/output -on test1 -c 2

clean:
	rm -rf build dist *.egg-info src/*.egg-info
	rm -rf .ruff_cache .venv
	rm -rf tests/output
	find . -type d -name __pycache__ -exec rm -rf {} +
	find . -type f -name "*.pyc" -delete