.PHONY: install lint format format-check type-check test test-cov verify fix build clean

PYTHON ?= 3.12

install:
	uv sync --group dev
	uv run --python $(PYTHON) python -c "import dftlearn; print(dftlearn.__version__)"

lint:
	uv run ruff check src

format:
	uv run ruff format src

format-check:
	uv run ruff format --check src

type-check:
	uvx ty check src

test:
	uv run --python $(PYTHON) pytest

test-cov:
	uv run --python $(PYTHON) pytest --cov=dftlearn --cov-report=xml --cov-report=term-missing

verify: lint format-check test

fix:
	uv run ruff check --fix src
	uv run ruff format src

build:
	uv build

clean:
	rm -rf dist build .pytest_cache .ruff_cache .coverage coverage.xml
	find . -type d -name __pycache__ -prune -exec rm -rf {} +
