# GitHub Actions Workflows

[![Compilation with Optimizations](https://github.com/michel-leonard/C-Quadratic-Sieve/actions/workflows/compilation-with-optimizations.yaml/badge.svg)](https://github.com/michel-leonard/C-Quadratic-Sieve/actions/workflows/compilation-with-optimizations.yaml)
[![Cross-Platform Testing](https://github.com/michel-leonard/C-Quadratic-Sieve/actions/workflows/cross-platform.yaml/badge.svg)](https://github.com/michel-leonard/C-Quadratic-Sieve/actions/workflows/cross-platform.yaml)
[![JSON and CSV Output Generation](https://github.com/michel-leonard/C-Quadratic-Sieve/actions/workflows/json-csv-outputs.yaml/badge.svg)](https://github.com/michel-leonard/C-Quadratic-Sieve/actions/workflows/json-csv-outputs.yaml)
[![Large Numbers Factorization](https://github.com/michel-leonard/C-Quadratic-Sieve/actions/workflows/large-numbers.yml/badge.svg)](https://github.com/michel-leonard/C-Quadratic-Sieve/actions/workflows/large-numbers.yml)
[![Medium-Sized Numbers Factorization](https://github.com/michel-leonard/C-Quadratic-Sieve/actions/workflows/medium-sized-numbers.yml/badge.svg)](https://github.com/michel-leonard/C-Quadratic-Sieve/actions/workflows/medium-sized-numbers.yml)
[![Memory Safety](https://github.com/michel-leonard/C-Quadratic-Sieve/actions/workflows/memory-safety.yaml/badge.svg)](https://github.com/michel-leonard/C-Quadratic-Sieve/actions/workflows/memory-safety.yaml)
[![Optimal Multiplier Selection](https://github.com/michel-leonard/C-Quadratic-Sieve/actions/workflows/optimal-multiplier.yaml/badge.svg)](https://github.com/michel-leonard/C-Quadratic-Sieve/actions/workflows/optimal-multiplier.yaml)
[![Small Numbers Factorization](https://github.com/michel-leonard/C-Quadratic-Sieve/actions/workflows/small-numbers.yml/badge.svg)](https://github.com/michel-leonard/C-Quadratic-Sieve/actions/workflows/small-numbers.yml)

## Introduction
The project is progressing well, and several GitHub Actions workflows have been implemented to automate essential tasks. These workflows ensure code quality, project portability, and result reliability.

Each workflow can be triggered manually via `workflow_dispatch`, but scheduled runs (`schedule`) are also defined. These scheduled executions help ensure, with a slight delay, that everything continues to function correctly over time.

## List of Workflows

### 1. Compilation with Optimizations
Compiles the program with different optimization levels (`-O0`, `-O1`, `-O2`, `-O3`, `-Ofast`) and compares execution times to assess their impact.

### 2. Cross-Platform Testing
Verifies the project's compatibility by compiling and testing on **Ubuntu, macOS, and Windows**. Integrity checks are performed by comparing SHA256 hashes of results across systems.

### 3. JSON and CSV Output Generation
Tests the generation of results in JSON and CSV formats while validating their accuracy using a Python script.

### 4. Large Numbers Factorization
Tests the factorization of large numbers while monitoring performance and ensuring validated results.

### 5. Medium-Sized Numbers Factorization
Similar to the previous workflow but applied to a smaller range of numbers.

### 6. Memory Safety
Uses **Valgrind** to detect potential memory leaks during program execution.

### 7. Optimal Multiplier Selection
Evaluates the impact of the `qs-multiplier` parameter by comparing performance across different number sizes.

### 8. Small Numbers Factorization
Runs the program on smaller numbers to validate its correct operation in this context.

## Using the Workflows
- To trigger a workflow manually, go to the **Actions** tab in the GitHub repository and select **Run workflow**.
- Execution logs are accessible directly from GitHub Actions.
- In case of failure, check the error messages and adjust the code or workflow parameters as needed.

These workflows provide continuous monitoring of the project's proper functioning and performance. They also help quickly identify potential regressions.
