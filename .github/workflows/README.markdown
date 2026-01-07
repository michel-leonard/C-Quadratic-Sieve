# GitHub Actions Workflows

[![Cross-Platform Testing](https://github.com/michel-leonard/C-Quadratic-Sieve/actions/workflows/cross-platform.yml/badge.svg)](https://github.com/michel-leonard/C-Quadratic-Sieve/actions/workflows/cross-platform.yml)
[![Compilation with Optimizations](https://github.com/michel-leonard/C-Quadratic-Sieve/actions/workflows/compilation-with-optimizations.yml/badge.svg)](https://github.com/michel-leonard/C-Quadratic-Sieve/actions/workflows/compilation-with-optimizations.yml)
[![Verbosity Functional Independence](https://github.com/michel-leonard/C-Quadratic-Sieve/actions/workflows/verbosity-functional-independence.yml/badge.svg)](https://github.com/michel-leonard/C-Quadratic-Sieve/actions/workflows/verbosity-functional-independence.yml)
[![JSON and CSV Output Generation](https://github.com/michel-leonard/C-Quadratic-Sieve/actions/workflows/json-csv-outputs.yml/badge.svg)](https://github.com/michel-leonard/C-Quadratic-Sieve/actions/workflows/json-csv-outputs.yml)
[![Factorize very large numbers](https://github.com/michel-leonard/C-Quadratic-Sieve/actions/workflows/very--large-numbers.yml/badge.svg)](https://github.com/michel-leonard/C-Quadratic-Sieve/actions/workflows/very--large-numbers.yml)
[![Large Numbers Factorization](https://github.com/michel-leonard/C-Quadratic-Sieve/actions/workflows/large-numbers.yml/badge.svg)](https://github.com/michel-leonard/C-Quadratic-Sieve/actions/workflows/large-numbers.yml)
[![Medium-Sized Numbers Factorization](https://github.com/michel-leonard/C-Quadratic-Sieve/actions/workflows/medium-sized-numbers.yml/badge.svg)](https://github.com/michel-leonard/C-Quadratic-Sieve/actions/workflows/medium-sized-numbers.yml)
[![Factorize OEIS Sequences 1](https://github.com/michel-leonard/C-Quadratic-Sieve/actions/workflows/factorize-oeis-1.yml/badge.svg)](https://github.com/michel-leonard/C-Quadratic-Sieve/actions/workflows/factorize-oeis-1.yml)
[![Factorize OEIS Sequences 2](https://github.com/michel-leonard/C-Quadratic-Sieve/actions/workflows/factorize-oeis-2.yml/badge.svg)](https://github.com/michel-leonard/C-Quadratic-Sieve/actions/workflows/factorize-oeis-2.yml)
[![Optimal Multiplier Selection](https://github.com/michel-leonard/C-Quadratic-Sieve/actions/workflows/optimal-multiplier.yml/badge.svg)](https://github.com/michel-leonard/C-Quadratic-Sieve/actions/workflows/optimal-multiplier.yml)
[![Small Numbers Factorization](https://github.com/michel-leonard/C-Quadratic-Sieve/actions/workflows/small-numbers.yml/badge.svg)](https://github.com/michel-leonard/C-Quadratic-Sieve/actions/workflows/small-numbers.yml)


## Introduction
The project is progressing well, and several workflows have been implemented to automate essential tasks. These workflows ensure code quality, project portability, and result reliability.

Each workflow can be triggered manually via `workflow_dispatch`, but scheduled runs (`schedule`) are also defined. These scheduled executions help ensure, with a slight delay, that everything continues to function correctly over time.

On average, workflows take about **1 minute** to execute, with some heavier workflows like **Cross-Platform Testing** and **Optimal Multiplier Selection** taking **2 to 3 minutes**. Originally, this project serves as a foundation for exploring GitHub Actions (version 4).

## Workflow Categories

### **Regression Prevention**
- **Cross-Platform Testing** - Ensures compatibility across different operating systems (Ubuntu, macOS, Windows).
- **Verbosity Functional Independence** - Prevents side effects through continuous testing.
- **JSON and CSV Output Generation** - Validates correctness of output formats.

### **Deployment Assurance**
- **Compilation with Optimizations** - Verifies compilation and execution performance using both GCC and Clang.
- **Small Numbers Factorization** - Ensures basic correctness with small inputs.

### **Performance Evaluation**
- **Very Large Numbers Factorization** - Measures performance in processing 260-bit numbers.
- **Large Numbers Factorization** - Measures performance with large inputs.
- **Medium-Sized Numbers Factorization** - Tests intermediate workloads.
- **OEIS Factorization** - Tests a factorization workflow based on the Online Encyclopedia of Integer Sequences.
- **Optimal Multiplier Selection** - Evaluates the impact of the `--qs-multiplier` parameter.

## Using the Workflows
- To trigger a workflow manually, go to the **Actions** tab in the GitHub repository and select **Run workflow**.
- Execution logs are accessible directly from GitHub Actions.
- In case of failure, check the error messages and adjust the code or workflow parameters as needed.

These workflows provide continuous monitoring of the project's proper functioning and performance. They also help quickly identify potential regressions.

## Production Readiness

As of now, this software is considered **production ready** as a **single-thread stateless factorizer**.
