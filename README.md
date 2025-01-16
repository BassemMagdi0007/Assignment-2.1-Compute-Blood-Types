# Assignment-2.1-Compute-Blood-Types

This repository contains a Python-based implementation of a Bayesian Network designed to predict blood types and solve familial relationship-based queries using probabilistic inference. The program leverages the pgmpy library for Bayesian Network construction and variable elimination.

## Table of Contents

- [Introduction](#introduction)
  - Key Features
- [Setup](#setup)
  - Repository content
  - How to run the code
  - Used libraries
- [Code Structure](#code-structure)
- [Self Evaluation and Design Decisions](#design-decision)
- [Output Format](#output-format)

## Introduction
This project focuses on creating a Bayesian Network to model genetic inheritance of blood types in families and analyze probabilistic queries about relationships. It incorporates data from JSON files to dynamically build family structures and apply predefined Conditional Probability Distributions (CPDs) for inference.

### Key Features 
- Dynamic construction of Bayesian Networks based on input family data.
- Implementation of genetic inheritance rules for blood type prediction.
- Use of Variable Elimination for probabilistic inference.
- JSON-based input and output for ease of data handling.

## Setup
### This repository contains:
1) main.py: Main script to execute the Bayesian Network creation and inference.
2) problems/: Contains problem JSON files.
3) p-solutions/: Generated solutions after executing the code

### How to run the code: 
1) Make sure that pgmpy library is intalled if not use **pip install pgmpy**
1) **`maint.py`** and **problems** folder must be on the same folder
2) Run the code from the used editor or from the cmd **python main.py**.

### Used libraries:
**_os_**
**_pgmpy_**


