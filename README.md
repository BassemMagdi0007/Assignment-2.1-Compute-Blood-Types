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
This project implements a Bayesian Network-based system to infer blood types and calculate probabilities in a family tree using genetic inheritance rules. The repository contains a modular implementation that reads family data, builds a probabilistic model, and queries it to answer specific questions about blood type distributions. It is designed to solve multiple problems stored in JSON format and outputs solutions in a structured format.

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
2) **`main.py`** and **problems** folder must be on the same folder
3) Run the code from the used editor or from the cmd **python main.py**.

### Used libraries:
**_pgmpy_**: Python library for probabilistic graphical models that provides tools for creating, manipulating, and performing inference on Bayesian and Markov networks.
  - **BayesianNetwork:** A directed acyclic graph (DAG) that represents probabilistic dependencies among variables using nodes and edges, allowing you to model complex systems.


  - **TabularCPD:** A representation of conditional probability distributions in tabular form, which defines the probabilities of a variable given its parent variables.
  - **VariableElimination:** An inference algorithm used to compute marginal or conditional probabilities efficiently by systematically eliminating variables from the graph.


## Code Structure
1) **Library imports**
```python
import json
import logging
import os
from pgmpy.models import BayesianNetwork
from pgmpy.factors.discrete import TabularCPD
from pgmpy.inference import VariableElimination
```

2) **Predefined Constants and CPDs**
```python
GENOTYPE_CPD = [
    #AA   AB   AO   BA   BB   BO   OA   OB   OO
    [1.0, 0.0, ...],  # AA
    ...
]
OFFSPIRING_CPD = [
    # AA   AO   BB   BO   OO   AB
    [1.0, 0.5, ...],  # A
    ...
]
```
The script defines CPDs and constants that represent genetic inheritance patterns and blood type probabilities:

- GENOTYPE_CPD: Maps combinations of alleles (A, B, O) to genotypes (AA, AO, BB, etc.).
- OFFSPRING_CPD: Models inheritance probabilities of alleles from parents.
- SUM_6_4: Maps genotypes to possible blood types (A, B, AB, O). <br>

The allele distributions for two fictional regions ("North Wumponia" and "South Wumponia") are also defined as constants.

3) **Helper Functions**
```python
def load_json(filepath):
      # ...
```
Reads a JSON file from the specified path and parses its content. If the file is missing or contains invalid JSON, an error message is printed.

```python
def extract_data(data):
      # ...
```
Extracts specific information from the JSON object:
- Family tree relationships.
- Blood type test results.
- Queries for blood type prediction.
- The country of origin.
