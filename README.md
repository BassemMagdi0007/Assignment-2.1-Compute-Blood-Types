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
This project implements a Bayesian Network-based system to infer blood types and calculate probabilities in a family tree using genetic inheritance rules. The provided code is designed to solve problems based on Bayesian networks by determining the probability distribution of blood types for individuals in a family tree. The task involves modeling the inheritance of blood types through alleles and using test results to refine these probabilities.

The main objectives are:

 - **Modeling:** Represent family relationships and blood type inheritance as a Bayesian network.
 - **Inference:** Use Bayesian reasoning to calculate the likelihood of each blood type for a queried individual.
 - **Automation:** Dynamically build models from problem-specific JSON files and solve queries programmatically.

This implementation uses pgmpy, a Python library for probabilistic graphical models, to handle Bayesian networks and perform inference. The repository contains a modular implementation that reads family data, builds a probabilistic model, and queries it to answer specific questions about blood type distributions. It is designed to solve multiple problems stored in JSON format and outputs solutions in a structured format.

### Key Features 
 - **Dynamic Bayesian Network Construction:** <br>
Automatically constructs the network based on family trees, allele distributions, and test results provided in the input JSON files.

 - **Handling Allele Distributions by Country:** <br>
Supports allele distributions specific to the fictional countries North Wumponia and South Wumponia, which influence the likelihood of each blood type.

 - **Inference Automation:** <br>
Uses test results and queries to perform variable elimination, calculating blood type distributions for queried individuals.

 - **JSON-based Problem Solving:** <br>
Parses JSON problem files to extract family relationships, test results, and queries, making it adaptable to a variety of input configurations.

 - **Result Formatting:** <br>
Outputs solutions in the prescribed JSON format, ensuring compatibility with the assignment requirements.

## Setup
### Repository Content:
1) **`main.py`:** Main python script to execute the Bayesian Network creation and inference.
2) **problems/:** Problems Directory contains the JSON problem files.
3) **p-solutions/:** Solutions directory Stores the output JSON files with results in the following format:
```python
[
    {
        "type": "bloodtype",
        "person": "Rory",
        "distribution": {
            "O": 0.25,
            "A": 0.5,
            "B": 0.25,
            "AB": 0.0
        }
    }
]

```

### To set up the environment and run the code, follow these steps: 
1) **Install Python Dependencies:** Ensure you have Python installed (recommended version: Python 3.8 or later) and install the required libraries using:
    ```bash
    pip install pgmpy
    ```
2) **`main.py`** and **problems** directory must be on the same folder
3) Run the code from the used editor or from the cmd
    ```python
    python main.py
    ```

### Used libraries:
**_pgmpy_**: Python library for probabilistic graphical models that provides tools for creating, manipulating, and performing inference on Bayesian and Markov networks.
  - **BayesianNetwork:** A directed acyclic graph (DAG) that represents probabilistic dependencies among variables using nodes and edges, allowing you to model complex systems.
  - **TabularCPD:** A representation of conditional probability distributions in tabular form, which defines the probabilities of a variable given its parent variables.
  - **VariableElimination:** An inference algorithm used to compute marginal or conditional probabilities efficiently by systematically eliminating variables from the graph.

**_logging_**: The logging module is used to suppress warnings from pgmpy, which might clutter the console output.

**_json_**: The json module is utilized for reading and writing JSON files. It parses the input files containing family trees, test results, and queries (json.load) and saves the computed results in the specified output format (json.dump). This module ensures seamless handling of the JSON-based problem files.

**_os_**: The os module facilitates file and directory management. It is used to dynamically generate paths (os.path.join), create directories if they don’t exist (os.makedirs)

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

- GENOTYPE_CPD: <br>
    -  Maps combinations of alleles (A, B, O) to genotypes (AA, AO, BB, etc.).
    -  Rows represent genotypes; columns correspond to all possible allele combinations.
- OFFSPRING_CPD: <br>
    -  Represents the probabilistic inheritance from parents Genotypes to offspring alleels.
    -  Rows represent the ABO gene (A,B,O); columns correspond to genotypes.
- SUM_6_4:  <br>
    -  Maps genotypes to possible blood types (A, B, AB, O).
    -  Includes the relationship between genotype frequency and blood type probability, considering co-dominance (e.g., A results from AA and AO).
 - Countries CPDs: <br>
     -  The allele distributions for two fictional regions ("North Wumponia" and "South Wumponia") are also defined as constants.
        ```python
        cpd_north_wumponia = [[0.5], [0.25], [0.25]]
        cpd_south_wumponia = [[0.15], [0.55], [0.30]]
        ```
        
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
- Extracts relevant sections (e.g., family-tree, test-results, queries, and country) from the parsed JSON file.
- Returns a structured dictionary for easier access.

These functions ensure that input data is processed systematically and robustly, even in the presence of errors.

4) **Processing Problems**
```python
def process_problem(problem_type, problem_number):
      # ...
```
The process_problem function is the core of the script. It handles loading problem-specific data, constructing a Bayesian network based on genetic inheritance, performing inference, and saving the results in the required format.

**The function handles each problems in the following way:**
1) **Load and Parse Data:** <br>
    The function constructs the file path based on the problem type and number (e.g., problem-A-01.json). It loads the data using load_json() and extracts relevant details with extract_data().
    ```python
    filename = f'example-problems/problem-{problem_type}-{problem_number:02d}.json'
    data = load_json(filename)
    ```
    If the file is missing or invalid, the function logs an error and skips further processing.
  Then we extract relevant information from the parsed JSON data using
    ```python
        def extract_data(data):
        return {
                  #...
               }
     ```

2) **Determine CPD Based on Country:** <br>
    Based on the "country" field, it selects the appropriate allele distribution CPD (cpd_north_wumponia or cpd_south_wumponia).

3) **Build Family Structure:**
    1) Initializing Family Members: Each individual in the family_tree is dynamically represented as a dictionary with three attributes:
        - **role** (father, mother, parent, offspring).
        - **bloodtype** (e.g., A, B, AB, or O if known, otherwise None).
        - **offspring** (list of children).
        ```python
        family_members[subject] = {"role": None, "bloodtype": None, "offspring": []}
        ```
        Leading to eventually constructing the dictionary of dictionaries "family_member" in such manner: 
        ```python
        family_members = {
                "Kim": 
                {
                    "role": "parent",
                    "bloodtype": None,
                    "offspring": ["Ahmed"]
                },
                "Ahmed": 
                {
                    "role": "parent",
                    "bloodtype": "B",
                    "offspring": ["Calvin", "Linda"]
                }
            }
        ```
    2) Updating Relationships: For each relationship (father-of, mother-of, parent-of):
        -  The role of the subject (e.g., father, mother, parent) is set.
        -  The object (child) is added to the subject's offspring list.
        -  The relationships are also maintained in a relations dictionary.
        (e.g. For a family tree like:
        ```python
        [
            {"relation": "mother-of", "subject": "Kim", "object": "Ahmed"},
            {"relation": "father-of", "subject": "Ahmed", "object": "Calvin"},
            {"relation": "father-of", "subject": "Ahmed", "object": "Linda"},
            {"relation": "mother-of", "subject": "Lindsay", "object": "Dana"}
        ]
        ```
        The resulting relations dictionary will look like this:
        ```python
        relations = {
            "Kim": ["Ahmed"],
            "Ahmed": ["Calvin", "Linda"],
            "Lindsay": ["Dana"]
        }
        ```

   

























  
