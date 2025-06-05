# Compute-Blood-Types

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

**The main objectives are:**

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
**The script defines CPDs and constants that represent genetic inheritance patterns and blood type probabilities:**

- **GENOTYPE_CPD:** <br>
    -  Maps combinations of alleles (A, B, O) to genotypes (AA, AO, BB, etc.).
    -  Rows represent genotypes; columns correspond to all possible allele combinations.
- **OFFSPRING_CPD:** <br>
    -  Represents the probabilistic inheritance from parents Genotypes to offspring alleels.
    -  Rows represent the ABO gene (A,B,O); columns correspond to genotypes.
- **SUM_6_4:**  <br>
    -  Maps genotypes to possible blood types (A, B, AB, O).
    -  Includes the relationship between genotype frequency and blood type probability, considering co-dominance (e.g., A results from AA and AO).
 - **Countries CPDs:** <br>
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
**Extracts specific information from the JSON object:**
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
    The function constructs the file path based on the problem type and number (e.g., `problem-a-01.json`). It loads the data using `load_json()` and extracts relevant details with `extract_data()`.
    ```python
    filename = f'example-problems/problem-{problem_type}-{problem_number:02d}.json'
    data = load_json(filename)
    ```
    ```python
    extracted_data = extract_data(data)
    family_tree = extracted_data["family_tree"]
    test_results = extracted_data["test_results"]
    queries = extracted_data["queries"]
    country = extracted_data["country"]
    ```
    If the file is missing or invalid, the function logs an error and skips further processing.


2) **Determine CPD Based on Country:** <br>
    Based on the "country" field, it selects the appropriate allele distribution CPD (cpd_north_wumponia or cpd_south_wumponia).

3) **Build Family Structure:**
    1) **Initializing Family Members:** <br>
    Each individual in the `family_tree` is dynamically represented as a dictionary with three attributes:
        - **role** (father, mother, parent, offspring).
        - **bloodtype** (e.g., A, B, AB, or O if known, otherwise None).
        - **offspring** (list of children).
        ```python
        family_members[subject] = {"role": None, "bloodtype": None, "offspring": []}
        ```
        Leading to eventually constructing the dictionary of dictionaries `family_member` in such manner:
       
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
    3) **Updating Relationships:** <br>
    **For each relationship (father-of, mother-of, parent-of):**
        -  The role of the subject (e.g., father, mother, parent) is set.
        -  The object (child) is added to the subject's offspring list.
        -  The relationships are also maintained in a `relations` dictionary.
        (e.g. For a family tree like:
        ```python
        [
            {"relation": "mother-of", "subject": "Kim", "object": "Ahmed"},
            {"relation": "father-of", "subject": "Ahmed", "object": "Calvin"},
            {"relation": "father-of", "subject": "Ahmed", "object": "Linda"},
            {"relation": "mother-of", "subject": "Lindsay", "object": "Dana"}
        ]
        ```
        The resulting `relations` dictionary will look like this:
        ```python
        relations = {
            "Kim": ["Ahmed"],
            "Ahmed": ["Calvin", "Linda"],
            "Lindsay": ["Dana"]
        }
        ```
        This structure is later used to connect nodes in the Bayesian network and model inheritance patterns.
       
   5) **Assigning Blood Types:** <br>
   **Known blood types from `test_results` are linked to the corresponding individuals in the `family_members` dictionary:**
        ```python
        for result in test_results:
            person = result.get("person")
            family_members[person]["bloodtype"] = result.get("result")
        ```
    
       The preceding processes (i, ii, iii) will transform the data into a more structured and visually coherent format, which will enhance both data processing and debugging efficiency.
   
      **For instance, consider the example from `problem-a-00.json`:** <br>
      
         ```python
              "family-tree": [
          {
              "relation": "father-of",
              "subject": "Ayansh",
              "object": "Dana"
          },
          {
              "relation": "mother-of",
              "subject": "Lyn",
              "object": "Dana"
          }
         "test-results": [
          {
              "type": "bloodtype-test",
              "person": "Lyn",
              "result": "A"
          }
         ```
      **The debugging outputs would be as follows:** <br>
      
        ```python
        FATHER: Ayansh ( )
        OFFSPRING: Dana ( )
        MOTHER: Lyn (A)
        ```
  
        ```python
        MEMBER:  Ayansh
        INFO:  {'role': 'father', 'bloodtype': None, 'offspring': ['Dana']}
        MEMBER:  Lyn
        INFO:  {'role': 'mother', 'bloodtype': 'A', 'offspring': ['Dana']}
        ```
  
        ```python
        member:  Dana
        MOTHER:  ['Lyn']
        FATHER:  ['Ayansh']
        PARENT:  []
        ```  
  4) **Bayesian Network Construction:**
      - **Initializing Nodes: Each family member is represented in the Bayesian network by:**
          - Two nodes for alleles (Allele1, Allele2).
          - One node for genotype (Genotype).
        ```python
        complete_model.add_nodes_from([allele1, allele2, f"{member}_Genotype"])
        ```  
       - **Assigning CPDs:**
          - **For individuals without known parentage:**
              - Their alleles are directly assigned CPDs based on the regional distribution (country_cpd).
              ```python
              cpd_allele1 = TabularCPD(variable=allele1, variable_card=3, values=country_cpd)
              ```
          - **For offsprings:**
               - Their alleles depend on the genotype of their parents (if known) or the country_cpd otherwise.
              ```python
              if father:
                cpd_allele1 = TabularCPD(variable=allele1, variable_card=3, evidence[f{father[0]}_Genotype"],...)
              else:
                cpd_allele1 = TabularCPD(variable=allele1, variable_card=3, values=country_cpd)
              ```
       - **Defining Relationships Between Nodes:**         
          - Alleles (Allele1, Allele2) determine the genotype:
            ```python
            complete_model.add_edges_from([(allele1, f"{member}_Genotype"), (allele2, f"{member}_Genotype")]))
            ```          
          - Parent genotypes influence child alleles:
            ```python
            complete_model.add_edge(f"{person}_Genotype", f"{s}_Allele1")
            ```
            
  5) **Performing Inference:**
      - **Blood Type Queries:**
          - For queried individuals or those with known blood types, a `Bloodtype` node is added, and its CPD is set based on genotype:
          ```python
          cpd = TabularCPD(variable=bloodtype_node, variable_card=4, evidence=[f"{member}_Genotype"], ...)
          ```
      - **Evidence Integration:**
          - Known blood types from the family are used as evidence to refine the probability distributions:
          ```python
          evidence[f"{member}_Bloodtype"] = ['A', 'B', 'O', 'AB'].index(info["bloodtype"])
          ```
      - **Inference Execution:**
          - The `VariableElimination` method calculates:
               - Blood type probabilities for queries.
               - Genotypes and alleles for individuals not constrained by evidence.
          ```python
          overall_distribution = inference_complete.query(variables=[inference_variable], evidence=evidence)
          ```
          
## Self Evaluation and Design Decisions

### 1) Different approaches trials to create the Bayesian Network:

### [Network v1 Structure]: 
**The network consisted of the following nodes:**
- **Father's Allele (Allele1):** Represents one allele contributed by the father.
- **Mother's Allele (Allele2):** Represents one allele contributed by the mother.
- Offspring's Genotype: Combines Allele1 and Allele2 to determine the offspring's blood type.
As for the CPDs, calculations was done manually for a given bloodtypes; <br>

  <img width="269" alt="image" src="https://github.com/user-attachments/assets/b5d523dc-8584-4ac3-89b0-239c436ee4fa" />

**For example in `problem-a-00.json`:** <br>

  <img width="514" alt="image" src="https://github.com/user-attachments/assets/6d5f534a-ee97-471f-a249-e342d96a560b" />

**Given the following information from the input JSON file:**
- **Country:** North Wumponia
- **Mother’s Blood Type:** A

**The CPDs (Conditional Probability Distributions) were assigned as follows:**
- **Mother's Allele (Allele2):** The blood type A CPD was directly assigned to Allele2 to represent the mother's contribution.

- **Father's Allele (Allele1):** Since the father’s blood type was not provided, the Allele1 node was assigned the CPD specific to North Wumponia (cpd_north_wumponia).

#### **Limitations:**
**While this approach worked for inferring the offspring’s blood type distribution, it presented several limitations:**
- **Restricted Query Capabilities:**
   - The network was designed to support queries only for the offspring’s genotype, as this node had a variable cardinality of 4 (representing the blood types A, B, O, AB).
   - Parental nodes could not support similar queries due to the lack of structure to infer their genotypes or alleles.
- **Challenges with Parent Queries:** 
   - To query a parent's genotype, that parent node would need to be redefined as the Genotype node. This required:
      - Adjusting the edges, as the Father's Genotype and Mother's Genotype are independent.
      - Adding a new node to incorporate either the offspring’s blood type (if available) or external country-specific CPD data.
     
These changes introduced unnecessary complexity, making the network rigid and less adaptable for future problem scenarios.

### [Network v2 Structure]:
**The network consisted of the following nodes (FOR EACH MEMBER):**
- Allele1
- Allele2
- Genotype
  
The same idea of v1 was applied but in this network every person has his own network (Allele1, Allele2, Genotype) but there's no relation/connection between the parents genotypes and the offspring alleles. <br>

<img width="527" alt="image" src="https://github.com/user-attachments/assets/96fb0e49-d3b7-4bb9-9c17-894dc3030dfa" />

#### **Limitations:**
- **Dependency on Offspring Blood Type:**
    - A parent is dependent to the offspring bloodtype, so to establish this connection the offspring’s blood type CPD was assigned in a hardcoded way to the parent’s allele node. This approach lacked generality and introduced rigid dependencies that could not easily adapt to more complex scenarios.
- **Scaling Issues with Complex Networks:**
    - While the network produced correct results for all `problem-a-*` cases, its performance degraded as the problem complexity increased. The hardcoded dependencies and limited flexibility of the network design led to inaccuracies and incorrect inferences in more advanced scenarios.. 

### [Network Final Structure]:
In this Bayesian network structure, it consists of three types of nodes: 
- Alleles
- Genotypes
- Bloodtypes.
  
<img width="533" alt="image" src="https://github.com/user-attachments/assets/713373d0-742b-47bd-b39e-81f5a1fcb70b" />

Each individual is represented by two allele nodes `Member_Allele1` and `Member_Allele2` that combine to form their genotype node `Member_Genotype` via directed edges. The genotype node determines the blood type node `Member_Bloodtype` using predefined conditional probability distributions. Parent-child relationships are represented by directed edges from a parent's genotype node to the child's allele nodes `Parent_Genotype` → `Member_Allele1`. This structure reflects the inheritance of genetic traits, with the network explicitly encoding how alleles from parents combine and influence offspring's genotypes and blood types.

### Adressing cheap-bloodtype-test problems:
The cheap-blood-test problems introduce a 20% probability of incorrect test results, which complicated the inference process. <br>
Two approaches were attempted: <br>
- **Approach 1: Applying the 20% Error to Test Results:**
  ```python
  if result.get("type") == "cheap-bloodtype-test":
      if random.random() < 0.2:  # 20% error
          other_person = random.choice([p for p in family_members.keys() if p != person])
          family_members[person]["bloodtype"] = family_members[other_person]["bloodtype"]
      else:
          family_members[person]["bloodtype"] = result.get("result")
  ```
    This approach did not yield the correct results. The error propagation seemed to disrupt the network’s inference, as introducing random incorrect blood types caused inconsistencies in the CPD values.

- **Approach 2: Applying the 20% Error to Final Results:** <br>
After running the Bayesian network normally and obtaining the results, a 20% error was applied to the computed distributions for queried individuals. A 20% chance was introduced to randomly adjust these probabilities. <br>
This approach was also incorrect. Modifying the results post-inference failed to reflect the dependencies and relationships inherent in the Bayesian network.

**Conclusion for Cheap-Blood-Test Problems:** <br>
Despite multiple attempts, I was unable to implement an accurate solution that accounted for the 20% error in test results while preserving the integrity of the Bayesian network.

### Adressing no country specidfied problems:
For problems without a specified country, the CPD for alleles is undefined, and each country (North Wumponia and South Wumponia) must be considered equally likely.  <br>
Several methods were explored: <br>

- **Approach 1: Randomly Choosing a Country:** <br>
For each problem, a single country (North or South Wumponia) was selected randomly with a 50% probability, and the Bayesian network was solved using that country’s CPD.
  ```python
  cpd_unspecified_country = random.choice([cpd_north_wumponia, cpd_south_wumponia])
  ```
  This method yielded inconsistent results, as the random selection led to significant variations and incorrect answers.

- **Approach 2: Averaging Results from Both Countries:** <br>
The network was run twice for each problem—once with North Wumponia’s CPD and once with South Wumponia’s CPD. The results were averaged to account for the equal likelihood of each country.
  ```python
  north_results = run_network(cpd_north_wumponia)
  south_results = run_network(cpd_south_wumponia)
  cpd_unspecified_country = (north_results + south_results) / 2
  ```
  The results were still incorrect. Averaging results seemed to oversimplify the relationships between alleles and blood types.

- **Approach 3: Combining CPDs:**
The allele distributions from both countries were combined and averaged to create a new CPD for unspecified countries:
  ```python
  cpd_unspecified_country = [[(0.5 * north + 0.5 * south) for north, south in zip(north_row, south_row)] 
                             for north_row, south_row in zip(cpd_north_wumponia, cpd_south_wumponia)]
  ```
  This approach produced results close to the correct answers to some problems like `problem-e-00.json` but was not entirely accurate:
  ```python
    Current Answer: {'O': 0.275, 'A': 0.325, 'B': 0.4, 'AB': 0.0})
    Correct Answer: {'O': 0.2795, 'A': 0.2934, 'B': 0.4270, 'AB': 0.0},
  ```
  Although the combined CPD approach came close to the correct results, it did not fully capture the nuances of the problem. As a result, the exact solution for this problem type remains unresolved.

## Output Format
The computed probabilities are formatted into a dictionary for each query and stored in a list: 
```python
result = {
    "type": "bloodtype",
    "person": person,
    "distribution": { "O": ..., "A": ..., "B": ..., "AB": ... }
}
```
The results are saved to a JSON file in the p-solutions/ directory:
```python
with open(output_filename, 'w') as outfile:
    json.dump(results, outfile, indent=4)
```
The final format of the output in the files: 
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

  
