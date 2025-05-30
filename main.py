import json
import logging
import os
import random
from pgmpy.models import DiscreteBayesianNetwork
from pgmpy.factors.discrete import TabularCPD
from pgmpy.inference import VariableElimination
import glob

'''------------------------------------------------------------------------------------------------'''
# Suppress pgmpy warnings
logging.getLogger("pgmpy").setLevel(logging.ERROR)
# Alleles, NORTH: |A: [0.75,0.0,0.25], |B: [0.0,0.6667,0.3333], |O: [0.0,0.0,1.0], |AB: [0.5,0.5,0.0]
# Alleles, SOUTH: |A: [0.6,0.0,0.4],   |B: [0.0,0.7392,0.2608], |O: [0.0,0.0,1.0], |AB: [0.5,0.5,0.0]
'''Pre-defined Conditional Probability Distributions (CPDs) for the alleles and genotypes'''
GENOTYPE_CPD = [
    #AA   AB   AO   BA   BB   BO   OA   OB   OO
    [1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],  # AA
    [0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0],  # AO
    [0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0],  # BB
    [0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 1.0, 0.0],  # BO
    [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0],  # OO
    [0.0, 1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0],  # AB
    ]

OFFSPIRING_CPD = [
    # AA   AO   BB  BO   OO   AB
    [1.0, 0.5, 0.0, 0.0, 0.0, 0.5],  # A
    [0.0, 0.0, 1.0, 0.5, 0.0, 0.5],  # B
    [0.0, 0.5, 0.0, 0.5, 1.0, 0.0],  # O
]

SUM_6_4 = [
    # AA   AO   BB  BO     OO   AB
    [1.0, 1.0, 0.0, 0.0, 0.0, 0.0],  # A
    [0.0, 0.0, 1.0, 1.0, 0.0, 0.0],  # B
    [0.0, 0.0, 0.0, 0.0, 1.0, 0.0],  # O
    [0.0, 0.0, 0.0, 0.0, 0.0, 1.0],  # AB
]

'''------------------------------------------------------------------------------------------------'''
'''Read the JSON file and return its data'''
# Load a JSON file and return its data
def load_json(filepath):
    try:
        with open(filepath, 'r') as f:
            return json.load(f)
    except FileNotFoundError:
        print(f"Error loading JSON file: [Errno 2] No such file or directory: '{filepath}'")
    except json.JSONDecodeError:
        print(f"Error decoding JSON file: '{filepath}'")
    return None

# Extract relevant information from the parsed JSON data
def extract_data(data):
    return {
        #List of dictoinaries containing the family tree
        "family_tree": data.get("family-tree", []),
        #List of dictionaries containing the test results
        "test_results": data.get("test-results", []),
        #List of dictionaries containing the queries
        "queries": data.get("queries", []),
        #String containing the country
        "country": data.get("country", None)
    }

'''------------------------------------------------------------------------------------------------'''
def process_problem(problem_type, problem_number):
    # Conditional Probability Distributions (CPDs) for the alleles and genotypes
    cpd_north_wumponia = [[0.5], [0.25], [0.25]]
    cpd_south_wumponia = [[0.15], [0.55], [0.30]]
    # TODO: WORNG
    cpd_unspecified_country = [[0.325], [0.4], [0.275]]

    '''--------------------------------------------------------------------------------------------'''
    ''''Load and extract data from JSON file, check the country and define the country CPD'''
    # Load and extract data from JSON file
    filename = f'example-problems/problem-{problem_type}-{problem_number:02d}.json'
    data = load_json(filename)
    if not data:
        print(f"Skipping problem {problem_number} due to missing data.")
        return

    extracted_data = extract_data(data)
    #List of dictoinaries containing the family tree
    family_tree = extracted_data["family_tree"]
    #List of dictionaries containing the test results
    test_results = extracted_data["test_results"]
    #List of dictionaries containing the queries
    queries = extracted_data["queries"]
    #String containing the country
    country = extracted_data["country"]

    # Define country_cpd based on the country
    use_country_node = False
    if country == "North Wumponia":
        country_cpd = cpd_north_wumponia
    elif country == "South Wumponia":
        country_cpd = cpd_south_wumponia
    elif country is None:
        use_country_node = True
        # 3x2 CPD for alleles: rows are A, B, O; columns are North, South
        founder_allele_cpd = [
            [0.5, 0.15],   # A
            [0.25, 0.55],  # B
            [0.25, 0.30],  # O
        ]
    else:
        print(f"Skipping problem {problem_number} due to invalid or missing country: {country}")
        return
    
    '''--------------------------------------------------------------------------------------------'''
    '''Define the family members  and their relations'''
    # Dynamically define the family members and their relations
    # Dictionary of dictionaries
    family_members = {}
    ''' 
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
        }'''
    # a Dictionary of lists
    relations = {}
    '''' OUTPUT
    relations = {
            "Kim": ["Ahmed"],
            "Ahmed": ["Calvin", "Linda"],
            "Lindsay": ["Dana"]
        }'''

    # Iterate through every object in the family tree and update the family members and relations
    for key in family_tree:
        subject = key["subject"]
        object_ = key["object"]
        relation_type = key["relation"]

        # 1) INITIALIZE THE FAMILY MEMBERS AND RELATIONS
        # CREATE A DICTIONARY FOR EACH SUBJECT AND OBJECT IN THE FAMILY TREE IF THEY DO NOT EXIST
        if subject not in family_members:
            # a Dictionary
            family_members[subject] = {"role": None, "bloodtype": None, "offspring": []}
        if object_ not in family_members:
            # a Dictionary
            family_members[object_] = {"role": None, "bloodtype": None, "offspring": []}

        # 2) UPDATE THE ROLES AND RELATIONS FOR THE SUBJECT AND OBJECT
        # For father set the role and add the offspring
        if relation_type == "father-of":
            family_members[subject]["role"] = "father"
            family_members[subject]["offspring"].append(object_)
            # If the subject(key) is not in the relations dictionary, add it and set the value to an empty list
            if subject not in relations:
                relations[subject] = []
            # Append the object to the list of relations for the subject as his offspring
            relations[subject].append(object_)
            # Set the role of the object to offspring in the family_members dictionary
            family_members[object_]["role"] = "offspring"

        # For mother set the role and add the offspring
        elif relation_type == "mother-of":
            family_members[subject]["role"] = "mother"
            family_members[subject]["offspring"].append(object_)
            if subject not in relations:
                relations[subject] = []
            relations[subject].append(object_)
            family_members[object_]["role"] = "offspring"

        # For parent set the role and add the offspring
        elif relation_type == "parent-of":
            family_members[subject]["role"] = "parent"
            family_members[subject]["offspring"].append(object_)
            if subject not in relations:
                relations[subject] = []
            relations[subject].append(object_)
            family_members[object_]["role"] = "offspring"

    # 3) FIND THE BLOOD TYPE OF EACH FAMILY MEMBER IF EXISTS IN THE TEST RESULTS
    for result in test_results:
        person = result.get("person")
        if person in family_members:
            if result.get("type") == "cheap-bloodtype-test":
                # 20% chance of incorrect result
                if random.random() < 0.2:
                    # Randomly select a different person's blood type from the same country
                    other_person = random.choice([p for p in family_members.keys() if p != person])
                    family_members[person]["bloodtype"] = family_members[other_person]["bloodtype"]
                else:
                    family_members[person]["bloodtype"] = result.get("result")
            else:
                family_members[person]["bloodtype"] = result.get("result")


    '''--------------------------------------------------------------------------------------------'''
    ''''PRINTS FOR DEBUGGING'''
    # Define the Bayesian Network structure
    complete_model = DiscreteBayesianNetwork()  
    if use_country_node:
        # Add Country node and CPD
        complete_model.add_node("Country")
        cpd_country = TabularCPD(variable="Country", variable_card=2, values=[[0.5], [0.5]])
        complete_model.add_cpds(cpd_country)

    # Print the family structure for debugging
    for member, info in family_members.items():
        role = info['role'].upper() if info['role'] else 'Unknown role'
        bloodtype = info['bloodtype'] if info['bloodtype'] else ' '
        print(f"{role}: {member} ({bloodtype})")

    print("\nFAMILY MEMBERS: ", family_members)
    print("\nRELATIONS: ", relations)
    offsprings = [offspring for member, info in family_members.items() for offspring in info["offspring"]]
    print("\nOFFSPRINGS: ", offsprings)

    '''--------------------------------------------------------------------------------------------'''
    ''''CREATE ALLELES AND GENOTYPE FOR EACH FAMILY MEMBER'''
    for member, info in family_members.items():
        print("\nMEMBER: ", member)
        print("INFO: ", info)
        # Add allele 1 and 2 and genotype nodes for each member
        allele1 = f"{member}_Allele1"
        allele2 = f"{member}_Allele2"
        # Add nodes to the model
        complete_model.add_nodes_from([allele1, allele2, f"{member}_Genotype"])

    '''--------------------------------------------------------------------------------------------'''
    '''CREATE CPDs FOR EACH FAMILY MEMBER'''
    # Now add CPDs for each family member
    for member, info in family_members.items():
        allele1 = f"{member}_Allele1"
        allele2 = f"{member}_Allele2"

        if member not in offsprings:
            # Founder: no parents
            if use_country_node:
                cpd_allele1 = TabularCPD(variable=allele1, variable_card=3, evidence=["Country"], evidence_card=[2], values=founder_allele_cpd)
                cpd_allele2 = TabularCPD(variable=allele2, variable_card=3, evidence=["Country"], evidence_card=[2], values=founder_allele_cpd)
                complete_model.add_edge("Country", allele1)
                complete_model.add_edge("Country", allele2)
            else:
                cpd_allele1 = TabularCPD(variable=allele1, variable_card=3, values=country_cpd)
                cpd_allele2 = TabularCPD(variable=allele2, variable_card=3, values=country_cpd)
        else:
            mother = [m for m in family_members.keys() if member in family_members[m]["offspring"] and family_members[m]["role"] == "mother"]
            father = [f for f in family_members.keys() if member in family_members[f]["offspring"] and family_members[f]["role"] == "father"]  
            parent = [p for p in family_members.keys() if member in family_members[p]["offspring"] and family_members[p]["role"] == "parent"]
            # FATHER
            if father:
                cpd_allele1 = TabularCPD(variable=allele1, variable_card=3, evidence=[f"{father[0]}_Genotype"], evidence_card=[6], values=OFFSPIRING_CPD)
            elif parent:
                cpd_allele1 = TabularCPD(variable=allele1, variable_card=3, evidence=[f"{parent[0]}_Genotype"], evidence_card=[6], values=OFFSPIRING_CPD)
            elif use_country_node:
                cpd_allele1 = TabularCPD(variable=allele1, variable_card=3, evidence=["Country"], evidence_card=[2], values=founder_allele_cpd)
                complete_model.add_edge("Country", allele1)
            else:
                cpd_allele1 = TabularCPD(variable=allele1, variable_card=3, values=country_cpd)
            # MOTHER
            if mother:
                cpd_allele2 = TabularCPD(variable=allele2, variable_card=3, evidence=[f"{mother[0]}_Genotype"], evidence_card=[6], values=OFFSPIRING_CPD)
            elif use_country_node:
                cpd_allele2 = TabularCPD(variable=allele2, variable_card=3, evidence=["Country"], evidence_card=[2], values=founder_allele_cpd)
                complete_model.add_edge("Country", allele2)
            else:
                cpd_allele2 = TabularCPD(variable=allele2, variable_card=3, values=country_cpd)

        complete_model.add_cpds(cpd_allele1, cpd_allele2)
        complete_model.add_edges_from([
            (allele1, f"{member}_Genotype"),
            (allele2, f"{member}_Genotype")
        ])

        # add genotype for each member
        cpd_genotype = TabularCPD(
            variable=f"{member}_Genotype",
            variable_card=6,
            evidence=[allele1, allele2],
            evidence_card=[3, 3],
            values=GENOTYPE_CPD
        )
        complete_model.add_cpds(cpd_genotype)

    # Add edges between parents and children
    for person in family_members.keys():
        role, offspring = family_members[person]["role"], family_members[person]["offspring"]
        for s in offspring:
            if role == "father":
                complete_model.add_edge(f"{person}_Genotype", f"{s}_Allele1")
            elif role == "mother":
                complete_model.add_edge(f"{person}_Genotype", f"{s}_Allele2")
            elif role == "parent":
                complete_model.add_edge(f"{person}_Genotype", f"{s}_Allele1")

    queries_persons = [query.get("person") for query in queries]
    # Perform inference
    for member, info in family_members.items():
        # check if has a bloodtype or in query list
        if info["bloodtype"] or member in queries_persons:
            bloodtype_node = f"{member}_Bloodtype"
            complete_model.add_node(bloodtype_node)

            cpd = TabularCPD(variable=bloodtype_node, variable_card=4, evidence=[f"{member}_Genotype"],
                             evidence_card=[6], values=SUM_6_4)
            complete_model.add_cpds(cpd)
            complete_model.add_edge(f"{member}_Genotype", bloodtype_node)

    print("\nNODES: ",complete_model.nodes())
    print("\nEDGES: ", complete_model.edges())

    inference_complete = VariableElimination(complete_model)

    results = []
    for query in queries:
        person = query.get("person")
        if person in family_members:
            inference_variable = f"{person}_Bloodtype"
            evidence = {}
            for member, info in family_members.items():
                if info["bloodtype"]:
                    evidence[f"{member}_Bloodtype"] = ['A', 'B', 'O', 'AB'].index(info["bloodtype"])
            print("\nevidence: ", evidence)
            overall_distribution = inference_complete.query(variables=[inference_variable], evidence=evidence)

            for member in family_members.keys():
                geno = f"{member}_Genotype"
                if geno not in evidence:
                    print(inference_complete.query(variables=[geno], evidence=evidence))
                alleles = [f"{member}_Allele1", f"{member}_Allele2"]
                for allele in alleles:
                    if allele not in evidence:
                        print(inference_complete.query(variables=[allele], evidence=evidence))

            genotype_mapping = {0: "A", 1: "B", 2: "O", 3: "AB"}
            named_result = {genotype_mapping[state]: prob for state, prob in enumerate(overall_distribution.values)}
            result = {
                "type": "bloodtype",
                "person": person,
                "distribution": {
                    "O": round(named_result["O"], 9),
                    "A": round(named_result["A"], 9),
                    "B": round(named_result["B"], 9),
                    "AB": round(named_result["AB"], 9)
                }
            }
            results.append(result)

    # Save results to a JSON file
    output_filename = os.path.join(os.getcwd(), f'p-solutions/solution-{problem_type}-{problem_number:02d}.json')
    os.makedirs(os.path.dirname(output_filename), exist_ok=True)
    with open(output_filename, 'w') as outfile:
        json.dump(results, outfile, indent=4)

def main():
    # Ensure the p-solutions directory exists
    os.makedirs(os.path.join(os.getcwd(), 'p-solutions'), exist_ok=True)

    # Set your desired pattern here (e.g., 'problem-a-*.json')
    pattern = 'problem-e-*.json'
    problem_files = glob.glob(os.path.join(os.getcwd(), 'problems', pattern))

    for problem_file in problem_files:
        try:
            # Extract problem type and number from the filename
            filename = os.path.basename(problem_file)
            problem_type, problem_number = filename.split('-')[1], int(filename.split('-')[2].split('.')[0])
            print(f"\nProcessing problem {problem_number} of type {problem_type}...")
            process_problem(problem_type, problem_number)
        except Exception as e:
            print(f"Error processing problem {problem_number} of type {problem_type}: {e}")
            continue

if __name__ == "__main__":
    main()



"""RUN ALL"""
# def main():
#     # Ensure the p-solutions directory exists
#     os.makedirs(os.path.join(os.getcwd(), 'p-solutions'), exist_ok=True)

#     # Get all problem files in the problems folder
#     problem_files = glob.glob(os.path.join(os.getcwd(), 'problems', '*.json'))

#     for problem_file in problem_files:
#         try:
#             # Extract problem type and number from the filename
#             filename = os.path.basename(problem_file)
#             problem_type, problem_number = filename.split('-')[1], int(filename.split('-')[2].split('.')[0])
#             print(f"\nProcessing problem {problem_number} of type {problem_type}...")
#             process_problem(problem_type, problem_number)
#         except Exception as e:
#             print(f"Error processing problem {problem_number} of type {problem_type}: {e}")
#             continue

# if __name__ == "__main__":
#     main()