import json
from pgmpy.factors.discrete import TabularCPD
from pgmpy.models import BayesianNetwork
from pgmpy.inference import VariableElimination

# Mapping from indices to allele names
allele_mapping = {"A": 0, "B": 1, "O": 2}

# Define the CPDs for different scenarios
cpd_north_wamponia = [[0.5], [0.25], [0.25]]
cpd_south_wamponia = [[0.15], [0.55], [0.30]]

cpd_bloodtype_a = [[0.75], [0.0], [0.25]]
cpd_bloodtype_b = [[0.0], [0.75], [0.25]]
cpd_bloodtype_o = [[0.0], [0.0], [1.0]]
cpd_bloodtype_ab = [[0.5], [0.5], [0.0]]

# Define the network structure
genetic_model = BayesianNetwork(
    [
        ("Allele1", "Genotype"),  # Genotype is dependent on Allele1 and Allele2
        ("Allele2", "Genotype"),  # The two Alleles are independent of each other
    ]
)

# Define the CPD for Genotype based on Allele1 and Allele2
cpd_genotype = TabularCPD(
    variable="Genotype",
    variable_card=4,  # A, B, O, AB
    evidence=["Allele1", "Allele2"],
    evidence_card=[3, 3],  # A, B, O (for each allele)
    values=[
        # AA ,  AB ,   AO ,  BA ,  BB   ,  BO   ,  OA  ,   OB  ,  OO
        [1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0],  # A
        [0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 1.0, 0.0],  # B
        [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0],  # O
        [0.0, 1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0],  # AB
    ],
)

# Add the CPD for Genotype to the network
genetic_model.add_cpds(cpd_genotype)

# Load a JSON file and return its data
def load_json(filename):
    try:
        with open(filename, 'r') as file:
            return json.load(file)
    except Exception as e:
        print(f"Error loading JSON file: {e}")
        return None

# Extract relevant information from the parsed JSON data
def extract_data(data):
    return {
        "family_tree": data.get("family-tree", []),
        "test_results": data.get("test-results", []),
        "queries": data.get("queries", []),
        "country": data.get("country", None)
    }

# Function to update CPDs based on the scenario
def update_cpds(scenario_allele1, scenario_allele2):
    if scenario_allele1 == "north_wamponia":
        cpd_allele1 = TabularCPD(variable="Allele1", variable_card=3, values=cpd_north_wamponia)
    elif scenario_allele1 == "south_wamponia":
        cpd_allele1 = TabularCPD(variable="Allele1", variable_card=3, values=cpd_south_wamponia)
    elif scenario_allele1 == "bloodtype_a":
        cpd_allele1 = TabularCPD(variable="Allele1", variable_card=3, values=cpd_bloodtype_a)
    elif scenario_allele1 == "bloodtype_b":
        cpd_allele1 = TabularCPD(variable="Allele1", variable_card=3, values=cpd_bloodtype_b)
    elif scenario_allele1 == "bloodtype_o":
        cpd_allele1 = TabularCPD(variable="Allele1", variable_card=3, values=cpd_bloodtype_o)
    elif scenario_allele1 == "bloodtype_ab":
        cpd_allele1 = TabularCPD(variable="Allele1", variable_card=3, values=cpd_bloodtype_ab)
    else:
        raise ValueError("Unknown scenario for Allele1")

    if scenario_allele2 == "north_wamponia":
        cpd_allele2 = TabularCPD(variable="Allele2", variable_card=3, values=cpd_north_wamponia)
    elif scenario_allele2 == "south_wamponia":
        cpd_allele2 = TabularCPD(variable="Allele2", variable_card=3, values=cpd_south_wamponia)
    elif scenario_allele2 == "bloodtype_a":
        cpd_allele2 = TabularCPD(variable="Allele2", variable_card=3, values=cpd_bloodtype_a)
    elif scenario_allele2 == "bloodtype_b":
        cpd_allele2 = TabularCPD(variable="Allele2", variable_card=3, values=cpd_bloodtype_b)
    elif scenario_allele2 == "bloodtype_o":
        cpd_allele2 = TabularCPD(variable="Allele2", variable_card=3, values=cpd_bloodtype_o)
    elif scenario_allele2 == "bloodtype_ab":
        cpd_allele2 = TabularCPD(variable="Allele2", variable_card=3, values=cpd_bloodtype_ab)
    else:
        raise ValueError("Unknown scenario for Allele2")

    genetic_model.add_cpds(cpd_allele1, cpd_allele2)
    assert genetic_model.check_model(), "The model is invalid!"

# Query functions
def calculate_genotype_distribution():
    """
    Calculate the overall genotype distribution given the allele probabilities.
    """
    result = inference.query(variables=["Genotype"])
    
    # Map the result to genotype names
    genotype_mapping = {0: "A", 1: "B", 2: "O", 3: "AB"}
    named_result = {
        genotype_mapping[state]: prob
        for state, prob in enumerate(result.values)
    }
    
    # Order the result as "O", "A", "B", "AB"
    ordered_result = {key: named_result[key] for key in ["O", "A", "B", "AB"]}
    return ordered_result

def calculate_genotype_given_alleles(allele1, allele2):
    """
    Calculate the probability of a genotype given specific alleles.
    """
    evidence = {"Allele1": allele_mapping[allele1], "Allele2": allele_mapping[allele2]}
    result = inference.query(variables=["Genotype"], evidence=evidence)
    
    # Map the result to genotype names
    genotype_mapping = {0: "A", 1: "B", 2: "O", 3: "AB"}
    named_result = {
        genotype_mapping[state]: prob
        for state, prob in enumerate(result.values)
    }
    
    # Order the result as "O", "A", "B", "AB"
    ordered_result = {key: named_result[key] for key in ["O", "A", "B", "AB"]}
    return ordered_result

# Main function to process the JSON file and calculate the genotype distribution
def main():
    json_file = "../example-problems/problem-a-02.json"  # Replace with your actual file
    data = load_json(json_file)
    if data is None:
        return

    extracted_data = extract_data(data)
    family_tree = extracted_data["family_tree"]
    test_results = extracted_data["test_results"]
    queries = extracted_data["queries"]
    country = extracted_data["country"]

    # Get all distinct names from the JSON file
    names = set()
    for relation in family_tree:
        names.add(relation["subject"])
        names.add(relation["object"])
    for test in test_results:
        names.add(test["person"])
    for query in queries:
        names.add(query["person"])

    # Determine the object and subjects based on the queries
    if queries:
        genotype_object = queries[0]["person"]
        names.remove(genotype_object)
        allele1_subject, allele2_subject = names
    else:
        print("No queries found in the JSON file.")
        return

    # Determine the scenarios for Allele1 and Allele2
    scenario_allele1 = None
    scenario_allele2 = None

    for test in test_results:
        if test["person"] == allele1_subject:
            if test["result"] == "A":
                scenario_allele1 = "bloodtype_a"
            elif test["result"] == "B":
                scenario_allele1 = "bloodtype_b"
            elif test["result"] == "O":
                scenario_allele1 = "bloodtype_o"
            elif test["result"] == "AB":
                scenario_allele1 = "bloodtype_ab"
        elif test["person"] == allele2_subject:
            if test["result"] == "A":
                scenario_allele2 = "bloodtype_a"
            elif test["result"] == "B":
                scenario_allele2 = "bloodtype_b"
            elif test["result"] == "O":
                scenario_allele2 = "bloodtype_o"
            elif test["result"] == "AB":
                scenario_allele2 = "bloodtype_ab"

    if scenario_allele1 is None:
        if country == "North Wumponia":
            scenario_allele1 = "north_wamponia"
        elif country == "South Wumponia":
            scenario_allele1 = "south_wamponia"

    if scenario_allele2 is None:
        if country == "North Wumponia":
            scenario_allele2 = "north_wamponia"
        elif country == "South Wumponia":
            scenario_allele2 = "south_wamponia"

    # Update the CPDs based on the scenarios
    update_cpds(scenario_allele1, scenario_allele2)

    # Perform inference on the model
    global inference
    inference = VariableElimination(genetic_model)

    # Print the names of the subjects and object
    print(f"\nSubject 1 (Allele1): {allele1_subject}")
    print(f"Subject 2 (Allele2): {allele2_subject}")
    print(f"Object (Genotype): {genotype_object}")

    # Calculate the overall genotype distribution
    overall_distribution = calculate_genotype_distribution()
    print("\nOverall Genotype Distribution:")
    for genotype, prob in overall_distribution.items():
        print(f"{genotype}: {prob:.4f}")

if __name__ == "__main__":
    main()