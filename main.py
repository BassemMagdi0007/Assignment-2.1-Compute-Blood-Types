import json
import logging
from pgmpy.models import BayesianNetwork
from pgmpy.factors.discrete import TabularCPD
from pgmpy.inference import VariableElimination
import networkx as nx
import matplotlib.pyplot as plt

# Suppress pgmpy warnings
logging.getLogger("pgmpy").setLevel(logging.ERROR)

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

# Calculate conditional CPD based on blood type and country CPD
def calculate_conditional_cpd(bloodtype, country_cpd):
    # Extract probabilities from the country CPD
    p_A = country_cpd[0][0]
    p_B = country_cpd[1][0]
    p_O = country_cpd[2][0]

    if bloodtype == "A":
        p_AA = p_A * p_A
        p_AO = 2 * (p_A * p_O)
        p_A_total = p_AA + p_AO
        p_AA_given_A = p_AA / p_A_total
        p_AO_given_A = p_AO / p_A_total
        return [[p_AA_given_A + 0.5 * p_AO_given_A], [0.0], [0.5 * p_AO_given_A]]
    elif bloodtype == "B":
        p_BB = p_B * p_B
        p_BO = 2 * (p_B * p_O)
        p_B_total = p_BB + p_BO
        p_BB_given_B = p_BB / p_B_total
        p_BO_given_B = p_BO / p_B_total
        return [[0.0], [p_BB_given_B + 0.5 * p_BO_given_B], [0.5 * p_BO_given_B]]
    elif bloodtype == "O":
        return [[0.0], [0.0], [1.0]]
    elif bloodtype == "AB":
        return [[0.5], [0.5], [0.0]]
    else:
        raise ValueError("Invalid blood type")

def process_problem(problem_type, problem_number):
    # Conditional Probability Distributions (CPDs) for the alleles and genotypes
    cpd_north_wamponia = [[0.5], [0.25], [0.25]]
    cpd_south_wamponia = [[0.15], [0.55], [0.30]]

    # Load and extract data from JSON file
    filename = f'example-problems/problem-{problem_type}-{problem_number:02d}.json'
    data = load_json(filename)
    if not data:
        print(f"Skipping problem {problem_number} due to loading error.")
        return

    extracted_data = extract_data(data)
    family_tree = extracted_data["family_tree"]
    test_results = extracted_data["test_results"]
    queries = extracted_data["queries"]
    country = extracted_data["country"]

    # Define country_cpd based on the country
    if country == "North Wamponia":
        country_cpd = cpd_north_wamponia
    elif country == "South Wamponia":
        country_cpd = cpd_south_wamponia
    else:
        country_cpd = cpd_north_wamponia  # Default to North Wamponia if the country is not recognized

    # Dynamically define the family members and their relations
    family_members = {}
    relations = {}

    for key in family_tree:
        subject = key["subject"]
        object_ = key["object"]
        relation_type = key["relation"]

        # Initialize family members if not already present
        if subject not in family_members:
            family_members[subject] = {"role": None, "bloodtype": None, "offspring": []}
        if object_ not in family_members:
            family_members[object_] = {"role": None, "bloodtype": None, "offspring": []}

        # Update roles and relationships
        if relation_type == "father-of":
            family_members[subject]["role"] = "father"
            family_members[subject]["offspring"].append(object_)
            if subject not in relations:
                relations[subject] = []
            relations[subject].append(object_)
            # Set object_ as an offspring
            family_members[object_]["role"] = "offspring"
        elif relation_type == "mother-of":
            family_members[subject]["role"] = "mother"
            family_members[subject]["offspring"].append(object_)
            if subject not in relations:
                relations[subject] = []
            relations[subject].append(object_)
            # Set object_ as an offspring
            family_members[object_]["role"] = "offspring"



    # Find the blood type for each family member
    for result in test_results:
        person = result.get("person")
        if person in family_members:
            family_members[person]["bloodtype"] = result.get("result")

    # Print the family structure for debugging
    for member, info in family_members.items():
        role = info['role'].capitalize() if info['role'] else 'Unknown role'
        bloodtype = info['bloodtype'] if info['bloodtype'] else 'Unknown blood type'
        print(f"{role}: {member} ({bloodtype})")

    # Define the Bayesian Network structure
    complete_model = BayesianNetwork()

    # Create nodes and CPDs for each family member
    for member, info in family_members.items():
        if info["role"] == "father":
            if any(family_members[child]["bloodtype"] for child in info["offspring"]):
                first_offspring_bloodtype = next(family_members[child]["bloodtype"] for child in info["offspring"] if family_members[child]["bloodtype"])
                cpd_allele1 = TabularCPD(variable=f"{member}_Allele1", variable_card=3, values=calculate_conditional_cpd(first_offspring_bloodtype, country_cpd))
            else:
                cpd_allele1 = TabularCPD(variable=f"{member}_Allele1", variable_card=3, values=country_cpd)
            cpd_allele2 = TabularCPD(variable=f"{member}_Allele2", variable_card=3, values=country_cpd)
            cpd_genotype = TabularCPD(
                variable=f"{member}_Genotype",
                variable_card=4,
                evidence=[f"{member}_Allele1", f"{member}_Allele2"],
                evidence_card=[3, 3],
                values=[
                    [1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0],  # A
                    [0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 1.0, 0.0],  # B
                    [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0],  # O
                    [0.0, 1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0],  # AB
                ],
            )
            complete_model.add_edges_from([
                (f"{member}_Allele1", f"{member}_Genotype"),
                (f"{member}_Allele2", f"{member}_Genotype")
            ])
            complete_model.add_cpds(cpd_allele1, cpd_allele2, cpd_genotype)

        elif info["role"] == "mother":
            if any(family_members[child]["bloodtype"] for child in info["offspring"]):
                first_offspring_bloodtype = next(family_members[child]["bloodtype"] for child in info["offspring"] if family_members[child]["bloodtype"])
                cpd_allele1 = TabularCPD(variable=f"{member}_Allele1", variable_card=3, values=calculate_conditional_cpd(first_offspring_bloodtype, country_cpd))
            else:
                cpd_allele1 = TabularCPD(variable=f"{member}_Allele1", variable_card=3, values=country_cpd)
            cpd_allele2 = TabularCPD(variable=f"{member}_Allele2", variable_card=3, values=country_cpd)
            cpd_genotype = TabularCPD(
                variable=f"{member}_Genotype",
                variable_card=4,
                evidence=[f"{member}_Allele1", f"{member}_Allele2"],
                evidence_card=[3, 3],
                values=[
                    [1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0],  # A
                    [0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 1.0, 0.0],  # B
                    [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0],  # O
                    [0.0, 1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0],  # AB
                ],
            )
            complete_model.add_edges_from([
                (f"{member}_Allele1", f"{member}_Genotype"),
                (f"{member}_Allele2", f"{member}_Genotype")
            ])
            complete_model.add_cpds(cpd_allele1, cpd_allele2, cpd_genotype)

    # Create nodes and CPDs for offspring based on their parents
    for father, children in relations.items():
        for child in children:
            mother = next((m for m, c in relations.items() if child in c), None)
            if mother:
                # Determine CPDs for offspring alleles based on parents' blood types
                father_bloodtype = family_members[father]["bloodtype"]
                mother_bloodtype = family_members[mother]["bloodtype"]

                if not father_bloodtype:
                    cpd_allele1 = TabularCPD(variable=f"{child}_Allele1", variable_card=3, values=country_cpd)
                else:
                    cpd_allele1 = TabularCPD(variable=f"{child}_Allele1", variable_card=3, values=calculate_conditional_cpd(father_bloodtype, country_cpd))

                if not mother_bloodtype:
                    cpd_allele2 = TabularCPD(variable=f"{child}_Allele2", variable_card=3, values=country_cpd)
                else:
                    cpd_allele2 = TabularCPD(variable=f"{child}_Allele2", variable_card=3, values=calculate_conditional_cpd(mother_bloodtype, country_cpd))

                cpd_genotype = TabularCPD(
                    variable=f"{child}_Genotype",
                    variable_card=4,
                    evidence=[f"{child}_Allele1", f"{child}_Allele2"],
                    evidence_card=[3, 3],
                    values=[
                        [1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0],  # A
                        [0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 1.0, 0.0],  # B
                        [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0],  # O
                        [0.0, 1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0],  # AB
                    ],
                )
                complete_model.add_edges_from([
                    (f"{child}_Allele1", f"{child}_Genotype"),
                    (f"{child}_Allele2", f"{child}_Genotype")
                ])
                complete_model.add_cpds(cpd_allele1, cpd_allele2, cpd_genotype)

    # Perform inference
    inference_complete = VariableElimination(complete_model)

    for query in queries:
        person = query.get("person")
        if person in family_members:
            inference_variable = f"{person}_Genotype"
            print(f"Genotype Distribution for {person}:")
            overall_distribution = inference_complete.query(variables=[inference_variable])
            genotype_mapping = {0: "A", 1: "B", 2: "O", 3: "AB"}
            named_result = {genotype_mapping[state]: prob for state, prob in enumerate(overall_distribution.values)}
            print(f"O: {named_result['O']:.4f}")
            print(f"A: {named_result['A']:.4f}")
            print(f"B: {named_result['B']:.4f}")
            print(f"AB: {named_result['AB']:.4f}")

    # # Visualization of the unified Bayesian Network
    # plt.figure(figsize=(12, 8))
    # G = nx.DiGraph()
    # G.add_edges_from(complete_model.edges())
    # pos = nx.spring_layout(G)
    # nx.draw(
    #     G,
    #     pos,
    #     with_labels=True,
    #     node_size=3000,
    #     node_color="lightgreen",
    #     font_size=10,
    #     font_weight="bold",
    #     arrowsize=20,
    # )
    # plt.title("Unified Bayesian Network for Parents and Offspring")
    # plt.show()

def main():
    # Define the type of problems
    problem_type = 'a'
    # problem_number = 0

    # # Process problems from 0 to 14
    for problem_number in range(15):
        print(f"\nProcessing problem {problem_number}...")
        process_problem(problem_type, problem_number)

if __name__ == "__main__":
    main()