import json
from pgmpy.models import BayesianNetwork
from pgmpy.factors.discrete import TabularCPD
from pgmpy.inference import VariableElimination
import networkx as nx
import matplotlib.pyplot as plt

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

    # Dynamically define the family members (father, mother, offspring)
    family_members = {"father": None, "mother": None, "offspring": []}

    # Process family tree to identify the relations
    for relation in family_tree:
        subject = relation["subject"]
        object_ = relation["object"]
        # Ensure no redundancy; ignore if the subject is already assigned as father or offspring
        if relation["relation"] == "father-of":
            if family_members["father"] is None and subject not in family_members["offspring"]:
                family_members["father"] = subject
            if object_ not in family_members["offspring"]:
                family_members["offspring"].append(object_)

        elif relation["relation"] == "mother-of":
            if family_members["mother"] is None and subject not in family_members["offspring"]:
                family_members["mother"] = subject
            if object_ not in family_members["offspring"]:
                family_members["offspring"].append(object_)

    # Find the blood type for the father
    father_bloodtype = None
    for result in test_results:
        if result.get("person") == family_members["father"]:
            # print("FATHER BLOODTYPE DETECTED")
            father_bloodtype = result.get("result")
            break

    # Find the blood type for the mother
    mother_bloodtype = None
    for result in test_results:
        if result.get("person") == family_members["mother"]:
            # print("MOTHER BLOODTYPE DETECTED")
            mother_bloodtype = result.get("result")
            break

    # Find the blood type for the offspring
    offspring_bloodtype = {}
    for result in test_results:
        for child in family_members["offspring"]:
            if result.get("person") == child:
                # print(f"OFFSPRING BLOODTYPE DETECTED for {child}")
                offspring_bloodtype[child] = result.get("result")
            break

    # Print the family structure for debugging
    father_info = f"{family_members['father']} ({father_bloodtype if father_bloodtype else 'Unknown'})"
    mother_info = f"{family_members['mother']} ({mother_bloodtype if mother_bloodtype else 'Unknown'})"
    offspring_info = [f"{child} ({offspring_bloodtype if offspring_bloodtype else 'Unknown'})" for child in family_members["offspring"]]
    print(f"Father: {father_info}")
    print(f"Mother: {mother_info}")
    print(f"Offspring: {offspring_info}")

    # Define the Bayesian Network structure
    complete_model = BayesianNetwork()

    # Check if any query is related to offspring
    query_related_to_offspring = False
    query_related_to_father = False
    query_related_to_mother = False
    for query in queries:
        if query.get("person") in family_members["offspring"]:
            query_related_to_offspring = True
        elif query.get("person") == family_members["father"]:
            query_related_to_father = True
        elif query.get("person") == family_members["mother"]:
            query_related_to_mother = True

    # If query is related to offspring, assign CPDs directly to offspring alleles
    if query_related_to_offspring:
        if not father_bloodtype:
            cpd_allele1_offspring = TabularCPD(variable="Offspring_Allele1", variable_card=3, values=country_cpd)
        else:
            cpd_allele1_offspring = TabularCPD(variable="Offspring_Allele1", variable_card=3, values=calculate_conditional_cpd(father_bloodtype, country_cpd))

        if not mother_bloodtype:
            cpd_allele2_offspring = TabularCPD(variable="Offspring_Allele2", variable_card=3, values=country_cpd)
        else:
            cpd_allele2_offspring = TabularCPD(variable="Offspring_Allele2", variable_card=3, values=calculate_conditional_cpd(mother_bloodtype, country_cpd))

        # Define CPDs for Offspring Genotype
        cpd_genotype_offspring = TabularCPD(
            variable="Offspring_Genotype",
            variable_card=4,
            values=[
                [1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0],  # A
                [0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 1.0, 0.0],  # B
                [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0],  # O
                [0.0, 1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0],  # AB
            ],
            evidence=["Offspring_Allele1", "Offspring_Allele2"],
            evidence_card=[3, 3],
        )

        # Add nodes and edges for offspring
        complete_model.add_edges_from([
            ("Offspring_Allele1", "Offspring_Genotype"),
            ("Offspring_Allele2", "Offspring_Genotype")
        ])
        # Add all CPDs to the unified model
        complete_model.add_cpds(cpd_allele1_offspring, cpd_allele2_offspring, cpd_genotype_offspring)

    # If query is related to father, assign CPDs for father alleles and genotype
    if query_related_to_father:
        if offspring_bloodtype:
            first_offspring_bloodtype = next(iter(offspring_bloodtype.values()))
            cpd_allele1_father = TabularCPD(variable="father_Allele1", variable_card=3, values=calculate_conditional_cpd(first_offspring_bloodtype, country_cpd))
        else:
            cpd_allele1_father = TabularCPD(variable="father_Allele1", variable_card=3, values=country_cpd)
        cpd_allele2_father = TabularCPD(variable="father_Allele2", variable_card=3, values=country_cpd)

        cpd_genotype_father = TabularCPD(
            variable="father_Genotype",
            variable_card=4,
            evidence=["father_Allele1", "father_Allele2"],
            evidence_card=[3, 3],
            values=[
                [1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0],  # A
                [0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 1.0, 0.0],  # B
                [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0],  # O
                [0.0, 1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0],  # AB
            ],
        )

        complete_model.add_edges_from([
            ("father_Allele1", "father_Genotype"),
            ("father_Allele2", "father_Genotype")
        ])

        complete_model.add_cpds(cpd_allele1_father, cpd_allele2_father, cpd_genotype_father)

    if query_related_to_mother:
        if offspring_bloodtype:
            first_offspring_bloodtype = next(iter(offspring_bloodtype.values()))
            cpd_allele1_mother = TabularCPD(variable="mother_Allele1", variable_card=3, values=calculate_conditional_cpd(first_offspring_bloodtype, country_cpd))
        else:
            cpd_allele1_mother = TabularCPD(variable="mother_Allele1", variable_card=3, values=country_cpd)
        cpd_allele2_mother = TabularCPD(variable="mother_Allele2", variable_card=3, values=country_cpd)

        cpd_genotype_mother = TabularCPD(
            variable="mother_Genotype",
            variable_card=4,
            evidence=["mother_Allele1", "mother_Allele2"],
            evidence_card=[3, 3],
            values=[
                [1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0],  # A
                [0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 1.0, 0.0],  # B
                [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0],  # O
                [0.0, 1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0],  # AB
            ],
        )

        complete_model.add_edges_from([
            ("mother_Allele1", "mother_Genotype"),
            ("mother_Allele2", "mother_Genotype")
        ])

        complete_model.add_cpds(cpd_allele1_mother, cpd_allele2_mother, cpd_genotype_mother)

    # for cpd in complete_model.get_cpds():
        # print(f"CPD of {cpd.variable}:")
        # print(cpd)

    inference_complete = VariableElimination(complete_model)

    inference_variable = None
    if query_related_to_father:
        inference_variable = "father_Genotype"
    elif query_related_to_mother:
        inference_variable = "mother_Genotype"
    elif query_related_to_offspring:
        inference_variable = "Offspring_Genotype"

    if inference_variable:
        # Print the genotype distribution in the specified order
        print(f"Genotype Distribution for {inference_variable.split('_')[0]}:")
        overall_distribution = inference_complete.query(variables=[inference_variable])
        genotype_mapping = {0: "A", 1: "B", 2: "O", 3: "AB"}
        named_result = {genotype_mapping[state]: prob for state, prob in enumerate(overall_distribution.values)}
        print(f"O: {named_result['O']:.4f}")
        print(f"A: {named_result['A']:.4f}")
        print(f"B: {named_result['B']:.4f}")
        print(f"AB: {named_result['AB']:.4f}")

def main():
    # Define the type of problems
    problem_type = 'a'

    # Process problems from 0 to 14
    for problem_number in range(15):
        print(f"\nProcessing problem {problem_number}...")
        process_problem(problem_type, problem_number)

if __name__ == "__main__":
    main()