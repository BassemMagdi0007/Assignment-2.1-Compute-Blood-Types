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

# Define the CPDs for Allele1 and Allele2 (using North Wamponia as default)
cpd_allele1 = TabularCPD(variable="Allele1", variable_card=3, values=cpd_north_wamponia)
cpd_allele2 = TabularCPD(variable="Allele2", variable_card=3, values=cpd_north_wamponia)

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

# Add CPDs to the network
genetic_model.add_cpds(cpd_allele1, cpd_allele2, cpd_genotype)

# Check if the model is valid
assert genetic_model.check_model(), "The model is invalid!"

# Perform inference on the model
inference = VariableElimination(genetic_model)

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
    return named_result

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
    return named_result

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

# Example Queries
# Query 1: Overall genotype distribution
print("Overall Genotype Distribution:")
overall_distribution = calculate_genotype_distribution()
for genotype, prob in overall_distribution.items():
    print(f"{genotype}: {prob:.4f}")

# Example: Querying genotype distribution given specific alleles
print("\nGenotype Distribution given Alleles A and B:")
distribution_given_alleles = calculate_genotype_given_alleles("A", "B")
for genotype, prob in distribution_given_alleles.items():
    print(f"{genotype}: {prob:.4f}")

# Example: Update CPDs for Allele1 as bloodtype A and Allele2 as North Wamponia and query the genotype distribution
update_cpds("bloodtype_a", "north_wamponia")
print("\nOverall Genotype Distribution for Allele1 as bloodtype A and Allele2 as North Wamponia:")
overall_distribution_mixed = calculate_genotype_distribution()
for genotype, prob in overall_distribution_mixed.items():
    print(f"{genotype}: {prob:.4f}")