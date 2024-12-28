from pgmpy.models import BayesianNetwork
from pgmpy.factors.discrete import TabularCPD
from pgmpy.inference import VariableElimination
import networkx as nx
import matplotlib.pyplot as plt

# Define the Bayesian Network structure
complete_model = BayesianNetwork(
    [
        # Parent1 structure
        ("Parent1_Allele1", "Parent1_Genotype"),
        ("Parent1_Allele2", "Parent1_Genotype"),
        
        # Parent2 structure
        ("Parent2_Allele1", "Parent2_Genotype"),
        ("Parent2_Allele2", "Parent2_Genotype"),
        
        # Offspring structure dependent on Parent1 and Parent2
        ("Parent1_Genotype", "Offspring_Allele1"),
        ("Parent2_Genotype", "Offspring_Allele2"),
        ("Offspring_Allele1", "Offspring_Genotype"),
        ("Offspring_Allele2", "Offspring_Genotype"),
    ]
)

# Define CPDs for Parent1
cpd_allele1_parent1 = TabularCPD(variable="Parent1_Allele1", variable_card=3, values=[[0.5], [0.25], [0.25]])
cpd_allele2_parent1 = TabularCPD(variable="Parent1_Allele2", variable_card=3, values=[[0.5], [0.25], [0.25]])
cpd_genotype_parent1 = TabularCPD(
    variable="Parent1_Genotype",
    variable_card=4,
    evidence=["Parent1_Allele1", "Parent1_Allele2"],
    evidence_card=[3, 3],
    values=[
        # AA,  AB,  AO,  BA,  BB,  BO,  OA,  OB,  OO
        [1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0],  # A
        [0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 1.0, 0.0],  # B
        [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0],  # O
        [0.0, 1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0],  # AB
    ],
)

# Define CPDs for Parent2
cpd_allele1_parent2 = TabularCPD(variable="Parent2_Allele1", variable_card=3, values=[[0.5], [0.25], [0.25]])
cpd_allele2_parent2 = TabularCPD(variable="Parent2_Allele2", variable_card=3, values=[[0.5], [0.25], [0.25]])
cpd_genotype_parent2 = TabularCPD(
    variable="Parent2_Genotype",
    variable_card=4,
    evidence=["Parent2_Allele1", "Parent2_Allele2"],
    evidence_card=[3, 3],
    values=[
        # AA,  AB,  AO,  BA,  BB,  BO,  OA,  OB,  OO
        [1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0],  # A
        [0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 1.0, 0.0],  # B
        [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0],  # O
        [0.0, 1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0],  # AB
    ],
)

# Define CPDs for Offspring
cpd_allele1_offspring = TabularCPD(
    variable="Offspring_Allele1",
    variable_card=3,
    values=[
        # A ,  B ,  O ,  AB
        [1.0, 0.0, 0.0, 0.5],  # A
        [0.0, 1.0, 0.0, 0.5],  # B
        [0.0, 0.0, 1.0, 0.0],  # O
        
    ],
    evidence=["Parent1_Genotype"],
    evidence_card=[4],
)

cpd_allele2_offspring = TabularCPD(
    variable="Offspring_Allele2",
    variable_card=3,
    values=[
        # A ,  B ,  O ,  AB
        [1.0, 0.0, 0.0, 0.5],  # A
        [0.0, 1.0, 0.0, 0.5],  # B
        [0.0, 0.0, 1.0, 0.0],  # O
    ],
    evidence=["Parent2_Genotype"],
    evidence_card=[4],
)

cpd_genotype_offspring = TabularCPD(
    variable="Offspring_Genotype",
    variable_card=4,
    values=[
        # AA,  AB,  AO,  BA,  BB,  BO,  OA,  OB,  OO
        [1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0],  # A
        [0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 1.0, 0.0],  # B
        [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0],  # O
        [0.0, 1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0],  # AB
    ],
    evidence=["Offspring_Allele1", "Offspring_Allele2"],
    evidence_card=[3, 3],
)

# Add all CPDs to the unified model
complete_model.add_cpds(
    cpd_allele1_parent1,
    cpd_allele2_parent1,
    cpd_genotype_parent1,
    cpd_allele1_parent2,
    cpd_allele2_parent2,
    cpd_genotype_parent2,
    cpd_allele1_offspring,
    cpd_allele2_offspring,
    cpd_genotype_offspring,
)

# Validate the combined model
assert complete_model.check_model(), "The combined Bayesian Network is invalid!"

# Inference on the unified model
inference_complete = VariableElimination(complete_model)

# Example: Querying Offspring's genotype distribution
print("Overall Genotype Distribution for Offspring (from unified model):")
overall_distribution_offspring = inference_complete.query(variables=["Offspring_Genotype"])
genotype_mapping = {0: "A", 1: "B", 2: "O", 3: "AB"}
named_result_offspring = {genotype_mapping[state]: prob for state, prob in enumerate(overall_distribution_offspring.values)}
for genotype, prob in named_result_offspring.items():
    print(f"{genotype}: {prob:.4f}")

# Visualization of the unified Bayesian Network
plt.figure(figsize=(12, 8))
G = nx.DiGraph()
G.add_edges_from(complete_model.edges())
pos = nx.spring_layout(G)
nx.draw(
    G,
    pos,
    with_labels=True,
    node_size=3000,
    node_color="lightgreen",
    font_size=10,
    font_weight="bold",
    arrowsize=20,
)
plt.title("Unified Bayesian Network for Parents and Offspring")
plt.show()
