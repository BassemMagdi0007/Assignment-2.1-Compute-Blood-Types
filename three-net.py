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

# Predefined values
# Mapping from indices to allele names for printing
allele_mapping = {"A": 0, "B": 1, "O": 2}
reverse_allele_mapping = {0: "A", 1: "B", 2: "O"}

# Define the CPDs for different scenarios
cpd_north_wamponia = [[0.5], [0.25], [0.25]]
cpd_south_wamponia = [[0.15], [0.55], [0.30]]

cpd_bloodtype_a = [[0.75], [0.0], [0.25]]
cpd_bloodtype_b = [[0.0], [0.75], [0.25]]
cpd_bloodtype_o = [[0.0], [0.0], [1.0]]
cpd_bloodtype_ab = [[0.5], [0.5], [0.0]]

# Load and extract data from JSON file
filename = "../../example-problems/problem-a-00.json"
data = load_json(filename)
if data:
    extracted_data = extract_data(data) #dict
    family_tree = extracted_data["family_tree"] #list
    test_results = extracted_data["test_results"] #list
    queries = extracted_data["queries"] #list
    country = extracted_data["country"] #str

# TODO: Should we ignore the mother entirely in case of father blood type requested? 
# Dynamically define the family members (father, mother, offspring)
family_members = {"father": None, "mother": None, "offspring": []}

# Process family tree to identify the relations
for relation in family_tree:
    subject = relation["subject"]
    object_ = relation["object"]
    
    if relation["relation"] == "father-of":
        # FIXME: is it right to igmnore father from offspring list?
        # Ensure no redundancy; ignore if the subject is already assigned as mother or offspring
        if family_members["father"] is None and subject not in family_members["offspring"]:
            family_members["father"] = subject
        if object_ not in family_members["offspring"]:
            family_members["offspring"].append(object_)

    elif relation["relation"] == "mother-of":
        # Ensure no redundancy; ignore if the subject is already assigned as father or offspring
        if family_members["mother"] is None and subject not in family_members["offspring"]:
            family_members["mother"] = subject
        if object_ not in family_members["offspring"]:
            family_members["offspring"].append(object_)

# Find the blood type for the father
father_bloodtype = None
for result in test_results:
    # Compare the name of the person in father 
    if result.get("person") == family_members["father"]:
        father_bloodtype = result.get("result")
        break

# Find the blood type for the mother
mother_bloodtype = None
for result in test_results:
    if result.get("person") == family_members["mother"]:
        mother_bloodtype = result.get("result")
        break

# Print the family structure for debugging
father_info = f"{family_members['father']} ({father_bloodtype if father_bloodtype else 'Unknown'})"
mother_info = f"{family_members['mother']} ({mother_bloodtype if mother_bloodtype else 'Unknown'})"
offspring_info = [f"{child} (Unknown)" for child in family_members["offspring"]]
print(f"Father: {father_info}")
print(f"Mother: {mother_info}")
print(f"Offspring: {offspring_info}")

# Define the Bayesian Network structure
complete_model = BayesianNetwork(
    [
        # Father structure
        ("father_Allele1", "father_Genotype"),
        ("father_Allele2", "father_Genotype"),

        # Mother structure
        ("mother_Allele1", "mother_Genotype"),
        ("mother_Allele2", "mother_Genotype"),

        # Offspring structure dependent on father and mother
        ("father_Genotype", "Offspring_Allele1"),
        ("mother_Genotype", "Offspring_Allele2"),
        ("Offspring_Allele1", "Offspring_Genotype"),
        ("Offspring_Allele2", "Offspring_Genotype"),
    ]
)

# if father bloodtype known assign its cpd value directly to offspring allele
# TODO: cpd values for father_allele1 and father_allele2 and father_genotype must be assigned? 
# Can they be ignored? 
# Can I store them in a place and ignore it? 
# Better: make the network dynamic to not include father_allele1 and father_allele2 and father_genotype in case of known bloodtype!
if father_bloodtype:
    if father_bloodtype == "A":
        cpd_allele1_offspring = TabularCPD(variable="Offspring_Allele1", variable_card=3, values=cpd_bloodtype_a)
    elif father_bloodtype == "B":
        cpd_allele1_offspring = TabularCPD(variable="Offspring_Allele1", variable_card=3, values=cpd_bloodtype_b)
    elif father_bloodtype == "O":
        cpd_allele1_offspring = TabularCPD(variable="Offspring_Allele1", variable_card=3, values=cpd_bloodtype_o)
    elif father_bloodtype == "AB":
        cpd_allele1_offspring = TabularCPD(variable="Offspring_Allele1", variable_card=3, values=cpd_bloodtype_ab)
else:
    if country == "North Wamponia":
        cpd_allele1_father = TabularCPD(variable="father_Allele1", variable_card=3, values=cpd_north_wamponia)
        cpd_allele2_father = TabularCPD(variable="father_Allele2", variable_card=3, values=cpd_north_wamponia)
    elif country == "South Wamponia":
        cpd_allele1_father = TabularCPD(variable="father_Allele1", variable_card=3, values=cpd_south_wamponia)
        cpd_allele2_father = TabularCPD(variable="father_Allele2", variable_card=3, values=cpd_south_wamponia)
    else:
        # Default to North Wamponia if country is not specified
        cpd_allele1_father = TabularCPD(variable="father_Allele1", variable_card=3, values=cpd_north_wamponia)
        cpd_allele2_father = TabularCPD(variable="father_Allele2", variable_card=3, values=cpd_north_wamponia)

    cpd_genotype_father = TabularCPD(
        variable="father_Genotype",
        variable_card=4,
        evidence=["father_Allele1", "father_Allele2"],
        evidence_card=[3, 3],
        values=[
            # AA, AB, AO, BA, BB, BO, OA, OB, OO
            [1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0],  # A
            [0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 1.0, 0.0],  # B
            [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0],  # O
            [0.0, 1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0],  # AB
        ],
    )

# Assign CPDs for mother based on country and blood type
if mother_bloodtype:
    if mother_bloodtype == "A":
        cpd_allele2_offspring = TabularCPD(variable="Offspring_Allele2", variable_card=3, values=cpd_bloodtype_a)
    elif mother_bloodtype == "B":
        cpd_allele2_offspring = TabularCPD(variable="Offspring_Allele2", variable_card=3, values=cpd_bloodtype_b)
    elif mother_bloodtype == "O":
        cpd_allele2_offspring = TabularCPD(variable="Offspring_Allele2", variable_card=3, values=cpd_bloodtype_o)
    elif mother_bloodtype == "AB":
        cpd_allele2_offspring = TabularCPD(variable="Offspring_Allele2", variable_card=3, values=cpd_bloodtype_ab)
else:
    if country == "North Wamponia":
        cpd_allele1_mother = TabularCPD(variable="mother_Allele1", variable_card=3, values=cpd_north_wamponia)
        cpd_allele2_mother = TabularCPD(variable="mother_Allele2", variable_card=3, values=cpd_north_wamponia)
    elif country == "South Wamponia":
        cpd_allele1_mother = TabularCPD(variable="mother_Allele1", variable_card=3, values=cpd_south_wamponia)
        cpd_allele2_mother = TabularCPD(variable="mother_Allele2", variable_card=3, values=cpd_south_wamponia)
    else:
        # Default to North Wamponia if country is not specified
        cpd_allele1_mother = TabularCPD(variable="mother_Allele1", variable_card=3, values=cpd_north_wamponia)
        cpd_allele2_mother = TabularCPD(variable="mother_Allele2", variable_card=3, values=cpd_north_wamponia)

    cpd_genotype_mother = TabularCPD(
        variable="mother_Genotype",
        variable_card=4,
        evidence=["mother_Allele1", "mother_Allele2"],
        evidence_card=[3, 3],
        values=[
            # AA, AB, AO, BA, BB, BO, OA, OB, OO
            [1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0],  # A
            [0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 1.0, 0.0],  # B
            [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0],  # O
            [0.0, 1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0],  # AB
        ],
    )

# Define CPDs for Offspring
if not father_bloodtype:
    cpd_allele1_offspring = TabularCPD(
        variable="Offspring_Allele1",
        variable_card=3,
        values=[
            # A ,  B ,  O ,  AB
            [1.0, 0.0, 0.0, 0.5],  # A
            [0.0, 1.0, 0.0, 0.5],  # B
            [0.0, 0.0, 1.0, 0.0],  # O
        ],
        evidence=["father_Genotype"],
        evidence_card=[4],
    )

if not mother_bloodtype:
    cpd_allele2_offspring = TabularCPD(
        variable="Offspring_Allele2",
        variable_card=3,
        values=[
            # A ,  B ,  O ,  AB
            [1.0, 0.0, 0.0, 0.5],  # A
            [0.0, 1.0, 0.0, 0.5],  # B
            [0.0, 0.0, 1.0, 0.0],  # O
        ],
        evidence=["mother_Genotype"],
        evidence_card=[4],
    )

cpd_genotype_offspring = TabularCPD(
    variable="Offspring_Genotype",
    variable_card=4,
    values=[
        # AA, AB, AO, BA, BB, BO, OA, OB, OO
        [1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0],  # A
        [0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 1.0, 0.0],  # B
        [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0],  # O
        [0.0, 1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0],  # AB
    ],
    evidence=["Offspring_Allele1", "Offspring_Allele2"],
    evidence_card=[3, 3],
)

# Add all CPDs to the unified model
if not father_bloodtype:
    complete_model.add_cpds(cpd_allele1_father, cpd_allele2_father, cpd_genotype_father)
if not mother_bloodtype:
    complete_model.add_cpds(cpd_allele1_mother, cpd_allele2_mother, cpd_genotype_mother)
complete_model.add_cpds(cpd_allele1_offspring, cpd_allele2_offspring, cpd_genotype_offspring)

# Validate the combined model
assert complete_model.check_model(), "The combined Bayesian Network is invalid!"

# Inference on the unified model
inference_complete = VariableElimination(complete_model)

# Example: Querying Offspring's genotype distribution
print("\nOverall Genotype Distribution for Offspring (from unified model):")
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