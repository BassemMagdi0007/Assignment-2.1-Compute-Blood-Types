import json
import logging
from constants import BloodTypeConstants  
from pgmpy.models import BayesianNetwork
from pgmpy.factors.discrete import TabularCPD
from pgmpy.inference import VariableElimination
import matplotlib.pyplot as plt
import networkx as nx

# ------------------------------
# Utility Functions
# ------------------------------

def load_json(filename):
    """Load a JSON file and return its data."""
    try:
        with open(filename, 'r') as file:
            return json.load(file)
    except Exception as e:
        print(f"Error loading JSON file: {e}")
        return None


def extract_data(data):
    """Extract relevant information from the parsed JSON data."""
    return {
        "family_tree": data.get("family-tree", []),
        "test_results": data.get("test-results", []),
        "queries": data.get("queries", []),
        "country": data.get("country", None)
    }


def add_cpds_to_network(network, country):
    """
    Adds Conditional Probability Distributions (CPDs) to the Bayesian Network.

    Parameters:
    - network (BayesianNetwork): The Bayesian network to which CPDs will be added.
    - country (str): The country for which the blood type distribution is used.

    Returns:
    - bool: True if the network is valid after adding the CPDs, otherwise False.
    """
    try:
        # Validate the country
        if country not in BloodTypeConstants.COUNTRY_BLOOD_TYPE_DISTRIBUTION:
            raise ValueError(f"Country {country} is not in the blood type distribution table.")

        blood_type_distribution = BloodTypeConstants.COUNTRY_BLOOD_TYPE_DISTRIBUTION[country]

        # Convert the distribution to allele probabilities
        allele_probabilities = [
            [blood_type_distribution["A"]], 
            [blood_type_distribution["B"]],
            [blood_type_distribution["O"]],
        ]
        # print("allele_probabilities: ", allele_probabilities)

        # Define CPDs for Father and Mother alleles
        try:
            father_allele1_cpd = TabularCPD("Father_Allele1", 3, allele_probabilities)
            # print("father_allele1_cpd: ", father_allele1_cpd)
        except Exception as e:
            raise ValueError(f"Error defining 'Father_Allele1' CPD: {e}")

        try:
            father_allele2_cpd = TabularCPD("Father_Allele2", 3, allele_probabilities)
            # print("father_allele2_cpd: ", father_allele2_cpd)
        except Exception as e:
            raise ValueError(f"Error defining 'Father_Allele2' CPD: {e}")

        try:
            mother_allele1_cpd = TabularCPD("Mother_Allele1", 3, allele_probabilities)
            # print("mother_allele1_cpd: ", mother_allele1_cpd)
        except Exception as e:
            raise ValueError(f"Error defining 'Mother_Allele1' CPD: {e}")

        try:
            mother_allele2_cpd = TabularCPD("Mother_Allele2", 3, allele_probabilities)
            # print("mother_allele2_cpd: ", mother_allele2_cpd)
        except Exception as e:
            raise ValueError(f"Error defining 'Mother_Allele2' CPD: {e}")

        # Define CPDs for Father and Mother blood types
        bloodtype_values = [
            list(BloodTypeConstants.PARENT_TO_OFFSPRING_GENE[gene].values())
            for gene in BloodTypeConstants.PARENT_TO_OFFSPRING_GENE.keys()
        ]
        # print("bloodtype_values: ", bloodtype_values)

        try:
            father_bloodtype_cpd = TabularCPD(
                variable="Father_BloodType",
                variable_card=9,
                values=bloodtype_values,
                evidence=["Father_Allele1", "Father_Allele2"],
                evidence_card=[3, 3],
            )
            print("father_bloodtype: ", father_bloodtype_cpd)
            
        except Exception as e:
            raise ValueError(f"Error defining 'Father_BloodType' CPD: {e}")

        try:
            mother_bloodtype_cpd = TabularCPD(
                variable="Mother_BloodType",
                variable_card=9,
                values=bloodtype_values,
                evidence=["Mother_Allele1", "Mother_Allele2"],
                evidence_card=[3, 3],
            )
        except Exception as e:
            raise ValueError(f"Error defining 'Mother_BloodType' CPD: {e}")

        # Define CPDs for Offspring Alleles
        try:
            # offspring_allele1_cpd = TabularCPD(
            #     variable="Offspring_Allele1",
            #     variable_card=3,
            #     values=[
            #         list(BloodTypeConstants.GENE_COMBINATION_TO_BLOODTYPE[combination].values())
            #         for combination in BloodTypeConstants.GENE_COMBINATION_TO_BLOODTYPE
            #     ],
            #     evidence=["Father_BloodType"],
            #     evidence_card=[9],
            # )
            offspring_allele1_cpd = TabularCPD(
            variable='Offspring_Allele1',
            variable_card=3,  # Possible alleles: A, B, O
            values=[
                [1, 0.5, 0.5, 0, 0, 0],  # Probabilities for 'A'
                [0, 0, 0, 1, 0.5, 0],    # Probabilities for 'B'
                [0, 0.5, 0.5, 0, 0.5, 1]  # Probabilities for 'O'
            ],
            evidence=['Father_Bloodtype'],
            evidence_card=6  # Blood types: AA, AB, AO, BB, BO, OO
        )
        except Exception as e:
            raise ValueError(f"Error defining 'Offspring_Allele1' CPD: {e}")

        try:
            offspring_allele2_cpd = TabularCPD(
                variable="Offspring_Allele2",
                variable_card=3,
                values=[
                    list(BloodTypeConstants.GENE_COMBINATION_TO_BLOODTYPE[combination].values())
                    for combination in BloodTypeConstants.GENE_COMBINATION_TO_BLOODTYPE
                ],
                evidence=["Mother_BloodType"],
                evidence_card=[9],
            )
        except Exception as e:
            raise ValueError(f"Error defining 'Offspring_Allele2' CPD: {e}")

        # Define CPD for Offspring Blood Type
        try:
            offspring_bloodtype_values = []

            for allele1 in ["A", "B", "O"]:
                for allele2 in ["A", "B", "O"]:
                    gene_combination = allele1 + allele2
                    if gene_combination not in BloodTypeConstants.PARENT_TO_OFFSPRING_GENE:
                        gene_combination = allele2 + allele1  # Handle reverse combinations (e.g., "OA" -> "AO")

                    probabilities = BloodTypeConstants.PARENT_TO_OFFSPRING_GENE.get(gene_combination, {})

                    # For each of the 6 possible offspring blood types, get the probability from the table
                    offspring_bloodtype_values.append(
                        [probabilities.get(bloodtype, 0) for bloodtype in ["AA", "AO", "BB", "BO", "AB", "OO"]]
                    )

            # Convert the bloodtype values to the correct shape (6, 9)
            offspring_bloodtype_values = list(zip(*offspring_bloodtype_values))  # Transpose the matrix to (6, 9)


            offspring_bloodtype_cpd = TabularCPD(
                variable="Offspring_BloodType",
                variable_card=9,  # 6 possible blood types
                values=offspring_bloodtype_values,
                evidence=["Offspring_Allele1", "Offspring_Allele2"],
                evidence_card=[3, 3],  # 3 possible alleles: A, B, O
            )
        except Exception as e:
            raise ValueError(f"Error defining 'Offspring_BloodType' CPD: {e}")

        # Add all CPDs to the network
        network.add_cpds(
            father_allele1_cpd, father_allele2_cpd,
            mother_allele1_cpd, mother_allele2_cpd,
            father_bloodtype_cpd, mother_bloodtype_cpd,
            offspring_allele1_cpd, offspring_allele2_cpd,
            offspring_bloodtype_cpd,
        )

        # Validate the model
        is_valid = network.check_model()
        print(f"Is the model valid after adding CPDs? {is_valid}")
        return is_valid

    except Exception as e:
        print(f"Error while adding CPDs: {e}")
        return False


def bayesian_network(country):
    """
    Constructs and returns the Bayesian Network for blood type inference.
    """
    network = BayesianNetwork()

    # Define nodes
    network.add_nodes_from(["Father_Allele1", "Father_Allele2", "Father_BloodType"])
    network.add_nodes_from(["Mother_Allele1", "Mother_Allele2", "Mother_BloodType"])
    network.add_nodes_from(["Offspring_Allele1", "Offspring_Allele2", "Offspring_BloodType"])

    # Define edges
    network.add_edges_from([("Father_Allele1", "Father_BloodType"), ("Father_Allele2", "Father_BloodType")])
    network.add_edges_from([("Mother_Allele1", "Mother_BloodType"), ("Mother_Allele2", "Mother_BloodType")])
    network.add_edges_from([("Father_BloodType", "Offspring_Allele1"), ("Mother_BloodType", "Offspring_Allele2")])
    network.add_edges_from([("Offspring_Allele1", "Offspring_BloodType"), ("Offspring_Allele2", "Offspring_BloodType")])

    # Add CPDs to the network
    if not add_cpds_to_network(network, country):
        print("Failed to create a valid Bayesian network.")
        return None

    print("Bayesian network created successfully.")
    return network


def visualize_bayesian_network(network):
    """
    Visualizes the Bayesian Network using networkx and matplotlib.

    Parameters:
    - network (BayesianNetwork): The Bayesian network to visualize.
    """
    try:
        # Convert pgmpy's BayesianNetwork to a networkx DiGraph
        graph = nx.DiGraph()
        graph.add_edges_from(network.edges())

        # Draw the graph
        plt.figure(figsize=(10, 8))
        pos = nx.spring_layout(graph, seed=42)  # Define node positions
        nx.draw(
            graph,
            pos,
            with_labels=True,
            node_size=3000,
            node_color="skyblue",
            font_size=10,
            font_weight="bold",
            edge_color="gray"
        )
        plt.title("Bayesian Network Structure")
        plt.show()

    except Exception as e:
        print(f"Error visualizing the Bayesian network: {e}")


# ------------------------------
# Main Execution
# ------------------------------

if __name__ == "__main__":
    json_filename = "example-problems/problem-a-15.json"
    data = load_json(json_filename)
    # print("data: ", data)

    if data:
        extracted_data = extract_data(data)
        country = extracted_data.get("country")

        if country:
            network = bayesian_network(country)
            if network:
                visualize_bayesian_network(network)
        else:
            print("No valid country found in JSON.")
    else:
        print("Failed to load JSON data.")
