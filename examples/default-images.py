"""A simple script for generating images of the default scaffolds and
substituents included by this framework.
"""
from constructure.scaffolds import SCAFFOLDS
from constructure.substituents import SUBSTITUENTS
from constructure.utilities.openeye import smiles_to_image_grid


def main():

    scaffold_smiles = []
    scaffold_labels = []

    for label, scaffold in SCAFFOLDS.items():
        scaffold_smiles.append(scaffold.smiles)
        scaffold_labels.append(label)

    smiles_to_image_grid(scaffold_smiles, "scaffolds.png", scaffold_labels, cols=6)

    substituent_smiles = []
    substituent_labels = []

    for label, substituent in SUBSTITUENTS.items():
        substituent_smiles.append(substituent.replace("[R]", "[R1]-"))
        substituent_labels.append(label)

    smiles_to_image_grid(
        substituent_smiles, "substituents.png", substituent_labels, cols=6
    )


if __name__ == "__main__":
    main()
