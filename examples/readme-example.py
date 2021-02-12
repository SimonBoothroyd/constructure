"""A runnable script containing the minimal code from the README example."""


def main():

    # Retrieve the default imine scaffold.
    from constructure.scaffolds import SCAFFOLDS

    scaffold = SCAFFOLDS["imine"]

    # Import a constructor and use it to enumerate all combinations of a set of
    # substituents attached to a scaffold.
    from constructure.constructors import RDKitConstructor as Constructor
    from constructure.substituents import SUBSTITUENTS

    smiles = Constructor.enumerate_combinations(
        scaffold,
        substituents={
            1: [SUBSTITUENTS["methyl"], SUBSTITUENTS["phenyl"]],
            2: [SUBSTITUENTS["isopropyl"], SUBSTITUENTS["ethyl"]],
            3: [SUBSTITUENTS["hydrogen"], SUBSTITUENTS["ethyl"]],
        },
    )

    # Save the 2D structures to an image file.
    from constructure.utilities.rdkit import smiles_to_image_grid

    smiles_to_image_grid(smiles, "imines.png", cols=4)


if __name__ == "__main__":
    main()
