from typing import List, Optional

from constructure.utilities import requires_package


@requires_package("rdkit")
def smiles_to_image_grid(
    smiles: List[str],
    output_path: str,
    labels: Optional[List[str]] = None,
    cols: int = 8,
    cell_width: int = 200,
    cell_height: int = 200,
):
    """Saves a list of SMILES patterns as an image of their corresponding 2D structures.

    Args:
        smiles: The list of SMILES patterns to include in the image.
        output_path: The path to save the image to.
        labels: A set of labels (one per smiles pattern) to show underneath each 2D
            structure.
        cols: The number of 2D structures per row.
        cell_width: The width (in px) to draw each 2D structure.
        cell_height: The height (in px) to draw each 2D structure.
    """

    from rdkit import Chem
    from rdkit.Chem import Draw

    molecules = [Chem.MolFromSmiles(pattern) for pattern in smiles]

    image = Draw.MolsToGridImage(molecules, cols, (cell_width, cell_height), labels)
    image.save(output_path)


@requires_package("rdkit")
def remove_duplicate_smiles(smiles: List[str]) -> List[str]:
    """Returns the list of unique SMILES patterns in a specified list.

    Args:
        smiles: The original list of SMILES patterns

    Returns:
        The unique SMILES patterns.
    """

    from rdkit import Chem

    # Use a separate closed list to retain the original ordering.
    unique_smiles = []
    closed_list = set()

    for smiles_pattern in smiles:

        rd_molecule = Chem.MolFromSmiles(smiles_pattern)
        unique_pattern = Chem.MolToSmiles(
            rd_molecule, isomericSmiles=True, allHsExplicit=True
        )

        if unique_pattern not in closed_list:
            unique_smiles.append(unique_pattern)

        closed_list.add(unique_pattern)

    return unique_smiles
