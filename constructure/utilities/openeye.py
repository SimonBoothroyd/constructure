"""A set of utilities to aid in interfacing with the OpenEye toolkits."""
import functools
import logging
import math
from typing import List, Optional, TypeVar

from typing_extensions import Literal

from constructure.utilities.utilities import MissingOptionalDependency, requires_package

logger = logging.getLogger(__name__)

T = TypeVar("T")


def requires_oe_package(
    package_name: Literal["oechem", "oedepict"]
):  # pragma: no cover
    """A decorator which checks that the required OpenEye package is installed
    and licensed.

    Args:
        package_name: The name of the required OpenEye package.
    """

    def inner_decorator(function):
        @requires_package(f"openeye.{package_name}")
        @functools.wraps(function)
        def wrapper(*args, **kwargs):

            from openeye import oechem, oedepict

            if (package_name == "oechem" and not oechem.OEChemIsLicensed()) or (
                package_name == "oedepict" and not oedepict.OEDepictIsLicensed()
            ):

                raise MissingOptionalDependency(f"openeye.{package_name}", True)

            return function(*args, **kwargs)

        return wrapper

    return inner_decorator


@requires_oe_package("oechem")
@requires_oe_package("oedepict")
def smiles_to_image_grid(
    smiles: List[str],
    output_path: str,
    labels: Optional[List[str]] = None,
    cols: int = 8,
    cell_width: int = 200,
    cell_height: int = 200,
):  # pragma: no cover
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

    from openeye import oechem, oedepict

    rows = math.ceil(len(smiles) / cols)

    image = oedepict.OEImage(cell_width * cols, cell_height * rows)
    grid = oedepict.OEImageGrid(image, rows, cols)

    opts = oedepict.OE2DMolDisplayOptions(
        grid.GetCellWidth(), grid.GetCellHeight(), oedepict.OEScale_AutoScale
    )
    opts.SetAromaticStyle(oedepict.OEAromaticStyle_Circle)
    opts.SetTitleLocation(oedepict.OETitleLocation_Bottom)

    for i, (smi, cell) in enumerate(zip(smiles, grid.GetCells())):

        smi = smi if labels is None else f"{smi} {labels[i]}"

        mol = oechem.OEGraphMol()
        oechem.OESmilesToMol(mol, smi)

        oedepict.OEPrepareDepiction(mol)

        disp = oedepict.OE2DMolDisplay(mol, opts)
        oedepict.OERenderMolecule(cell, disp)

    oedepict.OEWriteImage(output_path, image)


@requires_oe_package("oechem")
def remove_duplicate_smiles(smiles: List[str]) -> List[str]:  # pragma: no cover
    """Returns the list of unique SMILES patterns in a specified list.

    Args:
        smiles: The original list of SMILES patterns

    Returns:
        The unique SMILES patterns.
    """

    from openeye import oechem

    # Use a separate closed list to retain the original ordering.
    unique_smiles = []
    closed_list = set()

    smiles_options = (
        oechem.OESMILESFlag_Canonical
        | oechem.OESMILESFlag_Isotopes
        | oechem.OESMILESFlag_RGroups
        | oechem.OESMILESFlag_AtomStereo
        | oechem.OESMILESFlag_BondStereo
        | oechem.OESMILESFlag_Hydrogens
    )

    for smiles_pattern in smiles:

        oe_molecule = oechem.OEMol()
        oechem.OESmilesToMol(oe_molecule, smiles_pattern)

        unique_pattern = oechem.OECreateSmiString(oe_molecule, smiles_options)

        if unique_pattern not in closed_list:
            unique_smiles.append(unique_pattern)

        closed_list.add(unique_pattern)

    return unique_smiles
