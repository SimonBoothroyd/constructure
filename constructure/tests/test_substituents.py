import pytest
from rdkit import Chem

from constructure.substituents import SUBSTITUENTS


@pytest.mark.parametrize("substituent", SUBSTITUENTS.values())
def test_default_scaffolds(substituent):

    # Make sure the SMILES pattern is valid and contains a single R group.
    substituent = substituent.replace("[R]", "[1*]")
    substituent_molecule = Chem.MolFromSmiles(substituent)

    assert (
        sum(1 for atom in substituent_molecule.GetAtoms() if atom.HasProp("dummyLabel"))
        == 1
    )
