import os

from constructure.utilities.rdkit import remove_duplicate_smiles, smiles_to_image_grid


def test_smiles_to_image_grid(tmpdir):

    smiles_to_image_grid(["C"], os.path.join(tmpdir, "c.png"))
    assert os.path.isfile(os.path.join(tmpdir, "c.png"))


def test_remove_duplicate_smiles():

    unique_smiles = remove_duplicate_smiles(["C(Cl)(F)(Br)", "C(F)(Cl)(Br)"])

    assert len(unique_smiles) == 1
