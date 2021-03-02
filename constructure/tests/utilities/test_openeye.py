import os

import pytest

from constructure.utilities import MissingOptionalDependency
from constructure.utilities.openeye import (
    remove_duplicate_smiles,
    requires_oe_package,
    smiles_to_image_grid,
)

pytest.importorskip("openeye.oechem")


@requires_oe_package("oechem")
def dummy_oe_function():
    return 5


def test_requires_oe_package(monkeypatch):

    from openeye import oechem

    monkeypatch.setattr(oechem, "OEChemIsLicensed", lambda: False)

    with pytest.raises(MissingOptionalDependency) as error_info:
        dummy_oe_function()

    assert error_info.value.library_name == "openeye.oechem"
    assert error_info.value.license_issue is True


def test_smiles_to_image_grid(tmpdir):

    smiles_to_image_grid(["C"], os.path.join(tmpdir, "c.png"))
    assert os.path.isfile(os.path.join(tmpdir, "c.png"))


def test_remove_duplicate_smiles():

    unique_smiles = remove_duplicate_smiles(["C(Cl)(F)(Br)", "C(F)(Cl)(Br)"])

    assert len(unique_smiles) == 1
