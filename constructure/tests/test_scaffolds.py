import pytest

from constructure.constructors import RDKitConstructor
from constructure.scaffolds import SCAFFOLDS


@pytest.mark.parametrize("scaffold", SCAFFOLDS.values())
def test_default_scaffolds(scaffold):

    # Make sure the number of R groups in the `smiles` pattern matches the `r_groups`
    # attributes.
    assert RDKitConstructor.n_replaceable_groups(scaffold) == len(scaffold.r_groups)
    assert {*scaffold.r_groups} == {i + 1 for i in range(len(scaffold.r_groups))}
