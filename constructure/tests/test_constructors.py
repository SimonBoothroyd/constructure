import importlib

import pytest

from constructure.constructors import Constructor, OpenEyeConstructor, RDKitConstructor
from constructure.scaffolds import Scaffold
from constructure.tests import does_not_raise

try:
    importlib.import_module("openeye")
    CONSTRUCTORS = [OpenEyeConstructor, RDKitConstructor]
except ImportError:
    CONSTRUCTORS = [RDKitConstructor]


@pytest.mark.parametrize("constructor", CONSTRUCTORS)
@pytest.mark.parametrize(
    "scaffold, expected",
    [
        (Scaffold(smiles="C([R1])", r_groups={1: ["hydrogen"]}), 1),
        (
            Scaffold(
                smiles="C([R1])([R2])", r_groups={i: ["hydrogen"] for i in range(2)}
            ),
            2,
        ),
        (
            Scaffold(
                smiles="C([R1])([R2])([R3])",
                r_groups={i: ["hydrogen"] for i in range(3)},
            ),
            3,
        ),
        (
            Scaffold(
                smiles="C([R1])([R2])([R3])([R4])",
                r_groups={i: ["hydrogen"] for i in range(4)},
            ),
            4,
        ),
    ],
)
def test_n_replaceable_groups(constructor: Constructor, scaffold, expected):
    assert constructor.n_replaceable_groups(scaffold) == expected


@pytest.mark.parametrize("constructor", CONSTRUCTORS)
@pytest.mark.parametrize(
    "scaffold, expected",
    [
        (Scaffold(smiles="C([R1])", r_groups={1: ["hydrogen"]}), [1]),
        (
            Scaffold(
                smiles="C([R1])([R2])", r_groups={i: ["hydrogen"] for i in range(2)}
            ),
            [1, 2],
        ),
        (
            Scaffold(
                smiles="C([R1])([R2])([R3])",
                r_groups={i: ["hydrogen"] for i in range(3)},
            ),
            [1, 2, 3],
        ),
        (
            Scaffold(
                smiles="C([R1])([R2])([R3])([R4])",
                r_groups={i: ["hydrogen"] for i in range(4)},
            ),
            [1, 2, 3, 4],
        ),
        (
            Scaffold(
                smiles="C1(=C([C]([R4])=C2C(=[C]([R7])1)[N]([C]([R2])=[C]([R3])2)[H])O)O",
                r_groups={i: ["hydrogen"] for i in (2, 3, 4, 7)},
            ),
            [2, 3, 4, 7],
        ),
    ],
)
def test_get_replaceable_r_groups(constructor: Constructor, scaffold, expected):
    assert constructor.get_replaceable_r_groups(scaffold) == expected


@pytest.mark.parametrize("constructor", CONSTRUCTORS)
def test_attach_substituents(constructor: Constructor):

    from rdkit import Chem

    scaffold = Scaffold(
        smiles="C([R1])([R2])([R3])([R4])",
        r_groups={1: ["hydrogen"], 2: ["alkyl"], 3: ["aryl"], 4: ["halogen"]},
    )

    smiles = constructor.attach_substituents(
        scaffold, {1: "[R][H]", 2: "[R]C", 3: "[R]c1ccccc1", 4: "[R]Cl"}
    )
    smiles = Chem.MolToSmiles(Chem.MolFromSmiles(smiles))

    assert smiles == "CC(Cl)c1ccccc1"


@pytest.mark.parametrize("constructor", CONSTRUCTORS)
@pytest.mark.parametrize(
    "substituent, expected",
    [
        ("[R][H]", "hydrogen"),
        ("[H][R]", "hydrogen"),
        ("[R]C", "alkyl"),
        ("C[R]", "alkyl"),
        ("[R]c1ccccc1", "aryl"),
        ("c1ccc([R])cc1", "aryl"),
        ("[R]C(=O)", "acyl"),
        ("C([R])(=O)", "acyl"),
        ("[R][C+](=O)", "acyl"),
        ("[R]Cl", "halogen"),
        ("Cl[R]", "halogen"),
        ("[R]O", "hetero"),
        ("[R]N", "hetero"),
        ("[R][Li]", None),
    ],
)
def test_classify_substituent(constructor: Constructor, substituent, expected):
    assert constructor.classify_substituent(substituent) == expected


@pytest.mark.parametrize("constructor", CONSTRUCTORS)
def test_classify_substituent_error(constructor: Constructor):

    with pytest.raises(ValueError, match="The substituent C does not contain the "):
        constructor.classify_substituent("C")


@pytest.mark.parametrize(
    "substituents, expected_raises",
    [
        (
            {1: ["[R][H]"], 2: ["[R]C"], 3: ["[R]c1ccccc1"], 4: ["[R]Cl"]},
            does_not_raise(),
        ),
        (
            {1: ["[R][H]"], 3: ["[R]c1ccccc1"]},
            pytest.raises(
                ValueError, match="Please provide substituents for the R2, R4 groups"
            ),
        ),
        (
            {1: [], 2: ["[R]C"], 3: ["[R]c1ccccc1"], 4: ["[R]Cl"]},
            pytest.raises(
                ValueError, match="Please provide substituents for the R1 groups"
            ),
        ),
        (
            {
                1: ["[R][H]"],
                2: ["[R]C"],
                3: ["[R]c1ccccc1"],
                4: ["[R]Cl"],
                5: ["[R]Cl"],
            },
            pytest.raises(
                ValueError, match="The scaffold does not contain R5 groups for"
            ),
        ),
        (
            {1: ["[R][H]"], 2: ["[R]C=O"], 3: ["[R]c1ccccc1"], 4: ["[R]Cl"]},
            pytest.raises(ValueError, match="The R2 group only accepts "),
        ),
    ],
)
def test_validate_substituents(substituents, expected_raises):

    scaffold = Scaffold(
        smiles="C([R1])([R2])([R3])([R4])",
        r_groups={1: ["hydrogen"], 2: ["alkyl"], 3: ["aryl"], 4: ["halogen"]},
    )

    with expected_raises:
        RDKitConstructor.validate_substituents(scaffold, substituents)


def test_validate_replaceable_r_groups():
    scaffold = Scaffold(
        smiles="C([R1])([R1])([R3])([R4])",
        r_groups={1: ["hydrogen"], 2: ["alkyl"], 3: ["aryl"], 4: ["halogen"]},
    )
    err = "Duplicate R-group values found"
    with pytest.raises(ValueError, match=err):
        RDKitConstructor.validate_substituents(
            scaffold, {1: ["[R][H]"], 2: ["[R]C"], 3: ["[R]c1ccccc1"], 4: ["[R]Cl"]}
        )


@pytest.mark.parametrize("constructor", CONSTRUCTORS)
def test_remove_duplicate_smiles(constructor: Constructor):

    unique_smiles = constructor._remove_duplicate_smiles(
        ["C(Cl)(F)(Br)", "C(F)(Cl)(Br)"]
    )

    assert len(unique_smiles) == 1


@pytest.mark.parametrize("constructor", CONSTRUCTORS)
@pytest.mark.parametrize(
    "smiles, r_groups, substituents",
    [
        (
            "C([R1])C(O)CC([R2])",
            {1: ["alkyl"], 2: ["acyl"]},
            {1: ["[R]C", "[R]CC"], 2: ["[R]C=O", "[R]C(=O)C"]},
        ),
        (
            "C([R9])C(O)CC([R32])",
            {9: ["alkyl"], 32: ["acyl"]},
            {9: ["[R]C", "[R]CC"], 32: ["[R]C=O", "[R]C(=O)C"]},
        ),
    ],
)
def test_enumerate_combinations_combinatorial(
    constructor: Constructor, smiles: str, r_groups: dict, substituents: dict
):

    from rdkit import Chem

    scaffold = Scaffold(
        smiles=smiles,
        r_groups=r_groups,
    )

    enumerated_smiles = constructor.enumerate_combinations(
        scaffold=scaffold,
        substituents=substituents,
        mode="combinatorial",
    )

    assert len(enumerated_smiles) == 4

    assert {"CCCC(O)CCC=O", "CCC(O)CCC=O", "CCC(O)CCC(C)=O", "CCCC(O)CCC(C)=O"} == {
        Chem.MolToSmiles(Chem.MolFromSmiles(smiles)) for smiles in enumerated_smiles
    }


@pytest.mark.parametrize("constructor", CONSTRUCTORS)
@pytest.mark.parametrize(
    "validate, expected_raises",
    [(False, None), (True, pytest.raises(ValueError, match="group only accepts "))],
)
def test_enumerate_combinations_validation(
    constructor: Constructor, validate: bool, expected_raises
):

    if expected_raises is None:
        expected_raises = does_not_raise()

    scaffold = Scaffold(smiles="C([R1])", r_groups={1: ["alkyl"]})

    with expected_raises:

        constructor.enumerate_combinations(
            scaffold=scaffold, substituents={1: ["[R]O"]}, validate=validate
        )
