import pytest

from constructure.constructors import OpenEyeConstructor, RDKitConstructor
from constructure.scaffolds import Scaffold
from constructure.tests import does_not_raise


@pytest.mark.parametrize("constructor", [OpenEyeConstructor, RDKitConstructor])
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
def test_n_replaceable_groups(constructor, scaffold, expected):
    assert constructor.n_replaceable_groups(scaffold) == expected


@pytest.mark.parametrize("constructor", [OpenEyeConstructor, RDKitConstructor])
def test_attach_substituents(constructor):

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


@pytest.mark.parametrize("constructor", [OpenEyeConstructor, RDKitConstructor])
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
        ("[R]O", None),
    ],
)
def test_classify_substituent(constructor, substituent, expected):
    assert constructor.classify_substituent(substituent) == expected


@pytest.mark.parametrize("constructor", [OpenEyeConstructor, RDKitConstructor])
def test_classify_substituent_error(constructor):

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
        OpenEyeConstructor.validate_substituents(scaffold, substituents)


@pytest.mark.parametrize("constructor", [OpenEyeConstructor, RDKitConstructor])
def test_remove_duplicate_smiles(constructor):

    unique_smiles = constructor._remove_duplicate_smiles(
        ["C(Cl)(F)(Br)", "C(F)(Cl)(Br)"]
    )

    assert len(unique_smiles) == 1


@pytest.mark.parametrize("constructor", [OpenEyeConstructor, RDKitConstructor])
def test_enumerate_combinations_combinatorial(constructor):

    from rdkit import Chem

    scaffold = Scaffold(
        smiles="C([R1])C(O)CC([R2])", r_groups={1: ["alkyl"], 2: ["acyl"]}
    )

    enumerated_smiles = constructor.enumerate_combinations(
        scaffold=scaffold,
        substituents={1: ["[R]C", "[R]CC"], 2: ["[R]C=O", "[R]C(=O)C"]},
        mode="combinatorial",
    )

    assert len(enumerated_smiles) == 4

    assert {"CCCC(O)CCC=O", "CCC(O)CCC=O", "CCC(O)CCC(C)=O", "CCCC(O)CCC(C)=O"} == {
        Chem.MolToSmiles(Chem.MolFromSmiles(smiles)) for smiles in enumerated_smiles
    }
