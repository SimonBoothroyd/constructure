import abc
import itertools
import re
from typing import Dict, List, Literal, Optional

from constructure.scaffolds import Scaffold
from constructure.utilities import requires_package
from constructure.utilities.openeye import requires_oe_package


class Constructor(abc.ABC):
    """A class which provides methods for enumerating the possible attachments of
    substituents to different molecular scaffolds."""

    @classmethod
    @abc.abstractmethod
    def n_replaceable_groups(cls, scaffold: Scaffold) -> int:
        """Returns the number of R groups on a scaffold

        Args:
            scaffold: The scaffold with R groups defined.

        Returns:
            The number of found R groups.
        """
        raise NotImplementedError()

    @classmethod
    @abc.abstractmethod
    def attach_substituents(
        cls, scaffold: Scaffold, substituents: Dict[int, str]
    ) -> str:
        """Replaces the R groups on a specified scaffold with a set of substituents.

        Args:
            scaffold: The scaffold with R groups defined.
            substituents: A dictionary of substituents to add to the scaffold, where
                each key corresponds to a labelled R group (e.g. 1 corresponds to R1)
                and each value the SMILES pattern of the corresponding substituent to
                attach.

        Returns:
            The SMILES pattern of the substituted molecule.
        """
        raise NotImplementedError()

    @classmethod
    @abc.abstractmethod
    def classify_substituent(
        cls, substituent: str
    ) -> Optional[Literal["hydrogen", "aryl", "alkyl", "acyl", "halogen", "hetero"]]:
        """A function which will attempt to classify a substituent as being either
        a hydrogen, aryl, alkyl, acyl or a halogen.

        Args:
            substituent: The SMILES pattern representing the substituent.

        Returns:
            The classification if one could be determined, otherwise ``None``.
        """
        raise NotImplementedError()

    @classmethod
    def validate_substituents(
        cls,
        scaffold: Scaffold,
        substituents: Dict[int, List[str]],
    ):
        """Attempts to validate that a set of substituents are suitable for attachment
        to a given scaffold based off the the scaffolds defined R groups and the
        classification of each substituent.

        Args:
           scaffold: The scaffold with R groups defined.
           substituents: A dictionary of substituents to add to the scaffold, where
                each key corresponds to a labelled R group (e.g. 1 corresponds to R1)
                and each value a list of SMILES patterns corresponding to possible
                substituents to attach.
        """

        # Begin by checking that the user has specified a set of substitutes for each
        # R groups on the scaffold.
        n_expected_groups = cls.n_replaceable_groups(scaffold)

        missing_r_groups = {i + 1 for i in range(n_expected_groups)} - {*substituents}
        missing_r_groups.update(i for i in substituents if len(substituents[i]) == 0)

        if len(missing_r_groups) > 0:
            missing_r_groups_string = ", ".join(f"R{i}" for i in missing_r_groups)

            raise ValueError(
                f"Please provide substituents for the {missing_r_groups_string} "
                f"groups."
            )

        # Make sure the user didn't provide substituents for non-existent R substituents.
        extra_r_groups = {*substituents} - {i + 1 for i in range(n_expected_groups)}

        if len(extra_r_groups) > 0:
            extra_r_groups_string = ", ".join(f"R{i}" for i in extra_r_groups)

            raise ValueError(
                f"The scaffold does not contain {extra_r_groups_string} groups for "
                f"which a set of substituents were provided."
            )

        # Make sure each provided substituent is valid for the R group - i.e. check
        # the substituent is alkyl, aryl, acyl or a halogen.
        for i in substituents:

            for substituent in substituents[i]:

                substituent_type = cls.classify_substituent(substituent)

                if substituent_type in scaffold.r_groups[i]:
                    continue

                expected_substituents = ", ".join(scaffold.r_groups[i])

                raise ValueError(
                    f"The R{i} group only accepts {expected_substituents} "
                    f"substituents. {substituent} is an {substituent_type} "
                    f"substituent."
                )

    @classmethod
    @abc.abstractmethod
    def _remove_duplicate_smiles(cls, smiles: List[str]) -> List[str]:
        """Returns the list of unique SMILES patterns in a specified list.

        Args:
            smiles: The original list of SMILES patterns

        Returns:
            The unique SMILES patterns.
        """

        raise NotImplementedError()

    @classmethod
    def enumerate_combinations(
        cls,
        scaffold: Scaffold,
        substituents: Dict[int, List[str]],
        mode: Literal["combinatorial"] = "combinatorial",
    ) -> List[str]:
        """Attempts to enumerate the possible ways of attaching a set of substituents to
        a specified scaffold.

        Args:
            scaffold: The scaffold with R groups defined.
            substituents: A dictionary of substituents to add to the scaffold, where
                each key corresponds to a labelled R group (e.g. 1 corresponds to R1)
                and each value a list of SMILES patterns corresponding to possible
                substituents to attach.
            mode: The mode in which to attach the substituents.

                ``'combinatorial'``: Enumerates all combinatorial ways of attaching
                the substituents to the scaffold. The number of these can grow large
                for scaffolds with several R groups and multiple possible substituents
                per R group.

        Returns
            A list of SMILES patterns representing the possible decorated scaffolds.
        """

        # Make sure the substituents are appropriate for the specified scaffold.
        cls.validate_substituents(scaffold, substituents)

        if mode == "combinatorial":

            combinations = itertools.product(
                *(substituents[i + 1] for i in range(len(substituents)))
            )

            group_combinations = [
                {i + 1: substituent for i, substituent in enumerate(combination)}
                for combination in combinations
            ]

        else:
            raise NotImplementedError()

        smiles = [
            cls.attach_substituents(scaffold, group_combination)
            for group_combination in group_combinations
        ]

        return cls._remove_duplicate_smiles(smiles)


class RDKitConstructor(Constructor):
    """A class which provides methods for enumerating the possible attachments of
    substituents to different molecular scaffolds using RDKit."""

    @classmethod
    @requires_package("rdkit")
    def n_replaceable_groups(cls, scaffold: Scaffold) -> int:

        from rdkit import Chem

        scaffold_smiles = re.sub(r"\(\[R([1-9])+]\)", r"([\1*])", scaffold.smiles)
        scaffold_molecule = Chem.MolFromSmiles(scaffold_smiles)

        return sum(
            1 for atom in scaffold_molecule.GetAtoms() if atom.HasProp("dummyLabel")
        )

    @classmethod
    @requires_package("rdkit")
    def attach_substituents(
        cls, scaffold: Scaffold, substituents: Dict[int, str]
    ) -> str:

        from rdkit import Chem
        from rdkit.Chem import rdChemReactions

        scaffold_smiles = re.sub(r"\(\[R([1-9])+]\)", r"([\1*])", scaffold.smiles)
        reactant = (Chem.MolFromSmiles(scaffold_smiles),)

        for i in range(len(substituents)):

            substituent = substituents[i + 1].replace("[R]", "[*:1]")

            rxn = rdChemReactions.ReactionFromSmarts(f"[*:1]-[{i + 1}*]>>{substituent}")
            products = rxn.RunReactants(reactant)

            assert len(products) == 1 and len(products[0]) == 1
            reactant = products[0]

        assert len(reactant) == 1
        return Chem.MolToSmiles(reactant[0])

    @classmethod
    @requires_package("rdkit")
    def classify_substituent(
        cls, substituent: str
    ) -> Optional[Literal["hydrogen", "aryl", "alkyl", "acyl", "halogen", "hetero"]]:

        from rdkit import Chem

        if substituent.find("[R]") < 0 or substituent.find("[R]") != substituent.rfind(
            "[R]"
        ):

            raise ValueError(
                f"The substituent {substituent} does not contain the expected [R] "
                f"attachment point."
            )

        rdkit_molecule = Chem.MolFromSmiles(substituent.replace("[R]", "[1*]"))

        # Find the dummy atom.
        dummy_atom = [
            atom for atom in rdkit_molecule.GetAtoms() if atom.HasProp("dummyLabel")
        ][0]

        assert len(dummy_atom.GetNeighbors()) == 1

        attachment_atom = dummy_atom.GetNeighbors()[0]

        # Check for hydogren
        if attachment_atom.GetAtomicNum() == 1:
            return "hydrogen"

        # Check for halogen
        if attachment_atom.GetAtomicNum() in [9, 17, 35, 53]:
            return "halogen"

        # Check for hetero
        if attachment_atom.GetAtomicNum() in [7, 8, 15, 16]:
            return "hetero"

        # Check for alkyl
        if (
            attachment_atom.GetAtomicNum() == 6
            and attachment_atom.GetTotalDegree() == 4
        ):
            return "alkyl"

        # Check for aryl
        if attachment_atom.GetAtomicNum() in [6, 7] and attachment_atom.GetIsAromatic():
            return "aryl"

        # Check for acyl
        if attachment_atom.GetAtomicNum() == 6 and any(
            neighbor.GetAtomicNum() == 8 and neighbor.GetTotalDegree() == 1
            for neighbor in attachment_atom.GetNeighbors()
        ):
            return "acyl"

        return None

    @classmethod
    @requires_package("rdkit")
    def _remove_duplicate_smiles(cls, smiles: List[str]) -> List[str]:
        from constructure.utilities.rdkit import remove_duplicate_smiles

        return remove_duplicate_smiles(smiles)


class OpenEyeConstructor(Constructor):  # pragma: no cover
    """A class which provides methods for enumerating the possible attachments of
    substituents to different molecular scaffolds using the OpenEye toolkit."""

    @classmethod
    @requires_oe_package("oechem")
    def n_replaceable_groups(cls, scaffold: Scaffold) -> int:

        from openeye import oechem

        scaffold_molecule = oechem.OEMol()
        oechem.OESmilesToMol(scaffold_molecule, scaffold.smiles)

        return sum(
            1 for atom in scaffold_molecule.GetAtoms() if oechem.OEIsRGroup()(atom)
        )

    @classmethod
    @requires_oe_package("oechem")
    def attach_substituents(
        cls, scaffold: Scaffold, substituents: Dict[int, str]
    ) -> str:

        from openeye import oechem

        scaffold_smiles = re.sub(r"\(\[R([1-9])+]\)", r"([\1*])", scaffold.smiles)

        reactant = oechem.OEMol()
        oechem.OESmilesToMol(reactant, scaffold_smiles)

        for i in range(len(substituents)):

            substituent = substituents[i + 1].replace("[R]", "[*:1]")

            rxn = oechem.OEUniMolecularRxn(f"[*:1]-[{i + 1}*]>>{substituent}")
            rxn(reactant)

        return oechem.OEMolToSmiles(oechem.OEMol(reactant))

    @classmethod
    @requires_oe_package("oechem")
    def classify_substituent(
        cls, substituent: str
    ) -> Optional[Literal["hydrogen", "aryl", "alkyl", "acyl", "halogen", "hetero"]]:

        from openeye import oechem

        if substituent.find("[R]") < 0 or substituent.find("[R]") != substituent.rfind(
            "[R]"
        ):

            raise ValueError(
                f"The substituent {substituent} does not contain the expected [R] "
                f"attachment point."
            )

        oe_molecule = oechem.OEMol()
        oechem.OESmilesToMol(oe_molecule, substituent.replace("[R]", "[R1]"))

        # Find the dummy atom.
        dummy_atom = [
            atom for atom in oe_molecule.GetAtoms() if oechem.OEIsRGroup()(atom)
        ][0]

        assert dummy_atom.GetExplicitDegree() == 1

        attachment_atom = next(iter(dummy_atom.GetAtoms()))

        # Check for hydogren
        if attachment_atom.GetAtomicNum() == 1:
            return "hydrogen"

        # Check for halogen
        if attachment_atom.GetAtomicNum() in [9, 17, 35, 53]:
            return "halogen"

        # Check for hetero
        if attachment_atom.GetAtomicNum() in [7, 8, 15, 16]:
            return "hetero"

        # Check for alkyl
        if (
            attachment_atom.GetAtomicNum() == 6
            and (
                attachment_atom.GetExplicitDegree()
                + attachment_atom.GetImplicitHCount()
            )
            == 4
        ):
            return "alkyl"

        # Check for aryl
        if attachment_atom.GetAtomicNum() in [6, 7] and attachment_atom.IsAromatic():
            return "aryl"

        # Check for acyl
        if attachment_atom.GetAtomicNum() == 6 and any(
            neighbor.GetAtomicNum() == 8
            and (neighbor.GetExplicitDegree() + neighbor.GetImplicitHCount()) == 1
            for neighbor in attachment_atom.GetAtoms()
        ):
            return "acyl"

        return None

    @classmethod
    @requires_oe_package("oechem")
    def _remove_duplicate_smiles(cls, smiles: List[str]) -> List[str]:
        from constructure.utilities.openeye import remove_duplicate_smiles

        return remove_duplicate_smiles(smiles)
