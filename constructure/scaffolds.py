from typing import TYPE_CHECKING, Dict, List

from pydantic import BaseModel
from typing_extensions import Literal

RGroup = Literal["hydrogen", "acyl", "alkyl", "aryl", "halogen", "hetero"]

if TYPE_CHECKING:
    RGroup = str


class Scaffold(BaseModel):
    """A model describing a particular molecular scaffold.

    Examples:
        Defines a scaffold for a carbon atom which can have three different
        substituents, each with a different type::

            >>> Scaffold(
            >>>     smiles="C([R1])([R2])([R3])",
            >>>     r_groups={
            >>>         1: ["alkyl"], 2: ["aryl"], 3: ["halogen"]
            >>>     },
            >>> )
    )
    """

    smiles: str
    r_groups: Dict[int, List[RGroup]]


SCAFFOLDS = {
    # cation
    # anion
    # carbonyl compound
    "aldehyde": Scaffold(
        smiles="C([R1])(=O)", r_groups={1: ["hydrogen", "alkyl", "aryl"]}
    ),
    "ketone": Scaffold(
        smiles="C([R1])(=O)([R2])",
        r_groups={1: ["alkyl", "aryl"], 2: ["alkyl", "aryl"]},
    ),
    # thiocarbonyl compound
    "thioaldehyde": Scaffold(
        smiles="C([R1])(=S)", r_groups={1: ["hydrogen", "alkyl", "aryl"]}
    ),
    "thioketone": Scaffold(
        smiles="C([R1])(=S)([R2])",
        r_groups={1: ["hydrogen", "alkyl", "aryl"], 2: ["hydrogen", "alkyl", "aryl"]},
    ),
    "imine": Scaffold(
        smiles="C([R1])(=N([R3]))([R2])",
        r_groups={
            1: ["hydrogen", "alkyl", "aryl"],
            2: ["hydrogen", "alkyl", "aryl"],
            3: ["hydrogen", "alkyl", "aryl"],
        },
    ),
    "hydrazone": Scaffold(
        smiles="C([R1])(=N(N([R3])([R4])))([R2])",
        r_groups={
            1: ["hydrogen", "alkyl", "aryl"],
            2: ["hydrogen", "alkyl", "aryl"],
            3: ["hydrogen", "alkyl", "aryl"],
            4: ["hydrogen", "alkyl", "aryl"],
        },
    ),
    "semicarbazone": Scaffold(
        smiles="C([R1])(=N(N(C(=O)N([R4])([R5]))([R3])))([R2])",
        r_groups={
            1: ["hydrogen", "alkyl", "aryl"],
            2: ["hydrogen", "alkyl", "aryl"],
            3: ["hydrogen", "alkyl", "aryl"],
            4: ["hydrogen", "alkyl", "aryl"],
            5: ["hydrogen", "alkyl", "aryl"],
        },
    ),
    "thiosemicarbazone": Scaffold(
        smiles="C([R1])(=N(N(C(=S)N([R4])([R5]))([R3])))([R2])",
        r_groups={
            1: ["hydrogen", "alkyl", "aryl"],
            2: ["hydrogen", "alkyl", "aryl"],
            3: ["hydrogen", "alkyl", "aryl"],
            4: ["hydrogen", "alkyl", "aryl"],
            5: ["hydrogen", "alkyl", "aryl"],
        },
    ),
    "oxime": Scaffold(
        smiles="C([R1])(=NO)([R2])",
        r_groups={
            1: ["hydrogen", "alkyl", "aryl"],
            2: ["hydrogen", "alkyl", "aryl"],
        },
    ),
    "oxime ether": Scaffold(
        smiles="C([R1])(=NO([R3]))([R2])",
        r_groups={
            1: ["hydrogen", "alkyl", "aryl"],
            2: ["hydrogen", "alkyl", "aryl"],
            3: ["alkyl", "aryl"],
        },
    ),
    "ketene": Scaffold(
        smiles="C([R1])(=C(=O))([R2])",
        r_groups={
            1: ["hydrogen", "alkyl", "aryl"],
            2: ["hydrogen", "alkyl", "aryl"],
        },
    ),
    "ketene acetal derivative": Scaffold(
        smiles="C([R1])([R2])=C([R3])([R4])",
        r_groups={
            1: ["hydrogen", "alkyl", "aryl"],
            2: ["hydrogen", "alkyl", "aryl"],
            3: ["hetero"],
            4: ["hetero"],
        },
    ),
    "carbonyl hydrate": Scaffold(
        smiles="C([R1])(O)(O)([R2])",
        r_groups={
            1: ["hydrogen", "alkyl", "aryl"],
            2: ["hydrogen", "alkyl", "aryl"],
        },
    ),
    "hemiacetal": Scaffold(
        smiles="C([R1])(O([R3]))(O)([R2])",
        r_groups={
            1: ["hydrogen", "alkyl", "aryl"],
            2: ["hydrogen", "alkyl", "aryl"],
            3: ["alkyl", "aryl"],
        },
    ),
    "acetal": Scaffold(
        smiles="C([R1])(O([R3]))(O([R4]))([R2])",
        r_groups={
            1: ["hydrogen", "alkyl", "aryl"],
            2: ["hydrogen", "alkyl", "aryl"],
            3: ["alkyl", "aryl"],
            4: ["alkyl", "aryl"],
        },
    ),
    "hemiaminal": Scaffold(
        smiles="C([R1])(O([R3]))(N([R4])([R5]))([R2])",
        r_groups={
            1: ["hydrogen", "alkyl", "aryl"],
            2: ["hydrogen", "alkyl", "aryl"],
            3: ["hydrogen", "alkyl", "aryl"],
            4: ["hydrogen", "alkyl", "aryl"],
            5: ["hydrogen", "alkyl", "aryl"],
        },
    ),
    "aminal": Scaffold(
        smiles="C([R1])(N([R3])([R4]))(N([R5])([R6]))([R2])",
        r_groups={
            1: ["hydrogen", "alkyl", "aryl"],
            2: ["hydrogen", "alkyl", "aryl"],
            3: ["hydrogen", "alkyl", "aryl"],
            4: ["hydrogen", "alkyl", "aryl"],
            5: ["hydrogen", "alkyl", "aryl"],
            6: ["hydrogen", "alkyl", "aryl"],
        },
    ),
    "thiohemiaminal": Scaffold(
        smiles="C([R1])(S([R3]))(N([R4])([R5]))([R2])",
        r_groups={
            1: ["hydrogen", "alkyl", "aryl"],
            2: ["hydrogen", "alkyl", "aryl"],
            3: ["hydrogen", "alkyl", "aryl"],
            4: ["hydrogen", "alkyl", "aryl"],
            5: ["hydrogen", "alkyl", "aryl"],
        },
    ),
    "thioacetal": Scaffold(
        smiles="C([R1])(S([R3]))(S([R4]))([R2])",
        r_groups={
            1: ["hydrogen", "alkyl", "aryl"],
            2: ["hydrogen", "alkyl", "aryl"],
            3: ["alkyl", "aryl"],
            4: ["alkyl", "aryl"],
        },
    ),
    "enamine": Scaffold(
        smiles="C([R1])(=C([R3])(N([R4])([R5])))([R2])",
        r_groups={
            1: ["hydrogen", "acyl", "alkyl", "aryl"],
            2: ["hydrogen", "acyl", "alkyl", "aryl"],
            3: ["hydrogen", "acyl", "alkyl", "aryl"],
            4: ["hydrogen", "acyl", "alkyl", "aryl"],
            5: ["hydrogen", "acyl", "alkyl", "aryl"],
        },
    ),
    "enol": Scaffold(
        smiles="C([R1])(=C([R3])O)([R2])",
        r_groups={
            1: ["hydrogen", "acyl", "alkyl", "aryl"],
            2: ["hydrogen", "acyl", "alkyl", "aryl"],
            3: ["hydrogen", "acyl", "alkyl", "aryl"],
        },
    ),
    "enolether": Scaffold(
        smiles="C([R1])([R2])=C([R3])O([R4])",
        r_groups={
            1: ["hydrogen", "acyl", "alkyl", "aryl"],
            2: ["hydrogen", "acyl", "alkyl", "aryl"],
            3: ["hydrogen", "acyl", "alkyl", "aryl"],
            4: ["alkyl", "aryl"],
        },
    ),
    # hydroxy compound
    "alcohol": Scaffold(
        smiles="O([R1])",
        r_groups={
            1: ["alkyl", "aryl"],
        },
    ),
    "primary alcohol": Scaffold(
        smiles="C([R1])(O)",
        r_groups={
            1: ["alkyl", "aryl"],
        },
    ),
    "secondary alcohol": Scaffold(
        smiles="C([R1])(O)([R2])",
        r_groups={
            1: ["alkyl", "aryl"],
            2: ["alkyl", "aryl"],
        },
    ),
    "tertiary alcohol": Scaffold(
        smiles="C([R1])(O)([R2])([R3])",
        r_groups={
            1: ["alkyl", "aryl"],
            2: ["alkyl", "aryl"],
            3: ["alkyl", "aryl"],
        },
    ),
    "1,2-diol": Scaffold(
        smiles="C([R1])([R2])(O)(C([R3])([R4])(O))",
        r_groups={
            1: ["hydrogen", "alkyl", "aryl"],
            2: ["hydrogen", "alkyl", "aryl"],
            3: ["hydrogen", "alkyl", "aryl"],
            4: ["hydrogen", "alkyl", "aryl"],
        },
    ),
    "1,2-aminoalcohol": Scaffold(
        smiles="C([R1])([R2])(O)(C([R3])([R4])(N))",
        r_groups={
            1: ["hydrogen", "alkyl", "aryl"],
            2: ["hydrogen", "alkyl", "aryl"],
            3: ["hydrogen", "alkyl", "aryl"],
            4: ["hydrogen", "alkyl", "aryl"],
        },
    ),
    # phenol,
    # 1,2-diphenol
    "enediol": Scaffold(
        smiles="C([R1])(O)(=C([R2])(O))",
        r_groups={
            1: ["hydrogen", "alkyl", "aryl"],
            2: ["hydrogen", "alkyl", "aryl"],
        },
    ),
    "ether": Scaffold(
        smiles="O([R1])([R2])",
        r_groups={
            1: ["alkyl", "aryl"],
            2: ["alkyl", "aryl"],
        },
    ),
    # dialkylether
    # alkylarylether
    # diarylether
    "thioether": Scaffold(
        smiles="S([R1])([R2])",
        r_groups={
            1: ["alkyl", "aryl"],
            2: ["alkyl", "aryl"],
        },
    ),
    "disulfide": Scaffold(
        smiles="S([R1])(S([R2]))",
        r_groups={
            1: ["alkyl", "aryl"],
            2: ["alkyl", "aryl"],
        },
    ),
    "peroxide": Scaffold(
        smiles="O([R1])(O([R2]))",
        r_groups={
            1: ["hydrogen", "alkyl", "aryl"],
            2: ["alkyl", "aryl"],
        },
    ),
    "hydroperoxide": Scaffold(
        smiles="O([R1])O",
        r_groups={
            1: ["alkyl", "aryl"],
        },
    ),
    "hydrazine derivative": Scaffold(
        smiles="N([R1])([R2])(N([R3])([R4]))",
        r_groups={
            1: ["hydrogen", "acyl", "alkyl", "aryl"],
            2: ["hydrogen", "acyl", "alkyl", "aryl"],
            3: ["hydrogen", "acyl", "alkyl", "aryl"],
            4: ["hydrogen", "acyl", "alkyl", "aryl"],
        },
    ),
    "hydroxylamine": Scaffold(
        smiles="N([R1])([R2])(O([R3]))",
        r_groups={
            1: ["hydrogen", "alkyl", "aryl"],
            2: ["hydrogen", "alkyl", "aryl"],
            3: ["hydrogen", "alkyl", "aryl"],
        },
    ),
    "amine": Scaffold(
        smiles="N([R1])([R2])([R3])",
        r_groups={
            1: ["hydrogen", "alkyl", "aryl"],
            2: ["hydrogen", "alkyl", "aryl"],
            3: ["hydrogen", "alkyl", "aryl"],
        },
    ),
    # primary,
    # secondary,
    # tertiary
    "quaternary ammonium salt": Scaffold(
        smiles="[N+]([R1])([R2])([R3])([R4])",
        r_groups={
            1: ["hydrogen", "alkyl", "aryl"],
            2: ["hydrogen", "alkyl", "aryl"],
            3: ["hydrogen", "alkyl", "aryl"],
            4: ["hydrogen", "alkyl", "aryl"],
        },
    ),
    "N-oxide": Scaffold(
        smiles="[N+]([R1])([R2])([R3])([O-])",
        r_groups={
            1: ["alkyl", "aryl"],
            2: ["alkyl", "aryl"],
            3: ["alkyl", "aryl"],
        },
    ),
    # halogen deriv.
    # alkyl halide
    "alkyl flouride": Scaffold(
        smiles="F([R1])",
        r_groups={
            1: ["alkyl"],
        },
    ),
    "alkyl chloride": Scaffold(
        smiles="Cl([R1])",
        r_groups={
            1: ["alkyl"],
        },
    ),
    "alkyl bromide": Scaffold(
        smiles="Br([R1])",
        r_groups={
            1: ["alkyl"],
        },
    ),
    "alkyl iodide": Scaffold(
        smiles="I([R1])",
        r_groups={
            1: ["alkyl"],
        },
    ),
    # aryl halide
    "aryl flouride": Scaffold(
        smiles="F([R1])",
        r_groups={
            1: ["aryl"],
        },
    ),
    "aryl chloride": Scaffold(
        smiles="Cl([R1])",
        r_groups={
            1: ["aryl"],
        },
    ),
    "aryl bromide": Scaffold(
        smiles="Br([R1])",
        r_groups={
            1: ["aryl"],
        },
    ),
    "aryl iodide": Scaffold(
        smiles="I([R1])",
        r_groups={
            1: ["aryl"],
        },
    ),
    # organometallic
    # carboxylic acid deriv
    "carboxylic acid": Scaffold(
        smiles="C([R1])(=O)O",
        r_groups={
            1: ["hydrogen", "alkyl", "aryl"],
        },
    ),
    "carboxylic acid salt": Scaffold(
        smiles="C([R1])(=O)[O-]",
        r_groups={
            1: ["hydrogen", "alkyl", "aryl"],
        },
    ),
    "carboxylic acid ester": Scaffold(
        smiles="C([R1])(=O)O([R2])",
        r_groups={
            1: ["hydrogen", "alkyl", "aryl"],
            2: ["alkyl", "aryl"],
        },
    ),
    # lactone
    "carboxylic acid amide": Scaffold(
        smiles="C([R1])(=O)N([R2])([R3])",
        r_groups={
            1: ["hydrogen", "alkyl", "aryl"],
            2: ["hydrogen", "alkyl", "aryl"],
            3: ["hydrogen", "alkyl", "aryl"],
        },
    ),
    # primary,
    # secondary,
    # tertiary.
    # lactam
    "carboxylic acid hydrazide": Scaffold(
        smiles="C([R1])(=O)N([R2])(N([R3])([R4]))",
        r_groups={
            1: ["hydrogen", "alkyl", "aryl"],
            2: ["hydrogen", "alkyl", "aryl"],
            3: ["hydrogen", "alkyl", "aryl"],
            4: ["hydrogen", "alkyl", "aryl"],
        },
    ),
    "carboxylic acid azide": Scaffold(
        smiles="C([R1])(=O)N=[N+]=[N-]",
        r_groups={
            1: ["hydrogen", "alkyl", "aryl"],
        },
    ),
    "hydroxamic acid": Scaffold(
        smiles="C([R1])(=O)NO",
        r_groups={
            1: ["hydrogen", "alkyl", "aryl"],
        },
    ),
    "carboxylic acid amidine": Scaffold(
        smiles="C([R1])(=N([R4]))N([R2])([R3])",
        r_groups={
            1: ["hydrogen", "alkyl", "aryl"],
            2: ["hydrogen", "alkyl", "aryl"],
            3: ["hydrogen", "alkyl", "aryl"],
            4: ["hydrogen", "alkyl", "aryl"],
        },
    ),
    "carboxylic acid amidrazone": Scaffold(
        smiles="C([R1])(=N([R5]))(N([R2])(N([R3])([R4])))",
        r_groups={
            1: ["hydrogen", "alkyl", "aryl"],
            2: ["hydrogen", "alkyl", "aryl"],
            3: ["hydrogen", "alkyl", "aryl"],
            4: ["hydrogen", "alkyl", "aryl"],
            5: ["hydrogen", "alkyl", "aryl"],
        },
    ),
    "nitrile": Scaffold(
        smiles="C([R1])#N",
        r_groups={
            1: ["hydrogen", "alkyl", "aryl"],
        },
    ),
    # acyl halide
    "acyl flouride": Scaffold(
        smiles="FC(=O)([R1])",
        r_groups={
            1: ["hydrogen", "alkyl", "aryl"],
        },
    ),
    "acyl chloride": Scaffold(
        smiles="ClC(=O)([R1])",
        r_groups={
            1: ["hydrogen", "alkyl", "aryl"],
        },
    ),
    "acyl bromide": Scaffold(
        smiles="BrC(=O)([R1])",
        r_groups={
            1: ["hydrogen", "alkyl", "aryl"],
        },
    ),
    "acyl iodide": Scaffold(
        smiles="IC(=O)([R1])",
        r_groups={
            1: ["hydrogen", "alkyl", "aryl"],
        },
    ),
    "acyl cyanide": Scaffold(
        smiles="C([R1])(=O)C#N",
        r_groups={
            1: ["hydrogen", "alkyl", "aryl"],
        },
    ),
    "imido ester": Scaffold(
        smiles="C([R1])(=N([R3]))O([R2])",
        r_groups={
            1: ["hydrogen", "alkyl", "aryl"],
            2: ["alkyl", "aryl"],
            3: ["hydrogen", "alkyl", "aryl"],
        },
    ),
    "imidoyl halide": Scaffold(
        smiles="C([R1])(=N([R2]))([R3])",
        r_groups={
            1: ["hydrogen", "alkyl", "aryl"],
            2: ["hydrogen", "alkyl", "aryl"],
            3: ["halogen"],
        },
    ),
    # thiocarboxylic acid deriv
    "thiocarboxylic acid": Scaffold(
        smiles="C([R1])(=S)O",
        r_groups={
            1: ["hydrogen", "alkyl", "aryl"],
        },
    ),
    # SH variant
    "thiocarboxylic acid ester": Scaffold(
        smiles="C([R1])(=S)O([R2])",
        r_groups={
            1: ["hydrogen", "alkyl", "aryl"],
            2: ["alkyl", "aryl"],
        },
    ),
    # SH variant
    # thiolactone
    "thiocarboxylic acid amide": Scaffold(
        smiles="C([R1])(=S)N([R2])([R3])",
        r_groups={
            1: ["hydrogen", "alkyl", "aryl"],
            2: ["hydrogen", "alkyl", "aryl"],
            3: ["hydrogen", "alkyl", "aryl"],
        },
    ),
    # thiolactam
    "imidothioester": Scaffold(
        smiles="C([R1])(=N([R3]))S([R2])",
        r_groups={
            1: ["hydrogen", "alkyl", "aryl"],
            2: ["alkyl", "aryl"],
            3: ["hydrogen", "alkyl", "aryl"],
        },
    ),
    "oxohetarene": Scaffold(
        smiles="n1([R1])ccccc1=O",
        r_groups={
            1: ["hydrogen", "alkyl", "aryl"],
        },
    ),
    "thioxohetarene": Scaffold(
        smiles="n1([R1])ccccc1=S",
        r_groups={
            1: ["hydrogen", "alkyl", "aryl"],
        },
    ),
    "iminohetarene": Scaffold(
        smiles="N([R2])=C1C=CC=CN1([R1])",
        r_groups={
            1: ["hydrogen", "alkyl", "aryl"],
            2: ["hydrogen", "alkyl", "aryl"],
        },
    ),
    # orthocarboxylic acid deriv.
    "carboxylic acid orthoester": Scaffold(
        smiles="C([R1])(O([R2]))(O([R3]))(O([R4]))",
        r_groups={
            1: ["hydrogen", "alkyl", "aryl"],
            2: ["alkyl", "aryl"],
            3: ["alkyl", "aryl"],
            4: ["alkyl", "aryl"],
        },
    ),
    "carboxylic acid amide acetal": Scaffold(
        smiles="C([R1])(O([R2]))(O([R3]))(N([R4])([R5]))",
        r_groups={
            1: ["hydrogen", "alkyl", "aryl"],
            2: ["alkyl", "aryl"],
            3: ["alkyl", "aryl"],
            4: ["hydrogen", "alkyl", "aryl"],
            5: ["hydrogen", "alkyl", "aryl"],
        },
    ),
    "carboxylic acid anhydride": Scaffold(
        smiles="C([R1])(=O)OC(=O)([R2])",
        r_groups={
            1: ["hydrogen", "alkyl", "aryl"],
            2: ["hydrogen", "alkyl", "aryl"],
        },
    ),
    "carboxylic acid imide": Scaffold(
        smiles="C([R1])(=O)N([R3])(C(=O)([R2]))",
        r_groups={
            1: ["hydrogen", "alkyl", "aryl"],
            2: ["hydrogen", "alkyl", "aryl"],
            3: ["hydrogen", "alkyl", "aryl"],
        },
    ),
    # carboxylic acid unsubst imide
    # carboxylic acid subst imide
    # CO2 deriv
    # carbonic acid deriv
    "carbonic acid monoester": Scaffold(
        smiles="OC(=O)O([R1])",
        r_groups={
            1: ["alkyl", "aryl"],
        },
    ),
    "carbonic acid diester": Scaffold(
        smiles="C(=O)(O([R1]))(O([R2]))",
        r_groups={
            1: ["alkyl", "aryl"],
            2: ["alkyl", "aryl"],
        },
    ),
    "carbonic acid ester halide": Scaffold(
        smiles="C(=O)(O([R1]))([R2])",
        r_groups={
            1: ["alkyl", "aryl"],
            2: ["halogen"],
        },
    ),
    # thiocarbonic acid deriv
    "thiocarbonic acid monoester": Scaffold(
        smiles="OC(=S)O([R1])",
        r_groups={
            1: ["alkyl", "aryl"],
        },
    ),
    "thiocarbonic acid diester": Scaffold(
        smiles="C(=S)(O([R1]))(O([R2]))",
        r_groups={
            1: ["alkyl", "aryl"],
            2: ["alkyl", "aryl"],
        },
    ),
    "thiocarbonic acid ester halide": Scaffold(
        smiles="C(=S)(O([R1]))([R2])",
        r_groups={
            1: ["alkyl", "aryl"],
            2: ["halogen"],
        },
    ),
    # carbamic acid deriv
    "carbamic acid": Scaffold(
        smiles="OC(=O)N([R1])([R2])",
        r_groups={
            1: ["hydrogen", "alkyl", "aryl"],
            2: ["hydrogen", "alkyl", "aryl"],
        },
    ),
    "carbamic acid ester": Scaffold(
        smiles="C(=O)(N([R1])([R2]))(O([R3]))",
        r_groups={
            1: ["hydrogen", "alkyl", "aryl"],
            2: ["hydrogen", "alkyl", "aryl"],
            3: ["alkyl", "aryl"],
        },
    ),
    "carbamic acid halide": Scaffold(
        smiles="C([R3])(=O)(N([R1])([R2]))",
        r_groups={
            1: ["hydrogen", "alkyl", "aryl"],
            2: ["hydrogen", "alkyl", "aryl"],
            3: ["halogen"],
        },
    ),
    # thiocarbamic acid deriv
    "thiocarbamic acid": Scaffold(
        smiles="OC(=S)N([R1])([R2])",
        r_groups={
            1: ["hydrogen", "alkyl", "aryl"],
            2: ["hydrogen", "alkyl", "aryl"],
        },
    ),
    "thiocarbamic acid ester": Scaffold(
        smiles="C(=S)(N([R1])([R2]))(O([R3]))",
        r_groups={
            1: ["hydrogen", "alkyl", "aryl"],
            2: ["hydrogen", "alkyl", "aryl"],
            3: ["alkyl", "aryl"],
        },
    ),
    "thiocarbamic acid halide": Scaffold(
        smiles="C([R3])(=S)(N([R1])([R2]))",
        r_groups={
            1: ["hydrogen", "alkyl", "aryl"],
            2: ["hydrogen", "alkyl", "aryl"],
            3: ["halogen"],
        },
    ),
    "urea": Scaffold(
        smiles="C(=O)(N([R1])([R2]))(N([R3])([R4]))",
        r_groups={
            1: ["hydrogen", "alkyl", "aryl"],
            2: ["hydrogen", "alkyl", "aryl"],
            3: ["hydrogen", "alkyl", "aryl"],
            4: ["hydrogen", "alkyl", "aryl"],
        },
    ),
    "isourea": Scaffold(
        smiles="C(O([R4]))(N([R1])([R2]))(=N([R3]))",
        r_groups={
            1: ["hydrogen", "alkyl", "aryl"],
            2: ["hydrogen", "alkyl", "aryl"],
            3: ["hydrogen", "alkyl", "aryl"],
            4: ["hydrogen", "alkyl", "aryl"],
        },
    ),
    "thiourea": Scaffold(
        smiles="C(=S)(N([R1])([R2]))(N([R3])([R4]))",
        r_groups={
            1: ["hydrogen", "alkyl", "aryl"],
            2: ["hydrogen", "alkyl", "aryl"],
            3: ["hydrogen", "alkyl", "aryl"],
            4: ["hydrogen", "alkyl", "aryl"],
        },
    ),
    "isothiourea": Scaffold(
        smiles="C(S([R4]))(N([R1])([R2]))(=N([R3]))",
        r_groups={
            1: ["hydrogen", "alkyl", "aryl"],
            2: ["hydrogen", "alkyl", "aryl"],
            3: ["hydrogen", "alkyl", "aryl"],
            4: ["hydrogen", "alkyl", "aryl"],
        },
    ),
    "guanidine": Scaffold(
        smiles="C(=N([R5]))(N([R1])([R2]))(N([R3])([R4]))",
        r_groups={
            1: ["hydrogen", "alkyl", "aryl"],
            2: ["hydrogen", "alkyl", "aryl"],
            3: ["hydrogen", "alkyl", "aryl"],
            4: ["hydrogen", "alkyl", "aryl"],
            5: ["hydrogen", "alkyl", "aryl"],
        },
    ),
    "semicarbazide": Scaffold(
        smiles="C(=O)(N([R1])([R2]))(N([R3])N([R4])([R5]))",
        r_groups={
            1: ["hydrogen", "alkyl", "aryl"],
            2: ["hydrogen", "alkyl", "aryl"],
            3: ["hydrogen", "alkyl", "aryl"],
            4: ["hydrogen", "alkyl", "aryl"],
            5: ["hydrogen", "alkyl", "aryl"],
        },
    ),
    "thiosemicarbazide": Scaffold(
        smiles="C(=S)(N([R1])([R2]))(N([R3])N([R4])([R5]))",
        r_groups={
            1: ["hydrogen", "alkyl", "aryl"],
            2: ["hydrogen", "alkyl", "aryl"],
            3: ["hydrogen", "alkyl", "aryl"],
            4: ["hydrogen", "alkyl", "aryl"],
            5: ["hydrogen", "alkyl", "aryl"],
        },
    ),
    "azide": Scaffold(
        smiles="[N-]=[N+]=N([R1])",
        r_groups={
            1: ["alkyl", "aryl"],
        },
    ),
    "azo compound": Scaffold(
        smiles="N([R1])=N([R2])",
        r_groups={
            1: ["alkyl", "aryl"],
            2: ["alkyl", "aryl"],
        },
    ),
    "diazonium salt": Scaffold(
        smiles="N([R1])=[N+]",
        r_groups={
            1: ["alkyl", "aryl"],
        },
    ),
    "isonitrile": Scaffold(
        smiles="[C-]#[N+]([R1])",
        r_groups={
            1: ["alkyl", "aryl"],
        },
    ),
    "cyanate": Scaffold(
        smiles="C(#N)O([R1])",
        r_groups={
            1: ["alkyl", "aryl"],
        },
    ),
    "isocyanate": Scaffold(
        smiles="O=C=N([R1])",
        r_groups={
            1: ["alkyl", "aryl"],
        },
    ),
    "thiocyanate": Scaffold(
        smiles="C(#N)S([R1])",
        r_groups={
            1: ["alkyl", "aryl"],
        },
    ),
    "isothiocyanate": Scaffold(
        smiles="S=C=N([R1])",
        r_groups={
            1: ["alkyl", "aryl"],
        },
    ),
    "carbodiimide": Scaffold(
        smiles="N([R1])=C=N([R2])",
        r_groups={
            1: ["alkyl", "aryl"],
            2: ["alkyl", "aryl"],
        },
    ),
    "nitroso compound": Scaffold(
        smiles="O=N([R1])",
        r_groups={
            1: ["alkyl", "aryl"],
        },
    ),
    "nitro compound": Scaffold(
        smiles="N([R1])(=O)(=O)",
        r_groups={
            1: ["alkyl", "aryl"],
        },
    ),
    "nitrite": Scaffold(
        smiles="O=NO([R1])",
        r_groups={
            1: ["alkyl", "aryl"],
        },
    ),
    "nitrate": Scaffold(
        smiles="N(=O)(=O)O([R1])",
        r_groups={
            1: ["alkyl", "aryl"],
        },
    ),
    # sulfuric acid deriv.
    # sulfuric acid
    "sulfuric acid monoester": Scaffold(
        smiles="S(=O)(=O)(O)O([R1])",
        r_groups={
            1: ["alkyl", "aryl"],
        },
    ),
    "sulfuric acid diester": Scaffold(
        smiles="S(=O)(=O)(O([R2]))O([R1])",
        r_groups={
            1: ["alkyl", "aryl"],
            2: ["alkyl", "aryl"],
        },
    ),
    "sulfuric acid amide ester": Scaffold(
        smiles="S(=O)(=O)(N([R2])([R3]))O([R1])",
        r_groups={
            1: ["alkyl", "aryl"],
            2: ["hydrogen", "alkyl", "aryl"],
            3: ["hydrogen", "alkyl", "aryl"],
        },
    ),
    "sulfuric acid amide": Scaffold(
        smiles="S(=O)(=O)(N([R1])([R2]))O",
        r_groups={
            1: ["hydrogen", "alkyl", "aryl"],
            2: ["hydrogen", "alkyl", "aryl"],
        },
    ),
    "sulfuric acid diamide": Scaffold(
        smiles="S(=O)(=O)(N([R1])([R2]))N([R3])([R4])",
        r_groups={
            1: ["hydrogen", "alkyl", "aryl"],
            2: ["hydrogen", "alkyl", "aryl"],
            3: ["hydrogen", "alkyl", "aryl"],
            4: ["hydrogen", "alkyl", "aryl"],
        },
    ),
    # sulfuryl halide
    # sulfonic acid deriv
    # sulfonic acid
    "sulfonic acid ester": Scaffold(
        smiles="S(=O)(=O)(O)([R1])",
        r_groups={
            1: ["alkyl", "aryl"],
        },
    ),
    "sulfonamide": Scaffold(
        smiles="S(=O)(=O)(N([R2])([R3]))([R1])",
        r_groups={
            1: ["alkyl", "aryl"],
            2: ["hydrogen", "alkyl", "aryl"],
            3: ["hydrogen", "alkyl", "aryl"],
        },
    ),
    "sulfonyl halide": Scaffold(
        smiles="S([R1])(=O)(=O)([R2])",
        r_groups={
            1: ["alkyl", "aryl"],
            2: ["halogen"],
        },
    ),
    "sulfone": Scaffold(
        smiles="S([R1])(=O)(=O)([R2])",
        r_groups={
            1: ["alkyl", "aryl"],
            2: ["alkyl", "aryl"],
        },
    ),
    "sulfoxide": Scaffold(
        smiles="S([R1])(=O)([R2])",
        r_groups={
            1: ["alkyl", "aryl"],
            2: ["alkyl", "aryl"],
        },
    ),
    # sulfinic acid deriv.
    "sulfinic acid": Scaffold(
        smiles="S([R1])(=O)O",
        r_groups={
            1: ["alkyl", "aryl"],
        },
    ),
    "sulfinic acid ester": Scaffold(
        smiles="S([R1])(=O)O([R2])",
        r_groups={
            1: ["alkyl", "aryl"],
            2: ["alkyl", "aryl"],
        },
    ),
    "sulfinic acid halide": Scaffold(
        smiles="S([R1])(=O)([R2])",
        r_groups={
            1: ["alkyl", "aryl"],
            2: ["halogen"],
        },
    ),
    "sulfinic acid amide": Scaffold(
        smiles="S([R1])(=O)N([R2])([R3])",
        r_groups={
            1: ["alkyl", "aryl"],
            2: ["hydrogen", "alkyl", "aryl"],
            3: ["hydrogen", "alkyl", "aryl"],
        },
    ),
    # sulfenic acid deriv.
    "sulfenic acid": Scaffold(
        smiles="S([R1])O",
        r_groups={
            1: ["alkyl", "aryl"],
        },
    ),
    "sulfenic acid ester": Scaffold(
        smiles="S([R1])O([R2])",
        r_groups={
            1: ["alkyl", "aryl"],
            2: ["alkyl", "aryl"],
        },
    ),
    "sulfenic acid halide": Scaffold(
        smiles="S([R1])([R2])",
        r_groups={
            1: ["alkyl", "aryl"],
            2: ["halogen"],
        },
    ),
    "sulfenic acid amide": Scaffold(
        smiles="S([R1])N([R2])([R3])",
        r_groups={
            1: ["alkyl", "aryl"],
            2: ["hydrogen", "alkyl", "aryl"],
            3: ["hydrogen", "alkyl", "aryl"],
        },
    ),
    "thiol": Scaffold(
        smiles="S([R1])",
        r_groups={
            1: ["alkyl", "aryl"],
        },
    ),
    "alkylthiol": Scaffold(
        smiles="S([R1])",
        r_groups={
            1: ["alkyl"],
        },
    ),
    "arylthiol": Scaffold(
        smiles="S([R1])",
        r_groups={
            1: ["aryl"],
        },
    ),
    # phosphoric acid deriv.
    # phosphoric acid.
    "phosphoric acid ester": Scaffold(
        smiles="P([R2])([R3])(=O)O([R1])",
        r_groups={
            1: ["alkyl", "aryl"],
            2: ["hetero", "halogen"],
            3: ["hetero", "halogen"],
        },
    ),
    "phosphoric acid halide": Scaffold(
        smiles="P([R1])([R2])([R3])(=O)",
        r_groups={
            1: ["halogen"],
            2: ["hetero", "halogen"],
            3: ["hetero", "halogen"],
        },
    ),
    "phosphoric acid amide": Scaffold(
        smiles="P([R3])([R4])(=O)N([R1])([R2])",
        r_groups={
            1: ["hydrogen", "alkyl", "aryl"],
            2: ["hydrogen", "alkyl", "aryl"],
            3: ["hetero", "halogen"],
            4: ["hetero", "halogen"],
        },
    ),
    # thiophosphoric acid deriv.
    # thiophosphoric acid
    "thiophosphoric acid ester": Scaffold(
        smiles="P([R2])([R3])(=S)O([R1])",
        r_groups={
            1: ["alkyl", "aryl"],
            2: ["hetero", "halogen"],
            3: ["hetero", "halogen"],
        },
    ),
    "thiophosphoric acid halide": Scaffold(
        smiles="P([R1])([R2])([R3])(=S)",
        r_groups={
            1: ["halogen"],
            2: ["hetero", "halogen"],
            3: ["hetero", "halogen"],
        },
    ),
    "thiophosphoric acid amide": Scaffold(
        smiles="P([R3])([R4])(=S)N([R1])([R2])",
        r_groups={
            1: ["hydrogen", "alkyl", "aryl"],
            2: ["hydrogen", "alkyl", "aryl"],
            3: ["hetero", "halogen"],
            4: ["hetero", "halogen"],
        },
    ),
    # phosphonic acid deriv.
    "phosphonic acid": Scaffold(
        smiles="P([R1])(=O)(O)(O)",
        r_groups={
            1: ["alkyl", "aryl"],
        },
    ),
    # phosphonic acid ester
    "phosphine": Scaffold(
        smiles="P([R1])([R2])([R3])",
        r_groups={
            1: ["alkyl", "aryl"],
            2: ["alkyl", "aryl"],
            3: ["alkyl", "aryl"],
        },
    ),
    "phosphinoxide": Scaffold(
        smiles="P(=O)([R1])([R2])([R3])",
        r_groups={
            1: ["alkyl", "aryl"],
            2: ["alkyl", "aryl"],
            3: ["alkyl", "aryl"],
        },
    ),
    # boronic acid deriv.
    # boronic acid
    "boronic acid": Scaffold(
        smiles="B([R1])(O)(O)",
        r_groups={
            1: ["alkyl", "aryl"],
        },
    ),
    "boronic acid ester": Scaffold(
        smiles="B([R1])([R3])O([R2])",
        r_groups={
            1: ["alkyl", "aryl"],
            2: ["alkyl", "aryl"],
            3: ["halogen", "hetero"],
        },
    ),
    "alkene": Scaffold(
        smiles="C([R1])([R2])=C([R3])([R4])",
        r_groups={
            1: ["hydrogen", "alkyl", "aryl"],
            2: ["hydrogen", "alkyl", "aryl"],
            3: ["hydrogen", "alkyl", "aryl"],
            4: ["hydrogen", "alkyl", "aryl"],
        },
    ),
    "alkyne": Scaffold(
        smiles="C([R1])#C([R2])",
        r_groups={
            1: ["hydrogen", "alkyl", "aryl"],
            2: ["hydrogen", "alkyl", "aryl"],
        },
    ),
    # aromatic compound
    # heterocyclic compound
    "alpha-aminoacid": Scaffold(
        smiles="C([R1])(N([R2]))(C(=O)O)",
        r_groups={
            1: ["hydrogen", "alkyl", "aryl"],
            2: ["hydrogen", "alkyl", "aryl"],
        },
    ),
    "alpha-hydroxyacid": Scaffold(
        smiles="C([R1])(O)C(=O)O",
        r_groups={
            1: ["hydrogen", "alkyl", "aryl"],
        },
    ),
    # Non-checkmol
    "phenyl": Scaffold(
        smiles="c1([R1])ccccc1",
        r_groups={
            1: ["hydrogen", "alkyl", "aryl"],
        },
    ),
    "benzyl": Scaffold(
        smiles="C([R1])c1ccccc1",
        r_groups={
            1: ["hydrogen", "alkyl", "aryl"],
        },
    ),
    "o-phenylene": Scaffold(
        smiles="c1([R1])c([R2])cccc1",
        r_groups={
            1: ["hydrogen", "alkyl", "aryl", "acyl", "hetero", "halogen"],
            2: ["hydrogen", "alkyl", "aryl", "acyl", "hetero", "halogen"],
        },
    ),
    "m-phenylene": Scaffold(
        smiles="c1([R1])cc([R2])ccc1",
        r_groups={
            1: ["hydrogen", "alkyl", "aryl", "acyl", "hetero", "halogen"],
            2: ["hydrogen", "alkyl", "aryl", "acyl", "hetero", "halogen"],
        },
    ),
    "p-phenylene": Scaffold(
        smiles="c1([R1])ccc([R2])cc1",
        r_groups={
            1: ["hydrogen", "alkyl", "aryl", "acyl", "hetero", "halogen"],
            2: ["hydrogen", "alkyl", "aryl", "acyl", "hetero", "halogen"],
        },
    ),
    "o-pyridine": Scaffold(
        smiles="n1c([R1])cccc1",
        r_groups={
            1: ["hydrogen", "alkyl", "aryl", "acyl", "hetero", "halogen"],
        },
    ),
    "m-pyridine": Scaffold(
        smiles="n1cc([R1])ccc1",
        r_groups={
            1: ["hydrogen", "alkyl", "aryl", "acyl", "hetero", "halogen"],
        },
    ),
    "p-pyridine": Scaffold(
        smiles="n1ccc([R1])cc1",
        r_groups={
            1: ["hydrogen", "alkyl", "aryl", "acyl", "hetero", "halogen"],
        },
    ),
}
