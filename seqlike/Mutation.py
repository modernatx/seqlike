"""Mutation and MutationSet class definitions.

Placeholder for docstrings: Please see issue #57 on GitHub. https://github.com/modernatx/seqlike/issues/57

TODO: Flesh out this docstring better.
"""
from typing import Optional, Iterable
from copy import deepcopy


class Mutation:
    def __init__(
        self,
        mutation_string: Optional[str] = None,
        wt_letter: str = "",
        position: Optional[int] = None,
        mutant_letter: Optional[str] = None,
    ):
        """Initialize a Mutation object.

        We initialize a mutation object using the following order of priority:

        1. If a mutation_string in the form of "K15R" (wt-position-mutation)
           or "15-" (position-mutation)
           or "^15R" is provided,
           we use the mutation_string itself to initialize the object.
        2. Otherwise, both a position and a mutant letter must be provided.
        """
        if position is None and mutant_letter is None and mutation_string is None:
            raise ValueError("At least (position, mutant_letter) or (mutation_string) must be provided!")
        if mutation_string:
            wt_letter, position, mutant_letter = parse(mutation_string)

        self.wt_letter = wt_letter
        if position is None:
            raise ValueError("Mutations must have a position specified!")
        self.position = position

        if mutant_letter is None:
            raise ValueError("Mutations must have a mutant letter specified!")
        self.mutant_letter = mutant_letter
        self._one_indexed = False  # internal flag, not to be modified.

    def __repr__(self):
        return self.__str__()

    def __str__(self):
        wt_letter = self.wt_letter
        if self.wt_letter is None:
            wt_letter = ""
        return f"{wt_letter}{self.position}{self.mutant_letter}"

    def __add__(self, other: int):
        position = self.position + other
        obj = deepcopy(self)
        obj.position = position
        return obj

    def __sub__(self, other: int):
        return self.__add__(other * -1)

    def __gt__(self, other: "Mutation"):
        if self.position != other.position:
            return self.position > other.position
        else:
            return self.mutant_letter > other.mutant_letter

    def __lt__(self, other: "Mutation"):
        if self.position != other.position:
            return self.position < other.position
        else:
            return self.mutant_letter < other.mutant_letter

    def __eq__(self, other: "Mutation"):
        return (
            self.position == other.position
            and self.mutant_letter == other.mutant_letter
            and self.wt_letter == other.wt_letter
        )


class Substitution(Mutation):
    pass


class Deletion(Mutation):
    def __init__(self, mutation_string=None, position=None, mutant_letter=None):
        super().__init__(mutation_string=mutation_string, position=position, mutant_letter="-")
        # Simply override mutant_letter regardless of what users pass in.
        self.mutant_letter = "-"


class Insertion(Mutation):
    def __str__(self):
        return f"^{self.position}{self.mutant_letter}"


from seqlike.alphabets import STANDARD_AA_SET, STANDARD_NT_SET

STANDARD_LETTERS = STANDARD_AA_SET.union(STANDARD_NT_SET)


def parse(mutation_string: str):
    """Parse mutation string."""
    # Case 1: Mutation string begins with ^.
    if mutation_string[0] == "^" or mutation_string[0] in STANDARD_LETTERS:
        wt = mutation_string[0]
        pos = mutation_string[1:-1]
        mut = mutation_string[-1]

    else:
        wt = None
        pos = mutation_string[:-1]
        mut = mutation_string[-1]
    try:
        return wt, int(pos), mut
    except ValueError as e:
        raise ValueError(
            "It looks like you have passed in extraneous letters; "
            "please check that your mutation string conforms to the pattern "
            "<wt><pos><mut> (e.g. K15R) or <pos><mut> (e.g. 15R)!"
        ) from e


def magical_parse(mutation_string: str):
    """Magically parse a mutation string to return the correct type of mutation.

    :param mutation_string: The mutation string to parse.
    :returns: One of an Insertion, Deletion, or Substitution.
    """
    if mutation_string[0] == "^":
        return Insertion(mutation_string=mutation_string)
    elif mutation_string[-1] == "-":
        return Deletion(mutation_string=mutation_string)
    else:
        return Substitution(mutation_string=mutation_string)
