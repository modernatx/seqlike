"""Mutation and MutationSet class definitions.

Placeholder for docstrings: Please see issue #57 on GitHub. https://github.com/modernatx/seqlike/issues/57

TODO: Flesh out this docstring better.
"""
from copy import deepcopy
from typing import Optional

from seqlike.alphabets import STANDARD_AA_SET, STANDARD_NT_SET


class Mutation:
    def __new__(
        cls,
        mutation_string: str,
    ):
        """Magical constructor that dispatches to the individual classes."""

        if mutation_string is not None and mutation_string[0] == "^":
            return super().__new__(Insertion)
        elif mutation_string is not None and mutation_string[-1] == "-":
            return super().__new__(Deletion)
        else:
            return super().__new__(Substitution)

    def __init__(
        self,
        mutation_string: Optional[str] = None,
    ):
        """Initialize a Mutation object.

        If a mutation_string in the form of "K15R" (wt-position-mutation)
        or "15-" (position-mutation)
        or "^15R" is provided,
        we use the mutation_string itself to initialize the object.
        """
        if mutation_string is None:
            raise ValueError("A mutation string must be provided!")
        wt_letter, position, mutant_letter = parse_mutation(mutation_string)

        self.wt_letter = wt_letter
        self.position = position
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
        """Add an integer to a mutation to offset the position."""
        position = self.position + other
        obj = deepcopy(self)
        obj.position = position
        return obj

    def __sub__(self, other: int):
        """Subtract an integer from a mutation to offset the position."""
        return self.__add__(other * -1)

    def __gt__(self, other: "Mutation"):
        """Greater than comparison between two mutations."""
        if self.position != other.position:
            return self.position > other.position
        else:
            return self.mutant_letter > other.mutant_letter

    def __lt__(self, other: "Mutation"):
        """Less than comparison between two mutations."""
        if self.position != other.position:
            return self.position < other.position
        else:
            return self.mutant_letter < other.mutant_letter

    def __eq__(self, other: "Mutation"):
        """Equality comparison between two mutations."""
        return (
            self.position == other.position
            and self.mutant_letter == other.mutant_letter
            and self.wt_letter == other.wt_letter
        )

    def __deepcopy__(self, memo):
        """Deepcopy implementation.

        Necessary b/c of our use of __new__.

        <!-- #noqa: DAR101 -->
        <!-- #noqa: DAR201 -->
        """
        mutation_string = str(self)

        return Mutation(mutation_string)


class Substitution(Mutation):
    """Substitution class.

    This is an empty class, exists mostly for providing a Type for dispatch.
    """

    pass


class Deletion(Mutation):
    """Deletion class.

    This is an empty class, exists mostly for providing a Type for dispatch.
    """

    pass


class Insertion(Mutation):
    def __str__(self):
        return f"^{self.position}{self.mutant_letter}"


STANDARD_LETTERS = STANDARD_AA_SET.union(STANDARD_NT_SET)


def parse_mutation(mutation_string: str):
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
