"""MutationSet class definition."""

from copy import deepcopy
from multipledispatch import dispatch
from typing import Iterable, Union, List
from .Mutation import Mutation


class MutationSet(list):
    """MutationSet class definition.

    MutationSets are an iterable container around Mutations.
    The Mutations inside a MutationSet are orderable, hence this class inherits from lists.
    """

    def __init__(self, mutations: Iterable[Union[Mutation, str]]):
        """Initialize a MutationSet.

        MutationSets are initializable from an iterable of Mutations or strings.
        """
        # Parse a list of mutation strings:
        parsed_mutations = []
        for i, m in enumerate(mutations):
            if isinstance(m, str):
                parsed_mutations.append(Mutation(m))
            elif isinstance(m, Mutation):
                parsed_mutations.append(m)
            else:
                raise ValueError(
                    f"Mutations must be an iterable of mutation strings or Mutations. Element {m} at position {i} violates this assumption!"
                )
        self.mutations = sorted(parsed_mutations)

    @property
    def positions(self) -> List[int]:
        """Return the position of mutations in the MutationSet.

        :returns: A list of integers.
        """
        positions = []
        for mutation in self:
            positions.append(mutation.position)
        return positions

    def __add__(self, other):
        """Add stuff to mutations."""
        # self.mutations = [m + other for m in self.mutations]
        return _add(self, other)

    def __str__(self):
        return ";".join(str(i) for i in self.mutations)

    def __repr__(self):
        return str(self.mutations)

    def __iter__(self):
        return iter(self.mutations)

    def __next__(self):
        return next(self.mutations)

    def to_str(self):
        """Convert the MutationSet object into a string.

        This is a convenience method to help match the API of SeqLikes,
        which also have a `.to_str()` method.
        """
        return self.__str__()


@dispatch(MutationSet, int)
def _add(obj: MutationSet, other: int):
    """Offset all Mutations' positions in a MutationSet by an integer amount.

    This enables us to write:

    ```python
    ms = MutationSet(...)
    offset_ms = ms + 1
    ```

    instead of:

    ```python
    offset_ms = MutationSet([m + 1 for m in ms])
    """
    mutations = [deepcopy(m) + other for m in obj.mutations]
    return MutationSet(mutations=mutations)


@dispatch(MutationSet, MutationSet)
def _add(obj: MutationSet, other: MutationSet):
    """Add a MutationSet to another MutationSet.

    Add two MutationSets together.
    """
    new_mutations = deepcopy(obj.mutations)
    for mutation in other.mutations:
        new_mutations.append(mutation)

    new_mutations = sorted(new_mutations)

    return MutationSet(mutations=new_mutations)
