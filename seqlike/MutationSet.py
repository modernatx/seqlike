"""MutationSet class definition."""

from multipledispatch import dispatch
from typing import Iterable, Union
from .Mutation import Mutation


class MutationSet(list):
    def __init__(self, mutations: Iterable[Union[Mutation, str]]):
        self.mutations = mutations

    def __add__(self, other):
        """Only defined for integers for now, should be defined for other mutation sets!"""
        # self.mutations = [m + other for m in self.mutations]
        return _add(self, other)

    def __str__(self):
        return str(self.mutations)

    def __repr__(self):
        return str(self.mutations)

    def __iter__(self):
        return iter(self.mutations)

    def __next__(self):
        return next(self.mutations)


@dispatch(MutationSet, int)
def _add(obj: MutationSet, other: int):
    print(obj.mutations)
    obj.mutations = [m + other for m in obj.mutations]
    return obj


@dispatch(MutationSet, MutationSet)
def _add(obj: MutationSet, other: MutationSet):
    obj.mutations = obj.mutations + other.mutations
    return obj
