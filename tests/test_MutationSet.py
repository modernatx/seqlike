from seqlike import aaSeqLike
from seqlike.Mutation import Mutation
from seqlike.MutationSet import MutationSet


def test_MutationSet():
    sub1 = Mutation("3K")
    ins1 = Mutation("^4D")
    del1 = Mutation("2-")

    muts_list = [sub1, ins1, del1]
    muts = MutationSet(muts_list)
    assert len(muts) == len(muts_list) == 3
