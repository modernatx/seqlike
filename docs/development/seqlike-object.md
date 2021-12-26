# The SeqLike Object

## Purpose of this document

The purpose of this document is to accurately document the SeqLike object's intended behaviour.
It is not a replacement for docstrings,
which serve as a reference for the arguments of the SeqLike methods.

You should read this document if:

1. You are interested in understanding the internals of how the SeqLike object is organized.
2. You wish to leverage that knowledge to suggest improvements.

What is explicitly out-of-scope for this document are:

1. Places that we need help (please consult code for TODOs)

## SeqLike Object Overview

The SeqLike object aims to be an omnibus object
that can accept any kind of biological sequence-like Python data types.
This means strings, BioPython `Seq` or `SeqRecord` objects,
NumPy arrays (and their almost equivalent PyTorch tensors),
and more importantly, other SeqLike objects.

Internal to SeqLike objects are a few key attributes
that control their behaviour.

Firstly, we maintain a dual representation of nucleotides and amino acids
as BioPython SeqRecord objects. These are `._aa_record` and `._nt_record`,
as well as a pointer `._seqrecord` that points to only one of the two.

Secondly, we maintain state identifying whether a SeqLike object
is an _amino acid_ `_type` or a _nucleotide_ `_type`.
(Though we use the term "type" here,
please don't confuse it with Python object types!)

Thirdly, there is an `alphabet` attribute.

The `alphabet` identifies the set of valid characters
that can exist inside the sequence,
Order matters for the alphabet!

There is also an `_index_encoder` and `_onehot_encoder` attribute that gets set,
which allow us to easily take a string sequence
and convert it into a starter numerical representation for downstream purposes
(e.g. as inputs to ML models).

Finally, there is a `codon_map` attribute that gets set.
Setting this allows us to perform back translation
from amino acid sequences to nucleotid sequences easily.

The main arguments to the SeqLike object constructor are:

- `sequence`: A SeqLike Type object.
- `seq_type`: An optional argument that sets the primary `._type` of the SeqLike object.
- `alphabet`: The set of valid characters for the biological sequence.
We provide a standard set that users can use
in the top-level `seqlike` namespace and strongly recommend sticking to them.
- `codon_map`: An optional argument that enables back-translation.

In summary, the main attributes are:

- `._aa_record`: A BioPython SeqRecord object.
- `._nt_record`: A BioPython Seqrecord object.
- `._seqrecord`: A BioPython SeqRecord object,
which is always a pointer to one of `._aa_record` and `._nt_record`.
- `._type`: Primary sequence type. A string, one of `("NT", "AA")`
- `.alphabet`: A string of valid characters for the sequence.
- `._index_encoder`/`._onehot_encoder`: An encoder function
for converting the string into a numerical representation.
- `.codon_map`: A callable (i.e. function)
that accepts and returns an AA SeqLike object.

## How the constructor works

The `seqlike` package provides the `SeqLike` class,
which is used for constructing SeqLike objects.
In the constructor, the `seq_type` argument is required
in order to avoid incorrectly inferring whether a sequence
is a nucleotide or amino acid sequence.
For convenience, there are also `aaSeqLike` and `ntSeqLike` factory functions
that return a SeqLike constructed with the appropriate `seq_type` specified.
Inside the SeqLike object's constructor,
we use a dispatching pattern that returns a tuple of attributes
based on the Python object types that are passed in.
This design choice is highly inspired by the dispatching pattern
prevalent in the Julia language,
which makes reasoning about the behaviour of the constructor much easier.
It also makes the code flatter and easier to read.
(Zen of Python, friends!)
The main function that determines SeqLike attribute values
is called `_construct_seqlike`.
At a high level, this is what happens inside `_construct_seqlike`:

1. We determine the `_type` and `alphabet` of the sequence based on the passed-in `seq_type`, `alphabet`, and `sequence` arguments.
2. We then identify the appropriate `_index_encoder` and `_onehot_encoder` to attach to the object based on the `alphabet`.
3. We then construct a BioPython SeqRecord object based on the `sequence`, `_index_encoder`, and `_onehot_encoder`.
4. Finally, we return the SeqLike object attributes to be set.

Using this pattern ensures consistency in the way we reason about SeqLike objects.
By using the dispatching pattern,
we can also reason clearly about combinations of cases
and what their intended behaviour ought to be.

Coming up, let's see how some of these cases are handled.

### God-Mode: When the sequence is a SeqLike object

This mode takes top priority.
If a SeqLike object is passed into the SeqLike constructor,
such as what happens in the `.nt()` and `.aa()` methods,
then the dispatched constructor function simply copies all relevant attributes
and returns a tuple of them to be set in the object constructor.

Doing so helps with downstream logic;
we never have to worry about the SeqLike object type
in determining how to construct the object.

<!-- TODO: This has to be implemented! -->

### Determine the Sequence Type

The `._type` attribute semantically maps to the "sequence type",
i.e. an amino acid sequence or nucleotide sequence.
This attribute can never be `None`;
**it must be either "AA" or "NT".**
<!-- (NOTE: This should be tested!) -->
We determine `._type` using two pieces of information:
the `seq_type` argument and the `sequence` argument.

Here are the cases that we consider:

1. If `seq_type` is a string, then `._type` is set to `seq_type`.
2. If `seq_type` is `None`, then `._type` is inferred using the `sequence` argument.

### Determine the alphabet

Next up, we need to know some information about the sequence's alphabet
in order to determine the encoder mode and the relevant encoders.
Here are a few cases that we handle:

If the alphabet is user-supplied,
our primary assumption is that it is semantically valid for the user's applications.
_We don't bother validating the alphabet here!_
Otherwise, the alphabet is inferred using `seq_type` and `sequence`;
we do our best to map it to one of our supplied alphabets.
We think in most cases, you can simply _avoid specifying the alphabet_
and trust that our alphabet inference code will do the right thing.

### Create the encoder objects

Next up, we need to obtain the encoders!
We maintain dual encoders in SeqLike objects,
one for index encoding and one for one-hot encoding.
These allow us to get array representations of the sequence.
The alphabet order is super duper important here;
it will determine the ordering of one-hot and index encodings.
Here's an example:
if your alphabet is "AUGC" and you ask for an index encoding,
the mappings will be A->0, U->1, G->2 and C->3.
However, if your alphabet is "AGCU", then
you'll still have A->0, but G->1, C->2 and U->3.
Try to be sane and stick to the alphabets provided by us :).
The encoders, `_onehot_encoder` and `_index_encoder`,
are then both directly created off the alphabet.

### Handle seqrecord creation

Now, we finally handle the sequence object
that will become the _aa_record and/or _nt_record.
The first thing we need to check is if the object type of `sequence` passed in
is an array-like thing or not.
If it is, we need to convert it back into a `str` type by way of the encoders
and their `inverse_transform` methods.
Now, we're ready to turn the `sequence` into SeqRecord objects.

There is a hierarchy here:

- `ArrayType` objects are converted into `str` objects,
which are added to a `Seq` object inside a `SeqRecord`.
- `str` objects are deep copied and added to a `Seq`
that is inside a `SeqRecord` object.
- `Seq` objects are deep copied and added to a `SeqRecord` object.
- `SeqRecord` objects are directly deep copied.

We can decompose this into a logical information flow
that roughly follows our constructor implementation:

```
ArrayType -> str -> Seq -> SeqRecord
```

In essence, if we start with an `ArrayType`,
then we call the ArrayType-dispatched function,
and return the string representation.
We then take the string-dispatched function
and return the Seq representation.
And so on and so-forth.
All of this is handled by dispatching on a single function `record_from`.
(If you're not familiar with what dispatching is,
this is inspired from the Julia language.
We use the `multipledispatch` library to handle dispatching,
check out its repository [here](https://github.com/mrocklin/multipledispatch)
to read more about it.)

By using dispatching, we ensure that whatever data type is passed in for `sequence`
will be handled by the correct function at the correct step,
while also minimizing code duplication inside the constructor.
Now, while we haven't profiled this just yet,
there _may_ be some performance penalty for the nested calls.
<!-- TODO: Do a profiling exercise? -->
Canonically, though, you might be most used to `SeqLikes`
being one of `str`, `Seq`, or `SeqRecord`,
and that's perfectly logical.
Regardless, you can rest assured that if you pass in any of those three,
the dispatched function will get called at the right place.

Finally, there are other `kwargs` that are passed into the constructor -
apart from the `_index_encoder` and `_onehot_encoder`
that are necessary for converting `ArrayType`s,
these get fed directly into the `SeqRecord` constructor.

## Conversion between NT and AA

Converting between the nucleotide and amino acid representations of a sequence
is one of the core features that the `SeqLike` object brings to the table.
(We show in the tutorial notebooks the kind of situations in which
you might want to switch between the two representations.)
As you, the biosequence practitioner might know,
there are tricky situations that we have to worry about,
so we'll go through some of those in detail here.

### The semantics of `.aa()` and .`nt()` vs. `.translate()` and `.back_translate()`

You might be tempted to think that
`.aa()` and `.nt()` mean `.translate()` and `.back_translate()`.
Thankfully (and we mean it!), this is not the case.
Allow us to explain why.
What `.nt()` or `.aa()` are supposed to do
is return a deepcopy of the `SeqLike` object
with the appropriate `._type`,
`._index_encoder`,
`._onehot_encoder`,
and `.alphabet` attributes set correctly.
It doesn't explicitly do the `.translate()` or `.back_translate()` operation.
You can think of `.nt()` and `.aa()` as switching between views of the same sequence.

By contrast, `.translate()` and `.back_translate()`
are intended to actively (re)compute the internal representations
`._nt_record` and `._aa_record`.
Back translation is especially tricky - for the purposes of engineering mRNA,
you might want to have a _stochastic_ back-translation
so that you can generate multiple variants _in silico_.
Being Pythonistas,
we know from the Zen of Python that _explicit is better than implicit_.
As such, we believe that changing the internal representation
is something that should be _explicitly_ requested
rather than _implicitly_ done when switching representations.

With all of that said, there _is_ one asymmetry
that we would like to point out.
It is easy to go from a nucleotide representation to an amino acid representation,
if we assume a codon table is present.
It's _not_ easy to go backwards, because the process is probabilistic:
multiple NT representations could be valid for the same AA representation.
Regardless, `.aa()` has an argument, `auto_translate`,
which is set to `True` by default,
and `.nt()` has an argument, `auto_back_translate`,
which is set to `True` by default.
Consider it our way of providing magic to ourselves and to our users!

## SeqRecord-like behaviour with SeqLike attributes

For users who are used to BioPython SeqRecord objects,
we have implemented a custom `__getattr__`
that delegates attribute access to the underlying SeqRecord object
that is stored in either the `._nt_record` (for "NT" SeqLikes) attribute
or the `._aa_record` (for "AA" SeqLikes) attribute
when one asks for a SeqRecord attribute.

The overall effect here is that we can access SeqRecord attributes
_as if_ they were native to SeqLike objects even though their access is delegated,
such as accessing the SeqRecord's `description` or `dbxrefs`.

```python
s = SeqLike("AMPELQYTV", seq_type="AA", alphabet=STANDARD_AA)
s.description  # a SeqRecord attribute accessed on a SeqLike object.
s.dbxrefs      # also a SeqRecord attribute
```

## Slicing Behaviour

When we slice into a SeqLike object,
we basically are slicing into the underlying SeqRecord object.
To do so, we have implemented a custom `__getitem__` method
that guarantees correct slicing behaviour w.r.t.:

1. the underlying SeqRecord object(s) (`._nt_record` and `._aa_record`), and
2. the SeqRecord `.letter_annotations` dictionaries,

## Code Notes

### Validation functions

Our validation functions should _never, ever, ever_ return anything.
They should error loudly if the validation check errors out or else be silent.
We will be utterly confused reading the code
if a validation function ever returns anything,
so please don't do that!
