"""Data validation functions.

As far as possible, we try to put them here.
Unfortunately, sometimes we may run into circular import issues.
In that case, place the validation functions where they need to be
to avoid circular imports.
"""


def validate_seq_type(seq_type: str):
    """Validate the user provided seq_type argument.
    :param seq_type: User provided seq_type.
    :raises ValueError: Raises value error if the user provided seq_type is not one of "DNA", "RNA", "NT", "AA"
    """
    valid_seq_types = ("DNA", "RNA", "NT", "AA")
    err_msg = (
        f"seq_type must be one of the following: {valid_seq_types}."
        "The argument is case-insensitive, so you can pass in lower-case if you'd like."
    )
    if not seq_type in valid_seq_types:
        raise ValueError(err_msg)
