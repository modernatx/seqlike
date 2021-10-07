# Documentation Standards and Practices

This doc describes the practices surrounding documentation.

## Style

We use sphinx-style docstrings throughout the codebase.


## Standards

We use `interrogate` to help us maintain docstring coverage of our functions.

We also use `darglint` to help us check that our docstrings cover all arguments.
The only exception to this is functions that are `@dispatch`-ed across different types.
There, to avoid text duplication,
we use `#noqa: DAR101` and `#noqa: DAR201` inside the docstring
to tell `darglint` to ignore those functions.

!!! danger "Exceptions must be documented!"

    `DAR401`, which catches any undocumented exceptions raised, cannot be `#noqa`'d!

!!! note "noqa must be placed at the end of a docstring."

    The `#noqa: <ERRORCODE>` must be placed at the _end_ of a docstring,
    otherwise it will get rendered as Markdown inside the API docs.
