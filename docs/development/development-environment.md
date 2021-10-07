# Setting up your development environment manually

In order to hack on SeqLike, you'll want to have your development environment set up.
This document can serve as a tutorial document that shows you how to get setup.

## Assumed knowledge

We'll throw this point out first: this is not really a beginner-friendly document.
If you're here, we're assuming some knowledge:

1. You know about virtual environments and their flavours (`venv`, `conda`, etc.)
2. You're comfortable setting up virtual environments on your own.
3. The doc is `conda`-heavy, so familiarity with `conda` commands will be helpful.

There are usually multiple valid ways to set up an environment.
If you're looking for a pretty fool-proof way to set up development environments
without fiddling with all sorts of stuff,
try out [development containers][devcontainer].

[devcontainer]: /development/devcontainer/

## Steps

### tl;dr

Copy/paste this exactly _only if_ you feel comfortable with the commands
and know exactly what they're doing.

```bash
conda env create -f environment.yml
pip install -e .
pytest
```

### Set up conda environment

We provide an `environment.yml` file that should get you started.
Create the environment using the `conda env create` command:

```bash
conda env create -f environment.yml
```

### Install the package in development mode

To hack on the package, and, more importantly, to run tests,
you'll want to install SeqLike in development mode:

```bash
pip install -e .
```

### Run tests

The surest way of knowing whether your environment is installed correctly
is to run the test suite straight away.
You can do this by running:

```bash
pytest
```

### Build docs

This is another good way to test whether your environment is installed correctly or not.

```bash
mkdocs serve
```

You should be able to see the docs previewed at port 8000 on your localhost.

### Update conda environment

If you ever need to update the environment, say, to obtain new package versions:

```bash
conda env update -f environment.yml
```
