# Development Containers

## Why development containers?

This is our recommended way of getting set up for development,
because it'll eliminate the vast majority of development environment quriks
that are difficult to debug.
There are other reasons to use development containers
that are detailed in a variety of blog posts linked below.
If you're new to the concept,
we recommend that you give them a read.

- [An introduction to local development with containers](https://increment.com/development/an-introduction-to-local-development-with-containers/)
- [Container Driven Development](https://hackernoon.com/container-driven-development-1bd08c2f126d)

Microsoft has an awesome tutorial that you should walk through
before you attempt to use development containers.
The tutorial is available [here](https://docs.microsoft.com/en-us/learn/modules/use-docker-container-dev-env-vs-code/).

To use development containers, we recommend using Visual Studio Code,
which has the best support for development containers
amongst the IDEs that we have encountered.

## Getting started with SeqLike's development containers

The rough order of steps is as follows:

0. Ensure you have Docker running on your computer.
1. Clone the repository locally, preferrably using SSH.
2. Ensure that you have your [SSH agent running][ssh].
3. Open the repository inside VSCode.
4. (If you don't have it) Install the VSCode Remote extension.
5. When prompted, re-open the repository inside the Dev container we have specified.

[ssh]: https://code.visualstudio.com/docs/remote/containers#_using-ssh-keys

Once that is done, verify that you have everything correctly installed
by executing the following commands to run the test suite:

```bash
pytest
```

Also, you can build the docs:

```bash
mkdocs serve
```

You should then be able to go to `http://localhost:8000` to view the docs.

## Features of our development container

### Automated development install

SeqLike is automatically installed in development mode
thus allowing you to test code changes as soon as they are made.

### Python development tools

We installed the following VSCode extensions:

- Pylance (`ms-python.vscode-pylance`)
- Jupyter (`ms-toolsai.jupyter`)
- MS Python (`ms-python.python`)

### Pre-commit hooks

To ensure code quality, pre-commit hooks are installed into the environment.
You'll always be notified if your code needs to be modified
because they didn't pass our automated code checks.

### Nord theme

We give you a beautiful and aesthetically pleasing theme for VSCode.
Rest assured, this doesn't affect whatever your original theme settings are.

### Port 8000

We use `mkdocs` to build the docs.
`mkdocs` comes with a built-in preview server
which runs on port 8000 when you execute `mkdocs serve`.
That port is forwarded to your `localhost`.
