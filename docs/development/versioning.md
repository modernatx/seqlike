# Versioning

We use bump2version to help us manage versioning.

## But why?

Isn't versioning just about updating version numbers?

Well, kind of.

You see, GitHub has additional functionality
to help us archive versions of our software package.
This is an added benefit in addition to the releases that we cut on PyPI.
We essentially need to use git tags to help manage the process.
`bump2version` provides the tooling necessary
to automate the creation of git tags whenever we bump the version of SeqLike.

## What do we do when we need to cut a new release?

OK, so I'm guessing you're convinced of the need to have another tool in the toolbox.
:)

Let's talk about how to cut a new release, then.

Firstly, go to the GitHub repository's [actions].

[actions]: https://github.com/modernatx/seqlike/actions

Then, click on the workflow "Release".

Then, click on the dropdown named "Run workflow" (on the right side).

Finally, identify whether the bump is a major, minor, or patch release
and click "Run workflow"!

## Magical! How did you make that happen?

We essentially take advantage of GitHub actions' workflow triggers.

If you're curious to see how exactly the action is structured,
head over to the file `.github/workflows/release.yml`.
There are highly readable notes and section headers
that explain what logical order of steps are happening.
