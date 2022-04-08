# Optimizing for Import Speed

As SeqLike forms the basis of downstream tooling,
import timing becomes important as well.
As such, we have implemented a few tricks to enable fast importing of the library.

## Lazy Loading

In accordance with [SPEC-001][spec],
we use `lazy_loader` to lazily load Python packages into our modules.
Doing so ensures that we do not incur the time penalty
associated with importing large packages, such as:

- pandas
- NumPy
- matplotlib
- bokeh

[spec]: https://scientific-python.org/specs/spec-0001/

The general pattern we use here is to do:

```python
import lazy_loader as lazy

pd = lazy.load("pandas") # equivalent to `import pandas as pd`
```

This pattern is generally useful only for top-level packages.

## Nesting imports

The second strategy for deferring imports
is to import them within functions and class methods.
We do this when there are submodules
or specific functions/classes that we need to access.
One example is Bokeh:

```python
def some_bokeh_function():
    from bokeh.plotting import figure
    ...
    p1 = figure(...)
```

## Try/Excepts

The third pattern we use is the try/except pattern.
This one is particularly useful
for distinguishing between notebook and shell environments.
For example, in `draw_utils.py`, we use the following pattern:

```python
try:
    get_ipython
    from bokeh.io import output_notebook

    output_notebook()
```

This way, we avoid executing `output_notebook()`
if we are in a CLI/shell environment and not in a Jupyter notebook environment.
