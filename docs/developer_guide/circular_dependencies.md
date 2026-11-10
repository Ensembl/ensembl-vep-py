# Avoiding circular dependency

The heavy use of subfolders to group similar objects means that, depending on the way we structure our modules, there could be circular dependencies during import. The main source of circular dependencies here comes from using `__init__.py` files to expose a package's public API.

The circular import is avoided using the following approaches:

## When importing from same package
When importing within the same package, always import directly from the concrete module rather than through a package's public API.

**Prefer:**
```python
from ensembl.vep.models.genic.transcript import Transcript
```
**Over:**
```python
from ensembl.vep.models.genic import Transcript
```

We also prefer absolute imports over relative imports, even within the same package:
**Prefer:**
```python
from ensembl.vep.models.genic.transcript import Transcript
```
**Over:**
```python
from .transcript import Transcript
```
Absolute imports are more explicit, easier to trace, and consistent with [PEP 8](https://peps.python.org/pep-0008/#imports).

## Importing for type annotations
When we are importing some module that is solely used for type annotations, `typing.TYPE_CHECKING` should be used -
```python
if typing.TYPE_CHECKING:
    from some_module import some_function
```

It allows importing the module only when type checking is performed (e.g. by static analysis tools), but not during runtime.

## Prevent direct circular import
And, of course, we need to prevent direct circular imports (module A importing module B, and module B importing module A) by making sure we keep one-way dependency direction, i.e. the base class cannot inherit from any of the descendant classes.
