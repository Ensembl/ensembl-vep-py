# Style Conventions

## Typing Conventions

All code must use explicit type annotations for function signatures, class attributes, and variables where appropriate.

**Rules:**
- Use standard Python typing (PEP 484, PEP 604).
- For return types that may be a value or None, use `Type | None` instead of `Optional[Type]`. This is more cleaner. Also type checking does not happen at runtime, so does not affect users using python < 3.10 unless we are using type introspection dynamically at runtime.

  - Example:
    ```python
    def get_feature() -> Feature | None:
        ...
    ```
    Instead of:
    ```python
    def get_feature() -> Optional[Feature]:
        ...
    ```
- Use built-in types where possible (e.g., `list`, `dict`, `str`).
- Always annotate arguments and return values in public class, functions etc. Exceptions can be tests and private functions.

We will enforce these conventions through static analysis tools (e.g., mypy) as part of our CI pipeline to ensure consistency and catch type-related issues early in the development process.

## Logging Implementation

All logging in this project should use Python's built-in `logging` module. The recommended setup is as follows:

- `logging.basicConfig` has been used in `ensembl.vep.common.setup_logging` function to setup logging configuration. This function should be called at the entry point of the application to ensure consistent logging behavior across all modules.
- Log messages is sent to a file named `{output_file}_warnings.txt` consistent with old Ensembl VEP.
- The log format includes timestamp, log level, logger name, and message.
- For timestamps (`%Y-%m-%dT%H:%M:%S`) ISO 8601 format is used.
- Always use named loggers (`logging.getLogger(__name__)`) for module-level logging.

Example setup in a module:

```python
import logging

logger = logging.getLogger(__name__)

def some_function():
    logger.info("This is an info message.")
    logger.warning("This is a warning message.")
```

**Open question:**
Should we have a backward compatible option to have same manner of message formatting as old VEP?

## Docstring Style

All docstrings in this project must follow the [Google Python Style Guide](https://google.github.io/styleguide/pyguide.html#38-comments-and-docstrings) for docstrings.

**Rationale:**
- Google style docstrings provide clear, consistent documentation and improve readability.
- The docstring will facilitate automated documentation generation using mkdocs. 

> No docstring, no visibility in the public docs!

**Key points:**
- Use triple double quotes (`"""`) for all docstrings.
- Include sections for Args, Returns, Raises, and Examples as appropriate.
- Document types for all parameters and return values.
- Keep descriptions concise and informative.
