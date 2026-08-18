# Typing Conventions

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

