# Logging Implementation

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

