# See the NOTICE file distributed with this work for additional information
# regarding copyright ownership.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
"""Annotated-variant output writer registry."""

from __future__ import annotations

from collections.abc import Callable

from ensembl.vep.output.factory.write.base import OutputWriter
from ensembl.vep.output.factory.write.json import JSONOutputWriter
from ensembl.vep.output.factory.write.vcf import VCFOutputWriter

__all__ = [
    "available_output_writers",
    "get_output_writer",
    "register_output_writer",
]

OutputWriterFactory = Callable[[str], OutputWriter]
"""Build an output writer using the run's command line."""

_OUTPUT_WRITER_REGISTRY: dict[str, OutputWriterFactory] = {}


def register_output_writer(name: str) -> Callable[[OutputWriterFactory], OutputWriterFactory]:
    """Register an output writer factory under a format name."""

    def decorator(factory: OutputWriterFactory) -> OutputWriterFactory:
        _OUTPUT_WRITER_REGISTRY[name] = factory
        return factory

    return decorator


def available_output_writers() -> tuple[str, ...]:
    """Return the sorted set of registered output formats."""
    return tuple(sorted(_OUTPUT_WRITER_REGISTRY))


def get_output_writer(output_format: str, *, command_line: str = "") -> OutputWriter:
    """Return the writer for one output format.

    Raises:
        ValueError: If no writer supports ``output_format``.
    """
    try:
        factory = _OUTPUT_WRITER_REGISTRY[output_format]
    except KeyError as error:
        supported_formats = ", ".join(available_output_writers()) or "<none>"
        raise ValueError(f"no writer for {output_format!r}; supported: {supported_formats}") from error
    return factory(command_line)


@register_output_writer("json")
def _json_writer(_command_line: str) -> OutputWriter:
    """Build the JSON-lines writer."""
    return JSONOutputWriter()


@register_output_writer("vcf")
def _vcf_writer(command_line: str) -> OutputWriter:
    """Build the VCF writer with its command-line header."""
    return VCFOutputWriter(command_line=command_line)
