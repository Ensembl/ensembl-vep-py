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
"""The registered corpus and parameters to construct."""

from __future__ import annotations

import hashlib
import json
from collections.abc import Mapping, Sequence
from typing import TypeAlias

from attrs import frozen

__all__ = ["CorpusParameter", "CorpusSpec"]

CorpusParameter: TypeAlias = str | int | float | bool | None
"""A scalar parameter accepted by a registered corpus factory."""

_DEFAULT_PARAMETERS: dict[str, CorpusParameter] = {
    "fuzz": 4000,
    "step": 13,
    "genome_seed": 7,
    "boundary_seed": 11,
    "fuzz_seed": 23,
    "heavy_genome_seed": 31,
    "sweep_seed": 41,
}

_SYNTHETIC_KINDS = ("synthetic.default", "synthetic.heavy", "synthetic.all")


def _parameter_items(
    parameters: Mapping[str, CorpusParameter],
) -> tuple[tuple[str, CorpusParameter], ...]:
    """Freeze parameters in stable key order."""
    return tuple(sorted(parameters.items()))


def _parameter_value(raw: str) -> CorpusParameter:
    """Decode a command-line JSON scalar, retaining ordinary strings."""
    try:
        decoded = json.loads(raw)
    except json.JSONDecodeError:
        return raw
    if decoded is None or isinstance(decoded, (str, int, float, bool)):
        return decoded
    raise ValueError(f"corpus parameter values must be scalars, got {raw!r}")


def _option_parameters(options: object) -> dict[str, CorpusParameter]:
    """Return generic parameters declared by a config-like object."""
    configured = getattr(options, "parameters", None)
    if configured is None:
        return {}
    if isinstance(configured, Mapping):
        mapping_parameters: dict[str, CorpusParameter] = {}
        for name, value in configured.items():
            if not isinstance(name, str):
                raise TypeError("corpus parameter names must be strings")
            if value is not None and not isinstance(value, (str, int, float, bool)):
                raise TypeError(f"corpus parameter {name!r} must be a scalar")
            mapping_parameters[name] = value
        return mapping_parameters
    if isinstance(configured, Sequence) and not isinstance(configured, str):
        sequence_parameters: dict[str, CorpusParameter] = {}
        for item in configured:
            if not isinstance(item, str) or "=" not in item:
                raise ValueError(f"corpus parameter must be KEY=VALUE, got {item!r}")
            name, raw_value = item.split("=", 1)
            if not name:
                raise ValueError(f"corpus parameter needs a name, got {item!r}")
            sequence_parameters[name] = _parameter_value(raw_value)
        return sequence_parameters
    raise TypeError("corpus parameters must be a mapping or a sequence of KEY=VALUE strings")


@frozen(init=False, kw_only=True, auto_attribs=True)
class CorpusSpec:
    """Describe a registered corpus construction.

    ``parameters`` is deliberately open: each registered factory owns its parameter validation.
    Synthetic convenience arguments remain available alongside generic factory parameters.

    Attributes:
        kind (str): Registered corpus name.
        parameters (Mapping[str, CorpusParameter]): Factory-owned construction parameters.
    """

    kind: str
    _parameter_items: tuple[tuple[str, CorpusParameter], ...]

    def __init__(
        self,
        *,
        kind: str = "synthetic.all",
        parameters: Mapping[str, CorpusParameter] | None = None,
        fuzz: int | None = None,
        step: int | None = None,
        genome_seed: int | None = None,
        boundary_seed: int | None = None,
        fuzz_seed: int | None = None,
        heavy_genome_seed: int | None = None,
        sweep_seed: int | None = None,
    ) -> None:
        """Create a specification, accepting synthetic parameter overrides.

        Explicit synthetic arguments replace values from ``parameters``; unspecified synthetic
        parameters retain their defaults.

        Args:
            kind (str): Registered corpus name.
            parameters (Mapping[str, CorpusParameter] | None): Factory-owned scalar parameters.
            fuzz (int | None): Default synthetic random-variant count override.
            step (int | None): Heavy synthetic variant spacing override.
            genome_seed (int | None): Default synthetic genome seed override.
            boundary_seed (int | None): Boundary-variant seed override.
            fuzz_seed (int | None): Random-variant seed override.
            heavy_genome_seed (int | None): Heavy synthetic genome seed override.
            sweep_seed (int | None): Heavy synthetic variant seed override.

        Raises:
            TypeError: If a parameter name or value is not a supported scalar.
            ValueError: If the kind is empty or a built-in count or step is invalid.
        """
        if not kind:
            raise ValueError("corpus kind must not be empty")

        resolved = dict(_DEFAULT_PARAMETERS) if kind in _SYNTHETIC_KINDS else {}
        if parameters is not None:
            resolved.update(parameters)
        overrides = {
            "fuzz": fuzz,
            "step": step,
            "genome_seed": genome_seed,
            "boundary_seed": boundary_seed,
            "fuzz_seed": fuzz_seed,
            "heavy_genome_seed": heavy_genome_seed,
            "sweep_seed": sweep_seed,
        }
        for name, override_value in overrides.items():
            if override_value is not None:
                resolved[name] = override_value
        for name, parameter_value in resolved.items():
            if not isinstance(name, str):
                raise TypeError("corpus parameter names must be strings")
            if parameter_value is not None and not isinstance(parameter_value, (str, int, float, bool)):
                raise TypeError(f"corpus parameter {name!r} must be a scalar")

        if kind in _SYNTHETIC_KINDS:
            resolved_fuzz = resolved.get("fuzz")
            resolved_step = resolved.get("step")
            if not isinstance(resolved_fuzz, int) or isinstance(resolved_fuzz, bool) or resolved_fuzz < 0:
                raise ValueError("fuzz must be a non-negative integer")
            if not isinstance(resolved_step, int) or isinstance(resolved_step, bool) or resolved_step <= 0:
                raise ValueError("step must be a positive integer")

        object.__setattr__(self, "kind", kind)
        object.__setattr__(self, "_parameter_items", _parameter_items(resolved))

    @classmethod
    def from_options(cls, options: object) -> CorpusSpec:
        """Build a specification from corpus configuration attributes.

        Args:
            options (object): Configuration-like object with ``kind`` and corpus parameters.

        Returns:
            CorpusSpec: Normalised construction specification.

        Raises:
            TypeError: If configured generic parameters have invalid types.
            ValueError: If generic parameters are malformed or built-in values are invalid.
        """
        kind = getattr(options, "kind", "synthetic.all")
        parameters: dict[str, CorpusParameter] = {}
        corpus_directory = getattr(options, "corpus_dir", None)
        if corpus_directory is not None and kind == "loaded.files":
            parameters.update(
                {
                    "fasta": str(corpus_directory / "corpus.fa"),
                    "gff3": str(corpus_directory / "corpus.gff3"),
                    "vcf": str(corpus_directory / "corpus.vcf"),
                }
            )
        elif corpus_directory is not None and kind == "golden.file":
            parameters["directory"] = str(corpus_directory)
        parameters.update(_option_parameters(options))
        if kind in _SYNTHETIC_KINDS:
            for name in _DEFAULT_PARAMETERS:
                if hasattr(options, name):
                    parameters[name] = getattr(options, name)
        return cls(kind=kind, parameters=parameters)

    @classmethod
    def from_dict(cls, value: Mapping[str, object]) -> CorpusSpec:
        """Build a specification from its JSON-compatible representation.

        Args:
            value (Mapping[str, object]): Mapping produced by ``as_dict``.

        Returns:
            CorpusSpec: Decoded construction specification.

        Raises:
            TypeError: If the kind or parameters have invalid types.
            ValueError: If decoded built-in parameters are invalid.
        """
        kind = value.get("kind", "synthetic.all")
        if not isinstance(kind, str):
            raise TypeError("corpus kind must be a string")
        raw_parameters = value.get("parameters")
        if raw_parameters is None:
            parameters = {name: value[name] for name in _DEFAULT_PARAMETERS if name in value}
        elif isinstance(raw_parameters, Mapping):
            parameters = dict(raw_parameters)
        else:
            raise TypeError("corpus parameters must be a mapping")
        return cls(kind=kind, parameters=parameters)  # type: ignore[arg-type]

    @property
    def parameters(self) -> dict[str, CorpusParameter]:
        """Return a copy of the registered factory parameters."""
        return dict(self._parameter_items)

    def parameter(self, name: str, default: CorpusParameter = None) -> CorpusParameter:
        """Return one factory parameter or ``default`` when absent.

        Args:
            name (str): Parameter name.
            default (CorpusParameter): Value returned when the parameter is absent.

        Returns:
            CorpusParameter: Configured or default scalar value.
        """
        return self.parameters.get(name, default)

    def with_parameters(self, **parameters: CorpusParameter) -> CorpusSpec:
        """Return a specification with selected parameters replaced.

        Args:
            parameters (CorpusParameter): Scalar values to add or replace.

        Returns:
            CorpusSpec: New specification retaining the same corpus kind.
        """
        merged = self.parameters
        merged.update(parameters)
        return CorpusSpec(kind=self.kind, parameters=merged)

    def as_dict(self) -> dict[str, object]:
        """Return the JSON-compatible representation."""
        return {"kind": self.kind, "parameters": self.parameters}

    def _integer_parameter(self, name: str) -> int:
        """Return a synthetic integer parameter with its type checked."""
        value = self.parameters[name]
        if not isinstance(value, int) or isinstance(value, bool):
            raise TypeError(f"corpus parameter {name!r} must be an integer")
        return value

    @property
    def fuzz(self) -> int:
        """Return the default synthetic corpus random-variant count."""
        return self._integer_parameter("fuzz")

    @property
    def step(self) -> int:
        """Return the heavy-locus sweep step."""
        return self._integer_parameter("step")

    @property
    def genome_seed(self) -> int:
        """Return the default synthetic genome seed."""
        return self._integer_parameter("genome_seed")

    @property
    def boundary_seed(self) -> int:
        """Return the boundary-variant seed."""
        return self._integer_parameter("boundary_seed")

    @property
    def fuzz_seed(self) -> int:
        """Return the random-variant seed."""
        return self._integer_parameter("fuzz_seed")

    @property
    def heavy_genome_seed(self) -> int:
        """Return the heavy-locus genome seed."""
        return self._integer_parameter("heavy_genome_seed")

    @property
    def sweep_seed(self) -> int:
        """Return the heavy-locus sweep seed."""
        return self._integer_parameter("sweep_seed")

    @property
    def kinds(self) -> tuple[str, ...]:
        """Return the concrete corpus kinds represented by the specification."""
        if self.kind == "synthetic.all":
            return ("synthetic.default", "synthetic.heavy")
        return (self.kind,)

    def seeds_by_kind(self) -> dict[str, dict[str, int]]:
        """Return synthetic seeds for the selected built-in corpora."""
        seeds: dict[str, dict[str, int]] = {}
        if "synthetic.default" in self.kinds:
            seeds["synthetic.default"] = {
                "genome": self.genome_seed,
                "boundary": self.boundary_seed,
                "fuzz": self.fuzz_seed,
            }
        if "synthetic.heavy" in self.kinds:
            seeds["synthetic.heavy"] = {
                "genome": self.heavy_genome_seed,
                "sweep": self.sweep_seed,
            }
        return seeds

    def parameter_hash(self) -> str:
        """Return a short stable digest of the requested construction."""
        payload = json.dumps(self.as_dict(), sort_keys=True).encode("utf-8")
        return hashlib.sha256(payload).hexdigest()[:16]
