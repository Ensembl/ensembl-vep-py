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
"""Boolean expression tree for consequence rules."""

from __future__ import annotations

from typing import TYPE_CHECKING

from attrs import field, frozen

from ensembl.vep.annotators.consequence.facts.base import get_fact, total_cost

if TYPE_CHECKING:
    from ensembl.vep.annotators.consequence.context.fact import FactContext

__all__ = ["And", "Atom", "Expr", "Not", "Or", "all_of", "any_of", "expr_from_json", "no"]


class Expr:
    """Represent a Boolean expression over fact atoms."""

    def atoms(self) -> set[str]:
        """Return every fact atom referenced by this expression.

        Returns:
            set[str]: The referenced fact atom names.
        """
        raise NotImplementedError

    def evaluate(self, fact_context: FactContext) -> bool:
        """Evaluate this expression against a fact context, computing facts as needed.

        Args:
            fact_context (FactContext): The overlap to evaluate against.

        Returns:
            bool: The expression's truth value.
        """
        raise NotImplementedError

    def cost(self) -> int:
        """Return the cold cost of evaluating this expression.

        Returns:
            int: The total cost of every atom referenced, counting each shared dependency once.
        """
        return total_cost(self.atoms())

    def evaluate_partial(self, known: dict[str, bool]) -> bool | None:
        """Evaluate this expression against a partial (Kleene) assignment.

        Returns True/False when the atoms in ``known`` already settle the formula, and None when
        they do not. This is what lets an annotator decide a rule from cheap facts alone and
        never touch the expensive tiers.

        Args:
            known (dict[str, bool]): A partial assignment of fact atom names to boolean values.

        Returns:
            bool | None: The expression's truth value if settled by ``known``, else None.
        """
        raise NotImplementedError

    def __and__(self, other: Expr) -> Expr:
        return And(self, other)

    def __or__(self, other: Expr) -> Expr:
        return Or(self, other)

    def __invert__(self) -> Expr:
        return Not(child=self)


@frozen(kw_only=False, auto_attribs=True)
class Atom(Expr):
    """Represent one fact atom by name."""

    name: str

    def atoms(self) -> set[str]:
        """Return this atom's name."""
        return {self.name}

    def evaluate(self, fact_context: FactContext) -> bool:
        """Return this atom's value in a fact context."""
        return bool(get_fact(fact_context, self.name))

    def evaluate_partial(self, known: dict[str, bool]) -> bool | None:
        """Return this atom's known value, if assigned."""
        return known.get(self.name)

    def __repr__(self) -> str:
        return self.name


@frozen(kw_only=False, auto_attribs=True)
class Not(Expr):
    """Negate one child expression."""

    child: Expr

    def atoms(self) -> set[str]:
        """Return the atoms referenced by the child expression."""
        return self.child.atoms()

    def evaluate(self, fact_context: FactContext) -> bool:
        """Return the negated child value in a fact context."""
        return not self.child.evaluate(fact_context)

    def evaluate_partial(self, known: dict[str, bool]) -> bool | None:
        """Return the negated child value when it is known."""
        value = self.child.evaluate_partial(known)
        return None if value is None else not value

    def __repr__(self) -> str:
        return f"!{self.child!r}"


@frozen(kw_only=False, auto_attribs=True, init=False)
class And(Expr):
    """Combine child expressions with logical conjunction."""

    children: tuple[Expr, ...] = field()

    def __init__(self, *children: Expr) -> None:
        object.__setattr__(self, "children", tuple(children))

    def atoms(self) -> set[str]:
        """Return the atoms referenced by all child expressions."""
        return set().union(*(child.atoms() for child in self.children)) if self.children else set()

    def evaluate(self, fact_context: FactContext) -> bool:
        """Return whether every child is true in a fact context."""
        # Cheapest conjunct first: a false cheap test skips the expensive ones.
        # TBD: check by cost() or closure cost()? The latter is more expensive to compute but
        # should be more accurate.
        for child in sorted(self.children, key=lambda child: child.cost()):
            if not child.evaluate(fact_context):
                return False
        return True

    def evaluate_partial(self, known: dict[str, bool]) -> bool | None:
        """Return the conjunction when the known values settle it."""
        unknown = False
        for child in self.children:
            value = child.evaluate_partial(known)
            if value is False:
                return False
            if value is None:
                unknown = True
        return None if unknown else True

    def __repr__(self) -> str:
        return "(" + " & ".join(repr(child) for child in self.children) + ")"


@frozen(kw_only=False, auto_attribs=True, init=False)
class Or(Expr):
    """Combine child expressions with logical disjunction."""

    children: tuple[Expr, ...] = field()

    def __init__(self, *children: Expr) -> None:
        object.__setattr__(self, "children", tuple(children))

    def atoms(self) -> set[str]:
        """Return the atoms referenced by all child expressions."""
        return set().union(*(child.atoms() for child in self.children)) if self.children else set()

    def evaluate(self, fact_context: FactContext) -> bool:
        """Return whether any child is true in a fact context."""
        for child in sorted(self.children, key=lambda child: child.cost()):
            if child.evaluate(fact_context):
                return True
        return False

    def evaluate_partial(self, known: dict[str, bool]) -> bool | None:
        """Return the disjunction when the known values settle it."""
        unknown = False
        for child in self.children:
            value = child.evaluate_partial(known)
            if value is True:
                return True
            if value is None:
                unknown = True
        return None if unknown else False

    def __repr__(self) -> str:
        return "(" + " | ".join(repr(child) for child in self.children) + ")"


ExprLike = Expr | str


def _coerce(value: ExprLike) -> Expr:
    """Coerce a bare fact name into an Atom, passing an existing Expr through.

    Args:
        value (ExprLike): Either an already-built Expr, or a fact atom name.

    Returns:
        Expr: ``value`` itself, or ``Atom(value)`` if it was a string.
    """
    return Atom(value) if isinstance(value, str) else value


def all_of(*operands: ExprLike) -> Expr:
    """Build a conjunction from a mix of Expr nodes and bare fact atom names.

    Args:
        *operands (ExprLike): Expr nodes and/or fact atom names.

    Returns:
        Expr: An And of the coerced operands.
    """
    return And(*(_coerce(operand) for operand in operands))


def any_of(*operands: ExprLike) -> Expr:
    """Build a disjunction from a mix of Expr nodes and bare fact atom names.

    Args:
        *operands (ExprLike): Expr nodes and/or fact atom names.

    Returns:
        Expr: An Or of the coerced operands.
    """
    return Or(*(_coerce(operand) for operand in operands))


def no(operand: ExprLike) -> Expr:
    """Negate a single Expr node or bare fact atom name.

    Args:
        operand (ExprLike): An Expr node, or a fact atom name.

    Returns:
        Expr: A Not wrapping the coerced operand.
    """
    return Not(_coerce(operand))


def expr_from_json(payload: object) -> Expr:
    """Compile a JSON boolean expression payload (as used in the consequence catalogue) into an Expr.

    The payload grammar is:

    * a bare string -> Atom
    * ``{"not": <payload>}`` -> Not
    * ``{"all": [<payload>, ...]}`` -> And
    * ``{"any": [<payload>, ...]}`` -> Or

    Args:
        payload (object): A JSON-decoded boolean expression, in the grammar above.

    Returns:
        Expr: The compiled expression tree.

    Raises:
        ValueError: If ``payload`` does not match the grammar above.
    """
    if isinstance(payload, str):
        return Atom(payload)
    if isinstance(payload, dict):
        if "not" in payload and len(payload) == 1:
            return no(expr_from_json(payload["not"]))
        if "all" in payload and len(payload) == 1:
            return all_of(*(expr_from_json(item) for item in payload["all"]))
        if "any" in payload and len(payload) == 1:
            return any_of(*(expr_from_json(item) for item in payload["any"]))
    raise ValueError(f"Cannot compile boolean expression from payload: {payload!r}")
