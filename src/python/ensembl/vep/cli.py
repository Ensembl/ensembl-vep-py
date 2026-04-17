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
"""Hello world module"""

__all__ = ["hello"]

import argparse
from .ensembl_vep_core import sum_as_string  # noqa: F401,F403 pylint: disable=E0401,E0611


def hello(word: str) -> str:
    """Returns the concatenation of "Hello " with the given string.

    Args:
        word: The string you want to concatenate to "Hello ".

    """
    return sum_as_string("Hello ", word)


def main() -> None:
    """Module's entry-point."""

    parser = argparse.ArgumentParser(prog="Greetings program", description="Says 'hi' to the argument given")
    parser.add_argument("word")
    args = parser.parse_args()

    print(hello(args.word))
