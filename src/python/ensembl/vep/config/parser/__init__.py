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
"""Readers for each configuration source."""

from ensembl.vep.config.parser.cli import (
    add_section_arguments,
    build_parser,
    check_no_flag_collisions,
    parse_args,
)
from ensembl.vep.config.parser.common import flag_for, resolve_option_key
from ensembl.vep.config.parser.config import read_config_file, read_legacy_config_file
from ensembl.vep.config.parser.env import default_env_name, read_environment

__all__ = [
    "add_section_arguments",
    "build_parser",
    "check_no_flag_collisions",
    "default_env_name",
    "flag_for",
    "parse_args",
    "read_config_file",
    "read_environment",
    "read_legacy_config_file",
    "resolve_option_key",
]
