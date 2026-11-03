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
"""Registered corpus export."""

from __future__ import annotations

import logging
from typing import ClassVar

from ensembl.vep.config import Config
from ensembl.vep.corpus.exporter import get_exporter
from ensembl.vep.corpus.registry import CorpusOptionsError, build_corpus_from_options
from ensembl.vep.runner.base import ExitCode
from ensembl.vep.runner.registry import register_runner

__all__ = ["Runner"]

logger = logging.getLogger(__name__)


@register_runner("corpus.save")
class Runner:
    """Runner for exporting a corpus as FASTA, GFF3 and VCF.

    Attributes:
        name (ClassVar[str]): The registry name.
        sections (ClassVar[tuple[str, ...]]): The configuration sections read.
        description (ClassVar[str]): One-line description for the subcommand list.
    """

    name: ClassVar[str] = "corpus.save"
    sections: ClassVar[tuple[str, ...]] = ("logging", "corpus", "corpus_save")
    description: ClassVar[str] = "write a registered corpus as FASTA, GFF3 and VCF"

    def __init__(self, config: Config) -> None:
        """Hold the configuration."""
        self.config = config

    def run(self) -> ExitCode:
        """Write the configured corpus and report what was written."""
        corpus = self.config.corpus
        if corpus.corpus_dir is None:
            logger.error("give --corpus_dir DIR to say where the corpus should be written")
            return ExitCode.USAGE

        try:
            built_corpus = build_corpus_from_options(corpus)
        except CorpusOptionsError as unusable:
            logger.error("%s", unusable)
            return ExitCode.USAGE
        try:
            manifest = get_exporter(self.config.corpus_save.exporter).export(
                built_corpus,
                corpus.corpus_dir,
                force=self.config.corpus_save.force,
            )
        except KeyError as unknown_exporter:
            logger.error("%s", unknown_exporter.args[0])
            return ExitCode.USAGE
        except FileExistsError as present:
            logger.error("%s", present)
            return ExitCode.FAILURE
        finally:
            built_corpus.close()

        logger.info(
            "wrote %s: %d transcripts, %d variants, %d overlaps on %s (identity %s)",
            corpus.corpus_dir,
            manifest.transcripts,
            manifest.variants,
            manifest.overlaps,
            ", ".join(manifest.chromosomes),
            manifest.identity,
        )
        return ExitCode.OK
