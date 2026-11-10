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
"""Worker process count must not change what a run reports."""

from __future__ import annotations

from pathlib import Path

import ensembl.vep.annotators  # noqa: F401  pylint: disable=unused-import

# ^ registers the strategies the default annotator resolves.
from ensembl.vep.common.enums import StartMethod
from ensembl.vep.config import Config, ExecutionConfig, RunConfig, SourceConfig
from ensembl.vep.runner.base import ExitCode
from ensembl.vep.runner.run.default import Runner

GFF3 = """\
chr1\tensembl\tgene\t1000\t5000\t.\t+\t.\tID=gene:G1;Name=GENE1;biotype=protein_coding
chr1\tensembl\tmRNA\t1000\t5000\t.\t+\t.\tID=transcript:ENST001;Parent=gene:G1;biotype=protein_coding
chr1\tensembl\texon\t1000\t2000\t.\t+\t.\tID=exon:E1;Parent=transcript:ENST001
chr1\tensembl\tCDS\t1100\t2000\t.\t+\t0\tID=cds:P1;Parent=transcript:ENST001;protein_id=P1
chr1\tensembl\texon\t3000\t5000\t.\t+\t.\tID=exon:E2;Parent=transcript:ENST001
chr1\tensembl\tCDS\t3000\t4900\t.\t+\t0\tID=cds:P1;Parent=transcript:ENST001;protein_id=P1
chr2\tensembl\tgene\t1000\t5000\t.\t+\t.\tID=gene:G2;Name=GENE2;biotype=protein_coding
chr2\tensembl\tmRNA\t1000\t5000\t.\t+\t.\tID=transcript:ENST002;Parent=gene:G2;biotype=protein_coding
chr2\tensembl\texon\t1000\t2000\t.\t+\t.\tID=exon:E3;Parent=transcript:ENST002
chr2\tensembl\tCDS\t1100\t2000\t.\t+\t0\tID=cds:P2;Parent=transcript:ENST002;protein_id=P2
chr2\tensembl\texon\t3000\t5000\t.\t+\t.\tID=exon:E4;Parent=transcript:ENST002
chr2\tensembl\tCDS\t3000\t4900\t.\t+\t0\tID=cds:P2;Parent=transcript:ENST002;protein_id=P2
"""
"""Two contigs, each with one coding transcript spanning two exons."""

VCF = """\
##fileformat=VCFv4.3
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
chr1\t1500\trs_c1_coding\tA\tC\t.\tPASS\t.
chr1\t1700\trs_c1_multi\tA\tC,G\t.\tPASS\t.
chr1\t1200\trs_c1_indel\tC\tCAT\t.\tPASS\t.
chr1\t30000\trs_c1_intergenic\tA\tC\t.\tPASS\t.
chr2\t1500\trs_c2_coding\tA\tC\t.\tPASS\t.
chr2\t1700\trs_c2_multi\tA\tC,G\t.\tPASS\t.
chr2\t1200\trs_c2_indel\tC\tCAT\t.\tPASS\t.
"""
"""Four records on chr1 and three on chr2, so a ``buffer_size=3`` batch boundary falls both
mid-contig and exactly at the contig change: ``iter_vcf_buffers`` never mixes regions in one
batch, so this yields three buffers (3, 1, 3) split across two worker processes."""


def _fasta_and_fai(sequences: dict[str, str]) -> tuple[str, str]:
    """Return single-line-per-contig FASTA text and its matching ``.fai`` index text."""
    fasta_parts: list[str] = []
    fai_lines: list[str] = []
    offset = 0
    for name, sequence in sequences.items():
        header = f">{name}\n"
        fasta_parts.append(header)
        offset += len(header)
        fai_lines.append(f"{name}\t{len(sequence)}\t{offset}\t{len(sequence)}\t{len(sequence) + 1}\n")
        fasta_parts.append(sequence + "\n")
        offset += len(sequence) + 1
    return "".join(fasta_parts), "".join(fai_lines)


def _inputs(directory: Path) -> tuple[Path, Path, Path]:
    """Write the multi-contig VCF, GFF3 and FASTA fixture, returning their paths."""
    variants = directory / "variants.vcf"
    variants.write_text(VCF, encoding="utf-8")
    transcripts = directory / "transcripts.gff3"
    transcripts.write_text(GFF3, encoding="utf-8")
    fasta = directory / "ref.fa"
    fasta_text, fai_text = _fasta_and_fai({"chr1": "ACGT" * 10000, "chr2": "ACGT" * 10000})
    fasta.write_text(fasta_text, encoding="utf-8")
    fasta.with_suffix(".fa.fai").write_text(fai_text, encoding="utf-8")
    return variants, transcripts, fasta


def _data_lines(output: Path) -> list[str]:
    """The output's records, without its headers."""
    return [line for line in output.read_text(encoding="utf-8").splitlines() if not line.startswith("#")]


def _csq_sets(lines: list[str]) -> dict[str, frozenset[str]]:
    """Return each variant's CSQ entries as an order-independent set, keyed by its ID."""
    sets: dict[str, frozenset[str]] = {}
    for line in lines:
        columns = line.split("\t")
        variant_id = columns[2]
        csq = ""
        for field in columns[7].split(";"):
            if field.startswith("CSQ="):
                csq = field[len("CSQ=") :]
                break
        sets[variant_id] = frozenset(csq.split(",")) if csq else frozenset()
    return sets


class TestParallelismEquivalence:
    """Two worker processes must annotate exactly what one process would."""

    def test_two_processes_match_one_process_at_the_same_buffer_size(self, tmp_path: Path) -> None:
        """Order, per-variant CSQ sets, and full output bytes all agree across process counts.

        CSQ *entry* order within one variant may legitimately depend on transcript insertion order
        into the per-payload ``BinnedFeatureCache``, which a different ``buffer_size`` can disturb;
        comparing sets rather than tuples accepts that. Both runs here share ``buffer_size=3``, so
        that source of difference is absent too, and the two outputs are checked byte-identical.
        """
        # Inputs are written once and shared: the writer embeds the ``--input_file`` path in a
        # VCF header comment, and giving the two runs different input paths would make that
        # comment differ even though the annotation itself is identical.
        variants, transcripts, fasta = _inputs(tmp_path)

        serial_output = self._run(tmp_path / "serial.vcf", variants, transcripts, fasta, processes=1)
        parallel_output = self._run(tmp_path / "parallel.vcf", variants, transcripts, fasta, processes=2)

        serial_lines = _data_lines(serial_output)
        parallel_lines = _data_lines(parallel_output)
        assert [line.split("\t")[2] for line in parallel_lines] == [
            line.split("\t")[2] for line in serial_lines
        ]
        assert _csq_sets(parallel_lines) == _csq_sets(serial_lines)
        assert parallel_output.read_text(encoding="utf-8") == serial_output.read_text(encoding="utf-8")

    @staticmethod
    def _run(output: Path, variants: Path, transcripts: Path, fasta: Path, *, processes: int) -> Path:
        """Annotate the fixture through the real runner at a given process count."""
        config = Config(
            source=SourceConfig(input_file=variants, gff=transcripts, fasta=fasta),
            execution=ExecutionConfig(
                processes=processes,
                start_method=StartMethod.SPAWN,
                buffer_size=3,
            ),
            run=RunConfig(output_file=output, no_progress=True, no_stats=True),
        )
        assert Runner(config).run() == ExitCode.OK
        return output
