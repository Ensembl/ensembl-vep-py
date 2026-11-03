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
"""The worker contract: one payload in, one serialisable result out."""

from __future__ import annotations

from typing import Protocol, TypeVar, runtime_checkable

__all__ = ["Worker"]

PayloadT_contra = TypeVar("PayloadT_contra", contravariant=True)
ResultT_co = TypeVar("ResultT_co", covariant=True)


@runtime_checkable
class Worker(Protocol[PayloadT_contra, ResultT_co]):
    """Process-local runtime that runs one payload at a time.

    A worker is built once per process, after that process starts, and then reused: it holds open
    sources and warm caches that neither cross a process boundary nor survive a rebuild. So what
    the executor transports is the factory, never the worker, and the worker itself never has to
    be serialisable -- only its payloads and its results do.

    One method, because that is all an executor needs of it. What a payload holds, and what the
    result carries, belongs to the kind of work: ``AnnotationBufferWorker`` takes a payload of
    input records and answers with detached annotations.
    """

    def run(self, payload: PayloadT_contra) -> ResultT_co:
        """Return the result of this payload, against the sources this worker holds.

        Args:
            payload (PayloadT_contra): One unit of work, as its executor's payloads are grouped.

        Returns:
            ResultT_co: What the payload produced, serialisable enough to cross back.
        """
