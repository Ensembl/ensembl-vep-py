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

class BaseClassInstantiationError(TypeError):
    """Raised when attempting to instantiate a base class directly.
    """
    def __init__(self, msg='''Trying to instantiate a base class
                 Base classes in this codebase intentionally does not use __slot__ and should never 
    be instantiated — only subclassed. Concrete subclasses of such base classes must define 
    __slots__ to ensure memory efficiency.
                 ''', *args, **kwargs):
        super().__init__(msg, *args, **kwargs)