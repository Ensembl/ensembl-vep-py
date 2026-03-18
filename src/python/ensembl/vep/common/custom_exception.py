class BaseClassInstantiationError(TypeError):
    """Raised when attempting to instantiate a base class directly.
    """
    def __init__(self, msg='''Trying to instantiate a base class
                 Base classes in this codebase intentionally does not use __slot__ and should never 
    be instantiated — only subclassed. Concrete subclasses of such base classes must define 
    __slots__ to ensure memory efficiency.
                 ''', *args, **kwargs):
        super().__init__(msg, *args, **kwargs)