class CachedAttributeGuardMixin:
    __slots__ = ()
    
    def __setattr__(self, name, value):
        """Override __setattr__ to guard against setting attributes that are cached properties or invalidate cached properties when their dependencies change.

        TBD: these can be very expansive as it will be called on every attribute set, so we may want to 
        optimize the checks or use a different approach if performance becomes an issue. Current implemenatation
        is required because @cached_property attributes cannot be used by attrs.on_setattr. Can we instead 
        re-implement the caching mechanism using attrs features and avoid the need for @cached_property?  
        
        Args:            
            name (str): The name of the attribute being set.
            value: The value to set the attribute to.

        Raises:            
            AttributeError: If attempting to set a cached property directly. 
        """
        if name not in type(self)._DERIVED and name not in type(self)._INVALIDATES:
            return super().__setattr__(name, value)

        if name in type(self)._DERIVED and name in self.__slots__:
            raise AttributeError(f"'{name}' is a cached property and cannot be reassigned")
        
        for invalidation_candidate in type(self)._INVALIDATES.get(name, ()):
            try:
                delattr(self, invalidation_candidate)
            except AttributeError:
                pass

        return super().__setattr__(name, value)