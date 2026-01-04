import numpy as np
from functools import lru_cache, wraps

# https://stackoverflow.com/questions/52331944/cache-decorator-for-numpy-arrays
def np_cache(function):
    """Cache function results so that we don't need to rerun it if called with the same inputs.
    Numpy arrays are not hashable by themselves, so this is a workaround using tuple casting"""
    @lru_cache
    def cached_wrapper(*args, **kwargs):

        args = [np.array(a) if isinstance(a, tuple) and not isinstance(a, (int, float)) else a for a in args]
        kwargs = {
            k: np.array(v) if isinstance(v, tuple) and not isinstance(v, (int, float)) else v for k, v in kwargs.items()
        }

        # call the function now, with the array arguments
        return function(*args, **kwargs)

    # wrapper to convert array inputs to hashables, for caching
    @wraps(function)
    def wrapper(*args, **kwargs):
        args = [tuple(tuple(row) for row in a) if isinstance(a, np.ndarray) and a.ndim > 1 else tuple(a) for a in args]
        kwargs = {
            k: tuple(tuple(row) for row in v) if isinstance(v, np.ndarray) and v.ndim > 1 else tuple(v) for k, v in kwargs.items()
        }
        return cached_wrapper(*args, **kwargs)

    # copy lru_cache attributes over too
    wrapper.cache_info = cached_wrapper.cache_info
    wrapper.cache_clear = cached_wrapper.cache_clear

    return wrapper


def np_cache_minimal(function):
    """Same as np_cache, but for a function with only 1D array inputs or int/float inputs."""
    @lru_cache
    def cached_wrapper(*args, **kwargs):

        args = [np.array(a) if isinstance(a, tuple) else a for a in args]
        kwargs = {
            k: np.array(v) if isinstance(v, tuple) else v for k, v in kwargs.items()
        }

        # call the function now, with the array arguments
        return function(*args, **kwargs)

    # wrapper to convert array inputs to hashables, for caching
    @wraps(function)
    def wrapper(*args, **kwargs):
        args = [tuple(a) if isinstance(a, np.ndarray) else a for a in args]
        kwargs = {
            k: tuple(v) if isinstance(v, np.ndarray) else v for k, v in kwargs.items()
        }
        return cached_wrapper(*args, **kwargs)

    # copy lru_cache attributes over too
    wrapper.cache_info = cached_wrapper.cache_info
    wrapper.cache_clear = cached_wrapper.cache_clear

    return wrapper