
def lazy(f):
    def wrapper(*args, **kwargs):
        if hasattr(wrapper, 'val'):
            return wrapper.val
        res = f(*args, **kwargs)
        wrapper.val = res
        return res
    return wrapper
