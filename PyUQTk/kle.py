# Import the low-level C/C++ module
if __package__ or "." in __name__:
    from . import _kle
else:
    import _kle
