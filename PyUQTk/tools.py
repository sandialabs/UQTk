# Import the low-level C/C++ module
if __package__ or "." in __name__:
    from . import _tools
else:
    import _tools
