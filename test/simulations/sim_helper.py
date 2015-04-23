"""
Random helper functions for simulations
"""
import os

def get_path_str(path, pardir):
    """
    If absolute path, then returns the path as is.
    If relative path, then returns absolute path of concatenated pardir/path
    :param str path:  absolute or relative file or directory path
    :param str pardir: parent directory to concatenate to path if path is relative directory
    :return str: absolute resolved path
    """
    if not os.path.isabs(path):
        return os.path.join(pardir, path)
    else:
        return path
