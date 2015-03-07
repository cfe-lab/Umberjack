from argparse import ArgumentParser

class ConfigArgParser(ArgumentParser):
    """
    Override argparse.ArgumentParser to converts python config files to commandline arguments.
    """
    def convert_arg_line_to_args(self, arg_line):
        """
        Parse config file.
        :param str arg_line: a line in the file used for ArgumentParser arguments
        :return iterator : an iterator for each commandline argument equivalent to the config file
        """
        for arg in arg_line.split():
            if not arg.strip():
                continue
            yield arg
