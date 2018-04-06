
import sys


class Callback(object):
    """Callback base class
    """
    def __init__(self):
        pass

    def on_unfolding_begin(self):
        pass

    def on_unfolding_end(self):
        pass

    def on_iteration_begin(self, iteration, params):
        pass

    def on_iteration_end(self, iteration, params):
        pass


class Logger(Callback):
    """Logger callback

    Writes test statistic information for each iteration to sys.stdout.
    """
    def __init__(self):
        super(Callback, self).__init__()
        pass

    def on_iteration_end(self, iteration, params):
        """Writes to sys.stdout

        Parameters
        ----------
        iteration : int
            Unfolding iteration (i.e. iteration=1 after first unfolding
            iteration, etc.).
        params : dict
            Dictionary containing key value pairs for the current test
            statistic value (``'ts_iter'``) and the final test statistic stopping
            condition (``'ts_stopping'``).
        """
        output = ('Iteration {}: ts = {:0.4f}, ts_stopping ='
                  ' {}\n'.format(iteration,
                                 params['ts_iter'],
                                 params['ts_stopping']))
        sys.stdout.write(output)
