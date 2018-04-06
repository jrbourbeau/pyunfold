
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
    """
    def __init__(self):
        super(Callback, self).__init__()
        pass

    def on_iteration_end(self, iteration, params):
        output = ('Iteration {}: ts = {:0.4f}, ts_stopping ='
                  ' {}\n'.format(iteration + 1,
                                 params['ts_iter'],
                                 params['ts_stopping']))
        sys.stdout.write(output)
