
import pytest

from pyunfold.unfold import iterative_unfold
from pyunfold.callbacks import Logger


@pytest.mark.parametrize('callbacks', [[Logger()], Logger()])
def test_logger(capsys, callbacks):

    # Run example case
    data = [100, 150]
    data_err = [10, 12.2]
    response = [[0.9, 0.1],
                [0.1, 0.9]]
    response_err = [[0.01, 0.01],
                    [0.01, 0.01]]
    efficiencies = [0.4, 0.67]
    efficiencies_err = [0.01, 0.01]
    # Perform iterative unfolding
    unfolded_results = iterative_unfold(data, data_err,
                                        response, response_err,
                                        efficiencies, efficiencies_err,
                                        return_iterations=True,
                                        callbacks=callbacks)

    # Get stdout and std err from iterative_unfold
    out, err = capsys.readouterr()

    # Build expected output
    expected_output = ''
    for row_index, row in unfolded_results.iterrows():
        row_output = ('Iteration {}: ts = {:0.4f}, ts_stopping ='
                      ' {}\n'.format(row_index + 1,
                                     row['ts_iter'],
                                     row['ts_stopping']))
        expected_output += row_output

    assert expected_output == out
