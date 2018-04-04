
from collections import namedtuple
from .utils import safe_inverse, assert_kwargs_not_none


# Create path namedtuple object
mctable_names = ['eff',
                 'eff_err',
                 'pec',
                 'pec_err',
                 'NCmc',
                 'num_effects',
                 'num_causes',
                 ]

MCTables = namedtuple('MCTables', mctable_names)


@assert_kwargs_not_none('efficiencies',
                        'efficiencies_err',
                        'response',
                        'response_err')
def make_mctables(efficiencies=None, efficiencies_err=None, response=None,
                  response_err=None):

    efficiencies_err_inv = safe_inverse(efficiencies_err)
    NCmc = (efficiencies * efficiencies_err_inv)**2

    num_effects = response.shape[0]
    num_causes = response.shape[1]

    mctables = MCTables(eff=efficiencies,
                        eff_err=efficiencies_err,
                        pec=response,
                        pec_err=response_err,
                        NCmc=NCmc,
                        num_causes=num_causes,
                        num_effects=num_effects
                        )
    return mctables
