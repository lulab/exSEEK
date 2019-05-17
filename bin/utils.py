from collections import namedtuple
import numpy as np
import inspect
GFFRecord = namedtuple('GFFRecord', ['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attr'])

def read_gff(fin, zero_based=False):
    '''Read a GFF file
    Args:
        fin: file-like object
        zero_based: subtract start position by 1 if True
    Yields:
        feature: GFFRecord object
    '''
    for lineno, line in enumerate(fin):
        if line.startswith('#'):
            continue
        c = line.strip().split('\t')
        if len(c) < 9:
            raise ValueError('less than 9 columns found in GFF file at line{}'.format(lineno + 1))
        c[3] = int(c[3])
        c[4] = int(c[4])
        if zero_based:
            c[3] -= 1
        attrs = {}
        for a in c[8].split(';'):
            i = a.find('=')
            key = a[:i]
            val = a[(i + 1):]
            attrs[key] = val
        c[8] = attrs
        yield GFFRecord._make(c)

def search_dict(data, keys):
    '''Intersect given keys with a dict and return the subset

    Parameters:
        data: dict
        keys: list of keys to search
    Returns:
        dict
    '''
    return {key:data[key] for key in keys if key in data}

def get_feature_importances(estimator):
    from sklearn.feature_selection import RFE
    from sklearn.neural_network import MLPClassifier
    from sklearn.svm import SVC

    '''Get feature importance attribute of an estimator
    '''
    if hasattr(estimator, 'coef_'):
        return np.ravel(np.abs(estimator.coef_))
    elif hasattr(estimator, 'feature_importances_'):
        return np.ravel(estimator.feature_importances_)
    elif isinstance(estimator, RFE):
        ranking = estimator.ranking_.astype('float') - 1
        return np.ravel(1.0 - ranking/ranking.max())
    elif isinstance(estimator, MLPClassifier):
        return np.zeros(estimator.coefs_[0].shape[0])
    elif isinstance(estimator, SVC):
        # SVC does not have feature importance
        return np.zeros(estimator.support_vectors_.shape[1])
    else:
        raise ValueError('the estimator should have either coef_ or feature_importances_ attribute')

def function_has_arg(func, arg):
    return arg in inspect.getargspec(func).args

def python_args_to_r_args(args):
    def to_r_arg(s):
        if isinstance(arg, str):
            return "'" + arg + "'"
        elif isinstance(arg, bool):
            return 'TRUE' if arg else 'FALSE'
        elif isinstance(arg, int):
            return str(arg)
        elif isinstance(arg, float):
            return str(arg)
        else:
            raise TypeError('{}: unsupported type: {}'.format(arg, type(arg)))

    r_args = {}
    for name, arg in args.items():
        is_atom_type = False
        for atom_type in (str, int, bool, float):
            if isinstance(arg, atom_type):
                r_args[name] = to_r_arg(arg)
                is_atom_type = True
                break
        if not is_atom_type:
            if isinstance(arg, tuple) or isinstance(arg, list):
                r_args[name] = 'c(' + ','.join(map(to_r_arg, arg)) + ')'
            else:
                raise TypeError('{}: unsupported type: {}'.format(arg, type(arg)))
    return ', '.join(('%s=%s'%(k, v) for k, v in r_args.items()))