from nose.tools import assert_equal
from biotree.utilities import read_1567, linear_fit, model_1567
import numpy as np
import pandas as pd


def test_read_1567():
    BT1567 = read_1567()
    assert_equal(max(BT1567.index), 392.48)


def test_linear_fit():
    x = np.array([1, 2, 3])
    y = np.array([2, 4, 6])
    ser = pd.Series(data=y, index=x)
    fit = linear_fit(ser)
    assert_equal(fit, (0, 2))


def test_model_1567():
    model = model_1567()
    assert_equal(model[479], 678.16765451324272)
    assert_equal(model[40], 41.593525156541098)
    assert_equal(list(model[100]), [41.5935251565411, -86.98916604412273])
