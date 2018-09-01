from nose.tools import assert_equal
from biotree.pressure import Pressure


def test_read():
    test = Pressure()
    test.path = "./data/"
    test.read("perfusion_test.xlsx")
    assert_equal(test.data.shape, (13689, 8))
    assert_equal(max(test.data.index), 1368.8)
    assert_equal(list(test.columns),
                 ['BT1880_HL', 'BT1881_HL', 'BT1882_HL', 'BT1883_HL',
                  'BT1884_HL', 'BT1885_HL', 'BT1886_HL', 'BT1887_HL'])
