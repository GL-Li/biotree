from nose.tools import assert_equal
from biotree.sensor import Sensor


def test_instance():
    test = Sensor("this is a test.")
    assert_equal(test.description, "this is a test.")


def test_read():
    test = Sensor("test read")
    test.path = "./data/"
    test.read("foo.csv")
    test.read("bar.csv")
    assert_equal(list(test.data.columns), ["foo", "bar"])
    assert_equal(test.data.shape, (9001, 2))
    assert_equal(max(test.index), 900.0)
