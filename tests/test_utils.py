import numpy as np
from ts_cluster.zminer.utils import load_file, preprocess, get_interval_duration, hash_interval, compare_intervals, remove_intervals_from_database, make_constraints

def test_load_file() -> None:
    # CASE 1
    filepath_1 = "src/ts_cluster/zminer/data/ASL_BU_1.txt"
    loaded_file_1 = load_file(filepath_1)
    assert loaded_file_1.shape == (14802, 4)

    # CASE 2
    filepath_2 = "src/ts_cluster/zminer/data/ASL_BU_2.txt"
    loaded_file_2 = load_file(filepath_2)
    assert loaded_file_2.shape == (41761, 4)


def test_preprocess() -> None:
    # CASE 1 - SIMPLE CASE
    database_1 = np.array([[0, 12, 1, 300], [0, 12, 2, 444],  [1, 23, 56, 99],  [1, 54, 22, 100]])
    tseq, tdis, tintv, aintv, avgtime, eseqdb, unique_labels, initial_supports = preprocess.py_func(database_1)

    assert tseq == 2 #0, 1
    assert tdis == 3 #12, 23, 54
    assert tintv == 4 #four intervals
    assert aintv == 2.0 #average interval size per e-seq
    assert avgtime == 272.0
    assert eseqdb[0][0] == (12, 1, 300)
    assert unique_labels[0] == {12}
    assert initial_supports[12] == 1

    # CASE 2
    filepath_2 = "src/ts_cluster/zminer/data/ASL_BU_2.txt"
    loaded_file_2 = load_file(filepath_2)
    tseq, tdis, tintv, aintv, avgtime, eseqdb, unique_labels, initial_supports = preprocess.py_func(loaded_file_2)
    assert tseq == 1839
    assert tdis == 254
    assert tintv == 41761
    assert np.round(aintv, 2) == 22.71 #average interval size per e-seq
    assert np.round(avgtime, 2) == 3067.50
    assert eseqdb[0][0] == (126, 1, 3334)
    assert unique_labels[0] == {1,12,14,20,21,32,66,93,95,126,128,194,196,202,204,205,206,208,209,211}


def test_get_interval_duration() -> None:
    # CASE 1
    assert get_interval_duration.py_func((0, 1, 10)) == 9
    
    # CASE 2
    assert get_interval_duration.py_func((0, 20, 100)) == 80

def test_hash_interval() -> None:
    # CASE 1
    assert hash_interval.py_func((0, 1, 10)) == hash((0, 1, 10))

    # CASE 2
    assert hash_interval.py_func((10, 20, 100)) == hash((10, 20, 100))
    
def test_compare_intervals() -> None:
    assert compare_intervals.py_func((10, 100, 120), (10, 100, 140)) == True

    assert compare_intervals.py_func((10, 100, 120), (10, 60, 140)) == False

def test_remove_intervals_from_database() -> None:
    test_database = np.array([[0, 126, 1, 300], [0, 122, 2, 444],  [1, 23, 56, 99],  [1, 54, 22, 100]])
    _, _, _, _, _, eseqdb, _, initial_support = preprocess(test_database)

    assert len(remove_intervals_from_database.py_func(eseqdb, initial_support, 1)) == 2
    assert len(remove_intervals_from_database.py_func(eseqdb, initial_support, 1)[0]) == 2

    assert len(remove_intervals_from_database.py_func(eseqdb, initial_support, 2)) == 2
    assert len(remove_intervals_from_database.py_func(eseqdb, initial_support, 2)[0]) == 0


def test_make_constraints() -> None:
    # CASE 1
    test_database = np.array([[0, 126, 1, 300], [0, 122, 2, 444],  [1, 23, 56, 99],  [1, 54, 22, 100]])
    np.testing.assert_array_equal(make_constraints.py_func(np.array([0,0,0,0,0,0,0,0]), test_database), [0., 0., 0., 0., 0., 0., 0., 0.])
    
    # CASE 2
    # 30% minsUp, epsilon = 1, gap = inf, timeout = inf, level = 20, print = True
    constraints = make_constraints.py_func(np.array([0.3,1,-1,-1,20,1,1]), test_database)
    np.testing.assert_array_equal(constraints, [0.3, 1.2, 1., np.inf, np.inf, 20., 1., 1.])