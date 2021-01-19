import pytest

from constructure.utilities import MissingOptionalDependency, requires_package
from constructure.utilities.utilities import _CONDA_INSTALLATION_COMMANDS


def test_requires_package_found():
    @requires_package("constructure")
    def dummy_function():
        return 42

    assert dummy_function() == 42


def test_requires_package_unknown_missing():
    @requires_package("fake-package-42")
    def dummy_function():
        pass

    with pytest.raises(MissingOptionalDependency) as error_info:
        dummy_function()

    assert "The required fake-package-42 module could not be imported." in str(
        error_info.value
    )


def test_requires_package_known_missing(monkeypatch):

    monkeypatch.setitem(
        _CONDA_INSTALLATION_COMMANDS, "fake-package-42", "conda install ..."
    )

    @requires_package("fake-package-42")
    def dummy_function():
        pass

    with pytest.raises(MissingOptionalDependency) as error_info:
        dummy_function()

    assert "Try installing the package by running `conda install ...`" in str(
        error_info.value
    )
