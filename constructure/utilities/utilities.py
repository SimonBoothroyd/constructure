import functools
import importlib

_CONDA_INSTALLATION_COMMANDS = {
    "openeye": "conda install -c openeye openeye-toolkits",
    "rdkit": "conda install -c conda-forge rdkit",
}


class MissingOptionalDependency(ValueError):
    """An exception raised when an optional dependency is required
    but cannot be found.

    Attributes:
        library_name: The name of the missing library.
        license_issue: Whether the library was importable but was unusable due to a
            missing license.
    """

    def __init__(self, library_name: str, license_issue: bool = False):
        """
        Args:
            library_name: The name of the missing library.
            license_issue: Whether the library was importable but was unusable due to a
                missing license.
        """

        message = f"The required {library_name} module could not be imported."
        conda_command = _CONDA_INSTALLATION_COMMANDS.get(
            library_name.split(".")[0], None
        )

        if license_issue:
            message = f"{message} This is due to a missing license."
        elif conda_command is not None:
            message = (
                f"{message} Try installing the package by running `{conda_command}`."
            )

        super(MissingOptionalDependency, self).__init__(message)

        self.library_name = library_name
        self.license_issue = license_issue


def requires_package(import_path: str):
    """A decorator which checks that the required package / module is installed.

    Args:
        import_path: The import path of the required package / module.
    """

    def inner_decorator(function):
        @functools.wraps(function)
        def wrapper(*args, **kwargs):

            try:
                importlib.import_module(import_path)
            except (ImportError, ModuleNotFoundError):
                raise MissingOptionalDependency(import_path, False)

            return function(*args, **kwargs)

        return wrapper

    return inner_decorator
