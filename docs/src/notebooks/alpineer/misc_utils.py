import os
import warnings
from collections.abc import Iterable
from typing import Any

import matplotlib.pyplot as plt
import numpy as np

from alpineer import io_utils


def save_figure(save_dir, save_file, dpi=None):
    """Verify save_dir and save_file, then save to specified location

    Args:
        save_dir (str):
            the name of the directory we wish to save to
        save_file (str):
            the name of the file we wish to save to
        dpi (float):
            the resolution of the figure
    """

    # path validation
    io_utils.validate_paths(save_dir)

    # verify that if save_dir specified, save_file must also be specified
    if save_file is None:
        raise FileNotFoundError("save_dir specified but no save_file specified")

    plt.savefig(os.path.join(save_dir, save_file), dpi=dpi)


def create_invalid_data_str(invalid_data):
    """Creates a easy to read string for ValueError statements.

    Args:
        invalid_data (list[str]): A list of strings containing the invalid / missing data

    Returns:
        str: Returns a formatted string for more detailed ValueError outputs.
    """
    # Holder for the error string
    err_str_data = ""

    # Adding up to 10 invalid values to the err_str_data.
    for idx, data in enumerate(invalid_data[:10], start=1):
        err_msg = "{idx:{fill}{align}{width}} {message}\n".format(
            idx=idx,
            message=data,
            fill=" ",
            align="<",
            width=12,
        )
        err_str_data += err_msg

    return err_str_data


def make_iterable(a: Any, ignore_str: bool = True):
    """Convert noniterable type to singleton in list

    Args:
        a (T | Iterable[T]):
            value or iterable of type T
        ignore_str (bool):
            whether to ignore the iterability of the str type

    Returns:
        List[T]:
            a as singleton in list, or a if a was already iterable.
    """
    return (
        a
        if isinstance(a, Iterable)
        and not ((isinstance(a, str) and ignore_str) or isinstance(a, type))
        else [a]
    )


def verify_in_list(warn: bool = False, **kwargs) -> bool:
    """Verify at least whether the values in the first list exist in the second

    Args:
        warn (bool):
            Whether to issue warning instead of error, defaults to False
        **kwargs (list, list):
            Two lists, but will work for single elements as well.
            The first list specified will be tested to see
            if all its elements are contained in the second.

    Raises:
        ValueError:
            if not all values in the first list are found in the second
        Warning:
            if not all values are found and warn is True
    """

    if len(kwargs) != 2:
        raise ValueError("You must provide 2 arguments to verify_in_list")

    rhs_name, lhs_name = map(lambda s: s.replace("_", " "), kwargs.keys())

    lhs_list, rhs_list = map(lambda l: list(make_iterable(l)), kwargs.values())

    # Check if either list inputs are `None` (or both)
    if any(v == [None] for v in (lhs_list, rhs_list)):
        return True
    # If the Left Hand List is empty, the Right Hand List will always be contained in it.
    if len(lhs_list) == 0:
        return True
    # If the Right Hand List is empty, the Left Hand List will never be contained in it.
    # Throw a ValueError
    if len(rhs_list) == 0:
        raise ValueError(f"The list {lhs_name} is empty.")

    if not np.isin(lhs_list, rhs_list).all():
        # Calculate the difference between the `lhs_list` and the `rhs_list`
        difference = [str(val) for val in lhs_list if val not in rhs_list]

        # Only printing up to the first 10 invalid values.
        err_str = (
            "Not all values given in list {0:^} were found in list {1:^}.\n "
            "Displaying {2} of {3} invalid value(s) for list {4:^}\n"
        ).format(
            rhs_name,
            lhs_name,
            min(len(difference), 10),
            len(difference),
            rhs_name,
        )

        err_str += create_invalid_data_str(difference)

        if warn:
            warnings.warn(err_str)
        else:
            raise ValueError(err_str)
    return True


def verify_same_elements(enforce_order=False, warn=False, **kwargs) -> bool:
    """Verify if two lists contain the same elements regardless of count

    Args:
        enforce_order (bool):
            Whether to also check for the same ordering between the two lists
        warn (bool):
            Whether to issue warning instead of error, defaults to False
        **kwargs (list, list):
            Two lists

    Raises:
        ValueError:
            if the two lists don't contain the same elements
    """

    if len(kwargs) != 2:
        raise ValueError("You must provide 2 list arguments to verify_same_elements")

    list_one, list_two = kwargs.values()

    try:
        list_one_cast = list(list_one)
        list_two_cast = list(list_two)
    except TypeError:
        raise ValueError("Both arguments provided must be lists or list types")

    # Check if either list inputs are `None` (or both)
    if any(v == [None] for v in [list_one_cast, list_two_cast]):
        return True
    # If both lists are empty
    if all(len(v) == 0 for v in [list_one_cast, list_two_cast]):
        return True
    # If either list one, or two have a length of 0, but not both.
    if (len(list_one_cast) == 0) is not (len(list_two_cast) == 0):
        return False

    list_one_name, list_two_name = kwargs.keys()
    list_one_name = list_one_name.replace("_", " ")
    list_two_name = list_two_name.replace("_", " ")

    if not np.all(set(list_one_cast) == set(list_two_cast)):
        # Values in list one that are not in list two
        missing_vals_1 = [str(val) for val in (set(list_one_cast) - set(list_two_cast))]

        # Values in list two that are not in list one
        missing_vals_2 = [str(val) for val in (set(list_two_cast) - set(list_one_cast))]

        # Total missing values
        missing_vals_total = [str(val) for val in set(list_one_cast) ^ set(list_two_cast)]

        err_str = (
            "{0} value(s) provided for list {1:^} and list {2:^} are not found in both lists.\n"
        ).format(len(missing_vals_total), list_one_name, list_two_name)

        # Only printing up to the first 10 invalid values for list one.
        err_str += ("{0:>13} \n").format(
            "Displaying {0} of {1} value(s) in list {2} that are missing from list {3}\n".format(
                min(len(missing_vals_1), 10), len(missing_vals_1), list_one_name, list_two_name
            )
        )
        err_str += create_invalid_data_str(missing_vals_1) + "\n"

        # Only printing up to the first 10 invalid values for list two
        err_str += ("{0:>13} \n").format(
            "Displaying {0} of {1} value(s) in list {2} that are missing from list {3}\n".format(
                min(len(missing_vals_2), 10), len(missing_vals_2), list_two_name, list_one_name
            )
        )
        err_str += create_invalid_data_str(missing_vals_2) + "\n"

        if warn:
            warnings.warn(err_str)
        else:
            raise ValueError(err_str)
    elif enforce_order and list_one_cast != list_two_cast:
        first_bad_index = next(
            i for i, (l1, l2) in enumerate(zip(list_one_cast, list_two_cast)) if l1 != l2
        )

        err_str = "Lists %s and %s ordered differently: values %s and %s do not match at index %d"

        if warn:
            warnings.warn(
                err_str
                % (
                    list_one_name,
                    list_two_name,
                    list_one_cast[first_bad_index],
                    list_two_cast[first_bad_index],
                    first_bad_index,
                )
            )
        else:
            raise ValueError(
                err_str
                % (
                    list_one_name,
                    list_two_name,
                    list_one_cast[first_bad_index],
                    list_two_cast[first_bad_index],
                    first_bad_index,
                )
            )
    return True
