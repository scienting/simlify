# This file is licensed under the Prosperity Public License 3.0.0.
# You may use, copy, and share it for noncommercial purposes.
# Commercial use is allowed for a 30-day trial only.
#
# Contributor: Scienting Studio
# Source Code: https://github.com/scienting/simlify
#
# See the LICENSE.md file for full license terms.

from typing import Any

import importlib
from collections.abc import Iterable

from loguru import logger


def get_obj_from_string(import_string: str) -> object:
    """Retrieves an object based on an import string.

    The import string should specify the full path to the desired object,
    starting from the root module. This function dynamically imports the
    module and retrieves the specified object (which can be a class,
    function, or any other Python object).

    Args:
        import_string: The import string, which is a dot-separated string
            representing the module path and the object name. For example:
            `"os.path.join"` or `"my_module.MyClass"`.

    Returns:
        The object identified by the import string.

    Raises:
        ImportError: If the module specified in the import string cannot be found.
        AttributeError: If the object specified in the import string does not
            exist within the imported module.

    Examples:
        ```python
        import os

        join_func = get_obj_from_string("os.path.join")
        assert join_func is os.path.join
        ```
    """
    logger.debug("Importing {}", import_string)
    module_name, obj_name = import_string.rsplit(".", 1)
    module = importlib.import_module(module_name)
    obj = getattr(module, obj_name)
    return obj


def simple_generator(iterable: Iterable[Any]) -> Any:
    """Iterates through an iterable and yields each item.

    Args:
        iterable: The iterable to yield elements from.

    Yields:
        An element from the input `iterable`.

    Examples:
        ```python
        list(simple_generator([1, 2, 3]))
        # Expected output: [1, 2, 3]
        ```

        ```python
        list(simple_generator("abc"))
        # Expected output: ['a', 'b', 'c']
        ```
    """
    for i in iterable:
        yield i
