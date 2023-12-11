import importlib

from loguru import logger


def get_obj_from_string(import_string: str) -> object:
    """Retrieves an object based on an import string and object name.

    Args:
        import_string: The import string, starting from the root module, containing
            the desired object.

    Returns:
        The object identified by the import string.
    """
    logger.debug("Importing {}", import_string)
    module_name, obj_name = import_string.rsplit(".", 1)
    module = importlib.import_module(module_name)
    obj = getattr(module, obj_name)
    return obj
