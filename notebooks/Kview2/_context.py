import numpy as np

blacklist = ["try_all_threshold",
             "test",
             "minimum_of_touching_neighbors_map",
             "maximum_of_touching_neighbors_map",
             "standard_deviation_of_touching_neighbors_map",
             "mean_of_touching_neighbors_map",
             "z_position_of_minimum_z_projection",
             "z_position_of_maximum_z_projection",
             ]

def nop(image):
    """No operation, returns the image unchanged"""
    return image



class Context():
    def __init__(self, variables):
        # from ._utilities import logo

        self._functions = {"stackview.nop": nop}
        self._modules = { }
        # self._images = {"no_image": logo}
        #self._label_images = {}
        self._images = {} 

        self.parse(variables.items())

    def parse(self, items, prefix: str = None):

        from types import ModuleType
        from typing import Callable
        import stackview
        from ._utilities import is_image, is_label_image, count_image_parameters

        for name, value in items:
            if not name.startswith("_"):
                if isinstance(value, ModuleType) and value is not stackview:
                    if prefix is None:
                        self._modules[name] = value
                        self.parse({key: getattr(value, key) for key in dir(value)}.items(), prefix=name + ".")
                elif isinstance(value, Callable):
                    if count_image_parameters(value) > 0:
                        if name not in blacklist:
                            if prefix is None:
                                self._functions[name] = value
                            else:
                                self._functions[prefix + name] = value
                elif isinstance(value, dict) and ("marker_dict_cropped" in str(name)):
                    # If a dictionary of images is found, add them to the _images attribute.
                    if prefix is None:
                        self._images.update(value)
                    else:
                        self._images.update({f"{prefix}{k}": v for k, v in value.items()})
                elif is_image(value) and (("signal_sub1" == str(name)) or ("signal_sub2" == str(name)) or ("dn_signal" == str(name)) or ("signal_gauss" == str(name)) or ("signal_final" == str(name)) or ("signal_final_clahe" == str(name))):
                    if prefix is None:
                        self._images[name] = value
                    else:
                        self._images[prefix + name] = value


