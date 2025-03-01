
def interact(func,
             image = None,
             *args,
             continuous_update: bool = True,
             context:dict = None,
             zoom_factor:float = 1.0,
             zoom_spline_order:int = 0,
             colormap:str = None,
             display_min:float = None,
             display_max:float = None,
             min_value:int = 0, 
             max_value:int = 1000,
             min_value_float:float = 0.1, 
             max_value_float:float = 10.0, 
             step_float:float = 0.1, 
             step:int = 100,
             viewer = None,
             **kwargs):
    
    """Takes a function which has an image as first parameter and additional parameters.
    It will build a user interface consisting of sliders for numeric parameters and parameters
    that are called "footprint" or "selem".

    Parameters
    ----------
    func : function
    image : Image, optional
        If not provided, context must be provided instead.
    args
    continuous_update : bool, optional
        Update the image while dragging the mouse, default: False
    context:dict, optional
        A dictionary of (name:image), allows showing a pulldown of available images, e.g.: globals()
    zoom_factor: float, optional
        Allows showing the image larger (> 1) or smaller (<1)
    zoom_spline_order: int, optional
        Spline order used for interpolation (default=0, nearest-neighbor)
    colormap: str, optional
        Matplotlib colormap name or "pure_green", "pure_magenta", ...
    display_min: float, optional
        Lower bound of properly shown intensities
    display_max: float, optional
        Upper bound of properly shown intensities
    viewer: _SliceViewer, optional
        The viewer where the result image should be shown.
    kwargs

    """
    import inspect
    import ipywidgets
    from ._context import Context
    from skimage import morphology
    from ipywidgets import GridspecLayout

    if 'cupy.ndarray' in str(type(image)):
        image = image.get()

    if isinstance(func, tuple):
        func_name = func[1]
        func = func[0]
    else:
        func_name = func.__name__

    exposable_parameters = []
    footprint_parameters = []
    image_parameters = []

    if context is not None and not isinstance(context, Context):
        context = Context(context)

    image_passed = image is not None

    if context is not None and image is None:
        image = next(iter(context._images.values()))

    sig = inspect.signature(func)

    for key in sig.parameters.keys():
        exposable = False
        default_value = 0
        if isinstance(sig.parameters[key].default, int) or isinstance(sig.parameters[key].default, str):
            default_value = sig.parameters[key].default
            min_value, max_value, step = min_value, max_value, step

        elif isinstance(sig.parameters[key].default, float):
            min_value_float, max_value_float, step_float = min_value_float, max_value_float, step_float

        int_slider = ipywidgets.IntSlider
        float_slider = ipywidgets.FloatSlider

        if min_value is None:
            min_value = 0
        if max_value is None:
            max_value = 10

        default_slider_layout = ipywidgets.Layout(width='300px', description_width='150px')

        if sig.parameters[key].annotation is str:
            default_value = ipywidgets.Text(value=default_value,
                                       continuous_update=continuous_update)
            exposable = True
        elif sig.parameters[key].annotation is int:
            default_value = int_slider(min=min_value, max=max_value, step=step, value=default_value, continuous_update=continuous_update,layout=default_slider_layout, style={'description_width': 'initial'})
            exposable = True
        elif sig.parameters[key].annotation is float or 'sigma' in key or 'radius' in key:
            default_value = float_slider(min=min_value_float, max=max_value_float, step=step_float, value=default_value, continuous_update=continuous_update,layout=default_slider_layout, style={'description_width': 'initial'})
            exposable = True
        elif key.startswith("is_") or sig.parameters[key].annotation is bool:
            default_value = ipywidgets.Checkbox(value=default_value)
            exposable = True
        elif key == 'footprint' or key == 'selem' or key == 'structuring_element':
            footprint_parameters.append(key)
            default_value = int_slider(min=0, max=10, step=1, value=default_value,  continuous_update=continuous_update)
            exposable = True
        elif context is not None:

            image_parameters.append(key)
            default_value = ipywidgets.Dropdown(
                options=list(context._images.keys()))
            exposable = True

        if exposable:
            if key in kwargs.keys():
                default_value = kwargs[key]
          
            exposable_parameters.append(inspect.Parameter(key, inspect.Parameter.KEYWORD_ONLY, default=default_value))

    viewer_was_none = viewer is None
    if viewer_was_none:
        from ._slice_viewer import _SliceViewer
        from ._crop import _Cropper
        viewer = _SliceViewer(image, zoom_factor=zoom_factor, zoom_spline_order=zoom_spline_order, colormap=colormap, display_min=display_min, display_max=display_max)
        viewer_crop = _Cropper(image, zoom_factor=zoom_factor*7, zoom_spline_order=zoom_spline_order, colormap=colormap, display_min=display_min, display_max=display_max)
    viewer.slice_slider.continuous_update=continuous_update

    execution_blocked = True

    def worker_function(*args, **kwargs):

        command = func_name + "("
        if image_passed:
            command = command + "..."

        for key in [e.name for e in exposable_parameters]:

            if key in footprint_parameters:
                if len(image.shape) == 2:
                    command = command + ", " + key + "=disk(" + str(kwargs[key]) + ")"
                    kwargs[key] = morphology.disk(kwargs[key])
                elif len(image.shape) == 3:
                    command = command + ", " + key + "=ball(" + str(kwargs[key]) + ")"
                    kwargs[key] = morphology.ball(kwargs[key])
            elif key in image_parameters:
                command = command + ", " + key + "=" + str(kwargs[key])
                kwargs[key] = context._images[kwargs[key]]
            else:
                command = command + ", " + key + "=" + str(kwargs[key])
        command = command + ")"

        if not execution_blocked:
            if image_passed:
                result = func(image, *args, **kwargs)
            else:
                result = func(*args, **kwargs)
            viewer.image = result
            viewer_crop.image = result
            viewer_crop.update()

        viewer.slice_slider.max = viewer.image.shape[0] - 1
        viewer.configuration_updated(None)
  
    worker_function.__signature__ = inspect.Signature(exposable_parameters)
    inter = ipywidgets.interactive(worker_function, dict(manual=False, auto_display=False))
    execution_blocked = False
    inter.update()
    
    output_widgets = []
    output_widgets.append(inter)
    if viewer_was_none:
        full_view = ipywidgets.HBox([ipywidgets.VBox([viewer.view], layout=ipywidgets.Layout(object_fit='cover'))], layout=ipywidgets.Layout( object_fit='cover'))
        crop_view = ipywidgets.Box([viewer_crop], layout=ipywidgets.Layout(height='800px',width='1000px',object_fit='cover'))

    box_layout = ipywidgets.Layout(description_width='auto', 
                    border='solid', width='100%', align_items='center', justify_content='center')
    widget_result = ipywidgets.HBox(output_widgets, layout=box_layout)
 
    grid = GridspecLayout(4, 4, height='800px', width='1600px', grid_gap="5px 5px")
    grid[:, 1:] =  crop_view 
    grid[:2, 0:1] = widget_result
    grid[2:, 0:1] = full_view
    
    return grid


def guess_range(name, annotation):
    if name == 'footprint' or name == 'selem' or name == 'structuring_element':
        return 0, 100, 1
    if 'sigma' in name:
        return 0, 25, 1
    if 'radius' in name:
        return 0, 100, 1
    if 'factor' in name:
        return 0, 100, 1
    if name == 'angle' or "degrees" in name:
        return 0, 360, 15
    return None, None, None
