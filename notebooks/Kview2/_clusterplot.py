import numpy as np
from ._grid import grid
from ._curtain import curtain
from ._slice import slice
from ._scatterplot import scatterplot
import functools


def clusterplot(image, labels, df, column_x:str="x", column_y:str="y", column_selection:str="selection", 
                     zoom_factor:float=1, zoom_spline_order:int=0, figsize=(4, 4), alpha:float=0.5, markersize=4):


    print(f"Shape line 45: {labels.shape}")
    print(f"Labels shape: {labels.shape if labels is not None else None}")
    print(f"DataFrame shape: {df.shape if df is not None else None}")
    labels = np.asarray(labels)

    if column_selection in df.columns:
        selection = df["selection"].tolist()
        # Create a lookup array with shape (n_labels + 1,)
        lookup = np.asarray([-1] + list(selection)) * 1 + 1
        # Use where to maintain 2D structure
        selected_image = np.where(labels > 0, lookup[labels], 0).astype(np.uint32)
    else:
        selected_image = ((labels > 0) * 1).astype(np.uint32)

    if image is None:
        image_display = slice(selected_image, zoom_factor=zoom_factor, zoom_spline_order=zoom_spline_order)
    else:
        image_display = curtain(image, selected_image, zoom_factor=zoom_factor, zoom_spline_order=zoom_spline_order, alpha=alpha)

    def update(selection, label_image, selected_image, widget):
        temp = np.take(np.asarray([-1] + list(selection)) * 1 + 1, label_image)

        # overwrite the pixels in the given image
        np.copyto(selected_image, temp.astype(selected_image.dtype))

        # redraw the visualization
        widget.update()

    update_selection = functools.partial(update, label_image=labels, selected_image=selected_image,
                                         widget=image_display)

    scatterplot = scatterplot(df, column_x, column_y, column_selection, figsize=figsize,
                                        selection_changed_callback=update_selection, markersize=markersize)

    return grid([[
        image_display,
        scatterplot,

    ]])
