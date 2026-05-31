"""
Declarative tool registry for the KINTSUGI MCP server.

Each entry pairs a public MCP tool name with its handler (``module``/``attr`` in
``kintsugi.mcp.tools``), a top-level ``description``, and the ``input_schema`` that
is advertised to clients verbatim. Keeping this as plain data (no ``mcp`` import)
lets ``kintsugi mcp tools`` and tests read the registry without the optional
``[claude]`` dependency installed.

The input schemas preserve the per-parameter descriptions and ``enum`` constraints
that FastMCP's signature-derived schemas would otherwise drop. Argument *validation*
still runs from each handler's function signature (types) plus the semantic
``path_safety`` checks inside the handler body; this schema is the contract clients
read when deciding how to call a tool.
"""

from __future__ import annotations

from typing import Any

# name -> (handler module, handler attr, description, advertised input schema)
TOOL_SPECS: list[dict[str, Any]] = [
    # ---- Signal Isolation ----
    {
        "name": "load_channel",
        "module": "signal_isolation",
        "attr": "load_channel",
        "description": (
            "Load a channel image from a KINTSUGI project. Returns image metadata "
            "and prepares it for processing."
        ),
        "input_schema": {
            "type": "object",
            "properties": {
                "project_path": {
                    "type": "string",
                    "description": "Path to the KINTSUGI project directory",
                },
                "cycle": {
                    "type": "string",
                    "description": "Cycle name or number (e.g., 'cyc001' or '1')",
                },
                "channel": {
                    "type": "string",
                    "description": "Channel name (e.g., 'CD3', 'DAPI')",
                },
            },
            "required": ["project_path", "cycle", "channel"],
        },
    },
    {
        "name": "subtract_blank",
        "module": "signal_isolation",
        "attr": "subtract_blank",
        "description": (
            "Subtract autofluorescence/blank channel from signal. Uses optimized "
            "Dask operations for large images."
        ),
        "input_schema": {
            "type": "object",
            "properties": {
                "signal_channel": {
                    "type": "string",
                    "description": "Name of loaded signal channel",
                },
                "blank_channel": {
                    "type": "string",
                    "description": "Name of loaded blank channel",
                },
                "blank_clip_factor": {
                    "type": "integer",
                    "description": "Clip blank values below this threshold to 0 (default: 0)",
                    "default": 0,
                },
                "blank_scale_factor": {
                    "type": "number",
                    "description": "Scale factor for blank subtraction (default: 1.0)",
                    "default": 1.0,
                },
                "smooth_low": {
                    "type": "boolean",
                    "description": "Apply smoothing to low-intensity regions",
                    "default": False,
                },
                "smooth_high": {
                    "type": "boolean",
                    "description": "Apply smoothing to high-intensity regions",
                    "default": False,
                },
                "erosion": {
                    "type": "integer",
                    "description": "Erosion disk size for void removal (default: 0)",
                    "default": 0,
                },
            },
            "required": ["signal_channel", "blank_channel"],
        },
    },
    {
        "name": "denoise",
        "module": "signal_isolation",
        "attr": "denoise",
        "description": "Apply denoising filters (percentile, uniform, median) to reduce noise.",
        "input_schema": {
            "type": "object",
            "properties": {
                "channel": {
                    "type": "string",
                    "description": "Name of loaded channel to denoise",
                },
                "method": {
                    "type": "string",
                    "enum": ["percentile", "uniform", "median"],
                    "description": "Denoising method",
                },
                "filter_size": {
                    "type": "integer",
                    "description": "Filter kernel size (default: 3)",
                    "default": 3,
                },
                "percentile": {
                    "type": "integer",
                    "description": "Percentile value for percentile filter (default: 10)",
                    "default": 10,
                },
            },
            "required": ["channel", "method"],
        },
    },
    {
        "name": "apply_clahe",
        "module": "signal_isolation",
        "attr": "apply_clahe",
        "description": (
            "Apply Contrast Limited Adaptive Histogram Equalization to enhance local contrast."
        ),
        "input_schema": {
            "type": "object",
            "properties": {
                "channel": {
                    "type": "string",
                    "description": "Name of loaded channel",
                },
                "clip_limit": {
                    "type": "number",
                    "description": "Clip limit for contrast limiting (default: 0.01)",
                    "default": 0.01,
                },
                "tile_grid_size": {
                    "type": "integer",
                    "description": "Size of grid for histogram equalization (default: 70)",
                    "default": 70,
                },
                "nbins": {
                    "type": "integer",
                    "description": "Number of histogram bins (default: 128)",
                    "default": 128,
                },
            },
            "required": ["channel"],
        },
    },
    {
        "name": "clean_background",
        "module": "signal_isolation",
        "attr": "clean_background",
        "description": "Remove background and small objects from the image.",
        "input_schema": {
            "type": "object",
            "properties": {
                "channel": {
                    "type": "string",
                    "description": "Name of loaded channel",
                },
                "background_threshold": {
                    "type": "integer",
                    "description": "Intensity threshold below which pixels are set to 0 (default: 100)",
                    "default": 100,
                },
                "remove_small_objects": {
                    "type": "boolean",
                    "description": "Whether to remove small connected components",
                    "default": False,
                },
                "min_object_size": {
                    "type": "integer",
                    "description": "Minimum object size to keep (default: 30)",
                    "default": 30,
                },
            },
            "required": ["channel"],
        },
    },
    {
        "name": "gaussian_subtract",
        "module": "signal_isolation",
        "attr": "gaussian_subtract",
        "description": (
            "Subtract Gaussian-blurred version to remove structured background/autofluorescence."
        ),
        "input_schema": {
            "type": "object",
            "properties": {
                "channel": {
                    "type": "string",
                    "description": "Name of loaded channel",
                },
                "sigma": {
                    "type": "integer",
                    "description": "Gaussian sigma value (default: 10)",
                    "default": 10,
                },
                "scale_factor": {
                    "type": "number",
                    "description": "Scale factor for subtraction (default: 0.1)",
                    "default": 0.1,
                },
            },
            "required": ["channel"],
        },
    },
    {
        "name": "denoise_advanced",
        "module": "signal_isolation",
        "attr": "denoise_advanced",
        "description": (
            "Apply advanced denoising: adaptive (auto-select best filter), bilateral "
            "(edge-preserving), nlm (non-local means), bm3d (block matching 3D), n2v "
            "(self-supervised Noise2Void), or patch_svd."
        ),
        "input_schema": {
            "type": "object",
            "properties": {
                "channel": {
                    "type": "string",
                    "description": "Name of loaded channel to denoise",
                },
                "method": {
                    "type": "string",
                    "enum": ["adaptive", "bilateral", "nlm", "bm3d", "n2v", "patch_svd"],
                    "description": "Advanced denoising method (default: adaptive)",
                    "default": "adaptive",
                },
                "strength": {
                    "type": "string",
                    "enum": ["light", "medium", "strong", "auto"],
                    "description": "Denoising strength (default: auto)",
                    "default": "auto",
                },
                "preserve_edges": {
                    "type": "boolean",
                    "description": "Prioritize edge preservation (default: true)",
                    "default": True,
                },
                "n2v_epochs": {
                    "type": "integer",
                    "description": "Training epochs for N2V method (default: 50)",
                    "default": 50,
                },
                "model_path": {
                    "type": "string",
                    "description": "Path to pre-trained N2V model (optional, skip training if provided)",
                },
            },
            "required": ["channel"],
        },
    },
    # ---- Quality Assessment ----
    {
        "name": "assess_quality",
        "module": "quality_assessment",
        "attr": "assess_quality",
        "description": (
            "Assess channel quality using heuristic and/or deep learning metrics. "
            "Returns SNR, contrast, and confidence scores."
        ),
        "input_schema": {
            "type": "object",
            "properties": {
                "channel": {
                    "type": "string",
                    "description": "Name of loaded channel to assess",
                },
                "use_dl_model": {
                    "type": "boolean",
                    "description": "Use deep learning model if available (default: False)",
                    "default": False,
                },
                "confidence_threshold": {
                    "type": "number",
                    "description": "Threshold below which manual review is flagged (default: 0.85)",
                    "default": 0.85,
                },
            },
            "required": ["channel"],
        },
    },
    {
        "name": "compute_snr",
        "module": "quality_assessment",
        "attr": "compute_snr",
        "description": "Compute Signal-to-Noise Ratio for a channel.",
        "input_schema": {
            "type": "object",
            "properties": {
                "channel": {
                    "type": "string",
                    "description": "Name of loaded channel",
                },
            },
            "required": ["channel"],
        },
    },
    # ---- Visualization ----
    {
        "name": "get_image_stats",
        "module": "visualization",
        "attr": "get_image_stats",
        "description": "Get statistics for a loaded image (min, max, mean, std, histogram).",
        "input_schema": {
            "type": "object",
            "properties": {
                "channel": {
                    "type": "string",
                    "description": "Name of loaded channel",
                },
            },
            "required": ["channel"],
        },
    },
    {
        "name": "get_thumbnail",
        "module": "visualization",
        "attr": "get_thumbnail",
        "description": "Get a downsampled thumbnail of the image for preview.",
        "input_schema": {
            "type": "object",
            "properties": {
                "channel": {
                    "type": "string",
                    "description": "Name of loaded channel",
                },
                "max_size": {
                    "type": "integer",
                    "description": "Maximum dimension of thumbnail (default: 512)",
                    "default": 512,
                },
            },
            "required": ["channel"],
        },
    },
    # ---- Workflow ----
    {
        "name": "list_channels",
        "module": "workflow",
        "attr": "list_channels",
        "description": "List all available channels in the current project/cycle.",
        "input_schema": {
            "type": "object",
            "properties": {
                "project_path": {
                    "type": "string",
                    "description": "Path to KINTSUGI project",
                },
                "cycle": {
                    "type": "string",
                    "description": "Cycle name (optional - lists all if not specified)",
                },
            },
            "required": ["project_path"],
        },
    },
    {
        "name": "save_processed",
        "module": "workflow",
        "attr": "save_processed",
        "description": "Save the processed channel to the project output directory.",
        "input_schema": {
            "type": "object",
            "properties": {
                "channel": {
                    "type": "string",
                    "description": "Name of processed channel to save",
                },
                "output_name": {
                    "type": "string",
                    "description": "Output filename (optional - auto-generated if not provided)",
                },
                "format": {
                    "type": "string",
                    "enum": ["tiff", "zarr", "ome-tiff"],
                    "description": "Output format (default: tiff)",
                    "default": "tiff",
                },
            },
            "required": ["channel"],
        },
    },
    {
        "name": "get_processing_history",
        "module": "workflow",
        "attr": "get_processing_history",
        "description": (
            "Get the processing history for a channel, showing all operations applied."
        ),
        "input_schema": {
            "type": "object",
            "properties": {
                "channel": {
                    "type": "string",
                    "description": "Name of channel",
                },
            },
            "required": ["channel"],
        },
    },
    {
        "name": "suggest_parameters",
        "module": "workflow",
        "attr": "suggest_parameters",
        "description": (
            "Analyze channel and suggest processing parameters based on image characteristics."
        ),
        "input_schema": {
            "type": "object",
            "properties": {
                "channel": {
                    "type": "string",
                    "description": "Name of loaded channel",
                },
                "blank_channel": {
                    "type": "string",
                    "description": "Name of blank channel for AF estimation (optional)",
                },
            },
            "required": ["channel"],
        },
    },
    {
        "name": "generate_jupyter_cell",
        "module": "workflow",
        "attr": "generate_jupyter_cell",
        "description": (
            "Generate a Jupyter notebook cell for interactive parameter tuning with Kview2 widgets."
        ),
        "input_schema": {
            "type": "object",
            "properties": {
                "operation": {
                    "type": "string",
                    "enum": [
                        "blank_subtraction",
                        "denoise",
                        "clahe",
                        "clean",
                        "gaussian_subtract",
                        "full_pipeline",
                    ],
                    "description": "Operation to create interactive tuner for",
                },
                "channel": {
                    "type": "string",
                    "description": "Name of loaded channel",
                },
                "initial_params": {
                    "type": "object",
                    "description": "Initial parameter values (optional)",
                },
            },
            "required": ["operation", "channel"],
        },
    },
    # ---- Channel Clustering ----
    {
        "name": "cluster_channels",
        "module": "signal_isolation",
        "attr": "cluster_channels_tool",
        "description": (
            "Cluster channels by image feature similarity. Groups channels that share "
            "similar intensity profiles, SNR, and texture so parameters can be tuned "
            "once per cluster."
        ),
        "input_schema": {
            "type": "object",
            "properties": {
                "project_path": {
                    "type": "string",
                    "description": "Path to KINTSUGI project directory",
                },
                "n_clusters": {
                    "type": "integer",
                    "description": "Number of clusters (auto-select if omitted)",
                },
                "wavelength_aware": {
                    "type": "boolean",
                    "description": "Enforce same-wavelength clustering (default: true)",
                    "default": True,
                },
            },
            "required": ["project_path"],
        },
    },
    {
        "name": "propagate_parameters",
        "module": "signal_isolation",
        "attr": "propagate_parameters_tool",
        "description": (
            "Apply tuned signal isolation parameters from a representative channel to all "
            "members of a cluster."
        ),
        "input_schema": {
            "type": "object",
            "properties": {
                "project_path": {
                    "type": "string",
                    "description": "Path to KINTSUGI project directory",
                },
                "cluster_id": {
                    "type": "integer",
                    "description": "Cluster ID to propagate parameters to",
                },
                "params": {
                    "type": "object",
                    "description": "Processing parameters (blank_params and/or clean_params dicts)",
                },
                "output_dir": {
                    "type": "string",
                    "description": "Output directory (default: project signal_isolated dir)",
                },
            },
            "required": ["project_path", "cluster_id", "params"],
        },
    },
    # ---- Parameter Learning ----
    {
        "name": "get_learned_parameters",
        "module": "learning",
        "attr": "get_learned_parameters",
        "description": (
            "Get recommended parameters based on learned history from similar tissue/marker "
            "combinations. Returns parameters with confidence scores."
        ),
        "input_schema": {
            "type": "object",
            "properties": {
                "tissue_type": {
                    "type": "string",
                    "description": "Type of tissue (e.g., 'tonsil', 'lymph node', 'pancreas')",
                },
                "marker_name": {
                    "type": "string",
                    "description": "Marker name (e.g., 'CD3', 'CD20', 'DAPI')",
                },
                "operation": {
                    "type": "string",
                    "description": "Processing operation (e.g., 'blank_subtraction', 'denoise')",
                },
                "project_path": {
                    "type": "string",
                    "description": "Project path for project-specific learning (optional)",
                },
                "channel": {
                    "type": "string",
                    "description": "Channel name to get image characteristics from (optional)",
                },
            },
            "required": ["tissue_type", "marker_name", "operation"],
        },
    },
    {
        "name": "record_successful_parameters",
        "module": "learning",
        "attr": "record_successful_parameters",
        "description": (
            "Record a successful parameter set for future learning. Call after user approves "
            "processing results."
        ),
        "input_schema": {
            "type": "object",
            "properties": {
                "tissue_type": {
                    "type": "string",
                    "description": "Type of tissue processed",
                },
                "marker_name": {
                    "type": "string",
                    "description": "Marker name",
                },
                "operation": {
                    "type": "string",
                    "description": "Processing operation",
                },
                "parameters": {
                    "type": "object",
                    "description": "Parameter dictionary that was used",
                },
                "quality_before": {
                    "type": "number",
                    "description": "Quality score before processing",
                    "default": 0.0,
                },
                "quality_after": {
                    "type": "number",
                    "description": "Quality score after processing",
                    "default": 0.0,
                },
                "user_approved": {
                    "type": "boolean",
                    "description": "Whether user approved the result",
                    "default": True,
                },
                "user_notes": {
                    "type": "string",
                    "description": "Optional notes from user",
                },
                "channel": {
                    "type": "string",
                    "description": "Channel name (for characteristics)",
                },
                "project_path": {
                    "type": "string",
                    "description": "Project path",
                },
            },
            "required": ["tissue_type", "marker_name", "operation", "parameters"],
        },
    },
    {
        "name": "suggest_with_learning",
        "module": "learning",
        "attr": "suggest_with_learning",
        "description": (
            "Get comprehensive parameter suggestions combining image analysis AND learned "
            "history. This is the primary recommendation tool."
        ),
        "input_schema": {
            "type": "object",
            "properties": {
                "channel": {
                    "type": "string",
                    "description": "Loaded channel name",
                },
                "tissue_type": {
                    "type": "string",
                    "description": "Tissue type",
                },
                "marker_name": {
                    "type": "string",
                    "description": "Marker name",
                },
                "blank_channel": {
                    "type": "string",
                    "description": "Blank channel for AF estimation (optional)",
                },
                "project_path": {
                    "type": "string",
                    "description": "Project path (optional)",
                },
            },
            "required": ["channel", "tissue_type", "marker_name"],
        },
    },
    {
        "name": "approve_and_learn",
        "module": "learning",
        "attr": "approve_and_learn",
        "description": (
            "Approve processing results and record ALL parameters for learning. Call when "
            "user approves final processed result."
        ),
        "input_schema": {
            "type": "object",
            "properties": {
                "channel": {
                    "type": "string",
                    "description": "Channel name",
                },
                "tissue_type": {
                    "type": "string",
                    "description": "Tissue type",
                },
                "marker_name": {
                    "type": "string",
                    "description": "Marker name",
                },
                "operations_params": {
                    "type": "object",
                    "description": "Dict mapping operation names to their parameters",
                },
                "quality_before": {
                    "type": "number",
                    "description": "Quality score before all processing",
                },
                "quality_after": {
                    "type": "number",
                    "description": "Quality score after all processing",
                },
                "notes": {
                    "type": "string",
                    "description": "User notes",
                    "default": "",
                },
                "project_path": {
                    "type": "string",
                    "description": "Project path",
                },
            },
            "required": [
                "channel",
                "tissue_type",
                "marker_name",
                "operations_params",
                "quality_before",
                "quality_after",
            ],
        },
    },
    {
        "name": "get_learning_statistics",
        "module": "learning",
        "attr": "get_learning_statistics",
        "description": "Get statistics about the parameter learning database.",
        "input_schema": {
            "type": "object",
            "properties": {
                "project_path": {
                    "type": "string",
                    "description": "Project path (optional)",
                },
            },
            "required": [],
        },
    },
    # ---- Automated Optimization (Tier 2) ----
    {
        "name": "optimize_parameters",
        "module": "signal_isolation",
        "attr": "optimize_parameters",
        "description": (
            "Find optimal signal isolation parameters via Bayesian optimization (Optuna). "
            "Automatically searches the parameter space — no manual slider adjustment needed."
        ),
        "input_schema": {
            "type": "object",
            "properties": {
                "signal_channel": {
                    "type": "string",
                    "description": "Name of loaded signal channel",
                },
                "blank_channel": {
                    "type": "string",
                    "description": "Name of loaded blank channel",
                },
                "n_trials": {
                    "type": "integer",
                    "description": "Maximum optimization trials (default: 80)",
                    "default": 80,
                },
                "timeout": {
                    "type": "integer",
                    "description": "Maximum seconds for optimization (default: 300)",
                    "default": 300,
                },
                "optimize_clean": {
                    "type": "boolean",
                    "description": "Also optimize background cleaning parameters (default: true)",
                    "default": True,
                },
                "warm_start_params": {
                    "type": "object",
                    "description": "Initial parameter guess for faster convergence (optional)",
                },
                "project_path": {
                    "type": "string",
                    "description": "Project path for storing optimization study (optional)",
                },
            },
            "required": ["signal_channel", "blank_channel"],
        },
    },
    {
        "name": "predict_parameters",
        "module": "signal_isolation",
        "attr": "predict_parameters",
        "description": (
            "Predict signal isolation parameters from image features using a trained Random "
            "Forest model. Includes confidence scores and uncertainty estimates."
        ),
        "input_schema": {
            "type": "object",
            "properties": {
                "signal_channel": {
                    "type": "string",
                    "description": "Name of loaded signal channel",
                },
                "blank_channel": {
                    "type": "string",
                    "description": "Name of loaded blank channel (optional)",
                },
                "project_path": {
                    "type": "string",
                    "description": "Project path containing predictor model (optional)",
                },
                "model_path": {
                    "type": "string",
                    "description": "Explicit path to predictor model .joblib (optional)",
                },
            },
            "required": ["signal_channel"],
        },
    },
    {
        "name": "estimate_background",
        "module": "signal_isolation",
        "attr": "estimate_background",
        "description": (
            "Estimate and subtract background using SMO (Silver Mountain Operator). "
            "Parameter-free — no blank channel needed."
        ),
        "input_schema": {
            "type": "object",
            "properties": {
                "channel": {
                    "type": "string",
                    "description": "Name of loaded channel image",
                },
                "kernel_size": {
                    "type": "integer",
                    "description": "SMO kernel size (default: 7, robust across 5-15)",
                    "default": 7,
                },
                "return_background": {
                    "type": "boolean",
                    "description": "Store estimated background as a loaded image (default: false)",
                    "default": False,
                },
            },
            "required": ["channel"],
        },
    },
]
