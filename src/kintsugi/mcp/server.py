"""
KINTSUGI MCP Server

Exposes image processing tools for Claude Code integration.
"""

from __future__ import annotations

import asyncio
import json
import logging
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional

try:
    from mcp.server import Server
    from mcp.server.stdio import stdio_server
    from mcp.types import (
        Tool,
        TextContent,
        ImageContent,
        EmbeddedResource,
    )
    MCP_AVAILABLE = True
except ImportError:
    MCP_AVAILABLE = False
    Server = None

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
)
logger = logging.getLogger("kintsugi.mcp")


# Global state for the current project/session
_session_state: Dict[str, Any] = {
    "project_path": None,
    "current_channel": None,
    "loaded_images": {},
    "processing_history": [],
}


def create_server() -> "Server":
    """Create and configure the KINTSUGI MCP server."""
    if not MCP_AVAILABLE:
        raise ImportError(
            "MCP package not installed. Install with: pip install kintsugi[claude]"
        )

    server = Server("kintsugi")

    # Import tool modules
    from kintsugi.mcp.tools import signal_isolation, quality_assessment, visualization, workflow

    # Register tools from each module
    _register_signal_isolation_tools(server)
    _register_quality_assessment_tools(server)
    _register_visualization_tools(server)
    _register_workflow_tools(server)

    return server


def _register_signal_isolation_tools(server: "Server"):
    """Register signal isolation tools."""

    @server.list_tools()
    async def list_tools() -> List[Tool]:
        """List available tools."""
        return [
            # Signal Isolation Tools
            Tool(
                name="load_channel",
                description="Load a channel image from a KINTSUGI project. Returns image metadata and prepares it for processing.",
                inputSchema={
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
            ),
            Tool(
                name="subtract_blank",
                description="Subtract autofluorescence/blank channel from signal. Uses optimized Dask operations for large images.",
                inputSchema={
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
            ),
            Tool(
                name="denoise",
                description="Apply denoising filters (percentile, uniform, median) to reduce noise.",
                inputSchema={
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
            ),
            Tool(
                name="apply_clahe",
                description="Apply Contrast Limited Adaptive Histogram Equalization to enhance local contrast.",
                inputSchema={
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
            ),
            Tool(
                name="clean_background",
                description="Remove background and small objects from the image.",
                inputSchema={
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
            ),
            Tool(
                name="gaussian_subtract",
                description="Subtract Gaussian-blurred version to remove structured background/autofluorescence.",
                inputSchema={
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
            ),
            Tool(
                name="denoise_advanced",
                description="Apply advanced denoising: adaptive (auto-select best filter), bilateral (edge-preserving), nlm (non-local means), bm3d (block matching 3D), n2v (self-supervised Noise2Void), or patch_svd.",
                inputSchema={
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
            ),
            # Quality Assessment Tools
            Tool(
                name="assess_quality",
                description="Assess channel quality using heuristic and/or deep learning metrics. Returns SNR, contrast, and confidence scores.",
                inputSchema={
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
            ),
            Tool(
                name="compute_snr",
                description="Compute Signal-to-Noise Ratio for a channel.",
                inputSchema={
                    "type": "object",
                    "properties": {
                        "channel": {
                            "type": "string",
                            "description": "Name of loaded channel",
                        },
                    },
                    "required": ["channel"],
                },
            ),
            # Visualization Tools
            Tool(
                name="get_image_stats",
                description="Get statistics for a loaded image (min, max, mean, std, histogram).",
                inputSchema={
                    "type": "object",
                    "properties": {
                        "channel": {
                            "type": "string",
                            "description": "Name of loaded channel",
                        },
                    },
                    "required": ["channel"],
                },
            ),
            Tool(
                name="get_thumbnail",
                description="Get a downsampled thumbnail of the image for preview.",
                inputSchema={
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
            ),
            # Workflow Tools
            Tool(
                name="list_channels",
                description="List all available channels in the current project/cycle.",
                inputSchema={
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
            ),
            Tool(
                name="save_processed",
                description="Save the processed channel to the project output directory.",
                inputSchema={
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
            ),
            Tool(
                name="get_processing_history",
                description="Get the processing history for a channel, showing all operations applied.",
                inputSchema={
                    "type": "object",
                    "properties": {
                        "channel": {
                            "type": "string",
                            "description": "Name of channel",
                        },
                    },
                    "required": ["channel"],
                },
            ),
            Tool(
                name="suggest_parameters",
                description="Analyze channel and suggest processing parameters based on image characteristics.",
                inputSchema={
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
            ),
            Tool(
                name="generate_jupyter_cell",
                description="Generate a Jupyter notebook cell for interactive parameter tuning with Kview2 widgets.",
                inputSchema={
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
            ),
            # Parameter Learning Tools
            Tool(
                name="get_learned_parameters",
                description="Get recommended parameters based on learned history from similar tissue/marker combinations. Returns parameters with confidence scores.",
                inputSchema={
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
            ),
            Tool(
                name="record_successful_parameters",
                description="Record a successful parameter set for future learning. Call after user approves processing results.",
                inputSchema={
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
            ),
            Tool(
                name="suggest_with_learning",
                description="Get comprehensive parameter suggestions combining image analysis AND learned history. This is the primary recommendation tool.",
                inputSchema={
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
            ),
            Tool(
                name="approve_and_learn",
                description="Approve processing results and record ALL parameters for learning. Call when user approves final processed result.",
                inputSchema={
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
                    "required": ["channel", "tissue_type", "marker_name", "operations_params", "quality_before", "quality_after"],
                },
            ),
            Tool(
                name="get_learning_statistics",
                description="Get statistics about the parameter learning database.",
                inputSchema={
                    "type": "object",
                    "properties": {
                        "project_path": {
                            "type": "string",
                            "description": "Project path (optional)",
                        },
                    },
                    "required": [],
                },
            ),
        ]

    @server.call_tool()
    async def call_tool(name: str, arguments: Dict[str, Any]) -> List[TextContent]:
        """Handle tool calls."""
        try:
            # Import tool implementations
            from kintsugi.mcp.tools import signal_isolation, quality_assessment, visualization, workflow, learning

            # Route to appropriate handler
            if name == "load_channel":
                result = await signal_isolation.load_channel(**arguments)
            elif name == "subtract_blank":
                result = await signal_isolation.subtract_blank(**arguments)
            elif name == "denoise":
                result = await signal_isolation.denoise(**arguments)
            elif name == "apply_clahe":
                result = await signal_isolation.apply_clahe(**arguments)
            elif name == "clean_background":
                result = await signal_isolation.clean_background(**arguments)
            elif name == "gaussian_subtract":
                result = await signal_isolation.gaussian_subtract(**arguments)
            elif name == "denoise_advanced":
                result = await signal_isolation.denoise_advanced(**arguments)
            elif name == "assess_quality":
                result = await quality_assessment.assess_quality(**arguments)
            elif name == "compute_snr":
                result = await quality_assessment.compute_snr(**arguments)
            elif name == "get_image_stats":
                result = await visualization.get_image_stats(**arguments)
            elif name == "get_thumbnail":
                result = await visualization.get_thumbnail(**arguments)
            elif name == "list_channels":
                result = await workflow.list_channels(**arguments)
            elif name == "save_processed":
                result = await workflow.save_processed(**arguments)
            elif name == "get_processing_history":
                result = await workflow.get_processing_history(**arguments)
            elif name == "suggest_parameters":
                result = await workflow.suggest_parameters(**arguments)
            elif name == "generate_jupyter_cell":
                result = await workflow.generate_jupyter_cell(**arguments)
            # Learning tools
            elif name == "get_learned_parameters":
                result = await learning.get_learned_parameters(**arguments)
            elif name == "record_successful_parameters":
                result = await learning.record_successful_parameters(**arguments)
            elif name == "suggest_with_learning":
                result = await learning.suggest_with_learning(**arguments)
            elif name == "approve_and_learn":
                result = await learning.approve_and_learn(**arguments)
            elif name == "get_learning_statistics":
                result = await learning.get_learning_statistics(**arguments)
            else:
                result = {"error": f"Unknown tool: {name}"}

            return [TextContent(type="text", text=json.dumps(result, indent=2))]

        except Exception as e:
            logger.exception(f"Error in tool {name}")
            return [
                TextContent(
                    type="text",
                    text=json.dumps({"error": str(e), "tool": name}, indent=2),
                )
            ]

    return server


async def run_server():
    """Run the MCP server using stdio transport."""
    if not MCP_AVAILABLE:
        print("Error: MCP package not installed. Install with: pip install kintsugi[claude]")
        sys.exit(1)

    server = create_server()
    logger.info("Starting KINTSUGI MCP server...")

    async with stdio_server() as (read_stream, write_stream):
        await server.run(read_stream, write_stream, server.create_initialization_options())


def main():
    """Entry point for the MCP server."""
    asyncio.run(run_server())


if __name__ == "__main__":
    main()
