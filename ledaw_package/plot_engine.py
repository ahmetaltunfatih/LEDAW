import os
import shutil
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from .nbody_engine import normalize_path


def plot_heatmap_std_led(file_name, sheet_name, figsize, vmin, vmax, fig_format, set_dpi, save_as='default', cutoff_annot=None, submatrix_coords_to_be_highlighted=None, display_heatmap=False):
    """
    Draws a standard LED heat map with optional annotation settings and highlighted submatrices.

    Parameters:
    -----------
    file_name : str
        The name of the Excel file containing the data.
    sheet_name : str
        The sheet name within the Excel file to read the data from.
    figsize : tuple
        A tuple specifying the size of the figure (width, height).
    vmin : float
        The minimum value to map to the colormap. Data values lower than vmin will be clamped to the minimum color.
    vmax : float
        The maximum value to map to the colormap. Data values higher than vmax will be clamped to the maximum color.
    fig_format : str
        The format to save the figure in (e.g., 'png', 'jpg', 'tif', 'svg', 'eps').
	set_dpi : integer
		Dpi of the figures to be saved
	save_as : str
		The path where the figure will be saved plus file name. By default, it is written to the working directory as
		sheet_name.fig_format
    cutoff_annot : float, optional
        If specified, values between -cutoff_annot and cutoff_annot will not be annotated on the heatmap.
        Default is None, meaning all values will be annotated.
    submatrix_coords_to_be_highlighted : tuple of tuples, optional
        A tuple specifying the coordinates of the submatrix to be highlighted with a black border.
        Format: ((starting row, starting column), (ending row, ending column)).
        Default is None, meaning no submatrix will be highlighted.
	display_heatmap: boolean
		If True, it displays the plot in the environment where the code is executed.

    Details:
    --------
    - vmin and vmax control the color mapping of the heatmap. Values outside this range are clamped to the closest color.
    - The heatmap will have a #f5ecb3 grid separating the cells, and the color bar will display the range of values.

    Example Usage:
    --------------
    hm(file_name="data.xlsx", sheet_name="Sheet1", figsize=(12, 12), vmin=-100, vmax=100, fig_format="png", 
       cutoff_annot=1.0, submatrix_coords_to_be_highlighted=((0, 1), (0, 13))
    """

    plt.close('all')  # Close any existing figures to prevent interference

    df = pd.read_excel(file_name, sheet_name, index_col=0)
    df.sort_index(axis=0, ascending=True, inplace=True)
    df.sort_index(axis=1, ascending=True, inplace=True)
    df = df.apply(lambda x: x.map(lambda y: 0.0 if abs(y) < 0.05 else y))
    
    # Replace values between -cutoff_annot and +cutoff_annot with NaN
    if cutoff_annot is not None:
        df = df.where((df < -cutoff_annot) | (df > cutoff_annot))

    plt.figure(figsize=figsize)

    # Define custom colormap with blue-white-red gradient
    cmap_colors = [(1, 0, 0), (1, 1, 1), (0, 0, 1)]
    cmap = LinearSegmentedColormap.from_list("custom_cmap", cmap_colors)

    # Create the heatmap without border lines
    ax = sns.heatmap(df,
                     annot=True,
                     annot_kws={"size": 13, "color": "black", "fontfamily": "Arial"}, 
                     fmt=".1f",
                     linewidths=0,
                     cmap=cmap,
                     center=0.0,
                     vmin=vmin,
                     vmax=vmax,
                     cbar_kws={"shrink": 0.50, "orientation": "vertical", "pad": 0.05},
                     square=True,
                     mask=df.isna())

    # Adjust the color bar text size and font
    cbar = ax.collections[0].colorbar
    cbar.ax.tick_params(labelsize=12, width=2) 
    for label in cbar.ax.get_yticklabels():
        label.set_fontfamily('Arial')
 
    # Add a black box around the color bar
    for spine in cbar.ax.spines.values():
        spine.set_edgecolor('black')
        spine.set_linewidth(2)        

    # Set font family for x and y tick labels
    ax.set_xticklabels(ax.get_xticklabels(), fontfamily='Arial', fontsize=13)
    ax.set_yticklabels(ax.get_yticklabels(), fontfamily='Arial', fontsize=13)

    # Use a custom grid for the upper triangle and the diagonals
    for i in range(df.shape[0]):
        for j in range(i + 1, df.shape[1]):
            ax.add_patch(plt.Rectangle((j, i), 1, 1, fill=False, edgecolor='#f5ecb3', lw=3))

    # Add an extra box-shaped border for the diagonal cells
    for i in range(df.shape[1]):
        ax.add_patch(plt.Rectangle((i, i), 1, 1, fill=False, edgecolor='#f5ecb3', lw=3))

    # Adjust the outer frame at the top and right
    ax.axhline(y=0, color='#f5ecb3', linewidth=6)
    ax.axvline(x=df.shape[0], color='#f5ecb3', linewidth=7)
            
    # Adjustment of the lengths of lines to be added to the left of the first diagonal and to the bottom of the last diagonal
    linewidth = 5
    adjustment_factor = linewidth / 2.0 / 72.0
    
    # Add a line at the bottom of the last cell of the last row
    last_row_index = df.shape[0] - 1
    last_column_index = df.shape[1] - 1
    bottom_line_start = (last_column_index + adjustment_factor, last_row_index + 1 - adjustment_factor)
    bottom_line_end = (last_column_index + 1 - adjustment_factor, last_row_index + 1 - adjustment_factor)
    ax.add_line(plt.Line2D([bottom_line_start[0], bottom_line_end[0]], [bottom_line_start[1], bottom_line_end[1]], color='#f5ecb3', linewidth=linewidth))
    
    # Add a line at the left of the first cell of the first row
    first_cell = (0, 0)
    left_line_end = (0, 1 - adjustment_factor)
    left_line_start = (0, 0 - adjustment_factor)
    ax.add_line(plt.Line2D([left_line_start[0], left_line_end[0]], [left_line_start[1], left_line_end[1]], color='#f5ecb3', linewidth=linewidth))
    
    # Surround submatrix by a black box from the specified coordinates in submatrix_coords_to_be_highlighted
    if submatrix_coords_to_be_highlighted is not None:
        ul_row, ul_col = submatrix_coords_to_be_highlighted[0]  # Upper-left corner (row, col)
        lr_row, lr_col = submatrix_coords_to_be_highlighted[1]  # Lower-right corner (row, col)

        # Draw the left border
        if ul_col == 0:
            ax.add_line(plt.Line2D([ul_col, ul_col], [ul_row, lr_row + 0.9], color='black', linewidth=6))
        else:
            ax.add_line(plt.Line2D([ul_col, ul_col], [ul_row, lr_row + 1], color='black', linewidth=3))

        # Draw the top border
        if ul_row == 0:
            ax.add_line(plt.Line2D([ul_col + adjustment_factor + 0.02, lr_col + 0.9], [ul_row, ul_row], color='black', linewidth=6))
        else:
            ax.add_line(plt.Line2D([ul_col + adjustment_factor, lr_col + 0.99], [ul_row, ul_row], color='black', linewidth=2.8))

        # Draw the right border
        if lr_col == df.shape[1] - 1:
            ax.add_line(plt.Line2D([lr_col + 1 , lr_col + 1], [ul_row - 1.5*adjustment_factor, lr_row + 1 - 1.5*adjustment_factor], color='black', linewidth=7))
        else:
            ax.add_line(plt.Line2D([lr_col + 1 , lr_col + 1], [ul_row, lr_row + 1], color='black', linewidth=3))
     
        # Draw the bottom border
        ax.add_line(plt.Line2D([ul_col + 0.15*adjustment_factor, lr_col + 0.98], [lr_row + 1, lr_row + 1], color='black', linewidth=3.3))
        
    ax.tick_params(labelright=True, labeltop=True, labelleft=False, labelbottom=False, size=0)

    plt.xticks(fontsize=13)
    plt.yticks(fontsize=13, rotation=0)

    # Set default save_as if set to "default"
    if save_as == "default":
        save_as = f"{sheet_name}.{fig_format}"

    plt.savefig(save_as, dpi=set_dpi, format=fig_format, bbox_inches='tight')

    if display_heatmap:
        plt.show()
    else:
        plt.close()


def plot_heatmap_fp_led(file_name, sheet_name, figsize, vmin, vmax, fig_format, set_dpi, save_as='default', cutoff_annot=None, submatrix_coords_to_be_highlighted=None, display_heatmap=False):
    """
    Draws a fp-LED heat map with optional annotation settings and highlighted submatrices.

    Parameters:
    -----------
    file_name : str
        The name of the Excel file containing the data.
    sheet_name : str
        The sheet name within the Excel file to read the data from.
    figsize : tuple
        A tuple specifying the size of the figure (width, height).
    vmin : float
        The minimum value to map to the colormap. Data values lower than vmin will be clamped to the minimum color.
    vmax : float
        The maximum value to map to the colormap. Data values higher than vmax will be clamped to the maximum color.
    fig_format : str
        The format to save the figure in (e.g., 'png', 'jpg', 'tif', 'svg', 'eps').
	set_dpi : integer
		Dpi of the figures to be saved
	save_as : str
		The path where the figure will be saved plus file name. By default, it is written to the working directory as
		sheet_name.fig_format
    cutoff_annot : float, optional
        If specified, values between -cutoff_annot and cutoff_annot will not be annotated on the heatmap.
        Default is None, meaning all values will be annotated.
    submatrix_coords_to_be_highlighted : tuple of tuples, optional
        A tuple specifying the coordinates of the submatrix to be highlighted with a black border.
        Format: ((top left corner row index, top left corner column index), (bottom right cell's row index, bottom right column index)).
        Default is None, meaning no submatrix will be highlighted.
	display_heatmap: boolean
		If True, it displays the plot in the environment where the code is executed.

    Details:
    --------
    - vmin and vmax control the color mapping of the heatmap. Values outside this range are clamped to the closest color.
    - The heatmap will have a light yellow grid separating the cells, and the color bar will display the range of values.

    Example Usage:
    --------------
    hm(file_name="data.xlsx", sheet_name="Sheet1", figsize=(12, 12), vmin=-100, vmax=100, fig_format="png", 
       cutoff_annot=1.0, submatrix_coords_to_be_highlighted=((0, 1), (0, 13))
    """

    plt.close('all')  # Close any existing figures to prevent interference

    df = pd.read_excel(file_name, sheet_name, index_col=0)
    df.sort_index(axis=0, ascending=True, inplace=True)
    df.sort_index(axis=1, ascending=True, inplace=True)
    df = df.apply(lambda x: x.map(lambda y: 0.0 if abs(y) < 0.05 else y))
    
    # Replace values between -cutoff_annot and +cutoff_annot with NaN
    if cutoff_annot is not None:
        df = df.where((df < -cutoff_annot) | (df > cutoff_annot))
    
    plt.figure(figsize=figsize)

    # Define custom colormap with blue-white-red gradient
    cmap_colors = [(1, 0, 0), (1, 1, 1), (0, 0, 1)]
    cmap = LinearSegmentedColormap.from_list("custom_cmap", cmap_colors)

    # Create the heatmap without border lines
    ax = sns.heatmap(df,
                     annot=True,
                     annot_kws={"size": 13, "color": "black", "fontfamily": "Arial"},
                     fmt=".1f",
                     linewidths=0,
                     cmap=cmap,
                     center=0.0,
                     vmin=vmin,
                     vmax=vmax,
                     cbar_kws={"shrink": 0.50, "orientation": "vertical", "pad": 0.05},
                     square=True,
                     mask=df.isna()) 

    # Adjust the color bar text size and font
    cbar = ax.collections[0].colorbar
    cbar.ax.tick_params(labelsize=12, width=2) 
    for label in cbar.ax.get_yticklabels():
        label.set_fontfamily('Arial')
 
    # Add a black box around the color bar
    for spine in cbar.ax.spines.values():
        spine.set_edgecolor('black')
        spine.set_linewidth(2)        
    
    # Use a custom grid for the upper triangle
    for i in range(df.shape[0]):
        for j in range(i + 1, df.shape[1]):
            ax.add_patch(plt.Rectangle((j, i), 1, 1, fill=False, edgecolor='#f5ecb3', lw=3))

    # Adjustment of the lengths of lines to be added to the left of first diagonal and to the bottom of the last diagonal
    linewidth = 5
    adjustment_factor = linewidth / 2.0 / 72.0

    # Add a border line starting from column 1
    border_line_start_x = 1 + adjustment_factor
    border_line_end_x = df.shape[1] 
    ax.add_line(plt.Line2D([border_line_start_x, border_line_end_x], [0, 0], color='#f5ecb3', linewidth=linewidth))

    # Add a vertical border line at the right
    border_line_start_y = 0
    border_line_end_y = df.shape[0]-1 - adjustment_factor
    ax.add_line(plt.Line2D([df.shape[1], df.shape[1]], [border_line_start_y, border_line_end_y], color='#f5ecb3', linewidth=6))    

    # Surround submatrix by a black box from the specified coordinates in submatrix_coords_to_be_highlighted
    if submatrix_coords_to_be_highlighted is not None:
        ul_row, ul_col = submatrix_coords_to_be_highlighted[0]  # Upper-left corner (row, col)
        lr_row, lr_col = submatrix_coords_to_be_highlighted[1]  # Lower-right corner (row, col)

        # Draw the left border
        if ul_col == 0:
            ax.add_line(plt.Line2D([ul_col, ul_col], [ul_row, lr_row + 0.9], color='black', linewidth=6))
        else:
            ax.add_line(plt.Line2D([ul_col, ul_col], [ul_row, lr_row + 1], color='black', linewidth=3))

        # Draw the top border
        if ul_row == 0:
            ax.add_line(plt.Line2D([ul_col + adjustment_factor + 0.02, lr_col + 0.9], [ul_row, ul_row], color='black', linewidth=6))
        else:
            ax.add_line(plt.Line2D([ul_col + adjustment_factor, lr_col + 0.99], [ul_row, ul_row], color='black', linewidth=2.8))

        # Draw the right border
        if lr_col == df.shape[1] - 1:
            ax.add_line(plt.Line2D([lr_col + 1 , lr_col + 1], [ul_row - 1.5*adjustment_factor, lr_row + 1 - 1.5*adjustment_factor], color='black', linewidth=7))
        else:
            ax.add_line(plt.Line2D([lr_col + 1 , lr_col + 1], [ul_row, lr_row + 1], color='black', linewidth=3))
     
        # Draw the bottom border
        ax.add_line(plt.Line2D([ul_col + 0.15*adjustment_factor, lr_col + 0.98], [lr_row + 1, lr_row + 1], color='black', linewidth=3.3))
    
    # Set font family for x and y tick labels
    ax.set_xticklabels(ax.get_xticklabels(), fontfamily='Arial', fontsize=13)
    ax.set_yticklabels(ax.get_yticklabels(), fontfamily='Arial', fontsize=13)

    # Get rid of the first x_ticklabel and the last y_ticklabel
    ax.set_xticks(ax.get_xticks()[1:])
    ax.set_yticks(ax.get_yticks()[:-1])
        
    ax.tick_params(labelright=True, labeltop=True, labelleft=False, labelbottom=False, size=0)

    plt.xticks(fontsize=13)
    plt.yticks(fontsize=13, rotation=0)

    # Set default save_as if set to "default"
    if save_as == "default":
        save_as = f"{sheet_name}.{fig_format}"

    plt.savefig(save_as, dpi=set_dpi, format=fig_format, bbox_inches='tight')

    if display_heatmap:
        plt.show()
    else:
        plt.close()


def find_subdirs_under_base_path(base_path, depth_level=2):
    """Find all subdirectories under the given base path up to the specified depth level."""
    subdirectories = []
    for root, dirs, files in os.walk(base_path):
        if root.count(os.sep) == base_path.count(os.sep) + depth_level:
            # Only process directories containing the Excel files
            if "Summary_Standard_LED_matrices.xlsx" in files or "Summary_fp-LED_matrices.xlsx" in files:
                subdirectories.append(root)
    return subdirectories


def delete_existing_plot_directories_first(directory):
    """Delete the HEAT-MAP-STD-LED and HEAT-MAP-fp-LED directories if they exist."""
    std_output_dir = os.path.join(directory, "HEAT-MAP-STD-LED")
    fp_output_dir = os.path.join(directory, "HEAT-MAP-fp-LED")
    
    if os.path.exists(std_output_dir):
        shutil.rmtree(std_output_dir)
        #print(f"Deleted directory: {std_output_dir}")
    
    if os.path.exists(fp_output_dir):
        shutil.rmtree(fp_output_dir)
        #print(f"Deleted directory: {fp_output_dir}")


def process_std_led_heatmaps(directory, params):
    """Process standard LED heatmaps and save them."""
    std_output_dir = os.path.join(directory, "HEAT-MAP-STD-LED")
    os.makedirs(std_output_dir, exist_ok=True)

    std_file_name = os.path.join(directory, "Summary_Standard_LED_matrices.xlsx")
    if os.path.exists(std_file_name):
        with pd.ExcelFile(std_file_name) as xls:
            for sheet_name in xls.sheet_names:
                save_as = os.path.join(std_output_dir, f"{sheet_name}.{params['fig_format']}")
                plot_heatmap_std_led(
                    file_name=std_file_name,
                    sheet_name=sheet_name,
                    figsize=params["figsize"],
                    vmin=params["vmin"],
                    vmax=params["vmax"],
                    fig_format=params["fig_format"],
                    set_dpi=params["set_dpi"],
                    save_as=save_as,
                    cutoff_annot=params["cutoff_annot"],
                    submatrix_coords_to_be_highlighted=params["submatrix_coords_to_be_highlighted"]
                )


def process_fp_led_heatmaps(directory, params):
    """Process fp-LED heatmaps and save them."""
    fp_output_dir = os.path.join(directory, "HEAT-MAP-fp-LED")
    os.makedirs(fp_output_dir, exist_ok=True)

    fp_file_name = os.path.join(directory, "Summary_fp-LED_matrices.xlsx")
    if os.path.exists(fp_file_name):
        with pd.ExcelFile(fp_file_name) as xls:
            for sheet_name in xls.sheet_names:
                save_as = os.path.join(fp_output_dir, f"{sheet_name}.{params['fig_format']}")
                plot_heatmap_fp_led(
                    file_name=fp_file_name,
                    sheet_name=sheet_name,
                    figsize=params["figsize"],
                    vmin=params["vmin"],
                    vmax=params["vmax"],
                    fig_format=params["fig_format"],
                    set_dpi=params["set_dpi"],
                    save_as=save_as,
                    cutoff_annot=params["cutoff_annot"],
                    submatrix_coords_to_be_highlighted=params["submatrix_coords_to_be_highlighted"]
                )


def process_fp_led_heatmaps_alt(directory, params):
    """Process standard LED heatmaps and save them."""
    std_output_dir = os.path.join(directory, "HEAT-MAP-fp-LED")
    os.makedirs(std_output_dir, exist_ok=True)

    std_file_name = os.path.join(directory, "Summary_fp-LED_matrices.xlsx")
    if os.path.exists(std_file_name):
        with pd.ExcelFile(std_file_name) as xls:
            for sheet_name in xls.sheet_names:
                save_as = os.path.join(std_output_dir, f"{sheet_name}.{params['fig_format']}")
                plot_heatmap_std_led(
                    file_name=std_file_name,
                    sheet_name=sheet_name,
                    figsize=params["figsize"],
                    vmin=params["vmin"],
                    vmax=params["vmax"],
                    fig_format=params["fig_format"],
                    set_dpi=params["set_dpi"],
                    save_as=save_as,
                    cutoff_annot=params["cutoff_annot"],
                    submatrix_coords_to_be_highlighted=params["submatrix_coords_to_be_highlighted"]
                )


def heatmap_plot_engine(app_instance, base_path, plot_params_for_std_led_matrices, plot_params_for_fp_led_matrices, show_diag_cells_for_fp_led=False, delete_existing_heatmap_directories_first=True, directory_level=0):
    """
    Main function to process all subdirectories and generate heat maps.

    Parameters:
    -----------
    base_path : str
        The base directory where the subdirectories are located.
    
    show_diag_cells_for_fp_led : bool, optional, default=False
        If True, diagonal cells with NaN values are shown on fp-LED heat maps for consistency with standard LED.
        
    delete_existing_heatmap_directories_first : bool, optional, default=True
        If True, existing directories where the plots will be saved are deleted before generating new ones.

    directory_level : int, optional, default=2
        The depth level of directories to search for the Excel files.
        0 for base path, 1 for first level, 2 for second level.
    
    plot_params_for_std_led_matrices : dict
        Dictionary containing the plotting parameters for standard LED matrices.
        
    plot_params_for_fp_led_matrices : dict
        Dictionary containing the plotting parameters for fp-LED matrices.
    """
    
    # Normalize the base path
    normalized_base_path = normalize_path(base_path)

    # Find subdirectories based on the specified directory level
    if app_instance.worker.is_cancelled:
        return

    subdirectories = find_subdirs_under_base_path(normalized_base_path, depth_level=directory_level)

    for directory in subdirectories:
        normalized_directory = normalize_path(directory)
        print(f"Saving heat maps under: {normalized_directory}")

        # Optionally delete existing heatmap directories
        if delete_existing_heatmap_directories_first: 
            delete_existing_plot_directories_first(normalized_directory)

        if app_instance.worker.is_cancelled:
            return

        # Generate and save the heatmaps for standard LED matrices
        process_std_led_heatmaps(normalized_directory, plot_params_for_std_led_matrices)
        
        if app_instance.worker.is_cancelled:
            return

        print(f"   Standard LED heat maps were saved to HEAT-MAP-STD-LED")
        
        # Process fp-LED heatmaps based on the diagonal cell option
        if show_diag_cells_for_fp_led:
            if app_instance.worker.is_cancelled:
                return

            process_fp_led_heatmaps_alt(normalized_directory, plot_params_for_fp_led_matrices)

            if app_instance.worker.is_cancelled:
                return
        else:
            if app_instance.worker.is_cancelled:
                return

            process_fp_led_heatmaps(normalized_directory, plot_params_for_fp_led_matrices)

            if app_instance.worker.is_cancelled:
                return
            
        print(f"   fp-LED heat maps were saved to HEAT-MAP-fp-LED")
    
    print('\n')
    print('*' * 100)
    print("Heat map generation job was terminated NORMALLY")
    print('*' * 100)


def heatmap_plot_engine_no_gui(base_path, plot_params_for_std_led_matrices, plot_params_for_fp_led_matrices, show_diag_cells_for_fp_led=False, delete_existing_heatmap_directories_first=True, directory_level=0):
    """
    Main function to process all subdirectories and generate heat maps.

    Parameters:
    -----------
    base_path : str
        The base directory where the subdirectories are located.
    
    show_diag_cells_for_fp_led : bool, optional, default=False
        If True, diagonal cells with NaN values are shown on fp-LED heat maps for consistency with standard LED.
        
    delete_existing_heatmap_directories_first : bool, optional, default=True
        If True, existing directories where the plots will be saved are deleted before generating new ones.

    directory_level : int, optional, default=2
        The depth level of directories to search for the Excel files.
        0 for base path, 1 for first level, 2 for second level.
    
    plot_params_for_std_led_matrices : dict
        Dictionary containing the plotting parameters for standard LED matrices.
        
    plot_params_for_fp_led_matrices : dict
        Dictionary containing the plotting parameters for fp-LED matrices.
    """
    
    # Normalize the base path
    normalized_base_path = normalize_path(base_path)

    # Find subdirectories based on the specified directory level

    subdirectories = find_subdirs_under_base_path(normalized_base_path, depth_level=directory_level)

    for directory in subdirectories:
        normalized_directory = normalize_path(directory)
        print(f"Saving heat maps under: {normalized_directory}")

        # Optionally delete existing heatmap directories
        if delete_existing_heatmap_directories_first: 
            delete_existing_plot_directories_first(normalized_directory)

        # Generate and save the heatmaps for standard LED matrices
        process_std_led_heatmaps(normalized_directory, plot_params_for_std_led_matrices)
        
        print(f"   Standard LED heat maps were saved to HEAT-MAP-STD-LED")
        
        # Process fp-LED heatmaps based on the diagonal cell option
        if show_diag_cells_for_fp_led:
            process_fp_led_heatmaps_alt(normalized_directory, plot_params_for_fp_led_matrices)

        else:
            process_fp_led_heatmaps(normalized_directory, plot_params_for_fp_led_matrices)

        print(f"   fp-LED heat maps were saved to HEAT-MAP-fp-LED")
    
    print('\n')
    print('*' * 100)
    print("Heat map generation job was terminated NORMALLY")
    print('*' * 100)

