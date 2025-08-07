import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde
import csv
import os
import matplotlib.colors as mcolors
from scipy.spatial import ConvexHull
from matplotlib.patches import Polygon


data_path = '/Users/gabrielnathan/Tools/Galilean_moons/data/'


def compute_density(xy_data, grid_size=100):
    """
    Computes a 2D density estimate over a grid.
    
    Parameters:
    xy_data: List or array of (x, y) coordinates.
    grid_size: The number of points along each axis in the density grid.
    
    Returns:
    x_grid, y_grid: Meshgrid arrays for the density plot.
    density: 2D array of density values over the grid.
    """
    # Split xy_data into x and y components
    x, y = np.array(xy_data).T
    
    # Create a meshgrid for the background density plot
    x_grid, y_grid = np.meshgrid(
        np.linspace(x.min() - 1, x.max() + 1, grid_size), 
        np.linspace(y.min() - 1, y.max() + 1, grid_size)
    )
    xy_grid = np.vstack([x_grid.ravel(), y_grid.ravel()])
    density = gaussian_kde(np.vstack([x, y]))(xy_grid).reshape(x_grid.shape)
    return x_grid, y_grid, density

def load_xy_data_from_csv(filename):
    """
    Reads a CSV file of (x, y) coordinates and converts it into a list of (x, y) tuples.
    
    Parameters:
    filename: Path to the CSV file containing x, y coordinates.
    
    Returns:
    xy_data: List of (x, y) tuples suitable for use in plotting functions.
    """
    xy_data = []
    with open(filename, 'r') as file:
        csv_reader = csv.reader(file)
        # Skip the header if present
        next(csv_reader, None)  # Skip the first row if it contains headers
        for row in csv_reader:
            try:
                # Convert each row's x and y values to float
                x, y = float(row[0]), float(row[1])
                xy_data.append((x, y))
            except ValueError:
                # Skip rows that can't be converted to float
                print(f"Skipping invalid row: {row}")
    return xy_data
    
def create_fading_cmap(density_cmap='Blues', levels=10, start_fade_level=3):
    """
    Creates a colormap that starts opaque and fades to transparency at higher levels.
    
    Parameters:
    density_cmap: Name of the base colormap to start with.
    levels: Number of contour levels.
    start_fade_level: Level at which to start fading to transparency.
    
    Returns:
    transparent_cmap: Colormap with controlled fading transparency.
    """
    # Get the colormap based on the string name
    base = plt.cm.get_cmap(density_cmap, levels)
    colors = base(np.linspace(0, 1, levels))
    
    # Make the initial levels opaque and start fading after `start_fade_level`
    fade_alpha = np.ones(levels)
    fade_alpha[start_fade_level:] = np.exp(-np.linspace(0, 5, levels - start_fade_level))
    colors[:, -1] = fade_alpha  # Apply the fade to the alpha channel
    
    transparent_cmap = mcolors.ListedColormap(colors)
    return transparent_cmap

def plot_with_density(xy_data, scatter_color, density_cmap, grid_size=100, label=None, ax=None, density_threshold=0.01):
    """
    Plots a scatter plot with a density elevation profile that fades at higher levels.
    
    Parameters:
    xy_data: List or array of (x, y) coordinates.
    scatter_color: Color of the scatter points.
    density_cmap: Color map name (string) for the density elevation profile.
    grid_size: Controls resolution of the density background.
    label: Label for the scatter points in the legend.
    ax: Matplotlib axis to plot on. Creates new axis if None.
    density_threshold: Threshold below which density values are masked out.
    """
    if ax is None:
        fig, ax = plt.subplots(figsize=(8, 6))
    
    # Separate x and y for single and multiple points
    if len(xy_data) == 1:
        x, y = [xy_data[0][0]], [xy_data[0][1]]
    else:
        x, y = np.array(xy_data).T
    
    # Plot density if more than one point
    if len(xy_data) > 1:
        # Compute density
        x_grid, y_grid, density = compute_density(xy_data, grid_size=grid_size)
        
        # Mask out very low-density values to avoid the white square artifact
        density = np.ma.masked_where(density < density_threshold, density)
        
        # Apply custom colormap that fades at higher levels
        fading_cmap = create_fading_cmap(density_cmap=density_cmap, levels=10)
        
        # Plot density contour in background with fading colormap
        ax.contourf(x_grid, y_grid, density, levels=10, cmap=fading_cmap, zorder=0)
    
    # Plot scatter points on top with higher zorder
    ax.scatter(x, y, color=scatter_color, edgecolor='black', s=50, label=label, zorder=10)
    
# def plot_this_category(filename, label, scatter_color, density_cmap, data_path="", ax=None):
    """
    Wrapper function to load a CSV file of coordinates and plot with density elevation.
    
    Parameters:
    filename: Name of the file (without .csv) containing x, y coordinates.
    label: Label for the plot legend.
    scatter_color: Color of the scatter points.
    density_cmap: Color map for the density elevation profile.
    data_path: Path to the directory containing the CSV file (optional).
    ax: Matplotlib axis to plot on. Creates new axis if None.
    """
    # Construct full path to the CSV file
    full_path = os.path.join(data_path, f"{filename}.csv")
    
    # Load the xy data from the CSV file
    xy_data = load_xy_data_from_csv(full_path)
    
    # Plot using the loaded data and specified colors, passing the label and axis
    plot_with_density(xy_data, scatter_color=scatter_color, density_cmap=density_cmap, label=label, ax=ax)

# # Example usage to plot all on the same figure:
# fig, ax = plt.subplots(figsize=(10, 8))

# # Define a brown color gradient from light tan to dark brown
# brown_colors = ["#f4e1c6", "#d2a679", "#a0522d", "#8b4513", "#5d2e0a"]
# # Create a colormap from this brown gradient
# brown_cmap = mcolors.LinearSegmentedColormap.from_list("custom_brown", brown_colors)

# # usage:

# plot_this_category('CI_oxygen', "CI", scatter_color='blue', density_cmap='Blues', data_path=data_path, ax=ax)
# plot_this_category('CV_oxygen', "CV", scatter_color='red', density_cmap='Reds', data_path=data_path, ax=ax)
# plot_this_category('ordinary_oxygen', "Ordinary", scatter_color='brown', density_cmap=brown_cmap, data_path=data_path, ax=ax)
# plot_this_category('CM_oxygen', "CM", scatter_color='green', density_cmap='Greens', data_path=data_path, ax=ax)
# plot_this_category('enstatite_oxygen', "Enstatite", scatter_color='purple', density_cmap='Purples', data_path=data_path, ax=ax)
# plot_this_category('ryugu_oxygen', "Ryugu", scatter_color='pink', density_cmap='Pinks', data_path=data_path, ax=ax)


# # Customize the plot
# plt.xlabel(u'$\delta^{17}$O, ‰')
# plt.ylabel(u'$\delta^{18}$O, ‰')
# plt.title("Triple oxygen isotopes")
# plt.legend()
# plt.show()


################################################################################################
################################################################################################

# #PLOT WITH FILLED HULL, SHARP EDGES
# def plot_this_category(filename, label, scatter_color, hull_color, data_path="", ax=None):
#     """
#     Wrapper function to load a CSV file of coordinates and plot with convex hull outline.
    
#     Parameters:
#     filename: Name of the file (without .csv) containing x, y coordinates.
#     label: Label for the plot legend.
#     scatter_color: Color of the scatter points.
#     hull_color: Color of the convex hull outline.
#     data_path: Path to the directory containing the CSV file (optional).
#     ax: Matplotlib axis to plot on. Creates new axis if None.
#     """
#     # Construct full path to the CSV file
#     full_path = os.path.join(data_path, f"{filename}.csv")
    
#     # Load the xy data from the CSV file
#     xy_data = load_xy_data_from_csv(full_path)
    
#     # Plot using the loaded data, passing the label, scatter color, and hull color
#     plot_with_filled_hull(xy_data, scatter_color=scatter_color, hull_color=hull_color, label=label, ax=ax)

# def plot_with_filled_hull(xy_data, scatter_color, hull_color, label=None, ax=None, alpha=0.3):
#     """
#     Plots a scatter plot with a filled, transparent convex hull around the set of points.
    
#     Parameters:
#     xy_data: List or array of (x, y) coordinates.
#     scatter_color: Color of the scatter points.
#     hull_color: Color of the filled convex hull.
#     label: Label for the scatter points in the legend.
#     ax: Matplotlib axis to plot on. Creates new axis if None.
#     alpha: Transparency level for the hull fill (0 to 1).
#     """
#     if ax is None:
#         fig, ax = plt.subplots(figsize=(8, 6))
    
#     # Convert xy_data to a NumPy array for easier indexing
#     xy_data = np.array(xy_data)
#     x, y = xy_data.T
    
#     # Plot the scatter points
#     ax.scatter(x, y, color=scatter_color, edgecolor='black', s=50, label=label, zorder=10)
    
#     # Plot the filled convex hull around the points
#     if len(xy_data) > 2:  # Convex hull requires at least 3 points
#         hull = ConvexHull(xy_data)
#         hull_points = xy_data[hull.vertices]  # Get the coordinates of hull vertices
        
#         # Create and add a filled polygon for the convex hull
#         polygon = Polygon(hull_points, closed=True, facecolor=hull_color, edgecolor=None, alpha=alpha, zorder=5)
#         ax.add_patch(polygon)
    

# # Example usage to plot all categories on the same figure:
# fig, ax = plt.subplots(figsize=(10, 8))

# # data_path = "path/to/your/data/"
# plot_this_category('CI_oxygen', "CI", scatter_color='blue', hull_color='lightblue', data_path=data_path, ax=ax)
# plot_this_category('CV_oxygen', "CV", scatter_color='red', hull_color='lightcoral', data_path=data_path, ax=ax)
# plot_this_category('ordinary_oxygen', "Ordinary", scatter_color='brown', hull_color='sandybrown', data_path=data_path, ax=ax)
# plot_this_category('CM_oxygen', "CM", scatter_color='green', hull_color='lightgreen', data_path=data_path, ax=ax)
# plot_this_category('enstatite_oxygen', "Enstatite", scatter_color='purple', hull_color='plum', data_path=data_path, ax=ax)
# plot_this_category('ryugu_oxygen', "Ryugu", scatter_color='pink', hull_color='lightpink', data_path=data_path, ax=ax)

# plt.xlabel(u'$\delta^{17}$O, ‰')
# plt.ylabel(u'$\delta^{18}$O, ‰')
# plt.title("Triple oxygen isotopes")
# plt.legend()
# plt.show()
################################################################################################
################################################################################################



def plot_this_category_just_scatter(filename, label, scatter_color, hull_color, data_path="", ax=None):
    """
    Wrapper function to load a CSV file of coordinates and plot with convex hull outline.
    
    Parameters:
    filename: Name of the file (without .csv) containing x, y coordinates.
    label: Label for the plot legend.
    scatter_color: Color of the scatter points.
    hull_color: Color of the convex hull outline.
    data_path: Path to the directory containing the CSV file (optional).
    ax: Matplotlib axis to plot on. Creates new axis if None.
    """
    # Construct full path to the CSV file
    full_path = os.path.join(data_path, f"{filename}.csv")
    
    # Load the xy data from the CSV file
    xy_data = load_xy_data_from_csv(full_path)
    
    # Plot using the loaded data, passing the label, scatter color, and hull color
    plot_with_no_hull(xy_data, scatter_color=scatter_color, hull_color=hull_color, label=label, ax=ax)


def plot_with_no_hull(xy_data, scatter_color, hull_color, label=None, ax=None, alpha=0.5, bandwidth=0.35):
    """
    Plots a scatter plot with a filled, transparent, smooth hull around the set of points.
    
    Parameters:
    xy_data: List or array of (x, y) coordinates.
    scatter_color: Color of the scatter points.
    hull_color: Color of the filled hull.
    label: Label for the scatter points in the legend.
    ax: Matplotlib axis to plot on. Creates new axis if None.
    alpha: Transparency level for the hull fill (0 to 1).
    bandwidth: Controls the smoothness of the hull outline.
    """
    if ax is None:
        fig, ax = plt.subplots(figsize=(8, 6))
    
    # Convert xy_data to NumPy array and split into x and y for plotting
    xy_data = np.array(xy_data)
    x, y = xy_data.T
    

    # # Generate a smooth contour around the points using KDE
    # if len(xy_data) > 1:
    #     # Define a grid for KDE
    #     x_min, x_max = x.min() - 0.5, x.max() + 0.5
    #     y_min, y_max = y.min() - 0.5, y.max() + 0.5
    #     xx, yy = np.meshgrid(np.linspace(x_min, x_max, 100), np.linspace(y_min, y_max, 100))
        
    #     # Perform KDE on the grid
    #     xy = np.vstack([x, y])
    #     kde = gaussian_kde(xy, bw_method=bandwidth)
    #     zz = kde(np.vstack([xx.ravel(), yy.ravel()])).reshape(xx.shape)
        
    #     # Plot filled contour for a smooth hull
    #     ax.contourf(xx, yy, zz, levels=[0.01, zz.max()], colors=[hull_color], alpha=alpha, zorder=10)
    # Plot the scatter points
    ax.scatter(x, y, color=scatter_color, edgecolor='black', s=70, label=label, zorder=10)
    

    # Customize the plot
    ax.set_xlabel("X-axis")
    ax.set_ylabel("Y-axis")
    if label:
        ax.legend()

def plot_this_category(filename, label, scatter_color, hull_color, data_path="", ax=None):
    """
    Wrapper function to load a CSV file of coordinates and plot with convex hull outline.
    
    Parameters:
    filename: Name of the file (without .csv) containing x, y coordinates.
    label: Label for the plot legend.
    scatter_color: Color of the scatter points.
    hull_color: Color of the convex hull outline.
    data_path: Path to the directory containing the CSV file (optional).
    ax: Matplotlib axis to plot on. Creates new axis if None.
    """
    # Construct full path to the CSV file
    full_path = os.path.join(data_path, f"{filename}.csv")
    
    # Load the xy data from the CSV file
    xy_data = load_xy_data_from_csv(full_path)
    
    # Plot using the loaded data, passing the label, scatter color, and hull color
    plot_with_smooth_hull(xy_data, scatter_color=scatter_color, hull_color=hull_color, label=label, ax=ax)


def plot_with_smooth_hull(xy_data, scatter_color, hull_color, label=None, ax=None, alpha=0.5, bandwidth=0.35):
    """
    Plots a scatter plot with a filled, transparent, smooth hull around the set of points.
    
    Parameters:
    xy_data: List or array of (x, y) coordinates.
    scatter_color: Color of the scatter points.
    hull_color: Color of the filled hull.
    label: Label for the scatter points in the legend.
    ax: Matplotlib axis to plot on. Creates new axis if None.
    alpha: Transparency level for the hull fill (0 to 1).
    bandwidth: Controls the smoothness of the hull outline.
    """
    if ax is None:
        fig, ax = plt.subplots(figsize=(8, 6))
    
    # Convert xy_data to NumPy array and split into x and y for plotting
    xy_data = np.array(xy_data)
    x, y = xy_data.T
    

    # Generate a smooth contour around the points using KDE
    if len(xy_data) > 1:
        # Define a grid for KDE
        x_min, x_max = x.min() - 0.5, x.max() + 0.5
        y_min, y_max = y.min() - 0.5, y.max() + 0.5
        xx, yy = np.meshgrid(np.linspace(x_min, x_max, 100), np.linspace(y_min, y_max, 100))
        
        # Perform KDE on the grid
        xy = np.vstack([x, y])
        kde = gaussian_kde(xy, bw_method=bandwidth)
        zz = kde(np.vstack([xx.ravel(), yy.ravel()])).reshape(xx.shape)
        
        # Plot filled contour for a smooth hull
        ax.contourf(xx, yy, zz, levels=[0.01, zz.max()], colors=[hull_color], alpha=alpha, zorder=10)
    # Plot the scatter points
    ax.scatter(x, y, color=scatter_color, edgecolor='black', s=50, label=label, zorder=10)
    

    # Customize the plot
    # ax.set_xlabel(fontsize=16)
	# ax.set_ylabel(fontsize=16)
    # ax.set_title("Plot Title", fontsize=18)
    ax.legend(fontsize=14)
    if label:
        ax.legend(fontsize=14)


def plot_line(m, b, x_range=(-10, 10), line_color='blue', line_style='-', line_label=None, ax=None):
    """
    Plots a line with the equation y = mx + b.
    
    Parameters:
    m : float
        Slope of the line.
    b : float
        Intercept of the line.
    x_range : tuple, optional
        Range of x values as (xmin, xmax). Default is (-10, 10).
    line_color : str, optional
        Color of the line. Default is 'blue'.
    line_label : str, optional
        Label for the line, useful for legend. Default is None.
    """
    # Generate x values and corresponding y values
    x_values = np.linspace(x_range[0], x_range[1], 100)
    y_values = m * x_values + b
    
    # Plot the line
    ax.plot(x_values, y_values, color=line_color, label=line_label, linestyle=line_style, linewidth=2.0)
    # plt.xlabel("X-axis")
    # plt.ylabel("Y-axis")
    # plt.title(f"Plot of Line: y = {m}x + {b}")
    
    # Show legend if label is provided
    if line_label:
        plt.legend()
    





plt.rcParams.update({'font.size': 14})  # Adjust number as needed (e.g., 14 for larger font)
# Example usage to plot all categories on the same figure:
fig, ax = plt.subplots(figsize=(10, 8))

# data_path = "path/to/your/data/"
# plot_this_category('CI_oxygen', "CI", scatter_color='blue', hull_color='lightblue', data_path=data_path, ax=ax)
# plot_this_category('CV_oxygen', "CV", scatter_color='red', hull_color='lightcoral', data_path=data_path, ax=ax)
# plot_this_category('ordinary_oxygen', "Ordinary", scatter_color='brown', hull_color='sandybrown', data_path=data_path, ax=ax)
# plot_this_category('CM_oxygen', "CM", scatter_color='green', hull_color='lightgreen', data_path=data_path, ax=ax)
# plot_this_category('enstatite_oxygen', "Enstatite", scatter_color='purple', hull_color='plum', data_path=data_path, ax=ax)
# plot_this_category('ryugu_oxygen', "Ryugu", scatter_color='pink', hull_color='lightpink', data_path=data_path, ax=ax)

# plt.xlabel(u'$\delta^{17}$O, ‰')
# plt.ylabel(u'$\delta^{18}$O, ‰')
# plt.title("Triple oxygen isotopes")
# plt.legend()
# plt.show()
def plot_point_with_ranges(ax, x, y, x_range, y_range, point_color='blue', range_color='blue', point_label=None, alpha=0.3):
    """
    Plots a single point with an open circle, a horizontal range, and a vertical range with transparency on a given axis.
    
    Parameters:
    ax : matplotlib.axes.Axes
        Axis on which to plot the point and ranges.
    x : float
        X-coordinate of the point.
    y : float
        Y-coordinate of the point.
    x_range : tuple
        Range of x values as (xmin, xmax) for the horizontal bar.
    y_range : tuple
        Range of y values as (ymin, ymax) for the vertical bar.
    point_color : str, optional
        Color of the point. Default is 'blue'.
    range_color : str, optional
        Color of the range bars. Default is 'blue'.
    point_label : str, optional
        Label for the point, useful for legend. Default is None.
    alpha : float, optional
        Transparency level for the range bars (0 to 1). Default is 0.3.
    """
    # Plot the horizontal range with transparency
    ax.fill_betweenx([y - 0.1, y + 0.1], x_range[0], x_range[1], color=range_color, alpha=alpha)

    # Plot the vertical range with transparency
    ax.fill_between([x - 0.1, x + 0.1], y_range[0], y_range[1], color=range_color, alpha=alpha)

    # Plot the point with an open circle
    ax.errorbar(x, y, xerr=0, yerr=0, fmt='o', color=point_color, markerfacecolor='none', markersize=15, label=point_label, capsize=15)
    
    # # Labeling and grid
    # ax.set_xlabel("X-axis")
    # ax.set_ylabel("Y-axis")
    # ax.set_title("Point with Horizontal and Vertical Ranges")
    # ax.grid(True)
    
    # Show legend if label is provided
    if point_label:
        ax.legend()

# Example usage
# fig, ax = plt.subplots(figsize=(8, 6))
# plot_point_with_ranges(ax, x=2.3, y=2.3, x_range=(0.99, 7.55), y_range=(0.99, 7.55), point_color='purple', range_color='purple', point_label="Europa (Bierson & Nimmo)", alpha=0.2)
plot_point_with_ranges(ax, x=2.3, y=2.3, x_range=(2.2, 2.4), y_range=(2.2, 2.4), point_color='purple', range_color='purple', point_label="Europa (Bierson & Nimmo)", alpha=0.2)
plot_point_with_ranges(ax, x=0.5, y=0.5, x_range=(0.06, 0.39), y_range=(0.06, 0.39), point_color='orange', range_color='orange', point_label="Ganymede (Bierson & Nimmo)", alpha=0.2)
plot_point_with_ranges(ax, x=0.08, y=0.08, x_range=(0.06,0.1), y_range=(0.06,0.1), point_color='green', range_color='green', point_label="Callisto (Bierson & Nimmo)", alpha=0.2)

# plt.show()

# data_path = "path/to/your/data/"
plot_this_category_just_scatter('CI_oxygen', "CI", scatter_color='#1C01F0', hull_color='lightblue', data_path=data_path, ax=ax)
plot_this_category_just_scatter('CV_oxygen', "CV", scatter_color='#F00103', hull_color='lightcoral', data_path=data_path, ax=ax)
plot_this_category_just_scatter('ordinary_oxygen', "Ordinary", scatter_color='#9C9340', hull_color='sandybrown', data_path=data_path, ax=ax)
plot_this_category_just_scatter('CM_oxygen', "CM", scatter_color='#3ECF2B', hull_color='lightgreen', data_path=data_path, ax=ax)
plot_this_category_just_scatter('enstatite_oxygen', "Enstatite", scatter_color='#9C00C1', hull_color='plum', data_path=data_path, ax=ax)
plot_this_category_just_scatter('ryugu_oxygen', "Ryugu", scatter_color='#DB14D6', hull_color='lightpink', data_path=data_path, ax=ax)


ax.scatter(0, 0, color='blue', s=100, marker='s', label="Earth (SMOW)")
ax.scatter(3.2, 4.1, color='red', s=100, marker='s', label="Mars (SNC)")
ax.scatter(20.5, 17.7, color='teal', s=100, label="Asteroidal H2O")


# plot lines
plot_line(m=0.52, b=0, x_range=(-10, 30), line_color='black', line_style='solid', line_label='TFL', ax=ax)

plot_line(m=0.94, b=-4.47, x_range=(-10, 30), line_color='grey', line_style='solid', line_label='CCAM', ax=ax)

# plot_line(m=1.0, b=0, x_range=(-10, 30), line_color='grey', line_style='dotted', line_label='TFL', ax=ax)
    
# Arrow parameters
start_x, start_y = 12, 12
slope = 0.94
arrow_length = 10
end_x = start_x + arrow_length
end_y = start_y + arrow_length * slope

# Draw a thicker arrow
ax.annotate(
    '', xy=(end_x, end_y), xytext=(start_x, start_y),
    arrowprops=dict(facecolor='teal', edgecolor='teal', linewidth=5, arrowstyle='->')
)

# Add black, bold, and rotated text parallel to the arrow
text_x = (start_x + end_x) / 2  # midpoint x for positioning the text
text_y = (start_y + end_y) / 2 + 1  # slightly above the midpoint
rotation_angle = np.degrees(np.arctan(slope))

ax.text(
    text_x, text_y, "to cometary ices", color='black', fontsize=12, fontweight='bold',
    rotation=rotation_angle, rotation_mode='anchor'
)

plt.ylabel(u'$\delta^{17}$O, ‰', fontsize=20)
plt.xlabel(u'$\delta^{18}$O, ‰', fontsize=20)

plt.xlim(-10, 25)
plt.ylim(-10, 25)
# plt.title("Triple oxygen isotopes")
plt.legend(loc='upper left',fontsize=12)
plt.show()







