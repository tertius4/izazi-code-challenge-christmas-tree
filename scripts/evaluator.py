#!/usr/bin/env python3
"""
Christmas Tree Packing Challenge - Evaluator

This script validates and scores submissions for the Christmas Tree Packing Challenge.
It checks for:
- Correct submission format
- No overlapping trees
- Valid coordinate constraints

Usage:
    python evaluator.py submission.csv
    python evaluator.py submission.csv --visualize 10  # Visualize the 10-tree configuration
"""

import argparse
import sys
from decimal import Decimal, getcontext

import pandas as pd
from shapely import affinity
from shapely.geometry import Polygon
from shapely.ops import unary_union
from shapely.strtree import STRtree

# Set precision for Decimal calculations
getcontext().prec = 25
SCALE_FACTOR = Decimal('1e15')

# Coordinate constraints
MIN_COORD = -100
MAX_COORD = 100


def create_tree_polygon(center_x, center_y, angle):
    """
    Creates a Christmas tree polygon at the specified position and rotation.

    The tree shape consists of:
    - A trunk at the bottom (0.15 wide, 0.2 tall)
    - Three tiers of branches (base 0.7 wide, middle 0.4 wide, top 0.25 wide)
    - A tip at the top (y = 0.8)

    Args:
        center_x: X coordinate of tree center (as Decimal or string)
        center_y: Y coordinate of tree center (as Decimal or string)
        angle: Rotation angle in degrees (as Decimal or string)

    Returns:
        A Shapely Polygon representing the tree
    """
    center_x = Decimal(str(center_x))
    center_y = Decimal(str(center_y))
    angle = Decimal(str(angle))

    # Tree dimensions
    trunk_w = Decimal('0.15')
    trunk_h = Decimal('0.2')
    base_w = Decimal('0.7')
    mid_w = Decimal('0.4')
    top_w = Decimal('0.25')

    # Y coordinates for each tier
    tip_y = Decimal('0.8')
    tier_1_y = Decimal('0.5')
    tier_2_y = Decimal('0.25')
    base_y = Decimal('0.0')
    trunk_bottom_y = -trunk_h

    # Create the tree polygon vertices (starting from tip, going clockwise)
    initial_polygon = Polygon([
        # Tip
        (float(Decimal('0.0') * SCALE_FACTOR), float(tip_y * SCALE_FACTOR)),
        # Right side - Top Tier
        (float(top_w / Decimal('2') * SCALE_FACTOR), float(tier_1_y * SCALE_FACTOR)),
        (float(top_w / Decimal('4') * SCALE_FACTOR), float(tier_1_y * SCALE_FACTOR)),
        # Right side - Middle Tier
        (float(mid_w / Decimal('2') * SCALE_FACTOR), float(tier_2_y * SCALE_FACTOR)),
        (float(mid_w / Decimal('4') * SCALE_FACTOR), float(tier_2_y * SCALE_FACTOR)),
        # Right side - Bottom Tier
        (float(base_w / Decimal('2') * SCALE_FACTOR), float(base_y * SCALE_FACTOR)),
        # Right Trunk
        (float(trunk_w / Decimal('2') * SCALE_FACTOR), float(base_y * SCALE_FACTOR)),
        (float(trunk_w / Decimal('2') * SCALE_FACTOR), float(trunk_bottom_y * SCALE_FACTOR)),
        # Left Trunk
        (float(-(trunk_w / Decimal('2')) * SCALE_FACTOR), float(trunk_bottom_y * SCALE_FACTOR)),
        (float(-(trunk_w / Decimal('2')) * SCALE_FACTOR), float(base_y * SCALE_FACTOR)),
        # Left side - Bottom Tier
        (float(-(base_w / Decimal('2')) * SCALE_FACTOR), float(base_y * SCALE_FACTOR)),
        # Left side - Middle Tier
        (float(-(mid_w / Decimal('4')) * SCALE_FACTOR), float(tier_2_y * SCALE_FACTOR)),
        (float(-(mid_w / Decimal('2')) * SCALE_FACTOR), float(tier_2_y * SCALE_FACTOR)),
        # Left side - Top Tier
        (float(-(top_w / Decimal('4')) * SCALE_FACTOR), float(tier_1_y * SCALE_FACTOR)),
        (float(-(top_w / Decimal('2')) * SCALE_FACTOR), float(tier_1_y * SCALE_FACTOR)),
    ])

    # Rotate around origin, then translate to position
    rotated = affinity.rotate(initial_polygon, float(angle), origin=(0, 0))
    translated = affinity.translate(
        rotated,
        xoff=float(center_x * SCALE_FACTOR),
        yoff=float(center_y * SCALE_FACTOR)
    )
    return translated


def check_overlap(polygons):
    """
    Check if any polygons overlap (intersect but don't just touch).

    Args:
        polygons: List of Shapely polygons

    Returns:
        Tuple of (has_overlap: bool, overlapping_pairs: list of (i, j) tuples)
    """
    if len(polygons) <= 1:
        return False, []

    tree_index = STRtree(polygons)
    overlapping_pairs = []

    for i, poly in enumerate(polygons):
        possible_indices = tree_index.query(poly)
        for j in possible_indices:
            if j <= i:
                continue  # Skip self and already-checked pairs
            if polygons[i].intersects(polygons[j]) and not polygons[i].touches(polygons[j]):
                overlapping_pairs.append((i, j))

    return len(overlapping_pairs) > 0, overlapping_pairs


def calculate_bounding_box_side(polygons):
    """
    Calculate the side length of the smallest square bounding box for a set of polygons.

    Args:
        polygons: List of Shapely polygons

    Returns:
        Side length of the square bounding box (as Decimal)
    """
    if not polygons:
        return Decimal('0')

    union = unary_union(polygons)
    bounds = union.bounds  # (minx, miny, maxx, maxy)

    minx = Decimal(str(bounds[0])) / SCALE_FACTOR
    miny = Decimal(str(bounds[1])) / SCALE_FACTOR
    maxx = Decimal(str(bounds[2])) / SCALE_FACTOR
    maxy = Decimal(str(bounds[3])) / SCALE_FACTOR

    width = maxx - minx
    height = maxy - miny

    return max(width, height)


def parse_submission_value(value):
    """
    Parse a submission value (prefixed with 's') to a float.

    Args:
        value: String like 's0.123456' or 's-1.5'

    Returns:
        Float value
    """
    if isinstance(value, str) and value.startswith('s'):
        return float(value[1:])
    return float(value)


def load_submission(filepath):
    """
    Load and parse a submission CSV file.

    Args:
        filepath: Path to the submission CSV file

    Returns:
        DataFrame with columns: n (number of trees), tree_idx, x, y, deg
    """
    df = pd.read_csv(filepath, index_col='id')

    # Parse the id column to extract n and tree index
    records = []
    for idx, row in df.iterrows():
        parts = idx.split('_')
        n = int(parts[0])
        tree_idx = int(parts[1])
        x = parse_submission_value(row['x'])
        y = parse_submission_value(row['y'])
        deg = parse_submission_value(row['deg'])
        records.append({'n': n, 'tree_idx': tree_idx, 'x': x, 'y': y, 'deg': deg})

    return pd.DataFrame(records)


def validate_submission(df):
    """
    Validate submission format and constraints.

    Args:
        df: DataFrame from load_submission

    Returns:
        Tuple of (is_valid: bool, errors: list of error messages)
    """
    errors = []

    # Check we have all configurations from 1 to 200
    expected_ids = set()
    for n in range(1, 201):
        for t in range(n):
            expected_ids.add((n, t))

    actual_ids = set(zip(df['n'], df['tree_idx']))

    missing = expected_ids - actual_ids
    if missing:
        errors.append(f"Missing {len(missing)} tree entries. Examples: {list(missing)[:5]}")

    extra = actual_ids - expected_ids
    if extra:
        errors.append(f"Extra {len(extra)} tree entries. Examples: {list(extra)[:5]}")

    # Check coordinate constraints
    out_of_bounds = df[(df['x'] < MIN_COORD) | (df['x'] > MAX_COORD) |
                       (df['y'] < MIN_COORD) | (df['y'] > MAX_COORD)]
    if len(out_of_bounds) > 0:
        errors.append(f"{len(out_of_bounds)} trees have coordinates outside [{MIN_COORD}, {MAX_COORD}]")

    return len(errors) == 0, errors


def evaluate_submission(df, verbose=True):
    """
    Evaluate a submission and calculate the total score.

    Score = sum over all n from 1 to 200 of (s_n^2 / n)
    where s_n is the side length of the square bounding box for the n-tree configuration.

    Args:
        df: DataFrame from load_submission
        verbose: If True, print progress

    Returns:
        Tuple of (total_score: float, per_config_scores: dict, overlap_errors: list)
    """
    total_score = Decimal('0')
    per_config_scores = {}
    overlap_errors = []

    for n in range(1, 201):
        config_df = df[df['n'] == n]

        # Create tree polygons for this configuration
        polygons = []
        for _, row in config_df.iterrows():
            poly = create_tree_polygon(row['x'], row['y'], row['deg'])
            polygons.append(poly)

        # Check for overlaps
        has_overlap, pairs = check_overlap(polygons)
        if has_overlap:
            overlap_errors.append(f"Configuration n={n}: {len(pairs)} overlapping pairs")

        # Calculate bounding box and score
        side = calculate_bounding_box_side(polygons)
        config_score = (side * side) / Decimal(n)
        per_config_scores[n] = float(config_score)
        total_score += config_score

        if verbose and n % 20 == 0:
            print(f"  Evaluated configurations 1-{n}...")

    return float(total_score), per_config_scores, overlap_errors


def visualize_configuration(df, n):
    """
    Visualize a specific n-tree configuration.

    Args:
        df: DataFrame from load_submission
        n: Number of trees to visualize
    """
    try:
        import matplotlib.pyplot as plt
        from matplotlib.patches import Rectangle
    except ImportError:
        print("matplotlib is required for visualization. Install with: pip install matplotlib")
        return

    config_df = df[df['n'] == n]
    if len(config_df) == 0:
        print(f"No data for configuration n={n}")
        return

    # Create polygons
    polygons = []
    for _, row in config_df.iterrows():
        poly = create_tree_polygon(row['x'], row['y'], row['deg'])
        polygons.append(poly)

    # Calculate bounding box
    side = calculate_bounding_box_side(polygons)
    union = unary_union(polygons)
    bounds = union.bounds

    minx = Decimal(str(bounds[0])) / SCALE_FACTOR
    miny = Decimal(str(bounds[1])) / SCALE_FACTOR
    maxx = Decimal(str(bounds[2])) / SCALE_FACTOR
    maxy = Decimal(str(bounds[3])) / SCALE_FACTOR

    width = maxx - minx
    height = maxy - miny

    # Plot
    fig, ax = plt.subplots(figsize=(8, 8))
    colors = plt.cm.viridis([i / n for i in range(n)])

    for i, poly in enumerate(polygons):
        x_scaled, y_scaled = poly.exterior.xy
        x = [Decimal(str(val)) / SCALE_FACTOR for val in x_scaled]
        y = [Decimal(str(val)) / SCALE_FACTOR for val in y_scaled]
        ax.plot([float(xi) for xi in x], [float(yi) for yi in y], color=colors[i])
        ax.fill([float(xi) for xi in x], [float(yi) for yi in y], alpha=0.5, color=colors[i])

    # Draw bounding square
    square_x = minx if width >= height else minx - (side - width) / 2
    square_y = miny if height >= width else miny - (side - height) / 2

    bounding_square = Rectangle(
        (float(square_x), float(square_y)),
        float(side),
        float(side),
        fill=False,
        edgecolor='red',
        linewidth=2,
        linestyle='--',
    )
    ax.add_patch(bounding_square)

    padding = 0.5
    ax.set_xlim(float(square_x - Decimal(str(padding))),
                float(square_x + side + Decimal(str(padding))))
    ax.set_ylim(float(square_y - Decimal(str(padding))),
                float(square_y + side + Decimal(str(padding))))
    ax.set_aspect('equal', adjustable='box')
    ax.set_title(f'{n} Trees - Side Length: {float(side):.6f}')
    plt.show()


def main():
    parser = argparse.ArgumentParser(
        description='Evaluate a Christmas Tree Packing Challenge submission'
    )
    parser.add_argument('submission', help='Path to submission CSV file')
    parser.add_argument('--visualize', type=int, metavar='N',
                        help='Visualize the N-tree configuration')
    parser.add_argument('--quiet', action='store_true',
                        help='Suppress progress output')

    args = parser.parse_args()

    print(f"Loading submission from {args.submission}...")
    try:
        df = load_submission(args.submission)
    except Exception as e:
        print(f"Error loading submission: {e}")
        sys.exit(1)

    print("Validating submission format...")
    is_valid, errors = validate_submission(df)
    if not is_valid:
        print("Validation errors:")
        for error in errors:
            print(f"  - {error}")
        sys.exit(1)
    print("Format validation passed.")

    if args.visualize:
        visualize_configuration(df, args.visualize)
        return

    print("\nEvaluating submission...")
    total_score, per_config, overlap_errors = evaluate_submission(df, verbose=not args.quiet)

    if overlap_errors:
        print("\nOVERLAP ERRORS DETECTED:")
        for error in overlap_errors[:10]:  # Show first 10
            print(f"  - {error}")
        if len(overlap_errors) > 10:
            print(f"  ... and {len(overlap_errors) - 10} more")
        print("\nSubmission would be REJECTED due to overlapping trees.")
        sys.exit(1)

    print(f"\n{'='*50}")
    print(f"TOTAL SCORE: {total_score:.6f}")
    print(f"{'='*50}")
    print("\nLower scores are better.")


if __name__ == '__main__':
    main()
