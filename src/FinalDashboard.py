#!/usr/bin/env python
# coding: utf-8

# In[5]:


"""
=====================================================================================
 MAHONING COUNTY / YOUNGSTOWN FOOD ACCESS & VULNERABILITY INFORMATION SYSTEM
=====================================================================================

WHAT THIS PROGRAM DOES 
-----------------------------------------
Some neighborhoods have far fewer places to buy fresh, affordable food than others,
and some residents have no car or bus to reach the places that do exist. This program:

    1. Reads five data files (neighborhood boundaries, census statistics, food stores,
       free-food resources such as pantries, and the WRTA bus schedule).
    2. Measures, for every neighborhood, how far people must travel to reach a real
       grocery store, a bus stop, and a food pantry.
    3. Combines those distances with income, poverty and car ownership into one
       "vulnerability score" from 0 (least vulnerable) to 100 (most vulnerable).
    4. Builds an interactive web map that planners, residents and food pantries can click.

THE DATA FILES (all downloaded earlier; the program finds them automatically)
-----------------------------------------------------------------------------
    mah_census_tracts.shp   Boundaries of the 75 census tracts in Mahoning County
    acs_data.csv            Census (ACS) statistics for each tract: population, income ...
    grocery_stores.csv      Grocery, convenience and other food stores with coordinates
    resources.csv           Food pantries, free meals, food giveaways, with coordinates
    wrta_gfts.zip(.zip)     WRTA bus data in GTFS format (stops, routes, schedules)

HOW THE PROGRAM IS ORGANIZED 
---------------------------------------------------------------
    STEP 1   Load census tracts ............ Data acquisition and integration  (R1)
    STEP 2   Load census (ACS) statistics .. Data acquisition and integration
    STEP 3   Load food stores .............. Data acquisition and integration
    STEP 4   Load community food resources . Data acquisition and integration  (R5)
    STEP 5   Load WRTA bus data ............ Data acquisition and integration
    STEP 6   Save cleaned data ............. Spatial database
    STEP 7   Measure food access ........... Food access analysis              (R3)
    STEP 8   Compute vulnerability score ... Vulnerability assessment          (R4)
    STEP 9   Save results .................. Spatial database
    STEP 10  Build interactive map ......... Interactive GIS dashboard      (R1, R2, R5)

SHORT GLOSSARY OF TERMS USED IN THE COMMENTS
--------------------------------------------
    Census tract ... a small statistical neighborhood (about 1,200-8,000 people).
    ACS ............ American Community Survey: the Census Bureau's yearly statistics.
    GTFS ........... the standard file format transit agencies use to publish schedules.
    DataFrame ...... a table (rows and columns), like an Excel sheet inside Python.
    GeoDataFrame ... a DataFrame where each row also has a shape on a map
                     (a point, a line or a polygon).
    CRS ............ "coordinate reference system": the way positions are written down.
                     Latitude/longitude (EPSG:4326) is good for web maps but measures in
                     degrees. EPSG:3734 (Ohio State Plane) measures in FEET, so we use it
                     whenever we need real distances.
    Spatial join ... attaching information to a point based on which shape it falls in
                     (e.g. "which tract is this store in?").
    Buffer ......... a circle drawn around a point at a given radius.
    NaN ............ "not a number": how Python marks a missing value.
"""

# -------------------------------------------------------------------------------------
# IMPORTS: the toolboxes this program borrows from
# -------------------------------------------------------------------------------------
import os                      # read environment settings (used for the data folder)
import re                      # "regular expressions": pattern matching inside text
import sys                     # lets us stop the program with a friendly message
import zipfile                 # opens .zip files (the bus data comes zipped)
from pathlib import Path       # a safe way to work with file and folder paths

import numpy as np             # fast math on lists of numbers
import pandas as pd            # tables of data (DataFrames)

try:
    import geopandas as gpd                    # tables of data that also have map shapes
    from shapely.geometry import LineString    # lets us draw a bus route as a line
    import folium                              # builds the interactive web map
    from folium import plugins                 # extras for the map (e.g. a search box)
except ImportError as err:
    # If a package is not installed, explain exactly what to type, then stop.
    sys.exit(f"Missing package '{err.name}'. Install everything with:\n"
             "  pip install geopandas pandas numpy shapely pyproj folium")


# =====================================================================================
# CONFIGURATION
# =====================================================================================

# Folder that holds your downloaded files. The program also searches inside sub-folders,
# so it does not matter if the files sit in "Peter" or in "Peter\ForPython".
DATA_DIR = Path(os.environ.get("FOOD_DATA_DIR", r"C:\Users\Kibow\Documents\Youngstown_Food_Access"))
OUT_DIR = DATA_DIR / "output3"            # all results are written to this folder

# Coordinate systems (see glossary): one for web maps, one for measuring distance in feet.
CRS_WGS84 = "EPSG:4326"
CRS_PROJ = "EPSG:3734"
FT_PER_MILE = 5280

# Any point outside this rectangle (roughly Mahoning County plus a margin) is treated as
# a data-entry mistake. This is how a coordinate of 0,0 gets caught.
BBOX = {"lat": (40.6, 41.5), "lon": (-81.3, -80.2)}

# Which values of "Store_Type" count as a FULL-SERVICE grocery store (lowercase).
# Convenience stores, dollar stores, bakeries and farmers markets are NOT included
# because they usually lack a full range of fresh food. Edit this after reviewing the
# list the program prints in STEP 3.
FULL_SERVICE_TYPES = {"supermarket", "super store", "grocery store"}

TRANSIT_WALK_MI = 0.25     # a bus stop within a quarter mile counts as "walkable"
RADIUS_MI = 1.0            # "within reach" means within one mile

# --- Phase 2 unresolved design decision #1: how to measure distance -------------------
# This program uses STRAIGHT-LINE distance from each tract's center to the nearest place.
# (Distance along real roads would need extra road data; see the note in STEP 7.)

# --- Phase 2 unresolved design decision #2: the vulnerability formula -----------------
# Each ingredient is scaled 0-1 (1 = worst) and multiplied by its weight below.
# The weights must add up to 1.0. Change them and re-run to test other formulas.
WEIGHTS = {
    "dist_fullservice_mi": 0.30,   # farther from a full-service grocery  = more vulnerable
    "no_vehicle_rate": 0.25,       # more households without a car       = more vulnerable
    "low_income": 0.15,            # lower median household income       = more vulnerable
    "poverty_rate": 0.15,          # more people below the poverty line  = more vulnerable
    "dist_stop_mi": 0.15,          # farther from a WRTA bus stop        = more vulnerable
}
assert abs(sum(WEIGHTS.values()) - 1.0) < 1e-9, "WEIGHTS must add up to 1.0"

LOW_POP_THRESHOLD = 500    # tracts with fewer people get a "use caution" note
CATEGORY_LABELS = ["Low", "Moderate", "High", "Very High"]     # four equal-size groups
CATEGORY_COLORS = {"Low": "#2ca25f", "Moderate": "#fee08b",    # map colors for each group
                   "High": "#fc8d59", "Very High": "#d73027"}


# =====================================================================================
# SMALL HELPER FUNCTIONS used by several steps
# =====================================================================================
def banner(title, explanation=""):
    """Print a clear heading (and a one-line explanation) so you can follow progress."""
    print("\n" + "=" * 78)
    print(title)
    if explanation:
        print("  " + explanation)
    print("=" * 78)


def find_file(*names):
    """
    Look for a data file inside DATA_DIR (and all of its sub-folders).
    Pass several possible names; the first one found is used.
    """
    if not DATA_DIR.exists():
        raise FileNotFoundError(f"DATA_DIR does not exist: {DATA_DIR}\n"
                                "-> Edit DATA_DIR in the CONFIGURATION section.")
    for name in names:
        hits = [path for path in sorted(DATA_DIR.rglob(name))
                if path.is_file()
                and not path.name.endswith(".lock")      # ignore ArcGIS temporary lock files
                and OUT_DIR not in path.parents]         # ignore our own output folder
        if hits:
            return hits[0]
    raise FileNotFoundError(f"Could not find any of {names} under {DATA_DIR}")


def read_csv_any(path, **kwargs):
    """
    Read a CSV file. Files saved by Excel often use a Windows text encoding instead of
    UTF-8 (resources.csv does), so we try several encodings until one works.
    """
    for encoding in ("utf-8-sig", "cp1252", "latin-1"):
        try:
            table = pd.read_csv(path, encoding=encoding, **kwargs)
            table.columns = [c.strip() for c in table.columns]   # remove stray spaces
            return table
        except UnicodeDecodeError:
            continue                                             # try the next encoding
    raise UnicodeDecodeError("csv", b"", 0, 1, f"could not decode {path}")


def to_num(series):
    """
    Turn a column of text into numbers.
    Commas, $ and % signs are removed. Anything that is not a number (such as the "-"
    the Census uses for 'not available') becomes NaN, and so do negative values, because
    the Census uses negative codes like -666666666 to mean 'no data'.
    """
    cleaned = series.astype(str).str.replace(r"[,$%]", "", regex=True).str.strip()
    numbers = pd.to_numeric(cleaned, errors="coerce")     # "coerce" = bad text -> NaN
    return numbers.mask(numbers < 0)


def std_address(value):
    """
    Put an address in a standard form ("123 North Main Street" -> "123 N MAIN ST") so
    the same address typed two different ways can be recognized as a duplicate.
    """
    if pd.isna(value):
        return ""
    text = re.sub(r"[.,]", " ", str(value).upper())
    abbreviations = {r"\bSTREET\b": "ST", r"\bAVENUE\b": "AVE", r"\bROAD\b": "RD",
                     r"\bBOULEVARD\b": "BLVD", r"\bDRIVE\b": "DR", r"\bLANE\b": "LN",
                     r"\bNORTH\b": "N", r"\bSOUTH\b": "S", r"\bEAST\b": "E", r"\bWEST\b": "W"}
    for pattern, short in abbreviations.items():
        text = re.sub(pattern, short, text)
    return re.sub(r"\s+", " ", text).strip()              # collapse repeated spaces


def valid_coords(lat, lon):
    """True for rows whose latitude AND longitude fall inside the BBOX rectangle."""
    inside = lat.between(*BBOX["lat"]) & lon.between(*BBOX["lon"])
    return inside.fillna(False)                           # missing coordinates = invalid


def make_points(table, lon_col="lon", lat_col="lat"):
    """
    Turn a table that has longitude/latitude columns into a map layer of points,
    then convert it to the feet-based CRS so we can measure distances.
    """
    layer = gpd.GeoDataFrame(table.copy(),
                             geometry=gpd.points_from_xy(table[lon_col], table[lat_col]),
                             crs=CRS_WGS84)
    return layer.to_crs(CRS_PROJ)


def tag_with_tract(points, tracts):
    """
    SPATIAL JOIN: find out which census tract each point falls inside.
    Adds two columns:  Tract_ID      (the tract, or empty if outside the county)
                       in_study_area (True if the point is inside Mahoning County)
    """
    joined = gpd.sjoin(points, tracts[["Tract_ID", "geometry"]],
                       how="left", predicate="within")    # "left" keeps points outside too
    joined = joined.drop(columns=[c for c in ("index_right",) if c in joined.columns])
    joined = joined[~joined.index.duplicated(keep="first")]  # a point exactly on a border
    joined["in_study_area"] = joined["Tract_ID"].notna()     # could match two tracts
    return joined


# =====================================================================================
# STEP 1: CENSUS TRACTS  (the neighborhood boundaries everything else is attached to)
# =====================================================================================
def load_tracts():
    """Read the tract shapefile and return a clean table with one row per tract."""
    banner("STEP 1: Census tracts",
           "Loading neighborhood boundaries and preparing them for distance measuring.")
    path = find_file("mah_census_tracts.shp", "census_tracts.shp", "*tracts*.shp")
    tracts = gpd.read_file(path)
    print(f"Loaded {len(tracts)} tracts from {path.name} (coordinate system: {tracts.crs})")
    if tracts.crs is None:
        raise ValueError("The shapefile has no .prj file, so its coordinate system is unknown.")

    # The official 11-digit tract code (state + county + tract) becomes our unique ID.
    id_column = "GEOID" if "GEOID" in tracts.columns else "tract_id"
    tracts["Tract_ID"] = tracts[id_column].astype(str).str.extract(r"(\d{11})$")[0]
    tracts = tracts.dropna(subset=["Tract_ID"]).drop_duplicates("Tract_ID")
    tracts = tracts[tracts.geometry.notna() & ~tracts.geometry.is_empty].copy()

    # The Census Bureau publishes an "internal point" (a center) for every tract.
    # Distances are measured starting from this point.
    if {"INTPTLAT", "INTPTLON"} <= set(tracts.columns):
        tracts["int_lat"] = tracts["INTPTLAT"].astype(float)
        tracts["int_lon"] = tracts["INTPTLON"].astype(float)
    else:                                                  # fallback if columns are missing
        centers = tracts.geometry.representative_point().to_crs(CRS_WGS84)
        tracts["int_lat"], tracts["int_lon"] = centers.y.values, centers.x.values

    tracts["Tract_Name"] = (tracts["NAMELSAD"] if "NAMELSAD" in tracts.columns
                            else "Tract " + tracts["Tract_ID"])        # e.g. "Census Tract 8011"
    tracts = tracts.to_crs(CRS_PROJ)                       # switch to feet-based coordinates
    tracts["geometry"] = tracts.geometry.buffer(0)         # repairs any broken polygon edges
    tracts["area_sqmi"] = tracts.geometry.area / FT_PER_MILE ** 2    # square feet -> sq. miles
    keep = ["Tract_ID", "Tract_Name", "area_sqmi", "int_lat", "int_lon", "geometry"]
    return tracts[keep].reset_index(drop=True)


# =====================================================================================
# STEP 2: ACS DATA  (population, income, poverty, car ownership for each tract)
# =====================================================================================
def load_acs():
    """Read acs_data.csv and return clean numbers keyed by Tract_ID."""
    banner("STEP 2: Census (ACS) statistics",
           "Loading income, poverty and vehicle statistics for each tract.")
    raw = read_csv_any(find_file("acs_data.csv"), dtype=str)   # read everything as text first
    acs = pd.DataFrame({
        "Tract_ID": raw["tract_id"].str.strip().str.zfill(11),  # make sure it is 11 digits
        "population": to_num(raw["total_population"]),
        "median_income": to_num(raw["median_income"]),          # "-" becomes NaN
        "poverty_rate": to_num(raw["poverty_rate"]),            # % below poverty line
        "no_vehicle_rate": to_num(raw["no_vehicle_rate"]),      # % households with no car
    }).drop_duplicates("Tract_ID")
    print(f"Loaded {len(acs)} rows. Missing values per column:")
    print(acs.isna().sum().to_string())
    return acs


def merge_acs(tracts, acs):
    """
    Attach the ACS statistics to the tract boundaries (matching on Tract_ID), compute
    population density, and write a plain-language note for any tract with weak data.
    """
    tracts = tracts.merge(acs, on="Tract_ID", how="left")
    tracts["pop_density"] = tracts["population"] / tracts["area_sqmi"]   # people per sq. mile
    print(f"Tracts matched to ACS data: {tracts['population'].notna().sum()} of {len(tracts)}")

    # Build a "data_notes" column such as "median income missing; low population".
    problems = pd.DataFrame({
        "median income missing": tracts["median_income"].isna(),
        "poverty rate missing": tracts["poverty_rate"].isna(),
        "vehicle rate missing": tracts["no_vehicle_rate"].isna(),
        "low population": tracts["population"] < LOW_POP_THRESHOLD,
    })
    tracts["data_notes"] = problems.apply(lambda row: "; ".join(row.index[row]), axis=1)
    return tracts


# =====================================================================================
# STEP 3: FOOD RETAILERS  (grocery stores, convenience stores, farmers markets ...)
# =====================================================================================
# Group the many "Store_Type" labels into a few categories used on the map.
STORE_CATEGORY = {"supermarket": "Full-service grocery", "super store": "Full-service grocery",
                  "grocery store": "Full-service grocery", "farmers market": "Farmers market",
                  "convenience store": "Convenience store", "specialty store": "Specialty store",
                  "other": "Other (dollar/general)"}


def load_stores():
    """Read grocery_stores.csv, clean it, and return one row per unique store."""
    banner("STEP 3: Food retailers",
           "Cleaning store records: removing duplicates and bad coordinates.")
    raw = read_csv_any(find_file("grocery_stores.csv"), dtype={"Tract_ID": str, "Zip_Code": str})
    print(f"Loaded {len(raw)} rows")

    stores = pd.DataFrame({
        "Store_ID": raw["store_id"].astype(str).str.strip(),
        "Name": raw["store_name"].astype(str).str.strip(),
        "Address": raw["store_address"].astype(str).str.strip(),
        "City": raw["city"].astype(str).str.strip(),
        "Zip": raw["Zip_Code"].astype(str).str.strip(),
        "Type": raw["Store_Type"].fillna("Unknown").astype(str).str.strip()
                   .replace({"Farmers and Markets": "Farmers Market"}),   # tidier label
        "lat": to_num(raw["Latitude"]),
        "lon": pd.to_numeric(raw["Longitude"], errors="coerce"),   # longitude is negative here
        "Tract_ID_source": raw["Tract_ID"].astype(str).str.strip(),  # tract listed in the file
    })

    # Remove duplicates: same ID, or same name at the same (standardized) address.
    stores["addr_std"] = stores["Address"].map(std_address)
    before = len(stores)
    stores = stores.drop_duplicates("Store_ID").drop_duplicates(["Name", "addr_std"])
    print(f"Removed {before - len(stores)} duplicate stores")

    # Remove stores whose coordinates are missing or outside the study region.
    bad = ~valid_coords(stores["lat"], stores["lon"])
    if bad.any():
        stores[bad].to_csv(OUT_DIR / "stores_needs_fixing.csv", index=False)
        print(f"{int(bad.sum())} stores have invalid coordinates -> output/stores_needs_fixing.csv")
        stores = stores[~bad]

    # Label each store: is it a full-service grocery, and which map category is it in?
    stores["full_service"] = stores["Type"].str.lower().isin(FULL_SERVICE_TYPES)
    stores["Category"] = (stores["Type"].str.lower().map(STORE_CATEGORY)
                          .fillna("Other (dollar/general)"))

    print("\nStore types found (check that the full-service definition looks right):")
    print(stores.groupby(["Type", "full_service"]).size().to_string())

    # The file has no "open / closed" column, so give the user a checklist instead.
    print("\nNOTE: the file has no operating-status column, so open/closed is NOT verified.")
    stores.loc[stores["full_service"], ["Store_ID", "Name", "Address", "City", "Type"]] \
        .to_csv(OUT_DIR / "full_service_stores_to_verify.csv", index=False)
    print("Wrote output/full_service_stores_to_verify.csv (confirm each store is still open).")
    return stores.reset_index(drop=True)


# =====================================================================================
# STEP 4: COMMUNITY FOOD RESOURCES  (pantries, free meals, food giveaways)
# =====================================================================================
def classify_resource(details, name):
    """
    resources.csv has no "type" column, only free text in 'details'. This reads that text
    and picks ONE main type, checking the most important service first:
        pantry  >  giveaway/distribution  >  meals  >  other
    """
    text = f"{details} {name}".lower()
    if re.search(r"pantry|food bank", text):
        return "Food Pantry"
    if re.search(r"giveaway|distribution|food:|food and|food &|food assistance|"
                 r"emergency|feeding|basket", text):
        return "Food Giveaway"
    if re.search(r"meal|kitchen|breakfast|lunch|dinner|dining|soup|cafe|caf\u00e9", text):
        return "Free Meals"
    return "Other / Call Ahead"


def load_resources():
    """Read resources.csv, clean it, and return one row per resource with a type label."""
    banner("STEP 4: Community food resources",
           "Cleaning pantry / meal / giveaway records and deciding each one's type.")
    raw = read_csv_any(find_file("resources.csv"), dtype={"tract_id": str, "zip_code": str})
    print(f"Loaded {len(raw)} rows")

    resources = pd.DataFrame({
        "Resource_ID": raw["resource_id"].astype(str).str.strip(),
        "Name": raw["resource_name"].astype(str).str.strip(),
        "Address": raw["resource_address"].fillna("Address not listed").astype(str).str.strip(),
        "Phone": raw["phone"].fillna("").astype(str).str.strip(),
        # Some 'details' cells contain line breaks; show them as " | " instead.
        "Details": raw["details"].fillna("").astype(str)
                      .str.replace(r"\s*\n\s*", " | ", regex=True).str.strip(),
        "Requirements": raw["requirements"].fillna("").astype(str).str.strip(),
        "lat": pd.to_numeric(raw["latitude"], errors="coerce"),
        "lon": pd.to_numeric(raw["longitude"], errors="coerce"),
    }).drop_duplicates("Resource_ID")
    resources["Type"] = [classify_resource(d, n)
                         for d, n in zip(resources["Details"], resources["Name"])]

    # A coordinate of 0,0 (or anything outside the region) means the location is unknown.
    bad = ~valid_coords(resources["lat"], resources["lon"])
    if bad.any():
        resources[bad].to_csv(OUT_DIR / "resources_needs_fixing.csv", index=False)
        print(f"{int(bad.sum())} resource(s) lack valid coordinates and were not mapped: "
              f"{', '.join(resources.loc[bad, 'Name'])} -> output/resources_needs_fixing.csv")
        resources = resources[~bad]
    print(resources["Type"].value_counts().to_string())
    return resources.reset_index(drop=True)


# =====================================================================================
# STEP 5: WRTA BUS DATA (GTFS)
# =====================================================================================
def load_gtfs():
    """
    Open the zipped GTFS feed and read the five tables we need:
      stops (where buses stop) | routes (route names/colors) | trips (each scheduled run)
      shapes (the path each run follows) | stop_times (when each run reaches each stop)
    """
    banner("STEP 5: WRTA bus (GTFS) data", "Reading bus stops, routes and schedules.")
    path = find_file("wrta_gfts.zip", "wrta_gfts.zip.zip", "wrta*.zip")
    print("Reading", path)
    # For the big stop_times table we only load the two columns we actually use.
    wanted = {"stops": None, "routes": None, "trips": None, "shapes": None,
              "stop_times": ["trip_id", "stop_id"]}
    tables = {}
    with zipfile.ZipFile(path) as archive:
        inside = {Path(n).name: n for n in archive.namelist()}   # file name -> path in zip
        for table_name, columns in wanted.items():
            file_name = f"{table_name}.txt"
            if file_name not in inside:
                raise FileNotFoundError(f"{file_name} not found inside {path.name}")
            with archive.open(inside[file_name]) as handle:
                tables[table_name] = pd.read_csv(handle, dtype=str, usecols=columns)
    print({name: len(table) for name, table in tables.items()}, "rows per table")
    return tables


def build_transit(gtfs, tracts):
    """
    Turn the raw GTFS tables into map layers.
    Returns:  stops       - one point per bus stop (with a count of scheduled visits)
              routes      - one line per bus route (with its official color)
              stop_routes - a list of which routes serve which stops
    """
    # ---------- STOPS: tidy up, count scheduled visits, and place them on the map ------
    stops = gtfs["stops"].copy()
    stops["lat"] = pd.to_numeric(stops["stop_lat"], errors="coerce")
    stops["lon"] = pd.to_numeric(stops["stop_lon"], errors="coerce")
    stops = stops[valid_coords(stops["lat"], stops["lon"])]
    # "service_events" = how many times a bus is scheduled to stop here (a simple
    # measure of how frequently the stop is served).
    visits = gtfs["stop_times"].groupby("stop_id").size().rename("service_events").reset_index()
    stops = stops.merge(visits, on="stop_id", how="left").fillna({"service_events": 0})
    stops = stops.rename(columns={"stop_id": "Stop_ID", "stop_name": "Stop_Name"})
    stops = make_points(stops[["Stop_ID", "Stop_Name", "service_events", "lat", "lon"]])
    stops = tag_with_tract(stops, tracts)
    print(f"Stops: {len(stops)} total, {int(stops['in_study_area'].sum())} inside Mahoning County")

    # ---------- WHICH ROUTES SERVE WHICH STOP ------------------------------------------
    trips = gtfs["trips"][["trip_id", "route_id"]].drop_duplicates()
    stop_routes = (gtfs["stop_times"][["trip_id", "stop_id"]].drop_duplicates()
                   .merge(trips, on="trip_id")[["stop_id", "route_id"]].drop_duplicates()
                   .rename(columns={"stop_id": "Stop_ID", "route_id": "Route_ID"}))

    # ---------- ROUTES: connect the dots in shapes.txt to draw each route as a line ----
    shapes = gtfs["shapes"].copy()
    shapes["lat"] = pd.to_numeric(shapes["shape_pt_lat"])
    shapes["lon"] = pd.to_numeric(shapes["shape_pt_lon"])
    shapes["order"] = pd.to_numeric(shapes["shape_pt_sequence"])     # order of the dots
    line_rows = []
    for shape_id, dots in shapes.sort_values(["shape_id", "order"]).groupby("shape_id"):
        if len(dots) >= 2:                                            # a line needs 2+ points
            line_rows.append({"shape_id": shape_id,
                              "geometry": LineString(zip(dots["lon"], dots["lat"]))})
    lines = gpd.GeoDataFrame(pd.DataFrame(line_rows), geometry="geometry", crs=CRS_WGS84)

    # A route has several shapes (inbound, outbound, variations); merge them per route.
    shape_to_route = gtfs["trips"][["route_id", "shape_id"]].drop_duplicates()
    lines = lines.merge(shape_to_route, on="shape_id").to_crs(CRS_PROJ)
    routes = lines[["route_id", "geometry"]].dissolve(by="route_id", as_index=False)

    # Add the route name and WRTA's official color.
    names = gtfs["routes"].copy()
    for column in ("route_short_name", "route_long_name", "route_color"):
        if column not in names.columns:
            names[column] = np.nan
    routes = routes.merge(names[["route_id", "route_short_name", "route_long_name",
                                 "route_color"]], on="route_id", how="left")
    routes["Route_Name"] = (routes["route_long_name"].fillna(routes["route_short_name"])
                            .fillna(routes["route_id"]))
    routes["color"] = "#" + routes["route_color"].fillna("2b6cb0").astype(str).str.lstrip("#")
    routes = routes.rename(columns={"route_id": "Route_ID"})
    routes = routes[["Route_ID", "Route_Name", "color", "geometry"]]
    routes["geometry"] = routes.geometry.simplify(40)    # fewer points = smaller, faster map
    print(f"Routes: {len(routes)}")
    return stops, routes, stop_routes


# =====================================================================================
# STEPS 6 and 9: SAVE TO THE "SPATIAL DATABASE" (a GeoPackage file)
# =====================================================================================
def save_layers(gpkg_path, layers):
    """
    Save map layers into one GeoPackage (.gpkg) file. Think of it as a single-file
    database that GIS programs such as ArcGIS Pro and QGIS can open directly.
    `layers` is a dictionary of {layer name: table}.
    """
    for layer_name, table in layers.items():
        table.to_file(gpkg_path, layer=layer_name, driver="GPKG")
    print("Saved layers:", ", ".join(layers), "->", gpkg_path)


# =====================================================================================
# STEP 7: FOOD ACCESS INDICATORS
# =====================================================================================
def nearest(origins, id_column, targets, extra=()):
    """
    For every origin (e.g. a tract center), find the NEAREST target (e.g. a grocery
    store) and return the distance in MILES. `extra` lists target columns to carry along
    (for example, the name of the nearest store).
    """
    target_layer = targets[list(extra) + ["geometry"]]
    joined = gpd.sjoin_nearest(origins[[id_column, "geometry"]], target_layer,
                               how="left", distance_col="_d")
    joined = joined[~joined.index.duplicated(keep="first")]   # if two targets tie, keep one
    result = joined[[id_column, "_d"] + list(extra)].copy()
    result["_d"] = result["_d"] / FT_PER_MILE                 # feet -> miles
    return result


def count_within(origins, id_column, targets, miles, label):
    """
    Count how many targets lie within `miles` of each origin: draw a circle (buffer)
    around each origin and count the targets that fall inside it.
    """
    circles = origins[[id_column, "geometry"]].copy()
    circles["geometry"] = circles.geometry.buffer(miles * FT_PER_MILE)
    marked = targets[["geometry"]].assign(_one=1)             # every target counts as 1
    joined = gpd.sjoin(circles, marked, how="left", predicate="intersects")
    return joined.groupby(id_column)["_one"].sum().fillna(0).astype(int).rename(label).reset_index()


def compute_indicators(tracts, stores, resources, stops, stop_routes):
    """
    Calculate the food-access measures for every tract and add them as new columns.
    Also adds `near_transit` to each store (is it within walking distance of a bus stop?).
    """
    banner("STEP 7: Food access indicators",
           "Measuring distances to groceries, bus stops and pantries for every tract.")

    # Each tract is represented by one point (its Census internal point).
    origins = gpd.GeoDataFrame(
        tracts[["Tract_ID"]].copy(),
        geometry=gpd.points_from_xy(tracts["int_lon"], tracts["int_lat"]),
        crs=CRS_WGS84).to_crs(CRS_PROJ)

    full_service = stores[stores["full_service"]]
    result = tracts[["Tract_ID"]].copy()

    # (a) Nearest FULL-SERVICE grocery (and its name) and nearest food store of ANY kind.
    near = nearest(origins, "Tract_ID", full_service, extra=["Store_ID", "Name"]).rename(
        columns={"_d": "dist_fullservice_mi", "Store_ID": "nearest_store_id",
                 "Name": "nearest_store_name"})
    result = result.merge(near, on="Tract_ID", how="left")
    near = nearest(origins, "Tract_ID", stores).rename(columns={"_d": "dist_anystore_mi"})
    result = result.merge(near, on="Tract_ID", how="left")

    # (b) Nearest WRTA bus stop.
    near = nearest(origins, "Tract_ID", stops).rename(columns={"_d": "dist_stop_mi"})
    result = result.merge(near, on="Tract_ID", how="left")

    # (c) Nearest pantry/giveaway (food to take home) and nearest farmers market.
    take_home = resources[resources["Type"].isin(["Food Pantry", "Food Giveaway"])]
    near = nearest(origins, "Tract_ID", take_home).rename(columns={"_d": "dist_pantry_mi"})
    result = result.merge(near, on="Tract_ID", how="left")
    markets = stores[stores["Type"] == "Farmers Market"]
    if len(markets):
        near = nearest(origins, "Tract_ID", markets).rename(
            columns={"_d": "dist_farmers_market_mi"})
        result = result.merge(near, on="Tract_ID", how="left")

    # (d) How many places are within one mile of the tract center?
    for layer, label in [(full_service, "n_fullservice_1mi"), (stores, "n_stores_1mi"),
                         (resources, "n_resources_1mi")]:
        result = result.merge(count_within(origins, "Tract_ID", layer, RADIUS_MI, label),
                              on="Tract_ID", how="left")

    # (e) TRANSIT ACCESS TO STORES: is each store within walking distance of a bus stop?
    #     Then, for each tract, is its nearest full-service grocery reachable by bus?
    store_stop = nearest(stores, "Store_ID", stops).rename(columns={"_d": "store_to_stop_mi"})
    stores = stores.merge(store_stop, on="Store_ID", how="left")
    stores["near_transit"] = stores["store_to_stop_mi"] <= TRANSIT_WALK_MI
    on_transit = stores.set_index("Store_ID")["near_transit"]
    result["nearest_store_on_transit"] = (result["nearest_store_id"].map(on_transit)
                                          .map({True: "Yes", False: "No"}).fillna("Unknown"))

    # (f) Bus service INSIDE each tract: number of stops, scheduled visits, routes.
    in_county = stops[stops["in_study_area"]]
    per_tract = in_county.groupby("Tract_ID").agg(
        n_stops=("Stop_ID", "size"), service_events=("service_events", "sum")).reset_index()
    routes_per_tract = (in_county[["Stop_ID", "Tract_ID"]].merge(stop_routes, on="Stop_ID")
                        .groupby("Tract_ID")["Route_ID"].nunique().rename("n_routes")
                        .reset_index())
    result = (result.merge(per_tract, on="Tract_ID", how="left")
                    .merge(routes_per_tract, on="Tract_ID", how="left"))
    for column in ("n_stops", "service_events", "n_routes"):
        result[column] = result[column].fillna(0).astype(int)    # no stops = 0, not blank

    tracts = tracts.merge(result, on="Tract_ID", how="left")
    summary_cols = ["dist_fullservice_mi", "dist_anystore_mi", "dist_stop_mi",
                    "dist_pantry_mi", "n_stores_1mi", "n_resources_1mi"]
    print(tracts[summary_cols].describe().round(2).to_string())
    print(f"\nTracts with NO full-service grocery within {RADIUS_MI:g} mile: "
          f"{int((tracts['n_fullservice_1mi'] == 0).sum())} of {len(tracts)}")

    # NOTE on Phase 2 design decision #1: to measure distance along ROADS instead of
    # straight lines, download a road network with the 'osmnx' package and replace
    # nearest() with shortest-path lengths. Compare both and report the difference.
    return tracts, stores


# =====================================================================================
# STEP 8: VULNERABILITY SCORE
# =====================================================================================
def scale_0_to_1(series, low_pct=5, high_pct=95):
    """
    Put a column on a common 0-to-1 scale so different measures can be combined fairly
    (miles, dollars and percentages cannot be added directly).
      * Missing values are filled with the median (a neutral guess).
      * Extreme values are capped at the 5th / 95th percentile so one odd tract
        does not squash everyone else toward zero.
      * The lowest value becomes 0 and the highest becomes 1.
    """
    values = series.astype(float)
    values = values.fillna(values.median())
    low, high = np.nanpercentile(values, [low_pct, high_pct])
    values = values.clip(low, high)
    return (values - low) / (high - low) if high > low else values * 0.0


def score_vulnerability(tracts, weights=WEIGHTS):
    """
    Combine the ingredients into one 0-100 vulnerability score per tract:
        score = 100 x (weight1 x ingredient1 + weight2 x ingredient2 + ...)
    Also ranks the tracts (1 = most vulnerable) and splits them into four equal groups.
    """
    tracts = tracts.copy()
    parts = pd.DataFrame(index=tracts.index)       # the scaled ingredients, one column each
    parts["dist_fullservice_mi"] = scale_0_to_1(tracts["dist_fullservice_mi"])
    parts["no_vehicle_rate"] = scale_0_to_1(tracts["no_vehicle_rate"])
    parts["low_income"] = 1 - scale_0_to_1(tracts["median_income"])   # flipped: LOW income is bad
    parts["poverty_rate"] = scale_0_to_1(tracts["poverty_rate"])
    parts["dist_stop_mi"] = scale_0_to_1(tracts["dist_stop_mi"])

    tracts["vuln_score"] = (sum(parts[name] * w for name, w in weights.items()) * 100).round(1)
    tracts["vuln_rank"] = tracts["vuln_score"].rank(ascending=False, method="min").astype(int)
    # qcut splits the tracts into 4 equally sized groups (quartiles) by score.
    groups = pd.qcut(tracts["vuln_score"].rank(method="first"), 4, labels=CATEGORY_LABELS)
    tracts["vuln_category"] = groups.astype(str)
    return tracts


# =====================================================================================
# STEP 10: INTERACTIVE MAP (the dashboard)
# =====================================================================================
def fmt(value, pattern, missing="n/a"):
    """Format a number for display, e.g. fmt(2.345, '{:.1f} mi') -> '2.3 mi'."""
    return missing if pd.isna(value) else pattern.format(value)


def tract_profiles(tracts_wgs):
    """
    Build the text shown when someone clicks a neighborhood (R2: neighborhood profile).
    Returns the map layer plus the list of fields to display in the pop-up.
    """
    t = tracts_wgs
    profile = pd.DataFrame({
        "Tract_Label": t["Tract_Name"],
        "Vulnerability": (t["vuln_category"] + " (score " + t["vuln_score"].astype(str)
                          + ", rank " + t["vuln_rank"].astype(str) + " of " + str(len(t)) + ")"),
        "Population": t["population"].map(lambda v: fmt(v, "{:,.0f}")),
        "Median income": t["median_income"].map(lambda v: fmt(v, "${:,.0f}")),
        "Poverty rate": t["poverty_rate"].map(lambda v: fmt(v, "{:.1f}%")),
        "Households with no vehicle": t["no_vehicle_rate"].map(lambda v: fmt(v, "{:.1f}%")),
        "Nearest grocery": (t["nearest_store_name"].fillna("n/a") + " ("
                            + t["dist_fullservice_mi"].map(lambda v: fmt(v, "{:.2f} mi")) + ")"),
        "Grocery reachable by bus": t["nearest_store_on_transit"],
        "Nearest bus stop": t["dist_stop_mi"].map(lambda v: fmt(v, "{:.2f} mi")),
        "Nearest pantry": t["dist_pantry_mi"].map(lambda v: fmt(v, "{:.2f} mi")),
        "Food stores within 1 mile": (t["n_stores_1mi"].astype(str) + " ("
                                      + t["n_fullservice_1mi"].astype(str) + " full-service)"),
        "Food resources within 1 mile": t["n_resources_1mi"].astype(str),
        "WRTA routes in tract": t["n_routes"].astype(str),
        "Data notes": t["data_notes"].replace("", "none"),
        "vuln_category": t["vuln_category"],            # used for coloring, hidden in pop-up
    })
    layer = gpd.GeoDataFrame(profile, geometry=t.geometry.values, crs=CRS_WGS84)
    shown_fields = [c for c in profile.columns if c != "vuln_category"]
    return layer, shown_fields


def add_tract_layer(world_map, profile_layer, shown_fields):
    """Color every tract by vulnerability, add click pop-ups and a search box."""
    tract_layer = folium.GeoJson(
        profile_layer, name="Vulnerability by census tract",
        style_function=lambda feature: {
            "fillColor": CATEGORY_COLORS.get(feature["properties"]["vuln_category"], "#999999"),
            "color": "#555555", "weight": 1, "fillOpacity": 0.55},
        highlight_function=lambda feature: {"weight": 3, "color": "black"},   # on mouse-over
        tooltip=folium.GeoJsonTooltip(fields=["Tract_Label", "vuln_category"],
                                      aliases=["Tract", "Vulnerability"]),
        popup=folium.GeoJsonPopup(fields=shown_fields,
                                  aliases=[f + ":" for f in shown_fields], max_width=360),
    ).add_to(world_map)
    plugins.Search(layer=tract_layer, search_label="Tract_Label", collapsed=False,
                   placeholder="Search tract (e.g. Census Tract 8011)",
                   search_zoom=14, position="topright").add_to(world_map)


def add_transit_layers(world_map, routes_wgs, stops_wgs):
    """Draw WRTA routes (in their official colors) and, optionally, every bus stop."""
    folium.GeoJson(
        routes_wgs[["Route_ID", "Route_Name", "color", "geometry"]], name="WRTA routes",
        style_function=lambda feature: {"color": feature["properties"]["color"],
                                        "weight": 3, "opacity": 0.85},
        tooltip=folium.GeoJsonTooltip(fields=["Route_Name"], aliases=["Route"]),
    ).add_to(world_map)

    stop_group = folium.FeatureGroup(name="WRTA stops", show=False)   # hidden until ticked
    for stop in stops_wgs.itertuples():
        folium.CircleMarker([stop.geometry.y, stop.geometry.x], radius=2, color="#333333",
                            fill=True, fill_opacity=0.8,
                            tooltip=str(stop.Stop_Name)).add_to(stop_group)
    stop_group.add_to(world_map)


def add_store_layers(world_map, stores_wgs):
    """One switchable layer per store category. Full-service groceries get bigger dots."""
    colors = {"Full-service grocery": "#1b7837", "Farmers market": "#e08214",
              "Convenience store": "#7b3294", "Specialty store": "#2166ac",
              "Other (dollar/general)": "#999999"}
    shown_by_default = ("Full-service grocery", "Farmers market")
    for category, group in stores_wgs.groupby("Category"):
        layer = folium.FeatureGroup(name=f"Store: {category}", show=category in shown_by_default)
        for store in group.itertuples():
            by_bus = "yes" if store.near_transit else "no"
            popup_text = (f"<b>{store.Name}</b><br>{store.Type}<br>{store.Address}, {store.City}"
                          f"<br>Within {TRANSIT_WALK_MI:g} mi of a bus stop: {by_bus}")
            folium.CircleMarker(
                [store.geometry.y, store.geometry.x],
                radius=6 if store.full_service else 4,
                color="black" if store.full_service else colors.get(category, "#999999"),
                weight=1.5 if store.full_service else 1,
                fill=True, fill_color=colors.get(category, "#999999"), fill_opacity=0.9,
                popup=folium.Popup(popup_text, max_width=280)).add_to(layer)
        layer.add_to(world_map)


def add_resource_layers(world_map, resources_wgs):
    """One switchable layer per resource type (pantry, giveaway, meals, other)."""
    icon_colors = {"Food Pantry": "red", "Food Giveaway": "orange", "Free Meals": "blue",
                   "Other / Call Ahead": "gray"}
    for rtype, group in resources_wgs.groupby("Type"):
        layer = folium.FeatureGroup(name=f"Resource: {rtype}")
        for res in group.itertuples():
            popup_text = (f"<b>{res.Name}</b><br>{rtype}<br>{res.Address}<br>{res.Phone}"
                          f"<br><i>{res.Details}</i>")
            if res.Requirements:
                popup_text += f"<br>Requirements: {res.Requirements}"
            folium.Marker([res.geometry.y, res.geometry.x],
                          icon=folium.Icon(color=icon_colors.get(rtype, "gray"),
                                           icon="cutlery", prefix="glyphicon"),
                          popup=folium.Popup(popup_text, max_width=300)).add_to(layer)
        layer.add_to(world_map)


def add_legend_and_title(world_map, tracts):
    """Add the title, the color legend and a short list of the five highest-priority tracts."""
    top5 = tracts.nsmallest(5, "vuln_rank")
    top5_html = "".join(f"<li>{name} &ndash; {score:.0f}</li>"
                        for name, score in zip(top5["Tract_Name"], top5["vuln_score"]))
    swatches = "".join(f'<div><span style="background:{color};padding:0 10px;'
                       f'margin-right:6px"></span>{label}</div>'
                       for label, color in CATEGORY_COLORS.items())
    world_map.get_root().html.add_child(folium.Element(
        '<div style="position:fixed;bottom:25px;left:15px;z-index:9999;background:white;'
        'padding:10px 14px;border:1px solid #888;border-radius:6px;font-size:13px;'
        f'max-width:260px"><b>Food access vulnerability</b>{swatches}'
        '<hr style="margin:6px 0"><b>Top 5 priority tracts</b>'
        f'<ol style="margin:4px 0 0 18px;padding:0">{top5_html}</ol></div>'))
    world_map.get_root().html.add_child(folium.Element(
        '<h3 style="position:fixed;top:8px;left:55px;z-index:9999;background:white;'
        'padding:4px 10px;border-radius:6px;margin:0;font-family:sans-serif">'
        'Mahoning County / Youngstown Food Access Dashboard</h3>'))


def build_map(tracts, stores, resources, stops, routes, out_file):
    """Assemble all layers into one interactive HTML map and save it."""
    banner("STEP 10: Interactive dashboard",
           "Building the web map with neighborhood profiles, stores, pantries and bus routes.")
    # Web maps need latitude/longitude, so convert every layer back from feet-based CRS.
    tracts_wgs = tracts.to_crs(CRS_WGS84)
    min_x, min_y, max_x, max_y = tracts_wgs.total_bounds        # outline of the county

    # Base maps: OpenStreetMap is the default. The two extra choices appear in the
    # layer menu so the user can switch (none of these need an account or API key).
    world_map = folium.Map(location=[(min_y + max_y) / 2, (min_x + max_x) / 2],
                           zoom_start=11, tiles="OpenStreetMap")
    folium.TileLayer(
        tiles="https://server.arcgisonline.com/ArcGIS/rest/services/Canvas/"
              "World_Light_Gray_Base/MapServer/tile/{z}/{y}/{x}",
        attr="Tiles &copy; Esri &mdash; Esri, DeLorme, NAVTEQ",
        name="Esri Light Gray (low clutter)", max_zoom=16).add_to(world_map)
    folium.TileLayer(
        tiles="https://server.arcgisonline.com/ArcGIS/rest/services/"
              "World_Imagery/MapServer/tile/{z}/{y}/{x}",
        attr="Tiles &copy; Esri &mdash; Source: Esri, Maxar, Earthstar Geographics",
        name="Esri Satellite").add_to(world_map)
    world_map.fit_bounds([[min_y, min_x], [max_y, max_x]])      # zoom to show the whole county

    profile_layer, shown_fields = tract_profiles(tracts_wgs)
    add_tract_layer(world_map, profile_layer, shown_fields)
    add_transit_layers(world_map, routes.to_crs(CRS_WGS84), stops.to_crs(CRS_WGS84))
    add_store_layers(world_map, stores.to_crs(CRS_WGS84))
    add_resource_layers(world_map, resources.to_crs(CRS_WGS84))
    add_legend_and_title(world_map, tracts)
    folium.LayerControl(collapsed=False).add_to(world_map)      # the on/off layer menu
    world_map.save(str(out_file))
    print("Dashboard saved:", out_file)


# =====================================================================================
# MAIN: runs all ten steps in order
# =====================================================================================
def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)         # create the output folder if needed
    print("Data folder:  ", DATA_DIR)
    print("Output folder:", OUT_DIR)

    # ---- STEPS 1-5: load and clean every dataset ------------------------------------
    tracts = load_tracts()                                         # STEP 1
    tracts = merge_acs(tracts, load_acs())                         # STEP 2

    stores = tag_with_tract(make_points(load_stores()), tracts)    # STEP 3
    disagree = (stores["Tract_ID"].notna() &
                (stores["Tract_ID"] != stores["Tract_ID_source"])).sum()
    print(f"Stores inside Mahoning tracts: {int(stores['in_study_area'].sum())} of {len(stores)}"
          " (stores just outside the county are kept so border tracts are measured fairly)")
    print(f"Stores whose CSV Tract_ID disagrees with their map location: {int(disagree)}")

    resources = tag_with_tract(make_points(load_resources()), tracts)   # STEP 4
    stops, routes, stop_routes = build_transit(load_gtfs(), tracts)     # STEP 5

    # ---- STEP 6: save the cleaned layers --------------------------------------------
    banner("STEP 6: Saving cleaned layers", "Writing the cleaned data to a GeoPackage file.")
    gpkg = OUT_DIR / "food_access.gpkg"
    if gpkg.exists():
        gpkg.unlink()                                  # start fresh each run
    save_layers(gpkg, {"tracts_clean": tracts, "food_retailers": stores,
                       "community_resources": resources, "transit_stops": stops,
                       "transit_routes": routes})

    # ---- STEP 7: distances and counts -----------------------------------------------
    tracts, stores = compute_indicators(tracts, stores, resources, stops, stop_routes)

    # ---- STEP 8: the vulnerability score --------------------------------------------
    banner("STEP 8: Vulnerability index",
           "Combining distances, income, poverty and car ownership into one score.")
    tracts = score_vulnerability(tracts)
    print("Tracts per vulnerability group:")
    print(tracts["vuln_category"].value_counts().reindex(CATEGORY_LABELS).to_string())
    print("\nTop 10 most vulnerable tracts:")
    columns = ["Tract_Name", "vuln_score", "dist_fullservice_mi", "no_vehicle_rate",
               "median_income", "poverty_rate"]
    print(tracts.nsmallest(10, "vuln_rank")[columns].round(2).to_string(index=False))

    # ---- STEP 9: save the final results ---------------------------------------------
    banner("STEP 9: Saving results", "Writing scores to the GeoPackage and to a CSV file.")
    save_layers(gpkg, {"tract_vulnerability": tracts, "food_retailers_scored": stores})
    (tracts.drop(columns="geometry").sort_values("vuln_rank")
           .to_csv(OUT_DIR / "tract_vulnerability.csv", index=False))
    print("Wrote output/tract_vulnerability.csv (open it in Excel)")

    # ---- STEP 10: the interactive map -----------------------------------------------
    build_map(tracts, stores, resources, stops, routes, OUT_DIR / "food_access_dashboard.html")
    print("\nDONE. Open output/food_access_dashboard.html in your web browser.")


# This line means: "only run main() when the file is run directly" (not when imported).
if __name__ == "__main__":
    main()


# In[ ]:





# In[ ]:




