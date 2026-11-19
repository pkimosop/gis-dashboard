# GIS Dashboard Starter

A beginner-friendly Leaflet dashboard with a small local Python data-preparation script.

## Important

The included census polygon and grocery store locations are fictional demonstration data.
Replace them with verified data before using the dashboard for research or decisions.

## Project structure

- `data/input/`: source data
- `src/run_pipeline.py`: local Python processing
- `site/index.html`: dashboard webpage
- `site/style.css`: dashboard styling
- `site/script.js`: Leaflet map and data loading
- `site/data/`: GeoJSON files loaded by the dashboard
- `logs/`: local processing log
- `outputs/`: local processing report

## Run the Python script locally

From the project folder in PowerShell:

```powershell
py src/run_pipeline.py
```

The script uses only Python's standard library. No packages need to be installed.

## View the dashboard locally

From the project folder:

```powershell
py -m http.server 8000 --directory site
```

Open `http://localhost:8000` in your browser.

Press `Ctrl+C` in PowerShell to stop the local web server.

## GitHub

GitHub is used to store and version-control the project files. This starter project does not include a GitHub Actions workflow and does not run Python on GitHub.
