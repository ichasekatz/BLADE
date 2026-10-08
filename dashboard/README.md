# BLADE Dashboard

Live monitoring dashboard for BLADE pipeline runs. No base code modified — reads
the same TOML config as `full_framework.py` and watches the output filesystem.

## What it shows

| Panel | Content |
|---|---|
| Stage bar | Pipeline stage status (pending / running / done) with live pulse |
| Systems list | All compositions: queued / running (energy count) / TDB done |
| Center viewer | Live trajectory (`.xyz`) or CONTCAR for selected system |
| Energy chart | Energy convergence across relaxed structures |
| SQS tiles | 10 parallel `mcsqs` runs — structure + correlation score per run |
| Outputs strip | All PNG/GIF outputs as they appear; click to expand |
| Log panel | Last 80 lines of `nohup.out` |

## Setup

### 1. Start the dashboard server on HPC

```bash
cd ~/BLADE/BLADE

# Pass the same TOML you give full_framework.py
python3 dashboard/server.py examples/full_framework.toml --port 8080 --log nohup.out &
```

Multiple config files (same multi-toml syntax as `full_framework.py`):

```bash
python3 dashboard/server.py examples/full_framework.toml examples/configs/borides_hedb.toml --port 8080 --log nohup.out &
```

### 2. Open an SSH tunnel on your laptop

```bash
ssh -L 8080:localhost:8080 ichasekatz@DrALabSpark1
```

### 3. Open the dashboard

[http://localhost:8080](http://localhost:8080)

The page polls every 4 seconds automatically.

## Enable live trajectory

Set `track_trajectory = true` in your TOML's `[tdb_fit]` section (or in the advanced toml).
BLADE writes `.xyz` files into each composition directory during MLIP relaxation,
and the dashboard streams them frame by frame.

## Enable SQS structure preview

`bestsqs<N>.out` files appear automatically in `Files/SQS/<lattice>_<n>/` as `mcsqs -ip=N`
runs. The dashboard converts ATAT format to XYZ and renders each run's current best
structure in 3Dmol. The tile with the lowest correlation objective highlights in teal.

## Stop the server

```bash
pkill -f "dashboard/server.py"
```

## Notes

- No dependencies beyond Python 3.11+ stdlib (uses `tomllib`)
- The server serves only your `Files/` directory — no other paths are accessible
- Images are sent as base64; large GIF files may be slow to load in the strip
