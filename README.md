# A Practical Introduction to Polynomial Optimization using the Moment-SOS Hierarchy 

This repository contains the source code for both notes and experiments used in the talks.
You can view the rendered notes online at <https://lcwllmr.github.io/momsos/>.
In order to reproduce the experiments or play around with them, follow these instructions to set up the development environment.

```bash
git clone https://github.com/lcwllmr/momsos.git
cd momsos

# first option: make sure Python >= 3.12 is installed
python -m venv .venv
source .venv/bin/activate # if on linux with bash
.venv\Scripts\Activate.ps1 # if on windows with powershell
pip install -e .
python code/momsos/experiments/motzkin_minimize.py

# second option: install [uv](https://docs.astral.sh/uv)
uv sync
uv run code/momsos/experiments/motzkin_minimize.py
```

It's best to start by inspecting the experiments in `[code/momsos/experiments](https://github.com/lcwllmr/momsos/tree/main/code/momsos/experiments)` and go from there.
