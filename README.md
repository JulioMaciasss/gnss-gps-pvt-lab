# 🛰️ GNSS/GPS PVT Lab — MATLAB

Implementation of a complete GPS/GNSS positioning pipeline from raw RINEX files to a full Single Point Positioning (SPP) solver with error corrections.

Developed as part of the MIT (Ingeniería de Telecomunicaciones) GNSS lab at ETSIT, Universidad de Málaga.

---

## 📌 What this does

Starting from raw RINEX navigation and observation files, this project implements step by step:

1. **RINEX parsing** — reads navigation ephemeris and pseudorange observables
2. **Satellite orbit computation** — computes satellite ECEF position and velocity from broadcast ephemeris
3. **Least-squares PVT solver** — iterative nonlinear LS to estimate user position and clock bias
4. **Error corrections applied:**
   - Satellite clock bias (broadcast + relativistic correction)
   - Sagnac effect (Earth rotation correction)
   - Ionospheric delay (Klobuchar model)
   - Tropospheric delay
5. **DOP computation** — PDOP, TDOP, GDOP from the geometry matrix
6. **KML export** — visualise the position track in Google Earth

---

## 📁 Structure

```
├── Data/                    # RINEX files (nav + obs): parking, UMA lab, dynamic car
├── Functions/               # Helper functions (solveLS, cart2geod, global2localPos, ...)
├── Corrections/             # Ionospheric and tropospheric correction models
├── Rinex/                   # RINEX reader utilities
├── Ex1_Rinex_SatPos.m       # Exercise 1: RINEX reading + satellite position
├── Ex2_PVT.m                # Exercise 2: Full PVT solver
├── satellite_orbits_clock.m
├── solveLS.m
└── output.kml               # Example Google Earth output
```

---

## 📊 Datasets used

| Dataset | Scenario |
|---|---|
| `RinexPar*` | Static — parking lot (open sky) |
| `RinexUMA*` | Static — inside university lab (degraded signal) |
| `RinexCar*` | Dynamic — car trajectory |

---

## 🛠️ Requirements

- MATLAB (tested on R2023b+)
- No additional toolboxes required

---

## 🔑 Key concepts

`RINEX` · `GPS ephemeris` · `pseudorange` · `least-squares` · `Klobuchar ionosphere` · `Sagnac correction` · `ECEF → LLH` · `DOP` · `SPP`

---

## 📎 Related work

This lab is part of broader GNSS performance research. My BSc and MSc theses focus on **OTFS waveform assessment for high-mobility and aerial communications** — a scenario where robust positioning and timing are critical.

→ [OTFS vs OFDM — High Mobility repo](https://github.com/JulioMaciasss/OTFS-OFDM-High-Mobility)
