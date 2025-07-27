# DRT_Aberration_Compensation

This repository contains a demo MATLAB implementation of an aberration compensation algorithm for off-axis digital holographic microscopy (DHM). The approach leverages a **dimensionality reduction technique based on Singular Value Decomposition (SVD)** to remove aberrations and spurious fringes from reconstructed phase images.

---

## 📁 Repository Contents

| File | Description |
|------|-------------|
| `funcDRT.m` | Main function that implements the DRT-based aberration removal algorithm. |
| `zernfun2.m` | Utility to simulate Zernike-based aberration terms. |
| `Unwrap_TIE_DCT_Iter.m` | Phase unwrapping function using the Transport-of-Intensity Equation (TIE) with DCT-based iterative refinement. |
| `main.m` | Example script to demonstrate aberration removal from a simulated holographic phase object. |
| `recorded_hologram.png` | Sample input hologram image (optional – if included). |
| `output/` | Folder to save the corrected phase image output (create if needed). |

---

## ▶️ How to Use

1. **Open MATLAB** and navigate to the project folder.

2. **Add all files to the path**:
   ```matlab
   addpath(genpath(pwd))

