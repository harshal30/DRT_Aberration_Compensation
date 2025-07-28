# DRT_Aberration_Compensation

This repository contains a demo MATLAB implementation of an aberration compensation algorithm for off-axis digital holographic microscopy (DHM). The approach leverages a **dimensionality reduction technique based on Singular Value Decomposition (SVD)** to remove aberrations and spurious fringes from reconstructed phase images.

---

## 📁 Repository Contents

| File | Description |
|------|-------------|
| `funcDRT.m` | Main function that implements the DRT-based aberration removal algorithm. |
| `zernfun2.m` | Utility to simulate Zernike-based aberration terms. |
| `Unwrap_TIE_DCT_Iter.m` | Phase unwrapping function using the Transport-of-Intensity Equation (TIE) with DCT-based iterative refinement. |
| `Demo.m` | Example script to demonstrate aberration removal from a simulated cell as a phase object. |
---

## ▶️ How to Use

1. **Open MATLAB** and navigate to the project folder.

2. **Add all files to the path** and run Demo.m

## 📚 Citation

If you use this code in your research or publication, please cite:
```bibtex
@article{Harshal_2024,
  title = {Dimensionality reduction technique based phase aberration compensation and spurious fringe removal in off-axis digital holographic microscopy},
  journal = {Optics and Lasers in Engineering},
  volume = {172},
  pages = {107853},
  year = {2024},
  issn = {0143-8166},
  doi = {https://doi.org/10.1016/j.optlaseng.2023.107853},
  url = {https://www.sciencedirect.com/science/article/pii/S0143816623003822},
  author = {Harshal Chaudhari and Rishikesh Kulkarni and Pradeep Kumar Sundaravadivelu and Rajkumar P. Thummer and M.K. Bhuyan},
  keywords = {Digital holographic microscopy, Phase aberration, Phase aberration compensation, Dimensionality reduction},
}

