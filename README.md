# Mini Linear Algebra Library

Welcome to my custom linear algebra library! This repository serves as both a functional tool and a computational exploration of linear algebra concepts.

---

## About Me

I'm **Atharv**. I am currently pursuing my B.Tech in ECE + AI/ML at the **Netaji Subhash University of Technology** (Batch of 2029).

**Core Interests:**

* Scientific Computing
* Machine Learning
* Linear Algebra
* Numerical Methods

**Contact Information:**

* **Email:**  atharv.shrivastava.ug25@nsut.ac.in
* **LinkedIn:** https://www.linkedin.com/in/atharv-shrivastav-535871312?utm_source=share_via&utm_content=profile&utm_medium=member_android

---

## Project Overview

This project aims to build a lightweight, custom linear algebra library from scratch. It provides fundamental vector operations, matrix operations, advanced matrix decompositions, and practical mathematical applications like least squares.

The primary motivation behind this build is to develop a deep, computational understanding of Linear Algebra rather than just theoretical knowledge. Additionally, architecting this hybrid library demonstrates my ability to design, build, and document non-trivial software systems.

---

## Goals

* **Balanced Performance:** The objective is to create a fast, mini linear algebra library. While it is not intended to outpace heavily optimized libraries like NumPy, it provides highly competitive speeds for its scale.
* **Hybrid Architecture:** The system uses a two-language approach. The backend is written in **Julia** to leverage its exceptional execution speed, while the frontend is wrapped in **Python** to maintain an easy, intuitive syntax.
* **Accessibility:** It introduces a custom `Matrix` object to Python. This is specifically designed for students and developers who are comfortable with mathematical paradigms, allowing them to utilize Julia's high-speed processing through Python's highly readable syntax.

---

## Repository Structure

* **`README.md`:** This file. The high-level overview, installation, and progress.
* **`phase_wise_plan/`:** The phase-by-phase development roadmap.
* **`phase[1]_Vectors/`:** Phase 1, vector operations.
* **`phase_2_Matrix/`:** Phase 2, matrix architecture and basic operations.
* **`phase_3_{ref,rref,det,inverse}/`:** Phase 3, row reduction and linear solvers.
* **`phase_4_decompositions/`:** Phase 4, LU, QR, and SVD decompositions.
* **`package_file/`:** The installable `NLA_athsh` package (Python frontend, Julia backend), including the eigensolver module (Phase 5).

Each phase folder has its own `README.md` that lists the functions in that phase, what each one does, and usage notes. Documentation lives inside these READMEs rather than in separate PDF or Word files.

---

## Current Progress

| Phase | Topic | Status |
| :--- | :--- | :--- |
| **Phase 1** | Vector Space Foundations | Done |
| **Phase 2** | Matrix Architectures | Done |
| **Phase 3** | Row Reductions & Linear Solvers | Done |
| **Phase 4** | Matrix Decompositions | In progress: LU done, QR correction in progress, SVD pending |
| **Phase 5** | Spectral Theory & Eigensolvers | Done |
| **Phase 6** | Applied Machine Learning Models | Pending |
