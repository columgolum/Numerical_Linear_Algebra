# Vector Engine (`vec` Class) Documentation

Welcome to the **Vector Engine** module. This subfolder houses the high-performance hybrid architecture engineered for processing $N$-dimensional vector spaces. By combining Python’s elegant syntax with Julia’s lightning-fast compiled computational loops, this engine ensures near-native CPU speeds for vector mathematics.

---

##  Index
1. [User Manual: Detailed Method Breakdown](#1-user-manual-detailed-method-breakdown)
2. [Architectural Evolution History](#2-architectural-evolution-history)

---

## 1. User Manual: Detailed Method Breakdown

An operational overview of the functionalities provided by the latest hybrid Python-Julia `vec` class. 

>  **Prerequisite:** Ensure Python, Julia, and the `juliacall` library are fully configured on your local system before executing these methods.

### Core Class API Reference

| Method | Primary Function | Design & Safety Features / Usage Notes |
| :--- | :--- | :--- |
| **`__init__(*comp)`** | Instantiates an $N$-dimensional vector via optimized Julia arrays. | **Data Security:** Isolates components in a pseudo-private structure (`__comp`) to prevent external tampering. |
| **`len()`** | Returns the exact number of spatial dimensions ($N$). | **Usage Note:** Returns spatial coordinate count. Use `mag()` if you are seeking the geometric Euclidean length. |
| **`get_vec()`** | Returns a direct pointer to the compiled backend Julia array. | **Performance:** Passes memory references directly to ensure zero runtime Python translation overhead. |
| **`__repr__()`** | Provides a structured, cleanly aligned console printout. | **Readability:** Leverages Julia’s native `jl.sprint` for clean matrix-style formatting. |
| **`__add__(p1, p2)`** | Computes element-wise vector addition ($p_1 + p_2$) using optimized CPU threads. | **Validation:** Automatically halts execution with a native `ArgumentError` if vector dimensions mismatch. |
| **`__sub__(p1, p2)`** | Computes component-wise vector subtraction ($p_1 - p_2$). | **Data Protection:** Functional approach—preserves the original source vectors and avoids dangerous in-place corruption. |
| **`scale(v, a=1)`** | Multiplies the entire vector layout by a numerical scalar value $a$. | **Usage Note:** Leaving parameter $a$ empty defaults to `1`, effectively generating a safe clone of the vector. |
| **`mag(v)`** | Calculates the absolute Euclidean spatial distance ($\|v\|$) from the coordinate origin. | **Output Format:** Yields a single high-precision scalar float, making direct comparisons straightforward. |
| **`dot(u, v)`** | Computes the algebraic inner dot product ($u \cdot v$) across two vector fields. | **Validation:** Triggers an explicit structural error if the input
