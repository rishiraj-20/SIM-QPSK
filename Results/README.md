# Experimental Results

This folder contains the representative output figures generated from the experimental evaluation of the proposed receiver detection techniques for SIM-QPSK based Discrete-Modulated CV-QKD.

The figures are organized as follows:

| File | Description |
|------|-------------|
| `Circular.png` | Constellation showing the Circular (Radial) detector decision region. |
| `Rectangular.png` | Constellation showing the Rectangular detector decision region. |
| `Angular-Radial.png` | Constellation showing the proposed combined Angular-Radial detector decision region. |
| `All_Detectors.png` | Combined visualization comparing the Circular, Rectangular, and Angular-Radial detector decision regions on the same constellation. |
| `bob.jpeg` | Received constellation at Bob after completion of the entire DSP pipeline, including synchronization, CFO correction, timing recovery, channel estimation, equalization, and phase correction. |
| `eve.jpeg` | Received constellation at Eve after the complete DSP pipeline, illustrating the degraded symbol distribution compared to Bob. |
| `psift_t.jpeg` | Sifting Probability (Psift) obtained by varying the angular threshold (θ) of the Angular-Radial detector while keeping the radial threshold fixed. |
| `qber_t.jpeg` | Quantum Bit Error Rate (QBER) obtained by varying the angular threshold (θ) of the Angular-Radial detector. |
| `secret_fraction_t.jpeg` | Secret Fraction obtained by varying the angular threshold (θ). |
| `skr_t.jpeg` | Secret Key Rate (SKR) obtained by varying the angular threshold (θ). |
| `psift_r.jpeg` | Sifting Probability (Psift) obtained by varying the radial threshold (ρ), i.e., the minimum decision radius (*r<sub>min</sub>*), of the Angular-Radial detector while keeping the angular threshold fixed. |
| `qber_r.jpeg` | Quantum Bit Error Rate (QBER) obtained by varying the radial threshold (ρ). |
| `secret_fraction_r.jpeg` | Secret Fraction obtained by varying the radial threshold (ρ). |
| `skr_r.jpeg` | Secret Key Rate (SKR) obtained by varying the radial threshold (ρ). |

## Notes

- The suffix **`_t`** denotes experiments performed by sweeping the angular threshold **θ** of the proposed Angular-Radial detector.
- The suffix **`_r`** denotes experiments performed by sweeping the radial threshold **ρ** (implemented as the minimum decision radius **r<sub>min</sub>**) while keeping the angular threshold constant.
- All figures were generated using experimentally acquired data from the RF-over-Fiber SIM-QPSK testbed.
