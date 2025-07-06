# Cyclic-iMBV

## Abstract

Detecting intramyocardial blood volume (iMBV) dynamics from end-systole (ES) to end-diastole (ED) may help identify and assess ischemic heart disease. This study aims to demonstrate the feasibility of capturing cyclic ES-to-ED iMBV changes using a novel myocardial “T1-tracking” sequence with ferumoxytol-enhanced (FE) MRI.

To address the confounding effects of propagation of spin history, in-flow effects, and through-plane motion in 2D imaging, a new continuous 2D/3D “hybrid” steady-state sequence was developed, combining slice-exciting (2D) and slab-exciting (3D) RF pulses. Phase-specific myocardial T1 maps were generated using a Bloch equation-based lookup table. T1 values before and after ferumoxytol infusion were used to calculate iMBV in three coronary territories at ES and ED phases. 

## Overview

![Figure_3](https://github.com/user-attachments/assets/cb6210f2-c7b1-47d7-becb-69c1beb98a78)

Conceptual description of the new hybrid 2D/3D sequence and how it maintains the out-of-slice magnetization in steady-state. (a) End-systolic hybrid 2D/3D pulse sequence wherein 3D and 2D RF pulses are applied continuously with golden-angle radial RF-spoiled gradient-recalled echo (SPGR) readouts. The proposed sequence starts with a series of rapid 3D excitation pulses to drive all spins in the volume (covering the whole heart, shown as a blue stack) to steady-state, then switches between 2D and 3D pulses based on a user-defined “2D readout/imaging” interval (to image the 2D short-axis slice, shown in purple) that is synchronized with the R-wave to enable systolic (as shown) or diastolic (not shown) imaging. (b) The magnetization state in the spins outside of the 2D imaging slice is schematically described; a key feature of the proposed hybrid 2D/3D sequence is that the 3D excitation pulses ensure that the out-of-slice magnetization (which is disrupted during 2D excitation) returns to steady state before the next 2D-imaging interval. This in turn ensures robustness to confounding effects of blood in-flow and through-plane motion.

## Code

This repository provides an example code to reconstruct raw k-space data (from an animal study) acquired with hybrid 2D/3D pulse sequence.

- Download the zip file "demo_recon_hybrid2D3D.zip" and unzip
- Add all subfolders to your MATLAB path
- Navigate to the folder "DEMO_recon_hybrid2D3D" and run "demo_recon_Hybrid2D3D.m" 

## Contact

Hazar Benan Unal (hunal@purdue.edu)
